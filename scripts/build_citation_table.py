#!/usr/bin/env python
"""
Recover the originating publication for every genome assembly in a chrombase
genome list, and emit a citable table.

Genome databases assembled at tree-of-life scale routinely contain thousands of
assemblies, and it is not possible to cite each originating paper in the body of
a manuscript.  The usual fallback -- listing accessions and NCBI submitters in a
supplementary table -- credits institutions but not the people who did the
assembly, and it is invisible to citation indices.

This script closes as much of that gap as the public APIs allow.  For each
accession it collects the NCBI assembly metadata, then tries several routes to
find the paper the assembly came from, scores the candidates, and writes one row
per assembly with a resolved DOI/PMID and a ready-to-paste reference string.

Resolution routes, in order of trust:

  1. ``bioproject``     -- the ``<Publication>`` element of the assembly's
                           BioProject record.  This is asserted by the submitter
                           themselves, so it is treated as authoritative.
  2. ``assembly_name``  -- Europe PMC full-text search for the assembly name
                           (e.g. ``xbOstEdul1.1``).  Very high precision for
                           consortium genome notes, which quote the assembly
                           name verbatim.  Auto-generated NCBI names of the form
                           ``ASM123456v1`` carry no information and are skipped.
  3. ``accession``      -- Europe PMC full-text search for the bare accession
                           (``GCA_021613535``).
  4. ``bioproject_text``-- Europe PMC full-text search for the BioProject
                           accession (``PRJNA680543``).  Lowest precision: papers
                           that *use* an assembly cite the same BioProject as the
                           paper that produced it, so hits are scored down unless
                           other evidence agrees.

Candidates are scored on route, whether the organism name appears in the title,
whether the title looks like a genome announcement, and how close the
publication date is to the assembly release date.  A paper published well before
the assembly was released cannot be its source; a paper published years after it
is more likely a downstream user.

Nothing here invents a citation.  Assemblies with no recoverable publication are
written out with ``confidence=none`` and their submitter retained, so they can at
least be credited at the group level.

Results are cached on disk, so an interrupted run resumes cheaply and repeat
runs cost no API calls.

Example
-------
    python scripts/build_citation_table.py \\
        --genome-list genome_database/genome_list.tsv \\
        --out citation_table.tsv \\
        --bibtex citation_table.bib \\
        --email you@example.org
"""

import argparse
import csv
import html
import json
import os
import re
import sys
import threading
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor

import requests

NCBI_DATASETS_REPORT = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/dataset_report"
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EPMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

USER_AGENT = "chrombase-build_citation_table/1.0 (https://github.com/conchoecia/chrombase)"

# NCBI's own placeholder assembly names. They are derived from the accession and
# never appear in a paper, so searching for them only returns noise.
AUTO_ASSEMBLY_NAME = re.compile(r"^ASM\d+v\d+$", re.IGNORECASE)

# Titles of dedicated genome announcements / genome notes.
# Consortium assembly names built from a Tree of Life ID ("xbOstEdul1.1").
# Papers usually quote the ToLID without the assembly version suffix.
TOLID_ASSEMBLY_NAME = re.compile(r"^([a-z]{2}[A-Z][a-z]{2}[A-Z][a-z]{3}\d+)[.\-]")

GENOME_NOTE_TITLE = re.compile(
    r"(genome sequence of|genome assembly|chromosome[- ]level|chromosome[- ]scale|"
    r"reference genome|draft genome|genome of the|de novo genome|genome report)",
    re.IGNORECASE,
)

OUTPUT_COLUMNS = [
    "assembly_accession",
    "current_accession",
    "organism_name",
    "tax_id",
    "assembly_name",
    "assembly_level",
    "source_database",
    "submitter",
    "submitter_normalized",
    "bioproject_accession",
    "release_date",
    "pub_doi",
    "pub_pmid",
    "pub_pmcid",
    "pub_year",
    "pub_title",
    "pub_journal",
    "pub_authors",
    "citation",
    "evidence_route",
    "confidence",
    "score",
    "n_candidates",
    "notes",
]


# --------------------------------------------------------------------------
# plumbing
# --------------------------------------------------------------------------

class RateLimiter:
    """Crude global rate limiter; NCBI allows 3 req/s, or 10 with an API key."""

    def __init__(self, per_second):
        self.interval = 1.0 / float(per_second)
        self.lock = threading.Lock()
        self.next_ok = 0.0

    def wait(self):
        with self.lock:
            now = time.monotonic()
            if now < self.next_ok:
                time.sleep(self.next_ok - now)
                now = time.monotonic()
            self.next_ok = now + self.interval


class Cache:
    """One JSON file per namespace, flushed periodically and at exit."""

    def __init__(self, directory):
        self.dir = directory
        os.makedirs(directory, exist_ok=True)
        self.data = {}
        self.dirty = set()
        self.lock = threading.Lock()

    def _path(self, namespace):
        return os.path.join(self.dir, f"{namespace}.json")

    def _load(self, namespace):
        if namespace not in self.data:
            try:
                with open(self._path(namespace)) as handle:
                    self.data[namespace] = json.load(handle)
            except (OSError, ValueError):
                self.data[namespace] = {}
        return self.data[namespace]

    def get(self, namespace, key):
        with self.lock:
            return self._load(namespace).get(key)

    def put(self, namespace, key, value):
        with self.lock:
            self._load(namespace)[key] = value
            self.dirty.add(namespace)

    def flush(self):
        with self.lock:
            for namespace in sorted(self.dirty):
                tmp = self._path(namespace) + ".tmp"
                with open(tmp, "w") as handle:
                    json.dump(self.data[namespace], handle)
                os.replace(tmp, self._path(namespace))
            self.dirty.clear()


def clean_text(value):
    """Europe PMC marks up titles with escaped HTML ('&lt;i&gt;Adineta vaga&lt;/i&gt;').

    Unescape first, then strip the tags that unescaping reveals, so the result is
    plain text fit to paste into a table.
    """
    text = html.unescape(value or "")
    text = re.sub(r"<[^>]+>", "", text)
    return re.sub(r"\s+", " ", text).strip()


def make_session():
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})
    # Credentials are attached per-request by ncbi_params(), never set on the
    # session: a session-wide api_key would also be sent to Europe PMC.
    session.ncbi_params = {}
    return session


def ncbi_params(session):
    """Contact details and API key, for NCBI hosts only."""
    return dict(getattr(session, "ncbi_params", {}) or {})


def retry_after(response, attempt=0):
    """Seconds to wait before retrying, honouring the server's own advice.

    NCBI does not always send Retry-After, and once it has started refusing it
    stays unhappy for a while, so the fallback is deliberately patient.
    """
    header = response.headers.get("Retry-After")
    if header:
        try:
            return min(120.0, max(1.0, float(header)))
        except (TypeError, ValueError):
            pass
    if response.status_code == 429:
        return min(120.0, 10.0 * (2 ** attempt))
    return min(30.0, 2.0 * (2 ** attempt))


def request_json(session, url, limiter, method="GET", **kwargs):
    for attempt in range(5):
        limiter.wait()
        try:
            response = session.request(method, url, timeout=120, **kwargs)
            if response.status_code == 429 or response.status_code >= 500:
                time.sleep(retry_after(response, attempt))
                raise requests.RequestException(f"HTTP {response.status_code}")
            response.raise_for_status()
            # NCBI occasionally emits a raw tab inside a JSON string (in the
            # echoed query translation), which strict JSON rejects. The payload
            # is otherwise fine, so parse leniently rather than lose the record.
            return json.loads(response.text, strict=False)
        except (requests.RequestException, ValueError) as exc:
            if attempt == 4:
                sys.stderr.write(f"    give up on {url}: {exc}\n")
                return None
            time.sleep(2 ** attempt)
    return None


def request_text(session, url, limiter, method="GET", **kwargs):
    for attempt in range(5):
        limiter.wait()
        try:
            response = session.request(method, url, timeout=120, **kwargs)
            if response.status_code == 429 or response.status_code >= 500:
                time.sleep(retry_after(response, attempt))
                raise requests.RequestException(f"HTTP {response.status_code}")
            response.raise_for_status()
            return response.text
        except requests.RequestException as exc:
            if attempt == 4:
                sys.stderr.write(f"    give up on {url}: {exc}\n")
                return None
            time.sleep(2 ** attempt)
    return None


# --------------------------------------------------------------------------
# submitter names
# --------------------------------------------------------------------------

SUBMITTER_ALIASES = {
    "wellcome trust sanger institute": "Wellcome Sanger Institute",
    "wellcome sanger institute": "Wellcome Sanger Institute",
    "sanger institute": "Wellcome Sanger Institute",
    "genoscope": "Genoscope (CEA)",
    "genoscope cea": "Genoscope (CEA)",
    "genoscope - centre national de sequencage": "Genoscope (CEA)",
    "vertebrate genomes project": "Vertebrate Genomes Project",
    "the vertebrate genomes project": "Vertebrate Genomes Project",
    "vgp": "Vertebrate Genomes Project",
    "dna zoo": "DNA Zoo",
    "dnazoo": "DNA Zoo",
    "bgi": "BGI",
    "usda": "United States Department of Agriculture",
    "united states department of agriculture": "United States Department of Agriculture",
    "usda-ars": "United States Department of Agriculture",
}


# Words that stay lowercase inside a name, unless they lead it.
NAME_MINOR_WORDS = {"of", "the", "and", "for", "at", "in", "on", "de", "del", "della",
                    "di", "da", "do", "dos", "des", "du", "la", "le", "les", "van",
                    "von", "der", "den", "y", "e"}


def _titlecase_token(word, index):
    lowered = word.lower()
    if lowered in NAME_MINOR_WORDS:
        return lowered if index > 0 else lowered.capitalize()
    letters = re.sub(r"[^A-Za-z]", "", word)
    if word.isupper() and len(letters) <= 4:
        # Short and capitalised: probably an acronym (CEA, UPF, SC, USDA).
        return word
    return lowered[:1].upper() + lowered[1:]


def normalize_submitter(name):
    """Collapse the spelling variants NCBI accumulates for one institution.

    NCBI stores the submitter as free text, so a single institution appears
    under several strings ('WELLCOME SANGER INSTITUTE', 'Wellcome Sanger
    Institute').  Counting credit per institution requires collapsing these.
    """
    if not name:
        return ""
    key = re.sub(r"[^a-z0-9 ]+", " ", name.lower())
    key = re.sub(r"\s+", " ", key).strip()
    if key in SUBMITTER_ALIASES:
        return SUBMITTER_ALIASES[key]
    # No alias. NCBI stores submitters as free text, so the same institution can
    # arrive shouting ("WELLCOME SANGER INSTITUTE") or whispering ("university of
    # cologne"). Case-fold those two extremes to one spelling so credit counts
    # group correctly, and leave anything already mixed-case untouched -- it is
    # the submitter's own capitalisation, acronyms and all.
    stripped = name.strip()
    if not (stripped.isupper() or stripped.islower()):
        return stripped
    return " ".join(_titlecase_token(word, index)
                    for index, word in enumerate(stripped.split()))


# --------------------------------------------------------------------------
# stage 1: assembly metadata
# --------------------------------------------------------------------------

def read_genome_list(path):
    """Accept a chrombase genome_list.tsv, an NCBI datasets TSV, or a bare list."""
    accessions = []
    with open(path, newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        if "\t" in sample:
            reader = csv.DictReader(handle, delimiter="\t")
            field = None
            for candidate in ("accession", "Assembly Accession", "assembly_accession"):
                if reader.fieldnames and candidate in reader.fieldnames:
                    field = candidate
                    break
            if field is None:
                handle.seek(0)
                return [line.split("\t")[0].strip() for line in handle
                        if line.strip() and not line.startswith("#")]
            for row in reader:
                value = (row.get(field) or "").strip()
                if value:
                    accessions.append(value)
        else:
            accessions = [line.strip() for line in handle
                          if line.strip() and not line.startswith("#")]
    seen = set()
    unique = []
    for accession in accessions:
        if accession not in seen:
            seen.add(accession)
            unique.append(accession)
    return unique


def fetch_assembly_reports(session, limiter, cache, accessions, chunk_size=100):
    """Bulk-fetch NCBI Datasets assembly reports, one POST per chunk."""
    reports = {}
    todo = []
    for accession in accessions:
        cached = cache.get("assembly", accession)
        if cached is None:
            todo.append(accession)
        else:
            reports[accession] = cached

    for start in range(0, len(todo), chunk_size):
        chunk = todo[start:start + chunk_size]
        payload = {"accessions": chunk, "returned_content": "ASSEMBLY_ONLY", "page_size": 1000}
        result = request_json(session, NCBI_DATASETS_REPORT, limiter, method="POST",
                              json=payload, params=ncbi_params(session))
        returned = set()
        for report in (result or {}).get("reports", []):
            for key in {report.get("accession"), report.get("current_accession")}:
                if key in chunk:
                    reports[key] = report
                    cache.put("assembly", key, report)
                    returned.add(key)
        # Accessions NCBI no longer serves (suppressed, replaced, or withdrawn).
        for accession in chunk:
            if accession not in returned:
                reports[accession] = {}
                cache.put("assembly", accession, {})
        sys.stderr.write(f"  assembly metadata {min(start + chunk_size, len(todo))}/{len(todo)}\n")
        cache.flush()
    return reports


def flatten_report(accession, report):
    info = (report or {}).get("assembly_info", {}) or {}
    organism = (report or {}).get("organism", {}) or {}
    return {
        "assembly_accession": accession,
        "current_accession": report.get("current_accession", "") if report else "",
        "organism_name": organism.get("organism_name", ""),
        "tax_id": organism.get("tax_id", ""),
        "assembly_name": info.get("assembly_name", ""),
        "assembly_level": info.get("assembly_level", ""),
        "source_database": (report or {}).get("source_database", "").replace("SOURCE_DATABASE_", ""),
        "submitter": info.get("submitter", ""),
        "submitter_normalized": normalize_submitter(info.get("submitter", "")),
        "bioproject_accession": info.get("bioproject_accession", ""),
        "release_date": info.get("release_date", ""),
    }


# --------------------------------------------------------------------------
# stage 2: candidate publications
# --------------------------------------------------------------------------

def prefetch_bioproject_publications(session, limiter, cache, accessions, chunk_size=100):
    """Resolve many BioProject accessions to their publications in a few calls.

    E-utilities `id=` takes a numeric UID, not an accession. Passing an accession
    does not fail loudly: NCBI reads the digits out of "PRJEB90089" and returns
    UID 90089, an unrelated project, whose publications would then be credited to
    the wrong assembly. So each accession has to be resolved to a UID first.

    Done one accession at a time that is two requests each, which for a few
    thousand projects is enough traffic to get rate-limited. esearch accepts a
    disjunction of accessions and efetch accepts many UIDs at once, so a whole
    chunk costs two requests instead of two hundred. Records are matched back to
    the accession they report as their own, and an accession missing from the
    response is left uncached so the next run retries it.
    """
    todo = [accession for accession in dict.fromkeys(accessions)
            if accession and cache.get("bioproject", accession) is None]
    if not todo:
        return
    sys.stderr.write(f"  {len(todo)} BioProjects to look up\n")

    for start in range(0, len(todo), chunk_size):
        chunk = todo[start:start + chunk_size]
        term = " OR ".join(f"{accession}[Project Accession]" for accession in chunk)
        params = dict(ncbi_params(session), db="bioproject", term=term,
                      retmode="json", retmax=str(chunk_size * 2))
        found = request_json(session, f"{NCBI_EUTILS}/esearch.fcgi", limiter,
                             method="POST", data=params)
        uids = []
        try:
            uids = found["esearchresult"]["idlist"]
        except (KeyError, TypeError):
            uids = []
        if not uids:
            continue

        params = dict(ncbi_params(session), db="bioproject", id=",".join(uids), retmode="xml")
        text = request_text(session, f"{NCBI_EUTILS}/efetch.fcgi", limiter,
                            method="POST", data=params)
        if text is None:
            continue
        try:
            root = ET.fromstring(text)
        except ET.ParseError:
            continue

        for summary in root.iter("DocumentSummary"):
            own = [node.get("accession") for node in summary.iter("ArchiveID")
                   if node.get("accession")]
            publications = []
            for node in summary.iter("Publication"):
                identifier = (node.get("id") or "").strip()
                if not identifier:
                    continue
                if identifier.startswith("10.") or "/" in identifier:
                    publications.append({"doi": identifier, "pmid": ""})
                elif identifier.isdigit():
                    publications.append({"doi": "", "pmid": identifier})
            # Trust a record only for the accession it reports as its own.
            for accession in own:
                if accession in chunk:
                    cache.put("bioproject", accession, publications)

        sys.stderr.write(f"  BioProjects {min(start + chunk_size, len(todo))}/{len(todo)}\n")
        cache.flush()


def bioproject_publications(session, limiter, cache, bioproject):
    """Publications the submitter attached to the BioProject record.

    Populated by prefetch_bioproject_publications(); a miss here means the
    lookup did not come back, and is left uncached so a later run retries it.
    """
    if not bioproject:
        return []
    return cache.get("bioproject", bioproject) or []


def epmc_search(session, limiter, cache, query, page_size=10):
    cached = cache.get("epmc", query)
    if cached is not None:
        return cached
    params = {"query": query, "format": "json", "pageSize": page_size, "resultType": "lite"}
    result = request_json(session, EPMC_SEARCH, limiter, params=params)
    hits = []
    if result:
        for hit in result.get("resultList", {}).get("result", []):
            hits.append({
                "pmid": hit.get("pmid", ""),
                "pmcid": hit.get("pmcid", ""),
                "doi": hit.get("doi", ""),
                "title": clean_text(hit.get("title", "")),
                "journal": hit.get("journalTitle", "") or "",
                "authors": hit.get("authorString", "") or "",
                "year": hit.get("pubYear", "") or "",
                "date": hit.get("firstPublicationDate", "") or "",
                "source": hit.get("source", ""),
                "id": hit.get("id", ""),
            })
    cache.put("epmc", query, hits)
    return hits


def lookup_publication(session, limiter, cache, doi=None, pmid=None):
    """Turn a bare DOI or PMID into full bibliographic metadata."""
    if pmid:
        hits = epmc_search(session, limiter, cache, f"EXT_ID:{pmid} AND SRC:MED", page_size=1)
    elif doi:
        hits = epmc_search(session, limiter, cache, f'DOI:"{doi}"', page_size=1)
    else:
        return None
    if hits:
        return hits[0]
    if doi:
        return {"pmid": pmid or "", "pmcid": "", "doi": doi, "title": "", "journal": "",
                "authors": "", "year": "", "date": "", "source": "", "id": ""}
    return None


# --------------------------------------------------------------------------
# stage 3: scoring
# --------------------------------------------------------------------------

ROUTE_BASE = {
    "bioproject": 10,        # submitter-asserted; authoritative
    "assembly_name": 4,
    "accession": 4,
    "bioproject_text": 2,
    "organism_title": 1,
    "organism_genus_title": 1,
}

# Routes that identify the assembly itself. `organism_title` only identifies the
# *species*, so it can never be more than circumstantial: if two groups each
# published an assembly for the same species, it cannot tell them apart.
MAX_CONFIDENCE = {"organism_title": "medium", "organism_genus_title": "low"}
CONFIDENCE_ORDER = ["none", "low", "medium", "high", "authoritative"]


def release_year(row):
    match = re.match(r"(\d{4})", row.get("release_date", "") or "")
    return int(match.group(1)) if match else None


def candidate_year(candidate):
    match = re.match(r"(\d{4})", candidate.get("year") or candidate.get("date") or "")
    return int(match.group(1)) if match else None


def score_candidate(candidate, route, row, sole_hit=False):
    """Score how likely `candidate` is the paper that produced this assembly."""
    score = ROUTE_BASE.get(route, 0)
    reasons = [route]

    title = (candidate.get("title") or "").lower()
    organism = (row.get("organism_name") or "").strip()
    parts = organism.split()
    genus = parts[0].lower() if parts else ""
    binomial = " ".join(parts[:2]).lower() if len(parts) >= 2 else ""

    if binomial and binomial in title:
        score += 3
        reasons.append("binomial_in_title")
    elif genus and genus in title:
        score += 2
        reasons.append("genus_in_title")

    if GENOME_NOTE_TITLE.search(title):
        score += 2
        reasons.append("genome_note_title")

    # A BioProject that exactly one paper in the literature mentions is very
    # likely the paper that created it, rather than one of many reusing it.
    if sole_hit:
        score += 3
        reasons.append("sole_bioproject_hit")

    # Prefer the peer-reviewed version of a work over its preprint. Preprints are
    # still kept, and are often the only public record of a consortium assembly.
    if candidate.get("source") == "MED":
        score += 1
        reasons.append("peer_reviewed")

    published = candidate_year(candidate)
    released = release_year(row)
    if published and released:
        delta = published - released
        if -1 <= delta <= 2:
            score += 1
            reasons.append("date_consistent")
        elif delta < -1:
            # The paper predates the assembly's public release: it cannot be
            # describing this assembly.
            score -= 5
            reasons.append("predates_release")
        elif delta > 3:
            # Long after release; more likely a study reusing the assembly.
            score -= 3
            reasons.append("long_after_release")

    return score, reasons


def confidence_for(score, route, ambiguous=False):
    if route == "bioproject":
        return "authoritative"
    if score >= 8:
        level = "high"
    elif score >= 5:
        level = "medium"
    elif score >= 3:
        level = "low"
    else:
        return "none"

    cap = MAX_CONFIDENCE.get(route)
    if cap and CONFIDENCE_ORDER.index(level) > CONFIDENCE_ORDER.index(cap):
        level = cap
    if ambiguous and route in ("organism_title", "organism_genus_title"):
        # Several assemblies of this species, from different groups, are in the
        # input: a species-level match cannot say which one this paper describes.
        level = "low"
    return level


def format_citation(pub):
    """Compact plain-text reference, safe to paste into a supplementary table."""
    bits = []
    if pub.get("authors"):
        bits.append(pub["authors"].rstrip("."))
    if pub.get("year"):
        bits.append(f"({pub['year']})")
    if pub.get("title"):
        bits.append(pub["title"].rstrip(".") + ".")
    if pub.get("journal"):
        bits.append(pub["journal"].rstrip(".") + ".")
    elif pub.get("source") == "PPR":
        bits.append("Preprint.")
    if pub.get("doi"):
        bits.append(f"https://doi.org/{pub['doi']}")
    elif pub.get("pmid"):
        bits.append(f"PMID:{pub['pmid']}")
    return " ".join(bits).strip()


def organism_title_query(organism):
    """Species-level fallback: a genome paper whose title names this species."""
    parts = (organism or "").split()
    if len(parts) < 2 or not parts[1][:1].islower():
        return None
    binomial = " ".join(parts[:2])
    return (f'TITLE:"{binomial}" AND (TITLE:genome OR TITLE:assembly '
            f'OR TITLE:chromosome OR TITLE:genomic)')


def organism_genus_query(organism):
    """Weak fallback for species whose name differs between NCBI and the paper."""
    parts = (organism or "").split()
    if not parts:
        return None
    return (f'TITLE:"{parts[0]}" AND (TITLE:genome OR TITLE:assembly '
            f'OR TITLE:chromosome)')


def resolve_one(session, limiter_ncbi, limiter_epmc, cache, row, ambiguous_species=frozenset()):
    """Find the best originating-publication candidate for one assembly."""
    candidates = []

    for publication in bioproject_publications(session, limiter_ncbi, cache,
                                               row.get("bioproject_accession", "")):
        full = lookup_publication(session, limiter_epmc, cache,
                                  doi=publication.get("doi"), pmid=publication.get("pmid"))
        if full:
            if publication.get("doi") and not full.get("doi"):
                full["doi"] = publication["doi"]
            candidates.append((full, "bioproject", False))

    name = row.get("assembly_name", "")
    if name and not AUTO_ASSEMBLY_NAME.match(name):
        queries = [f'"{name}"']
        tolid = TOLID_ASSEMBLY_NAME.match(name)
        if tolid:
            queries.append(f'"{tolid.group(1)}"')
        for query in queries:
            for hit in epmc_search(session, limiter_epmc, cache, query, page_size=5):
                candidates.append((hit, "assembly_name", False))

    accession = (row.get("assembly_accession") or "").split(".")[0]
    if accession:
        for hit in epmc_search(session, limiter_epmc, cache, f'"{accession}"', page_size=5):
            candidates.append((hit, "accession", False))

    bioproject = row.get("bioproject_accession", "")
    if bioproject:
        hits = epmc_search(session, limiter_epmc, cache, f'"{bioproject}"', page_size=5)
        for hit in hits:
            candidates.append((hit, "bioproject_text", len(hits) == 1))

    query = organism_title_query(row.get("organism_name", ""))
    if query:
        hits = epmc_search(session, limiter_epmc, cache, query, page_size=5)
        for hit in hits:
            candidates.append((hit, "organism_title", False))
        if not hits:
            # Species-level miss is often just a naming mismatch between NCBI and
            # the paper ("Knipowitschia caucasica" vs "Knipowitschia cf.
            # caucasica"). Retry on the genus, flagged as weak evidence.
            genus_query = organism_genus_query(row.get("organism_name", ""))
            if genus_query:
                for hit in epmc_search(session, limiter_epmc, cache, genus_query, page_size=5):
                    candidates.append((hit, "organism_genus_title", False))

    # Collapse duplicates across routes, keeping the highest-scoring route.
    best_by_id = {}
    for candidate, route, sole_hit in candidates:
        key = candidate.get("doi") or candidate.get("pmid") or candidate.get("id") or candidate.get("title")
        if not key:
            continue
        score, reasons = score_candidate(candidate, route, row, sole_hit=sole_hit)
        if key not in best_by_id or score > best_by_id[key][1]:
            best_by_id[key] = (candidate, score, route, reasons)

    row["n_candidates"] = len(best_by_id)
    if not best_by_id:
        row["confidence"] = "none"
        row["score"] = ""
        row["evidence_route"] = ""
        row["notes"] = "no candidate publication found"
        return row

    # Highest score wins; ties go to the earliest paper, since the assembly's own
    # description precedes the studies that reuse it.
    candidate, score, route, reasons = min(
        best_by_id.values(),
        key=lambda item: (-item[1], candidate_year(item[0]) or 9999))

    ambiguous = row.get("organism_name", "") in ambiguous_species
    confidence = confidence_for(score, route, ambiguous=ambiguous)
    if ambiguous and route in ("organism_title", "organism_genus_title"):
        reasons.append("ambiguous_multiassembly_species")
    if confidence == "none":
        row["confidence"] = "none"
        row["score"] = score
        row["evidence_route"] = ""
        row["notes"] = "candidates found but none scored above threshold"
        return row

    candidate = dict(candidate, title=clean_text(candidate.get("title", "")))
    row.update({
        "pub_doi": candidate.get("doi", ""),
        "pub_pmid": candidate.get("pmid", ""),
        "pub_pmcid": candidate.get("pmcid", ""),
        "pub_year": candidate.get("year", ""),
        "pub_title": candidate.get("title", ""),
        "pub_journal": candidate.get("journal", "") or ("preprint" if candidate.get("source") == "PPR" else ""),
        "pub_authors": candidate.get("authors", ""),
        "citation": format_citation(candidate),
        "evidence_route": route,
        "confidence": confidence,
        "score": score,
        "notes": ",".join(reasons),
    })
    return row


# --------------------------------------------------------------------------
# output
# --------------------------------------------------------------------------

def bibtex_key(row, used):
    author = (row.get("pub_authors") or "anon").split(",")[0].split()[0]
    author = re.sub(r"[^A-Za-z]", "", author) or "anon"
    year = row.get("pub_year") or "0000"
    base = f"{author}{year}"
    key = base
    suffix = 0
    while key in used:
        suffix += 1
        key = f"{base}{chr(ord('a') + suffix - 1)}"
    used.add(key)
    return key


BIBTEX_ESCAPE = {"&": r"\&", "%": r"\%", "#": r"\#", "$": r"\$", "_": r"\_"}


def bibtex_escape(value):
    """Escape the characters LaTeX treats specially, and drop stray braces."""
    value = (value or "").replace("{", "").replace("}", "")
    return "".join(BIBTEX_ESCAPE.get(character, character) for character in value)


def bibtex_authors(author_string):
    """Europe PMC returns 'Simion P, Narayan J, ...'; BibTeX wants ' and ' between authors."""
    cleaned = (author_string or "").strip().rstrip(".")
    names = [name.strip() for name in cleaned.split(",") if name.strip()]
    return " and ".join(names)


def write_bibtex(rows, path):
    """One entry per distinct publication, noting the assemblies it accounts for."""
    entries = {}
    order = []
    for row in rows:
        identifier = row.get("pub_doi") or row.get("pub_pmid")
        if not identifier:
            continue
        if identifier not in entries:
            entries[identifier] = (row, [])
            order.append(identifier)
        entries[identifier][1].append(row["assembly_accession"])

    used = set()
    with open(path, "w") as handle:
        for identifier in order:
            row, accessions = entries[identifier]
            key = bibtex_key(row, used)
            handle.write(f"@article{{{key},\n")
            for field, value in (("title", bibtex_escape(row.get("pub_title"))),
                                 ("author", bibtex_escape(bibtex_authors(row.get("pub_authors")))),
                                 ("journal", bibtex_escape(row.get("pub_journal"))),
                                 ("year", row.get("pub_year")),
                                 ("doi", row.get("pub_doi"))):
                if value:
                    handle.write(f"  {field} = {{{value}}},\n")
            shown = ", ".join(accessions[:5])
            if len(accessions) > 5:
                shown += f", and {len(accessions) - 5} more"
            label = "Assembly" if len(accessions) == 1 else f"{len(accessions)} assemblies"
            handle.write(f"  note = {{{label}: {shown}}},\n")
            handle.write("}\n\n")
    return len(entries)


def write_unlinked_report(rows, path):
    """Assemblies whose paper we found, but which NCBI does not link to it.

    These are the actionable cases: the publication exists and is public, but the
    BioProject record carries no <Publication> element, so no automated tool --
    including this one, on a first pass -- can connect the two. Submitters can fix
    this by adding the publication to their own BioProject record.
    """
    columns = ["assembly_accession", "organism_name", "submitter", "bioproject_accession",
               "pub_doi", "pub_pmid", "pub_title", "confidence", "evidence_route"]
    count = 0
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            found = row.get("pub_doi") or row.get("pub_pmid")
            if found and row.get("evidence_route") != "bioproject":
                writer.writerow(row)
                count += 1
    return count


def summarize(rows, handle=sys.stderr):
    from collections import Counter
    confidence = Counter(row.get("confidence") or "none" for row in rows)
    route = Counter(row.get("evidence_route") or "-" for row in rows)
    resolved = sum(1 for row in rows if row.get("pub_doi") or row.get("pub_pmid"))
    handle.write("\n=== summary ===\n")
    handle.write(f"assemblies:            {len(rows)}\n")
    handle.write(f"with a publication:    {resolved} ({100.0 * resolved / max(len(rows), 1):.1f}%)\n")
    handle.write("by confidence:\n")
    for level in ("authoritative", "high", "medium", "low", "none"):
        if confidence.get(level):
            handle.write(f"  {level:<14} {confidence[level]}\n")
    handle.write("by route:\n")
    for name, count in route.most_common():
        handle.write(f"  {name:<16} {count}\n")
    groups = {row.get("submitter_normalized") for row in rows if row.get("submitter_normalized")}
    handle.write(f"distinct submitting groups: {len(groups)}\n")


# --------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Recover originating publications for genome assemblies and write a citable table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-g", "--genome-list", required=True,
                        help="chrombase genome_list.tsv, an NCBI datasets TSV, or a file of accessions.")
    parser.add_argument("-o", "--out", required=True, help="Output TSV.")
    parser.add_argument("-b", "--bibtex", help="Also write a deduplicated BibTeX file.")
    parser.add_argument("-u", "--unlinked-report",
                        help="Also write a TSV of assemblies whose publication we recovered but "
                             "whose BioProject record does not link to it.")
    parser.add_argument("-c", "--cache-dir", default=".citation_cache",
                        help="Directory for cached API responses.")
    parser.add_argument("-e", "--email", default=os.environ.get("NCBI_EMAIL", ""),
                        help="Contact email sent to NCBI, as their usage policy asks.")
    parser.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY", ""),
                        help="NCBI API key; raises the E-utilities rate limit from 3/s to 10/s.")
    parser.add_argument("-t", "--threads", type=int, default=4,
                        help="Concurrent publication lookups.")
    parser.add_argument("-n", "--limit", type=int, help="Only process the first N accessions (for testing).")
    return parser.parse_args()


def main():
    args = parse_args()

    session = make_session()
    if args.email:
        session.ncbi_params["email"] = args.email
    if args.api_key:
        session.ncbi_params["api_key"] = args.api_key

    limiter_ncbi = RateLimiter(9 if args.api_key else 1.5)
    limiter_epmc = RateLimiter(6)
    cache = Cache(args.cache_dir)

    accessions = read_genome_list(args.genome_list)
    if args.limit:
        accessions = accessions[:args.limit]
    sys.stderr.write(f"{len(accessions)} accessions from {args.genome_list}\n")

    sys.stderr.write("stage 1: assembly metadata from NCBI Datasets\n")
    reports = fetch_assembly_reports(session, limiter_ncbi, cache, accessions)

    rows = []
    for accession in accessions:
        row = {column: "" for column in OUTPUT_COLUMNS}
        row.update(flatten_report(accession, reports.get(accession) or {}))
        if not row["organism_name"]:
            row["notes"] = "no NCBI assembly report (suppressed, replaced, or withdrawn)"
            row["confidence"] = "none"
        rows.append(row)

    sys.stderr.write("stage 2: BioProject publication links\n")
    prefetch_bioproject_publications(
        session, limiter_ncbi, cache,
        [row["bioproject_accession"] for row in rows if row["bioproject_accession"]])

    sys.stderr.write("stage 3: publication recovery\n")
    live = [row for row in rows if row["organism_name"]]
    done = [0]
    lock = threading.Lock()

    # Species represented by more than one assembly from different groups: a
    # species-level literature match cannot be attributed to a single assembly.
    by_species = {}
    for row in live:
        by_species.setdefault(row["organism_name"], set()).add(row["submitter_normalized"])
    ambiguous_species = frozenset(name for name, groups in by_species.items() if len(groups) > 1)
    sys.stderr.write(f"  {len(ambiguous_species)} species have assemblies from more than one group\n")

    def work(row):
        try:
            resolve_one(session, limiter_ncbi, limiter_epmc, cache, row,
                        ambiguous_species=ambiguous_species)
        except Exception as exc:                       # keep one bad row from killing the run
            row["notes"] = f"error: {exc}"
            row["confidence"] = "none"
        with lock:
            done[0] += 1
            if done[0] % 100 == 0:
                sys.stderr.write(f"  {done[0]}/{len(live)}\n")
                cache.flush()

    with ThreadPoolExecutor(max_workers=max(1, args.threads)) as pool:
        list(pool.map(work, live))
    cache.flush()

    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    sys.stderr.write(f"\nwrote {args.out}\n")

    if args.bibtex:
        count = write_bibtex(rows, args.bibtex)
        sys.stderr.write(f"wrote {args.bibtex} ({count} unique references)\n")

    if args.unlinked_report:
        count = write_unlinked_report(rows, args.unlinked_report)
        sys.stderr.write(f"wrote {args.unlinked_report} ({count} assemblies NCBI does not link)\n")

    summarize(rows)


if __name__ == "__main__":
    main()
