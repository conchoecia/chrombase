#!/usr/bin/env python3
"""Build a copyright-minimized public release from local validation outputs.

The release stores each short attributed evidence fragment once. Assembly
claims reference fragment IDs instead of repeating quotations. Full text,
abstracts, local paths, email addresses, and correspondence are never copied.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import re
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from typing import Iterable


MAX_EXCERPT_TOKENS = 25
SPACE_RE = re.compile(r"\s+")
TAG_RE = re.compile(r"<[^>]+>")
PMC_RE = re.compile(r"(PMC\d+)", re.I)
EMAIL_RE = re.compile(r"[A-Z0-9._%+\-]+@[A-Z0-9.\-]+\.[A-Z]{2,}", re.I)

PAPER_FIELDS = [
    "paper_id", "primary_doi", "primary_pmid", "primary_title", "pub_year",
    "pub_journal", "pub_authors", "n_assembly_claims", "paper_overall_status",
    "best_support_level", "worst_support_level", "n_exact_supported",
    "n_composite_supported", "n_unresolved", "n_downstream_flagged",
    "source_type", "source_document_id", "source_url", "source_license",
]

FRAGMENT_FIELDS = [
    "fragment_id", "paper_id", "excerpt_type", "excerpt", "locator",
    "source_block_sha256", "source_document_id", "source_url",
    "source_license", "excerpt_treatment", "rights_basis",
]

CLAIM_FIELDS = [
    "assembly_accession", "organism_name", "assembly_name", "paper_id",
    "primary_doi", "primary_pmid", "text_evidence_status", "support_level",
    "text_supports_current_primary", "production_fragment_id",
    "identity_fragment_id", "identifier_type", "identifier_value",
    "identity_context", "source_document_id", "ncbi_current_accession",
    "paired_accession", "biosample", "wgs_project", "secondary_pub_doi",
    "secondary_pub_pmid", "secondary_pub_title", "secondary_relationship",
    "legacy_verdict_frozen", "recommended_next_action",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1].lower()


def unique(values: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(x.strip() for x in values if x and x.strip()))


def source_document_id(source_path: str, pmid: str, doi: str) -> str:
    if source_path:
        return Path(source_path).name
    if pmid:
        return f"ABSTRACT_PMID_{pmid}"
    digest = hashlib.sha256(doi.lower().encode()).hexdigest()[:12]
    return f"ABSTRACT_DOI_{digest}"


def source_url(source_path: str, doi: str, pmid: str) -> str:
    match = PMC_RE.search(Path(source_path).name) if source_path else None
    if match:
        return f"https://pmc.ncbi.nlm.nih.gov/articles/{match.group(1).upper()}/"
    if doi:
        return f"https://doi.org/{doi}"
    if pmid:
        return f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
    return ""


def normalize_license(value: str) -> str:
    value = value.strip().replace("http://creativecommons.org", "https://creativecommons.org")
    return value.rstrip("/") if value else "UNKNOWN"


def xml_license(path: Path) -> str:
    if not path.exists() or path.suffix.lower() != ".xml":
        return ""
    try:
        root = ET.parse(path).getroot()
    except (ET.ParseError, OSError):
        return ""
    urls = []
    for license_node in (x for x in root.iter() if local_name(x.tag) == "license"):
        for node in license_node.iter():
            for key, value in node.attrib.items():
                if local_name(key) == "href" and value:
                    urls.append(value)
    cc = [x for x in urls if "creativecommons.org" in x.lower()]
    return normalize_license((cc or urls or [""])[0]) if (cc or urls) else ""


def load_license_maps(citation_dir: Path) -> tuple[dict[str, str], dict[str, str]]:
    crossref: dict[str, str] = {}
    crossref_path = citation_dir / "crossref_access_metadata.json"
    if crossref_path.exists():
        for doi, record in json.loads(crossref_path.read_text(encoding="utf-8")).items():
            urls = unique(x.get("URL", "") for x in record.get("licenses", []))
            if urls:
                crossref[doi.lower()] = normalize_license(urls[0])
    f1000: dict[str, str] = {}
    f1000_path = citation_dir / "f1000_open_text_metadata.json"
    if f1000_path.exists():
        for doi, record in json.loads(f1000_path.read_text(encoding="utf-8")).items():
            label = record.get("text_license", "")
            if label == "CC_BY_4":
                f1000[doi.lower()] = "https://creativecommons.org/licenses/by/4.0"
    return crossref, f1000


def article_license(source_path: str, doi: str, crossref: dict[str, str], f1000: dict[str, str]) -> str:
    if source_path:
        value = xml_license(Path(source_path))
        if value:
            return value
    key = doi.lower()
    return f1000.get(key) or crossref.get(key) or "UNKNOWN"


def sanitize_excerpt(value: str, anchor: str = "") -> tuple[str, str]:
    value = html.unescape(value or "")
    value = TAG_RE.sub(" ", value)
    value = SPACE_RE.sub(" ", value).strip()
    if not value:
        return "", ""
    tokens = value.split()
    treatment = "whitespace_normalized;html_removed"
    if len(tokens) > MAX_EXCERPT_TOKENS:
        index = len(tokens) // 2
        if anchor:
            anchor_lower = anchor.lower()
            for i, token in enumerate(tokens):
                if anchor_lower in token.lower():
                    index = i
                    break
        start = max(0, index - MAX_EXCERPT_TOKENS // 2)
        end = min(len(tokens), start + MAX_EXCERPT_TOKENS)
        start = max(0, end - MAX_EXCERPT_TOKENS)
        tokens = tokens[start:end]
        value = ("… " if start else "") + " ".join(tokens) + (" …" if end < len(value.split()) else "")
        treatment += ";windowed_to_25_tokens"
    return value, treatment


def fragment_id(paper_id: str, excerpt_type: str, excerpt: str, locator: str, block_hash: str) -> str:
    payload = "\t".join((paper_id, excerpt_type, excerpt, locator, block_hash))
    return "EF-" + hashlib.sha256(payload.encode()).hexdigest()[:16].upper()


def bibliographic_map(citation_dir: Path) -> dict[str, dict[str, str]]:
    path = citation_dir / "citation_table_verified.tsv"
    if not path.exists():
        return {}
    result: dict[str, dict[str, str]] = {}
    for row in read_tsv(path):
        doi = (row.get("originating_pub_doi_verified") or row.get("pub_doi") or "").lower()
        if doi and doi not in result:
            result[doi] = row
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--citation-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    packets = read_tsv(args.citation_dir / "announcement_evidence_packets.tsv")
    summaries = read_tsv(args.citation_dir / "paper_text_evidence_summary.tsv")
    summary_by_paper = {x["paper_id"]: x for x in summaries}
    bib = bibliographic_map(args.citation_dir)
    crossref, f1000 = load_license_maps(args.citation_dir)

    paper_rows = []
    paper_meta: dict[str, dict[str, str]] = {}
    for row in summaries:
        doi = row["primary_doi"]
        source_path_value = row.get("source_path", "")
        document = source_document_id(source_path_value, row["primary_pmid"], doi)
        url = source_url(source_path_value, doi, row["primary_pmid"])
        license_value = article_license(source_path_value, doi, crossref, f1000)
        biblio = bib.get(doi.lower(), {})
        public = {
            **row,
            "pub_year": biblio.get("pub_year", ""),
            "pub_journal": biblio.get("pub_journal", ""),
            "pub_authors": biblio.get("pub_authors", ""),
            "source_document_id": document,
            "source_url": url,
            "source_license": license_value,
        }
        paper_rows.append(public)
        paper_meta[row["paper_id"]] = public

    fragments: dict[str, dict[str, str]] = {}
    claim_rows = []
    for row in packets:
        paper = paper_meta[row["paper_id"]]
        production, production_treatment = sanitize_excerpt(row.get("production_excerpt", ""))
        identity, identity_treatment = sanitize_excerpt(
            row.get("identity_excerpt", ""), row.get("identifier_value", "")
        )
        production_id = ""
        identity_id = ""
        for excerpt_type, excerpt, treatment, locator, block_hash in (
            (
                "production", production, production_treatment,
                row.get("production_locator", ""), row.get("production_block_sha256", ""),
            ),
            (
                "identity", identity, identity_treatment,
                row.get("identity_locator", ""), row.get("identity_block_sha256", ""),
            ),
        ):
            if not excerpt:
                continue
            fid = fragment_id(row["paper_id"], excerpt_type, excerpt, locator, block_hash)
            fragments.setdefault(
                fid,
                {
                    "fragment_id": fid,
                    "paper_id": row["paper_id"],
                    "excerpt_type": excerpt_type,
                    "excerpt": excerpt,
                    "locator": locator,
                    "source_block_sha256": block_hash,
                    "source_document_id": paper["source_document_id"],
                    "source_url": paper["source_url"],
                    "source_license": paper["source_license"],
                    "excerpt_treatment": treatment,
                    "rights_basis": (
                        "source_open_license"
                        if "creativecommons.org" in paper["source_license"].lower()
                        else "short_quotation_for_scholarly_verification"
                    ),
                },
            )
            if excerpt_type == "production":
                production_id = fid
            else:
                identity_id = fid

        claim_rows.append(
            {
                **row,
                "support_level": row["announcement_support_level"],
                "production_fragment_id": production_id,
                "identity_fragment_id": identity_id,
                "source_document_id": paper["source_document_id"],
            }
        )

    claim_rows.sort(key=lambda x: (x["assembly_accession"], x["paper_id"]))
    fragment_rows = sorted(fragments.values(), key=lambda x: x["fragment_id"])
    review_rows = [x for x in claim_rows if x["support_level"] in {"3", "4", "5"}]

    write_tsv(args.output_dir / "papers.tsv", PAPER_FIELDS, sorted(paper_rows, key=lambda x: x["paper_id"]))
    write_tsv(args.output_dir / "evidence_fragments.tsv", FRAGMENT_FIELDS, fragment_rows)
    write_tsv(args.output_dir / "assembly_claims.tsv", CLAIM_FIELDS, claim_rows)
    write_tsv(args.output_dir / "review_queue.tsv", CLAIM_FIELDS, review_rows)

    print(f"papers={len(paper_rows)}")
    print(f"claims={len(claim_rows)}")
    print(f"unique_fragments={len(fragment_rows)}")
    print(f"review_claims={len(review_rows)}")
    if any(EMAIL_RE.search(str(value)) for row in paper_rows + fragment_rows + claim_rows for value in row.values()):
        print("refusing release: email-like content detected", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
