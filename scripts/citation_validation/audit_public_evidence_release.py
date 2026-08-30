#!/usr/bin/env python3
"""Audit the public, copyright-minimized citation evidence snapshot."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path


EMAIL_RE = re.compile(r"[A-Z0-9._%+\-]+@[A-Z0-9.\-]+\.[A-Z]{2,}", re.I)
ABSOLUTE_PATH_RE = re.compile(r"(?:^|[\t ])/(?:home|Users|media|tmp)/")
TAG_RE = re.compile(r"<[^>]+>")
EXPECTED = {"papers": 2600, "claims": 3016, "review": 1070}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release-dir", required=True, type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()

    errors = []
    papers = read_tsv(args.release_dir / "papers.tsv")
    fragments = read_tsv(args.release_dir / "evidence_fragments.tsv")
    claims = read_tsv(args.release_dir / "assembly_claims.tsv")
    review = read_tsv(args.release_dir / "review_queue.tsv")

    if len(papers) != EXPECTED["papers"]:
        errors.append(f"papers expected {EXPECTED['papers']}, found {len(papers)}")
    if len(claims) != EXPECTED["claims"]:
        errors.append(f"claims expected {EXPECTED['claims']}, found {len(claims)}")
    if len(review) != EXPECTED["review"]:
        errors.append(f"review expected {EXPECTED['review']}, found {len(review)}")

    paper_ids = [x["paper_id"] for x in papers]
    fragment_ids = [x["fragment_id"] for x in fragments]
    if len(paper_ids) != len(set(paper_ids)):
        errors.append("duplicate paper_id")
    if len(fragment_ids) != len(set(fragment_ids)):
        errors.append("duplicate fragment_id")
    paper_set = set(paper_ids)
    fragment_set = set(fragment_ids)

    claim_keys = [(x["assembly_accession"], x["paper_id"]) for x in claims]
    if len(claim_keys) != len(set(claim_keys)):
        errors.append("duplicate assembly-paper claim")
    for row in claims:
        if row["paper_id"] not in paper_set:
            errors.append(f"unknown paper {row['paper_id']}")
        for column in ("production_fragment_id", "identity_fragment_id"):
            if row[column] and row[column] not in fragment_set:
                errors.append(f"unknown fragment {row[column]}")

    for row in fragments:
        if row["paper_id"] not in paper_set:
            errors.append(f"fragment has unknown paper {row['paper_id']}")
        if len(row["excerpt"].split()) > 27:  # 25 words plus up to two ellipsis tokens
            errors.append(f"overlong fragment {row['fragment_id']}")
        if TAG_RE.search(row["excerpt"]):
            errors.append(f"HTML remains in {row['fragment_id']}")
        if not row["source_url"]:
            errors.append(f"missing source URL for {row['fragment_id']}")

    all_paths = [
        args.release_dir / "papers.tsv",
        args.release_dir / "evidence_fragments.tsv",
        args.release_dir / "assembly_claims.tsv",
        args.release_dir / "review_queue.tsv",
    ]
    for path in all_paths:
        text = path.read_text(encoding="utf-8")
        if EMAIL_RE.search(text):
            errors.append(f"email-like content in {path.name}")
        if ABSOLUTE_PATH_RE.search(text):
            errors.append(f"absolute local path in {path.name}")

    support = Counter(x["support_level"] for x in claims)
    statuses = Counter(x["text_evidence_status"] for x in claims)
    if support != Counter({"1": 19, "2": 1927, "3": 432, "4": 402, "5": 236}):
        errors.append(f"unexpected support counts: {dict(support)}")
    expected_review_keys = {
        (x["assembly_accession"], x["paper_id"])
        for x in claims
        if x["support_level"] in {"3", "4", "5"}
    }
    if expected_review_keys != {(x["assembly_accession"], x["paper_id"]) for x in review}:
        errors.append("review queue does not equal support levels 3-5")

    report_lines = [
        f"papers={len(papers)}",
        f"claims={len(claims)}",
        f"unique_fragments={len(fragments)}",
        f"review_claims={len(review)}",
        "support_counts=" + ";".join(f"{k}:{v}" for k, v in sorted(support.items())),
        "status_counts=" + ";".join(f"{k}:{v}" for k, v in sorted(statuses.items())),
        f"audit={'FAIL' if errors else 'PASS'}",
    ]
    report = "\n".join(report_lines) + "\n"
    if args.report:
        args.report.write_text(report, encoding="utf-8")
    print(report, end="")
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
