# Validating originating-paper candidates

`build_citation_table.py` retrieves and ranks candidate publications. This
separate validation layer asks whether each candidate actually says that its
authors produced the exact assembly build.

The 2026 snapshot covers 3,016 current assembly–paper claims across 2,600 unique
papers:

- 1,946 claims have exact text support: production language plus the exact
  assembly accession/name in either the same statement or a labelled
  deposition/build context.
- 432 combine paper production language with a project/sample identifier that
  NCBI maps to the build. These are composite support, not exact paper-text
  proof.
- 542 remain unresolved because production or exact-build identity is
  incomplete.
- 96 contain the matched identifier in downstream-use language and must not be
  accepted as originating citations on that basis.

The snapshot preserves the earlier verdict in `legacy_verdict_frozen`; a legacy
`CONFIRMED` value is not itself evidence that the paper announces the exact
build.

## Evidence hierarchy

1. Production statement and exact build identifier in the same statement.
2. Production statement plus exact accession/name in a labelled build table or
   deposition context.
3. Production statement plus a project/sample identifier linked to the build by
   NCBI metadata; manually validate.
4. Production proof or identity proof incomplete; unresolved.
5. Downstream context, absent announcement language, or no support; do not
   accept without correction.

## Public snapshot

[`validation/schultz_2026/`](../validation/schultz_2026/) contains:

- `papers.tsv`: one row per candidate paper, including source and license data.
- `evidence_fragments.tsv`: deduplicated short attributed fragments with source
  locators and SHA-256 hashes.
- `assembly_claims.tsv`: one row per assembly claim, referring to fragment IDs.
- `review_queue.tsv`: the support-level 3–5 claims.
- `audit.txt`: coverage and privacy/copyright-minimization checks.

No full text, abstract, figure, table, supplement, author email address, response
token, or local absolute path is included. Short excerpts are normalized for
whitespace, stripped of HTML, and limited to a small evidence window. Each is
linked to its paper and locator. Open-license excerpts retain their license URL;
otherwise the fragment is included only as a short quotation for scholarly
verification and commentary.

The repository GPL-3.0 license does not relicense third-party excerpts.
Copyright in those fragments remains with the respective rights holders.

## Rebuilding and auditing

Given the private/local validation outputs:

```bash
python3 scripts/citation_validation/build_public_evidence_release.py \
  --citation-dir /path/to/citation_verification \
  --output-dir validation/schultz_2026

python3 scripts/citation_validation/audit_public_evidence_release.py \
  --release-dir validation/schultz_2026 \
  --report validation/schultz_2026/audit.txt
```

Automatic extraction nominates evidence. Levels 3–5 remain review work, and a
wrong attribution remains worse than leaving the primary citation unresolved.
