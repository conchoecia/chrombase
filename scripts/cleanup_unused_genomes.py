"""Remove genome directories that are no longer referenced in accession TSVs.

This utility compares the genome directories produced by the Chrombase
workflows against the current accession TSV files.  Any directories under
``odp_ncbi_genome_db/output/source_data/annotated_genomes`` or
``odp_ncbi_genome_db/output/source_data/unannotated_genomes`` that are not
listed in the corresponding TSV will be deleted.  Use ``--dry-run`` to preview
the directories that would be removed.

The script is intended to be run manually whenever the TSV files change and
stale genome data should be purged.  It is not invoked by any Snakemake
workflow.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from typing import Iterable, Set

import yaml
import types

# Import GenDB from the repository's src directory
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
SRC_DIR = os.path.join(REPO_ROOT, "src")
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

# GenDB expects the third-party ``fasta`` package.  We only need a very small
# portion of GenDB, so stub out the dependency if it is missing.
if "fasta" not in sys.modules:
    sys.modules["fasta"] = types.ModuleType("fasta")

import GenDB


def _load_accessions(config_path: str, annotated: bool) -> Set[str]:
    """Return the set of assembly accessions from the config/TSV files."""

    with open(config_path, "r", encoding="utf-8") as fh:
        config = yaml.safe_load(fh) or {}

    # Ensure tool name for directory layout
    config["tool"] = config.get("tool", "odp_ncbi_genome_db")

    # For unannotated genomes the workflow uses a temporary ignore list.
    tempignore = None
    if not annotated:
        tempignore = os.path.join(REPO_ROOT, "data", "temporary_ignore_list.txt")

    config = GenDB.opening_logic_GenDB_build_db(
        config, chr_scale=True, annotated=annotated, tempignore_path=tempignore
    )

    return set(config.get("assemAnn", []))


def _cleanup(genome_dir: str, valid: Iterable[str], dry_run: bool) -> None:
    """Remove subdirectories of ``genome_dir`` not present in ``valid``.

    Parameters
    ----------
    genome_dir:
        Path to the directory containing genome subdirectories.
    valid:
        Iterable of accession strings that should be retained.
    dry_run:
        If ``True``, only report which directories would be removed.
    """

    if not os.path.isdir(genome_dir):
        return

    valid_set = set(valid)
    for entry in sorted(os.listdir(genome_dir)):
        path = os.path.join(genome_dir, entry)
        if not os.path.isdir(path):
            continue
        if not entry.startswith(("GCA_", "GCF_")):
            continue
        if entry not in valid_set:
            if dry_run:
                print(f"Would remove {path}")
            else:
                shutil.rmtree(path)
                print(f"Removed {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Remove stale genome directories")
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to the configuration file used for database building",
    )
    parser.add_argument(
        "--annotated",
        action="store_true",
        help="Clean annotated genome directories",
    )
    parser.add_argument(
        "--unannotated",
        action="store_true",
        help="Clean unannotated genome directories",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show directories that would be removed without deleting them",
    )
    args = parser.parse_args()

    if not args.annotated and not args.unannotated:
        # Default to cleaning both
        args.annotated = True
        args.unannotated = True

    if args.annotated:
        acc = _load_accessions(args.config, annotated=True)
        ann_dir = os.path.join(
            REPO_ROOT, "odp_ncbi_genome_db", "output", "source_data", "annotated_genomes"
        )
        _cleanup(ann_dir, acc, args.dry_run)

    if args.unannotated:
        acc = _load_accessions(args.config, annotated=False)
        unann_dir = os.path.join(
            REPO_ROOT, "odp_ncbi_genome_db", "output", "source_data", "unannotated_genomes"
        )
        _cleanup(unann_dir, acc, args.dry_run)


if __name__ == "__main__":
    main()
