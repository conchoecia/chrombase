# chrombase
Tools to build a database of chromosome-scale genomes.
This is developed for and intended to be used on Linux systems.

## Prerequisites

Developed and tested on Ubuntu and variants of Red Hat Linux.
Your working environment should have the following tools installed:
- `git`
- `make`
- `snakemake`
- `python3`
- `ete4`

Before using this software, build a local NCBI Taxonomy database with
`ete4` using [this guide](https://etetoolkit.github.io/ete/tutorial/tutorial_taxonomy.html#setting-up-local-copies-of-the-ncbi-and-gtdb-taxonomy-databases).

## Installation

Clone the repository and run `make` to fetch required dependencies,
including the [fastapy](https://github.com/aziele/fastapy) parser, and NCBI [datasets](https://github.com/ncbi/datasets).

If your network environment lacks a proper certificate store, you can skip
SSL verification (insecure) by setting `NCBI_CLI_SKIP_SSL=1` when invoking
`make`:

```bash
git clone https://github.com/aziele/chrombase
cd chrombase
make
NCBI_CLI_SKIP_SSL=1 make # if you encounter SSL issues

```

After installation, Chrombase provides Python utilities and Snakemake
workflows for assembling a curated collection of chromosome-scale genome
FASTA files from the NCBI datasets service. The container includes the
dependencies needed to run the download and processing pipelines.

## Repository layout

- `src/GenDB.py` – helper functions for downloading assemblies and
  preparing FASTA files.
- `src/GenDB_build_db_annotated_chr.snakefile` – workflow for annotated
  chromosome-scale genomes.
- `src/GenDB_build_db_unannotated_chr.snakefile` – workflow for
  unannotated chromosome-scale genomes.
- `src/GenDB_scrape_genomes_NCBI.snakefile` – utility rules for fetching
  data from NCBI.

## Running a workflow

Each Snakemake file can be run directly. For example, to build a
database of annotated or unannotated chromosome-scale genomes:

```bash
# annotated assemblies
snakemake -s src/GenDB_build_db_annotated_chr.snakefile --cores <n>

# unannotated assemblies
snakemake -s src/GenDB_build_db_unannotated_chr.snakefile --cores <n>
```

The workflows rely on the NCBI `datasets` and `dataformat` command-line
tools, which must be installed and discoverable in your `PATH`.

## Running a build on SLURM

`scripts/run_build_slurm.sh` is a controller job for the two build workflows.
Submit it from the database directory, the one with `config.yaml`:

```bash
export CHROMBASE=/path/to/chrombase
export NCBI_API_KEY=...   # optional, raises the NCBI request limit
sbatch $CHROMBASE/scripts/run_build_slurm.sh annotated     # or: unannotated
```

The wrapper runs Snakemake with two flags that matter for a database that is
rebuilt over time:
- `--rerun-triggers mtime`: editing a snakefile does not rebuild genomes that
  are already finished.
- `--drop-metadata`: `.snakemake/metadata` does not grow by one file per output.

A finished genome directory holds only these files:

- annotated: `.chr.fasta.gz`, `.chrFilt.pep.gz`, `.chrFilt.chrom.gz`,
  `.chrFilt.report.txt`, `.scaffold_df.all.tsv`, `.scaffold_df.chr.tsv`,
  `.yaml.part`
- unannotated: `.chr.fasta.gz`, `_annotated_with_<LG>.pep.gz`,
  `_annotated_with_<LG>.chrom.gz`, the two `.scaffold_df` tables and
  `.yaml.part`, plus `mapped_reads/<acc>/<LG>_to_<acc>.filt.paf`

NCBI downloads are staged in `$TMPDIR` and removed when the job ends, so an
interrupted job does not leave data packages or uncompressed FASTA files in
the database. Raw miniprot alignments are deleted once they have been
filtered.

## Cleaning stale genome directories

When accession TSV files change you may want to remove genomes that are no
longer referenced.  The repository provides a small helper script that
compares the existing genome directories against the current TSVs and deletes
any extras:

```bash
python scripts/cleanup_unused_genomes.py --config config.yaml
```

By default the script processes both annotated and unannotated genome
directories.  Use `--annotated` or `--unannotated` to restrict the cleanup to
one category.  Add `--dry-run` to preview the directories that would be
removed without deleting them:

```bash
python scripts/cleanup_unused_genomes.py --dry-run --config config.yaml
```

## Recovering citations for the genomes in a database

A database built by chrombase can hold thousands of assemblies, each one the
product of someone's work. NCBI records the *submitting institution* for an
assembly but not the paper it came from, so a study using the database has no
practical way to cite the people who produced the data.

`scripts/build_citation_table.py` nominates an originating-publication candidate
for as many assemblies as the public APIs allow, and writes a reviewable table:

```bash
python scripts/build_citation_table.py \
    --genome-list genome_database/genome_list.tsv \
    --out citation_table.tsv \
    --bibtex citation_table.bib \
    --unlinked-report ncbi_missing_publication_links.tsv \
    --email you@example.org
```

It queries the NCBI Datasets and BioProject APIs and Europe PMC, scores the
candidate publications, and records how each one was found along with a
confidence level. The score measures retrieval evidence; it is not proof that a
paper produced the exact assembly build. Assemblies with no candidate keep their
submitter so they can still be credited at the group level, and
`--unlinked-report` lists candidate publication links absent from NCBI for
review.

See [`docs/citation_recovery.md`](docs/citation_recovery.md) for the resolution
routes, the scoring, and the output columns. Candidate retrieval is followed by
the exact-build validation described in
[`docs/citation_validation.md`](docs/citation_validation.md).

## Citing chrombase

If you use `chrombase` in your work, please cite the following paper:

> Schultz, D.T., Blümel, A., Destanović, D., Sarigol, F., & Simakov, O. (2026).
> Topological mixing and irreversibility in animal chromosome evolution.
> *Science Advances*, **12**(34), eadz5561.
> [https://doi.org/10.1126/sciadv.adz5561](https://doi.org/10.1126/sciadv.adz5561)

See also [`CITATION.cff`](CITATION.cff).
