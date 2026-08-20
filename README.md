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

## Citing chrombase

If you use `chrombase` in your work, please cite the following paper:

> Schultz, D.T., Blümel, A., Destanović, D., Sarigol, F., & Simakov, O. (2026).
> Topological mixing and irreversibility in animal chromosome evolution.
> *Science Advances*, **12**(34), eadz5561.
> [https://doi.org/10.1126/sciadv.adz5561](https://doi.org/10.1126/sciadv.adz5561)

See also [`CITATION.cff`](CITATION.cff).
