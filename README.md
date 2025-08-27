# chrombase
Tools to build a database of chromosome-scale genomes

## Installation

Clone the repository and run `make` to fetch required dependencies,
including the [fastapy](https://github.com/aziele/fastapy) parser that is
placed under `dependencies/fasta-parser`.

```bash
git clone https://github.com/aziele/chrombase
cd chrombase
make
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
=======
## Installing required tools

This project relies on the NCBI `datasets` and `dataformat` command line
utilities.  The repository includes a helper script and make target that
downloads the most recent release of these binaries from GitHub and places
them in the `bin/` directory.  Simply run:

```
make
```

to fetch the tools along with other dependencies.
