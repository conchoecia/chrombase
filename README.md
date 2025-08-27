# chrombase
Tools to build a database of chromosome-scale genomes

## Installing required tools

This project relies on the NCBI `datasets` and `dataformat` command line
utilities.  The repository includes a helper script and make target that
downloads the most recent release of these binaries from GitHub and places
them in the `bin/` directory.  Simply run:

```
make
```

to fetch the tools along with other dependencies.

If your network environment lacks a proper certificate store, you can skip
SSL verification (insecure) by setting `NCBI_CLI_SKIP_SSL=1` when invoking
`make`:

```
NCBI_CLI_SKIP_SSL=1 make
```
