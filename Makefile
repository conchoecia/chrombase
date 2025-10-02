FASTAPY_URL = https://raw.githubusercontent.com/aziele/fastapy/main/fastapy.py
NCBI_GFF2CHROM_URL = https://raw.githubusercontent.com/conchoecia/odp/refs/heads/decay-branch/scripts/NCBIgff2chrom.py
GENBARGO_URL = https://github.com/conchoecia/genbargo

.PHONY: ncbi-tools clean
.ONESHELL:

all: fasta ncbi-tools genbargo

genbargo:
	mkdir -p dependencies
	git clone $(GENBARGO_URL) dependencies/genbargo

fasta:
	mkdir -p dependencies
	curl -L $(FASTAPY_URL) -o dependencies/fasta.py

ncbi-tools:
	python scripts/install_ncbi_cli.py

ncbigff2chrom:
	mkdir -p scripts
	curl -L $(NCBI_GFF2CHROM_URL) -o scripts/NCBIgff2chrom.py
	chmod +x scripts/NCBIgff2chrom.py

clean:
	rm -f dependencies/fasta.py bin/datasets bin/dataformat scripts/NCBIgff2chrom.py dependencies/genbargo
