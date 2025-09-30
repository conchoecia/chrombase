FASTAPY_URL = https://raw.githubusercontent.com/aziele/fastapy/main/fastapy.py
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

clean:
	rm -f dependencies/fasta.py bin/datasets bin/dataformat scripts/NCBIgff2chrom.py dependencies/genbargo
