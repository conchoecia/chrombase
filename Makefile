FASTAPY_URL = https://raw.githubusercontent.com/aziele/fastapy/main/fastapy.py


.PHONY: all fasta ncbi-tools clean gff2chrom

all: fasta ncbi-tools gff2chrom

fasta:
	mkdir -p dependencies
	curl -L $(FASTAPY_URL) -o dependencies/fasta.py

ncbi-tools:
	python scripts/install_ncbi_cli.py

clean:
	rm -f dependencies/fasta.py bin/datasets bin/dataformat scripts/NCBIgff2chrom.py
