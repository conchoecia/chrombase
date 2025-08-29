FASTAPY_URL = https://raw.githubusercontent.com/aziele/fastapy/main/fastapy.py
GFF2CHROM_URL = https://raw.githubusercontent.com/conchoecia/odp/refs/heads/main/scripts/NCBIgff2chrom.py


.PHONY: all fasta ncbi-tools clean gff2chrom

all: fasta ncbi-tools gff2chrom

fasta:
	mkdir -p dependencies
	curl -L $(FASTAPY_URL) -o dependencies/fasta.py

ncbi-tools:
	python scripts/install_ncbi_cli.py

gff2chrom:
	curl -L $(GFF2CHROM_URL) -o scripts/NCBIgff2chrom.py

clean:
	rm -f dependencies/fasta.py bin/datasets bin/dataformat scripts/NCBIgff2chrom.py
