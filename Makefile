FASTAPY_URL = https://raw.githubusercontent.com/aziele/fastapy/main/fastapy.py


.PHONY: all fasta ncbi-tools clean

all: fasta ncbi-tools

fasta:
	curl -L $(FASTAPY_URL) -o dependencies/fasta.py

ncbi-tools:
	python scripts/install_ncbi_cli.py

clean:
	rm -f dependencies/fasta.py dependencies/datasets dependencies/dataformat
