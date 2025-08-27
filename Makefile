FASTAPY_URL = https://raw.githubusercontent.com/aziele/fastapy/main/fastapy.py

.PHONY: fasta clean

fasta:
	curl -L $(FASTAPY_URL) -o src/fasta.py

clean:
	rm -f src/fasta.py
