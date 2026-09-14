#!/usr/bin/env python

"""
This python program contains functions that are used to download genomes from NCBI.
  - These tools are common to the scripts:
    - GenDB_build_db_unannotated_chr.snakefile
    - GenDB_build_db_unannotated_nonchr.snakefile
    - GenDB_build_db_annotated_chr.snakefile
    - GenDB_build_db_annotated_nonchr.snakefile
"""

import gzip
import io
from io import StringIO
import math
import os
import pandas as pd
from random import randint
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import warnings
import zipfile

# use the fasta package included with this repo
snakefile_path = os.path.dirname(os.path.realpath(__file__))
dependencies_path = os.path.join(snakefile_path, "../dependencies/fasta-parser")
sys.path.insert(1, dependencies_path)
import fasta

def create_directories_recursive_notouch(path):
    """
    Unlike os.makedirs, this function will not touch a directory if it already exists.
    This is useful for snakemake because it will not re-run a rule if the output already exists.
    """
    parts = os.path.normpath(path).split(os.path.sep)
    # An absolute path splits into an empty first part, and os.mkdir("") fails.
    current_path = os.path.sep if os.path.isabs(path) else ""

    for part in parts:
        if not part:
            continue
        current_path = os.path.join(current_path, part)
        if not os.path.exists(current_path):
            os.mkdir(current_path)

def _double_per_attempt(first, attempt, cap):
    """
    A resource request that doubles on every Snakemake retry and stops at `cap`.
    Unlike the lookup tables used before, this is defined for any attempt number. A table with
      fewer entries than the rule has attempts raised KeyError, which crashed the controller.
    """
    return int(min(cap, first * 2 ** (attempt - 1)))

def dlChrs_get_mem_mb(wildcards, attempt):
    """
    The genome is streamed out of the NCBI zip line by line, so the RAM needed no longer grows
      with chromosome size. The scaffold table is the largest thing held in memory.
    """
    return _double_per_attempt(2000, attempt, 16000)

def dlChrs_get_runtime(basepairs, attempt) -> int:
    """
    Minutes to download an assembly and write its chromosome-scale scaffolds to .chr.fasta.gz.
    30 minutes (this includes up to 4 minutes of randomized sleeps) plus 6 minutes per Gbp of
      assembly, doubled on each retry. The fixed 20 minute limit used before killed the largest
      downloads part way through. Download speed from NCBI varies more than genome size does:
      in a 34 genome test build a 2.3 Gbp genome took 22 minutes and a 19.8 Gbp genome 37 minutes.

    Returns an int, because SLURM only takes ints for runtime.
    """
    return int(math.ceil((30 + 6.0 * (basepairs / 1e9)) * 2 ** (attempt - 1)))

def miniprot_get_mem_mb(basepairs, attempt):
    """
    The amount of RAM needed for miniprot is dominated by indexing the genome, which grows with
      the assembly size. The first request is 14 GB per Gbp of assembly (at least 16 GB), and it
      doubles on each retry up to 1.5 TB.

    UPDATES:
      - 202609 - In the September 2025 build every genome started at 16 GB. That ran out of memory
                 for 400 genomes of 1.5 to 5 Gbp, and the largest requests that still ran out were
                 about 11 GB per Gbp. No genome under 1.4 Gbp ran out of memory at 16 GB.
      - 202609 - In a test build of 19 genomes of 0.1 to 22 Gbp, genomes of 1.7 Gbp and larger used
                 6.5 to 10 GB per Gbp.
    """
    first = max(16000, 14000 * (basepairs / 1e9))
    return _double_per_attempt(first, attempt, 1536000)

def prep_chrom_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome, the number of proteins, and the gff size.
    4 GB, doubling on each retry up to 64 GB.
    """
    return _double_per_attempt(4000, attempt, 64000)

def prep_chrom_get_time(wildcards, attempt):
    """
    The amount of minutes needed varies depending on the input size.
    16 minutes, doubling on each retry up to 256 minutes.
    """
    return _double_per_attempt(16, attempt, 256)


def print_gendb_config_summary(config, chr_scale=True, annotated=True, sample_n=5):
    tsv_keys = [
        "annotated_genome_chr_tsv",
        "annotated_genome_nonchr_tsv",
        "unannotated_genome_chr_tsv",
        "unannotated_genome_nonchr_tsv",
    ]

    print("\n[GenDB] Accessions TSV paths (from config):", flush=True)
    for k in tsv_keys:
        v = config.get(k)
        if v:
            abs_v = os.path.abspath(v)
            exists = os.path.exists(v)
            print(f"  {k}: {abs_v}  (exists={exists})", flush=True)
        else:
            print(f"  {k}: MISSING", flush=True)

    # Determine chosen_file using same logic as opening_logic_GenDB_build_db
    if annotated:
        chosen_key = "annotated_genome_chr_tsv" if chr_scale else "annotated_genome_nonchr_tsv"
    else:
        chosen_key = "unannotated_genome_chr_tsv" if chr_scale else "unannotated_genome_nonchr_tsv"

    chosen_file = config.get(chosen_key)
    print("", flush=True)
    if chosen_file:
        print(f"[GenDB] chosen_file key: {chosen_key}", flush=True)
        print(f"[GenDB] chosen_file path: {os.path.abspath(chosen_file)}  (exists={os.path.exists(chosen_file)})", flush=True)
    else:
        print(f"[GenDB][WARN] No chosen_file found in config for key: {chosen_key}", flush=True)

    # Summary of assemblies and scaffold-length mapping from config
    assem_list = config.get("assemAnn")
    scaflen_map = config.get("assemAnn_to_scaflen")

    if assem_list is None:
        print("[GenDB][WARN] config['assemAnn'] is missing or None", flush=True)
    else:
        print(f"[GenDB] Total assemblies (assemAnn): {len(assem_list)}", flush=True)
        if len(assem_list) > 0:
            print(f"[GenDB] First {min(sample_n, len(assem_list))} assembly accessions: {assem_list[:sample_n]}", flush=True)

    if scaflen_map is None:
        print("[GenDB][WARN] config['assemAnn_to_scaflen'] is missing or None", flush=True)
    else:
        # Show a few sample mappings
        items = list(scaflen_map.items())
        print(f"[GenDB] Total assembly->scaflen entries: {len(items)}", flush=True)
        if items:
            print(f"[GenDB] Sample assembly->scaflen mappings (first {min(sample_n, len(items))}):", flush=True)
            for asm, slen in items[:sample_n]:
                print(f"    {asm} -> {slen}", flush=True)

    # Optional integrity checks
    if assem_list is not None and scaflen_map is not None:
        missing_in_map = [a for a in assem_list if a not in scaflen_map]
        if missing_in_map:
            print(f"[GenDB][WARN] {len(missing_in_map)} assembly accessions present in 'assemAnn' but missing from 'assemAnn_to_scaflen' (showing up to 5): {missing_in_map[:5]}", flush=True)
    print("", flush=True)

def remove_entries_from_df(df, colname, remove_list):
    """
    This function takes in a pandas dataframe, a column name, and a list of entries to remove from that column.
    It returns a new dataframe with the entries removed.

    Input:
        - df:          a pandas dataframe
        - colname:     the name of the column to check
        - remove_list: a list of entries to remove from the column
    Output:
        - df:          a pandas dataframe with the entries removed
    """
    start_len = len(df)
    df = df[~df[colname].isin(remove_list)]
    end_len = len(df)
    print("Removed {} entries from the dataframe based on the remove_list.".format(start_len - end_len), file = sys.stderr)
    return df


def opening_logic_GenDB_build_db(config, chr_scale = None, annotated = None, tempignore_path = None):
    """
    This is common logic for all of the GenDB_build_db_*.snakefile scripts.
    Basically just pertains to parsing which files will go into this category.

    Takes the config object as input, and returns it modified.
    The user must provide boolean values (True, False) for chr_scale and annotated.

    Opens a directory containing the accesstion tsv files, finds the latest accessions.
    """
    # Do some logic to see if the user has procided enough information for us to analyse the genomes
    if ("directory" not in config) and ("accession_tsvs" not in config):
        raise IOerror("You must provide either a directory of the annotated and unannotated genome lists, or a list of the paths to those tsv files. Read the config file.")

    if "directory" in config:
        # ensure that the user also hasn't specified the accession tsvs
        if "accession_tsvs" in config:
            raise IOerror("You cannot provide both a directory of tsv files ('directory') and specify the exact path of the TSVs ('accession_tsvs'). Read the config file.")
        # ensure that the directory exists
        if not os.path.isdir(config["directory"]):
            raise IOerror("The directory of TSV files you provided does not exist. {}".format(config["directory"]))
        # get the paths to the tsv files
        latest_accesions = return_latest_accession_tsvs(config["directory"])
        # This is the output of the above function:
        #return_dict = {"annotated_chr":      os.path.join(directory_path, most_recent_annotated_chr_file),
        #               "annotated_nonchr":   os.path.join(directory_path, most_recent_annotated_nonchr_file),
        #               "unannotated_chr":    os.path.join(directory_path, most_recent_unannotated_chr_file),
        #               "unannotated_nonchr": os.path.join(directory_path, most_recent_unannotated_nonchr_file)
        #               }
        config["annotated_genome_chr_tsv"]       = latest_accesions["annotated_chr"]
        config["annotated_genome_nonchr_tsv"]    = latest_accesions["annotated_nonchr"]
        config["unannotated_genome_chr_tsv"]     = latest_accesions["unannotated_chr"]
        config["unannotated_genome_nonchr_tsv"]  = latest_accesions["unannotated_nonchr"]

        # Which file we use depends on whether we are doing annotated or unannotated genomes,
        #  and whether we are doing chromosome or non-chromosome scale genomes.
        chosen_file = None
        if annotated:
            if chr_scale:
                chosen_file = config["annotated_genome_chr_tsv"]
            else:
                chosen_file = config["annotated_genome_nonchr_tsv"]
        else:
            if chr_scale:
                chosen_file = config["unannotated_genome_chr_tsv"]
            else:
                chosen_file = config["unannotated_genome_nonchr_tsv"]
        #if chosen_file is None:
        #    raise IOError("You provided a directory for us to find the lists of genomes {}, but")
        # Use the chosen file, and add the entries to the config file so we can download them or not
        df = pd.read_csv(chosen_file, sep="\t", index_col=False)
        # strip leading and trailing whitespace from the column names because pandas can screw up sometimes
        df.columns = df.columns.str.strip()

        config["assemAnn"]            = df["Assembly Accession"].tolist()
        # make a dict where assemAnn is the key and the value is the "Assembly Stats Total Ungapped Length" column
        config["assemAnn_to_scaflen"] = dict(zip(df["Assembly Accession"].tolist(), df["Assembly Stats Total Ungapped Length"].tolist()))

        if tempignore_path is not None:
            # ensure that the file exists
            if not os.path.isfile(tempignore_path):
                raise IOError("The temporary ignore list file you provided does not exist. {}".format(tempignore_path))
            # read in the file, and make a list of the entries to remove
            remove_list = []
            with open(tempignore_path, "r") as f:
                for line in f:
                    line = line.strip()
                    if line and line[0] != "#":
                        remove_list.append(line)
            print("Removing {} entries from the dataframe based on the temporary ignore list.".format(len(remove_list)), file = sys.stderr)
            for entry in remove_list:
                config["assemAnn"].remove(entry)

    elif "accession_tsvs" in config:
        # ensure that the user also hasn't specified the directory
        if "directory" in config:
            raise IOerror("You cannot provide both a directory of tsv files ('directory') and specify the exact path of the TSVs ('accession_tsvs'). Read the config file.")
        # we haven't implemented this yet. I haven't found a use case where I would want to specifially pick a file path rather than just get the most recent one.
        raise NotImplementedError("We haven't implemented this yet. I haven't found a use case where I would want to specifially pick a file path rather than just get the most recent one.")

    return config

def download_assembly_scaffold_df(assembly_accession, datasetsEx, chrscale = True) -> pd.DataFrame:
    """
    Handles downloading the assembly scaffold df from NCBI.
    Uses the NCBI datasets executable tool to do so.

    Input:
        - assembly_accession: the assembly accession number of the genome to download
        - datasetsEx: path to the NCBI datasets executable
    Output:
        - df: a pandas dataframe containing the assembly scaffold df

    df output example:
                assembly_accession     assembly_unit assigned_molecule_location_type chr_name    gc_count  gc_percent  genbank_accession    length                role
        0      GCA_026151205.1  Primary Assembly                      Chromosome        1  10152078.0        38.5         CM047403.1  26206493  assembled-molecule
        1      GCA_026151205.1  Primary Assembly                      Chromosome        2   8914563.0        38.5         CM047404.1  23028209  assembled-molecule
        2      GCA_026151205.1  Primary Assembly                      Chromosome        3   8250754.0        39.0         CM047405.1  21284423  assembled-molecule
        3      GCA_026151205.1  Primary Assembly                      Chromosome        4   6284100.0        38.5         CM047406.1  16326010  assembled-molecule
        4      GCA_026151205.1  Primary Assembly                      Chromosome        5   9163041.0        39.0         CM047407.1  23637599  assembled-molecule
        ..                 ...               ...                             ...      ...         ...         ...                ...       ...                 ...
        342    GCA_026151205.1  Primary Assembly                      Chromosome       Un         NaN         NaN  JAIOUN010000343.1      3450   unplaced-scaffold
        343    GCA_026151205.1  Primary Assembly                      Chromosome       Un         NaN         NaN  JAIOUN010000344.1      3041   unplaced-scaffold
        344    GCA_026151205.1  Primary Assembly                      Chromosome       Un         NaN         NaN  JAIOUN010000345.1      2526   unplaced-scaffold
        345    GCA_026151205.1  Primary Assembly                      Chromosome       Un         NaN         NaN  JAIOUN010000346.1      1772   unplaced-scaffold
        346    GCA_026151205.1       non-nuclear                   Mitochondrion       MT      2109.0        20.0         CM047416.1     10476  assembled-molecule
    """
    # wait a random amount of time up to 120 seconds to space out requests to the NCBI servers.
    _polite_sleep("download the assembly df for {}".format(assembly_accession))

    # first use the datasets executable to get a list of scaffolds to download
    # launch it as a process from bash and collect the output to sys.stdout
    # ./datasets summary genome accession GCF_000001405.40 --report sequence --as-json-lines
    cmd = [datasetsEx, "summary", "genome", "accession", assembly_accession, "--report", "sequence", "--as-json-lines"] + _datasets_api_key_args()
    p = ""
    if chrscale:
        # if we are downloading the chromosome-scale scaffolds, then we can pre-filter the non-chromosome-scale scaffolds
        # add this to the command: | grep -v 'unplaced-scaffold'
        cmd = " ".join(shlex.quote(x) for x in cmd) + " | grep -v 'unplaced-scaffold' | grep -v 'chr_name\":\"Un'"
        print("Trying the following command: {}".format(_redact_api_key(cmd)), file = sys.stderr)
        p = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    else:
        print("Trying the following command: {}".format(_redact_api_key(" ".join(cmd))), file = sys.stderr)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    out, err = p.communicate()
    err = _redact_api_key(err)
    sys.stderr.write(err)
    if p.returncode != 0:
        raise ValueError("datasets summary for {} returned a non-zero exit code. Something went wrong. {}".format(assembly_accession, err))
    # the output is json lines. Convert it to a pandas dataframe
    json_buffer = StringIO(out)
    df = pd.read_json(json_buffer, lines=True)
    # return an error if the dataframe is empty
    if df.empty:
        raise ValueError("The dataframe is empty. This is probably because the assembly accession number you provided is not valid. {}".format(assembly_accession))
    # return a dataframe if there are any values in the "assembly_accession" column that are not the same as the one we provided
    if not df["assembly_accession"].unique()[0] == assembly_accession:
        raise ValueError("The assembly accession number you provided is not valid. {}".format(assembly_accession))
    # return the dataframe
    return df

def filter_scaffold_df_keep_chrs(df) -> pd.DataFrame:
    """
    Takes in an NCBI dataset of the scaffold entries and filters it down to a df that only contains the chromosome-scale scaffolds.
    The input dataframe looks like the one detailed in the docstring of download_assembly_scaffold_df(assembly_accession, datasetsEx)
    The strategy to get the chromosome-scale scaffolds is mostly to remove the things that we know aren't chromosome-scale.

    Input:
        - scaffold_df: a pandas dataframe containing the assembly scaffold df
    Output:
        - df: a pandas dataframe containing only the chromosome-scale scaffolds

    Notes:
        - DO NOT filter on "Chromosome" in the "assigned_molecule_location_type" column. This will remove linkage groups, which are also chromosome-scale scaffolds.
          Instead remove "Mitochondrion" entries from the "assigned_molecule_location_type" column.
    """
    startdf = df.copy()
    # get the most common assembly accession number in the df
    assembly_accession = df["assembly_accession"].mode()[0]
    start_length = len(df)

    # Remove mitochondria
    if "assigned_molecule_location_type" in df.columns:
        df = df[df["assigned_molecule_location_type"] != "Mitochondrion"]
    if "role" in df.columns:
        # remove things called "unplaced-scaffold" - these are non-chromosome-scale scaffolds for which the chromosome is uncertain
        df = df[df["role"] != "unplaced-scaffold"]
        # remove things called "alt-scaffold" - these are not the main scaffolds for the annotation
        df = df[df["role"] != "alt-scaffold"]
        # remove things that are called "unlocalized-scaffold", the chromosome to which these belond is known, but the location on the chromosome is not known
        df = df[df["role"] != "unlocalized-scaffold"]
    if "chr_name" in df.columns:
        # remove chr_name Un
        df = df[df["chr_name"] != "Un"]
    # raise an error if the dataframe is empty - we expect there to be at least one chromosome-scale scaffold
    if df.empty:
        print("This was the startdf: \n{}".format(startdf), file = sys.stderr)
        raise ValueError("The dataframe is empty after filtering for chromosome-scale scaffolds. The start length of the df was {}. The assembly_accession is {}".format(start_length, assembly_accession))
    return df

def _polite_sleep(label):
    """
    Wait a random amount of time before contacting NCBI, to space out the requests of many concurrent jobs.
    The maximum is 120 seconds by default. Set the environment variable CHROMBASE_MAX_SLEEP to lower it,
      or to 0 to skip the wait.
    """
    max_sleep = int(os.environ.get("CHROMBASE_MAX_SLEEP", "120"))
    if max_sleep <= 0:
        return
    sleeptime = randint(1, max_sleep)
    print("About to {}.\nSleeping for {} seconds to avoid overloading the NCBI servers.".format(label, sleeptime), file = sys.stderr)
    time.sleep(sleeptime)

def _datasets_api_key_args():
    """
    The --api-key arguments for the NCBI datasets tool, if NCBI_API_KEY is set in the environment.
    """
    key = os.environ.get("NCBI_API_KEY", "")
    return ["--api-key", key] if key else []

def _redact_api_key(text):
    """
    Removes the NCBI API key from a string before it is printed to a log.
    """
    key = os.environ.get("NCBI_API_KEY", "")
    return text.replace(key, "<NCBI_API_KEY>") if key else text

class _GzipTextWriter:
    """
    A text handle that writes a gzip file, through pigz when it is installed and through the
      gzip module otherwise. Both use compression level 6, like `gzip` on the command line.
    If the block inside `with` raises, the partial file is abandoned instead of finalized.
    """
    def __init__(self, path, threads = 1):
        self.path = path
        self.fileobj = open(path, "wb")
        self.proc = None
        if shutil.which("pigz"):
            self.proc = subprocess.Popen(["pigz", "-p", str(max(1, threads)), "-6", "-c"],
                                         stdin = subprocess.PIPE, stdout = self.fileobj, bufsize = 1 << 20)
            self.handle = io.TextIOWrapper(self.proc.stdin, encoding = "utf-8", newline = "\n")
        else:
            gz = gzip.GzipFile(fileobj = self.fileobj, mode = "wb", compresslevel = 6)
            self.handle = io.TextIOWrapper(gz, encoding = "utf-8", newline = "\n")

    def write(self, text):
        self.handle.write(text)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        if exc_type is None:
            self.handle.close()
            if self.proc is not None and self.proc.wait() != 0:
                self.fileobj.close()
                raise IOError("pigz exited with code {} while writing {}".format(self.proc.returncode, self.path))
        else:
            try:
                self.handle.close()
            except Exception:
                pass
            if self.proc is not None:
                self.proc.kill()
                self.proc.wait()
        self.fileobj.close()
        return False

def _write_selected_scaffolds(handle, out, lengths_by_col, assembly_accession, wrap = 80):
    """
    Streams a fasta file handle and writes only the chromosome-scale scaffolds to `out`.
    The output of each record is the same as fasta.Record.format(wrap=80): the header as
      '>id description' and the sequence wrapped to 80 columns.

    Inputs:
      - handle:         a text handle of the NCBI genomic fasta file
      - out:            a text handle to write to
      - lengths_by_col: {"genbank_accession": {scaffold: length}, "refseq_accession": {...}}. The first
                        header found in one of these decides which naming the fasta file uses.
      - assembly_accession: only used in messages
    Returns:
      - the dictionary of {scaffold: expected length} that was used, and the set of scaffolds written
    """
    wanted  = None    # {scaffold: expected length} of the column that matches the fasta headers
    seen    = set()
    name    = None    # the scaffold currently being written, None while skipping a record
    carry   = ""      # sequence not yet written because it is shorter than a full line
    seqlen  = 0

    def finish_record(name, carry, seqlen):
        if carry:
            out.write(carry + "\n")
        # The lengths in the sequence report do not always match, see https://github.com/ncbi/datasets/issues/301
        if seqlen != wanted[name]:
            print("WARNING: scaffold {} of {} is {} bp in the fasta file but {} bp in the sequence report.".format(
                  name, assembly_accession, seqlen, wanted[name]), file = sys.stderr)

    for line in handle:
        if line.startswith(">"):
            if name is not None:
                finish_record(name, carry, seqlen)
            thisid = line.split()[0][1:]
            name = None
            if wanted is None:
                for thiscol, lengths in lengths_by_col.items():
                    if thisid in lengths:
                        wanted = lengths
                        break
            if (wanted is not None) and (thisid in wanted):
                if thisid in seen:
                    raise ValueError("The scaffold {} appears twice in the fasta file of {}. Something went wrong.".format(thisid, assembly_accession))
                seen.add(thisid)
                name   = thisid
                desc   = line[len(thisid) + 1:].strip()
                out.write(">{} {}\n".format(thisid, desc) if desc else ">{}\n".format(thisid))
                carry  = ""
                seqlen = 0
        elif name is not None:
            seq = line.rstrip()
            seqlen += len(seq)
            carry  += seq
            if len(carry) >= wrap:
                full = len(carry) - (len(carry) % wrap)
                for i in range(0, full, wrap):
                    out.write(carry[i:i + wrap] + "\n")
                carry = carry[full:]
    if name is not None:
        finish_record(name, carry, seqlen)
    return wanted, seen

def download_chr_scale_genome(chr_df, datasetsEx, out_fasta_gz, threads = 1):
    """
    NOTE: until https://github.com/ncbi/datasets/issues/298 IS FIXED we download the entire genome and filter it ourselves.
    This function takes in a dataframe of known chromosome-scale scaffolds for a specific assembly accession number,
      downloads the genome from NCBI, and writes only those scaffolds to a gzipped fasta file.
    It replaces download_chr_scale_genome_from_df_WORKING, which unpacked the zip and wrote an uncompressed
      .chr.fasta next to the database files. When a job was killed, those were left behind.

    How it works:
      - The zip is downloaded into a temporary directory under $TMPDIR (job-local on SLURM). The directory
        is removed when the function exits, whether it succeeds or fails.
      - The genomic fasta is read straight out of the zip one line at a time. Nothing is unpacked, and the
        RAM needed does not grow with chromosome length.
      - The scaffolds are compressed as they are written, into the temporary directory. The finished file
        is copied to out_fasta_gz only after all of the checks pass.

    Special checks that are made along the way to make sure we get exactly what is needed:
      - Every chromosome-scale scaffold in chr_df is written exactly once, and nothing else is written.
      - Scaffolds whose length differs from the "length" column are reported on stderr.

    Inputs:
        - chr_df:       a pandas dataframe containing only the chromosome-scale scaffolds
        - datasetsEx:   path to the NCBI datasets executable
        - out_fasta_gz: path of the gzipped fasta file to write
        - threads:      number of compression threads
    Outputs:
        - Outputs 0 if the function completes successfully.
    """
    # get the most common value of the assembly_accession column
    assembly_accession = chr_df["assembly_accession"].mode()[0]
    _polite_sleep("download the genome assembly for {}".format(assembly_accession))

    lengths_by_col = {}
    for thiscol in ["genbank_accession", "refseq_accession"]:
        if thiscol in chr_df.columns:
            lengths_by_col[thiscol] = dict(zip(chr_df[thiscol].tolist(), chr_df["length"].tolist()))

    workdir = tempfile.mkdtemp(prefix = "chrombase_{}_".format(assembly_accession))
    try:
        # DOWNLOAD THE ASSEMBLY
        outzip = os.path.join(workdir, "{}.zip".format(assembly_accession))
        cmd = [datasetsEx, "download", "genome", "accession", assembly_accession,
               "--include", "genome", "--no-progressbar", "--filename", outzip] + _datasets_api_key_args()
        print("the command is {}".format(_redact_api_key(" ".join(cmd))), file = sys.stderr)
        p = subprocess.run(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
        sys.stdout.write(p.stdout)
        sys.stderr.write(_redact_api_key(p.stderr))
        if p.returncode != 0:
            raise ValueError("The datasets executable for {} returned a non-zero exit code. Something went wrong. {}".format(assembly_accession, _redact_api_key(p.stderr)))
        if not os.path.exists(outzip):
            raise ValueError("The output zip file does not exist. Something went wrong with the download. {}".format(outzip))

        # WRITE THE CHROMOSOME-SCALE SCAFFOLDS. The main assembly file is called GCF_905220415.1_{something}_genomic.fna
        tmp_gz = os.path.join(workdir, os.path.basename(out_fasta_gz))
        with zipfile.ZipFile(outzip) as zf:
            prefix  = "ncbi_dataset/data/{}/".format(assembly_accession)
            members = [x for x in zf.namelist() if x.startswith(prefix) and x.endswith("_genomic.fna")]
            if len(members) != 1:
                raise IOError("For Assembly Accession {}, expected one {}*_genomic.fna in the download but found {}.".format(assembly_accession, prefix, members))
            with zf.open(members[0]) as raw, _GzipTextWriter(tmp_gz, threads) as out:
                handle = io.TextIOWrapper(io.BufferedReader(raw, 1 << 20), encoding = "utf-8")
                wanted, seen = _write_selected_scaffolds(handle, out, lengths_by_col, assembly_accession)

        if wanted is None:
            raise ValueError("We didn't find any scaffolds in the fasta file that were in the dataframe. Something went wrong. Assembly accession: {}".format(assembly_accession))
        missing = set(wanted) - seen
        if missing:
            raise ValueError("These chromosome-scale scaffolds of {} were not in the downloaded fasta file: {}".format(assembly_accession, sorted(missing)))
        if len(seen) != len(chr_df):
            raise ValueError("The number of scaffolds in the output fasta file is not the same as the number of scaffolds in the input dataframe. Something went wrong. Assembly acession: {}".format(assembly_accession))

        # MOVE THE FINISHED FILE INTO THE DATABASE
        outdir = os.path.dirname(out_fasta_gz)
        if outdir:
            os.makedirs(outdir, exist_ok = True)
        shutil.copyfile(tmp_gz, out_fasta_gz)
    finally:
        shutil.rmtree(workdir, ignore_errors = True)
    return 0

def download_chr_scale_genome_from_df(chr_df, datasetsEx, output_dir):
    """
    NOTE: THIS WILL ONLY WORK WHEN https://github.com/ncbi/datasets/issues/298 IS FIXED
    This function takes in a dataframe of known chromosome-scale scaffolds for a specific assembly accession number,
      and downloads the genome from NCBI.

    Special checks that are made along the way to make sure we get exactly what is needed:
      - We check that we only add each needed scaffold exactly once. This avoids a problem of using a bash-based solution
        where we cannot finely control which fasta files are concatenated together.
      - We check that each fasta entry is exactly the length that is expected. This will help find cases where the
        download or decompression did not work for some reason.

    Other features:
      - This function will delete all of the unneeded intermediate files.

    Inputs:
        - chr_df:      a pandas dataframe containing only the chromosome-scale scaffolds
        - datasetsEx:  path to the NCBI datasets executable
        - output_dir:  The directory to which additional files will be saved.
                         Everything will be prefixed with the assembly_accession value.
    Outputs:
        - Outputs 0 if the function completes successfully.
        - Saves a fasta file to the output_dir that contains all of the chromosome-scale scaffolds.
          The file will be saved at {output_dir}/{assembly_accession}.chr.fasta
    """
    # get the most common value of the assembly_accession column
    assembly_accession = chr_df["assembly_accession"].mode()[0]
    # First, we need to check that the output directory exists. If it doesn't, then we need to create it.
    create_directories_recursive_notouch(output_dir)
    # cast the chr_name column to strings so we can work with them as strings later
    chr_df["chr_name"] = chr_df["chr_name"].astype(str)

    # DOWNLOAD THE SCAFFOLDS
    # Now we need to figure out which scaffolds we need to download. To do this, we can use the chromosome names and prompt the datasets tool to download them
    chrom_download_string = ",".join(chr_df["chr_name"].tolist())
    print("the chrom download string is {}".format(chrom_download_string), file = sys.stderr)
    # use the assembly accession number to name the output zip file
    outzip = os.path.join(output_dir, "{}.zip".format(assembly_accession))
    # set up the download process using the datasets executable
    cmd = [datasetsEx, "download", "genome", "accession", assembly_accession, "--chromosomes", chrom_download_string, "--filename", outzip]
    print("the command is {}".format(" ".join(cmd)), file = sys.stderr)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
    # wait for the process to finish
    out, err = p.communicate()
    sys.stdout.write(out)
    sys.stderr.write(err)
    if p.returncode != 0:
        raise ValueError("The datasets executable for {} returned a non-zero exit code. Something went wrong. {}".format(assembly_accession, err))
    # check that the output zip file exists
    if not os.path.exists(outzip):
        raise ValueError("The output zip file does not exist. Something went wrong with the download. {}".format(outzip))

    # UNZIP THE SCAFFOLDS
    # unzip the file, forcing it to overwrite any existing files, excluding the file called README.md
    unzip_dir = os.path.join(output_dir, assembly_accession)
    cmd = ["unzip", "-o", outzip, "-d", unzip_dir]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
    # wait for the process to finish
    out, err = p.communicate()
    sys.stdout.write(out)
    sys.stderr.write(err)
    if p.returncode != 0:
        raise ValueError("Unzip for {} returned a non-zero exit code. Something went wrong. {}".format(assembly_accession, err))
    # check that the output directory exists
    if not os.path.isdir(unzip_dir):
        raise ValueError("The output directory does not exist. Something went wrong with the unzip. {}".format(output_dir))

    # CHECK THAT ALL THE SCAFFOLD FILES ARE THERE
    for thischrom in chr_df["chr_name"].tolist():
        thisfile = os.path.join(unzip_dir, "ncbi_dataset/data/{}/chr{}.fna".format(assembly_accession,thischrom))
        if not os.path.exists(thisfile):
            # The file does not exist. Something went wrong with the download of one of the chromosomes
            raise IOError("For Assembly Accession {}, the file {} does not exist. Something went wrong with the download of this chromosome.".format(assembly_accession, thisfile))

    # CHECK THAT THE SCAFFOLDS ARE THE CORRECT LENGTH AND SAVE TO NEW FILE
    # To keep track of which scaffolds we need to include in the output file, use a copy of the chr_df and delete rows as we go.
    outfile = os.path.join(output_dir, "{}.chr.fasta".format(assembly_accession))
    outhandle = open(outfile, "w")
    scafs_to_parse_df = chr_df.copy()
    # Iterate through the files in the directory and:
    #  - go through each entry in the fasta file
    #    - check if the scaffold is in scafs_to_parse_df
    #    - if it is, then check that the length is correct
    #    - print it out to the outhandle
    #    - remove the corresponding row from scafs_to_parse_df
    for thischrom in chr_df["chr_name"].tolist():
        thisfile = os.path.join(unzip_dir, "ncbi_dataset/data/{}/chr{}.fna".format(assembly_accession,thischrom))
        for record in fasta.parse(thisfile):
            # get the corresponding row from scafs_to_parse_df. could be either "genbank_accession" or "refseq_accession"
            if record.id in scafs_to_parse_df["genbank_accession"].tolist() + scafs_to_parse_df["refseq_accession"].tolist():
                # there is a row in the dataframe that corresponds to this scaffold
                # Get all the indices that correspond to this scaffold, using either the "genbank_accession" or "refseq_accession" column
                matching_indices = scafs_to_parse_df.index[scafs_to_parse_df["genbank_accession"] == record.id].tolist() + scafs_to_parse_df.index[scafs_to_parse_df["refseq_accession"] == record.id].tolist()
                if len(matching_indices) > 1:
                    # There should only be one row that matches this scaffold. If there are more, then something went wrong.
                    raise ValueError("There are multiple rows in the dataframe that match this scaffold. Something went wrong (Presumably with NCBI datasets, as this dataframe was generated with 'datasets summary genome accession ACC --report sequence --as-json-lines'). {}".format(record.id))
                # If we get here, then there is only one row that matches this scaffold. Get the index of that row.
                rowix = matching_indices[0]
                # Get the expected length of this scaffold, based on this row's "length" column
                expected_length = scafs_to_parse_df.loc[rowix, "length"]
                # if the expected length doesn't match the actual length, then something went wrong.
                # THIS IS COMMENTED OUT UNTIL https://github.com/ncbi/datasets/issues/301 IS FIXED
                #if not expected_length == len(record.seq):
                #    raise ValueError("The expected length of scaffold {} from assembly accession {} is {}, but the actual length is {}. Something went wrong.".format(record.id, assembly_accession, expected_length, len(record.seq)))
                # After passing these checks, safe to print to the output file
                print(record.format(wrap=80), file=outhandle, end="")
                # remove the row from scafs_to_parse_df, reset the index
                scafs_to_parse_df = scafs_to_parse_df.drop(rowix).reset_index(drop=True)
            else:
                # This scaffold is not in the dataframe. We shouldn't be here, but maybe one unrequested scaffold was downloaded by the NCBI datasets tool.
                # One instance in which this happens often is when we tell datasets to download specific chromosomes, but it also downloads separate fasta files for the unplaced scaffolds.
                pass
    # close the output file
    outhandle.close()
    # check if there are any rows left in scafs_to_parse_df. If there are, then we didn't see one chromosome's sequence for some reason.
    if not scafs_to_parse_df.empty:
        # remove the output fasta file, since it is broken
        os.remove(outfile)
        # remove the unzipped directory and the zip file
        os.remove(outzip)
        shutil.rmtree(unzip_dir)
        raise ValueError("There are still rows in scafs_to_parse_df. This means that we didn't see one of the chromosomes. Something went wrong. {}".format(scafs_to_parse_df))
    # Just to be completely sure that we don't have any doubled scaffolds appearing in the output fasta file, read through it again with the fasta package
    seen_once = []
    for record in fasta.parse(outfile):
        if record.id in seen_once:
            raise ValueError("The scaffold {} appears twice in the output fasta file. Something went wrong.".format(record.id))
        else:
            seen_once.append(record.id)
    # Another excessive check, just make sure that the length of seen_once is the same length as the number of rows in the chr_df
    if not len(seen_once) == len(chr_df):
        raise ValueError("The number of scaffolds in the output fasta file is not the same as the number of scaffolds in the input dataframe. Something went wrong. Assembly acession: {}".format(assembly_accession))
    # If we get here, then everything probably worked as we intended. Remove the unneeded files.
    os.remove(outzip)
    shutil.rmtree(unzip_dir)
    return 0

def download_unzip_genome(assembly_accession, output_dir, datasetsEx,
                          chrscale = True, threads = 1):
    """
    Given an assembly accession number, download the genome from NCBI.
    If the assembly is chromosome-scale, then we:
      - save the NCBI sequence report of the assembly, and the chromosome-scale rows of it
      - write the chromosome-scale scaffolds to a gzipped fasta file
      - check that the fasta file contains all of the chromosome-scale scaffolds, each exactly once

    Inputs:
      - assembly_accession:   The assembly accession number of the genome to download
      - output_dir:           The directory to which the files will be saved.
                              Everything will be prefixed with the assembly_accession value.
      - datasetsEx:           path to the NCBI datasets executable
      - chrscale:             boolean, whether the genome is chromosome-scale or not. Invokes different processing.
      - threads:              number of compression threads for the fasta file

    Returns:
      - If a chromosome-scale genome, these files will be in the output_dir:
          {assembly_accession}.scaffold_df.all.tsv, {assembly_accession}.scaffold_df.chr.tsv, {assembly_accession}.chr.fasta.gz
      - Returns 0 if the function completes successfully.
    """
    # First, strip the white space surrounding the assembly accession number
    assembly_accession = assembly_accession.strip()

    if chrscale:
        scaffold_df = download_assembly_scaffold_df(assembly_accession, datasetsEx, chrscale = True)
        # we're going to save a file to a directory. First, make sure that the directory exists.
        create_directories_recursive_notouch(output_dir)
        # Originally, I was going to first check if the scaffold_df was already saved to output_dir.
        #  Then, I thought this would be a bad idea because the entry could change on NCBI's end, in which case it is best to download the scaffold info every time.
        scaffold_df_filepath = os.path.join(output_dir, assembly_accession + ".scaffold_df.all.tsv")
        # save the scaffold_df to a file
        scaffold_df.to_csv(scaffold_df_filepath, sep="\t", index=False)
        # Now we just get the chromosome-scale scaffolds
        chr_df = filter_scaffold_df_keep_chrs(scaffold_df)
        chr_df_filepath = os.path.join(output_dir, assembly_accession + ".scaffold_df.chr.tsv")
        chr_df.to_csv(chr_df_filepath, sep="\t", index=False)
        # Download the genome and keep only the chromosome-scale scaffolds
        out_fasta_gz = os.path.join(output_dir, assembly_accession + ".chr.fasta.gz")
        download_chr_scale_genome(chr_df, datasetsEx, out_fasta_gz, threads = threads)
        # check that there is a file in the output
        if not os.path.exists(out_fasta_gz):
            raise IOError("The output fasta file does not exist. Something went wrong with the download. {}".format(out_fasta_gz))
    else:
        # The genome is not chromsome-scale, or we want to download every scaffold.
        pass
    return 0

def contains_date(string_to_check):
    """
    From the string in question, just see if it contains a date in the format YYYYMMDDHHMM or YYYYMMDD
    """
    # split the string on the underscore
    split_string = string_to_check.replace(".tsv","").split("_")
    # check if any of the elements are a date
    string_contains_date = False
    for element in split_string:
        if element.isdigit() and (len(element) == 12 or len(element) == 8):
            # check that the years, months, days are legal
            if not int(element[0:4]) > 1990:
                raise IOError("The year in the date is not legal. {}".format(string_to_check))
            if not (int(element[4:6]) >= 1 and int(element[4:6]) <= 12):
                raise IOError("The month in the date is not legal. {}".format(string_to_check))
            if not (int(element[6:8]) >= 1 and int(element[6:8]) <= 31):
                raise IOError("The day in the date is not legal. {}".format(string_to_check))
            string_contains_date = True
    return string_contains_date

def _best_date_key_from_name(name: str):
    """
    Extract the latest date-like token from the filename and convert to a sortable key.
    Supports YYYYMMDD (8) and YYYYMMDDHHMM (12). For 8-digit dates, assume 23:59 so
    they correctly sort after any same-day times.
    Returns a tuple (YYYY, MM, DD, HH, MM) or None if no date is found.
    """
    tokens = re.findall(r'(?<!\d)((?:19|20)\d{6}(?:\d{4})?)(?!\d)', name)
    best = None
    for t in tokens:
        if len(t) == 12:  # YYYYMMDDHHMM
            key = (int(t[0:4]), int(t[4:6]), int(t[6:8]), int(t[8:10]), int(t[10:12]))
        elif len(t) == 8:  # YYYYMMDD
            key = (int(t[0:4]), int(t[4:6]), int(t[6:8]), 23, 59)
        else:
            continue
        if best is None or key > best:
            best = key
    return best

def return_latest_accession_tsvs(directory_path):
    """
    Return full paths to the most recent annotated/unannotated, chr/nonchr TSVs.
    - Files are detected by flexible name checks (e.g., 'annotated_chr', 'unannotated_nonchr').
    - Presence of a date is required (via contains_date(f)).
    - Recency is determined by the latest date token found in the filename.
    """
    if not os.path.isdir(directory_path):
        raise IOError(f"The directory of TSV files you provided does not exist. {directory_path}")

    files = [f for f in os.listdir(directory_path) if f.endswith(".tsv")]

    # Category filters (robust to extra suffixes like _fewerColumns, etc.)
    annotated_chr_files      = [f for f in files if f.startswith("annotated_")   and "_chr_"    in f and contains_date(f)]
    annotated_nonchr_files   = [f for f in files if f.startswith("annotated_")   and "_nonchr_" in f and contains_date(f)]
    unannotated_chr_files    = [f for f in files if f.startswith("unannotated_") and "_chr_"    in f and contains_date(f)]
    unannotated_nonchr_files = [f for f in files if f.startswith("unannotated_") and "_nonchr_" in f and contains_date(f)]

    cycledict = {
        "annotated_chr_files":      annotated_chr_files,
        "annotated_nonchr_files":   annotated_nonchr_files,
        "unannotated_chr_files":    unannotated_chr_files,
        "unannotated_nonchr_files": unannotated_nonchr_files,
    }

    # Warn (non-fatal) if a category is missing
    for label, flist in cycledict.items():
        if not flist:
            warnings.warn(f"There are no accession TSV files for this filetype: {label}\n"
                          f"Directory checked: {directory_path}")

    # Sort by latest embedded date key, newest first
    def sort_key(fname):
        key = _best_date_key_from_name(fname)
        # Put files without a parsable key at the very beginning (then reversed -> end)
        return key if key is not None else (0, 0, 0, 0, 0)

    for key in cycledict:
        cycledict[key].sort(key=sort_key, reverse=True)

    # Pick most recent or None
    most_recent_annotated_chr_file      = cycledict["annotated_chr_files"][0]      if cycledict["annotated_chr_files"]      else None
    most_recent_annotated_nonchr_file   = cycledict["annotated_nonchr_files"][0]   if cycledict["annotated_nonchr_files"]   else None
    most_recent_unannotated_chr_file    = cycledict["unannotated_chr_files"][0]    if cycledict["unannotated_chr_files"]    else None
    most_recent_unannotated_nonchr_file = cycledict["unannotated_nonchr_files"][0] if cycledict["unannotated_nonchr_files"] else None

    return {
        "annotated_chr":      os.path.join(directory_path, most_recent_annotated_chr_file)      if most_recent_annotated_chr_file      else None,
        "annotated_nonchr":   os.path.join(directory_path, most_recent_annotated_nonchr_file)   if most_recent_annotated_nonchr_file   else None,
        "unannotated_chr":    os.path.join(directory_path, most_recent_unannotated_chr_file)    if most_recent_unannotated_chr_file    else None,
        "unannotated_nonchr": os.path.join(directory_path, most_recent_unannotated_nonchr_file) if most_recent_unannotated_nonchr_file else None,
    }