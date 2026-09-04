import pickle
import os
from Bio import SeqIO
from tqdm import tqdm
import logging
import sys
import shutil
import subprocess
import tempfile
from pykmertools import OligoComputer

from mbcclr_utils import scan_dsk

logger = logging.getLogger('MetaBCC-LR')

# Depprecated
def run_filter(reads_path, output, ids=None):
    extension = reads_path.split(".")[-1]

    if extension in ["fq", "fastq"]:
        fmt = "fastq"
    else:
        fmt = "fasta"
    
    records = 0
    for record in SeqIO.parse(reads_path, fmt):
        records += 1
    
    logger.debug(f"Total of {records} reads to filter")
    output_fasta_file = open(f"{output}/misc/filtered_reads.fasta", "w+")

    if ids:
        output_truth_file = open(f"{output}/misc/filtered_truth.txt", "w+")
        for record, truth in tqdm(zip(SeqIO.parse(reads_path, fmt), open(ids)), total=records, desc="Filtering reads longer than 1000bp"):
            if len(record.seq) >= 1000:
                output_truth_file.write(truth)
                output_fasta_file.write(f">{record.id}\n{str(record.seq)}\n")
        output_truth_file.close()
    else:
        for record in tqdm(SeqIO.parse(reads_path, fmt), total=records, desc="Filtering reads longer than 1000bp"):
            if len(record.seq) >= 1000:
                output_fasta_file.write(f">{record.id}\n{str(record.seq)}\n")
    output_fasta_file.close()
    
def run_assign(output, threads):
    cmd = f"""{os.path.dirname(__file__)}/bin/assign "{output}/profiles/3mers" "{output}/profiles/15mers" "{output}/misc/cluster-stats.txt" {threads} {output}/final.txt """
    o = os.system(cmd)
    check_proc(o, "Assigning reads")    

def run_kmers(reads_path, output, k_size, threads):
    if not os.path.isdir(f"{output}/profiles"):
        os.makedirs(f"{output}/profiles")

    if reads_path.split(".")[-1].lower() in ["fq", "fastq"]:
        fmt = "fastq"
    else:
        fmt = "fasta"

    output_path = f"{output}/profiles/3mers"
    computer = OligoComputer(k_size)
    fallback = [0.0 for _ in computer.get_header(mins=True)]

    with open(reads_path, "r") as input_file, open(output_path, "w+") as output_file:
        for record in SeqIO.parse(input_file, fmt):
            seq = str(record.seq).upper()
            if len(seq) < k_size:
                profile = fallback
            else:
                try:
                    profile = computer.vectorise_one(seq, norm=True, mins=True)
                except ValueError:
                    seq = "".join([s for s in seq if s in "ACGT"])
                    if len(seq) < k_size:
                        profile = fallback
                    else:
                        try:
                            profile = computer.vectorise_one(seq, norm=True, mins=True)
                        except ValueError:
                            profile = fallback
            output_file.write(" ".join([f"{v:.6f}" for v in profile]) + "\n")

def run_15mer_counts(reads_path, output, threads):
    if not os.path.isdir(f"{output}/profiles"):
        os.makedirs(f"{output}/profiles")
    logger.debug("15-mer counts are generated with kmertools cov in run_15mer_vecs")

def run_15mer_vecs(reads_path, output, bin_size, bins, threads):
    if not os.path.isdir(f"{output}/profiles"):
        os.makedirs(f"{output}/profiles")

    kmertools = shutil.which("kmertools")
    if kmertools is None:
        local_cli = os.path.join(os.path.dirname(sys.executable), "kmertools")
        if os.path.isfile(local_cli) and os.access(local_cli, os.X_OK):
            kmertools = local_cli

    if kmertools is None:
        logger.error("Unable to locate kmertools CLI. Please install pykmertools in this environment.")
        sys.exit(1)

    profiles_path = f"{output}/profiles"
    with tempfile.TemporaryDirectory(prefix="kmertools-cov-", dir=profiles_path) as tmpdir:
        cmd = [
            kmertools,
            "cov",
            "-i", reads_path,
            "-o", tmpdir,
            "-k", "15",
            "-s", str(bin_size),
            "-c", str(bins),
            "-p", "spc",
            "-t", str(threads),
        ]
        logger.debug("CMD::" + " ".join(cmd))
        proc = subprocess.run(cmd)
        if proc.returncode != 0:
            check_proc(proc.returncode, "Counting 15-mer profiles")

        counts_path = os.path.join(tmpdir, "kmers.counts")
        vectors_path = os.path.join(tmpdir, "kmers.vectors")

        if not os.path.isfile(counts_path) or not os.path.isfile(vectors_path):
            logger.error("kmertools cov did not produce expected outputs: kmers.counts and kmers.vectors")
            sys.exit(1)

        shutil.move(counts_path, f"{profiles_path}/15mers-counts")
        shutil.move(vectors_path, f"{profiles_path}/15mers")

# depprecated
def run_dsk(output, max_memory, threads):
    logger.debug("Running DSK")
    cmdDSK = f"""dsk -verbose 0 -file "{output}/misc/filtered_reads.fasta" -kmer-size 15 -abundance-min 10 -out-dir "{output}/misc/DSK" -max-memory {max_memory} -nb-cores {threads}"""
    logger.debug("CMD::" + cmdDSK)
    o = os.system(cmdDSK)
    check_proc(o, "Running DSK")
    scan_dsk.scan_dsk(f"{output}/misc/DSK/filtered_reads.h5", threads, f"{output}/misc/DSK/")

def checkpoint(new_checkpoints, path):
    pickle.dump(new_checkpoints, open(path, "wb+"))

def load_checkpoints(path):
    if not os.path.isfile(path):
        data = {}
        data['completed'] = set()
        return data
    else:
        return pickle.load(open(path, "rb"))

def check_proc(ret, name=""):
    if ret != 0:
        if name!= "": logger.error(f"Error in step: {name}")
        logger.error("Failed due to an error. Please check the log. Good Bye!")
        sys.exit(ret)
