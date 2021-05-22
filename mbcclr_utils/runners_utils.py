import pickle
import os
from Bio import SeqIO
from tqdm import tqdm
import logging
import sys

from mbcclr_utils import scan_dsk

logger = logging.getLogger('MetaBCC-LR')

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

    cmd = f""""{os.path.dirname(__file__)}/bin/count-kmers" "{reads_path}" "{output}/profiles/3mers" {k_size} {threads}"""
    logger.debug("CMD::" + cmd)
    o = os.system(cmd)
    check_proc(o, "Counting Trimers")

def run_15mers(reads_path, output, bin_size, bins, threads):
    if not os.path.isdir(f"{output}/profiles"):
        os.makedirs(f"{output}/profiles")

    cmd = f""""{os.path.dirname(__file__)}/bin/count-15mers" "{reads_path}" "{output}/profiles/15mers-counts" {threads}"""
    logger.debug("CMD::" + cmd)
    o = os.system(cmd)
    check_proc(o, "Counting 15-mers")
    
    cmd = f""""{os.path.dirname(__file__)}/bin/search-15mers" "{output}/profiles/15mers-counts" "{reads_path}" "{output}/profiles/15mers" {bin_size} {bins} {threads}"""
    logger.debug("CMD::" + cmd)
    o = os.system(cmd)
    check_proc(o, "Counting 15-mer profiles")

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