import argparse
import os
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Separate reads in to bins.')

parser.add_argument('--reads', '-r', type=str, required=True)
parser.add_argument('--bins', '-b', type=str, required=True)
parser.add_argument('--output', '-o' ,type=str, required=True)

args = parser.parse_args()

readsPath = args.reads
readBinsPath = args.bins
readsType = readsPath.split(".")[-1].lower()

if readsType in ["fasta", "fna", "fa"]:
    readsType = "fasta"
else:
    readsType = "fastq"

if not os.path.exists(args.output):
    os.makedirs(args.output)

outputDir = args.output + "/"
outFiles = {}
bins = open(readBinsPath)

for record, bin_id in zip(SeqIO.parse(readsPath), open(readBinsPath)):
    bin_id = bin_id.strip().split("\n")[-1]
    bpath = f"{outputDir}/{bin_id}.fasta"

    if bpath not in outFiles:
        outFiles[bpath] = open(bpath, "w+")

    SeqIO.write(record, outFiles[bpath], "fasta")

for k, f in outFiles.items():
    f.close()
bins.close()