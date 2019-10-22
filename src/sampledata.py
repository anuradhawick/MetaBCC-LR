import sys
import random
from datetime import datetime
import re
import argparse

parser = argparse.ArgumentParser(description='Sample profiles for analysis.')

parser.add_argument('-p3',
                    help="Path to trimer profiles",
                    type=str,
                    required=True)
parser.add_argument('-p15',
                    help="Path to 15-mer coverage histogram profiles",
                    type=str,
                    required=True)
parser.add_argument('-ids',
                    help="Read ids of reads (For dry runs with ground truth)",
                    type=str,
                    required=False,
                    default=None)                    

parser.add_argument('-c',
                    help="Number of profiles to sample",
                    type=int,
                    required=True)  

args = parser.parse_args()


prof3 = args.p3
prof15 = args.p15
prof3Out = prof3 + "_sampled"
prof15Out = prof15 + "_sampled"
idsFile = args.ids
COUNT = args.c

o3 = open(prof3Out, "w+")
o15 = open(prof15Out, "w+")

random_linenos = []
ids = []
sampledIds = []

if idsFile:
    idsOut = idsFile + "_sampled"
    with open(idsFile) as f:
        ids = f.read().split()

with open(prof3) as f:
    linecount = sum(1 for line in f)
    f.seek(0)
    lst = list(range(0, linecount))
    random.shuffle(lst)
    random_linenos = sorted(random.sample(lst, COUNT), reverse = True)
    lines_copy = list(random_linenos)
    lineno = lines_copy.pop()

    for n, line in enumerate(f):
        if n == lineno:
            o3.write(line)
            if idsFile:
                sampledIds.append(re.sub(r'_Read|@|_[0-9]+', '', ids[n]))
            if len(lines_copy) > 0:
                lineno = lines_copy.pop()
        elif len(lines_copy) == 0:
            break


    with open(prof15) as f:
        lines_copy = list(random_linenos)
        lineno = lines_copy.pop()
        for n, line in enumerate(f):
            if n == lineno:
                o15.write(line)
                if len(lines_copy) > 0:
                    lineno = lines_copy.pop()
            elif len(lines_copy) == 0:
                break
if idsFile:
    with open(idsOut, "w+") as ff:
        for x in sampledIds:
            ff.write(x + "\n")