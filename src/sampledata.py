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
parser.add_argument('-rl',
                    help="Read length of each read",
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
lengthsFile = args.rl
COUNT = args.c

o3 = open(prof3Out, "w+")
o15 = open(prof15Out, "w+")

random_linenos = []
ids = []
sampledIds = []
lengths = []

if idsFile:
    idsOut = idsFile + "_sampled"
    with open(idsFile) as f:
        ids = f.read().strip().split("\n")

with open(lengthsFile) as f:
    lengths = list(map(int, f.read().strip().split("\n")))

if COUNT <= 0:
    COUNT = int(len(lengths)/100)

# zip order of reads and lengths
r_l = list(zip([x for x in range(len(lengths))], lengths))
# reverse sort by lengths; longest at first
srtd = sorted(r_l, key=lambda x: x[1], reverse=True)
# get longest 10%
longestSet = srtd[0:int(len(lengths)/10)]
# sample COUNT amount of reads and resort by their order
longestSet = random.sample(longestSet, COUNT)
srtByRid = sorted(longestSet, key=lambda x: x[0])

with open(prof3) as f:
    lines_copy = list(srtByRid)

    lineno = lines_copy.pop(0)[0]

    for n, line in enumerate(f):
        if n == lineno:
            o3.write(line)
            if idsFile:
                sampledIds.append(re.sub(r'_Read|@|_[0-9]+', '', ids[n]))
            if len(lines_copy) > 0:
                lineno = lines_copy.pop(0)[0]
        elif len(lines_copy) == 0:
            break


    with open(prof15) as f:
        lines_copy = list(srtByRid)
        lineno = lines_copy.pop(0)[0]
        for n, line in enumerate(f):
            if n == lineno:
                o15.write(line)
                if len(lines_copy) > 0:
                    lineno = lines_copy.pop(0)[0]
            elif len(lines_copy) == 0:
                break
if idsFile:
    with open(idsOut, "w+") as ff:
        for x in sampledIds:
            ff.write(x + "\n")

o3.close()
o15.close()