#!/bin/bash

readsPath=$1
outputPath=$2
threads=$3
memory=$4

if [ -z "$3" ] 
then
    threads=8
fi

if [ -z "$4" ] 
then
    memory=5000
fi

if [ -f "$readsPath" ]; then
    echo "Reading reads from path = " $readsPath 
else
    echo "Unable to load the reads file = " $readsPath
    exit 1
fi

[ -d $outputPath ] && rm -r $outputPath
mkdir -p $outputPath

echo "Filtering Reads"
#FILTER READS
`dirname $0`/bin/filter "$readsPath" "$outputPath"/filteredReads.fq

echo "Running DSK"
#DSK OPERATIONS
dsk -file "$outputPath"/filteredReads.fq -kmer-size 15 -abundance-min 10 -out-dir "$outputPath"/DSK -max-memory $memory -nb-cores $threads
python `dirname $0`/src/scan-dsk.py "$outputPath"/DSK/filteredReads.h5 "$outputPath"/DSK/ $threads
cat "$outputPath"/DSK/*.chunk > "$outputPath"/DSK/15mersCounts

echo "Counting Trimers"
#GET TRIMERS
mkdir -p "$outputPath"/profiles

`dirname $0`/bin/countTrimers "$outputPath"/filteredReads.fq "$outputPath"/profiles/3mers $threads

echo "Counting 15-mer profiles"
#GET 15-MER PROFILES
`dirname $0`/bin/search15mers "$outputPath"/DSK/15mersCounts "$outputPath"/filteredReads.fq "$outputPath"/profiles/15mers $threads

echo "Obtaining Read Ids"
#EXTRACT READ IDS
awk "NR%4==1" "$outputPath"/filteredReads.fq > "$outputPath"/readIdsFiltered


echo "Sampling Reads"
python `dirname $0`/src/sampledata.py -p3 "$outputPath"/profiles/3mers -p15 "$outputPath"/profiles/15mers -c 10000

echo "Binning sampled reads"
python `dirname $0`/src/Binner.py -p3 "$outputPath"/profiles/3mers_sampled -p15 "$outputPath"/profiles/15mers_sampled -o "$outputPath"/Binning

echo "Assigning reads"
`dirname $0`/bin/assign "$outputPath"/profiles/3mers "$outputPath"/profiles/15mers "$outputPath"/Binning/cluster-stats.txt 8 "$outputPath"/Binning/final.txt 

echo "Final Result = $outputPath/Binning/final.txt"