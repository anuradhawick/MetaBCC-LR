#!/bin/bash

echo "STARTED BUILDING C++ COMPONENTS"
echo "USING COMPILER"
echo `which g++`
echo `g++ -v`

[ -d bin ] && rm -r bin
mkdir bin

echo "BUILDING THE MetaBCC-LR 15 MER COMPUTATIONS"
g++ src/search-15mers.cpp -fopenmp -Wall -o bin/search15mers
echo "BUILDING THE MetaBCC-LR 3 MER COMPUTATIONS"
g++ src/count-tri.cpp -Wall -fopenmp -o bin/countTrimers
echo "BUILDING THE MetaBCC-LR READ FILTER"
g++ src/filter_reads.cpp -Wall -fopenmp -o bin/filter
echo "BUILDING READ ASSIGNER"
g++ src/assignBins.cpp -fopenmp -o bin/assign

echo "BUILD FINISHED"