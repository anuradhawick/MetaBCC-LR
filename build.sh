#!/bin/bash

[ -d mbcclr_utils/bin ] && rm -r mbcclr_utils/bin
mkdir mbcclr_utils/bin

echo "BUILDING THE MetaBCC-LR 15 MER COMPUTATIONS"
g++ mbcclr_utils/search-15mers.cpp -fopenmp -Wall -o mbcclr_utils/bin/search15mers
echo "BUILDING THE MetaBCC-LR 3 MER COMPUTATIONS"
g++ mbcclr_utils/count-tri.cpp -Wall -fopenmp -o mbcclr_utils/bin/countTrimers
echo "BUILDING READ ASSIGNER"
g++ mbcclr_utils/assign_bins.cpp -fopenmp -o mbcclr_utils/bin/assign
echo "BUILD FINISHED"