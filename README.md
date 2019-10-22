# MetaBCC-LR
**MetaBCC-LR**: MetaBCC-LR: Reference-free Binning of Long Reads using k-mer Coverage Histograms and Oligonucleotide Composition in Metagenomics Analysis

## Dependencies
MetaBCC-LR is coded purely using C++ (v9) and Python 3.6. To run MetaBCC-LR, you will need to install the following python and C++ modules.

### Python dependencies
* numpy 1.16.4 
* scipy 1.3.0 
* kneed 0.4.2
* seaborn 0.9.0
* h5py 2.9.0

### C++ requirements
* GCC version 9.1.0
* OpenMP 4.5 for multi processing

### Third party programs
* DSK: https://github.com/GATB/dsk
    * Add DSK binaries to your PATH variable

## Downloading MetaBCC-LR
To download MetaBCC-LR, you have to clone the MetaBCC-LR repository to your machine.

```
git clone https://github.com/anuradhawick/MetaBCC-LR.git
```

## Compiling the source code
```
cd MetaBCC-LR
./build.sh
```

## Running the MetaBCC-LR
In order to run MetaBCC-LR you are required to provide the reads in FASTQ format (Only format supported in current version).

```
cd MetaBCC-LR
./MetaBCC-LR -h

usage: MetaBCC-LR [-h] -r <READS PATH> [-t THREADS] [-i IDS] -o DEST

MetaBCC-LR Help. A tool developed for binning of metagenomics long reads
(PacBio/ONT). Tool utilizes composition and coverage profiles of reads based
on k-mer frequencies to perform dimension reduction. dimension reduced reads
are then clustered using DB-SCAN. Minimum RAM requirement is 9GB.

optional arguments:
  -h, --help       show this help message and exit
  -r <READS PATH>  Reads path (FASTQ)
  -t THREADS       Thread limit
  -i IDS           Read ids of reads (For dry runs with ground truth)
  -o DEST          Output directory
```
* Reads path should contain FASTQ reads in the standard FASTQ format.
* Output path is the foldername that you wish the results to be in.
* Specify the number of threads
* The program requires a minimum of 9GB to run. This is because we have optimized the coverage histogram generation process to accommodate all 15mers in RAM for faster lookup of counts.
