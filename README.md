# MetaBCC-LR

![GitHub](https://img.shields.io/github/license/anuradhawick/MetaBCC-LR)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/anuradhawick/MetaBCC-LR)

**MetaBCC-LR**: Metagenomics Binning by Coverage and Composition for Long Reads

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

usage: MetaBCC-LR [-h] -r <READS PATH> [-t <THREADS>] [-i <IDS>]
                  [-c <No of reads to sample>] [-s <Sensitivity>]
                  [-b <bin width for coverage histograms>] [--resume] -o
                  <DEST>

MetaBCC-LR Help. A tool developed for binning of metagenomics long reads
(PacBio/ONT). Tool utilizes composition and coverage profiles of reads based
on k-mer frequencies to perform dimension reduction. dimension reduced reads
are then clustered using DB-SCAN. Minimum RAM requirement is 9GB.

optional arguments:
  -h, --help            show this help message and exit
  -r <READS PATH>       Reads path (FASTQ)
  -t <THREADS>          Thread limit
  -i <IDS>              Read ids of reads (For dry runs with ground truth)
  -c <No of reads to sample>
                        Number of reads to sample in order to determine the
                        number of bins. Set to 1% of reads by default.
                        Changing this parameter will affect whether low
                        coverage species are separated or not.
  -s <Sensitivity>      Value between 1 and 10, Higher helps recovering low
                        abundant species (No. of species > 100)
  -b <bin width for coverage histograms>
                        Value of bx32 will be the total coverage of k-mers in
                        the coverage histograms. Usually k-mers are shifted
                        towards y-axis due to errors. By defaul b=10;
                        coverages upto 320X
  --resume              Continue from the last step or the binning step (which
                        ever comes first). Can save time needed to run DSK and
                        obtain k-mers. Ideal for sensitivity tuning
  -o <DEST>             Output directory
```
* Reads path should contain FASTQ reads in the standard FASTQ format.
* Output path is the foldername that you wish the results to be in.
* Specify the number of threads
* The program requires a minimum of 9GB to run. This is because we have optimized the coverage histogram generation process to accommodate all 15mers in RAM for faster lookup of counts.
