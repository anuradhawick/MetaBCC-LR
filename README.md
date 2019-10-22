# MetaBCC-LR
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
TBD
```
* Reads path should contain FASTQ reads in the standard FASTQ format.
* Output path is the folder name that you wish the results to be saved in.
* Specify the number of threads.
* Max memory for the DSK to run.
