<p align="center">
  <img src="MetaBCC-LR_logo.png" width="600" title="Final Labelling" alt="Final Labelling">
</p>

# MetaBCC-LR: Metagenomics Binning by Coverage and Composition for Long Reads

![GitHub](https://img.shields.io/github/license/anuradhawick/MetaBCC-LR)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/anuradhawick/MetaBCC-LR)

## Dependencies
MetaBCC-LR is coded purely using C++ (v9) and Python 3.6. To run MetaBCC-LR, you will need to install the following python and C++ modules.

### Python dependencies
* numpy 1.16.4 
* scipy 1.3.0 
* kneed 0.4.2
* seaborn 0.9.0
* h5py 2.9.0
* tabulate 0.8.7

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
* Build the binaries
```
cd MetaBCC-LR
python setup.py build
```
OR
```
sh build.sh
```    
* To install the program 
```
pip install .
```
OR add the program path to your $PATH variable.

## Running the MetaBCC-LR
In order to run MetaBCC-LR you are required to provide the reads in FASTQ or FASTA format.

```
cd MetaBCC-LR
./MetaBCC-LR -h

usage: MetaBCC-LR [-h] --reads-path READS_PATH [--threads THREADS]
                  [--ground-truth GROUND_TRUTH] [--sample-count SAMPLE_COUNT]
                  [--sensitivity SENSITIVITY] [--bin-size BIN_SIZE]
                  [--max-memory MAX_MEMORY] [--resume] --output <DEST>
                  [--version]

MetaBCC-LR Help. A tool developed for binning of metagenomics long reads
(PacBio/ONT). Tool utilizes composition and coverage profiles of reads based
on k-mer frequencies to perform dimension reduction. dimension reduced reads
are then clustered using DB-SCAN. Minimum RAM requirement is 9GB.

optional arguments:
  -h, --help            show this help message and exit
  --reads-path READS_PATH, -r READS_PATH
                        Reads path for binning
  --threads THREADS, -t THREADS
                        Thread count for computation
  --ground-truth GROUND_TRUTH, -g GROUND_TRUTH
                        Ground truth of reads for dry runs and sensitivity
                        tuning
  --sample-count SAMPLE_COUNT, -c SAMPLE_COUNT
                        Number of reads to sample in order to determine the
                        number of bins. Set to 1% of reads by default.
                        Changing this parameter will affect whether low
                        coverage species are separated or not.
  --sensitivity SENSITIVITY, -s SENSITIVITY
                        Value between 1 and 10, Higher helps recovering low
                        abundant species (No. of species > 100)
  --bin-size BIN_SIZE, -b BIN_SIZE
                        Value of bx32 will be the total coverage of k-mers in
                        the coverage histograms. Usually k-mers are shifted
                        towards y-axis due to errors. By defaul b=10;
                        coverages upto 320X
  --max-memory MAX_MEMORY, -m MAX_MEMORY
                        Default 5000. DSK k-mer counter accepts a max memory
                        parameter. However, the complete pipeline requires
                        5GB+ RAM. This is only to make DSK step faster, should
                        you have more RAM.
  --resume              Continue from the last step or the binning step (which
                        ever comes first). Can save time needed to run DSK and
                        obtain k-mers. Ideal for sensitivity tuning
  --output <DEST>, -o <DEST>
                        Output directory
  --version, -v         Show version.
```
* Output path is the foldername that you wish the results to be in.
* Specify the number of threads
* The program requires a minimum of 5GB to run. This is because we have optimized the coverage histogram generation process to accommodate all 15mers in RAM for faster lookup of counts.

## Citation

```bibtex
@article{10.1093/bioinformatics/btaa441,
    author = {Wickramarachchi, Anuradha and Mallawaarachchi, Vijini and Rajan, Vaibhav and Lin, Yu},
    title = "{MetaBCC-LR: metagenomics binning by coverage and composition for long reads}",
    journal = {Bioinformatics},
    volume = {36},
    number = {Supplement_1},
    pages = {i3-i11},
    year = {2020},
    month = {07},
    abstract = "{Metagenomics studies have provided key insights into the composition and structure of microbial communities found in different environments. Among the techniques used to analyse metagenomic data, binning is considered a crucial step to characterize the different species of micro-organisms present. The use of short-read data in most binning tools poses several limitations, such as insufficient species-specific signal, and the emergence of long-read sequencing technologies offers us opportunities to surmount them. However, most current metagenomic binning tools have been developed for short reads. The few tools that can process long reads either do not scale with increasing input size or require a database with reference genomes that are often unknown. In this article, we present MetaBCC-LR, a scalable reference-free binning method which clusters long reads directly based on their k-mer coverage histograms and oligonucleotide composition.We evaluate MetaBCC-LR on multiple simulated and real metagenomic long-read datasets with varying coverages and error rates. Our experiments demonstrate that MetaBCC-LR substantially outperforms state-of-the-art reference-free binning tools, achieving ∼13\\% improvement in F1-score and ∼30\\% improvement in ARI compared to the best previous tools. Moreover, we show that using MetaBCC-LR before long-read assembly helps to enhance the assembly quality while significantly reducing the assembly cost in terms of time and memory usage. The efficiency and accuracy of MetaBCC-LR pave the way for more effective long-read-based metagenomics analyses to support a wide range of applications.The source code is freely available at: https://github.com/anuradhawick/MetaBCC-LR.Supplementary data are available at Bioinformatics online.}",
    issn = {1367-4803},
    doi = {10.1093/bioinformatics/btaa441},
    url = {https://doi.org/10.1093/bioinformatics/btaa441},
    eprint = {https://academic.oup.com/bioinformatics/article-pdf/36/Supplement\_1/i3/33488763/btaa441.pdf},
}
```

## Notes

A setup and a bundler with few performance updates will be released in due course.
