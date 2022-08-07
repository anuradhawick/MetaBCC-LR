<p align="center">
  <img src="MetaBCC-LR_logo.png" width="600" title="Final Labelling" alt="Final Labelling">
</p>

# MetaBCC-LR: Metagenomics Binning by Coverage and Composition for Long Reads

![GitHub](https://img.shields.io/github/license/metagentools/MetaBCC-LR)
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
* umap-learn 0.5.1
* song-vis (latest version from github)

### C++ requirements
* GCC version 9.1.0
* OpenMP 4.5 for multi processing
* PThreads (any version should work)

<!-- ### Third party programs
* DSK: https://github.com/GATB/dsk
    * Add DSK binaries to your PATH variable -->

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

### Test run data

Extract test data from [here](https://anu365-my.sharepoint.com/:f:/g/personal/u6776114_anu_edu_au/EnV-rUq01pRHl1lH4Y8SaSwBwVVMKNAptbA6YW8RWX6Pqw?e=tDgy9v);

In order to run MetaBCC-LR you are required to provide the reads in FASTQ or FASTA format.

```
python mbcclr --resume -r test_data/data/reads.fasta -g test_data/data/ids.txt -o test_output -e umap -c 25000 -bs 10 -bc 10 -k 4
```

### Separate reads into Bins

You can use the script `reads2bins.py` to separate reads into bins. This is included in a separate script as you might want to play around with clustering sensitivity and sampling reads count to get a good final binning. You can look into images generated in `Output/images` directory to see if you have a good clustering of reads. Finally you can use the script `reads2bins.py` to separate reads.

Inputs:
* -r path to reads file used for binning
* -b output/final.txt (the file containing bin of each read)
* -o a destination directory to place final fasta files

```
usage: reads2bins.py [-h] --reads READS --bins BINS --output OUTPUT

Separate reads in to bins.

optional arguments:
  -h, --help            show this help message and exit
  --reads READS, -r READS
  --bins BINS, -b BINS
  --output OUTPUT, -o OUTPUT
```

### Usage and Help
```
cd MetaBCC-LR
./mbcclr -h

usage: mbcclr [-h] --reads-path READS_PATH [--embedding {tsne,umap,song}]
              [--k-size {3,4,5,6,7}] [--sample-count SAMPLE_COUNT]
              [--sensitivity SENSITIVITY] [--bin-size BIN_SIZE]
              [--bin-count BIN_COUNT] [--threads THREADS]
              [--ground-truth GROUND_TRUTH] [--resume] --output OUTPUT
              [--version]

MetaBCC-LR Help. A tool developed for binning of metagenomics long reads
(PacBio/ONT). Tool utilizes composition and coverage profiles of reads based
on k-mer frequencies to perform dimension reduction. dimension reduced reads
are then clustered using DB-SCAN. Minimum RAM requirement is 9GB.

optional arguments:
  -h, --help            show this help message and exit
  --reads-path READS_PATH, -r READS_PATH
                        Reads path for binning
  --embedding {tsne,umap,song}, -e {tsne,umap,song}
                        Embedding tool to be used for clustering
  --k-size {3,4,5,6,7}, -k {3,4,5,6,7}
                        Choice of k-mer for oligonucleotide frequency vector.
  --sample-count SAMPLE_COUNT, -c SAMPLE_COUNT
                        Number of reads to sample in order to determine the
                        number of bins. Set to 1% of reads by default.
                        Changing this parameter will affect whether low
                        coverage species are separated or not.
  --sensitivity SENSITIVITY, -s SENSITIVITY
                        Value between 1 and 10, Higher helps recovering low
                        abundant species (No. of species > 100)
  --bin-size BIN_SIZE, -bs BIN_SIZE
                        Size of each bin in coverage histogram.
  --bin-count BIN_COUNT, -bc BIN_COUNT
                        Number of bins in the coverage histogram.
  --threads THREADS, -t THREADS
                        Thread count for computation
  --ground-truth GROUND_TRUTH, -g GROUND_TRUTH
                        Ground truth of reads for dry runs and sensitivity
                        tuning
  --resume              Continue from the last step or the binning step (which
                        ever comes first). Can save time needed to run DSK and
                        obtain k-mers. Ideal for sensitivity tuning
  --output OUTPUT, -o OUTPUT
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

## Update

Program can be built and installed with `sh build` and `pip install .` 
We recommend using `sh build` and using the program without installing. Thus making it easier to fetch future upadates and run.

## New in v-2.X

* No need to have DSK, we have implemented a consice k-mer counting strategy using compare and swap (CAS).
* Supports UMAP and SONG embeddings. Please note that UMAP and SONG are still being improved. Needs more work from our side. But usable!
* Supports any input format **fasta, fastq** or **gzipped** formats of either. (Thanks for Klib by Attractive Chaos [blog](http://attractivechaos.github.io/klib/))
