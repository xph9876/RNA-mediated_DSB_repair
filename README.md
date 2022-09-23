# RNA-mediated DSB repair

## Description

Contains the code for the computational analysis presented in the [publication](#citation).

## Citation:

Jeon Y. *et al*. RNA-mediated double-strand break repair in human cells. (2022). *In preparation*.

## Files

* *scripts*: Python scripts for microhomology mediated end joining (MMEJ) analysis, read trimming, and read categorization.

* *refseq*: Reference sequences for MMEJ analysis.

* *analyze_MMEJ.sh* : Main script for generating figures for MMEJ analysis.

* *analyze_MMEJ_antisense.sh*: Main script for generating figures for MMEJ analysis (antisense libraries).

* *libinfo.tsv*: Library metadata for MMEJ analysis.

* *libinfo_antisense.tsv*: Library metadata for MMEJ analysis (antisense libraries).

* *analyze_nhej*: Scripts for nonhomologous end joining analysis (NHEJ). Plotting variation-distance graphs and variation-position histograms. See *analyze_nhej/README.md* for more details.

## Dependencies

The trimming and read-categorizing scripts do not have any dependencies except for Python VERSION.

MMEJ pipeline:
* Python (VERSION)
* scipy (VERSION)
* seaborn (VERSION)
* pandas (VERSION)
* xlsxwriter (VERSION)
* matplotlib (VERSION)
* statannotations (VERSION)

NHEJ pipeline: please see *analyze_nhej/README.txt*.

## Steps

### Trimming

The trimming script in *scripts/trim_F_tag.py* (forward strand file) and *scripts/trim_R_tag.py* (reverse strand files) should be run on the raw Illumina FASTQ files before all other steps.

### Read categorization

TODO

### MMEJ pipeline

1) Place the trimmed FASTQ files into the appropriate directory with the appropriate name (see comments in *analyze_MMEJ.sh* and *analyze_MMEJ_antisense.sh*).

2) Run *analyze_MMEJ.sh* and *analyze_MMEJ_antisense.sh*.

### NHEJ pipeline

1) Align the trimmed FASTQ files with the appropriate reference sequence using Bowtie2 (see *analyze_nhej/README.md* for examples).

2) Place SAM file output of alignment in the appropriate directory with appropriate name (see *analyze_nhej/README.md*).

3) Run either *run_all.ps1* (Windows) or *run_all.sh* (Unix). Note, the working directory must be *analyze_nhej* for this to work correctly.

See *analyze_nhej/README.md* for more details.