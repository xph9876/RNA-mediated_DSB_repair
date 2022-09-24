# RNA-mediated DSB repair

## Description

Contains the code for the computational analysis presented in the [publication](#citation).

## Citation

Jeon Y. *et al*. RNA-mediated double-strand break repair in human cells. (2022). *In preparation*.

## Files/Directories

* `scripts`: Python scripts for microhomology mediated end joining (MMEJ) analysis, read trimming, and read categorization.

* `refseq`: Reference sequences for MMEJ analysis.

* `analyze_MMEJ.sh` : Main script for generating figures for MMEJ analysis.

* `analyze_MMEJ_antisense.sh`: Main script for generating figures for MMEJ analysis (antisense libraries).

* `libinfo.tsv`: Library metadata for MMEJ analysis.

* `libinfo_antisense.tsv`: Library metadata for MMEJ analysis (antisense libraries).

* `analyze_nhej`: Scripts for nonhomologous end joining analysis (NHEJ). Plotting variation-distance graphs and variation-position histograms. See `analyze_nhej/README.md` for more details.

## Dependencies

The trimming and read-categorizing scripts do not have any dependencies other than Python 3.8.13.

MMEJ pipeline:

* Python 3.8.13
* Numpy 1.22.3
* Matplotlib 3.4.3
* Pandas 1.4.2
* Seaborn 0.11.2
* Xlsxwriter 3.0.3
* Scipy 1.7.3
* Statannotations 0.4.3

NHEJ pipeline: please see `analyze_nhej/README.txt`.

## Steps

### Trimming

The trimming scripts `scripts/trim_F_tag.py` (forward strand file) and `scripts/trim_R_tag.py` (reverse strand files) should be run on the raw Illumina FASTQ files before all other steps.

### Read Categorization

The reads are categorized into different repair mechanisms by `scripts/FIXME`. The input should be the trimmed FASTQ files from [trimming](#trimming).

### MMEJ Pipeline

1) Place the [trimmed](#trimming) FASTQ files into the appropriate directory with the appropriate name (see comments in `analyze_MMEJ.sh` and `analyze_MMEJ_antisense.sh`).

2) Run `analyze_MMEJ.sh` and `analyze_MMEJ_antisense.sh`.

### NHEJ Pipeline

1) Align the [trimmed](#trimming) FASTQ files with the appropriate reference sequence using Bowtie2 (see `analyze_nhej/README.md` for examples).

2) Place SAM file output of alignment in the appropriate directory with appropriate name (see `analyze_nhej/README.md`).

3) Run either `analyze_nhej/run_all.ps1` (Windows) or `analyze_nhej/run_all.sh` (Unix). Note, the working directory must be `analyze_nhej` for this to work correctly.

See `analyze_nhej/README.md` for more details.

## Contact

Webpages:
* [Francesca Storici Lab](https://storicilab.gatech.edu/)
* [Natasha Jonoska Group](https://knot.math.usf.edu/)

Code maintainers:
* [Youngkyu Jeon](mailto:yjeon39@gatech.edu)
* [Penghao Xu](mailto:pxu64@gatech.edu)
* [Tejasvi Channagiri](mailto:tchannagri@usf.edu)
