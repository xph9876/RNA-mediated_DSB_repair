# RNA-mediated DSB repair

## Description

Contains the code for the computational analysis presented in the [publication](#citation). Please mee the Methods section of this publication for a more detailed description of the code.

## Citation

Jeon Y. *et al*. RNA-mediated double-strand break repair in human cells. (2022). *In submission*.

## Files/Directories

* `trimming`: Python3 scripts for the [trimming](#trimming) stage.

* `flipped_intron`: Python3 scripts for the [frequency of flipped intron](#frequency-of-flipped-intron) stage.

* `intron`: Python3 scripts for the [frequency of intron](#frequency-of-intron) stage.

* `RTDR`: Python3 scripts for the [frequency of RTDR](#frequency-of-rtdr) stage.

* `RNA_seq`: Python3 scripts for the [categorization of RNA-seq](#categorization-of-rna-seq) stage.

* `MMEJ`: Python3 scripts for the [MMEJ pipeline](#mmej-pipeline).

* `permutation_test`: Python3 scripts for the [permutation test](#permutation-test) stage.

* `NHEJ`: Python3 scripts for the [NHEJ pipeline](#nhej-pipeline). See `NHEJ/README.md` for more details.

* `refseq`: Reference sequences for [MMEJ pipeline](#mmej-pipeline).

* `analyze_MMEJ.sh`: Main script for running the [MMEJ pipeline](#mmej-pipeline) for sense libraries.

* `analyze_MMEJ_antisense.sh`: Main script for running the [MMEJ pipeline](#mmej-pipeline) for antisense libraries.

* `libinfo.tsv`: Library metadata for the [MMEJ pipeline](#mmej-pipeline) for sense libraries.

* `libinfo_antisense.tsv`: Library metadata for the [MMEJ pipeline](#mmej-pipeline) for antisense libraries.

## Dependencies

The  [trimming](#trimming), [frequency of intron](#frequency-of-intron), [frequency of flipped intron](#frequency-of-flipped-intron), [frequency of RTDR](#frequency-of-rtdr), and [categorization of RNA-seq](#categorization-of-rna-seq) stages do not have any dependencies other than Python 3.8.13. The other stages' dependencies are below.

[Permutation test](#permutation-test):

* Python 3.8.13
* Xlsxwriter 3.0.3

[MMEJ pipeline](#mmej-pipeline):

* Python 3.8.13
* Numpy 1.22.3
* Matplotlib 3.4.3
* Pandas 1.4.2
* Seaborn 0.11.2
* Xlsxwriter 3.0.3
* Scipy 1.7.3
* Statannotations 0.4.3

[NHEJ pipeline](#nhej-pipeline): please see `NHEJ/README.md`.

## Stages

### Trimming

Trimming must be performed on all raw reads, DNA-seq and RNA-seq, before the other stages. First, the reads must first be processed with [cutadapt](https://cutadapt.readthedocs.io/en/stable/) and then with the trimming scripts `trimming/trim_F_tag.py` (forward strand reads) and `trimming/trim_R_tag.py` (reverse strand reads). See Methods of the [publication](#citation) for more details.

### Frequency of intron

The script `intron/freq_intron.py` calculates the frequency of reads with the intron. The input should be the trimmed FASTQ files from [trimming](#trimming).

### Frequency of flipped intron

The scripts `flipped_intron/flipped_intron*.py` calculate the frequency of the reads with the flipped intron in the `sense` and `antisense` constructs, for forward (`F`) and reverse (`R`) reads. The input should be the trimmed FASTQ files from [trimming](#trimming).

### Frequency of RTDR

The script `RTDR/RTDR_R.py` calculates the frequency of the reads with RNA-templated DNA repair (RTDR) for reverse strand reads. The input should be the trimmed FASTQ files from [trimming](#trimming).

### Categorization of RNA-seq

The scripts `RNA_seq/*.py` categorize RNA-seq reads which have been aligned with the [hisat2](http://daehwankimlab.github.io/hisat2/) software. The input should be the aligned SAM files. There are separate scripts for each of the different constructs, 5'-SplicingΔ, Antisense, BranchΔ, Sense, and pCMVΔ, as indicated by the file names. Note, the raw RNA-seq reads do not require [trimming](#trimming).

### MMEJ pipeline

Pipeline for the microhomology-mediated end joining (MMEJ) analysis. Computes microhomology pairs on the references sequences, calculates MMEJ frequencies in trimmed reads, and generates figures.

1) Place the FASTQ files into the appropriate directory with the appropriate name (see comments in `analyze_MMEJ.sh` and `analyze_MMEJ_antisense.sh`).

2) Run `analyze_MMEJ.sh` and `analyze_MMEJ_antisense.sh`.

### Permutation test

Performs permutations tests to compare the ratio BranchΔ/Sense of repair in wild-type vs. KO cells. The input frequencies should be the TSV files output from the [MMEJ pipeline](#mmej-pipeline). The output of the [frequency of intron](#frequency-of-intron), [frequency of flipped intron](#frequency-of-flipped-intron), and [frequency of RTDR](#frequency-of-rtdr) stages can also be used as long they are formatted in a TSV file with the proper columns: `Sample`, `Read`, `Cell_line`, `Genotype`, `Breaks`, `Count`, and, `Frequency`. See source code of `MMEJ/calc_mh_freq.py` or the output of the [MMEJ pipeline](#mmej-pipeline) for more details.

### NHEJ pipeline

The non-homologous end joining (NHEJ) analysis pipline. Extracts DSB-sequence windows from aligned reads, and plots variation-distance graphs and variation-position histograms.

1) Align the [trimmed](#trimming) FASTQ files with the appropriate reference sequence using Bowtie2 (see `NHEJ/README.md` for examples).

2) Place the SAM file output of alignment in the appropriate directory with appropriate name (see `NHEJ/README.md`).

3) Run either `NHEJ/run_all.ps1` (Windows) or `NHEJ/run_all.sh` (Unix). Note, the working directory of the terminal must be `NHEJ` for this to work correctly.

See `NHEJ/README.md` for more details.

## Contact

Webpages:
* [Francesca Storici Lab](https://storicilab.gatech.edu/)
* [Natasha Jonoska Group](https://knot.math.usf.edu/)

Code maintainers:
* [Youngkyu Jeon](mailto:yjeon39@gatech.edu)
* [Penghao Xu](mailto:pxu64@gatech.edu)
* [Tejasvi Channagiri](mailto:tchannagri@usf.edu)
