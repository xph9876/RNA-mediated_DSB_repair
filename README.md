# RNA-mediated DSB repair

## Description

Contains the code for the computational analysis presented in the [publication](#citation). Please see the Methods section of this publication for a more detailed description of the code.

## Citation

Jeon Y. *et al*. RNA-mediated double-strand break repair in human cells. (2022). bioRxiv. https://doi.org/10.1101/2022.11.01.514688.

## Files/Directories

* `trimming`: Python 3 scripts for the [trimming](#trimming) stage.
* `flipped_intron`: Python 3 scripts for the [frequency of flipped intron](#frequency-of-flipped-intron) stage.
* `intron`: Python 3 scripts for the [frequency of intron](#frequency-of-intron) stage.
* `intron_yeast`: Python 3 scripts for the [frequency of intron](#frequency-of-intron) stage (for yeast data).
* `RTDR`: Python3 scripts for the [frequency of R-TDR](#frequency-of-rtdr) stage.
* `RNA_seq`: Python 3 scripts for the [categorization of RNA-seq](#categorization-of-rna-seq) stage.
* `MMEJ`: Python 3 scripts for the [MMEJ pipeline](#mmej-pipeline).
* `refseq`: Reference sequences for [MMEJ pipeline](#mmej-pipeline).
* `analyze_MMEJ.sh`: Main script for running the [MMEJ pipeline](#mmej-pipeline) for sense libraries.
* `analyze_MMEJ_antisense.sh`: Main script for running the [MMEJ pipeline](#mmej-pipeline) for antisense libraries.
* `analyze_MMEJ_OX.sh`: Main script for running the [MMEJ pipeline](#mmej-pipeline) for RNase H1 overexpression/control libraries.
* `libinfo.tsv`: Library metadata for the [MMEJ pipeline](#mmej-pipeline) for sense libraries.
* `libinfo_antisense.tsv`: Library metadata for the [MMEJ pipeline](#mmej-pipeline) for antisense libraries.
* `libinfo_OX.tsv`: Library metadata for the [MMEJ pipeline](#mmej-pipeline) for RNase H1 overexpression/control libraries.
* `test_ratios`: Python 3 scripts for the [test ratios](#test-ratios) stage.
* `NHEJ`: Python 3 scripts for the [NHEJ pipeline](#nhej-pipeline). See `NHEJ/README.md` for more details.
* `demo`: Demonstrations of the various analyses and pipelines with example data.
* `unknown`: Python 3 scripts for analyzing the unclassified sequences in 1-DSB experiments.

## Dependencies

For stages [trimming](#trimming), [frequency of intron](#frequency-of-intron), [frequency of flipped intron](#frequency-of-flipped-intron),[frequency of RTDR](#frequency-of-rtdr), [categorization of RNA-seq](#categorization-of-rna-seq), [test ratios](#test-ratios), [MMEJ pipeline](#mmej-pipeline), please use `MMEJ/mmej_environment.yml` to create the conda environment.

```bash
# After conda installation, run in repository folder
conda env create -f MMEJ/mmej_environment.yml
# Before each use, activate conda environment
conda activate microhomology
```

For stage [NHEJ pipeline](#nhej-pipeline) and [unknown analysis](#unknown-analysis), please see `NHEJ/README.md`. The unknown analysis additionally requires the Python package BioPython 1.80.

## Installation
After installing all dependencies, just clone the repository to install.
```bash
git clone https://github.com/xph9876/RNA-mediated_DSB_repair.git
```

## Stages

### Trimming

Trimming must be performed on all raw reads, DNA-seq and RNA-seq, before the other stages. First, the reads must first be processed with [cutadapt](https://cutadapt.readthedocs.io/en/stable/) and then with the trimming scripts `trimming/trim_F_tag.py` (forward strand reads) and `trimming/trim_R_tag.py` (reverse strand reads). See the Methods of the [publication](#citation) for more details.

### Frequency of intron

The script `intron/freq_intron.py` calculates the frequency of reads with the intron. The input should be the trimmed FASTQ files from [trimming](#trimming). 

For the yeast data a different method of detection is used, which are implemented in `intron_yeast/freq_intron_with_seq_F.py` and `intron_yeast/freq_intron_with_seq_R.py`. These scripts should be run as in `intron_yeast/run.sh`.

### Frequency of flipped intron

The scripts `flipped_intron/flipped_intron*.py` calculate the frequency of the reads with the flipped intron in the `sense` and `antisense` constructs, for forward (`F`) and reverse (`R`) reads. The input should be the trimmed FASTQ files from [trimming](#trimming).

### Frequency of RTDR

The scripts `RTDR/RTDR_*.py` calculate the frequency of the reads with RNA-templated DNA repair (RTDR) for forward (`F`) and reverse (`R`) strand reads. The input should be the trimmed FASTQ files from [trimming](#trimming).

### Categorization of RNA-seq

The scripts `RNA_seq/*.py` categorize RNA-seq reads which have been aligned with the [hisat2](http://daehwankimlab.github.io/hisat2/) software. See the Methods of the [publication](#citation) for more details. The input should be the aligned SAM files. There are separate scripts for each of the different constructs, 5'-SplicingΔ, Antisense, BranchΔ, Sense, and pCMVΔ, as indicated by the file names.

### MMEJ pipeline

Pipeline for the microhomology-mediated end joining (MMEJ) analysis. Computes microhomology pairs on the references sequences, calculates MMEJ frequencies in trimmed reads, and generates figures.

1. Place the FASTQ files into the appropriate directory with the appropriate name (see comments in `analyze_MMEJ.sh`, `analyze_MMEJ_antisense.sh`, and `analyze_MMEJ_OX.sh`).
2. Set the variables in `analyze_MMEJ.sh`, `analyze_MMEJ_antisense.sh`, and `analyze_MMEJ_OX.sh` (e.g., the input and output directories).
3. Run `analyze_MMEJ.sh`, `analyze_MMEJ_antisense.sh`, and `analyze_MMEJ_OX.sh`.

### Test ratios

Performs Mann-Whitney U tests to compare the ratio Sense/BranchΔ (or other constructs) of repair in wild-type and RNase H2A (or other) cells. The input frequencies should be the TSV files output from the [MMEJ pipeline](#mmej-pipeline). The output of the [frequency of intron](#frequency-of-intron), [frequency of flipped intron](#frequency-of-flipped-intron), and [frequency of RTDR](#frequency-of-rtdr) stages can also be used as long they are formatted in a TSV file with the proper columns: `Sample`, `Read`, `Cell_line`, `Genotype`, `Breaks`, `Count`, and, `Frequency`. See source code of `MMEJ/calc_mh_freq.py` or the output of the [MMEJ pipeline](#mmej-pipeline) for more details.

The input frequencies used in the publication are included in the `test_ratios/input` subdirectory. To run the analysis use the `test_ratios/run.sh` script.

### NHEJ pipeline

The non-homologous end joining (NHEJ) analysis pipeline. Extracts DSB-sequence windows from aligned reads, and plots variation-distance graphs and variation-position histograms. See `NHEJ/README.md` for more details.

### Unknown analysis

An analysis pipline for examining the unclassified reads of the 1-DSB experiments. To run the analysis use the `unknown/run_all.sh` or `unknown/run_all.ps1` script. Note that the NHEJ pipeline must already have been run (particularly, the tables in `NHEJ/1_data_filter_nhej` must have been generated). Some of the input to this pipline in `unknown/input` has been derived from other pipelines such as `MMEJ`. The output tables will be written to `unknown/output`.

## Contact

Webpages:
* [Francesca Storici Lab](https://storicilab.gatech.edu/)
* [Natasha Jonoska Group](https://knot.math.usf.edu/)

Code maintainers:
* [Youngkyu Jeon](mailto:yjeon39@gatech.edu)
* [Penghao Xu](mailto:pxu64@gatech.edu)
* [Tejasvi Channagiri](mailto:tchannagri@usf.edu)
