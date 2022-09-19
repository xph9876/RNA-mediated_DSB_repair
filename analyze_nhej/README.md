# Analyze NHEJ

## Description

This folder contains the processing pipeline for the NHEJ (nonhomologous end joining) DNA repair. The main functionality such as filtering and plotting is implemented in Python3. The run*.[sh|ps1] scripts in the main directory are for reproducing the computations with datasets used in the [publication](#citation). Note *.sh scripts are intended to be run with Unis bash, while *.ps1 scripts are intended to be run with Windows PowerShell (although there is no semantic difference is how the Python scripts are called from either file type).

## Prerequisties

TODO: mention the package names/versions in the requirement.txt

## Citation

TODO

## Pipeline Stages

### 0_generate_scripts

Meta-scripts used to generate the \*.[sh|ps1] scripts for running the data processing pipeline on the data sets for the [publication](#citation). These are not required if using any of the pipline stages individually. All generate_\*.py scripts are run without arguments. 0_generate_scripts/run_all.[sh|ps1] runs them all.

### 1_process_nhej

Initial processing of raw SAM files which have been aligned to a reference sequence.

#### filter_nhej.py

Takes as input a SAM file and filters only alignments that are defined as being produced via NHEJ (see methods of the [publication](#citation) for the exact definition).

Arguments:

* --ref_seq_file: FASTA file with the reference sequence used to align the INPUT SAM file. Should contain a single nucleotide sequence in FASTA format.

* --sam_file: Aligned SAM input file. Must be created with [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (specific flags from Bowtie2 are used). Every read must be aligned with exactly the same reference sequence.

* --output: Output file. The result is a tab-separated file with columns: Sequence: nucleotide sequence of the read; CIGAR: CIGAR string for the alignment produced by Bowtie2 (see Bowtie2 documentation); Count: number of reads with the same sequence; Num_Subst: number of substitutions (AKA mismatches) in the alignment. Note, if multiple reads with identical nucleotide sequences and different CIGARs occur in SAM_FILE, then they will be a single row in the output representing the sequence with the Count column as the number of such reads and the CIGAR column arbitrarily selected from one of the CIGARs.

* --min_length: Minimum length of a read to pass filtering. Reads shorter than this are discarded.

* --dsb_pos: Position on reference sequence immediately upstream of DSB site (i.e., the DSB is between 1-based position DSB_POS and DSB_POS + 1.).

* --quiet: If present, do not output extra log message showing how many reads were discarded due to each filter criteria.

#### combine_repeats.py

Combines tables produced with [filter_nhej.py](#filternhejpy) representing biological repeats into a single file. The input must be tab-separate file produced with the  script and the ouput.

Arguments:

* --input: List of tab-separated files output from script [filter_nhej.py](#filternhejpy). Must have columns: Sequence, CIGAR, Count, Num_Subst. The files names must be of the form &lt;name&gt;&lowbar;&ast;, where &lt;name&gt; is the name of the library. 

* --output: Output tab-separated file. The output Count columns will be of the form Count_&lt;name&gt;, where &lt;name&gt; is the name of the library (see INPUT).

* --quiet: If present, do not output extra log messages.

### 2_get_window_data

Scripts for further processing the raw NHEJ data in order to extract window around the DSB, convert raw counts to frequencies, merge experiments which have been sequenced twice, and precompute other data used for downstream visulizations.

### get_window.py

Extract DSB-sequence for the NHEJ variation-distance graphs while discarding sequences that do not have proper anchors flanking the window.

Arguments:

* --input: TSV file output of combine_repeat.py. The columns must be: "Sequence", "CIGAR", "Count_&lt;X1&gt;", "Count_&lt;X2&gt;", etc. All the columns after "CIGAR" should be the counts for each repeat where &lt;Xi&gt; denotes the name of the library.

* --ref_seq_file: Reference sequence FASTA. Should contain a single nucleotide sequence in FASTA format.

* --output: Output directory. The output file will be TSV format and named like "window_*.py", with the suffix after "window_" indicating the characteristics of the file (e.g., with/without subsitutions, mean/count, or filtered). The library metadata is also output in "data_info.tsv".

* --dsb_pos: Position on reference sequence immediately upstream of DSB site (i.e., the DSB is between 1-based position DSB_POS and DSB_POS + 1.).

* --window_size: Size of window around DSB site to extract. The nucleotides at the positions {DSB_POS - WINDOW_SIZE + 1, ..., DSB_POS + WINDOW_SIZE} are extracted. The actual number of nucleotides extracted may vary depending on how many insertions/deletion the alignment of the sequence has.

* --anchor_size: Size of anchor on left/right of the window to check for mismatches. Reads with more than the allowed number of mismatches on the left/right anchor will be discarded. The mismatches on the left/right are counted separately.

* --anchor_mismatches: Maximum number of mismatches allowed on the left/right anchor sequences. Reads with more than the allowed number of mismatches on the left/right anchor will be discarded. The mismatches on the left/right are counted separately.

* --subst_type: Whether to keep or ignore alignment substitutions. If ignoring aligment substitutions, the corresponding nucleotide on the read is replaced with the reference sequenec nucleotide.

* --construct: Name of construct (e.g., Sense, Antisense, etc.,). Used only in the library metadata.

* --control_type: Whether this data is a control experiment and what type (e.g., noDSB, etc.,). Used only in the library metadata.

* --dsb_type: Whether this data is 1 DSB, 2 DSB, or 2 DSB antisense. Used only in the library metadata.

* --guide_rna: Type of guide RNA used in this experiment (e.g., sgRNA A, sgRNA B, etc.,). Used only in the library metadata.

* --strand: Strand of the reads in this library (e.g., R1 for forward, R2 for reverse). Used only in the library metadata.

* --cell_line: Cell in this library (e.g., WT or KO).

* --version: Version indicator since some libraries were sequenced twice (e.g., old, new, merged, none).

#### get_freq.py

Convert the raw read counts in the input data into frequencies using the input total reads. Outputs 3 files: (1) windows_freq.tsv: contains the all the sequences with the counts converted to frequencies. (2) windows_freq_filter.tsv: the previous file with the sequences removed whose frequency is <= FREQ_MIN in any of the repeats. (3) windows_freq_filter_mean.tsv: contains the means of the frequencies in the previous file (over all repeats).

Arguments:

* --input: Directory with output tables from get_window.py or get_merged.py.

* --total_reads: Total reads for each file. Must be the same number of arguments as the number of Count columns in the tables in INPUT.

* --output: Output directory.

* --subst_type: Whether to process the files with/without substitutions.

* --freq_min: Minimum frequency for output in windows_freq_filter_mean.tsv. Sequences with frequences <= FREQ_MIN are discarded.

#### get_merged.py (NEXT)

### 3_get_graph_data

#### get_graph_data.py

### 4_get_histogram_data

### 5_get_histogram_data

### 6_plot_histogram

### 7_get_pptx

* filte
* **run_01_process_nhej** - Filter the raw SAM files to retain the NHEJ patterns.
* run_02_combine_repeats - Combine repeat experiments and retain only common sequences.
3. run_03_make_main_data_withSubst, run_03_make_main_data_withoutSubst - Make processed/filtered data files for the variation-position graphs and 3D histograms. The "withSubst" script keeps the substitutions, and "withoutSubst" script ignores them.
4. run_04_make_graph_data_withSubst, run_04_make_graph_data_withoutSubst - Make processed data files for generating the variation-position graphs. The "withSubst"/"withoutSubst" script uses the respective output of the previous stage.
5. run_05_make_main_data_combined_withSubst, run_05_make_main_data_combined_withoutSubst - Makes the combined data files for comparing 2 experiments.
6. run_05_make_graph_data_combined_withSubst, run_05_make_graph_data_combined_withoutSubst - Makes the combined data files for comparing 2 experiments.
7. run_07_make_historam_3d - Plots the 3D variation-position histograms.
8. run_08_common_layout - Create the common layouts for the plotting the variation-position graphs. This insures that comparable experiments have vertices with the same variation placed in the same position.
9. run_09_plot_graph - Plots the variation-position graphs.
10. run_10_make_pptx_graph - Arranges the variation-position graphs in a grid in PPTX files.
11. run_11_make_pptx_histograms - Arranges the variation-position histograms in a grid in PPTX files.
12. run_12_plot_graph_main - Plots some of the variation-position graphs (slightly modified to reduce white space).

Running
-------
To run all stages of the pipeline use the "run_all.bat" file (Windows cmd) or "run_all.sh" file (Unix bash). To run one of the individual stages above use the corresponding "run*.bat" file (Windows cmd) or "run*.sh" file (Unix bash). To run individual Python files in the pipeline use Python 3 and refer to the .bat or .sh scripts for the usage, or run the files without any arguments to print a help message. These scripts have been tested on Windows 11, Python 3.10.4. Note, the steps in the pipeline must be run in proper order as later stages depend on output from previous stages.