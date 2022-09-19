# Analyze NHEJ

## Description

TODO: italicize all file names.

This folder contains the processing pipeline for the NHEJ (nonhomologous end joining) DNA repair. The main functionality such as filtering and plotting is implemented in Python3. The *run&ast;.sh* and *run&ast;.ps* scripts in the main directory are for reproducing the computations with datasets used in the [publication](#citation). Note *.sh scripts are intended to be run with Unix bash, while *.ps1 scripts are intended to be run with Windows PowerShell (although there is no semantic difference is how the Python scripts are called from either file type).

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

#### get_merged.py

Merge together library files from samples that have beem sequenced twice. Operates on the output of get_window.py.

Arguments:

* --input: Input directories with "windows_*.tsv" output of get_window.py. There must be the same number of columns in all input data sets. The "count" columns must be in the same order they should be merged in.

* --output: Output directory. Files will be named in the same manner as get_window.py.

* --subst_type: Whether to operate on file with or with substitutions.

* --version: Version ID of the merged output (e.g., "merged"). Used only for the metadata output data_info.tsv.

* --new_library_names: Names of the new merged libraries. Must be the same number of arguments as the number of count columns in the data.

#### get_freq_comparison.py

Combine two individual experiment directories to make a comparison experiment directory for comparison graphs. The experiments must be compatible: have all the same attributes except for the constructs which must be different. Operates of the output of get_freq.py. The tables from the two libraries are merged by joining on the alignment columns.

Arguments:

* --input: Directory of individual experiment frequency data produced with get_freq.py.

* --output: Output directory. Output files parallel the structure of get_freq.py output, except they have two frequencies columns for the two compared libraries.

* --subst_type: Whether to operate on files with/without substitutions.

### 3_get_graph_data

#### *get_graph_data.py*

Precompute the necessary data for the variation-distance graphs. Uses as input the output of the [2_get_window_data](#2getwindowdata) stage and outputs the following files:

* sequence_data_*.tsv: Table of the vertex data for the variation distance graphs. Each row represents a single vertex.

* edge_data_*.tsv: Table of edge data for variation-distance graphs. Each row represents a single edge.

* distance_matrix_*.tsv: Table of pairwise Levenshtein distances of the vertices.

* graph_stats_*.tsv: Summary statistics of the graph (e.g., number of vertices).

Arguments

* --input: Directory with output from [2_get_window_data](#2getwindowdata) stage.

* --output: Output directory.

* --subst_type: Whether to process the files with/without substitutions.

### 4_get_histogram_data

#### *get_histogram_data.py*

Precomputes the necessary data for the variation-position histgrams. Uses as input the output of the [3_get_graph_data](#3getgraphdata) stage and outputs the following files:

* *variation_&ast;.tsv*: Contains data on individual variations of each sequence window in the *sequence_&ast;.tsv* files. Each row represents a single variation (insertion, deletion, or substitution) in a sequence alignment.

* *variation_grouped_&ast;.tsv*: The same as data as the corresponding *variation_&ast;.tsv* grouped by (1) position of variation on the reference sequence, (2) total number of variations in the parent sequence, and (3) the type of variation (insertion, deletion, or substitution). The frequencies have have the aggregate of all individual variations with the same grouping characteristics.

Arguments:

* --input: Directory with output from [3_get_graph_data](#3getgraphdata).

* --output: Output directory.

* --subst_type: Whether to process the files with/without substitutions.

### 5_plot_graph

Code for laying out and plotting variation-distance graphs. Uses output from [3_get_graph_data](#3getgraphdata).

#### *get_precomputed_layout.py*

Precomputes the layout for groups of experiments by taking the union of all sequences in all input experiments and laying them out. This allows experiments with the same reference sequence to have the same coordinates for the same sequence. All experiments should have the same window reference sequence and be individual (not comparison) experiments. There are multiple possible layouts:

* radial: Arranges vetices in a radial grid around the reference sequence. Insertion (deletion) vertices are placed above (below) the reference sequence. Heuristics are used to improve aesthesics of vertex placement, though this may not easily egeneralize to new data.

* mds: Uses multi-dimensional scaling (MDS) to embed the pairwise Levenshtein distance matrix into 2D. 

* kamada: Uses the [*NetworkX*](http://networkx.org) package's implementation of the Kamada-Kawaii algorithm to layout the graph. Reference: T. Kamada and S. Kawai. An algorithm for drawing general undirected graphs. *Inform. Process. Lett.*, 31:7–15, 1989.

* universal: Uses a deterministic layout descibed in the Methods of the [publicaton](#citation).

Arguments:

* --input: List of data directories created with [*get_graph_data.py*](#getgraphdatapy).

* --output: Output directory. Two files are created within this directory paralleling the output of [*get_graph_data.py*](#getgraphdatapy): *sequence_data_&ast;.py* and *edge_data_&ast;.py*.

* --reverse_complement: Whether to reverse complement the sequence in the corresponding input. If present, the number of values must be the same as the number of input directories. "1" mean reverse complement the sequences and "0" means do not. Used for making a common layout for data sets that have window reference sequences that are the reverse complements of each other.

* --subst_type: Whether to process files with/without substitutions.

* --layout: One of "radial", "mds", "kamada", or "universal". See above for a description of each layout.

#### *plot_graph.py*

Does layout and plotting of the processed graph data from [3_get_graph_data](#3getgraphdata).

* --input: Directory with the data files produced with [3_get_graph_data](#3getgraphdata).

* --output: Output directory. Optional, if not given no output will be written.

* --ext
    '--ext',
    choices = ['png', 'html'],
    default = 'png',
    help = (
      'Which types of file to generate from the Plotly library:' +
      ' static PNG or interactive HTML.'
    ),
  )
  parser.add_argument(
    '--title',
    action = 'store_true',
    help = (
      'If present, adds a title to the plot showing the type of'
      ' and the name of the data set.'
    )
  )
  parser.add_argument(
    '--layout',
    choices = ['kamada', 'radial', 'mds', 'universal', 'fractal'],
    default = 'radial',
    help = 'The algorithm to use for laying out the graph.'
  )
  parser.add_argument(
    '--universal_layout_y_axis_x_pos',
    type = float,
    help = (
      'If present, shows a y-axis at the given x position' +
      ' on the universal layout showing the distances to the reference.'
    )
  )
  parser.add_argument(
    '--universal_layout_x_axis_deletion_y_pos',
    type = float,
    help = (
      'If present, shows a x-axis for deletions at the given y position' +
      ' on the universal layout showing the midpoints of the deleted ranges.'
    )
  )
  parser.add_argument(
    '--universal_layout_x_axis_insertion_y_pos',
    type = float,
    help = (
      'If present, shows a x-axis for insertions at the given y position' +
      ' on the universal layout showing the first nucleotide of inserted sequences.'
    )
  )
  parser.add_argument(
    '--universal_layout_y_axis_y_range',
    nargs = '+',
    type = float,
    help = (
      'If showing an y-axis for the universal layout,' +
      ' the min and max y-position of the line.'
    )
  )
  parser.add_argument(
    '--universal_layout_x_axis_x_range',
    nargs = '+',
    type = float,
    help = (
      'If showing an x-axis for the universal layout,' +
      ' the min and max x-position of the line.'
    )
  )
  parser.add_argument(
    '--universal_layout_y_axis_deletion_max_tick',
    type = int,
    help = (
      'If showing an y-axis for the universal layout,' +
        ' the max tick value for the deletion side.'
    )
  )
  parser.add_argument(
    '--universal_layout_y_axis_insertion_max_tick',
    type = int,
    help = (
      'If showing an y-axis for the universal layout,' +
        ' the max tick value for the insertion side.'
    )
  )
  parser.add_argument(
    '--subst_type',
    choices = library_constants.SUBST_TYPES,
    help = 'Whether to plot data with or without substitutions.',
    default = library_constants.SUBST_WITHOUT,
  )
  parser.add_argument(
    '--node_max_freq',
    type = float,
    help = (
      'Max frequency to determine node size.' +
      'Higher frequencies are clipped to this value.'
    ),
    default = library_constants.GRAPH_NODE_SIZE_MAX_FREQ,
  )
  parser.add_argument(
    '--node_min_freq',
    type = float,
    help = (
      'Min frequency to determine node size.' +
      'Lower frequencies are clipped to this value.'
    ),
    default = library_constants.GRAPH_NODE_SIZE_MIN_FREQ,
  )
  parser.add_argument(
    '--node_max_px',
    type = float,
    help = 'Largest node size as determined by the frequency.',
    default = library_constants.GRAPH_NODE_SIZE_MAX_PX,
  )
  parser.add_argument(
    '--node_min_px',
    type = float,
    help = 'Smallest node size as determined by the frequency.',
    default = library_constants.GRAPH_NODE_SIZE_MIN_PX,
  )
  parser.add_argument(
    '--node_outline_scale',
    type = float,
    default = library_constants.GRAPH_NODE_OUTLINE_WIDTH_SCALE,
    help = (
      'How much to scale the node outline width (thickness).' +
      ' Values > 1 increase the width; values < 1 decrease the width.'
    ),
  )
  parser.add_argument(
    '--variation_types',
    nargs = '+',
    help = (
      'The variation types that should be included in the graph.'
      ' This should be a list of the types:'
      ' "insertion", "deletion", "substitution", "none".' +
      ' Default value: "insertion", "deletion", "none".' +
      ' "none" means the reference sequence.',
    ),
  )
  parser.add_argument(
    '--edge_scale',
    type = float,
    help = (
      'How much to scale the edges width (thickness).' +
      ' Values > 1 increase the width; values < 1 decrease the width.'
    )
  )
  parser.add_argument(
    '--width_px',
    type = int,
    default = library_constants.GRAPH_WIDTH_PX,
    help = 'The width of the plot in pixels.'
  )
  parser.add_argument(
    '--height_px',
    type = int,
    default = library_constants.GRAPH_HEIGHT_PX,
    help = 'The height of the plot in pixels.'
  )
  parser.add_argument(
    '--line_width_scale',
    type = float,
    default = library_constants.GRAPH_LINE_WIDTH_SCALE,
    help = (
      'How much to scale the line widths (aka thickness).' +
      ' Values > 1 increase the width; values < 1 decrease the width.'
    ),
  )
  parser.add_argument(
    '--font_size_scale',
    type = float,
    default = library_constants.GRAPH_FONT_SIZE_SCALE,
    help = (
      'How much to scale the font size.' +
      ' Values > 1 increase the font size; values < 1 decrease it.'
    ),
  )
  parser.add_argument(
    '--precomputed_layout_dir',
    type = common_utils.check_dir,
    default = None,
    help = (
      'If present, gives the directory where the precomputed layouts are.' +
      ' If not present the layout is computed newly.'
    )
  )
  parser.add_argument(
    '--reverse_complement',
    action = 'store_true',
    help = (
      'If present, uses the reverse complement of sequences when determining the'
      ' node positions, and displaying labels and hover text.' +
      ' This affects the precomputed layout, universal layout, and fractal layout.'
    )
  )
  parser.add_argument(
    '--crop_x',
    nargs = '+',
    type = float,
    help = (
      'Range of the horizontal dimension to crop.' +
      ' Specified with normalized coords in range [0, 1].'
    ),
  )
  parser.add_argument(
    '--crop_y',
    nargs = '+',
    type = float,
    help = (
      'Range of the vertical dimension to crop.' +
      ' Specified in normalized coords in range [0, 1].'
    ),
  )
  parser.add_argument(
    '--range_x',
    type = float,
    nargs = '*',
    help = (
      'Range of x-axis for plotting.'
      'If not specified chosen automatically to either show all nodes or a preset value'
      ' for the layout.'
    ),
  )
  parser.add_argument(
    '--range_y',
    type = float,
    nargs = '*',
    help = (
      'Range of y-axis for plotting.'
      'If not specified chosen automatically to either show all nodes or a preset value'
      ' for the layout.'
    ),
  )
  parser.add_argument(
    '--legend',
    action = 'store_true',
    help = 'Whether to show a legend on the figure.'
  )
  parser.add_argument(
    '--legend_color_bar_scale',
    type = float,
    default = library_constants.GRAPH_LEGEND_COLORBAR_SCALE,
    help = 'How much to scale the legend color bar (for freq ratio coloring).'
  )
  parser.add_argument(
    '--separate_components',
    action = 'store_true',
    help = 'If present separate the connected components of the graph.'
  )

### 6_plot_histogram

### 7_get_pptx

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