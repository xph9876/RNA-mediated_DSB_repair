# NHEJ

## Description

This folder contains the processing pipeline for the NHEJ (non-homologous end joining) DNA repair. The main functionality such as filtering and plotting is implemented in Python3. The `run*.sh` and `run*.ps` scripts in the main directory are for reproducing the computations with data sets used in the [publication](#citation). Note `*.sh` scripts are intended to be run with Unix bash, while `*.ps1` scripts are intended to be run with Windows PowerShell (although there is no difference is how the Python scripts are called from either file type).

## Dependencies

Tested OS: Windows 10 Home.

* Software:
    * Python 3.10.6
    * Bowtie2 2.4.1
* Python packages:
    * Kaleido 0.2.1
    * Matplotlib 3.5.3
    * NetworkX 2.8.6
    * Numpy 1.23.2
    * Pandas 1.4.3
    * Pillow 9.2.0
    * Plotly 5.10.0
    * Psutil 5.9.1
    * Python-Levenshtein 0.12.2
    * Python-pptx 0.6.21
    * Requests 2.28.1
    * Scikit-Learn 1.1.2
    * Scipy 1.9.1
    * XlsxWriter 3.0.3

The full output of ```pip freeze``` is given in `python_packages.txt`. For PNG image output from the [Plotly](https://plotly.com/) library, a backend known as [Orca](https://github.com/plotly/orca) had to be used due to to problems with the [Kaleido](https://github.com/plotly/Kaleido) backend (at least on Windows 10). Please see [here](https://plotly.com/python/static-image-export/) for details on installation/setup of Orca.

## Citation

Jeon Y. *et al.* RNA-mediated double-strand break repair in human cells. (2022). *In submission*.

## Reproducing Analyses

To reproduce the NHEJ analyses of the [publication](#citation), the scripts with names of the form `*.sh` and `*.ps1` (such as `run_01_process_nhej.ps1` and `run_01_process_nhej.sh`) must be run in the order indicated by their numbering. These scripts also serve as usage examples for the Python scripts. The scripts `run_all.sh` and `run_all.ps1` run all stages. All scripts must be run with the `NHEJ` directory as the current working directory of the terminal.

## Pipeline Stages

### Alignment

The raw input for the NHEJ pipeline are the trimmed FASTQ files from the trimming stage (see parent directory). The trimmed FASTQ files must be alignment against the appropriate reference sequence in the directory `ref_seq`. Bowtie2 should be used with the default settings:

```
  bowtie2-build <ref_name>.fa <ref_name>
  bowtie2 -x <ref_name> <fastq_name>.fa -S <sam_name>.sam
```

Where `<ref_name>` is the name of reference sequence, `<fastq_name>` is the name of the trimmed FASTQ file, and `<sam_name>` is the name of the output SAM file.

For the pre-made `*.sh` and `*.ps1` scripts to be run successfully, all SAM files produced as described above must be placed in the `data_0_sam` directory. Each SAM file must also be named in the format: `<library>_<cell_line>_<guide_rna>_<construct>.sam` where

* `<library>`: Name of the library (alphanumeric).
* `<cell_line>`: Cell line. `WT` for wild type and `KO` for knock-out.
* `<guide_rna>`: Guide RNA. Choices: `sgA`, `sgB`, `sgAB`, or `sgCD`.
* `<strand>`: Strand. `R1` for forward and `R2` for reverse.
* `<construct>`: Plasmid construct. `sense` for Sense, `branch` for BranchΔ, `cmv` for pCMVΔ, `antisense` for Antisense, `splicing` for 5'-SplicingΔ.

Example: `yjl89_WT_sgCD_R1_antisense.sam`.

Please see `bowtie2_arguments.tsv` for a list of the libraries used in this analysis, their corresponding reference-sequence files, and the proper SAM file names.

### `0_generate_scripts`

Meta-scripts used to generate the `*.sh` and `*.ps1` scripts for running the data processing pipeline on the data sets for the [publication](#citation). These are not required if using any of the pipeline stages individually. All `generate_*.py` scripts are run without arguments. 

### `1_process_nhej`

Initial processing of the raw SAM files obtained from aligning the trimmed reads to the proper reference sequence.

#### `filter_nhej.py`

Takes as input a SAM file and filters only alignments that are defined as being produced via NHEJ (see methods of the [publication](#citation) for the exact definition).

Arguments:

* `--ref_seq_file`: FASTA file with the reference sequence used to align the input `SAM_FILE`. Should contain a single nucleotide sequence in FASTA format.

* `--sam_file`: Aligned SAM input file. Must be created with [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (specific flags from Bowtie2 are used). Every read must be aligned with exactly the same reference sequence.

* `--output`: Output file. The result is a tab-separated file with columns: `Sequence`: nucleotide sequence of the read; `CIGAR`: CIGAR string for the alignment produced by Bowtie2 (see [Bowtie2 documentation](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)); `Count`: number of reads with the same sequence; `Num_Subst`: number of substitutions (AKA mismatches) in the alignment. If multiple reads with identical nucleotide sequences and different CIGARs occur in `SAM_FILE`, then they will be a single row in the output representing the sequence with the `Count` column as the number of such reads and the `CIGAR` column arbitrarily selected from one of the CIGARs.

* `--min_length`: Minimum length of a read to pass filtering. Reads shorter than this are discarded.

* `--dsb_pos`: Position on reference sequence immediately upstream of DSB site (i.e., the DSB is between 1-based position `DSB_POS` and `DSB_POS + 1`.).

* `--quiet`: If present, do not output extra log message showing how many reads were discarded due to each filtering criteria.

#### `combine_repeat.py`

Combines multiple tables of biological repeats into a single file. Uses as input the TSV output of [`filter_nhej.py`](#filternhejpy).

Arguments:

* `--input`: List of TSV files output from script [`filter_nhej.py`](#filternhejpy). Must have columns: `Sequence`, `CIGAR`, `Count`, `Num_Subst`. The files names must be of the form `<name>_*.tsv`, where `<name>` is the name of the library. 

* `--output`: Output TSV file. The output `Count` columns will be of the form `Count_<name>` where `<name>` is the name of the library (see `INPUT`).

* `--quiet`: If present, do not output extra log messages.

### `2_get_window_data`

Scripts for further processing the raw NHEJ data in order to extract windows around the DSB, convert raw counts to frequencies, merge experiments which have been sequenced twice, and pre-compute other data used for downstream visualizations.

### `get_window.py`

Extract DSB-sequence for the NHEJ variation-distance graphs while discarding sequences that do not have proper anchors flanking the window.

Arguments:

* `--input`: TSV file output of [`combine_repeat.py`](#combinerepeatpy). The columns must be: `Sequence`, `CIGAR`, `Count_<X1>`, `Count_<X2>`, etc. All the columns after `CIGAR` should be the counts for each repeat where `<Xi>` denotes the name of the library.

* `--ref_seq_file`: Reference sequence FASTA. Should contain a single nucleotide sequence in FASTA format.

* `--output`: Output directory. The output file will be TSV format and named like `window_*.py`, with the suffix after `window_` indicating the characteristics of the file (e.g., with/without substitutions, frequency/count, filtered, etc.,). The library metadata is output in `data_info.tsv`.

* `--dsb_pos`: Position on reference sequence immediately upstream of DSB site (i.e., the DSB is between 1-based position `DSB_POS` and `DSB_POS + 1`.).

* `--window_size`: Size of window around DSB site to extract. The nucleotides at the positions `{DSB_POS - WINDOW_SIZE + 1, ..., DSB_POS + WINDOW_SIZE}` are extracted. The actual number of nucleotides extracted may vary depending on how many insertions/deletion the alignment of the sequence has.

* `--anchor_size`: Size of anchor on left/right of the window to check for mismatches. Reads with more than the allowed number of mismatches on the left/right anchor will be discarded. The mismatches on the left/right are counted separately.

* `--anchor_mismatches`: Maximum number of mismatches allowed on the left/right anchor sequences. Reads with more than the allowed number of mismatches on the left/right anchor will be discarded. The mismatches on the left/right are counted separately.

* `--subst_type`: Whether to keep or ignore alignment substitutions. If ignoring alignment substitutions, the corresponding nucleotide on the read is replaced with the reference sequence nucleotide.

* `--construct`: Name of construct (e.g., Sense, Antisense, etc.,). Used only in the library metadata.

* `--control_type`: Whether this data is a control experiment and what type (e.g., noDSB, etc.,). Used only in the library metadata.

* `--dsb_type`: Whether this data is 1-DSB, 20DSB, or 2-DSB-antisense. Used only in the library metadata.

* `--guide_rna`: Type of guide RNA used in this experiment (e.g., sgRNA A, sgRNA B, etc.,). Used only in the library metadata.

* `--strand`: Strand of the reads in this library (e.g., `R1` for forward, `R2` for reverse). Used only in the library metadata.

* `--cell_line`: Cell line in this library (e.g., `WT` or `KO`).

* `--version`: Version indicator since some libraries were sequenced twice (e.g., `old`, `new`, `merged`, or `none`).

#### `get_freq.py`

Convert the raw read counts in the input data into frequencies using the input total reads. Outputs 3 files: (1) `windows_freq.tsv`: contains the all the sequences with the counts converted to frequencies. (2) `windows_freq_filter.tsv`: the previous file with the sequences removed whose frequency is <= FREQ_MIN in any of the repeats. (3) `windows_freq_filter_mean.tsv`: contains the means of the frequencies in the previous file (over all repeats).

Arguments:

* `--input`: Directory with output tables from `get_window.py` or `get_merged.py`.

* `--total_reads`: Total reads for each file. The number of arguments must be the same as the number of `count_*` columns in the tables in `INPUT`.

* `--output`: Output directory.

* `--subst_type`: Whether to process the files with/without substitutions.

* `--freq_min`: Minimum frequency for output in `windows_freq_filter_mean.tsv`. Sequences with frequencies <= `FREQ_MIN` are discarded.

#### `get_merged.py`

Merge together library files from samples that have been sequenced twice. Operates on the output of [`get_window.py`](#getwindowpy).

Arguments:

* `--input`: Input directories with `windows_*.tsv` output of [`get_window.py`](#getwindowpy). There must be the same number of columns in all input data sets. The `count_*` columns must be in the same order they should be merged in.

* `--output`: Output directory. Files will be named in the same manner as [`get_window.py`](#getwindowpy).

* `--subst_type`: Whether to operate on files with or with substitutions.

* `--version`: Version ID of the merged output (e.g., `merged`). Used only to output metadata in the file `data_info.tsv`.

* `--new_library_names`: Names of the new merged libraries. Must be the same number of arguments as the number of `count_*` columns in the data.

#### `get_freq_comparison.py`

Combine two individual experiment directories to make a comparison experiment directory for comparison graphs. The experiments must be compatible: have all the same attributes except for the constructs, which must be different. Operates of the output of [`get_freq.py`](#getfreqpy). The tables from the two libraries are merged by joining on the alignment columns.

Arguments:

* `--input`: Directory of individual experiment frequency data produced with [`get_freq.py`](#getfreqpy).

* `--output`: Output directory. Output files parallel the structure of output from [`get_freq.py`](#getfreqpy), except they have two frequency columns, one for each compared library.

* `--subst_type`: Whether to operate on files with/without substitutions.

### `3_get_graph_data`

#### `get_graph_data.py`

Pre-compute the necessary data for the variation-distance graphs. Uses as input the output of the [`2_get_window_data`](#2getwindowdata) stage and outputs the following files:

* `sequence_data_*.tsv`: Table of the vertex data for the variation distance graphs. Each row represents a single vertex.

* `edge_data_*.tsv`: Table of edge data for variation-distance graphs. Each row represents a single edge.

* `distance_matrix_*.tsv`: Table of pairwise Levenshtein distances of the vertices.

* `graph_stats_*.tsv`: Summary statistics of the graph (e.g., number of vertices).

Arguments:

* `--input`: Directory with output from [`2_get_window_data`](#2getwindowdata) stage.

* `--output`: Output directory.

* `--subst_type`: Whether to process the files with/without substitutions.

### `4_get_histogram_data`

#### `get_histogram_data.py`

Pre-computes the necessary data for the variation-position histograms. Uses as input the output of the [`3_get_graph_data`](#3getgraphdata) stage and outputs the following files:

* `variation_*.tsv`: Contains data on individual variations of each sequence window in the `sequence_*.tsv` files. Each row represents a single variation (insertion, deletion, or substitution) in a sequence alignment.

* `variation_grouped_*.tsv`: The same as data as the corresponding `variation_*.tsv` grouped by (1) position of variation on the reference sequence, (2) total number of variations in the parent sequence, and (3) the type of variation (insertion, deletion, or substitution). The resulting frequencies are the sum of all individual variations with the same grouping characteristics.

Arguments:

* `--input`: Directory with output from [`3_get_graph_data`](#3getgraphdata).

* `--output`: Output directory.

* `--subst_type`: Whether to process the files with/without substitutions.

### `5_plot_graph`

Code for laying out and plotting variation-distance graphs. Uses the output from [`3_get_graph_data`](#3getgraphdata) as input. For more information about how the vertex sizes, vertex colors, edges, and layouts are computed please see the [publication](#citation).

A brief description of the layouts:

* `radial`: Arranges vertices in a radial grid around the reference sequence. Insertion (deletion) vertices are placed above (below) the reference sequence. Heuristics are used to improve aesthetics of vertex placement, though this may not easily generalize to new data.

* `mds`: Uses multi-dimensional scaling (MDS) to embed the pairwise Levenshtein distance matrix into 2D. 

* `kamada`: Uses the [NetworkX](http://networkx.org) package's implementation of the Kamada-Kawaii algorithm to layout the graph. Reference: T. Kamada and S. Kawaii. An algorithm for drawing general undirected graphs. *Inform. Process. Lett.*, 31:7–15, 1989.

* `universal`: Uses a deterministic layout described in the Methods of the [publication](#citation).

* `fractal`: Lays out only intsertion vertices. Places them in a fractal-like grid by repeatedly subdividing a square into 4 smaller squares. Each successive nucleotide of the inserted sequence determines which square is selected. `A` goes top-left, `C` goes top-right, `G` goes bottom-left, and `T` goes bottom-right.

#### `get_precomputed_layout.py`

Pre-computes the layout for groups of experiments by taking the union of all sequences in the specified experiments and laying them out. This allows experiments with the same reference sequence to have the same coordinates for the same sequence. All experiments should have the same windowed reference sequence and be individual (not comparison) experiments. 

See code of [`plot_graph.py`](#plotgraphpy) for the implementation of each layout.

Arguments:

* `--input`: List of data directories created with [`get_graph_data.py`](#getgraphdatapy). The sequencing in these data directories will be combined to produce the layout.

* `--output`: Output directory. Two files are created within this directory, paralleling the output of [`get_graph_data.py`](#getgraphdatapy): `sequence_data_*.py` and `edge_data_*.py`.

* `--reverse_complement`: Whether to reverse complement the sequence in the corresponding input. If present, the number of values must be the same as the number of input directories. `1` means reverse complement the sequences and `0` means do not. Used for making a common layout for data sets that have window reference sequences that are the reverse complements of each other (e.g., 2-DSB forward- and reverse-strand libraries).

* `--subst_type`: Whether to process files with/without substitutions.

* `--layout`: Choices: `radial`, `mds`, `kamada`, or `universal`. See the introduction of [`5_plot_graph`](#5plotgraph) for a description of each layout.

#### `plot_graph.py`

Does layout and plotting of the processed graph data from [`3_get_graph_data`](#3getgraphdata).

Arguments:

* `--input`: Directory with the data files produced with [`3_get_graph_data`](#3getgraphdata).

* `--output`: Output directory. Optional. If not given, no output will be written (useful only when using `--interactive`). If given, the output file has the same name as the `INPUT` directory with the appropriate files extension (e.g., `.png` or `.html`).

* `--ext`: File type of output. Choices: `png` or `html`. The library used to produce the output is [Plotly](https://plotly.com/). PNG output is a static image. HTML output is interactive and allows zooming/panning and provides details on vertices/edges through hover boxes.

* `--title`: If present, shows a title describing the experiment using the metadata stored in the `INPUT` directory.

* `--layout`: Layout algorithms to use for representing graph vertices in 2D. Choices: `kamada`, `radial`, `mds`, `universal`, `fractal`. See [`get_precomputed_layout.py`](#getprecomputedlayoutpy) for more details on each layout.

* `--universal_layout_y_axis_x_pos`: If present, shows a y-axis at the given x-position on the universal layout showing the distances to the reference. This value is in vertex coordinate space which has been arbitrarily set to be typically be in the region [-20, 20] x [-20, 20] for sequences in our data set. The value must thus be tweaked by hand to look right.

* `--universal_layout_x_axis_deletion_y_pos`: If present, shows an x-axis for deletions at the given y-position on the universal layout showing the approximate position of the deleted ranges. This value is in vertex coordinate space which has been arbitrarily set to be typically be in the region [-20, 20] x [-20, 20] for sequences in our data set. The value must thus be tweaked by hand to look right.

* `--universal_layout_x_axis_insertion_y_pos`: If present, shows an x-axis for insertions at the given y-position on the universal layout showing the first nucleotide of inserted sequences. This value is in vertex coordinate space which has been arbitrarily set to be typically be in the region [-20, 20] x [-20, 20] for sequences in our data set. The value must thus be tweaked by hand to look right.

* `--universal_layout_y_axis_y_range`: If showing an y-axis for the universal layout, the min and max y-position of the line. This value is in vertex coordinate space which has been arbitrarily set to be typically be in the region [-20, 20] x [-20, 20] for sequences in our data set. The value must thus be tweaked by hand to look right.

* `--universal_layout_x_axis_x_range`: If showing an x-axis for the universal layout, the min and max x-position of the line. This value is in vertex coordinate space which has been arbitrarily set to be typically be in the region [-20, 20] x [-20, 20] for sequences in our data set. The value must thus be tweaked by hand to look right.

* `--universal_layout_y_axis_deletion_max_tick`: If showing a y-axis for the universal layout, the max tick value for the deletion side. The ticks correspond to the number of deletions in the sequences in the corresponding row.

* `--universal_layout_y_axis_insertion_max_tick`: If showing a y-axis for the universal layout, the max tick value for the insertions side. The ticks correspond to the number of insertions in the sequences in the corresponding row.

* `--subst_type`: Whether to use the data files with or without substitutions.

* `--node_max_freq`: Max frequency to determine vertex size. All vertices with frequency >= `NODE_MAX_FREQ` get the largest possible vertex size.

* `--node_min_freq`: Max frequency to determine vertex size. All vertices with frequency <= `NODE_MIN_FREQ` get the smallest possible vertex size.

* `--node_max_px`: Largest possible vertex size in pixels.

* `--node_min_px`: Smallest possible vertex size in pixels.

* `--node_outline_scale`: How much to scale the node outline width (thickness). Values > 1 increase the width; values < 1 decrease the width.

* `--variation_types`: The variation types that should be included in the graph. This should be a sub-list of the types: `insertion`, `deletion`, `substitution`, `none`. Default value: [`insertion`, `deletion`, `none`]. `none` means the reference sequence.

* `--edge_scale`: How much to scale the edges width (thickness). Values > 1 increase the width; values < 1 decrease the width.

* `--width_px`: The width of the plot in pixels.

* `--height_px`: The height of the plot in pixels.

* `--line_width_scale`: How much to scale the line widths (AKA thickness). Values > 1 increase the width; values < 1 decrease the width.

* `--font_size_scale`: How much to scale the font size. Values > 1 increase the font size; values < 1 decrease the font size.

* `--precomputed_layout_dir`: If present, gives the directory where the precomputed layouts are located. This directory must contain the output of [`get_precomputed_layout.py`](#getprecomputedlayoutpy). If not present, the layout is computed newly.

* `--reverse_complement`: If present, uses the reverse complement of sequences when determining the node positions, and displaying labels and hover-text. This affects the precomputed layout, universal layout, and fractal layout.

* `--crop_x`: Range of the horizontal dimension to crop. Specified with normalized coords in range [0, 1]. Cropping is applied as the last step, so does not affect other coordinates specified in the arguments.

* `--crop_y`: Range of the vertical dimension to crop. Specified in normalized coords in range [0, 1]. Cropping is applied as the last step, so does not affect other coordinates specified in the arguments.

* `--range_x`: Range of x-axis for plotting. If not specified will be chosen to either show all nodes or a preset value for the layout, depending on the layout. This value is in vertex coordinate space which has been arbitrarily set to be typically be in the region [-20, 20] x [-20, 20] for sequences in our data set. The value must thus be tweaked by hand to look right.

* `--range_y`: Range of y-axis for plotting. If not specified will be chosen to either show all nodes or a preset value for the layout, depending on the layout. This value is in vertex coordinate space which has been arbitrarily set to be typically be in the region [-20, 20] x [-20, 20] for sequences in our data set. The value must thus be tweaked by hand to look right.

* `--legend`: If present, show a legend on the figure.

* `--legend_color_bar_scale`: How much to scale the legend color bar (for frequency ratio coloring). Values > 1 increase the size and values < 1 decrease the size.

* `--separate_components`: If present, separate the connected components of the graph. Currently, this means that the plot is partitioned into a grid with the reference sequence component getting the largest space and other components getting smaller spaces along the border.

### `6_plot_histogram`

#### `plot_histogram.py`

Plots the 3D histograms showing the variation-type/position/frequency of the DSB-sequence windows of experiments. Uses the output of [`4_get_histogram_data`](#4gethistogramdata). Uses [Matplotlib](https://matplotlib.org/) for the plotting.

Arguments:

* `--input`: Directory with the data files which are output from [`4_get_histogram_data`](#4gethistogramdata).

* `--output`: Output directory. Three images are output for each of the 3 variations: insertion, deletion, substitution. The files are named with the name of the `INPUT` directory and suffixed by the variation type.

* `--reverse_pos`: If present, reverse the x-axis positions. Useful if the input is from reverse-strand data and you want to compare it with forward-strand data.

* `--label_type`: Whether to index the x-axis by `absolute` positions on the reference sequence from `1` to `ref_length`, or `relative` positions from `-ref_length/2` to `ref_length/2` (skipping 0). The relative labeling assumes that the DSB is between `ref_length/2` and `ref_length/2 + 1`. Here `ref_length` refers to the length of windowed reference sequence, which in current analyses is always 20.

### `7_get_pptx`

#### `get_pptx.py`

Arrange images into a grid with labels and legend in a PPTX file. Uses the [Python-pptx](https://python-pptx.readthedocs.io/) package.

Arguments:

* `--input`: List of images to include in the grid.

* `--top_margin_labels`: Labels in the top margins of the grids. The number of arguments must be the same as the sum of the `NUM_COLS` arguments.

* `--left_margin_labels`: Labels in the left margins of the grids. The number of arguments must be the same as the sum of the `NUM_ROWS` arguments.

* `--output`: Output PPTX file.

* `--num_grids`: Number of separate grids to create. Each grid is stacked vertically below the previous grid.

* `--num_rows`: Number of rows in each grid. The number of arguments should equal `NUM_GRIDS`:.

* `--num_cols`: Number of columns in each grid. The number of arguments should equal `NUM_GRIDS`:.

* `--total_width`: Fraction of the page width to use for each grid. The number of arguments should equal `NUM_GRIDS`.

* `--node_max_freq`: The corresponding argument from [`plot_graph.py`](#plotgraphpy). Used to create the PPTX vertex size legends using PPTX shapes with the proper sizes. The number of arguments should equal `NUM_GRIDS`.

* `--node_min_freq`: The corresponding argument from [`plot_graph.py`](#plotgraphpy). Used to create the PPTX vertex size legends using PPTX shapes with the proper sizes. The number of arguments should equal `NUM_GRIDS`.

* `--node_max_px`: The corresponding argument from [`plot_graph.py`](#plotgraphpy). Used to create the PPTX vertex size legends using PPTX shapes with the proper sizes. The number of arguments should equal `NUM_GRIDS`.

* `--node_min_px`: The corresponding argument from [`plot_graph.py`](#plotgraphpy). Used to create the PPTX vertex size legends using PPTX shapes with the proper sizes. The number of arguments should equal `NUM_GRIDS`.

* `--title`: Page title displayed at the top of the page. Shows the directory name of `INPUT`.

* `--legends`: Which legends to draw in the file. These legends include the vertex size legend, vertex variation type color legend, vertex frequency ratio color legend, and edge legend. The legends are drawn outside the page limits so they will later have to be positioned by hand manually.

* `--template`: The PPTX file to use as a template. Currently, this is only used for determining the page size. The default template should be located at `7_get_pptx/template.pptx`.