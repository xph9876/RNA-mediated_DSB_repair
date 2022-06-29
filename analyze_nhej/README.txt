Description
-----------
This folder contains the processing pipeline for the NHEJ (nonhomologous end joining) DNA repair.

Pipeline stages
---------------
1. run_01_process_nhej - Filter the raw SAM files to retain the NHEJ patterns.
2. run_02_combine_repeats - Combine repeat experiments and retain only common sequences.
3. run_03_make_main_data_withSubst, run_03_make_main_data_withoutSubst - Make processed/filtered data files for the variation-position graphs and 3D histograms. The "withSubst" script keeps the substitutions, and "withoutSubst" script ignores them.
4. run_04_make_graph_data_withSubst, run_04_make_graph_data_withoutSubst - Make processed data files for generating the variation-position graphs. The "withSubst"/"withoutSubst" script uses the respective output of the previous stage.
5. run_05_make_main_data_combined_withSubst, run_05_make_main_data_combined_withoutSubst - Makes the combined data files for comparing 2 experiments.
6. run_05_make_graph_data_combined_withSubst, run_05_make_graph_data_combined_withoutSubst - Makes the combined data files for comparing 2 experiments.
7. run_07_make_historam_3d - Plots the 3D variation-position histograms.
8. run_08_common_layout - Create the common layouts for the plotting the variation-position graphs. This insure that comparable experiments have vertices with the same variation placed in the same position.
9. run_09_plot_graph - Plots the variation-position graphs.
10. run_10_make_pptx_graph - Arranges the variation-position graphs in a grid in PPTX files.
11. run_11_make_ppts_histograms - Arranges the variation-position histograms in a grid in PPTX files.

Running
-------
To run all stages of the pipeline use the "run_all.bat" file (Windows cmd) or "run_all.sh" file (Unix bash). To run one of the individual stages above use the corresponding "run*.bat" file (Windows cmd) or "run*.sh" file (Unix bash). To run individual Python files in the pipeline use Python 3 and refer to the .bat or .sh scripts for the usage, or run the files without any arguments to print a help message. These scripts have been tested on Windows 11, Python 3.10.4. Note, the steps in the pipeline must be run in proper order as later stages depend on output from previous stages.