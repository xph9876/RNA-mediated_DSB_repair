# 1. Alignment.
bowtie2-build ref_seq/test.fa test_bowtie
bowtie2 -x test_bowtie data_demo_fastq/test_1.fastq -S data_0_sam/test_1.sam
bowtie2 -x test_bowtie data_demo_fastq/test_2.fastq -S data_0_sam/test_2.sam
bowtie2 -x test_bowtie data_demo_fastq/test_3.fastq -S data_0_sam/test_3.sam
bowtie2 -x test_bowtie data_demo_fastq/test_4.fastq -S data_0_sam/test_4.sam

# 2. NHEJ filtering.
python 1_process_nhej\filter_nhej.py --sam_file data_0_sam/test_1.sam --ref_seq_file ref_seq/test.fa --output data_1_filter_nhej/test_1.tsv --min_length 50 --dsb_pos 50 --quiet
python 1_process_nhej/filter_nhej.py --sam_file data_0_sam/test_2.sam --ref_seq_file ref_seq/test.fa --output data_1_filter_nhej/test_2.tsv --min_length 50 --dsb_pos 50 --quiet
python 1_process_nhej/filter_nhej.py --sam_file data_0_sam/test_3.sam --ref_seq_file ref_seq/test.fa --output data_1_filter_nhej/test_3.tsv --min_length 50 --dsb_pos 50 --quiet
python 1_process_nhej/filter_nhej.py --sam_file data_0_sam/test_4.sam --ref_seq_file ref_seq/test.fa --output data_1_filter_nhej/test_4.tsv --min_length 50 --dsb_pos 50 --quiet

# 3. Combine repeats.
python 1_process_nhej/combine_repeat.py --input data_1_filter_nhej/test_1.tsv data_1_filter_nhej/test_2.tsv data_1_filter_nhej/test_3.tsv data_1_filter_nhej/test_4.tsv --output data_2_combine_repeat/test.tsv --quiet

# 4. Extract windows.
python 2_get_window_data/get_window.py --input data_2_combine_repeat/test.tsv --ref_seq_file ref_seq/test.fa --output data_3_window/test --dsb_pos 50 --dsb_type 2DSB --strand R1 --guide_rna sgTest --cell_line WT --construct Test --subst_type withSubst --control_type notControl --version versionNone
python 2_get_window_data/get_window.py --input data_2_combine_repeat/test.tsv --ref_seq_file ref_seq/test.fa --output data_3_window/test --dsb_pos 50 --dsb_type 2DSB --strand R1 --guide_rna sgTest --cell_line WT --construct Test --subst_type withoutSubst --control_type notControl --version versionNone

# 5. Get frequencies.
python 2_get_window_data/get_freq.py --input data_3_window/test --output data_3_window/test --subst_type withSubst --total_reads 3000 3000 3000 3000 --freq_min 1e-5
python 2_get_window_data/get_freq.py --input data_3_window/test --output data_3_window/test --subst_type withoutSubst --total_reads 3000 3000 3000 3000 --freq_min 1e-5

# 6. Get graph data.
python 3_get_graph_data/get_graph_data.py --input data_3_window/test --output data_4_graph/test --subst_type withSubst
python 3_get_graph_data/get_graph_data.py --input data_3_window/test --output data_4_graph/test --subst_type withoutSubst

# 7. Get histogram data.
python 4_get_histogram_data/get_histogram_data.py --input data_4_graph/test --output data_5_histogram/test --subst_type withSubst
python 4_get_histogram_data/get_histogram_data.py --input data_4_graph/test --output data_5_histogram/test --subst_type withoutSubst

# 8. Plot variation-distance graphs.
python 5_plot_graph/plot_graph.py --input data_4_graph/test --output plot/demo/html --ext html --layout universal  --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19
python 5_plot_graph/plot_graph.py --input data_4_graph/test --output plot/demo/ --ext png --layout universal  --width 2400 --height 1800 --range_x -12 13 --range_y -23 20 --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -20.5 18.5 --universal_layout_y_axis_insertion_max_tick 6 --universal_layout_y_axis_deletion_max_tick 19

# 9. Plot variation-position histograms.
python 6_plot_histogram\plot_histogram.py --input data_5_histogram\test --output plot\histogram\demo --label_type relative