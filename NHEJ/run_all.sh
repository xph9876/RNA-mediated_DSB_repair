./run/run_00_make_dirs.sh
./run/run_00_bowtie2_align.sh
./run/run_01_process_nhej.sh
./run/run_02_combine_repeat.sh
./run/run_03_get_window.sh
./run/run_04_graph_data.sh
./run/run_05_histogram_data.sh
./run/run_06_precomputed_layout.sh
./run/run_07_plot_graph_html.sh
./run/run_07_plot_graph_png.sh
./run/run_08_plot_histogram.sh
./run/run_09_pptx_graph.sh
./run/run_10_pptx_histogram.sh
./run/run_11_plot_graph_main_png.sh
python ./8_freq_analyze/freq_analyze.py --output ./data_7_freq_analyze/
