dirs=(
  "data_0_sam"
  "data_1_filter_nhej"
  "data_2_combine_repeat"
  "data_3_window"
  "data_4_graph"
  "data_5_histogram"
  "data_6_precomputed_layout"
  "data_bowtie2_build"
  "data_fastq"
)
for dir in $dirs; do
    if [ ! -d $dir ]; then
        mkdir $dir
    fi
done
./run_00_bowtie2_align.sh
./run_01_process_nhej.sh
./run_02_combine_repeat.sh
./run_03_get_window.sh
./run_04_graph_data.sh
./run_05_histogram_data.sh
./run_06_precomputed_layout.sh
./run_07_plot_graph_html.sh
./run_07_plot_graph_png.sh
./run_08_plot_histogram.sh
./run_09_pptx_graph.sh
./run_10_pptx_histogram.sh
./run_11_plot_graph_main_png.sh
