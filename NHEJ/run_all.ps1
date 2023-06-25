foreach (
  $dir in @(
    "data_0_sam",
    "data_1_filter_nhej",
    "data_2_combine_repeat",
    "data_3_window",
    "data_4_graph",
    "data_5_histogram",
    "data_6_precomputed_layout",
    "data_bowtie2_build",
    "data_fastq"
  )
) {
  if (!(Test-Path $dir)) {
    mkdir $dir
  }
}
./run_00_bowtie2_align.ps1
./run_01_process_nhej.ps1
./run_02_combine_repeat.ps1
./run_03_get_window.ps1
./run_04_graph_data.ps1
./run_05_histogram_data.ps1
./run_06_precomputed_layout.ps1
./run_07_plot_graph_html.ps1
./run_07_plot_graph_png.ps1
./run_08_plot_histogram.ps1
./run_09_pptx_graph.ps1
./run_10_pptx_histogram.ps1
./run_11_plot_graph_main_png.ps1
