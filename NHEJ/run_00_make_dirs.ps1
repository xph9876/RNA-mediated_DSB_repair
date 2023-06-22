foreach (
  $dir in @(
    "data_0_sam",
    "data_1_filter_nhej",
    "data_2_combine_repeat",
    "data_3_window",
    "data_4_graph",
    "data_5_histogram",
    "data_6_precomputed_layout",
    "data_bowtie2_build"
  )
) {
  if (!(Test-Path $dir)) {
    New-Item -ItemType Directory $dir
  }
}
