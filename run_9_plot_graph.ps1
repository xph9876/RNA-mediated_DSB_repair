param (
  [string]$layout = "radial"
)

if ($layout -eq "radial") {
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgAB_R1_sense -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgAB_R1_branch -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgAB_R1_cmv -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgAB_R1_sense -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgAB_R1_branch -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgAB_R1_cmv -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgAB_R2_sense -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout radial --reverse_complement  --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgAB_R2_branch -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout radial --reverse_complement  --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgAB_R2_cmv -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout radial --reverse_complement  --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgAB_R2_sense -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout radial --reverse_complement  --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgAB_R2_branch -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout radial --reverse_complement  --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgAB_R2_cmv -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout radial --reverse_complement  --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgA_R1_sense -o plots/graphs/individual/png --layout_dir layouts/1DSB_A -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgA_R1_branch -o plots/graphs/individual/png --layout_dir layouts/1DSB_A -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgA_R1_cmv -o plots/graphs/individual/png --layout_dir layouts/1DSB_A -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgA_R1_sense -o plots/graphs/individual/png --layout_dir layouts/1DSB_A -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgA_R1_branch -o plots/graphs/individual/png --layout_dir layouts/1DSB_A -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgA_R1_cmv -o plots/graphs/individual/png --layout_dir layouts/1DSB_A -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgB_R2_sense -o plots/graphs/individual/png --layout_dir layouts/1DSB_B -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgB_R2_branch -o plots/graphs/individual/png --layout_dir layouts/1DSB_B -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgB_R2_cmv -o plots/graphs/individual/png --layout_dir layouts/1DSB_B -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgB_R2_sense -o plots/graphs/individual/png --layout_dir layouts/1DSB_B -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgB_R2_branch -o plots/graphs/individual/png --layout_dir layouts/1DSB_B -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgB_R2_cmv -o plots/graphs/individual/png --layout_dir layouts/1DSB_B -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgCD_R1_antisense -o plots/graphs/individual/png --layout_dir layouts/2DSBanti_CD -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgCD_R1_splicing -o plots/graphs/individual/png --layout_dir layouts/2DSBanti_CD -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgCD_R2_antisense -o plots/graphs/individual/png --layout_dir layouts/2DSBanti_CD -ext png --layout radial --reverse_complement  --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgCD_R2_splicing -o plots/graphs/individual/png --layout_dir layouts/2DSBanti_CD -ext png --layout radial --reverse_complement  --width 2400 --height 1863

  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgAB_R1_sense_branch -o plots/graphs/combined/png --layout_dir layouts/2DSB_AB -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgAB_R1_sense_cmv -o plots/graphs/combined/png --layout_dir layouts/2DSB_AB -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgAB_R1_sense_branch -o plots/graphs/combined/png --layout_dir layouts/2DSB_AB -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgAB_R1_sense_cmv -o plots/graphs/combined/png --layout_dir layouts/2DSB_AB -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgAB_R2_sense_branch -o plots/graphs/combined/png --layout_dir layouts/2DSB_AB -ext png --layout radial --reverse_complement  --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgAB_R2_sense_cmv -o plots/graphs/combined/png --layout_dir layouts/2DSB_AB -ext png --layout radial --reverse_complement  --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgAB_R2_sense_branch -o plots/graphs/combined/png --layout_dir layouts/2DSB_AB -ext png --layout radial --reverse_complement  --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgAB_R2_sense_cmv -o plots/graphs/combined/png --layout_dir layouts/2DSB_AB -ext png --layout radial --reverse_complement  --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgA_R1_sense_branch -o plots/graphs/combined/png --layout_dir layouts/1DSB_A -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgA_R1_sense_cmv -o plots/graphs/combined/png --layout_dir layouts/1DSB_A -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgA_R1_sense_branch -o plots/graphs/combined/png --layout_dir layouts/1DSB_A -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgA_R1_sense_cmv -o plots/graphs/combined/png --layout_dir layouts/1DSB_A -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgB_R2_sense_branch -o plots/graphs/combined/png --layout_dir layouts/1DSB_B -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgB_R2_sense_cmv -o plots/graphs/combined/png --layout_dir layouts/1DSB_B -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgB_R2_sense_branch -o plots/graphs/combined/png --layout_dir layouts/1DSB_B -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgB_R2_sense_cmv -o plots/graphs/combined/png --layout_dir layouts/1DSB_B -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgCD_R1_antisense_splicing -o plots/graphs/combined/png --layout_dir layouts/2DSBanti_CD -ext png --layout radial --width 2400 --height 1863
  # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgCD_R2_antisense_splicing -o plots/graphs/combined/png --layout_dir layouts/2DSBanti_CD -ext png --layout radial --reverse_complement  --width 2400 --height 1863
} elseif ($layout -eq "universal") {
  # RERUN THE UNIVERSAL!!!
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgAB_R1_sense -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -22 20 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgAB_R1_branch -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -22 20 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgAB_R1_cmv -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -22 20 
  python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgAB_R1_sense -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -22 20 
  python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgAB_R1_branch -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -22 20 
  python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgAB_R1_cmv -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -22 20 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgAB_R2_sense -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -22 20 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgAB_R2_branch -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -22 20 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgAB_R2_cmv -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -22 20 
  python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgAB_R2_sense -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -22 20 
  python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgAB_R2_branch -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -22 20 
  python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgAB_R2_cmv -o plots/graphs/individual/png --layout_dir layouts/2DSB_AB -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -22 20 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgA_R1_sense -o plots/graphs/individual/png --layout_dir layouts/1DSB_A -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -22 18 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgA_R1_branch -o plots/graphs/individual/png --layout_dir layouts/1DSB_A -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -22 18 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgA_R1_cmv -o plots/graphs/individual/png --layout_dir layouts/1DSB_A -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -22 18 
  python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgA_R1_sense -o plots/graphs/individual/png --layout_dir layouts/1DSB_A -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -22 18 
  python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgA_R1_branch -o plots/graphs/individual/png --layout_dir layouts/1DSB_A -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -22 18 
  python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgA_R1_cmv -o plots/graphs/individual/png --layout_dir layouts/1DSB_A -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -22 18 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgB_R2_sense -o plots/graphs/individual/png --layout_dir layouts/1DSB_B -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -20 10 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgB_R2_branch -o plots/graphs/individual/png --layout_dir layouts/1DSB_B -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -20 10 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgB_R2_cmv -o plots/graphs/individual/png --layout_dir layouts/1DSB_B -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -20 10 
  python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgB_R2_sense -o plots/graphs/individual/png --layout_dir layouts/1DSB_B -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -20 10 
  python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgB_R2_branch -o plots/graphs/individual/png --layout_dir layouts/1DSB_B -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -20 10 
  python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgB_R2_cmv -o plots/graphs/individual/png --layout_dir layouts/1DSB_B -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -20 10 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgCD_R1_antisense -o plots/graphs/individual/png --layout_dir layouts/2DSBanti_CD -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -18 24 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgCD_R1_splicing -o plots/graphs/individual/png --layout_dir layouts/2DSBanti_CD -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -18 24 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgCD_R2_antisense -o plots/graphs/individual/png --layout_dir layouts/2DSBanti_CD -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -18 24 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgCD_R2_splicing -o plots/graphs/individual/png --layout_dir layouts/2DSBanti_CD -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -18 24
    
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgAB_R1_sense_branch -o plots/graphs/combined/png --layout_dir layouts/2DSB_AB -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -22 20 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgAB_R1_sense_cmv -o plots/graphs/combined/png --layout_dir layouts/2DSB_AB -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -22 20 
  python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgAB_R1_sense_branch -o plots/graphs/combined/png --layout_dir layouts/2DSB_AB -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -22 20 
  python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgAB_R1_sense_cmv -o plots/graphs/combined/png --layout_dir layouts/2DSB_AB -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -22 20 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgAB_R2_sense_branch -o plots/graphs/combined/png --layout_dir layouts/2DSB_AB -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -22 20 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgAB_R2_sense_cmv -o plots/graphs/combined/png --layout_dir layouts/2DSB_AB -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -22 20 
  python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgAB_R2_sense_branch -o plots/graphs/combined/png --layout_dir layouts/2DSB_AB -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -22 20 
  python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgAB_R2_sense_cmv -o plots/graphs/combined/png --layout_dir layouts/2DSB_AB -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -22 20 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgA_R1_sense_branch -o plots/graphs/combined/png --layout_dir layouts/1DSB_A -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -22 18 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgA_R1_sense_cmv -o plots/graphs/combined/png --layout_dir layouts/1DSB_A -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -22 18 
  python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgA_R1_sense_branch -o plots/graphs/combined/png --layout_dir layouts/1DSB_A -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -22 18 
  python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgA_R1_sense_cmv -o plots/graphs/combined/png --layout_dir layouts/1DSB_A -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -22 18 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgB_R2_sense_branch -o plots/graphs/combined/png --layout_dir layouts/1DSB_B -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -20 10 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgB_R2_sense_cmv -o plots/graphs/combined/png --layout_dir layouts/1DSB_B -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -20 10 
  python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgB_R2_sense_branch -o plots/graphs/combined/png --layout_dir layouts/1DSB_B -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -20 10 
  python 2_graph_processing/plot_graph.py -i libraries_4/KO_sgB_R2_sense_cmv -o plots/graphs/combined/png --layout_dir layouts/1DSB_B -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -20 10 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgCD_R1_antisense_splicing -o plots/graphs/combined/png --layout_dir layouts/2DSBanti_CD -ext png --layout universal --width 2400 --height 1863  --range_x -12 12 --range_y -18 24 
  python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgCD_R2_antisense_splicing -o plots/graphs/combined/png --layout_dir layouts/2DSBanti_CD -ext png --layout universal --reverse_complement  --width 2400 --height 1863  --range_x -12 12 --range_y -18 24
}