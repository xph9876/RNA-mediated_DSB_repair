param (
  [string]$layout = "radial"
)

$layout = "universal"

foreach ($cell in "WT", "KO") {
  # 1 DSB graphs
  python 3_make_pptx/make_pptx.py -i "plots/graphs/${layout}/individual/png/${cell}_sgA_R1_sense.png" "plots/graphs/${layout}/individual/png/${cell}_sgA_R1_branch.png" "plots/graphs/${layout}/individual/png/${cell}_sgA_R1_cmv.png" "plots/graphs/${layout}/individual/png/${cell}_sgB_R2_sense.png" "plots/graphs/${layout}/individual/png/${cell}_sgB_R2_branch.png" "plots/graphs/${layout}/individual/png/${cell}_sgB_R2_cmv.png" "plots/graphs/${layout}/combined/png/${cell}_sgA_R1_sense_branch.png" "plots/graphs/${layout}/combined/png/${cell}_sgA_R1_sense_cmv.png" "plots/graphs/${layout}/combined/png/${cell}_sgB_R2_sense_branch.png" "plots/graphs/${layout}/combined/png/${cell}_sgB_R2_sense_cmv.png" -lab "sgRNA A\nSense" "sgRNA A\nBranchΔ" "sgRNA A\nCMVΔ" "sgRNA B\nSense" "sgRNA B\nBranchΔ" "sgRNA B\nCMVΔ" "sgRNA A\nSense & BranchΔ" "sgRNA A\nSense & CMVΔ" "sgRNA B\nSense & BranchΔ" "sgRNA B\nSense & CMVΔ" -ng 2 -nr 2 2 -nc 3 2 -o "pptx/${cell}_1DSB_graphs.pptx" --legends node_size freq_ratio_sense_branch freq_ratio_sense_cmv node_type edge_type

  # 2 DSB graphs
  python 3_make_pptx/make_pptx.py -i "plots/graphs/${layout}/individual/png/${cell}_sgAB_R1_sense.png" "plots/graphs/${layout}/individual/png/${cell}_sgAB_R1_branch.png" "plots/graphs/${layout}/individual/png/${cell}_sgAB_R1_cmv.png" "plots/graphs/${layout}/individual/png/${cell}_sgAB_R2_sense.png" "plots/graphs/${layout}/individual/png/${cell}_sgAB_R2_branch.png" "plots/graphs/${layout}/individual/png/${cell}_sgAB_R2_cmv.png" "plots/graphs/${layout}/combined/png/${cell}_sgAB_R1_sense_branch.png" "plots/graphs/${layout}/combined/png/${cell}_sgAB_R1_sense_cmv.png" "plots/graphs/${layout}/combined/png/${cell}_sgAB_R2_sense_branch.png" "plots/graphs/${layout}/combined/png/${cell}_sgAB_R2_sense_cmv.png" -lab "sgRNA A & B\nForward strand\nSense" "sgRNA A & B\nForward strand\nBranchΔ" "sgRNA A & B\nForward strand\nCMVΔ" "sgRNA A & B\nReverse strand\nSense" "sgRNA A & B\nReverse strand\nBranchΔ" "sgRNA A & B\nReverse strand\nCMVΔ" "sgRNA A & B\nForward strand\nSense & BranchΔ" "sgRNA A & B\nForward strand\nSense & CMVΔ" "sgRNA A & B\nReverse strand\nSense & BranchΔ" "sgRNA A & B\nReverse strand\nSense & CMVΔ" -ng 2 -nr 2 2 -nc 3 2 -o "pptx/${cell}_2DSB_graphs.pptx" --legends node_size freq_ratio_sense_branch freq_ratio_sense_cmv node_type edge_type

  # Antisense graphs
  if ($cell -eq "WT") {
    python 3_make_pptx/make_pptx.py -i "plots/graphs/${layout}/individual/png/${cell}_sgCD_R1_antisense.png" "plots/graphs/${layout}/individual/png/${cell}_sgCD_R1_splicing.png" "plots/graphs/${layout}/individual/png/${cell}_sgCD_R1_antisense.png" "plots/graphs/${layout}/individual/png/${cell}_sgCD_R2_splicing.png" "plots/graphs/${layout}/combined/png/${cell}_sgCD_R1_antisense_splicing.png" "plots/graphs/${layout}/combined/png/${cell}_sgCD_R2_antisense_splicing.png" -lab "sgRNA C & D\nForward strand\nAntisense" "sgRNA C' & D\nForward strand\n5' splicingΔ" "sgRNA C & D\nForward strand\nAntisense" "sgRNA C' & D\nForward strand\n5' splicingΔ" "sgRNA C/C' & D\nForward strand\nAntisense & 5' splicingΔ" "sgRNA C/C' & D\nReverse strand\nAntisense & 5' splicingΔ" -ng 2 -nr 2 1 -nc 2 2 -o "pptx/${cell}_2DSBanti_graphs.pptx" --legends node_size freq_ratio_antisense_splicing node_type edge_type
  }
}
