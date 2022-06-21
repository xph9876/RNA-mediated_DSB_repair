param (
  [string]$layout = "radial"
)

if ($layout -eq "radial") {
  python 2_graph_processing/make_common_layout.py -i libraries_4/WT_sgAB_R1_sense libraries_4/WT_sgAB_R1_branch libraries_4/WT_sgAB_R1_cmv libraries_4/KO_sgAB_R1_sense libraries_4/KO_sgAB_R1_branch libraries_4/KO_sgAB_R1_cmv libraries_4/WT_sgAB_R2_sense libraries_4/WT_sgAB_R2_branch libraries_4/WT_sgAB_R2_cmv libraries_4/KO_sgAB_R2_sense libraries_4/KO_sgAB_R2_branch libraries_4/KO_sgAB_R2_cmv -rc 0 0 0 0 0 0 1 1 1 1 1 1 -o layouts/2DSB_AB --subst_type without --layout radial
  python 2_graph_processing/make_common_layout.py -i libraries_4/WT_sgA_R1_sense libraries_4/WT_sgA_R1_branch libraries_4/WT_sgA_R1_cmv libraries_4/KO_sgA_R1_sense libraries_4/KO_sgA_R1_branch libraries_4/KO_sgA_R1_cmv -o layouts/1DSB_A --subst_type without --layout radial
  python 2_graph_processing/make_common_layout.py -i libraries_4/WT_sgB_R2_sense libraries_4/WT_sgB_R2_branch libraries_4/WT_sgB_R2_cmv libraries_4/KO_sgB_R2_sense libraries_4/KO_sgB_R2_branch libraries_4/KO_sgB_R2_cmv -o layouts/1DSB_B --subst_type without --layout radial
  python 2_graph_processing/make_common_layout.py -i libraries_4/WT_sgCD_R1_antisense libraries_4/WT_sgCD_R1_splicing libraries_4/WT_sgCD_R2_antisense libraries_4/WT_sgCD_R2_splicing -rc 0 0 1 1 -o layouts/2DSBanti_CD --subst_type without --layout radial
} elseif ($layout -eq "universal") {
  python 2_graph_processing/make_common_layout.py -i libraries_4/WT_sgAB_R1_sense libraries_4/WT_sgAB_R1_branch libraries_4/WT_sgAB_R1_cmv libraries_4/KO_sgAB_R1_sense libraries_4/KO_sgAB_R1_branch libraries_4/KO_sgAB_R1_cmv libraries_4/WT_sgAB_R2_sense libraries_4/WT_sgAB_R2_branch libraries_4/WT_sgAB_R2_cmv libraries_4/KO_sgAB_R2_sense libraries_4/KO_sgAB_R2_branch libraries_4/KO_sgAB_R2_cmv -rc 0 0 0 0 0 0 1 1 1 1 1 1 -o layouts/2DSB_AB --subst_type without --layout universal
  python 2_graph_processing/make_common_layout.py -i libraries_4/WT_sgA_R1_sense libraries_4/WT_sgA_R1_branch libraries_4/WT_sgA_R1_cmv libraries_4/KO_sgA_R1_sense libraries_4/KO_sgA_R1_branch libraries_4/KO_sgA_R1_cmv -o layouts/1DSB_A --subst_type without --layout universal
  python 2_graph_processing/make_common_layout.py -i libraries_4/WT_sgB_R2_sense libraries_4/WT_sgB_R2_branch libraries_4/WT_sgB_R2_cmv libraries_4/KO_sgB_R2_sense libraries_4/KO_sgB_R2_branch libraries_4/KO_sgB_R2_cmv -rc 1 1 1 1 1 1 -o layouts/1DSB_B --subst_type without --layout universal
  python 2_graph_processing/make_common_layout.py -i libraries_4/WT_sgCD_R1_antisense libraries_4/WT_sgCD_R1_splicing libraries_4/WT_sgCD_R2_antisense libraries_4/WT_sgCD_R2_splicing -rc 0 0 1 1 -o layouts/2DSBanti_CD --subst_type without --layout universal
}