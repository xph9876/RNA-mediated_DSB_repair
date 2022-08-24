python 3_graphs/get_common_layout.py --input data_4_graphs/KO_sgAB_R1_branch data_4_graphs/KO_sgAB_R2_branch data_4_graphs/KO_sgAB_R1_cmv data_4_graphs/KO_sgAB_R2_cmv data_4_graphs/KO_sgAB_R1_sense data_4_graphs/KO_sgAB_R2_sense data_4_graphs/WT_sgAB_R1_branch data_4_graphs/WT_sgAB_R2_branch data_4_graphs/WT_sgAB_R1_cmv data_4_graphs/WT_sgAB_R2_cmv data_4_graphs/WT_sgAB_R1_sense data_4_graphs/WT_sgAB_R2_sense --output data_6_layouts/universal/2DSB --reverse_complement 0 1 0 1 0 1 0 1 0 1 0 1 --subst_type withoutSubst --layout universal
python 3_graphs/get_common_layout.py --input data_4_graphs/KO_sgA_R1_branch_30bpDown data_4_graphs/KO_sgA_R1_cmv_30bpDown data_4_graphs/KO_sgA_R1_sense_30bpDown data_4_graphs/KO_sgA_R1_branch_noDSB data_4_graphs/KO_sgA_R1_cmv_noDSB data_4_graphs/KO_sgA_R1_sense_noDSB data_4_graphs/KO_sgA_R1_branch data_4_graphs/KO_sgA_R1_cmv data_4_graphs/KO_sgA_R1_sense data_4_graphs/WT_sgA_R1_branch_30bpDown data_4_graphs/WT_sgA_R1_cmv_30bpDown data_4_graphs/WT_sgA_R1_sense_30bpDown data_4_graphs/WT_sgA_R1_branch_noDSB data_4_graphs/WT_sgA_R1_cmv_noDSB data_4_graphs/WT_sgA_R1_sense_noDSB data_4_graphs/WT_sgA_R1_branch data_4_graphs/WT_sgA_R1_cmv data_4_graphs/WT_sgA_R1_sense --output data_6_layouts/universal/1DSB_A --reverse_complement 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 --subst_type withoutSubst --layout universal
python 3_graphs/get_common_layout.py --input data_4_graphs/KO_sgB_R2_branch_30bpDown data_4_graphs/KO_sgB_R2_cmv_30bpDown data_4_graphs/KO_sgB_R2_sense_30bpDown data_4_graphs/KO_sgB_R2_branch_noDSB data_4_graphs/KO_sgB_R2_cmv_noDSB data_4_graphs/KO_sgB_R2_sense_noDSB data_4_graphs/KO_sgB_R2_branch data_4_graphs/KO_sgB_R2_cmv data_4_graphs/KO_sgB_R2_sense data_4_graphs/WT_sgB_R2_branch_30bpDown data_4_graphs/WT_sgB_R2_cmv_30bpDown data_4_graphs/WT_sgB_R2_sense_30bpDown data_4_graphs/WT_sgB_R2_branch_noDSB data_4_graphs/WT_sgB_R2_cmv_noDSB data_4_graphs/WT_sgB_R2_sense_noDSB data_4_graphs/WT_sgB_R2_branch data_4_graphs/WT_sgB_R2_cmv data_4_graphs/WT_sgB_R2_sense --output data_6_layouts/universal/1DSB_B --reverse_complement 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 --subst_type withoutSubst --layout universal
python 3_graphs/get_common_layout.py --input data_4_graphs/WT_sgCD_R1_antisense_merged data_4_graphs/WT_sgCD_R1_antisense_new data_4_graphs/WT_sgCD_R1_antisense_old data_4_graphs/WT_sgCD_R2_antisense_merged data_4_graphs/WT_sgCD_R2_antisense_new data_4_graphs/WT_sgCD_R2_antisense_old data_4_graphs/WT_sgCD_R1_splicing_merged data_4_graphs/WT_sgCD_R1_splicing_new data_4_graphs/WT_sgCD_R1_splicing_old data_4_graphs/WT_sgCD_R2_splicing_merged data_4_graphs/WT_sgCD_R2_splicing_new data_4_graphs/WT_sgCD_R2_splicing_old --output data_6_layouts/universal/2DSBanti --reverse_complement 0 0 0 1 1 1 0 0 0 1 1 1 --subst_type withoutSubst --layout universal
