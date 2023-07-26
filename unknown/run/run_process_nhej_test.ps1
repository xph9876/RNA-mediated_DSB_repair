# Code for running the 'filter_nhej.py' script on the test data
# Working directory must be 'NHEJ' to run this script
python 1_process_nhej/filter_nhej.py --sam_file data_0_sam/yjl255_WT_sgA_R1_sense.sam --ref_seq_file ref_seq/1DSB_R1_sense.fa --output data_1_filter_nhej/yjl255_WT_sgA_R1_sense.tsv --output_rejected data_1_filter_nhej/yjl255_WT_sgA_R1_sense_rejected.tsv --min_length 130 --dsb_pos 67 > ..\unknown\output\logs\yjl255_WT_sgA_R1_sense.txt
python 1_process_nhej/filter_nhej.py --sam_file data_0_sam/yjl259_WT_sgA_R1_branch.sam --ref_seq_file ref_seq/1DSB_R1_branch.fa --output data_1_filter_nhej/yjl259_WT_sgA_R1_branch.tsv --output_rejected data_1_filter_nhej/yjl259_WT_sgA_R1_branch_rejected.tsv --min_length 130 --dsb_pos 67 > ..\unknown\output\logs\yjl259_WT_sgA_R1_branch.txt
python 1_process_nhej/filter_nhej.py --sam_file data_0_sam/yjl267_WT_sgB_R2_sense.sam --ref_seq_file ref_seq/1DSB_R2_sense.fa --output data_1_filter_nhej/yjl267_WT_sgB_R2_sense.tsv --output_rejected data_1_filter_nhej/yjl267_WT_sgB_R2_sense_rejected.tsv --min_length 130 --dsb_pos 46 > ..\unknown\output\logs\yjl267_WT_sgB_R2_sense.txt
python 1_process_nhej/filter_nhej.py --sam_file data_0_sam/yjl271_WT_sgB_R2_branch.sam --ref_seq_file ref_seq/1DSB_R2_branch.fa --output data_1_filter_nhej/yjl271_WT_sgB_R2_branch.tsv --output_rejected data_1_filter_nhej/yjl271_WT_sgB_R2_branch_rejected.tsv --min_length 130 --dsb_pos 46 > ..\unknown\output\logs\yjl271_WT_sgB_R2_branch.txt
python 1_process_nhej/filter_nhej.py --sam_file data_0_sam/yjl244_WT_sgA_R1_sense_noDSB.sam --ref_seq_file ref_seq/1DSB_R1_sense.fa --output data_1_filter_nhej/yjl244_WT_sgA_R1_sense_noDSB.tsv --output_rejected data_1_filter_nhej/yjl244_WT_sgA_R1_sense_noDSB_rejected.tsv --min_length 130 --dsb_pos 67 > ..\unknown\output\logs\yjl244_WT_sgA_R1_sense_noDSB.txt
python 1_process_nhej/filter_nhej.py --sam_file data_0_sam/yjl244_WT_sgB_R2_sense_noDSB.sam --ref_seq_file ref_seq/1DSB_R2_sense.fa --output data_1_filter_nhej/yjl244_WT_sgB_R2_sense_noDSB.tsv --output_rejected data_1_filter_nhej/yjl244_WT_sgB_R2_sense_noDSB_rejected.tsv --min_length 130 --dsb_pos 46 > ..\unknown\output\logs\yjl244_WT_sgB_R2_sense_noDSB.txt
python 1_process_nhej/filter_nhej.py --sam_file data_0_sam/yjl245_WT_sgA_R1_branch_noDSB.sam --ref_seq_file ref_seq/1DSB_R1_branch.fa --output data_1_filter_nhej/yjl245_WT_sgA_R1_branch_noDSB.tsv --output_rejected data_1_filter_nhej/yjl245_WT_sgA_R1_branch_noDSB_rejected.tsv --min_length 130 --dsb_pos 67 > ..\unknown\output\logs\yjl245_WT_sgA_R1_branch_noDSB.txt
python 1_process_nhej/filter_nhej.py --sam_file data_0_sam/yjl245_WT_sgB_R2_branch_noDSB.sam --ref_seq_file ref_seq/1DSB_R2_branch.fa --output data_1_filter_nhej/yjl245_WT_sgB_R2_branch_noDSB.tsv --output_rejected data_1_filter_nhej/yjl245_WT_sgB_R2_branch_noDSB_rejected.tsv --min_length 130 --dsb_pos 46 --quiet > ..\unknown\output\logs\yjl245_WT_sgB_R2_branch_noDSB.txt
