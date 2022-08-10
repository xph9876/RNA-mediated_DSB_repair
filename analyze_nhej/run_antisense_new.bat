:: 02_combine_repeats
python 1_process_nhej/combine_repeats.py -i libraries_2/yjl349_WT_sgCD_R1_antisense.tsv libraries_2/yjl350_WT_sgCD_R1_antisense.tsv libraries_2/yjl351_WT_sgCD_R1_antisense.tsv libraries_2/yjl352_WT_sgCD_R1_antisense.tsv --total_reads 6302370 5988345 5741048 5405939 -o libraries_3_new/WT_sgCD_R1_antisense.tsv --quiet 
python 1_process_nhej/combine_repeats.py -i libraries_2/yjl353_WT_sgCD_R1_splicing.tsv libraries_2/yjl354_WT_sgCD_R1_splicing.tsv libraries_2/yjl355_WT_sgCD_R1_splicing.tsv libraries_2/yjl356_WT_sgCD_R1_splicing.tsv --total_reads 4177291 4647120 4405688 5054324 -o libraries_3_new/WT_sgCD_R1_splicing.tsv --quiet 
python 1_process_nhej/combine_repeats.py -i libraries_2/yjl349_WT_sgCD_R2_antisense.tsv libraries_2/yjl350_WT_sgCD_R2_antisense.tsv libraries_2/yjl351_WT_sgCD_R2_antisense.tsv libraries_2/yjl352_WT_sgCD_R2_antisense.tsv --total_reads 6161112 5856132 5613583 5281487 -o libraries_3_new/WT_sgCD_R2_antisense.tsv --quiet 
python 1_process_nhej/combine_repeats.py -i libraries_2/yjl353_WT_sgCD_R2_splicing.tsv libraries_2/yjl354_WT_sgCD_R2_splicing.tsv libraries_2/yjl355_WT_sgCD_R2_splicing.tsv libraries_2/yjl356_WT_sgCD_R2_splicing.tsv --total_reads 4062558 4535123 4280724 4926338 -o libraries_3_new/WT_sgCD_R2_splicing.tsv --quiet 

:: 03_make_main_data_withoutSubst
python 2_graph_processing/make_main_data.py -i libraries_3_new/WT_sgCD_R1_antisense.tsv -o libraries_4_new/WT_sgCD_R1_antisense -ref ref_seq/2DSBanti_R1_antisense.fa -dsb 66 --dsb_type 2a --strand R1 --hguide CD --cell WT --treatment antisense --subst_type without --control none
python 2_graph_processing/make_main_data.py -i libraries_3_new/WT_sgCD_R1_splicing.tsv -o libraries_4_new/WT_sgCD_R1_splicing -ref ref_seq/2DSBanti_R1_splicing.fa -dsb 66 --dsb_type 2a --strand R1 --hguide CD --cell WT --treatment splicing --subst_type without --control none
python 2_graph_processing/make_main_data.py -i libraries_3_new/WT_sgCD_R2_antisense.tsv -o libraries_4_new/WT_sgCD_R2_antisense -ref ref_seq/2DSBanti_R2_antisense.fa -dsb 47 --dsb_type 2a --strand R2 --hguide CD --cell WT --treatment antisense --subst_type without --control none
python 2_graph_processing/make_main_data.py -i libraries_3_new/WT_sgCD_R2_splicing.tsv -o libraries_4_new/WT_sgCD_R2_splicing -ref ref_seq/2DSBanti_R2_splicing.fa -dsb 47 --dsb_type 2a --strand R2 --hguide CD --cell WT --treatment splicing --subst_type without --control none

:: 03_make_main_data_withSubst
python 2_graph_processing/make_main_data.py -i libraries_3_new/WT_sgCD_R1_antisense.tsv -o libraries_4_new/WT_sgCD_R1_antisense -ref ref_seq/2DSBanti_R1_antisense.fa -dsb 66 --dsb_type 2a --strand R1 --hguide CD --cell WT --treatment antisense --subst_type with --control none
python 2_graph_processing/make_main_data.py -i libraries_3_new/WT_sgCD_R1_splicing.tsv -o libraries_4_new/WT_sgCD_R1_splicing -ref ref_seq/2DSBanti_R1_splicing.fa -dsb 66 --dsb_type 2a --strand R1 --hguide CD --cell WT --treatment splicing --subst_type with --control none
python 2_graph_processing/make_main_data.py -i libraries_3_new/WT_sgCD_R2_antisense.tsv -o libraries_4_new/WT_sgCD_R2_antisense -ref ref_seq/2DSBanti_R2_antisense.fa -dsb 47 --dsb_type 2a --strand R2 --hguide CD --cell WT --treatment antisense --subst_type with --control none
python 2_graph_processing/make_main_data.py -i libraries_3_new/WT_sgCD_R2_splicing.tsv -o libraries_4_new/WT_sgCD_R2_splicing -ref ref_seq/2DSBanti_R2_splicing.fa -dsb 47 --dsb_type 2a --strand R2 --hguide CD --cell WT --treatment splicing --subst_type with --control none

:: 04_make_graph_data_withoutSubst
python 2_graph_processing/make_graph_data.py --subst_type without -dir libraries_4_new/WT_sgCD_R1_antisense
python 2_graph_processing/make_graph_data.py --subst_type without -dir libraries_4_new/WT_sgCD_R1_splicing
python 2_graph_processing/make_graph_data.py --subst_type without -dir libraries_4_new/WT_sgCD_R2_antisense
python 2_graph_processing/make_graph_data.py --subst_type without -dir libraries_4_new/WT_sgCD_R2_splicing

:: 04_make_graph_data_withSubst
python 2_graph_processing/make_graph_data.py --subst_type with -dir libraries_4_new/WT_sgCD_R1_antisense
python 2_graph_processing/make_graph_data.py --subst_type with -dir libraries_4_new/WT_sgCD_R1_splicing
python 2_graph_processing/make_graph_data.py --subst_type with -dir libraries_4_new/WT_sgCD_R2_antisense
python 2_graph_processing/make_graph_data.py --subst_type with -dir libraries_4_new/WT_sgCD_R2_splicing

:: 05_make_main_data_combined_withSubst
python 2_graph_processing/make_main_data_combined.py -i libraries_4_new/WT_sgCD_R1_antisense libraries_4_new/WT_sgCD_R1_splicing -o libraries_4_new/WT_sgCD_R1_antisense_splicing --subst_type without
python 2_graph_processing/make_main_data_combined.py -i libraries_4_new/WT_sgCD_R2_antisense libraries_4_new/WT_sgCD_R2_splicing -o libraries_4_new/WT_sgCD_R2_antisense_splicing --subst_type without

:: 05_make_main_data_combined_withoutSubst
python 2_graph_processing/make_main_data_combined.py -i libraries_4_new/WT_sgCD_R1_antisense libraries_4_new/WT_sgCD_R1_splicing -o libraries_4_new/WT_sgCD_R1_antisense_splicing --subst_type with
python 2_graph_processing/make_main_data_combined.py -i libraries_4_new/WT_sgCD_R2_antisense libraries_4_new/WT_sgCD_R2_splicing -o libraries_4_new/WT_sgCD_R2_antisense_splicing --subst_type with

:: 06_make_graph_data_combined_withoutSubst
python 2_graph_processing/make_graph_data.py --subst_type without -dir libraries_4_new/WT_sgCD_R1_antisense_splicing
python 2_graph_processing/make_graph_data.py --subst_type without -dir libraries_4_new/WT_sgCD_R2_antisense_splicing

:: 06_make_graph_data_combined_withSubst
python 2_graph_processing/make_graph_data.py --subst_type with -dir libraries_4_new/WT_sgCD_R1_antisense_splicing
python 2_graph_processing/make_graph_data.py --subst_type with -dir libraries_4_new/WT_sgCD_R2_antisense_splicing
