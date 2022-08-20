:: 02_combine_repeats
python 1_process_nhej/combine_repeats.py -i libraries_2/yjl89_WT_sgCD_R1_antisense.tsv libraries_2/yjl90_WT_sgCD_R1_antisense.tsv libraries_2/yjl91_WT_sgCD_R1_antisense.tsv libraries_2/yjl92_WT_sgCD_R1_antisense.tsv --total_reads 7813781 7486736 8114856 6586577 -o libraries_3_old/WT_sgCD_R1_antisense.tsv --quiet 
python 1_process_nhej/combine_repeats.py -i libraries_2/yjl93_WT_sgCD_R1_splicing.tsv libraries_2/yjl94_WT_sgCD_R1_splicing.tsv libraries_2/yjl95_WT_sgCD_R1_splicing.tsv libraries_2/yjl96_WT_sgCD_R1_splicing.tsv --total_reads 9334492 9142001 8497185 8346242 -o libraries_3_old/WT_sgCD_R1_splicing.tsv --quiet 
python 1_process_nhej/combine_repeats.py -i libraries_2/yjl89_WT_sgCD_R2_antisense.tsv libraries_2/yjl90_WT_sgCD_R2_antisense.tsv libraries_2/yjl91_WT_sgCD_R2_antisense.tsv libraries_2/yjl92_WT_sgCD_R2_antisense.tsv --total_reads 7813781 7486736 8114856 6586577 -o libraries_3_old/WT_sgCD_R2_antisense.tsv --quiet 
python 1_process_nhej/combine_repeats.py -i libraries_2/yjl93_WT_sgCD_R2_splicing.tsv libraries_2/yjl94_WT_sgCD_R2_splicing.tsv libraries_2/yjl95_WT_sgCD_R2_splicing.tsv libraries_2/yjl96_WT_sgCD_R2_splicing.tsv --total_reads 9334492 9142001 8497185 8346242 -o libraries_3_old/WT_sgCD_R2_splicing.tsv --quiet 

:: 03_make_main_data_withoutSubst
python 2_graph_processing/make_main_data.py -i libraries_3_old/WT_sgCD_R1_antisense.tsv -o libraries_4_old/WT_sgCD_R1_antisense -ref ref_seq/2DSBanti_R1_antisense_old.fa -dsb 50 --dsb_type 2a --strand R1 --hguide CD --cell WT --treatment antisense --subst_type without --control none -f 0
python 2_graph_processing/make_main_data.py -i libraries_3_old/WT_sgCD_R1_splicing.tsv -o libraries_4_old/WT_sgCD_R1_splicing -ref ref_seq/2DSBanti_R1_splicing_old.fa -dsb 50 --dsb_type 2a --strand R1 --hguide CD --cell WT --treatment splicing --subst_type without --control none -f 0
python 2_graph_processing/make_main_data.py -i libraries_3_old/WT_sgCD_R2_antisense.tsv -o libraries_4_old/WT_sgCD_R2_antisense -ref ref_seq/2DSBanti_R2_antisense_old.fa -dsb 47 --dsb_type 2a --strand R2 --hguide CD --cell WT --treatment antisense --subst_type without --control none -f 0
python 2_graph_processing/make_main_data.py -i libraries_3_old/WT_sgCD_R2_splicing.tsv -o libraries_4_old/WT_sgCD_R2_splicing -ref ref_seq/2DSBanti_R2_splicing_old.fa -dsb 47 --dsb_type 2a --strand R2 --hguide CD --cell WT --treatment splicing --subst_type without --control none -f 0

:: 03_make_main_data_withSubst
python 2_graph_processing/make_main_data.py -i libraries_3_old/WT_sgCD_R1_antisense.tsv -o libraries_4_old/WT_sgCD_R1_antisense -ref ref_seq/2DSBanti_R1_antisense_old.fa -dsb 50 --dsb_type 2a --strand R1 --hguide CD --cell WT --treatment antisense --subst_type with --control none -f 0
python 2_graph_processing/make_main_data.py -i libraries_3_old/WT_sgCD_R1_splicing.tsv -o libraries_4_old/WT_sgCD_R1_splicing -ref ref_seq/2DSBanti_R1_splicing_old.fa -dsb 50 --dsb_type 2a --strand R1 --hguide CD --cell WT --treatment splicing --subst_type with --control none -f 0
python 2_graph_processing/make_main_data.py -i libraries_3_old/WT_sgCD_R2_antisense.tsv -o libraries_4_old/WT_sgCD_R2_antisense -ref ref_seq/2DSBanti_R2_antisense_old.fa -dsb 47 --dsb_type 2a --strand R2 --hguide CD --cell WT --treatment antisense --subst_type with --control none -f 0
python 2_graph_processing/make_main_data.py -i libraries_3_old/WT_sgCD_R2_splicing.tsv -o libraries_4_old/WT_sgCD_R2_splicing -ref ref_seq/2DSBanti_R2_splicing_old.fa -dsb 47 --dsb_type 2a --strand R2 --hguide CD --cell WT --treatment splicing --subst_type with --control none -f 0

:: 04_make_graph_data_withoutSubst
python 2_graph_processing/make_graph_data.py --subst_type without -dir libraries_4_old/WT_sgCD_R1_antisense
python 2_graph_processing/make_graph_data.py --subst_type without -dir libraries_4_old/WT_sgCD_R1_splicing
python 2_graph_processing/make_graph_data.py --subst_type without -dir libraries_4_old/WT_sgCD_R2_antisense
python 2_graph_processing/make_graph_data.py --subst_type without -dir libraries_4_old/WT_sgCD_R2_splicing

:: 04_make_graph_data_withSubst
python 2_graph_processing/make_graph_data.py --subst_type with -dir libraries_4_old/WT_sgCD_R1_antisense
python 2_graph_processing/make_graph_data.py --subst_type with -dir libraries_4_old/WT_sgCD_R1_splicing
python 2_graph_processing/make_graph_data.py --subst_type with -dir libraries_4_old/WT_sgCD_R2_antisense
python 2_graph_processing/make_graph_data.py --subst_type with -dir libraries_4_old/WT_sgCD_R2_splicing

:: 05_make_main_data_combined_withSubst
python 2_graph_processing/make_main_data_combined.py -i libraries_4_old/WT_sgCD_R1_antisense libraries_4_old/WT_sgCD_R1_splicing -o libraries_4_old/WT_sgCD_R1_antisense_splicing --subst_type without
python 2_graph_processing/make_main_data_combined.py -i libraries_4_old/WT_sgCD_R2_antisense libraries_4_old/WT_sgCD_R2_splicing -o libraries_4_old/WT_sgCD_R2_antisense_splicing --subst_type without

:: 05_make_main_data_combined_withoutSubst
python 2_graph_processing/make_main_data_combined.py -i libraries_4_old/WT_sgCD_R1_antisense libraries_4_old/WT_sgCD_R1_splicing -o libraries_4_old/WT_sgCD_R1_antisense_splicing --subst_type with
python 2_graph_processing/make_main_data_combined.py -i libraries_4_old/WT_sgCD_R2_antisense libraries_4_old/WT_sgCD_R2_splicing -o libraries_4_old/WT_sgCD_R2_antisense_splicing --subst_type with

:: 06_make_graph_data_combined_withoutSubst
python 2_graph_processing/make_graph_data.py --subst_type without -dir libraries_4_old/WT_sgCD_R1_antisense_splicing
python 2_graph_processing/make_graph_data.py --subst_type without -dir libraries_4_old/WT_sgCD_R2_antisense_splicing

:: 06_make_graph_data_combined_withSubst
python 2_graph_processing/make_graph_data.py --subst_type with -dir libraries_4_old/WT_sgCD_R1_antisense_splicing
python 2_graph_processing/make_graph_data.py --subst_type with -dir libraries_4_old/WT_sgCD_R2_antisense_splicing
