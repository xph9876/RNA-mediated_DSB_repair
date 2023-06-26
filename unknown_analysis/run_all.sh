python ./python/make_search_tables.py \
  -rs ./input/ref_seq/1DSB_R1_sense.fa \
  -rb input/ref_seq/1DSB_R1_branch.fa \
  -o ./input/search/ -b 119 173

python ./python/make_run_scripts.py \
  -i ./input/info \
  -o ./output \
  -p ./python \
  -f C:/Users/tchan/Code/DSB_Project/1DSB_alignment_analysis/input/fasta/1_trimmed \
  -n ../NHEJ/data_1_filter_nhej \
  -r ./run

./run/run_analyze_alignments.sh

./run/run_nhej_mmej.sh

./run/run_full_output.sh

./run/run_summary_output.sh

./run/run_compare_libraries.sh

./run/run_mean_tables.sh
