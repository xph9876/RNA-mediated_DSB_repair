python .\python\make_run_scripts.py `
  -i .\input\info `
  -o .\output `
  -p .\python `
  -f C:\Users\tchan\Code\DSB_Project\1DSB_alignment_analysis\input\fasta\1_trimmed `
  -n ..\NHEJ\data_1_filter_nhej `
  -r .\run
python .\python\make_search_tables.py -r .\input\ref_seq\1DSB_R1_sense.fa input\ref_seq\1DSB_R1_branch.fa -o .\input\search\