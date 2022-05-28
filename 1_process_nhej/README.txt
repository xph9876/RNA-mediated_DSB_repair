The example file is at:
  knot.math.usf.edu:/home/public/www_safe/jeon/SCMB/Scripts/yjl255_NHEJ_R1_aligned.sam
  https://knot.math.usf.edu/safe/users/jeon/SCMB/Scripts/yjl255_NHEJ_R1_aligned.sam

Example commands to run scripts:
python filter_nhej.py -fa 1DSB_ref.fa -sam test.sam -o output1.tsv -dsb "67" --min_length "130"
python filter_nhej.py -fa 1DSB_ref.fa -sam test.sam -o output2.tsv -dsb "67" --min_length "130"
python filter_nhej.py -fa 1DSB_ref.fa -sam test.sam -o output3.tsv -dsb "67" --min_length "130"
python filter_nhej.py -fa 1DSB_ref.fa -sam test.sam -o output4.tsv -dsb "67" --min_length "130"

python combine_repeats.py output1.tsv output2.tsv output3.tsv output4.tsv --total_reads 100 100 100 100 -o output_combined.tsv