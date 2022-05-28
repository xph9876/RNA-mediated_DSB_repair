The example file is at:
  knot.math.usf.edu:/home/public/www_safe/jeon/SCMB/Scripts/yjl255_NHEJ_R1_aligned.sam
  https://knot.math.usf.edu/safe/users/jeon/SCMB/Scripts/yjl255_NHEJ_R1_aligned.sam

* Instead of NM - XM to compute in/dels in Bowtie2 alignment, I used XG since they are the same
* I parsed the optional tags in the Bowtie2 output since according to the spec the order is not specified.
* I checked the (FLAG & 4) to filter out reads that did not align.
* I filtered out reads that did not align with position 1 of reference.
* I checked that the min length was achieved before all the other processing

Example of running script:
python mut_middle_indel.py -fa 1DSB_ref.fa -sam test.sam -o output.tsv -dsb "67" --min_length "140"