The example file is at:
  knot.math.usf.edu:/home/public/www_safe/jeon/SCMB/Scripts/yjl255_NHEJ_R1_aligned.sam
  https://knot.math.usf.edu/safe/users/jeon/SCMB/Scripts/yjl255_NHEJ_R1_aligned.sam

* Please see algorithm.pdf for the steps in the filtering (not 100% sure if it is the sam as the original algorithm)
* Instead of NM - XM to compute in/dels in Bowtie2 alignment, I used XG since they are the same
* I parsed the optional tags in the Bowtie2 output since according to the spec the order is not specified.
* I checked the (FLAG & 4) to filter out reads that did not align.
* I filtered out reads that did not align with position 1 of reference.