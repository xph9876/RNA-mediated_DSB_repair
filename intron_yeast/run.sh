# Set the directories to the correct path
R1_fastq_file_directory=../R1_fastq
R2_fastq_file_directory=../R2_fastq
python freq_intron_with_seq_F.py $R1_fastq_file_directory --exon CCAGAGCATGTATCATATGGTCCAGAAACCCTATACCTGTGTGGACGTTAATCACTTGCGATTG -o output.tsv
python freq_intron_with_seq_R.py $R2_fastq_file_directory --exon CCAAGCATTCCGGCTGGTCGCTAATCGTTGAGTGCATTGGTGACTTACACATAGACGACCATCACACCACTGAAGACTGCGGGATTGCTCTC -o output.tsv