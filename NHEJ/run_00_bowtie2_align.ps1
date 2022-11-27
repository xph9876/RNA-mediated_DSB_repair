bowtie2-build-s.exe ref_seq/2DSBanti_R1_antisense_old.fa data_bowtie2_build/2DSBanti_R1_antisense_old
bowtie2-build-s.exe ref_seq/2DSBanti_R2_antisense_old.fa data_bowtie2_build/2DSBanti_R2_antisense_old
bowtie2-build-s.exe ref_seq/2DSBanti_R1_splicing_old.fa data_bowtie2_build/2DSBanti_R1_splicing_old
bowtie2-build-s.exe ref_seq/2DSBanti_R2_splicing_old.fa data_bowtie2_build/2DSBanti_R2_splicing_old
bowtie2-build-s.exe ref_seq/2DSB_R1_sense.fa data_bowtie2_build/2DSB_R1_sense
bowtie2-build-s.exe ref_seq/2DSB_R2_sense.fa data_bowtie2_build/2DSB_R2_sense
bowtie2-build-s.exe ref_seq/2DSB_R1_branch.fa data_bowtie2_build/2DSB_R1_branch
bowtie2-build-s.exe ref_seq/2DSB_R2_branch.fa data_bowtie2_build/2DSB_R2_branch
bowtie2-build-s.exe ref_seq/2DSB_R1_cmv.fa data_bowtie2_build/2DSB_R1_cmv
bowtie2-build-s.exe ref_seq/2DSB_R2_cmv.fa data_bowtie2_build/2DSB_R2_cmv
bowtie2-build-s.exe ref_seq/1DSB_R1_sense.fa data_bowtie2_build/1DSB_R1_sense
bowtie2-build-s.exe ref_seq/1DSB_R2_sense.fa data_bowtie2_build/1DSB_R2_sense
bowtie2-build-s.exe ref_seq/1DSB_R1_branch.fa data_bowtie2_build/1DSB_R1_branch
bowtie2-build-s.exe ref_seq/1DSB_R2_branch.fa data_bowtie2_build/1DSB_R2_branch
bowtie2-build-s.exe ref_seq/1DSB_R1_cmv.fa data_bowtie2_build/1DSB_R1_cmv
bowtie2-build-s.exe ref_seq/1DSB_R2_cmv.fa data_bowtie2_build/1DSB_R2_cmv
bowtie2-build-s.exe ref_seq/2DSBanti_R1_antisense_new.fa data_bowtie2_build/2DSBanti_R1_antisense_new
bowtie2-build-s.exe ref_seq/2DSBanti_R2_antisense_new.fa data_bowtie2_build/2DSBanti_R2_antisense_new
bowtie2-build-s.exe ref_seq/2DSBanti_R1_splicing_new.fa data_bowtie2_build/2DSBanti_R1_splicing_new
bowtie2-build-s.exe ref_seq/2DSBanti_R2_splicing_new.fa data_bowtie2_build/2DSBanti_R2_splicing_new
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R1_antisense_old data_fastq/yjl89_R1.fastq -S data_0_sam/yjl89_WT_sgCD_R1_antisense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R2_antisense_old data_fastq/yjl89_R2.fastq -S data_0_sam/yjl89_WT_sgCD_R2_antisense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R1_antisense_old data_fastq/yjl90_R1.fastq -S data_0_sam/yjl90_WT_sgCD_R1_antisense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R2_antisense_old data_fastq/yjl90_R2.fastq -S data_0_sam/yjl90_WT_sgCD_R2_antisense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R1_antisense_old data_fastq/yjl91_R1.fastq -S data_0_sam/yjl91_WT_sgCD_R1_antisense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R2_antisense_old data_fastq/yjl91_R2.fastq -S data_0_sam/yjl91_WT_sgCD_R2_antisense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R1_antisense_old data_fastq/yjl92_R1.fastq -S data_0_sam/yjl92_WT_sgCD_R1_antisense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R2_antisense_old data_fastq/yjl92_R2.fastq -S data_0_sam/yjl92_WT_sgCD_R2_antisense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R1_splicing_old data_fastq/yjl93_R1.fastq -S data_0_sam/yjl93_WT_sgCD_R1_splicing.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R2_splicing_old data_fastq/yjl93_R2.fastq -S data_0_sam/yjl93_WT_sgCD_R2_splicing.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R1_splicing_old data_fastq/yjl94_R1.fastq -S data_0_sam/yjl94_WT_sgCD_R1_splicing.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R2_splicing_old data_fastq/yjl94_R2.fastq -S data_0_sam/yjl94_WT_sgCD_R2_splicing.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R1_splicing_old data_fastq/yjl95_R1.fastq -S data_0_sam/yjl95_WT_sgCD_R1_splicing.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R2_splicing_old data_fastq/yjl95_R2.fastq -S data_0_sam/yjl95_WT_sgCD_R2_splicing.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R1_splicing_old data_fastq/yjl96_R1.fastq -S data_0_sam/yjl96_WT_sgCD_R1_splicing.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R2_splicing_old data_fastq/yjl96_R2.fastq -S data_0_sam/yjl96_WT_sgCD_R2_splicing.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_sense data_fastq/yjl217_R1.fastq -S data_0_sam/yjl217_WT_sgAB_R1_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_sense data_fastq/yjl217_R2.fastq -S data_0_sam/yjl217_WT_sgAB_R2_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_sense data_fastq/yjl218_R1.fastq -S data_0_sam/yjl218_WT_sgAB_R1_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_sense data_fastq/yjl218_R2.fastq -S data_0_sam/yjl218_WT_sgAB_R2_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_sense data_fastq/yjl219_R1.fastq -S data_0_sam/yjl219_WT_sgAB_R1_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_sense data_fastq/yjl219_R2.fastq -S data_0_sam/yjl219_WT_sgAB_R2_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_sense data_fastq/yjl220_R1.fastq -S data_0_sam/yjl220_WT_sgAB_R1_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_sense data_fastq/yjl220_R2.fastq -S data_0_sam/yjl220_WT_sgAB_R2_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_branch data_fastq/yjl221_R1.fastq -S data_0_sam/yjl221_WT_sgAB_R1_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_branch data_fastq/yjl221_R2.fastq -S data_0_sam/yjl221_WT_sgAB_R2_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_branch data_fastq/yjl222_R1.fastq -S data_0_sam/yjl222_WT_sgAB_R1_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_branch data_fastq/yjl222_R2.fastq -S data_0_sam/yjl222_WT_sgAB_R2_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_branch data_fastq/yjl223_R1.fastq -S data_0_sam/yjl223_WT_sgAB_R1_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_branch data_fastq/yjl223_R2.fastq -S data_0_sam/yjl223_WT_sgAB_R2_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_branch data_fastq/yjl224_R1.fastq -S data_0_sam/yjl224_WT_sgAB_R1_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_branch data_fastq/yjl224_R2.fastq -S data_0_sam/yjl224_WT_sgAB_R2_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_cmv data_fastq/yjl225_R1.fastq -S data_0_sam/yjl225_WT_sgAB_R1_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_cmv data_fastq/yjl225_R2.fastq -S data_0_sam/yjl225_WT_sgAB_R2_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_cmv data_fastq/yjl226_R1.fastq -S data_0_sam/yjl226_WT_sgAB_R1_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_cmv data_fastq/yjl226_R2.fastq -S data_0_sam/yjl226_WT_sgAB_R2_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_cmv data_fastq/yjl227_R1.fastq -S data_0_sam/yjl227_WT_sgAB_R1_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_cmv data_fastq/yjl227_R2.fastq -S data_0_sam/yjl227_WT_sgAB_R2_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_cmv data_fastq/yjl228_R1.fastq -S data_0_sam/yjl228_WT_sgAB_R1_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_cmv data_fastq/yjl228_R2.fastq -S data_0_sam/yjl228_WT_sgAB_R2_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_sense data_fastq/yjl229_R1.fastq -S data_0_sam/yjl229_KO_sgAB_R1_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_sense data_fastq/yjl229_R2.fastq -S data_0_sam/yjl229_KO_sgAB_R2_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_sense data_fastq/yjl230_R1.fastq -S data_0_sam/yjl230_KO_sgAB_R1_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_sense data_fastq/yjl230_R2.fastq -S data_0_sam/yjl230_KO_sgAB_R2_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_sense data_fastq/yjl231_R1.fastq -S data_0_sam/yjl231_KO_sgAB_R1_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_sense data_fastq/yjl231_R2.fastq -S data_0_sam/yjl231_KO_sgAB_R2_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_sense data_fastq/yjl232_R1.fastq -S data_0_sam/yjl232_KO_sgAB_R1_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_sense data_fastq/yjl232_R2.fastq -S data_0_sam/yjl232_KO_sgAB_R2_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_branch data_fastq/yjl233_R1.fastq -S data_0_sam/yjl233_KO_sgAB_R1_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_branch data_fastq/yjl233_R2.fastq -S data_0_sam/yjl233_KO_sgAB_R2_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_branch data_fastq/yjl234_R1.fastq -S data_0_sam/yjl234_KO_sgAB_R1_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_branch data_fastq/yjl234_R2.fastq -S data_0_sam/yjl234_KO_sgAB_R2_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_branch data_fastq/yjl235_R1.fastq -S data_0_sam/yjl235_KO_sgAB_R1_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_branch data_fastq/yjl235_R2.fastq -S data_0_sam/yjl235_KO_sgAB_R2_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_branch data_fastq/yjl236_R1.fastq -S data_0_sam/yjl236_KO_sgAB_R1_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_branch data_fastq/yjl236_R2.fastq -S data_0_sam/yjl236_KO_sgAB_R2_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_cmv data_fastq/yjl237_R1.fastq -S data_0_sam/yjl237_KO_sgAB_R1_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_cmv data_fastq/yjl237_R2.fastq -S data_0_sam/yjl237_KO_sgAB_R2_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_cmv data_fastq/yjl238_R1.fastq -S data_0_sam/yjl238_KO_sgAB_R1_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_cmv data_fastq/yjl238_R2.fastq -S data_0_sam/yjl238_KO_sgAB_R2_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_cmv data_fastq/yjl239_R1.fastq -S data_0_sam/yjl239_KO_sgAB_R1_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_cmv data_fastq/yjl239_R2.fastq -S data_0_sam/yjl239_KO_sgAB_R2_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R1_cmv data_fastq/yjl240_R1.fastq -S data_0_sam/yjl240_KO_sgAB_R1_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSB_R2_cmv data_fastq/yjl240_R2.fastq -S data_0_sam/yjl240_KO_sgAB_R2_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_sense data_fastq/yjl244_R1.fastq -S data_0_sam/yjl244_WT_sgA_R1_sense_noDSB.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_sense data_fastq/yjl244_R2.fastq -S data_0_sam/yjl244_WT_sgB_R2_sense_noDSB.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_branch data_fastq/yjl245_R1.fastq -S data_0_sam/yjl245_WT_sgA_R1_branch_noDSB.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_branch data_fastq/yjl245_R2.fastq -S data_0_sam/yjl245_WT_sgB_R2_branch_noDSB.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl246_R1.fastq -S data_0_sam/yjl246_WT_sgA_R1_cmv_noDSB.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl246_R2.fastq -S data_0_sam/yjl246_WT_sgB_R2_cmv_noDSB.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_sense data_fastq/yjl255_R1.fastq -S data_0_sam/yjl255_WT_sgA_R1_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_sense data_fastq/yjl255_R1.fastq -S data_0_sam/yjl255_WT_sgA_R1_sense_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_sense data_fastq/yjl256_R1.fastq -S data_0_sam/yjl256_WT_sgA_R1_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_sense data_fastq/yjl256_R1.fastq -S data_0_sam/yjl256_WT_sgA_R1_sense_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_sense data_fastq/yjl257_R1.fastq -S data_0_sam/yjl257_WT_sgA_R1_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_sense data_fastq/yjl257_R1.fastq -S data_0_sam/yjl257_WT_sgA_R1_sense_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_sense data_fastq/yjl258_R1.fastq -S data_0_sam/yjl258_WT_sgA_R1_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_sense data_fastq/yjl258_R1.fastq -S data_0_sam/yjl258_WT_sgA_R1_sense_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_branch data_fastq/yjl259_R1.fastq -S data_0_sam/yjl259_WT_sgA_R1_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_branch data_fastq/yjl259_R1.fastq -S data_0_sam/yjl259_WT_sgA_R1_branch_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_branch data_fastq/yjl260_R1.fastq -S data_0_sam/yjl260_WT_sgA_R1_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_branch data_fastq/yjl260_R1.fastq -S data_0_sam/yjl260_WT_sgA_R1_branch_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_branch data_fastq/yjl261_R1.fastq -S data_0_sam/yjl261_WT_sgA_R1_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_branch data_fastq/yjl261_R1.fastq -S data_0_sam/yjl261_WT_sgA_R1_branch_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_branch data_fastq/yjl262_R1.fastq -S data_0_sam/yjl262_WT_sgA_R1_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_branch data_fastq/yjl262_R1.fastq -S data_0_sam/yjl262_WT_sgA_R1_branch_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl263_R1.fastq -S data_0_sam/yjl263_WT_sgA_R1_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl263_R1.fastq -S data_0_sam/yjl263_WT_sgA_R1_cmv_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl264_R1.fastq -S data_0_sam/yjl264_WT_sgA_R1_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl264_R1.fastq -S data_0_sam/yjl264_WT_sgA_R1_cmv_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl265_R1.fastq -S data_0_sam/yjl265_WT_sgA_R1_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl265_R1.fastq -S data_0_sam/yjl265_WT_sgA_R1_cmv_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl266_R1.fastq -S data_0_sam/yjl266_WT_sgA_R1_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl266_R1.fastq -S data_0_sam/yjl266_WT_sgA_R1_cmv_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_sense data_fastq/yjl267_R2.fastq -S data_0_sam/yjl267_WT_sgB_R2_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_sense data_fastq/yjl267_R2.fastq -S data_0_sam/yjl267_WT_sgB_R2_sense_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_sense data_fastq/yjl268_R2.fastq -S data_0_sam/yjl268_WT_sgB_R2_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_sense data_fastq/yjl268_R2.fastq -S data_0_sam/yjl268_WT_sgB_R2_sense_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_sense data_fastq/yjl269_R2.fastq -S data_0_sam/yjl269_WT_sgB_R2_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_sense data_fastq/yjl269_R2.fastq -S data_0_sam/yjl269_WT_sgB_R2_sense_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_sense data_fastq/yjl270_R2.fastq -S data_0_sam/yjl270_WT_sgB_R2_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_sense data_fastq/yjl270_R2.fastq -S data_0_sam/yjl270_WT_sgB_R2_sense_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_branch data_fastq/yjl271_R2.fastq -S data_0_sam/yjl271_WT_sgB_R2_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_branch data_fastq/yjl271_R2.fastq -S data_0_sam/yjl271_WT_sgB_R2_branch_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_branch data_fastq/yjl272_R2.fastq -S data_0_sam/yjl272_WT_sgB_R2_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_branch data_fastq/yjl272_R2.fastq -S data_0_sam/yjl272_WT_sgB_R2_branch_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_branch data_fastq/yjl273_R2.fastq -S data_0_sam/yjl273_WT_sgB_R2_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_branch data_fastq/yjl273_R2.fastq -S data_0_sam/yjl273_WT_sgB_R2_branch_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_branch data_fastq/yjl274_R2.fastq -S data_0_sam/yjl274_WT_sgB_R2_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_branch data_fastq/yjl274_R2.fastq -S data_0_sam/yjl274_WT_sgB_R2_branch_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl275_R2.fastq -S data_0_sam/yjl275_WT_sgB_R2_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl275_R2.fastq -S data_0_sam/yjl275_WT_sgB_R2_cmv_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl276_R2.fastq -S data_0_sam/yjl276_WT_sgB_R2_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl276_R2.fastq -S data_0_sam/yjl276_WT_sgB_R2_cmv_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl277_R2.fastq -S data_0_sam/yjl277_WT_sgB_R2_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl277_R2.fastq -S data_0_sam/yjl277_WT_sgB_R2_cmv_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl278_R2.fastq -S data_0_sam/yjl278_WT_sgB_R2_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl278_R2.fastq -S data_0_sam/yjl278_WT_sgB_R2_cmv_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_sense data_fastq/yjl281_R1.fastq -S data_0_sam/yjl281_KO_sgA_R1_sense_noDSB.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_sense data_fastq/yjl281_R2.fastq -S data_0_sam/yjl281_KO_sgB_R2_sense_noDSB.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_branch data_fastq/yjl282_R1.fastq -S data_0_sam/yjl282_KO_sgA_R1_branch_noDSB.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_branch data_fastq/yjl282_R2.fastq -S data_0_sam/yjl282_KO_sgB_R2_branch_noDSB.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl283_R1.fastq -S data_0_sam/yjl283_KO_sgA_R1_cmv_noDSB.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl283_R2.fastq -S data_0_sam/yjl283_KO_sgB_R2_cmv_noDSB.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_sense data_fastq/yjl292_R1.fastq -S data_0_sam/yjl292_KO_sgA_R1_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_sense data_fastq/yjl292_R1.fastq -S data_0_sam/yjl292_KO_sgA_R1_sense_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_sense data_fastq/yjl293_R1.fastq -S data_0_sam/yjl293_KO_sgA_R1_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_sense data_fastq/yjl293_R1.fastq -S data_0_sam/yjl293_KO_sgA_R1_sense_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_sense data_fastq/yjl294_R1.fastq -S data_0_sam/yjl294_KO_sgA_R1_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_sense data_fastq/yjl294_R1.fastq -S data_0_sam/yjl294_KO_sgA_R1_sense_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_sense data_fastq/yjl295_R1.fastq -S data_0_sam/yjl295_KO_sgA_R1_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_sense data_fastq/yjl295_R1.fastq -S data_0_sam/yjl295_KO_sgA_R1_sense_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_branch data_fastq/yjl296_R1.fastq -S data_0_sam/yjl296_KO_sgA_R1_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_branch data_fastq/yjl296_R1.fastq -S data_0_sam/yjl296_KO_sgA_R1_branch_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_branch data_fastq/yjl297_R1.fastq -S data_0_sam/yjl297_KO_sgA_R1_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_branch data_fastq/yjl297_R1.fastq -S data_0_sam/yjl297_KO_sgA_R1_branch_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_branch data_fastq/yjl298_R1.fastq -S data_0_sam/yjl298_KO_sgA_R1_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_branch data_fastq/yjl298_R1.fastq -S data_0_sam/yjl298_KO_sgA_R1_branch_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_branch data_fastq/yjl299_R1.fastq -S data_0_sam/yjl299_KO_sgA_R1_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_branch data_fastq/yjl299_R1.fastq -S data_0_sam/yjl299_KO_sgA_R1_branch_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl300_R1.fastq -S data_0_sam/yjl300_KO_sgA_R1_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl300_R1.fastq -S data_0_sam/yjl300_KO_sgA_R1_cmv_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl301_R1.fastq -S data_0_sam/yjl301_KO_sgA_R1_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl301_R1.fastq -S data_0_sam/yjl301_KO_sgA_R1_cmv_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl302_R1.fastq -S data_0_sam/yjl302_KO_sgA_R1_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl302_R1.fastq -S data_0_sam/yjl302_KO_sgA_R1_cmv_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl303_R1.fastq -S data_0_sam/yjl303_KO_sgA_R1_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl303_R1.fastq -S data_0_sam/yjl303_KO_sgA_R1_cmv_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_sense data_fastq/yjl304_R2.fastq -S data_0_sam/yjl304_KO_sgB_R2_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_sense data_fastq/yjl304_R2.fastq -S data_0_sam/yjl304_KO_sgB_R2_sense_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_sense data_fastq/yjl305_R2.fastq -S data_0_sam/yjl305_KO_sgB_R2_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_sense data_fastq/yjl305_R2.fastq -S data_0_sam/yjl305_KO_sgB_R2_sense_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_sense data_fastq/yjl306_R2.fastq -S data_0_sam/yjl306_KO_sgB_R2_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_sense data_fastq/yjl306_R2.fastq -S data_0_sam/yjl306_KO_sgB_R2_sense_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_sense data_fastq/yjl307_R2.fastq -S data_0_sam/yjl307_KO_sgB_R2_sense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_sense data_fastq/yjl307_R2.fastq -S data_0_sam/yjl307_KO_sgB_R2_sense_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_branch data_fastq/yjl308_R2.fastq -S data_0_sam/yjl308_KO_sgB_R2_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_branch data_fastq/yjl308_R2.fastq -S data_0_sam/yjl308_KO_sgB_R2_branch_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_branch data_fastq/yjl309_R2.fastq -S data_0_sam/yjl309_KO_sgB_R2_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_branch data_fastq/yjl309_R2.fastq -S data_0_sam/yjl309_KO_sgB_R2_branch_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_branch data_fastq/yjl310_R2.fastq -S data_0_sam/yjl310_KO_sgB_R2_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_branch data_fastq/yjl310_R2.fastq -S data_0_sam/yjl310_KO_sgB_R2_branch_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_branch data_fastq/yjl311_R2.fastq -S data_0_sam/yjl311_KO_sgB_R2_branch.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_branch data_fastq/yjl311_R2.fastq -S data_0_sam/yjl311_KO_sgB_R2_branch_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl312_R2.fastq -S data_0_sam/yjl312_KO_sgB_R2_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl312_R2.fastq -S data_0_sam/yjl312_KO_sgB_R2_cmv_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl313_R2.fastq -S data_0_sam/yjl313_KO_sgB_R2_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl313_R2.fastq -S data_0_sam/yjl313_KO_sgB_R2_cmv_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl314_R2.fastq -S data_0_sam/yjl314_KO_sgB_R2_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl314_R2.fastq -S data_0_sam/yjl314_KO_sgB_R2_cmv_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl315_R2.fastq -S data_0_sam/yjl315_KO_sgB_R2_cmv.sam
bowtie2-align-s.exe -x  data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl315_R2.fastq -S data_0_sam/yjl315_KO_sgB_R2_cmv_30bpDown.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R1_antisense_new data_fastq/yjl349_R1.fastq -S data_0_sam/yjl349_WT_sgCD_R1_antisense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R2_antisense_new data_fastq/yjl349_R2.fastq -S data_0_sam/yjl349_WT_sgCD_R2_antisense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R1_antisense_new data_fastq/yjl350_R1.fastq -S data_0_sam/yjl350_WT_sgCD_R1_antisense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R2_antisense_new data_fastq/yjl350_R2.fastq -S data_0_sam/yjl350_WT_sgCD_R2_antisense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R1_antisense_new data_fastq/yjl351_R1.fastq -S data_0_sam/yjl351_WT_sgCD_R1_antisense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R2_antisense_new data_fastq/yjl351_R2.fastq -S data_0_sam/yjl351_WT_sgCD_R2_antisense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R1_antisense_new data_fastq/yjl352_R1.fastq -S data_0_sam/yjl352_WT_sgCD_R1_antisense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R2_antisense_new data_fastq/yjl352_R2.fastq -S data_0_sam/yjl352_WT_sgCD_R2_antisense.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R1_splicing_new data_fastq/yjl353_R1.fastq -S data_0_sam/yjl353_WT_sgCD_R1_splicing.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R2_splicing_new data_fastq/yjl353_R2.fastq -S data_0_sam/yjl353_WT_sgCD_R2_splicing.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R1_splicing_new data_fastq/yjl354_R1.fastq -S data_0_sam/yjl354_WT_sgCD_R1_splicing.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R2_splicing_new data_fastq/yjl354_R2.fastq -S data_0_sam/yjl354_WT_sgCD_R2_splicing.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R1_splicing_new data_fastq/yjl355_R1.fastq -S data_0_sam/yjl355_WT_sgCD_R1_splicing.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R2_splicing_new data_fastq/yjl355_R2.fastq -S data_0_sam/yjl355_WT_sgCD_R2_splicing.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R1_splicing_new data_fastq/yjl356_R1.fastq -S data_0_sam/yjl356_WT_sgCD_R1_splicing.sam
bowtie2-align-s.exe -x  data_bowtie2_build/2DSBanti_R2_splicing_new data_fastq/yjl356_R2.fastq -S data_0_sam/yjl356_WT_sgCD_R2_splicing.sam
