bowtie2-build ref_seq/1DSB_R1_branch.fa data_bowtie2_build/1DSB_R1_branch
bowtie2-build ref_seq/1DSB_R1_cmv.fa data_bowtie2_build/1DSB_R1_cmv
bowtie2-build ref_seq/1DSB_R1_sense.fa data_bowtie2_build/1DSB_R1_sense
bowtie2-build ref_seq/1DSB_R2_branch.fa data_bowtie2_build/1DSB_R2_branch
bowtie2-build ref_seq/1DSB_R2_cmv.fa data_bowtie2_build/1DSB_R2_cmv
bowtie2-build ref_seq/1DSB_R2_sense.fa data_bowtie2_build/1DSB_R2_sense
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_branch data_fastq/yjl282_R1.fastq -S data_0_sam/yjl282_KO_sgA_R1_branch_noDSB.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_branch data_fastq/yjl296_R1.fastq -S data_0_sam/yjl296_KO_sgA_R1_branch.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_branch data_fastq/yjl297_R1.fastq -S data_0_sam/yjl297_KO_sgA_R1_branch.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_branch data_fastq/yjl298_R1.fastq -S data_0_sam/yjl298_KO_sgA_R1_branch.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_branch data_fastq/yjl299_R1.fastq -S data_0_sam/yjl299_KO_sgA_R1_branch.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl283_R1.fastq -S data_0_sam/yjl283_KO_sgA_R1_cmv_noDSB.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl300_R1.fastq -S data_0_sam/yjl300_KO_sgA_R1_cmv.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl301_R1.fastq -S data_0_sam/yjl301_KO_sgA_R1_cmv.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl302_R1.fastq -S data_0_sam/yjl302_KO_sgA_R1_cmv.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl303_R1.fastq -S data_0_sam/yjl303_KO_sgA_R1_cmv.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_sense data_fastq/yjl281_R1.fastq -S data_0_sam/yjl281_KO_sgA_R1_sense_noDSB.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_sense data_fastq/yjl292_R1.fastq -S data_0_sam/yjl292_KO_sgA_R1_sense.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_sense data_fastq/yjl293_R1.fastq -S data_0_sam/yjl293_KO_sgA_R1_sense.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_sense data_fastq/yjl294_R1.fastq -S data_0_sam/yjl294_KO_sgA_R1_sense.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_sense data_fastq/yjl295_R1.fastq -S data_0_sam/yjl295_KO_sgA_R1_sense.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_branch data_fastq/yjl245_R1.fastq -S data_0_sam/yjl245_WT_sgA_R1_branch_noDSB.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_branch data_fastq/yjl259_R1.fastq -S data_0_sam/yjl259_WT_sgA_R1_branch.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_branch data_fastq/yjl260_R1.fastq -S data_0_sam/yjl260_WT_sgA_R1_branch.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_branch data_fastq/yjl261_R1.fastq -S data_0_sam/yjl261_WT_sgA_R1_branch.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_branch data_fastq/yjl262_R1.fastq -S data_0_sam/yjl262_WT_sgA_R1_branch.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl246_R1.fastq -S data_0_sam/yjl246_WT_sgA_R1_cmv_noDSB.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl263_R1.fastq -S data_0_sam/yjl263_WT_sgA_R1_cmv.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl264_R1.fastq -S data_0_sam/yjl264_WT_sgA_R1_cmv.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl265_R1.fastq -S data_0_sam/yjl265_WT_sgA_R1_cmv.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_cmv data_fastq/yjl266_R1.fastq -S data_0_sam/yjl266_WT_sgA_R1_cmv.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_sense data_fastq/yjl244_R1.fastq -S data_0_sam/yjl244_WT_sgA_R1_sense_noDSB.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_sense data_fastq/yjl255_R1.fastq -S data_0_sam/yjl255_WT_sgA_R1_sense.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_sense data_fastq/yjl256_R1.fastq -S data_0_sam/yjl256_WT_sgA_R1_sense.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_sense data_fastq/yjl257_R1.fastq -S data_0_sam/yjl257_WT_sgA_R1_sense.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R1_sense data_fastq/yjl258_R1.fastq -S data_0_sam/yjl258_WT_sgA_R1_sense.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_branch data_fastq/yjl282_R2.fastq -S data_0_sam/yjl282_KO_sgB_R2_branch_noDSB.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_branch data_fastq/yjl308_R2.fastq -S data_0_sam/yjl308_KO_sgB_R2_branch.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_branch data_fastq/yjl309_R2.fastq -S data_0_sam/yjl309_KO_sgB_R2_branch.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_branch data_fastq/yjl310_R2.fastq -S data_0_sam/yjl310_KO_sgB_R2_branch.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_branch data_fastq/yjl311_R2.fastq -S data_0_sam/yjl311_KO_sgB_R2_branch.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl283_R2.fastq -S data_0_sam/yjl283_KO_sgB_R2_cmv_noDSB.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl312_R2.fastq -S data_0_sam/yjl312_KO_sgB_R2_cmv.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl313_R2.fastq -S data_0_sam/yjl313_KO_sgB_R2_cmv.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl314_R2.fastq -S data_0_sam/yjl314_KO_sgB_R2_cmv.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl315_R2.fastq -S data_0_sam/yjl315_KO_sgB_R2_cmv.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_sense data_fastq/yjl281_R2.fastq -S data_0_sam/yjl281_KO_sgB_R2_sense_noDSB.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_sense data_fastq/yjl304_R2.fastq -S data_0_sam/yjl304_KO_sgB_R2_sense.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_sense data_fastq/yjl305_R2.fastq -S data_0_sam/yjl305_KO_sgB_R2_sense.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_sense data_fastq/yjl306_R2.fastq -S data_0_sam/yjl306_KO_sgB_R2_sense.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_sense data_fastq/yjl307_R2.fastq -S data_0_sam/yjl307_KO_sgB_R2_sense.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_branch data_fastq/yjl245_R2.fastq -S data_0_sam/yjl245_WT_sgB_R2_branch_noDSB.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_branch data_fastq/yjl271_R2.fastq -S data_0_sam/yjl271_WT_sgB_R2_branch.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_branch data_fastq/yjl272_R2.fastq -S data_0_sam/yjl272_WT_sgB_R2_branch.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_branch data_fastq/yjl273_R2.fastq -S data_0_sam/yjl273_WT_sgB_R2_branch.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_branch data_fastq/yjl274_R2.fastq -S data_0_sam/yjl274_WT_sgB_R2_branch.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl246_R2.fastq -S data_0_sam/yjl246_WT_sgB_R2_cmv_noDSB.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl275_R2.fastq -S data_0_sam/yjl275_WT_sgB_R2_cmv.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl276_R2.fastq -S data_0_sam/yjl276_WT_sgB_R2_cmv.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl277_R2.fastq -S data_0_sam/yjl277_WT_sgB_R2_cmv.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_cmv data_fastq/yjl278_R2.fastq -S data_0_sam/yjl278_WT_sgB_R2_cmv.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_sense data_fastq/yjl244_R2.fastq -S data_0_sam/yjl244_WT_sgB_R2_sense_noDSB.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_sense data_fastq/yjl267_R2.fastq -S data_0_sam/yjl267_WT_sgB_R2_sense.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_sense data_fastq/yjl268_R2.fastq -S data_0_sam/yjl268_WT_sgB_R2_sense.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_sense data_fastq/yjl269_R2.fastq -S data_0_sam/yjl269_WT_sgB_R2_sense.sam
bowtie2 -x --no-hd --no-rc data_bowtie2_build/1DSB_R2_sense data_fastq/yjl270_R2.fastq -S data_0_sam/yjl270_WT_sgB_R2_sense.sam
