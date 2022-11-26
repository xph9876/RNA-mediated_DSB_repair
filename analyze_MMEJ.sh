#!/usr/bin/bash

# Author: Penghao
# Institute: Georgia Institute of Technology
# Description: Generate figures for RNA-mediated DSB repair study

# Parameters
# Path to GitHub repository
rep="."
scripts="${rep}/MMEJ/"
refseq="${rep}/refseq/"
# Path to trimmed reads folder. All the fastq files should be named as {sample name}_{R1/R2}.fq (eg. yjl217_R1.fq)
reads="/storage/home/hcoda1/0/pxu64/bio-storici/microhomology/trimmed_reads/"
# Path to output folder
output="./output/"
# Path to libinfo.tsv, should be 4 columns: Sample_name\tGenotype\tCelltype\tDSB
libinfo="libinfo.tsv"

# mkdir folder structures
if ! [ -e $output ]
then
    mkdir $output
fi

for aa in plots mh_sites freqs
do
    if ! [ -e $output/$aa ]
    then
        mkdir $output/$aa
    fi
done

for aa in barplots barplots_annot barplots_sum barplots_sum_annot boxplots_ee_ei_ratio_annot boxplots_ee_ei_ratio boxplots_wt_db_ratio scatterplot
do 
    if ! [ -e $output/plots/${aa} ]
    then
        mkdir $output/plots/$aa
    fi
done


# Step 1: Find MMEJ pairs from reference sequences
# find MMEJ pairs of each type
eval $scripts/find_microhomology.py $refseq/wt.fa 21 67 184 209 -o $output/mh_sites/exon_exon_wt_hg39.tsv &
eval $scripts/find_microhomology.py $refseq/db.fa 21 67 129 154 -o $output/mh_sites/exon_exon_db_hg39.tsv &
eval $scripts/find_microhomology.py $refseq/wt.fa 21 72 184 209 -o $output/mh_sites/exon_exon_wt_hg42.tsv &
eval $scripts/find_microhomology.py $refseq/db.fa 21 72 129 154 -o $output/mh_sites/exon_exon_db_hg42.tsv &
eval $scripts/find_microhomology.py $refseq/db.fa 21 67 129 154 -o $output/mh_sites/exon_exon_db_2dsb.tsv &
eval $scripts/find_microhomology.py $refseq/wt.fa 21 67 184 209 -o $output/mh_sites/exon_exon_wt_2dsb.tsv &
eval $scripts/find_microhomology.py $refseq/wt.fa 21 67 119 173 -o $output/mh_sites/exon_branch_wt_hg39.tsv &
eval $scripts/find_microhomology.py $refseq/wt.fa 119 173 184 209 -o $output/mh_sites/exon_branch_wt_hg42.tsv &
eval $scripts/find_microhomology.py $refseq/wt.fa 21 67 73 118 -o $output/mh_sites/exon_intron_wt_hg39.tsv &
eval $scripts/find_microhomology.py $refseq/db.fa 21 67 73 118 -o $output/mh_sites/exon_intron_db_hg39.tsv &
eval $scripts/find_microhomology.py $refseq/wt.fa 73 118 184 209 -o $output/mh_sites/exon_intron_wt_hg42.tsv &
eval $scripts/find_microhomology.py $refseq/db.fa 73 118 129 154 -o $output/mh_sites/exon_intron_db_hg42.tsv &
wait
# Merge all MMEJ sites
eval $scripts/merge_microhomology_sites.py $output/mh_sites/exon*.tsv -o $output/mh_sites/mmej.tsv
# Generate colored Excel file for analysis
eval $scripts/color_mmej.py $output/mh_sites/mmej.tsv -o $output/mmej_color.xlsx &
echo "MMEJ pairs generated!"


# Step 2: Calculate the frequency of each MMEJ pair
# Calculate the frequency for each reads, then sum up
for fq in $(eval ls $reads)
do
    eval $scripts/calc_mh_freq.py $output/mh_sites/mmej.tsv $libinfo $reads/$fq -o $output/freqs/${fq}.tsv &
done
wait
# merge
eval ls $output/freqs/*.tsv | head -n 1| xargs -I aa head -n 1 aa > $output/mmej_freqs.tsv
for tsv in $(eval ls $output/freqs | grep tsv)
do
    tail -n +2 $output/freqs/$tsv >> $output/mmej_freqs.tsv
done
# generate colored excel files
eval $scripts/color_mmej.py $output/mmej_freqs.tsv -o $output/mmej_freqs_color.xlsx &
echo "MMEJ frequency calculated!"

# Step 3: Generate figures
# Generate barplots for each mmej
eval $scripts/draw_barplot.py $output/mmej_freqs.tsv  -o $output/plots/barplots/MMEJ &
eval $scripts/draw_barplot.py $output/mmej_freqs.tsv  -o $output/plots/barplots_annot/MMEJ --annot &
eval $scripts/draw_barplot_sum.py $output/mmej_freqs.tsv  -o $output/plots/barplots_sum/MMEJ &
eval $scripts/draw_barplot_sum.py $output/mmej_freqs.tsv  -o $output/plots/barplots_sum_annot/MMEJ --annot &
wait

# Categorize the mmej barplots into overlapped data
cp $output/plots/barplots_annot $output/plots/barplots_grouped -r
eval $scripts/move_mmej_subfolder.py $output/mh_sites/mmej.tsv $output/plots/barplots_grouped

# move barplots
for aa in barplots barplots_annot
do
    mkdir $output/plots/${aa}/exon_exon $output/plots/${aa}/exon_branch $output/plots/${aa}/exon_intron
    mv $output/plots/${aa}/*EE* $output/plots/${aa}/exon_exon
    mv $output/plots/${aa}/*EB* $output/plots/${aa}/exon_branch
    mv $output/plots/${aa}/*EI* $output/plots/${aa}/exon_intron
done
echo "Barplots for MMEJ generated!"


# boxplots
eval $scripts/draw_ee_ei_ratio.py $output/mmej_freqs.tsv -o $output/plots/boxplots_ee_ei_ratio/MMEJ &
eval $scripts/draw_ee_ei_ratio.py $output/mmej_freqs.tsv -o $output/plots/boxplots_ee_ei_ratio_annot/MMEJ --annot &
eval $scripts/draw_wt_db_ratio.py $output/mmej_freqs.tsv -o $output/plots/boxplots_wt_db_ratio/MMEJ &
eval $scripts/draw_wt_db_ratio_KO_WT.py $output/mmej_freqs.tsv -o $output/plots/boxplots_wt_db_ratio/MMEJ &
wait
echo "Boxplots for cell type comparison generated!"

# scatterplot
for aa in ratios freqs ratios_left freqs_left ratios_right freqs_right
do
    if ! [ -e $output/plots/scatterplot/${aa} ]
    then
        mkdir $output/plots/scatterplot/${aa}
    fi
done

eval $scripts/draw_scatter.py $output/mmej_freqs.tsv -o $output/plots/scatterplot/MMEJ --remove_6bp
mv $output/plots/scatterplot/*ratios.png $output/plots/scatterplot/ratios 
mv $output/plots/scatterplot/*freqs.png $output/plots/scatterplot/freqs 
mv $output/plots/scatterplot/*ratios_left.png $output/plots/scatterplot/ratios_left 
mv $output/plots/scatterplot/*freqs_left.png $output/plots/scatterplot/freqs_left 
mv $output/plots/scatterplot/*ratios_right.png $output/plots/scatterplot/ratios_right 
mv $output/plots/scatterplot/*freqs_right.png $output/plots/scatterplot/freqs_right 

wait
echo "Scatter plots for the distance vs. frequency & ratios relationship generated!"
echo "Done!"
