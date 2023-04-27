mkdir -p output/flipped_intron
mkdir -p output/mmej
mkdir -p output/NHEJ
mkdir -p output/R-TDR
mkdir -p output/OX_WT
mkdir -p output/yeast
python test_ratios.py input/flipped_intron.tsv output/flipped_intron.xlsx -ca WT -cb KO -na Sense -nb BranchD -all
python test_ratios.py input/mmej.tsv output/mmej.xlsx -ca WT -cb KO -na Sense -nb BranchD -all
python test_ratios.py input/NHEJ.tsv output/NHEJ.xlsx -ca WT -cb KO -na Sense -nb BranchD -all
python test_ratios.py input/R-TDR.tsv output/R-TDR.xlsx -ca WT -cb KO -na Sense -nb BranchD -all
python test_ratios.py input/OX_WT.tsv output/OX_WT.xlsx -ca WT -cb OX -na Sense -nb BranchD -all
python test_ratios.py input/yeast.tsv output/yeast.xlsx -ca SPT3 -cb RHSPT3 -na Antisense -nb BranchD -all
