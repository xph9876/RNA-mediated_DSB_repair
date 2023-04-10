python test_ratios.py input/flipped_intron.tsv output/flipped_intron.xlsx -draw output/flipped_intron -ca WT -cb KO -na Sense -nb BranchD
python test_ratios.py input/mmej.tsv output/mmej.xlsx -draw output/mmej -ca WT -cb KO -na Sense -nb BranchD
python test_ratios.py input/NHEJ.tsv output/NHEJ.xlsx -draw output/NHEJ -ca WT -cb KO -na Sense -nb BranchD
python test_ratios.py input/R-TDR.tsv output/R-TDR.xlsx -draw output/R-TDR -ca WT -cb KO -na Sense -nb BranchD
python test_ratios.py input/OX_WT.tsv output/OX_WT.xlsx -draw output/OX_WT -ca WT -cb OX -na Sense -nb BranchD
python test_ratios.py input/yeast.tsv output/yeast.xlsx -draw output/yeast -ca SPT3 -cb RHSPT3 -na Antisense -nb BranchD
