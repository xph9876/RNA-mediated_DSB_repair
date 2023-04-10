python bootstrap.py input/flipped_intron.tsv output_diff/flipped_intron.xlsx -draw output_diff/flipped_intron -ca WT -cb KO -na Sense -nb BranchD -stat diff
python bootstrap.py input/mmej.tsv output_diff/mmej.xlsx -draw output_diff/mmej -ca WT -cb KO -na Sense -nb BranchD -stat diff
python bootstrap.py input/NHEJ.tsv output_diff/NHEJ.xlsx -draw output_diff/NHEJ -ca WT -cb KO -na Sense -nb BranchD -stat diff
python bootstrap.py input/R-TDR.tsv output_diff/R-TDR.xlsx -draw output_diff/R-TDR -ca WT -cb KO -na Sense -nb BranchD -stat diff
python bootstrap.py input/OX_WT.tsv output_diff/OX_WT.xlsx -draw output_diff/OX_WT -ca WT -cb OX -na Sense -nb BranchD -stat diff
python bootstrap.py input/yeast.tsv output_diff/yeast.xlsx -draw output_diff/yeast -ca SPT3 -cb RHSPT3 -na Sense -nb BranchD -stat diff

python bootstrap.py input/flipped_intron.tsv output_ratio/flipped_intron.xlsx -draw output_ratio/flipped_intron -ca WT -cb KO -na Sense -nb BranchD -stat ratio
python bootstrap.py input/mmej.tsv output_ratio/mmej.xlsx -draw output_ratio/mmej -ca WT -cb KO -na Sense -nb BranchD -stat ratio
python bootstrap.py input/NHEJ.tsv output_ratio/NHEJ.xlsx -draw output_ratio/NHEJ -ca WT -cb KO -na Sense -nb BranchD -stat ratio
python bootstrap.py input/R-TDR.tsv output_ratio/R-TDR.xlsx -draw output_ratio/R-TDR -ca WT -cb KO -na Sense -nb BranchD -stat ratio
python bootstrap.py input/OX_WT.tsv output_ratio/OX_WT.xlsx -draw output_ratio/OX_WT -ca WT -cb OX -na Sense -nb BranchD -stat ratio
python bootstrap.py input/yeast.tsv output_ratio/yeast.xlsx -draw output_ratio/yeast -ca SPT3 -cb RHSPT3 -na Sense -nb BranchD -stat ratio