if (!(Test-Path output)) {
    New-Item -ItemType Directory output
}
if (!(Test-Path output/alignment)) {
  New-Item -ItemType Directory output/alignment
}
if (!(Test-Path output/nhej_mmej)) {
  New-Item -ItemType Directory output/nhej_mmej
}
if (!(Test-Path output/alignment_full)) {
  New-Item -ItemType Directory output/alignment_full
}
if (!(Test-Path output/alignment_full/yjl271_WT_sgB_R2_branch)) {
  New-Item -ItemType Directory output/alignment_full/yjl271_WT_sgB_R2_branch
}

python .\make_full_output.py -i .\output\alignment\yjl271_WT_sgB_R2_branch.csv -o .\output\alignment_full\yjl271_WT_sgB_R2_branch

python .\make_summary_output.py -i .\output\alignment_full\yjl271_WT_sgB_R2_branch\unknown.csv -o .\output\summary\yjl271_WT_sgB_R2_branch -m unknown
python .\make_summary_output.py -i .\output\alignment_full\yjl271_WT_sgB_R2_branch\mmej.csv -o .\output\summary\yjl271_WT_sgB_R2_branch -m mmej
python .\make_summary_output.py -i .\output\nhej_mmej\yjl271_WT_sgB_R2_branch.csv -o .\output\summary\yjl271_WT_sgB_R2_branch -m nhej_mmej