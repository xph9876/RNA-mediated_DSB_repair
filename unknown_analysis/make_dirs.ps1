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