if (!(Test-Path output)) {
    New-Item -ItemType Directory output
}
if (!(Test-Path output/alignment)) {
  New-Item -ItemType Directory output/alignment
}
if (!(Test-Path output/nhej_mmej)) {
  New-Item -ItemType Directory output/nhej_mmej
}
