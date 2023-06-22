if (!(Test-Path output)) {
    New-Item -ItemType Directory output
}
if (!(Test-Path output/collected)) {
  New-Item -ItemType Directory output/collected
}
