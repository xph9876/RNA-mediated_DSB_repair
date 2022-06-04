# Code for remove the blank lines and headers of the SAM files
# and rename files slightly
$Libraries = Get-ChildItem libraries
if (-not(Test-Path "libraries2")) {
  New-Item "libraries2" -ItemType Directory
}
foreach ($File in $Libraries) {
  echo $File
  $Content = Get-Content $File.FullName
  $Content = $Content | Where-Object {($_ -ne "") -and ($_ -notmatch "^@")}
  $NewName = $File.Name -replace "2DSBs","2DSB"
  $Content | Out-File -FilePath ("libraries2\\" + $NewName) -Encoding utf8
}

# Code for running the NHEJ scripts

