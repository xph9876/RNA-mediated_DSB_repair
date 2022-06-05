# Code for remove the blank lines and headers of the SAM files
# and rename files slightly
$DirFrom = "libraries_0"
$DirTo = "libraries_1"
$Libraries = Get-ChildItem $DirFrom
if (-not(Test-Path $DirTo)) {
  New-Item $DirTo -ItemType Directory
}
foreach ($File in $Libraries) {
  $Content = Get-Content $File.FullName
  $Content = $Content | Where-Object {($_ -ne "") -and ($_ -notmatch "^@")}
  $NewName = $File.Name -replace "2DSBs","2DSB"
  $NewName = $NewName -replace "_100lines",""
  echo ($File.Name + " " + $NewName)
  $Content | Out-File -FilePath ($DirTo + "\\" + $NewName) -Encoding utf8
}

# Code for running the NHEJ scripts

