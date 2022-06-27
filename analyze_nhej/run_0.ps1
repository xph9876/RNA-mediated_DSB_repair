# Code for remove the blank lines from SAM files
# and rename files slightly
$DirFrom = "libraries_0"
$DirTo = "libraries_1"
$Libraries = Get-ChildItem $DirFrom
if (-not(Test-Path $DirTo)) {
  New-Item $DirTo -ItemType Directory
}
foreach ($File in $Libraries) {
  $Content = Get-Content $File.FullName
  $Content = $Content | Where-Object {($_ -ne "")}
  $NewName = $File.Name -replace "_100lines",""
  $Utf8NoBomEncoding = New-Object System.Text.UTF8Encoding $False
  Write-Output ($File.Name + " " + $NewName)
  [System.IO.File]::WriteAllLines(($DirTo + "\\" + $NewName), $Content, $Utf8NoBomEncoding)
}
