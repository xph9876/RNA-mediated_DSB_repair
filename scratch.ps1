param (
    [string]$server = "http://defaultserver"
)

if ($server -eq 10) {
  Write-Output "ten" hello
} else {
  Write-Output $server HELLO
}