# flatten_pcs.ps1
$outputFile = "flattened_pcs_repo.txt"
$files = Get-ChildItem -Recurse -Include *.cpp, *.hpp, *.h, *.c

if (Test-Path $outputFile) { Remove-Item $outputFile }

foreach ($file in $files) {
    "--------------------------------------------------" | Out-File -Append -FilePath $outputFile
    "FILE: $($file.FullName)" | Out-File -Append -FilePath $outputFile
    "--------------------------------------------------" | Out-File -Append -FilePath $outputFile
    Get-Content $file.FullName | Out-File -Append -FilePath $outputFile
    "" | Out-File -Append -FilePath $outputFile
}

Write-Host "Done! Flattened into $outputFile"