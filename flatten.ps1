# flatten.ps1
$outputFile = "flattened_pcs_repo.txt"

# Define the file extensions we care about (Source, Headers, Makefile)
$extensions = @("*.cpp", "*.c", "*.hpp", "*.h", "Makefile")

# Remove the output file if it already exists so we don't append to old data
if (Test-Path $outputFile) {
    Remove-Item $outputFile
}

Write-Host "Flattening repository..."

# Find all matching files in the current directory
Get-ChildItem -Path . -Include $extensions -File | ForEach-Object {
    Write-Host "Adding $($_.Name)..."
    
    # Create the separator header
    $header = "`n" + ("=" * 80) + "`n"
    $header += "File: $($_.Name)`n"
    $header += ("=" * 80) + "`n"
    
    # Read the file content
    $content = Get-Content $_.FullName -Raw
    
    # Append to the output file
    Add-Content -Path $outputFile -Value ($header + $content)
}

Write-Host "Done! Please upload: $outputFile"