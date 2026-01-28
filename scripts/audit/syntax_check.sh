#!/bin/bash
set -e

# Find all .R files in scripts/ (recursive)
echo "Scanning scripts/ for .R files..."
files=$(find scripts -name "*.R")

count=0
errors=0
failed_files=()

for file in $files; do
    count=$((count + 1))
    echo -n "Checking $file ... "
    
    if Rscript -e "tryCatch(parse(file='$file'), error=function(e) quit(status=1))" > /dev/null 2>&1; then
        echo "OK"
    else
        echo "FAIL"
        errors=$((errors + 1))
        failed_files+=("$file")
    fi
done

echo "----------------------------------------"
echo "Checked $count files."
echo "Found $errors syntax errors."

if [ $errors -gt 0 ]; then
    echo "Failed files:"
    for f in "${failed_files[@]}"; do
        echo "- $f"
    done
    exit 1
else
    echo "All files passed syntax check."
    exit 0
fi
