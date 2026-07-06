#!/bin/bash

TARGET="Wei-Yu.Lin@newcastle.ac.uk"

find . -type f -name "*.slurm" | while read -r file; do
    echo "Processing: $file"
    sed -i "/$TARGET/d" "$file"
done
