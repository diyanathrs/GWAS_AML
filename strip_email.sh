#!/bin/bash

TARGET="Wei-Yu.Lin@newcastle.ac.uk"


find . -type f -name "*.sh" | while read -r file; do
    echo "Processing: $file"
    sed -i "/$TARGET/d" "$file"
done
