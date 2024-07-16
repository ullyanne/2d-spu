#!/bin/bash

dir="instances/bke/T20"

find "$dir" -type f -name "*.txt" | while read file; do
    src/samplecode -f "$file"
done
