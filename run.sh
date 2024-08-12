#!/bin/bash

dir="instances/bke"

find "$dir" -type f -iname "*.txt" | while read file; do
    src/samplecode -f "$file"
done
