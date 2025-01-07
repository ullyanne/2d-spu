#!/bin/bash

dir="/home/ullyanne/Documents/2spp/instances/bke/T40/N10Burke.txt"

find "$dir" -type f -iname "*.txt" | while read file; do
    src/samplecode -f "$file" -k 8 -m 0.3 -e 0.4 -o 0.7 -p 2500 -t 14 -b 1000 -g 1200
done
