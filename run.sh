#!/bin/bash

dir="/home/ullyanne/Documents/2spp/instances/bke/T20/N1Burke.txt"

find "$dir" -type f -iname "*.txt" | while read file; do
    src/samplecode -f "$file" -k 3 -m 0.7 -e 0.3 -o 0.4 -p 2500 -t 14 -b 100 -g 1000
done
