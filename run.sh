#!/bin/bash

dir="/home/ullyanne/Documents/2spp/instances/bke/T20/N1Burke.txt"

find "$dir" -type f -iname "*.txt" | while read file; do
    src/samplecode -f "$file" -k 5 -m 0.3 -e 0.7 -o 0.8 -p 2500 -t 14 -b 1000 -g 1000
done
