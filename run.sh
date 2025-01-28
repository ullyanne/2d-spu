#!/bin/bash

dir="/home/ullyanne/Documents/2spp/instances/bke/T20/N3Burke.txt"

find "$dir" -type f -iname "*.txt" | while read file; do
    src/samplecode -f "$file" -k 4 -m 0.8 -e 0.2 -o 0.4 -p 300 -t 14 -b 300 -g 100
done