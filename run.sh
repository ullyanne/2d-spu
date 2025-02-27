#!/bin/bash

cd src/
make clean
make
cd ..

dir="/home/ullyanne/Documents/2spp/instances/bke/T20/N3Burke.txt"

find "$dir" -type f -iname "*.txt" | while read file; do
    src/samplecode -f "$file" -k 8 -m 0.7 -e 0.3 -o 0.8 -p 500 -t 1 -b 600 -g 400 -d
done