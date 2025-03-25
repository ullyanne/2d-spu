#!/bin/bash

cd src/
make clean
make
cd ..

dir="/home/ullyanne/Documents/2spp/instances/bke/T20/N11Burke.txt"

find "$dir" -type f -iname "*.txt" | while read file; do
    for i in {1..1}; do
        src/samplecode -f "$file" -k 2 -m 0.7 -e 0.3 -o 0.6 -p 60 -t 14 -b 100 -g 100 -d
    done
done