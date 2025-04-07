#!/bin/bash

cd src/
make clean
make
cd ..

dir="/home/ullyanne/Documents/2spp/instances/bke/T20/"

find "$dir" -type f -iname "*.txt" | while read file; do
    for i in {1..1}; do
        src/samplecode -f "$file" -k 4 -m 0.3 -e 0.15 -o 0.65 -p 250 -t 14 -b 200 -g 200 -d
    done
done