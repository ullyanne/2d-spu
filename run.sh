#!/bin/bash

cd src/
make clean
make
cd ..

dir="/home/ullyanne/Documents/2spp/instances/bke/T20/"

find "$dir" -type f -iname "*.txt" | while read file; do
    src/samplecode -f "$file" -k 4 -m 0.7 -e 0.3 -o 0.6 -p 50 -t 1 -b 600 -g 100 -d
done