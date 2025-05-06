#!/bin/bash

# cd src/
# make clean
# make
# cd ..

dir="/home/ullyanne/Documents/2spp/instances/bke/T20/N12Burke.txt"
find "$dir" -type f -iname "*.txt" | while read file; do
    echo "$file"
    for i in {1..6}; do
        src/samplecode -f "$file" -k 1 -p 332 -e 0.5256 -m 0.0104 -o 0.7411 -g 300 -a 14 -x 0.4440 -y 0.4856 -z 0.0704 -i 0.1 -l 0.2003 -v 0.2566 -s "$i" -d
    done
done