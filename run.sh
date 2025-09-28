#!/bin/bash

cd src/
make clean
make
cd ..

dir="./instances/"
find "$dir" -type f -iname "*.txt" | while read file; do
    echo "$file"

    count=0
    for i in {1..20}; do

        echo "Starting thread with seed $i for file $file" &

        src/brkgarls -f "$file" -k 3 -g 241 -n 4 -p 475 -e 0.2121 -m 0.0442 -o 0.6514 -a 19 -x 0.2436 -y 0.0216 -i 0.4687 -l 0.1032 -v 0.4462 -z 0.7348 -s "$i" -t 60 &

        ((count++))
        if (( count % 3 == 0 )); then
            wait 
        fi
    done

    wait 
done
