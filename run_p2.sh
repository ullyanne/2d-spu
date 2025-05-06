#!/bin/bash

# cd src/
# make clean
# make
# cd ..

dir="/home/ullyanne/Documents/2spp/instances/"
find "$dir" -type f -iname "*.txt" | while read file; do
    echo "Processing: $file"

    count=0
    for i in {1..20}; do
        # Print a seed to debug
        echo "Starting thread with seed $i for file $file" &

        # Executa o comando com o seed
        # src/samplecode -f "$file" -k 1 -p 279 -e 0.6760 -m 0.0039 -o 0.8688 -g 254 -a 14 -x 0.1260 -y 0.2503 -z 0.6237 -i 0.3508 -l 0.1204 -v 0.3320 -r 74 -s "$i" &

        src/samplecode -f "$file" -k 1 -p 392 -e 0.2567 -m 0.0202 -o 0.6001 -g 364 -a 20 -x 0.3186 -y 0.1852 -i 0.2365 -l 0.1399 -v 0.4838 -z 0.4962 -s "$i" &
        #src/samplecode -f "$file" -k 1 -p 279 -e 0.6760 -m 0.0039 -o 0.8688 -g 254 -a 14 -x 0.1260 -y 0.2503 -z 0.6237 -i 0.3508 -l 0.1204 -v 0.3320 -r 20 -s "$i" &


        # src/samplecode -f "$file" -k 1 -p 271 -e 0.4562 -m 0.0223 -o 0.7359 -g 400 -a 15 -x 0.3404 -y 0.13 -i 0.3338 -l 0.2665 -v 0.3653 -z 0.5296 -s 19 &

        # src/samplecode -f "$file" -k 1 -p 300 -e 0.5592 -m 0.0103 -o 0.6842 -g 327 -a 8 -x 0.2027 -y 0.2589 -i 0.3308 -l 0.2776 -v 0.4578 -s "$i" -z 0.5384 &

        # src/samplecode -f "$file" -k 1 -p 338 -e 0.611 -m 0.0135 -o 0.7978 -g 387 -a 14 -x 0.2204 -y 0.2544 -i 0.2809 -l 0.1597 -v 0.4525 -s "$i" -z 0.5252 &

        ((count++))
        if (( count % 5 == 0 )); then
            wait  # Espera as 3 execuções paralelas terminarem
        fi
    done

    wait  # Espera o restante se sobrar menos de 3 no final
done
