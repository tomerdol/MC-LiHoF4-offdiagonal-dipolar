#!/bin/bash

while read -r seed; do
    echo "deleting $seed"
    rm -f ../LiHoF4/data/results/high_jex_dilution_0.46/*${seed}.txt
    rm -f ../LiHoF4/data/lattice_output/high_jex_dilution_0.46/*${seed}.txt
    rm -f ../LiHoF4/checkpoints/high_jex_dilution_0.46/*${seed}.txt
    # rm "$path"
done < Seedlist
