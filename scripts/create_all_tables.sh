#!/bin/sh

arrayH=( 0.6 1.0 )

for H in "${arrayH[@]}"
do
qsub -V -S /bin/bash -cwd -N table"$H" -q smoshe.q scripts/create_one_table.sh "$H"
sleep 30
done
