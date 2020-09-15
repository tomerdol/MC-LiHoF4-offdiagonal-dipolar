#!/bin/sh

arrayH=( 0.0 0.3 )

for H in "${arrayH[@]}"
do
qsub -V -S /bin/bash -cwd -N table"$H" -q smoshe.q scripts/create_one_table.sh "$H"
sleep 30
done