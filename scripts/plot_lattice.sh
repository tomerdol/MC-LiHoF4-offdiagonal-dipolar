#!/bin/bash
#$ -cwd

L=$1
name="temp_Fe8_test8"
arrayMech=( "true" )
arrayH=( 0.0 )

for H in "${arrayH[@]}"
do

for mech in "${arrayMech[@]}"
do

IFS=$'\n,' read -d '' -ra temp_schedule <../"$SYS_NAME"/temperature_schedules/temp_schedule_"$L"_"$L"_"$name"_"$H"_"$mech".txt
for T in "${temp_schedule[@]}"; do
mkdir -p ../"$SYS_NAME"/figures/lattice_plot_"$L"_"$L"_"$name"_"$H"_"$mech"/T="$T"
for file in ../"$SYS_NAME"/data/lattice_output/"$name"/table_"$L"_"$L"_"$H"_"$T"_"$mech"_*.txt
do
seed=$(echo ${file##*_} | cut -f 1 -d '.')
awk -f plot_lattice.awk -v L=$L $file > temp_lattice.txt
nohup nice /gpfs0/system/Mathematica/12.0.0/Executables/math -noinit -batchinput -batchoutput < plot_lattice.txt > /dev/null 2>&1 &
sleep 10
mv lattice_output.png ../"$SYS_NAME"/figures/lattice_plot_"$L"_"$L"_"$name"_"$H"_"$mech"/T="$T"/"$seed".png
mv temp_lattice.txt ../"$SYS_NAME"/figures/lattice_plot_"$L"_"$L"_"$name"_"$H"_"$mech"/T="$T"/"$seed".txt

done
done

done
done

rm temp_lattice.txt

