#!/bin/bash
#$ -cwd
# Plot all final simulation states for all temperatures using Mathematica by running the script plot_lattice.txt.
# need to set the simulation parameters below.

# Check input
if [ $# -ne 1 ]
then
echo Usage $0 [L]
exit
fi

L=$1
name="temp_Fe8_test8"
arrayMech=( "true" )
arrayH=( 0.0 )

for H in "${arrayH[@]}"
do

for mech in "${arrayMech[@]}"
do
# get all T's
IFS=$'\n,' read -d '' -ra temp_schedule <../"$SYS_NAME"/temperature_schedules/temp_schedule_"$L"_"$L"_"$name"_"$H"_"$mech".txt
for T in "${temp_schedule[@]}"; do
# make a new dir in /figures with the project name and a sub-directory with the current T
mkdir -p ../"$SYS_NAME"/figures/lattice_plot_"$L"_"$L"_"$name"_"$H"_"$mech"/T="$T"
for file in ../"$SYS_NAME"/data/lattice_output/"$name"/table_"$L"_"$L"_"$H"_"$T"_"$mech"_*.txt
do
seed=$(echo ${file##*_} | cut -f 1 -d '.')
# translate the spin index from the lattice output file to real space
# coordinates that can be used for plotting the spin in place
awk -f plot_lattice.awk -v L=$L $file > temp_lattice.txt
# run Mathematica with the script plot_lattice.txt
nohup nice /gpfs0/system/Mathematica/12.0.0/Executables/math -noinit -batchinput -batchoutput < plot_lattice.txt > /dev/null 2>&1 &
sleep 10
# move the resulting figure to its place in the figures dir
mv lattice_output.png ../"$SYS_NAME"/figures/lattice_plot_"$L"_"$L"_"$name"_"$H"_"$mech"/T="$T"/"$seed".png
mv temp_lattice.txt ../"$SYS_NAME"/figures/lattice_plot_"$L"_"$L"_"$name"_"$H"_"$mech"/T="$T"/"$seed".txt

done
done

done
done

rm temp_lattice.txt

