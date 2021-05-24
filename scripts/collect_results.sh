#!/bin/bash

# remove previous "collected_results.txt" file
if [ -f collected_results.txt ]
then
rm collected_results.txt
fi

# loop over all files and collect results to "collected_results.txt"
for file in /gpfs0/smoshe/users/matityas/test.sh.o*
do
res=$(awk '{if(NR>9 && NF>0 && $1~/[0-9]+/) {print $1}}' ORS=',' "$file")
echo ${res::-1} >> collected_results.txt
done