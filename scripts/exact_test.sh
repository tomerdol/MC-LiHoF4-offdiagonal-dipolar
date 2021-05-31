#!/bin/bash

arrayH=( 0.0 0.04081633 0.08163265 0.12244898 0.16326531 0.20408163 0.24489796 0.28571429 0.32653061 0.36734694 0.40816327 0.44897959 0.48979592 0.53061224 0.57142857 0.6122449 0.65306122 0.69387755 0.73469388 0.7755102 0.81632653 0.85714286 0.89795918 0.93877551 0.97959184 1.02040816 1.06122449 1.10204082 1.14285714 1.18367347 1.2244898 1.26530612 1.30612245 1.34693878 1.3877551 1.42857143 1.46938776 1.51020408 1.55102041 1.59183673 1.63265306 1.67346939 1.71428571 1.75510204 1.79591837 1.83673469 1.87755102 1.91836735 1.95918367 2.0 )

name="exact_test"
L=1
T=1.53
mech="false"

# Do the analysis
tmp_file="sample_energy.txt"
if [ -f $tmp_file ]
then
rm $tmp_file
fi
echo H e >> $tmp_file

for H in "${arrayH[@]}"
do

file="./"$SYS_NAME"/data/results/""$name""/table_""$L""_""$L""_""$H""_""$T""_""$mech""_*.txt"

num_of_samples=$(tail -q -n 1 ${file} | awk 'END{print $1}')
if [ $num_of_samples != "index" ]
then
start_point=$((num_of_samples / 2))
energy=$(tail -q -n ${start_point} ${file} | awk '!/^ *#/{ total += $3/4; count++ } END { print total/count }')
magnetization=$(tail -q -n ${start_point} ${file} | awk '{ total += $2 } END { print total/NR }')
else
energy=""
magnetization=""
fi

echo $H $energy $magnetization >> $tmp_file
echo $H $energy $magnetization
done
