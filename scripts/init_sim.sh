#!/bin/bash

source ./scripts/sub.sh

# this should override the definitions from sub.sh
name="res_test_const"
arrayMech=( "false" )
arrayL=( 6 5 4 )
arrayH=( 0.0 0.3 )
minT_false=1.4
maxT_false=1.85
minT_true=1.55
maxT_true=2.0
delta=0.05

gen_temp_schedules() {
# generate temporary temperature schedules
for L in "${arrayL[@]}"
do
for H in "${arrayH[@]}"
do

yes n | /gpfs0/smoshe/projects/Python-3.8_old/bin/python3 ./scripts/gen_temp_schedule.py "$1" "$2" 24 | head -n1 > ./temperature_schedules/temp_schedule_"$L"_"$L"_"$4""$name"_"$H"_"$3".txt
echo "generated ./temperature_schedules/temp_schedule_${L}_${L}_$4${name}_${H}_$3.txt"

done
done
}
if true; then
# generate temporary temperature schedules for "true"
gen_temp_schedules $minT_true $maxT_true "true" "temp_"
# generate temporary temperature schedules for "false"
gen_temp_schedules $minT_false $maxT_false "false" "temp_"

# run initial temporary run
sub "temp_" 513 20
fi
# wait until temporary runs finish
running_jobs=$(qstat -u tomerdol | awk '$3 ~ /^tr/' | awk 'END {print NR}')
while [ $running_jobs -gt 0 ]
do
echo "there are still $running_jobs jobs running. waiting 1 hour."
sleep 1h
running_jobs=$(qstat -u tomerdol | awk '$3 ~ /^tr/' | awk 'END {print NR}')
done

#*************************************************************************
# check that there are results
#*************************************************************************
echo "checking that there are results from temporary runs"

for mech in "${arrayMech[@]}"
do
for L in "${arrayL[@]}"
do
for H in "${arrayH[@]}"
do

IFS=$'\n,' read -d '' -ra temp_schedule < ./temperature_schedules/temp_schedule_"$L"_"$L"_temp_"$name"_"$H"_"$mech".txt

for T in "${temp_schedule[@]}"
do


for file in ./data/results/temp_"$name"/table_"$L"_"$L"_"$H"_"$T"_"$mech"_*.txt
do
# find if file exists and is larger than 5000 bytes = 5 KB
if ! [[ $(find "$file" -type f -size +5k 2>/dev/null) ]]; then
    echo "File ${file} does not exist or is smaller than 5 KB. Exiting."
    exit 1
else
    # counts how many lines have a number in the first column
    num_of_result_lines=$(awk '{if ($1 ~ /^[0-9]*$/) {count++}} END {print count}' ${file})
    if [ $num_of_result_lines -eq 0 ]; then
        echo "File ${file} exists but does not seem to contain results. Exiting."
        exit 1
    fi
fi
done

done
done
done
done
#*************************************************************************


# run extract_data_EO.sh for each simulation
#*************************************************************************
echo "extracting data from temporary runs"
for mech in "${arrayMech[@]}"
do
for L in "${arrayL[@]}"
do
for H in "${arrayH[@]}"
do

IFS=$'\n,' read -d '' -ra temp_schedule < ./temperature_schedules/temp_schedule_"$L"_"$L"_temp_"$name"_"$H"_"$mech".txt

echo "reading temperatures from ./temperature_schedules/temp_schedule_${L}_${L}_temp_${name}_${H}_${mech}.txt"
echo "T e corr_length"

# Do the analysis
tmp_file="sample_energy.txt"
if [ -f $tmp_file ]
then
rm $tmp_file
fi
echo T e corr_length >> $tmp_file

for T in "${temp_schedule[@]}"
do

file="./data/results/temp_""$name""/table_""$L""_""$L""_""$H""_""$T""_""$mech""_*.txt"

num_of_samples=$(tail -q -n 1 ${file} | awk 'END{print $1}')
start_point=$((num_of_samples / 4))
energy=$(tail -q -n ${start_point} ${file} | awk '{ total += $3 } END { print total/NR }')
m2=$(tail -q -n ${start_point} ${file} | awk '{ total += ($2)^2 } END { print total/NR }')
mk2=$(tail -q -n ${start_point} ${file} | awk '{ total += $15 } END { print total/NR }')

corr_length=$(bc -l <<< "scale=5; sqrt(($m2/$mk2)-1)/(2*$L*s(4*a(1)/$L))")

echo $T $energy $corr_length >> $tmp_file
echo $T $energy $corr_length
done

sort -t" " -nk1 $tmp_file > ./data/analysis/sample_energy_"$L"_"$H"_temp_"$name"_"$mech".txt
rm $tmp_file

done
done
done

#*************************************************************************


# find initial Tc and create new temperature schedules
#*************************************************************************
echo "looking for initial Tc and creating new temperature schedules"
for mech in "${arrayMech[@]}"
do
for H in "${arrayH[@]}"
do

all_L="${arrayL[@]}"
initial_tc=$(/gpfs0/smoshe/projects/Python-3.8_old/bin/python3 ./scripts/find_initial_tc.py --h_ex_list ${H} --mech ${mech} --folder_list temp_${name} -L ${all_L})
echo "Hex=$H mech=$mech : initial Tc = $initial_tc"

for L in "${arrayL[@]}"
do
bash ./scripts/temp_set_script.sh ./data/analysis/sample_energy_"$L"_"$H"_temp_"$name"_"$mech".txt 24 $(bc <<< "$initial_tc-$delta") $(bc <<< "$initial_tc+$delta") | awk -vORS=, '{ print $1 }' | sed 's/,$/\n/' > ./temperature_schedules/temp_schedule_"$L"_"$L"_"$name"_"$H"_"$mech".txt
echo "created ./temperature_schedules/temp_schedule_${L}_${L}_${name}_${H}_${mech}.txt"
done

done
done
#*************************************************************************


# run simulations
echo "running simulations"
sub "" 4096 50
echo "done"
exit 0
