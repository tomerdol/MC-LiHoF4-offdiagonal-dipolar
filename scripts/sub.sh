#!/bin/bash

name="parallel_test2"
arrayMech=( "false" "true" )
arrayL=( 8 7 6 )
arrayH=( 0.0 0.3 )
minT_false=1.4
maxT_false=1.85
minT_true=1.55
maxT_true=2.0

sub() {
temp="$1"
max_sweeps="$2"
runs="$3"

temp_initial="$(echo $temp | head -c 1)"

local_name="${temp}${name}"


total=$((${#arrayMech[@]} * ${#arrayL[@]} * ${#arrayH[@]} * $runs))
echo "total num of simulations to run is $total"
# generate the required number of seeds (len(arrayL)*len(arrayH)) and store in seeds array. these are guaranteed to have no duplicates
seeds=()
echo "generating seeds"
while IFS= read -r line; do
	seeds+=( "$line" )
done < <( java -classpath "./bin/production/LiHoF4_transverse_field_MC/:./lib/*" utilities.GenerateSeeds "${JOB_ID}" "$total" )
i=0

for H in "${arrayH[@]}"
do

for mech in "${arrayMech[@]}"
do
mech_initial="$(echo $mech | head -c 1)"

for L in "${arrayL[@]}"
do

COUNT=0
while [ $COUNT -lt $runs ]; do
    used_slots=`free_slot smoshe.q | grep sge1081 | cut -d' ' -f 2 | cut -d'/' -f 1`
    #exclude sge1081 to leave at least 30 cores available for other users
    queues="lublin.q,smoshe.q@sge1082,smoshe.q@sge190,fairshare.q"
    
    if [ "$used_slots" != "" ]; then
    if [ $used_slots -le 24 ]; then
    queues="lublin.q,smoshe.q,fairshare.q"
    fi
    fi

    qsub -pe shared 24 -l mem_free=40G -V -S /bin/bash -cwd -N "$temp_initial"r"$L"_"$H"_"$COUNT"_"$mech_initial" -o ./output/ -e ./output/ -q "$queues" met_with_t.sh "$L" "$L" "$max_sweeps" "$H" "$mech" "$local_name" "${seeds[$i]}"
    #echo "${seeds[$i]}"
    ((i++))
    ((COUNT++))
    sleep 120
done
done
done
done
}

