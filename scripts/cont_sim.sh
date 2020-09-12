#!/bin/bash
unset module
name="temp_parallel_test2"
arrayMech=( "false" "true" )
arrayL=( 6 )
arrayH=( 0.0 0.3 )


max_sweeps=$1

i=0

for H in "${arrayH[@]}"
do

for mech in "${arrayMech[@]}"
do
mech_initial="$(echo $mech | head -c 1)"

for L in "${arrayL[@]}"
do

COUNT=0
for file in ./checkpoints/"$name"/save_state_"$L"_"$L"_"$H"_"$mech"_*.txt
do
    seed=$(echo ${file##*_} | cut -f 1 -d '.')
    
    used_slots=`free_slot smoshe.q | grep sge1081 | cut -d' ' -f 2 | cut -d'/' -f 1`
    #exclude sge1081 to leave at least 30 cores available for other users
    queues="lublin.q,smoshe.q@sge1082,smoshe.q@sge190,fairshare.q"
    
    if [ "$used_slots" != "" ]; then
    if [ $used_slots -le 24 ]; then
    queues="lublin.q,smoshe.q,fairshare.q"
    fi
    fi

    qsub -pe shared 24 -l mem_free=40G -V -S /bin/bash -cwd -N tr"$L"_"$H"_"$COUNT"_"$mech_initial" -o ./temp_output/ -e ./temp_output/ -q "$queues" met_with_t.sh "$L" "$L" "$max_sweeps" "$H" "$mech" "$name" "$seed"
    #echo "${seeds[$i]}"
    
    ((COUNT++))
    sleep 120
done
done
done
done

