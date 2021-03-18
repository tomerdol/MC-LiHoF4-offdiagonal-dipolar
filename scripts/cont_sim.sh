#!/bin/bash
unset module
name="temp_Fe8_test"
arrayMech=( "true" )
arrayL=( 3 4 5 )
arrayH=( 0.0 )

# Check input
if [[ ! $# > 0 || ! $1 =~ ^[0-9]+$ ]]
then
echo "Usage: $0 [# of sweeps] [extra parameter]"
exit
fi

max_sweeps=$1
extra=$2
i=0

for H in "${arrayH[@]}"
do

for mech in "${arrayMech[@]}"
do
mech_initial="$(echo $mech | head -c 1)"

for L in "${arrayL[@]}"
do

COUNT=0
for file in ./"$SYS_NAME"/checkpoints/"$name"/save_state_"$L"_"$L"_"$H"_"$mech"_*.txt
do
    seed=$(echo ${file##*_} | cut -f 1 -d '.')
    # check that the seed does not belong to an already running job
    qstat -u tomerdol | awk 'NR>2 {print $1}' | xargs -n1 qstat -j | grep -q "$seed" &> /dev/null
    if [ $? == 0 ]; then
      jobid=$(qstat -u tomerdol | awk 'NR>2 {print $1}' | xargs -n1 qstat -j | grep -B26 ${seed} | awk 'NR==1 {print $2}')
      echo "Seed ${seed} already running with job id ${jobid}."
    else
      used_slots=`free_slot smoshe.q | grep sge1081 | cut -d' ' -f 2 | cut -d'/' -f 1`
      #exclude sge1081 to leave at least 30 cores available for other users
      queues="lublin.q,smoshe.q@sge1082,smoshe.q@sge190,fairshare.q,smoshe.q@sge247,smoshe.q@sge249"

      if [ "$used_slots" != "" ]; then
      if [ $used_slots -le 24 ]; then
      queues="lublin.q,smoshe.q,fairshare.q"
      fi
      fi

      qsub -pe shared 24 -l mem_free=40G -V -S /bin/bash -cwd -N tr"$L"_"$H"_"$COUNT"_"$mech_initial" -o ./"$SYS_NAME"/output/ -e ./"$SYS_NAME"/output/ -q "$queues" ./scripts/met_with_t.sh "$L" "$L" "$max_sweeps" "$H" "$mech" "$name" "$seed" "$extra"
      #echo "${seeds[$i]}"

      ((COUNT++))
      sleep 120
    fi
done
done
done
done

