#!/bin/bash
# Continue a simulation that was previously stopped (either finished its run or aborted)
# Uncomment either parallel submission of serial submission below

unset module

# simulation parameters
name="high_jex_dilution_0.67"
arrayMech=( "true" )
arrayL=( 8 )
arrayH=( 0.0 )
nT=24
x=0.67
#jex="1.16e-3"
jex="3.91e-3"

#arrayH=( 0.0 0.04081633 0.08163265 0.12244898 0.16326531 0.20408163 0.24489796 0.28571429 0.32653061 0.36734694 0.40816327 0.44897959 0.48979592 0.53061224 0.57142857 0.6122449 0.65306122 0.69387755 0.73469388 0.7755102 0.81632653 0.85714286 0.89795918 0.93877551 0.97959184 1.02040816 1.06122449 1.10204082 1.14285714 1.18367347 1.2244898 1.26530612 1.30612245 1.34693878 1.3877551 1.42857143 1.46938776 1.51020408 1.55102041 1.59183673 1.63265306 1.67346939 1.71428571 1.75510204 1.79591837 1.83673469 1.87755102 1.91836735 1.95918367 2.0 )

# Check input
if [[ ! $# > 0 || ! $1 =~ ^[0-9]+$ ]]
then
echo "Usage: $0 [# of sweeps] [extra parameter]"
exit
fi

max_sweeps=$1
extra=$2
i=0

# loop through the defined simulations
for H in "${arrayH[@]}"
do

for mech in "${arrayMech[@]}"
do
# get the first character of the mechanism definition
mech_initial="$(echo $mech | head -c 1)"

for L in "${arrayL[@]}"
do

COUNT=0

# loop through all checkpoint files that correspond to the current simulation params
for file in ./"$SYS_NAME"/checkpoints/"$name"/save_state_"$L"_"$L"_"$H"_"$mech"_*.txt
do
    seed=$(echo ${file##*_} | cut -f 1 -d '.')
    # check that the seed does not belong to an already running job
    qstat -u tomerdol | awk 'NR>2 {print $1}' | xargs -n1 qstat -j | grep -q "$seed" &> /dev/null
    if [ $? == 0 ]; then
      jobid=$(qstat -u tomerdol | awk 'NR>2 {print $1}' | xargs -n1 qstat -j | grep -B26 ${seed} | awk 'NR==1 {print $2}')
      echo "Seed ${seed} already running with job id ${jobid}."
    else
      queues="lublin.q,smoshe.q@sge1082,smoshe.q@sge190,smoshe.q@sge247,smoshe.q@sge1081"

      # parallel submission
      qsub -pe shared "$nT" -l mem_free=40G -V -S /bin/bash -cwd -N r"$L"_"$H"_"$COUNT"_"$mech_initial" -o ./"$SYS_NAME"/output/ -e ./"$SYS_NAME"/output/ -q "$queues" ./scripts/met_with_t.sh "$L" "$L" "$max_sweeps" "$H" "$mech" "$name" "$seed" "$extra" "$x" "$jex"

      # single submission
      #qsub -l mem_free=4G -V -S /bin/bash -cwd -N r"$L"_"$H"_"$COUNT"_"$mech_initial" -o ./"$SYS_NAME"/output/ -e ./"$SYS_NAME"/output/ -q smoshe.q@sge1082,smoshe.q@sge190,smoshe.q@sge247,smoshe.q@sge1081,lublin.q,intel_all.q ./scripts/met_with_t_single.sh "$L" "$L" "$max_sweeps" "$H" "$mech" "$name" "$seed" "$extra" "$x" "$jex"

      ((COUNT++))
      sleep 10
    fi
done
done
done
done

