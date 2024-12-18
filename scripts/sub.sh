#!/bin/bash
# script for submitting new simulations according to the defined parameters.
# remember to set x and extra_par in the function itself.
# Usage: source ./scripts/sub.sh; sub [job name initial: "" | temp_] [# of sweeps] [# of independent runs] [# of temperatures]
# E.g., source ./scripts/sub.sh; sub "" 8200 1000 20

name="low_jex_dilution_1.0_2"
arrayMech=( "false" "true" )
arrayL=( 6 7 8 )
arrayH=( 0.0 )

#name="exact_test_const"
#arrayH=( 0.0 0.04081633 0.08163265 0.12244898 0.16326531 0.20408163 0.24489796 0.28571429 0.32653061 0.36734694 0.40816327 0.44897959 0.48979592 0.53061224 0.57142857 0.6122449 0.65306122 0.69387755 0.73469388 0.7755102 0.81632653 0.85714286 0.89795918 0.93877551 0.97959184 1.02040816 1.06122449 1.10204082 1.14285714 1.18367347 1.2244898 1.26530612 1.30612245 1.34693878 1.3877551 1.42857143 1.46938776 1.51020408 1.55102041 1.59183673 1.63265306 1.67346939 1.71428571 1.75510204 1.79591837 1.83673469 1.87755102 1.91836735 1.95918367 2.0 )

#name="Fe8_test8"
#arrayMech=( "false" "true" )
#arrayL=( 5 6 7 8 9 )
#arrayLexclude=()
#arrayH=( 0.0 )
#minT_false=0.3
#maxT_false=0.9
#minT_true=0.3
#maxT_true=0.9
#delta=0.08

sub() {
temp="$1"
max_sweeps="$2"
runs="$3"
nT="$4"
extra_par="_const"
x="1.0"
jex="1.16e-3"
#jex="3.91e-3"

temp_initial="$(echo $temp | head -c 1)"

local_name="${temp}${name}"

num_of_Ls=$((${#arrayL[@]} - ${#arrayLexclude[@]}))
total=$((${#arrayMech[@]} * ${num_of_Ls} * ${#arrayH[@]} * $runs))
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
if [[ ! " ${arrayLexclude[@]} " =~ " ${L} " ]]; then
COUNT=0
while [ $COUNT -lt $runs ]; do
    #used_slots=`free_slot smoshe.q | grep sge1081 | cut -d' ' -f 2 | cut -d'/' -f 1`
    #exclude sge1081 to leave at least 30 cores available for other users
    #queues="lublin.q,smoshe.q@sge1081,smoshe.q@sge1082,smoshe.q@sge190,smoshe.q@sge247,smoshe.q@sge249,fairshare.q"
    
    #if [ "$used_slots" != "" ]; then
    #if [ $used_slots -le $nT ]; then
    queues="lublin.q,smoshe.q@sge1081,smoshe.q@sge1082,smoshe.q@sge190,smoshe.q@sge247"
    #fi
    #fi

    # parallel submission
    qsub -pe shared $nT -l mem_free=40G -V -S /bin/bash -cwd -N "$temp_initial"r"$L"_"$H"_"$COUNT"_"$mech_initial" -o ./"$SYS_NAME"/output/ -e ./"$SYS_NAME"/output/ -q "$queues" scripts/met_with_t.sh "$L" "$L" "$max_sweeps" "$H" "$mech" "$local_name" "${seeds[$i]}" "$extra_par" "$x" "$jex"
    # sequential submission
    #qsub -l mem_free=2G -V -S /bin/bash -cwd -N "$temp_initial"r"$L"_"$H"_"$COUNT"_"$mech_initial" -o ./"$SYS_NAME"/output/ -e ./"$SYS_NAME"/output/ -q smoshe.q@sge247,smoshe.q@sge1081,smoshe.q@sge1082,lublin.q,intel_all.q scripts/met_with_t_single.sh "$L" "$L" "$max_sweeps" "$H" "$mech" "$local_name" "${seeds[$i]}" "$extra_par" "$x" "$jex"
    #echo "${seeds[$i]}"
    ((i++))
    ((COUNT++))
    sleep 30
done
fi
done
done
done
}

