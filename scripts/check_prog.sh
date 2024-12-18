#!/bin/bash
# Check the progress of a currently running job.
# Reads the last index from the job's results file and compares to the total
# number of sweeps it was given as a parameter.
# Run from the root project directory as bash ./scripts/check_prog.sh

# Check input
if [ $# -ne 1 ]
then
echo "Usage $0 [jobid | job name]"
exit
fi

# check is a job id was given or a job name, and the the other
if [[ $1 =~ ^[[:digit:]] ]]; then
    jobid=$1
    jobname=$(qstat -j $1 | grep job_name | awk '{print $2}')
else
    jobid=$(qstat -u tomerdol | tail -n +3 | awk '$3==name {print $1; exit}' name="$1")
    jobname=$1
fi

# get the args the job was run with
args=$(qstat -j ${jobid} | grep job_args | awk '{print $2}')
IFS=',' read -r -a args_array <<< "$args"

# find the job's results file
# array: $0=Lx $1=Lz $2=max_sweeps $3=extBx $4=mech $5=name $6=seed
files="./${SYS_NAME}/data/results/${args_array[5]}/table_${args_array[0]}_${args_array[1]}_${args_array[3]}_*_${args_array[4]}_${args_array[6]}.txt"

for f in ${files}; do

    ## Check if the glob gets expanded to existing files.
    ## If not, f here will be exactly the pattern above
    ## and the exists test will evaluate to false.
    if [ -e "$f" ]; then
    :
    else
    # if no files exist, it's the same as them being empty
    echo "0/${args_array[2]}"
    exit 0
    fi

    ## This is all we needed to know, so we can break after the first iteration
    break
done

temperature_array=()
last_step_array=()
# read last line from each temperature
while IFS= read -r line; do
    temperature_array+=( $( echo "$line" | cut -d: -f1 ) )
    # check that step is a number
    if [[ "$( echo "$line" | cut -d: -f2 )" =~ ^[0-9]+$ ]]; then
    last_step_array+=( $( echo "$line" | cut -d: -f2 ) )
    else
    last_step_array+=( 0 )
    fi
done < <( tail -n 1 ${files} | awk '$1=="==>" {n=split(substr($0, 5, length-8),a,"_"); next} $1!="" {print a[n-2]":"$1}' )

# check if all temperature are at the same last step
if [ "${#last_step_array[@]}" -gt 0 ] && [ $(printf "%s\000" "${last_step_array[@]}" |
       LC_ALL=C sort -z -u |
       grep -z -c .) -eq 1 ] ; then
  echo "${last_step_array[0]}/${args_array[2]}"
else
    echo -e "Temperature\tstep"
    count=0
    for t in "${temperature_array[@]}"; do
        echo -e "${t}\t${last_step_array[$count]}/${args_array[2]}"
        ((count++))
    done
fi

