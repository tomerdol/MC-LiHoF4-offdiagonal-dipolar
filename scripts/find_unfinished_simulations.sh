#!/bin/bash

# Check input: number of arguments is 1 or 2 and if 2 arguments are given then 2nd argument is numeric
if [[ $# -ne 1 && $# -ne 2 ]] || [[ $# -eq 2 && ! $2 =~ ^[0-9]+$ ]]
then
  echo Usage
  echo $0 [directory] [expected last line]
  echo Or
  echo echo $0 [directory]
exit
fi

name=$1
expected_last_line=$2

# Check if the directory exists
if [ ! -d ./$SYS_NAME/data/results/$name ]
then
echo Directory $name not found.
exit
fi

# Draw progress bar
function redraw_progress_bar { # int barsize, int base, int i, int top
    local barsize=$1
    local base=$2
    local current=$3
    local top=$4
    local j=0
    local progress=$(( ($barsize * ( $current - $base )) / ($top - $base ) ))
    echo -n "["
    for ((j=0; j < $progress; j++)) ; do echo -n '='; done
    echo -n '=>'
    for ((j=$progress; j < $barsize ; j++)) ; do echo -n ' '; done
    echo -n "] $(( $current )) / $top " $'\r'
}

function find_last_line {
  local seed=$1
  for file in ./"$SYS_NAME"/data/results/$name/table_*"$seed".txt
  do
  last_line=$(sed '/^#/d' "$file" | tail -n1 | awk '{print $1}')
  break
  done
  echo "$last_line"
}

res=$(find ./$SYS_NAME/data/results/$name/table_* -maxdepth 0 | wc -l)
echo "found $res results"

i=1
# find simulation with max written steps
if [[ $2 == "" ]]; then
    echo "Looking for simulation with maximum written steps..."
    expected_last_line=0
    for file in ./"$SYS_NAME"/data/results/$name/table_*
    do
      redraw_progress_bar 50 0 $i $res
      ((i++))
      last_line=$(sed '/^#/d' "$file" | tail -n1 | awk '{print $1}')
      if [ "$last_line" == "index" ]; then
        last_line=0
      fi
      if [ "$last_line" -gt "$expected_last_line" ]
      then
        expected_last_line=$last_line
      fi
    done
fi

i=1
all_files_ok=true
seeds_unfinished=()
printf "\n"
echo "Maximum # of steps: ${expected_last_line}"
printf "\n"
echo "Looking for unfinished simulations..."
for file in ./"$SYS_NAME"/data/results/$name/table_*
do
  redraw_progress_bar 50 0 $i $res
  ((i++))
  last_line=$(sed '/^#/d' "$file" | tail -n1 | awk '{print $1}')
  if [ "$last_line" == "index" ]; then
    last_line=0
  fi
  if [ "$last_line" -ne "$expected_last_line" ]
  then
    echo "File ${file} has ${last_line} lines!"
    all_files_ok=false
    seed=$(echo ${file##*_} | cut -f 1 -d '.')
    if [[ ! " ${seeds_unfinished[@]} " =~ " ${seed} " ]]; then
    seeds_unfinished+=("$seed")
    fi
  fi
done

printf "\n"

#last_bin=$(bc -l <<<"2*((2^$last_line)-1)")
#if [ "$all_files_ok" = true ]; then
#  echo "All files are ok!"
#  echo "All simulations have ${last_line} lines written to output (${last_bin} bins)."
#else
#  echo "Some unfinished simulations were found:"
#  for seed in "${seeds_unfinished[@]}"
#  do
#    last_line=$(find_last_line "$seed")
#    # check whether the seed belongs to an already running job
#    qstat -u tomerdol | awk 'NR>2 {print $1}' | xargs -n1 qstat -j | grep -q "$seed" &> /dev/null
#    if [ $? == 0 ]; then
#      jobid=$(qstat -u tomerdol | awk 'NR>2 {print $1}' | xargs -n1 qstat -j | grep -B26 ${seed} | awk 'NR==1 {print $2}')
#      echo "* Seed ${seed} already running with job id ${jobid}. (last line=${last_line})"
#    else
#      echo "* Seed ${seed} not currently running. (last line=${last_line})"
#    fi
#  done
#fi

declare -A seed_jobid
last_bin=$(bc -l <<<"2*((2^$last_line)-1)")
if [ "$all_files_ok" = true ]; then
  echo "All files are ok!"
  echo "All simulations have ${last_line} lines written to output (${last_bin} bins)."
else
  echo "Some unfinished simulations were found:"

  # create hash table of seed:jobid
  while read -r jobid ; do
    seed=$(qstat -j "$jobid" | grep job_args: | awk -F"," '{print $7}')
    if [[ "$seed" != "" ]]; then
      seed_jobid["$seed"]="$jobid"
    fi
  done < <(qstat -u tomerdol | awk 'NR>2 {print $1}')

#  for x in "${!seed_jobid[@]}"; do printf "[%s]=%s\n" "$x" "${seed_jobid[$x]}" ; done

  # go over unfinished seeds and check hash table
  for seed in "${seeds_unfinished[@]}"
  do
    last_line=$(find_last_line "$seed")
    # check whether the seed belongs to an already running job
    if [ ${seed_jobid[${seed}]+_} ]; then
      jobid=seed_jobid[${seed}]
      echo "* Seed ${seed} already running with job id ${jobid}. (last line=${last_line} / ${expected_last_line})"
    else
      echo "* Seed ${seed} not currently running. (last line=${last_line} / ${expected_last_line})"
    fi
  done
fi