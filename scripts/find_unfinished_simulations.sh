#!/bin/bash

name=$1
expected_last_line=$2

# Check input: number of arguments is 2 and 2nd argument is numeric
if [[ $# -ne 2 || ! $expected_last_line =~ ^[0-9]+$ ]]
then
  echo Usage $0 [directory] [expected last line]
exit
fi
# Check if the directory exists
if [ ! -d ./data/results/$name ]
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

res=$(find ./data/results/$name/table_* -maxdepth 0 | wc -l)
echo "found $res results"
i=1

all_files_ok=true
seeds_unfinished=()

for file in ./data/results/$name/table_*
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
if [ "$all_files_ok" = true ]; then
  echo "All files are ok!"
else
  echo "Some unfinished simulations were found:"
  for seed in "${seeds_unfinished[@]}"
  do
    # check whether the seed belongs to an already running job
    qstat -u tomerdol | awk 'NR>2 {print $1}' | xargs -n1 qstat -j | grep -q "$seed" &> /dev/null
    if [ $? == 0 ]; then
      jobid=$(qstat -u tomerdol | awk 'NR>2 {print $1}' | xargs -n1 qstat -j | grep -B26 ${seed} | awk 'NR==1 {print $2}')
      echo "* Seed ${seed} already running with job id ${jobid}."
    else
      echo "* Seed ${seed} not currently running."
    fi
  done
fi