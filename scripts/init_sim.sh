#!/bin/bash
# Initiates a simulation with the defined parmeters.
# Starts by generating a (geometric) temperature schedule over a large range,
# then runs a short simulation and analyzes the preliminary results to get a
# reasonable estimate of T_c. Uses that T_c to run a longer simulation with
# optimized temperature schedule close to that T_c.

source ./scripts/sub.sh

# number of temperatures to run (also affects the number of CPUs requested when running parallel jobs)
nT=24

# this should override the definitions from sub.sh
if [ $SYS_NAME == "LiHoF4" ]; then
name="low_jex_dilution_1.0"
arrayMech=( "true" "false" )
arrayL=( 6 7 8 )
# if adding a new L to existing simulations, we may want to address the
# previous L's to find an initial Tc but not actually re-run those simulation.
# those can be but in the arrayLexclude array
arrayLexclude=()
arrayH=( 0.3 )
minT_false=1.3
maxT_false=1.8
minT_true=1.6
maxT_true=2.0
delta=0.1
elif [ $SYS_NAME == "Fe8" ]; then
name="Fe8_test4"
arrayMech=("true")
arrayL=( 5 6 7 8 )
arrayLexclude=()
arrayH=( 0.0 )
minT_false=1.4
maxT_false=2.0
minT_true=0.3
maxT_true=0.9
delta=0.08
fi


gen_temp_schedules() {
  # generate temporary temperature schedules
  for L in "${arrayL[@]}"; do
    if [[ ! " ${arrayLexclude[@]} " =~ " ${L} " ]]; then
      for H in "${arrayH[@]}"; do

        yes n | /gpfs0/smoshe/projects/Python-3.8_old/bin/python3 ./scripts/gen_temp_schedule.py "$1" "$2" $nT | head -n1 >./"$SYS_NAME"/temperature_schedules/temp_schedule_"$L"_"$L"_"$4""$name"_"$H"_"$3".txt
        echo "generated ./${SYS_NAME}/temperature_schedules/temp_schedule_${L}_${L}_$4${name}_${H}_$3.txt"

      done
    fi
  done
}

# enables one to temporarily deactivate parts of this script when set to false.
# notice where the terminating "fi" is.
if true; then
# generate temporary temperature schedules for "true"
gen_temp_schedules $minT_true $maxT_true "true" "temp_"
# generate temporary temperature schedules for "false"
gen_temp_schedules $minT_false $maxT_false "false" "temp_"

# run initial temporary run
sub "temp_" 600 20 $nT


# wait until temporary runs finish
running_jobs=$(qstat -u tomerdol | awk '$3 ~ /^tr/' | awk 'END {print NR}')
while [ $running_jobs -gt 0 ]; do
  echo "there are still $running_jobs jobs running. waiting 1 hour."
  sleep 1h
  running_jobs=$(qstat -u tomerdol | awk '$3 ~ /^tr/' | awk 'END {print NR}')
done

#*************************************************************************
# check that there are results
#*************************************************************************
echo "checking that there are results from temporary runs"

for mech in "${arrayMech[@]}"; do
  for L in "${arrayL[@]}"; do
    for H in "${arrayH[@]}"; do
      # read back temperatures into temp_schedule array
      IFS=$'\n,' read -d '' -ra temp_schedule <./"$SYS_NAME"/temperature_schedules/temp_schedule_"$L"_"$L"_temp_"$name"_"$H"_"$mech".txt

      for T in "${temp_schedule[@]}"; do
        for file in ./"$SYS_NAME"/data/results/temp_"$name"/table_"$L"_"$L"_"$H"_"$T"_"$mech"_*.txt; do
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

# extract energy and magnetization data from each simulation
#*************************************************************************
echo "extracting data from temporary runs"
for mech in "${arrayMech[@]}"; do
  for L in "${arrayL[@]}"; do
    if [[ ! " ${arrayLexclude[@]} " =~ " ${L} " ]]; then
      for H in "${arrayH[@]}"; do

        IFS=$'\n,' read -d '' -ra temp_schedule <./"$SYS_NAME"/temperature_schedules/temp_schedule_"$L"_"$L"_temp_"$name"_"$H"_"$mech".txt

        echo "reading temperatures from ./${SYS_NAME}/temperature_schedules/temp_schedule_${L}_${L}_temp_${name}_${H}_${mech}.txt"
        echo "T e corr_length"

        tmp_file="sample_energy.txt"
        if [ -f $tmp_file ]; then
          rm $tmp_file
        fi
        echo T e corr_length >>$tmp_file

        for T in "${temp_schedule[@]}"; do

          file="./"$SYS_NAME"/data/results/temp_""$name""/table_""$L""_""$L""_""$H""_""$T""_""$mech""_*.txt"

          # read the data from all seeds
          num_of_samples=$(tail -q -n 1 ${file} | awk 'END{print $1}')
          # use just the last quarter of the results
          start_point=$((num_of_samples / 4))
          energy=$(tail -q -n ${start_point} ${file} | awk '!/^ *#/{ total += $3; count++ } END { print total/count }')
          m2=$(tail -q -n ${start_point} ${file} | awk '!/^ *#/{ total += ($2)^2; count++ } END { print total/count }')
          mk2=$(tail -q -n ${start_point} ${file} | awk '!/^ *#/{ total += $15; count++ } END { print total/count }')
          # calculate the correlation length
          corr_length=$(bc -l <<<"scale=5; sqrt(($m2/$mk2)-1)/(2*$L*s(4*a(1)/$L))")

          echo $T $energy $corr_length >>$tmp_file
          echo $T $energy $corr_length
        done

        # print the table of T, Energy, correlation length to a file in /data/analysis/
        sort -t" " -nk1 $tmp_file >./"$SYS_NAME"/data/analysis/sample_energy_"$L"_"$H"_temp_"$name"_"$mech".txt
        rm $tmp_file

      done
    fi
  done
done

#*************************************************************************

# find initial Tc and create new temperature schedules
#*************************************************************************
echo "looking for initial Tc and creating new temperature schedules"
for mech in "${arrayMech[@]}"; do
  for H in "${arrayH[@]}"; do

    all_L="${arrayL[@]}"
    # find the an initial guess for Tc based on the short temporary runs
    initial_tc=$(/gpfs0/smoshe/projects/Python-3.8_old/bin/python3 ./scripts/find_initial_tc.py -sys ${SYS_NAME} -h_ex ${H} --mech ${mech} --folder_list temp_${name} -L ${all_L})
    echo "Hex=$H mech=$mech : initial Tc = $initial_tc"

    # find an optimal temperature schedule based on the energies from the short temporary run
    for L in "${arrayL[@]}"; do
      if [[ ! " ${arrayLexclude[@]} " =~ " ${L} " ]]; then
        bash ./scripts/temp_set_script.sh ./"$SYS_NAME"/data/analysis/sample_energy_"$L"_"$H"_temp_"$name"_"$mech".txt $nT $(bc <<<"$initial_tc-$delta") $(bc <<<"$initial_tc+$delta") | awk -vORS=, '{ print $1 }' | sed 's/,$/\n/' >./"$SYS_NAME"/temperature_schedules/temp_schedule_"$L"_"$L"_"$name"_"$H"_"$mech".txt
        echo "created ./${SYS_NAME}/temperature_schedules/temp_schedule_${L}_${L}_${name}_${H}_${mech}.txt"
      fi
    done

  done
done
#*************************************************************************

fi
# run simulations
#echo "running simulations"
#sub "" 2048 50 $nT
echo "done"
exit 0
