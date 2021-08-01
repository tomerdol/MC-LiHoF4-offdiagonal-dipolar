#!/bin/bash
# Generate a temperature schedule based on the results of an old project
# run from /scripts

# Check input
if [ $# -lt 8 ]
then
echo "Usage $0 [Tc] [delta] [# of temperatures] [dev|beta|stable] [project_name] [mech] [Bex] {L (sorted)}"
exit
fi

initial_tc="$1"
delta="$2"
nT="$3"
root_dir="$4"
name="$5"
mech="$6"
H="$7"
arrayL=(${@:8})

# extract energy and magnetization data from each simulation
#*************************************************************************
echo "extracting data from temporary runs"
for L in "${arrayL[@]}"; do
  if [ -f "../../"$root_dir"/"$SYS_NAME"/temperature_schedules/temp_schedule_"$L"_"$L"_"$name"_"$H"_"$mech".txt" ]; then
    maxL="$L"
    IFS=$'\n,' read -d '' -ra temp_schedule < ../../"$root_dir"/"$SYS_NAME"/temperature_schedules/temp_schedule_"$L"_"$L"_"$name"_"$H"_"$mech".txt

    echo "reading temperatures from ../../${root_dir}/${SYS_NAME}/temperature_schedules/temp_schedule_${L}_${L}_temp_${name}_${H}_${mech}.txt"
    echo "T e corr_length"

    tmp_file="sample_energy.txt"
    if [ -f $tmp_file ]; then
      rm $tmp_file
    fi
    echo T e corr_length >>$tmp_file

    for T in "${temp_schedule[@]}"; do
      # remove trailing zeros in T
      T=$(echo ${T} | sed -e 's/[0]*$//g')
      file="../../""$root_dir""/""$SYS_NAME""/data/results/""$name""/binned_data/table_""$L""_""$L""_""$H""_""$T""_""$mech""_*.txt"

      if compgen -G "$file" > /dev/null; then
        # files of the form table_""$L""_""$L""_""$H""_""$T""_""$mech""_*.txt" exist
        # this means that the binned data is in regular text files

        # read the data from all seeds
        energy=$(tail -q -n 1 ${file} | awk '!/^ *#/{ total += $2; count++ } END { print total/count }')
        m2=$(tail -q -n 1 ${file} | awk '!/^ *#/{ total += $4; count++ } END { print total/count }')
        mk2=$(tail -q -n 1 ${file} | awk '!/^ *#/{ total += $17; count++ } END { print total/count }')
        # calculate the correlation length
        corr_length=$(bc -l <<<"scale=5; sqrt(($m2/$mk2)-1)/(2*$L*s(4*a(1)/$L))")
      else
        # this means that the binned data is in hdf5 files

        # save original working directory
        original_dir=$(pwd)
        # change working dir to correct stage (dev|stable)
        cd "../../""$root_dir""/scripts"
        # indentation is important here. run python script to
output=$(
python3 - <<END_SCRIPT
import bin_data
import analysis_tools
sim = analysis_tools.get_simulation("$L", "$name", "$H", "$mech", "$T")
df = bin_data.read_binned_data(sim, use_latest=False)
print(df.Energy.mean(),df['Magnetization^2'].mean(),df['mk2x'].mean())
END_SCRIPT
)
        energy="${output[1]}"
        m2="${output[2]}"
        mk2="${output[3]}"
        corr_length=$(bc -l <<<"scale=5; sqrt(($m2/$mk2)-1)/(2*$L*s(4*a(1)/$L))")

        cd "$original_dir"
      fi

      echo $T $energy $corr_length >>$tmp_file
      echo $T $energy $corr_length
    done

    # print the table of T, Energy, correlation length to a file in /data/analysis/
    sort -t" " -nk1 $tmp_file >../"$SYS_NAME"/data/analysis/sample_energy_"$L"_"$H"_"$name"_"$mech".txt
    rm $tmp_file

    echo "L: ${L}"
    new_temp_schedule=$(bash ./temp_set_script.sh ../"$SYS_NAME"/data/analysis/sample_energy_"$L"_"$H"_"$name"_"$mech".txt $nT $(bc <<<"$initial_tc-$delta") $(bc <<<"$initial_tc+$delta") | awk -vORS=, '{ print $1 }' | sed 's/,$/\n/')
    echo "$new_temp_schedule"
    read -p "Write to temperature schedule file? (Y/N): " confirm
    if [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]]; then
      read -p "new project name: " new_name
      echo "$new_temp_schedule" > ../"$SYS_NAME"/temperature_schedules/temp_schedule_"$L"_"$L"_"$new_name"_"$H"_"$mech".txt
      echo "Writing ../"$SYS_NAME"/temperature_schedules/temp_schedule_"$L"_"$L"_"$new_name"_"$H"_"$mech".txt"
    fi
  else
    echo "no temperature schedule exists for this L: ${L}"
    if [ -z ${maxL+x} ]; then
      echo "Lowest given L does not exist in given project. Make sure the L's are sorted from lowest to highest."
    else
      read -p "Copy L=${maxL} temperature schedule? (Y/N): " confirm
      # copy the temperature schedule from the maximum L to use with this one as well
      if [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]]; then
        cp ../"$SYS_NAME"/temperature_schedules/temp_schedule_"$maxL"_"$maxL"_"$new_name"_"$H"_"$mech".txt ../"$SYS_NAME"/temperature_schedules/temp_schedule_"$L"_"$L"_"$new_name"_"$H"_"$mech".txt
      fi
    fi
  fi
done