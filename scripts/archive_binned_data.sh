#!/bin/bash

if [ $# -ne 1 ]
then
echo Usage $0 [directory]
exit
fi

dir="../${SYS_NAME}/data/results/""$1"

tmp_data_file="$dir""/files_to_archive.tmp"
if [ -f $tmp_data_file ]
then
rm $tmp_data_file
fi

checkpoint_dir="../${SYS_NAME}/checkpoints/""$1"
tmp_checkpoints_file="$checkpoint_dir""/files_to_archive.tmp"
if [ -f $tmp_checkpoints_file ]
then
rm $tmp_checkpoints_file
touch $tmp_checkpoints_file
fi

# list all files eligible to be archived
for file in ${dir}/*
do
echo -n "Checking ${file} ... "
if [ ! -d "$file" ]; then
    # find matching checkpoint to also archive
    filename=${file##*/}
    name_no_ext="${filename%.txt}"
    Lx=$(echo "${name_no_ext}" | awk -F_ '{print $2}')
    Lz=$(echo ${name_no_ext} | awk -F_ '{print $3}')
    H=$(echo ${name_no_ext} | awk -F_ '{print $4}')
    T=$(echo ${name_no_ext} | awk -F_ '{print $5}' | tr . _)
    mech=$(echo ${name_no_ext} | awk -F_ '{print $6}')
    seed=$(echo ${name_no_ext} | awk -F_ '{print $7}')
    # remove tamperature and seed to get the .h5 file name
    binned_file="${dir}/binned_data/table_${Lx}_${Lz}_${H}_${mech}.h5"
    if [ -f "$binned_file" ]; then
        # first try to find the last sample in the last 100 lines
        # this is usually the case, unless a finished simulation
        # was re-started many time and the headers fill > 100 lines
        last_sample=$(tail -n100 "$file" | sed '/^#/d' | tail -n1 | awk '{print $1}')
        # only if last_sample turns out empty, look in the entire file
        if [ -z "$last_sample" ]; then
            last_sample=$(sed '/^#/d' "$file" | tail -n1 | awk '{print $1}')
        fi
        # if h5ls fails this means that the dataset does not exist so just put 0
        last_bin=$(h5ls "${binned_file}/s${seed}/T${T}" | grep -oP '{\K.*?(?=/)' || echo "0")
        last_bin=$(($last_bin-1))
        # both are numeric
        if ([[ $last_sample =~ ^[[:digit:]] ]] && [[ $last_bin =~ ^[[:digit:]] ]]); then
            last_sample_for_next_bin=$(echo "2*(2^($last_bin + 1) - 1)" | bc)
            # last written data is already binned
            if [ "$last_sample_for_next_bin" -gt "$last_sample" ]; then
                # we can add this file to the archive
                basename "$file" >> "$tmp_data_file"
                echo " binned."
                checkpoint_file="save_state_${Lx}_${Lz}_${H}_${mech}_${seed}.txt"
                if [ -f "$checkpoint_dir"/"$checkpoint_file" ]; then
                    grep -qxF "$checkpoint_file" "$tmp_checkpoints_file" || echo "$checkpoint_file" >> $tmp_checkpoints_file
                    #echo "$checkpoint_file" >> "$tmp_data_file"
                else
                    echo "checkpoint file does not exist."
                    exit 1
                fi
            else
                echo "has unbinned data: $last_sample_for_next_bin, $last_sample" 
                exit 1
            fi
        else
            echo "cound not find last indices (non-numeric)."
            exit 1
        fi
    else
        echo "no binned data file."
        exit 1
    fi
else
    echo "not a file."
fi
done

echo "archiving files..."
# archive eligible files
tar -czvf archived_results.tar.gz -C "$dir" -T "$tmp_data_file" --remove-files
echo "archiving checkpoints..."
tar -czvf archived_checkpoints.tar.gz -C "$checkpoint_dir" -T "$tmp_checkpoints_file" --remove-files

rm $tmp_data_file
rm $tmp_checkpoints_file

echo "done"
