#!/bin/bash
# Check the progress of all jobs currently running under user tomerdol.
# Reads the last index from each job's results file and compares to the total
# number of sweeps it was given as a parameter.
# The output is formatted like the output of qstat with an additional column "progress"
# Uses the check_prog.sh script.
# Run from the root project directory as bash ./scripts/check_overall_prog.sh

tmp_file="progress_tracking.txt"
if [ -f $tmp_file ]; then
rm $tmp_file
fi

total=0
count=0
echo "job-ID  prior   name       user         state submit/start at     queue                          slots  progress
----------------------------------------------------------------------------------------------------------------------"
qstat -u tomerdol | tail -n +3 | (while read line; do
jobid=$(echo $line | awk '{print $1}' )
output=$(bash scripts/check_prog.sh $jobid | tail -n1 | awk '{print $NF}' )
((count+=$(echo $output | cut -d/ -f1 )))
((total+=$(echo $output | cut -d/ -f2 )))
echo "$line""    ""$output"
echo "$line""  ""$(echo $output | cut -d/ -f1 )" >> $tmp_file
done
echo "$(echo "100*${count}/${total}" | bc -l)%")
