#!/bin/bash
#parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

#cd "$parent_path"

total=0
count=0
echo "job-ID  prior   name       user         state submit/start at     queue                          slots  progress
----------------------------------------------------------------------------------------------------------------------"
qstat -u tomerdol | tail -n +3 | (while read line; do
jobid=$(echo $line | awk '{print $1}' )
output=$(bash scripts/check_prog.sh $jobid | tail -n1 | awk '{print $NF}')
((count+=$(echo $output | cut -d/ -f1 )))
((total+=$(echo $output | cut -d/ -f2 )))
echo "$line""    ""$output"
done
echo "$(echo "100*${count}/${total}" | bc -l)%")
