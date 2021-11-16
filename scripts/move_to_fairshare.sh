#!/bin/bash

while :
do
num_of_running_jobs=$(qstat -u tomerdol | grep fairshare.q | awk 'BEGIN {sum=0} ; {sum+=$9} ; END {print sum}')
while [[ $num_of_running_jobs -lt 460 ]]
do
    job_id=$(qstat -u tomerdol | awk '{if ($5=="qw" && $8==24){print $1;exit;}}')
    qalter -q smoshe.q@sge247,smoshe.q@sge1081,smoshe.q@sge1082,lublin.q,fairshare.q "$job_id"
    sleep 20
    # make sure the altered job is running before moving on
    while [ "$(qstat -u tomerdol | awk -v jid="$job_id" '$1 == jid {print $5}')" != "r" ]; do
        sleep 1h
    done

    num_of_running_jobs=$(qstat -u tomerdol | grep fairshare.q | awk 'BEGIN {sum=0} ; {sum+=$9} ; END {print sum}')
done
sleep 10
done
