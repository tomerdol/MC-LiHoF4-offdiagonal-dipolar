#!/bin/bash
# Removes queued parallel jobs and resubmits them to intel_all.q in serial mode

# Check input
if [ $# -ne 2 ]
then
echo "Usage $0 [jobs name] [# of jobs to resubmit]"
exit
fi

search_key=$1
num_of_jobs_to_resubmit=$2
num_of_queued_jobs=$(qstat -u tomerdol | grep "$search_key" | awk '$5=="qw" {count++} END {print count}')

count=0
while [[ $num_of_queued_jobs -gt 0 && count -lt "$num_of_jobs_to_resubmit" ]]
do
    jobid=$(qstat -u tomerdol | grep qw | grep "$search_key" | head -n1 | awk '{print $1}')
    args=$(qstat -j ${jobid} | grep job_args | awk '{print $2}')
    jobname=$(qstat -j ${jobid} | grep job_name | awk '{print $2}')

    if [[ ${jobname:0:1} == "t" ]]; then temp="temp_"; else temp=""; fi

    # resubmit to intel_all.q in single mode
    new_jobid=$(qsub -l mem_free=40G -V -S /bin/bash -cwd -N "$jobname" -o ./"$SYS_NAME"/output/ -e ./"$SYS_NAME"/output/ -q smoshe.q@sge1082,lublin.q,intel_all.q scripts/met_with_t_single.sh "${args//,/ }")
    echo "$new_jobid"  # this echos the submission message
    # get the job id number from the output message
    new_jobid=$(echo $new_jobid | awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}')
    qdel "$jobid"
    
    sleep 8

    # make sure the new submitted job is running before grabbing the next queued job
    while [ "$(qstat -u tomerdol | awk -v jid="$new_jobid" '$1 == jid {print $5}')" != "r" ]; do
      echo "New job ($new_jobid) is not yet running. Sleeping for 10 seconds."
      sleep 10
    done

    ((count++))

    num_of_queued_jobs=$(qstat -u tomerdol | grep "$search_key" | awk '$5=="qw" {count++} END {print count}')
done
