#!/bin/bash

num_of_queued_jobs=$(qstat -u tomerdol | awk '$5=="qw" {count++} END {print count}')

count=0
while [[ $num_of_queued_jobs -gt 20  && count -lt 22 ]]
do
    jobid=$(qstat -u tomerdol | grep qw | head -n1 | awk '{print $1}')
    args=$(qstat -j ${jobid} | grep job_args | awk '{print $2}')
    jobname=$(qstat -j ${jobid} | grep job_name | awk '{print $2}')

    if [[ ${jobname:0:1} == "t" ]]; then temp="temp_"; else temp=""; fi

    # resubmit to intel_all.q in single mode
    qsub -l mem_free=40G -V -S /bin/bash -cwd -N "$jobname" -o ./output/ -e ./output/ -q smoshe.q,lublin.q,intel_all.q met_with_t_single.sh "${args//,/ }"
    qdel "$jobid"

    ((count++))
    sleep 30
    num_of_queued_jobs=$(qstat -u tomerdol | awk '$5=="qw" {count++} END {print count}')
done
