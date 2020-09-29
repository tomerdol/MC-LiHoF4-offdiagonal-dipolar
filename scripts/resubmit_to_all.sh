#!/bin/bash

num_of_queued_jobs=$(qstat -u tomerdol | grep r4_ | awk '$5=="qw" {count++} END {print count}')

count=0
while [[ $num_of_queued_jobs -gt 0 && count -lt 40 ]]
do
    jobid=$(qstat -u tomerdol | grep qw | grep r4_ | head -n1 | awk '{print $1}')
    args=$(qstat -j ${jobid} | grep job_args | awk '{print $2}')
    jobname=$(qstat -j ${jobid} | grep job_name | awk '{print $2}')

    if [[ ${jobname:0:1} == "t" ]]; then temp="temp_"; else temp=""; fi

    # resubmit to intel_all.q in single mode
    qsub -l mem_free=40G -V -S /bin/bash -cwd -N "$jobname" -o ./output/ -e ./output/ -q smoshe.q@sge1082,lublin.q,intel_all.q scripts/met_with_t_single.sh "${args//,/ }"
    qdel "$jobid"

    ((count++))
    sleep 10

    num_of_queued_jobs=$(qstat -u tomerdol | grep r4_ | awk '$5=="qw" {count++} END {print count}')
done
