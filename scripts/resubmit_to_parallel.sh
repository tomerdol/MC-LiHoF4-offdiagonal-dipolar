#!/bin/bash

num_of_queued_jobs=$(qstat -u tomerdol | awk '$9=="1" {count++} END {print count}')

count=0
while [[ $num_of_queued_jobs -gt 0 && count -lt 100 ]]
do
jobid=$(qstat -u tomerdol | awk '$9=="1" {print $1; exit}')
args=$(qstat -j ${jobid} | grep job_args | awk '{print $2}')
jobname=$(qstat -j ${jobid} | grep job_name | awk '{print $2}')

# resubmit to smoshe.q,lublin.q in parallel mode
qdel "$jobid"
qsub -pe shared 24 -l mem_free=40G -V -S /bin/bash -cwd -N "$jobname" -o ./"$temp"output/ -e ./"$temp"output/ -q smoshe.q,lublin.q met_with_t.sh "${args//,/ }"
((count++))
sleep 30
num_of_queued_jobs=$(qstat -u tomerdol | awk '$9=="1" {count++} END {print count}')
done
