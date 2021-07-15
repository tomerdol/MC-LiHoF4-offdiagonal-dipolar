#!/bin/bash
# Deletes serial jobs and resubmits them to smoshe.q,lublin.q in parallel mode

# Check input
if [ $# -ne 3 ]
then
echo "Usage $0 [job name] [# of jobs to resubmit] [# of temperatures]"
exit
fi

search_key=$1
num_of_jobs_to_resubmit=$2
nT=$3
num_of_queued_jobs=$(qstat -u tomerdol | grep "$search_key" | awk '$9=="1" {count++} END {print count}')

count=0
while [[ $num_of_queued_jobs -gt 0 && count -lt "$num_of_jobs_to_resubmit" ]]
do
jobid=$(qstat -u tomerdol | grep tr9_ | awk '$9=="1" {print $1; exit}')
args=$(qstat -j ${jobid} | grep job_args | awk '{print $2}')
jobname=$(qstat -j ${jobid} | grep job_name | awk '{print $2}')

# resubmit to smoshe.q,lublin.q in parallel mode
qdel "$jobid"
qsub -pe shared "$nT" -l mem_free=40G -V -S /bin/bash -cwd -N "$jobname" -o ./"$SYS_NAME"/output/ -e ./"$SYS_NAME"/output/ -q smoshe.q,lublin.q scripts/met_with_t.sh "${args//,/ }"
((count++))
sleep 10
num_of_queued_jobs=$(qstat -u tomerdol | grep "$search_key" | awk '$9=="1" {count++} END {print count}')
done
