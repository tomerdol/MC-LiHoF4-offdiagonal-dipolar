#!/bin/bash
# Continuously monitors the queues and periodically moves the slowest serial job to a parallel run on the private queues

nT=24

#bash ./scripts/check_overall_prog.sh

while :
do
while read i
do
    echo "moving job ${i}"
    job_id=$(echo $i | awk '{print $1}')
    
    if qstat -u tomerdol | grep -q -w "$job_id"; then
        # job is in the queue
        args=$(qstat -j ${job_id} | grep job_args | awk '{print $2}')
        jobname=$(qstat -j ${job_id} | grep job_name | awk '{print $2}')
         # resubmit to smoshe.q,lublin.q in parallel mode
        qdel "$job_id"
        queues="lublin.q,smoshe.q@sge1081,smoshe.q@sge1082,smoshe.q@sge190,smoshe.q@sge247,smoshe.q@sge249"
        new_jobid=$(qsub -pe shared "$nT" -l mem_free=40G -V -S /bin/bash -cwd -N "$jobname" -o ./"$SYS_NAME"/output/ -e ./"$SYS_NAME"/output/ -q "$queues" scripts/met_with_t.sh "${args//,/ }")

        echo "$new_jobid"  # this echos the submission message
        # get the job id number from the output message
        new_jobid=$(echo $new_jobid | awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}')

        sleep 20

        # make sure the new submitted job is running before grabbing the next queued job
        while [ "$(qstat -u tomerdol | awk -v jid="$new_jobid" '$1 == jid {print $5}')" != "r" ]; do
            echo "New job ($new_jobid) is not yet running. Sleeping for 1 hour."
            sleep 1h
        done

    else
        echo "job no longer running"
    fi

done < <(cat progress_tracking.txt | awk '{if($3 ~ /^r[0-9]/ &&  $5=="r" && $9==1){print $1" "$10} else if($3 ~ /^r[0-9]/ && $5="qw" && $8==1){print $1" "$9}}' | sort -k 2,2 | head -n 5)

qsub -cwd -q smoshe.q,lublin.q -V -S /bin/bash -N "track_prog" -o /dev/null ./scripts/check_overall_prog.sh
sleep 3h
done
