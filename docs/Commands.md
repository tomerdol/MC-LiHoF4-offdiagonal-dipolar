## Some useful commands for working with the SGE cluster

* Monitor running and queued jobs:
  ```
  qstat -u username
  ```
* Submit single job:
  ```
  qsub -cwd -q smoshe.q@sge247 -V -S /bin/bash -N "phase_diagram" -hold_jid 9918867
  ```
  -cwd: current working directory\
  -q: queue list\
  -V: pass environment variables to job\
  -S: shell to use\
  -N: job name\
  -hold_jid: wait for another job to finish before starting
* Start interactive session:
  ```
  qrsh -q smoshe.q
  ```
* Kill all jobs by some criteria:
  ```
  qstat -u tomerdol | awk '{if($3!="QRLOGIN" && NR>2 && $3 ~ /^tr7_0.0_/){print $1}}' > Worklist
  cat Worklist | xargs qdel
  rm Worklist
  ```
* Standard run:
  ```
  java -Dsystem=LiHoF4 -classpath "./bin/production/LiHoF4_transverse_field_MC/:./lib/*" simulation.montecarlo.Main -max_iter 80 -L 5,5 -d 0.6 -extBx 0.0 -continue_from_save yes -max_sweeps 513 -name dilution_test -temp_schedule ./"$SYS_NAME"/temperature_schedules/temp_schedule_5_5_dilution_test_0.0_false.txt -alpha 0.95 -verbose -mode s -seed 123456789 -Jex 1.16e-3 -interpolation_table_name _0.014 -p -print_sweep_num 1
  ```
* Format T<sub>c</sub> output for copy-paste into numpy array:
  ```
  cat phase_diagram_true_4_5_6_7_res_test_ex_0.014_res.txt | awk 'NR>1{print "["$1", " $2", " $3"],"}'
  ```
* Find smallest absolute value Bz fields in lattice_output:
  ```
  awk 'NR==1 {min=100} {if (/^[^#]/ && FNR>10 && ($6<0?$6*-1:$6)<min) {min=($6<0?$6*-1:$6); print min}} END {print min}' *
  ```
* Create phase diagram:
  ```
  python3 phase_diagram_bin.py --boot_num 500 -h_ex 0.0 --mech true --folder_list res_test_ex_0.014 -L 5 6 7 -sys LiHoF4
  ```
* Continue existing simulation:
  ```
  echo 'bash ./scripts/cont_sim.sh 16400 _0.014' | qsub -cwd -q smoshe.q,lublin.q -V -S /bin/bash -N "cont_sim"
  ```
* Delete (empty) SGE output files:
  ```
  find . -maxdepth 1 -type f -empty -name '*.[o,e][0-9][0-9][0-9][0-9][0-9][0-9]*' -delete
  find . -maxdepth 1 -type f -empty -name '*.p[o,e][0-9][0-9][0-9][0-9][0-9][0-9]*' -delete
  ```
  * Use -mtime +1 to only delete files that are 1 day old or older. 
* Archive lattice output (run from ./LiHoF4/data/lattice_output/project_directory/)
  ```
  tar -czvf archived_lattice_output.tar.gz --remove-files *.txt
  ```
* Re-pack .h5 file (save space and possibly improve i/o speed):
  ```
  ptrepack -o --chunkshape=auto --propindexes  table_5_5_0.0_false.h5 table_5_5_0.0_false2.h5
  ```
    * Ptrepack is installed with pytables and in located in the python /bin
    * After finished, we can remove table_5_5_0.0_false.h5 and rename table_5_5_0.0_false2.h5 to table_5_5_0.0_false.h5
* Sum statistics on methods used for self-consistent calculation:
  ```
  grep "# Used methods statistics:" *.o* | awk -F "[|]|," '{ for (i=1; i<=NF; i++) total[i] += $i; } END { for (i=1; i<=NF; i++) print i "=" total[i]+0 }'
  ```
* Restore from ESS backup (overwrite existing files):
  ```
  echo 'rsync -vah --ignore-times /gpfs0/smoshe/.snapshots/daily-2021.07.04-22.30.54/users/tomerdol/home/LiHoF4_transverse_field_MC/dev/ /gpfs0/smoshe/users/tomerdol/home/LiHoF4_transverse_field_MC/dev/' | qsub -cwd -q smoshe.q -V -S /bin/bash -N "restr_bck"
  ```
* kill vncserver in one line:
  ```
  kill `ps -ef | grep -i tomerdol | grep vnc | awk '$8=="/usr/bin/Xvnc"' | awk '{print $2}' `
  ```