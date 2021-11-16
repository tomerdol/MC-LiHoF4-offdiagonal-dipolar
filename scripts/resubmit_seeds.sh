#!/bin/bash

input="scripts/simulations.txt"
i=0
while IFS= read -r line
do
  qsub -l mem_free=2G -V -S /bin/bash -cwd -N r7_0.0_"$i"_f -o ./"$SYS_NAME"/output/ -e ./"$SYS_NAME"/output/ -q smoshe.q@sge247,smoshe.q@sge1081,smoshe.q@sge1082,lublin.q,intel_all.q scripts/met_with_t_single.sh 7 7 8200 0.0 false low_jex_dilution_0.46 "$line" _const 0.46 1.16e-3
  ((i++))
  sleep 10
done < "$input"
