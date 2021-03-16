#!/bin/sh

#arrayH=( 1.5 2.0 )
#arrayH=( 0.0 0.04081633 0.08163265 0.12244898 0.16326531 0.20408163 0.24489796 0.28571429 0.32653061 0.36734694 0.40816327 0.44897959 0.48979592 0.53061224 0.57142857 0.6122449 0.65306122 0.69387755 0.73469388 0.7755102 0.81632653 0.85714286 0.89795918 0.93877551 0.97959184 1.02040816 1.06122449 1.10204082 1.14285714 1.18367347 1.2244898 1.26530612 1.30612245 1.34693878 1.3877551 1.42857143 1.46938776 1.51020408 1.55102041 1.59183673 1.63265306 1.67346939 1.71428571 1.75510204 1.79591837 1.83673469 1.87755102 1.91836735 1.95918367 2.0 )
arrayH=( 0.0 0.5 1.0 1.5 2.0 )

for H in "${arrayH[@]}"
do
if [ $SYS_NAME == "LiHoF4" ]; then
    cosphi=1
    sinphi=0
elif [ $SYS_NAME == "Fe8" ]; then
    cosphi=0.375
    sinphi=0.927
fi

#qsub -V -S /bin/bash -cwd -N table"$H" -q smoshe.q@sge1081,smoshe.q@sge1082,smoshe.q@sge247,smoshe.q@sge249 scripts/create_one_table.sh "$H"
qsub -V -S /bin/bash -cwd -N table"$H" -q smoshe.q@sge1081,smoshe.q@sge1082,smoshe.q@sge247,smoshe.q@sge249 scripts/create_one_table.sh $(bc <<< "$H*$cosphi") $(bc <<< "$H*$sinphi")
sleep 30
done
