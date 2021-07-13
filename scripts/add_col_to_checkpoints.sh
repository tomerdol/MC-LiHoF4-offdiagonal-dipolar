#!/bin/bash
#
# Add a new column to an existing checkpoint file

if [ $# -ne 1 ]
then
echo Usage $0 [directory]
exit
fi

name=$1
# loop through all checkpoints in the given directory
for file in ./"$SYS_NAME"/checkpoints/"$name"/save_state_*.txt
do
# extract parameters from file names
filename=${file##*/}
name_no_ext="${filename%.txt}"
Lx=$(echo "${name_no_ext}" | awk -F_ '{print $3}')
Lz=$(echo ${name_no_ext} | awk -F_ '{print $4}')
H=$(echo ${name_no_ext} | awk -F_ '{print $5}')
mech=$(echo ${name_no_ext} | awk -F_ '{print $6}')
seed=$(echo ${name_no_ext} | awk -F_ '{print $7}')
if [ $mech == "true" ]; then
suppress_tag="-suppress"
else
suppress_tag=""
fi

# run fixBackwardsCompatibility with the new number of columns (here: 36)
echo "Fixing ${filename}."
java -Dsystem="$SYS_NAME" -classpath "./bin/production/LiHoF4_transverse_field_MC/:./lib/*" utilities.fixBackwardsCompatibility -max_iter 0 -max_sweeps 0 -s yes -mode s -L "$Lx","$Lz" -extBx "$H" -name "$name" -seed "$seed" -new_col_num 36 "$suppress_tag"

done
