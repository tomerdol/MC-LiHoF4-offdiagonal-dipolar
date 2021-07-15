#!/bin/bash
# Run all simulations under the given project name with the same number of sweeps they were originally ran with.
# this was originally used to plot the final simulation states for finished simulations when this feature was
# introduced.
# Can still be used to introduce new features to finished simulations.

# Check input
if [ $# -ne 1 ]
then
echo Usage $0 [directory]
exit
fi

name=$1
# loop through all the checkpoint files
for file in ./"$SYS_NAME"/checkpoints/"$name"/save_state_*.txt
do
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

echo "Fixing ${filename}."
# get the number of sweeps that were run
maxSweeps=$(ls ./"$SYS_NAME"/data/results/"$name"/table_"$Lx"_"$Lz"_"$H"_*_"$mech"_"$seed".txt | head -n1 | xargs grep -v '^#' | tail -n1 | awk '{print $1}')
# run simulation with same maxSweeps as before so that nothing happens except printing the lattice
java -Dsystem="$SYS_NAME" -classpath "./bin/production/LiHoF4_transverse_field_MC/:./lib/*" simulation.montecarlo.Main -max_iter 80 -L "$Lx","$Lz" -d 1.0 -extBx "$H" -continue_from_save yes -max_sweeps "$maxSweeps" -name "$name" -temp_schedule ./"$SYS_NAME"/temperature_schedules/temp_schedule_"$Lx"_"$Lz"_"$name"_"$H"_"$mech".txt -alpha 0.95 -verbose -mode s -seed "$seed" -Jex 1.16e-3 "$suppress_tag"

done
