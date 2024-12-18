# Run a serial MC simulation
# Usage: $0 [Lx] [Lz] [# of sweeps] [external Bx] [mechanism flag] [project name] [seed] [interpolation table file name extension] [dilution]
# $1=Lx $2=Lz $3=max_sweeps $4=extBx $5=mech $6=name $7=seed $8=interpolation_table_name $9=x

if [ $SYS_NAME == "LiHoF4" ]; then
#jex="1.16e-3"
#jex="3.91e-3"
#jex="0.0"
jex="${10}"
x=$9
elif [ $SYS_NAME == "Fe8" ]; then
jex="0.0"
x="1.0"
fi
if [ $5 == "true" ]; then
suppress_tag="-suppress"
else
suppress_tag=""
fi
# echo the seed
echo $7
# run the simulation
java -Dsystem="$SYS_NAME" -classpath "./bin/production/LiHoF4_transverse_field_MC/:./lib/*" simulation.montecarlo.Main -max_iter 80 -L "$1","$2" -d "$9" -extBx "$4" -continue_from_save yes -max_sweeps "$3" -name "$6" -temp_schedule ./"$SYS_NAME"/temperature_schedules/temp_schedule_"$1"_"$2"_"$6"_"$4"_"$5".txt -alpha 0.95 -verbose -mode s -seed "$7" -Jex "$jex" "$suppress_tag" -interpolation_table_name "$8"

echo "done"
