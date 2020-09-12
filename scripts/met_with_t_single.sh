# $1=Lx $2=Lz $3=max_sweeps $4=extBx $5=mech $6=name $7=seed
if [ $5 == "true" ]; then
suppress_tag="-suppress"
else
suppress_tag=""
fi
echo $7
java -classpath "./out/production/project7/:./lib/*" simulation.montecarlo.Main -max_iter 80 -L "$1","$2" -d 1.0 -extBx "$4" -continue_from_save yes -max_sweeps "$3" -name "$6" -temp_schedule ./temperature_schedules/temp_schedule_"$1"_"$2"_"$6"_"$4"_"$5".txt -alpha 0.95 -verbose -mode s -seed "$7" -Jex 1.16e-3 "$suppress_tag" -dummy 0

echo "done"
