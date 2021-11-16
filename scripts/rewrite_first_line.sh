#!/bin/bash

jex="1.16e-3"
name="low_jex_dilution_1.0"
x="1.0"
H="0.0"
L=6
mech="false"
files="../LiHoF4/data/results/""$name""/table_${L}_${L}_${H}_*_${mech}_*"

if [ "$mech" == "true" ]; then
suppress_tag="-suppress"
else
suppress_tag=""
fi

awk '{if ($0 !~ /^[[:space:]]*#/ && $1 != "index") {if ($1==1) print FILENAME; nextfile}}' ../LiHoF4/data/results/"$name"/table_"$L"_"$L"_* | awk -F "[_.]" '{print $15}' > tmp_file.txt
sort tmp_file.txt | uniq > seeds.txt
rm tmp_file.txt

echo "found all seeds"
while read -r seed; do
    echo "Seed: ${seed}"
    cd ..
    java -Dsystem="$SYS_NAME" -classpath "./bin/production/LiHoF4_transverse_field_MC/:./lib/*" simulation.montecarlo.Main -max_iter 80 -L "$L","$L" -d "$x" -extBx "$H" -continue_from_save no -max_sweeps 3 -name default -temp_schedule ./"$SYS_NAME"/temperature_schedules/temp_schedule_"$L"_"$L"_"$name"_"$H"_"$mech".txt -alpha 0.95 -verbose -mode s -seed "$seed" -Jex "$jex" "$suppress_tag" -interpolation_table_name _const
    cd scripts
    echo "finished running simulation"
    for temperature_file in ../LiHoF4/data/results/default/table_*"$seed".txt; do
        original_file=$(echo "$temperature_file" | sed 's/default/'"$name"'/g')
        text1=$(awk '{if ($0 !~ /^[[:space:]]*#/ && $1 != "index" && $1!=0) {print $0}}' "$temperature_file")
        text2=$(awk '{if ($0 !~ /^[[:space:]]*#/ && $1 != "index" && $1!=0 && $1 < 3) {print $0}}' "$original_file")
        echo "$temperature_file"
        if [ "$text1" = "$text2" ]; then
            echo "very good"
            missing_line=$(awk '{if ($0 !~ /^[[:space:]]*#/ && $1 != "index" && $1==0) {print $0}}' "$temperature_file")
            sed -i '/^[[:space:]]*index/a \'"$missing_line" "$original_file"
        else
            echo "not compatible"
            echo "$text1"
            echo "$text2"
        fi
    done
done < seeds.txt

