#!/bin/bash

# get list of L's received as cmdline args to find the result file names at the end
fname=""
for l in "$@"
do
fname="$fname""$l"_
done


# wait until fastup is low
thresh=`/storage/scripts/fastup | grep load | awk '{print $5}' | awk -F '.' '{print $1}'`
while [ $thresh -ge 4 ]
do
        echo "Fastup load average is high >= 4 !!!"
        sleep 30
        thresh=`/storage/scripts/fastup | grep load | awk '{print $5}' | awk -F '.' '{print $1}'`
done

# create phase diagram w/o mech
echo "~~~~~~~~~~~~~~~~ without mechanism ~~~~~~~~~~~~~~~~ "
python3 phase_diagram.py --boot_num 500 --h_ex_list 0.0 --mech true --folder_list short_range -L "$@"

# wait again for low gastup
thresh=`/storage/scripts/fastup | grep load | awk '{print $5}' | awk -F '.' '{print $1}'`
while [ $thresh -ge 4 ]
do
        echo "Fastup load average is high >= 4 !!!"
        sleep 30
        thresh=`/storage/scripts/fastup | grep load | awk '{print $5}' | awk -F '.' '{print $1}'`
done

# create phase diagram w/ mech
echo "~~~~~~~~~~~~~~~~~ with mechanism ~~~~~~~~~~~~~~~~~~ "
python3 phase_diagram.py --boot_num 500 --h_ex_list 0.0 --mech false --folder_list short_range -L "$@"

# wait again for low gastup
thresh=`/storage/scripts/fastup | grep load | awk '{print $5}' | awk -F '.' '{print $1}'`
while [ $thresh -ge 4 ]
do
        echo "Fastup load average is high >= 4 !!!"
        sleep 30
        thresh=`/storage/scripts/fastup | grep load | awk '{print $5}' | awk -F '.' '{print $1}'`
done

# plot the two curves
python3 plot_phase_diagram.py phase_diagram_true_"$fname"res.txt phase_diagram_false_"$fname"res.txt

#transfer the figures
rsync -avzhe ssh ../figures/ tomerdol@newphysnet1:~/graphs/
