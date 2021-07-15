#!/bin/bash
# plot histograms of multiple observables for multiple projects,
# close to the found Tc of each project.
# uncomment the desired project block -- colors refer to the marker
# colors in the phase diagram.

# tcs array contains the found Tc values, around which we want to plot the histograms
declare -A tcs
: '
#green
arrayH=( 0.3 0.6 1.0 1.5 )
name="res_test_ex_0.014_re"
mech="true"
tcs[0.3]=1.77
tcs[0.6]=1.745
tcs[1.0]=1.72
tcs[1.5]=1.64

#orange
arrayH=( 0.0 0.3 0.6 1.0 1.5 )
name="res_test_ex_0.014"
mech="false"
tcs[0.0]=1.61
tcs[0.3]=1.62
tcs[0.6]=1.651
tcs[1.0]=1.712
tcs[1.5]=1.716
'
#red
arrayH=( 1.5 )
name="res_test_ex_0.014"
mech="true"
tcs[0.0]=1.51
tcs[0.3]=1.507
tcs[0.6]=1.48
tcs[1.0]=1.45
tcs[1.5]=1.26

: '
#blue
#arrayH=( 0.0 0.3 0.6 1.0 1.5 )
#name="res_test_const2"
mech="false"
tcs[0.0]=1.556
tcs[0.3]=1.548
tcs[0.6]=1.534
tcs[1.0]=1.505
tcs[1.5]=1.452
tcs[2.0]=1.401
arrayH=( 2.0 )
name="res_test_const3"
'
for H in "${arrayH[@]}"
do
    python3 plot_histogram.py -L 6 --h_ex "$H" -m "$mech" -T 5.0 ${tcs[$H]} 0.0 -f "$name" --to_plot localBx
    python3 plot_histogram.py -L 6 --h_ex "$H" -m "$mech" -T 5.0 ${tcs[$H]} 0.0 -f "$name" --to_plot localBy
    python3 plot_histogram.py -L 6 --h_ex "$H" -m "$mech" -T 5.0 ${tcs[$H]} 0.0 -f "$name" --to_plot localBz
    python3 plot_histogram.py -L 6 --h_ex "$H" -m "$mech" -T 5.0 ${tcs[$H]} 0.0 -f "$name" --to_plot spinSize
done

