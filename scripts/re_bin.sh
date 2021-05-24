#!/bin/bash

arrayH=( 0.0 0.3 0.6 1.0 1.5 2.0 )
arrayname=( "res_test_const2" "res_test_const3" )
arrayMech=( "true" "false" )

for H in "${arrayH[@]}"
do
for name in "${arrayname[@]}"
do
for mech in "${arrayMech[@]}"
do
    python3 bin_data.py "$H" "$name" "$mech" 4 5 6
done
done
done


