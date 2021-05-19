#!/bin/bash

name="res_test"
arrayMech=( "false" "true" )
arrayL=( 6 5 4 )
arrayH=( 0.0 0.3 0.6 )

for L in "${arrayL[@]}"
do
for H in "${arrayH[@]}"
do
for mech in "${arrayMech[@]}"
do

echo "Binning data"
python3 bin_data.py "$H" "$name" "$mech" "$L"
echo "Removing raw data"
rm -f ../data/results/"$name"/table_"$L"_"$L"_"$H"_*_"$mech".txt

done
done
done



