#!/bin/bash

name="temp_Fe8_test8"
arrayMech=( "false" "true" )
arrayL=( 5 6 7 8 9 )
arrayH=( 0.0 )

for L in "${arrayL[@]}"
do
for H in "${arrayH[@]}"
do
for mech in "${arrayMech[@]}"
do

echo "Binning data"
python3 bin_data.py "$SYS_NAME" "$H" "$name" "$mech" "$L"
echo "Removing raw data"
rm -f ../"$SYS_NAME"/data/results/"$name"/table_"$L"_"$L"_"$H"_*_"$mech".txt

done
done
done



