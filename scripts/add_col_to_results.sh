#!/bin/bash
name=$1
for file in ./data/results/"$name"/table_*
do
gawk -i inplace -f ./scripts/add_col.awk $file 
done
