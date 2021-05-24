#!/bin/bash
name=$1
for file in ./"$SYS_NAME"/data/results/"$name"/table_*
do
gawk -i inplace -f ./scripts/add_col.awk $file 
done
