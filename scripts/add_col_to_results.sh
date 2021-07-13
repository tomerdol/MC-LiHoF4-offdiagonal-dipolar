#!/bin/bash
# Add new columns to existing results files
# by calling add_col.awk
# Added columns are currently: "mk2y","mk2z","tMaxConfigs" (can be changed in add_col.awk)

if [ $# -ne 1 ]
then
echo Usage $0 [directory]
exit
fi

name=$1
for file in ./"$SYS_NAME"/data/results/"$name"/table_*
do
gawk -i inplace -f ./scripts/add_col.awk $file 
done
