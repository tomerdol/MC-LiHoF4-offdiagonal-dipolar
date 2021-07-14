# run the Java program to create and save the table od dipolar interactions between spins, using the Ewald summation method
# run from the root project directory

# Check input
if [ $# -ne 3 ]
then
echo Usage $0 [L] [k_cutoff] [real_cutoff]
exit
fi

# run the java program which performs the calculations and saves the data to ./${SYS_NAME}/data/interaction/intTable_Lx_Lz.txt
java -Dsystem=${SYS_NAME} -classpath "./bin/production/LiHoF4_transverse_field_MC/:./lib/*" utilities.ewaldSum "$1" "$1" "$2" "$3"
