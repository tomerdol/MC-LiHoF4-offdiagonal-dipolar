import numpy as np
import sys
import os

start = sys.argv[1]
end = sys.argv[2]
num = sys.argv[3]

xdata = np.geomspace(float(start),float(end),num=int(num))

print(','.join(map(str,xdata)))
ans=input("Overwrite current temperature_schedule file [y/n]? ")
if (ans=='y'):
    script_dir = os.path.dirname(os.path.abspath(__file__)) #<-- absolute dir the script is in
    name=input('project name: ')
    bx=input('Bx: ')
    mech=input('mechanism (true/false): ')
    num_of_Ls=int(input('Number of system sizes: '))
    for i in range(0,num_of_Ls):
        L=input('L: ')
        rel_path = "/../temperature_schedules/temp_schedule_%s_%s_%s_%s_%s.txt"%(L,L,name,bx,mech)
        abs_file_path = script_dir + rel_path
        with open(abs_file_path, 'w') as outfile:
            outfile.write(','.join(map(str,xdata)))


