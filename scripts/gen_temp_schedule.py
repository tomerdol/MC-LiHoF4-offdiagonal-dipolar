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
	name=input('project name: ')
	bx=input('Bx: ')
	mech=input('mechanism (true/false): ')
	script_dir = os.path.dirname(os.path.abspath(__file__)) #<-- absolute dir the script is in
	rel_path = "/../temp_schedule_%s_%s_%s.txt"%(name,bx,mech)
	abs_file_path = script_dir + rel_path
	with open(abs_file_path, 'w') as outfile:
		outfile.write(','.join(map(str,xdata)))


