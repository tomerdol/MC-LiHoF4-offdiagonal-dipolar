import numpy as np
from numpy import linalg as LA
import math
import pandas as pd
import sys
import os

def write_to_file(name, data, Bx, By, Bz):
	script_dir = os.path.dirname(os.path.abspath(__file__)) #<-- absolute dir the script is in
	rel_path = "/../data/interactions/" + name + '.txt'
	print(script_dir)
	abs_file_path = script_dir + rel_path
	with open(abs_file_path, 'w') as outfile:
		# I'm writing a header here just for the sake of readability
		# Any line starting with "#" will be ignored by numpy.loadtxt
		outfile.write('# {0}\n'.format(data.shape))
		outfile.write('# Bz: %s\n'%np.array2string(Bz,threshold=sys.maxsize,max_line_width=np.nan,separator=','))
		outfile.write('# By: %s\n'%np.array2string(By,threshold=sys.maxsize,max_line_width=np.nan,separator=','))
		outfile.write('# Bx: %s\n'%np.array2string(Bx,threshold=sys.maxsize,max_line_width=np.nan,separator=','))
		# Iterating through a ndimensional array produces slices along
		# the last axis. This is equivalent to data[i,:,:] in this case
		for data_slice in data:

			# The formatting string indicates that I'm writing out
			# the values in left-justified columns 7 characters in width
			# with 2 decimal places.  
			np.savetxt(outfile, data_slice, fmt='%-7.8f')

			# Writing out a break to indicate different slices...
			outfile.write('# New Bz slice\n')

#constants
hbar = 1
J=8
deg_J = 2 * J + 1
g_L = 5/4 # Lande g-factor
u_B = 0.6717 # Bohr magneton [K/T]


# J matrices
jplus = hbar * np.diag(np.array( [ math.sqrt(J*(J+1) - m*(m+1)) for m in np.arange(-J,J) ] ), 1)
jminus = hbar * np.diag(np.array( [ math.sqrt(J*(J+1) - m*(m-1)) for m in np.arange(-J+1,J+1) ] ),-1)
jx = (jplus + jminus) * 0.5
jy = (jplus - jminus) * (-0.5j)
jz = hbar * np.diag(np.arange(-J,J+1))
I_J = np.diag(np.ones(int(round(deg_J))))


# crystal field equivalent operators
O02 = 3 * LA.matrix_power(jz,2) - J*(J+1)*I_J
O04 = 35 * LA.matrix_power(jz,4) - 30 * J * (J+1) * LA.matrix_power(jz,2) + 25*LA.matrix_power(jz,2) - 6 * J * (J+1) * I_J + 3 * J**2 * (J+1)**2 * I_J
O44C = 0.5 * (LA.matrix_power(jplus,4) + LA.matrix_power(jminus,4))
O06 = 231 * LA.matrix_power(jz,6) - 315*J*(J+1)*LA.matrix_power(jz,4) + 735*LA.matrix_power(jz,4) + 105 * J**2 * (J+1)**2 * LA.matrix_power(jz,2) - 525*J*(J+1)*LA.matrix_power(jz,2) + 294*LA.matrix_power(jz,2) - 5 * J**3 * (J+1)**3 * I_J + 40 * J**2 *(J+1)**2 * I_J - 60*J*(J+1)*I_J
O46C1 = 0.25 * (LA.matrix_power(jplus,4) + LA.matrix_power(jminus,4)) @ (11*LA.matrix_power(jz,2) - J*(J+1)*I_J - 38*I_J)
O46S1 = -0.25j * (LA.matrix_power(jplus,4) - LA.matrix_power(jminus,4)) @ (11 * LA.matrix_power(jz,2) - J*(J+1)*I_J - 38*I_J)
O46C = O46C1 + np.transpose(np.conj(O46C1))
O46S = O46S1 + np.transpose(np.conj(O46S1))


# crystal field parameters
B02 = -0.696
B04 = 4.06e-3
B06 =  4.64e-6
B44C = 0.0418
B46C = 8.12e-4
B46S = 1.137e-4


# crystal field Hamiltonian
H_cf = B02*O02 + B04*O04 + B06*O06 + B44C*O44C + B46C*O46C + B46S*O46S


# Magnetic field Zeeman term
meanBx = float(sys.argv[1])
"""
Bx = np.linspace(meanBx-2.,meanBx+2.,num=50)
By = np.linspace(-2.,2.,num=50)
Bz = np.linspace(-2.,2.,num=51)
"""
maxBx=3.0	# the real max is one less than this
maxBz=3.0	# the real max is one less than this
number=40	# the real number is twice this
#min_bz=0.0193
# next:
min_bz=0.058
Bx = np.geomspace(1,maxBx,num=number) - 1 + np.geomspace(1,maxBx,num=number)[1] - np.geomspace(1,maxBx,num=number)[0]
Bx = np.concatenate((np.flip(-1*Bx),Bx),axis=0)
Bx += meanBx
By = np.geomspace(1,maxBx,num=number) - 1 + np.geomspace(1,maxBx,num=number)[1] - np.geomspace(1,maxBx,num=number)[0]
By = np.concatenate((np.flip(-1*By),By),axis=0)
Bz = np.geomspace(1,maxBz,num=number) - 1 + min_bz
Bz = np.concatenate((np.flip(-1*Bz),Bz),axis=0)
"""
Bx = np.geomspace(1,maxBx,num=number) - 1 + np.geomspace(1,maxBx,num=number)[1] - np.geomspace(1,maxBx,num=number)[0]
Bx = np.concatenate((np.flip(-1*Bx),Bx),axis=0)
Bx += meanBx
By = np.geomspace(1,maxBx,num=number) - 1 + np.geomspace(1,maxBx,num=number)[1] - np.geomspace(1,maxBx,num=number)[0]
By = np.concatenate((np.flip(-1*By),By),axis=0)
Bz = np.geomspace(1,maxBz,num=number) - 1 + np.geomspace(1,maxBz,num=number)[1] - np.geomspace(1,maxBz,num=number)[0]
Bz = np.concatenate((np.flip(-1*Bz),Bz),axis=0)
"""
# create zeeman term for different combination of Bx,By,Bz
res_energy_up=[]
res_energy_down=[]
res_magnetic_moment_up=[]
res_magnetic_moment_down=[]
for i, bz in enumerate(Bz):
	res_energy_y_up=[]
	res_energy_y_down=[]
	res_magnetic_moment_y_up=[]
	res_magnetic_moment_y_down=[]
	for by in By:
		res_energy_x_up=[]
		res_energy_x_down=[]
		res_magnetic_moment_x_up=[]
		res_magnetic_moment_x_down=[]
		for bx in Bx:
			H_zeeman = u_B*g_L*(bx*jx + by*jy + bz*jz)    # zeeman term
			H = H_cf - H_zeeman                 # full hamiltonian
			w,v = LA.eigh(H)
			# initially assume lower level is up and upper is down
			energy_up = w[0]
			energy_down = w[1]
			magnetic_moment_up = np.real(np.diagonal(np.conj(v.T)@jz@v)[0])
			magnetic_moment_down = np.real(np.diagonal(np.conj(v.T)@jz@v)[1])
			
			
			
			if (magnetic_moment_up < magnetic_moment_down):
				# switch moments
				temp = magnetic_moment_up
				magnetic_moment_up = magnetic_moment_down
				magnetic_moment_down = temp
				# switch energies
				temp = energy_up
				energy_up = energy_down
				energy_down = temp
			
			if  bz==0 and (bx!=0 or by!=0): #math.isclose(bz,0,abs_tol=1e-15)
				print('(%s,%s,%s): %s'%(bx,by,bz,magnetic_moment_up))
				print('(%s,%s,%s): %s'%(bx,by,bz,magnetic_moment_down))
				print('(%s,%s,%s): %s'%(bx,by,bz,energy_up))
				print('(%s,%s,%s): %s'%(bx,by,bz,energy_down))
			res_energy_x_up.append(energy_up)
			res_energy_x_down.append(energy_down)
			
			res_magnetic_moment_x_up.append(magnetic_moment_up)
			res_magnetic_moment_x_down.append(magnetic_moment_down)
		res_energy_y_up.append(res_energy_x_up)
		res_energy_y_down.append(res_energy_x_down)
		res_magnetic_moment_y_up.append(res_magnetic_moment_x_up)
		res_magnetic_moment_y_down.append(res_magnetic_moment_x_down)
	res_energy_up.append(res_energy_y_up)
	res_energy_down.append(res_energy_y_down)
	res_magnetic_moment_up.append(res_magnetic_moment_y_up)
	res_magnetic_moment_down.append(res_magnetic_moment_y_down)
	print(str(100*i/Bz.size) + "%")
		
energy_up_arr = np.array(res_energy_up)
energy_down_arr = np.array(res_energy_down)
magnetic_moment_up_arr = np.array(res_magnetic_moment_up)
magnetic_moment_down_arr = np.array(res_magnetic_moment_down)
print("check up and down magnetic moment arrays are the same: " + str(np.allclose((-1)*(np.flip(magnetic_moment_down_arr, 0)), magnetic_moment_up_arr, atol=1e-15)))
print("check up and down energy arrays are the same: " + str(np.allclose(np.flip(energy_down_arr, 0), energy_up_arr, atol=1e-15)))
if (np.allclose((-1)*(np.flip(magnetic_moment_down_arr, 0)), magnetic_moment_up_arr, atol=1e-15)):
	write_to_file('magnetic_moment_up_arr_%1.2f_0.058'%meanBx,magnetic_moment_up_arr, Bx, By, Bz)
 	#pass
else:
    print('magnetic moment table not transposable!')
if (np.allclose(np.flip(energy_down_arr, 0), energy_up_arr, atol=1e-15)):
	write_to_file('energy_up_arr_%1.2f_0.058'%meanBx,energy_up_arr, Bx, By, Bz)
 	#pass
else:
    print('energy table not transposable!')


