"""
Create a magnetic moment and energy interpolation table, i.e. calculates the
energy and magnetic moment of the "up" and "down" states of a Fe8 molecule under the
Fe8 crystal field potential and an applied magnetic field on a grid
of such applied field: (Bx,By,Bz) and saves the results so that they can be read
and used by interpolating within the grid during the MC simulation.
See further details in the documentation of the Java class simulation.montecarlo.CrystalField.
This one allows for magnetic moments that might be positive or negative
for both "up" and "down", given a large enough applied magnetic field.
It was originally introduced when using the Broyden method for the self-
consistent calculation of the magnetic moments since that required large
(possibly unphysical) fields be tested on the way to the solution.
"""
import numpy as np
from numpy import linalg as LA
import math
import sys
import os

def write_to_file(name, data, Bx, By, Bz):
	"""
	Writes the calculated table to a txt file in /Fe8/data/interactions/.
	:param name: name of the txt file
	:param data: data (3D NumPy array) to save in the file
	:param Bx: grid Bx values
	:param By: grid By values
	:param Bz: grid Bz values
	:return:
	"""
	script_dir = os.path.dirname(os.path.abspath(__file__)) #<-- absolute dir the script is in
	rel_path = "/../Fe8/data/interactions/" + name + '.txt'
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

# constants
hbar=1
D=0.294
E=0.046
u_B=0.6717 # Bohr magneton [K/T]
g_L=2 # Lande g-factor
J=10
deg_J = 2 * J + 1

# J matrices
jplus = hbar * np.diag(np.array( [ math.sqrt(J*(J+1) - m*(m+1)) for m in np.arange(-J,J) ] ), 1)
jminus = hbar * np.diag(np.array( [ math.sqrt(J*(J+1) - m*(m-1)) for m in np.arange(-J+1,J+1) ] ),-1)
jx = (jplus + jminus) * 0.5
jy = (jplus - jminus) * (-0.5j)
jz = hbar * np.diag(np.arange(-J,J+1))
I_J = np.diag(np.ones(int(round(deg_J))))

# "crystal field" Hamiltonian
H_cf = -D*LA.matrix_power(jz,2) + E*(LA.matrix_power(jx,2)-LA.matrix_power(jy,2))

# Magnetic field Zeeman term
# we use a geometric sequence so that the grid is
# denser closer to zero where we want better resolution.
meanBx = float(sys.argv[1])
meanBy = float(sys.argv[2])
maxBx=3.0       # the real max is one less than this
maxBz=1.2       # the real max is one less than this
number=40       # the real number is twice this

min_bz=0.0014
Bx = np.geomspace(1,maxBx,num=number) - 1 + np.geomspace(1,maxBx,num=number)[1]- np.geomspace(1,maxBx,num=number)[0]
Bx = np.concatenate((np.flip(-1*Bx),Bx),axis=0)
Bx += meanBx
By = np.geomspace(1,maxBx,num=number) - 1 + np.geomspace(1,maxBx,num=number)[1] - np.geomspace(1,maxBx,num=number)[0]
By = np.concatenate((np.flip(-1*By),By),axis=0)
By += meanBy
Bz = np.geomspace(1,maxBz,num=number) - 1 + min_bz
Bz = np.concatenate((np.flip(-1*Bz),Bz),axis=0)

# Create full standard table
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
			j=0
			if np.allclose(w[j],w[j:j+2],atol=10**-15):
				# for degenerate pairs of states, lift the degeneracy by diagonalizing the J_z block
				change_of_basis_matrix_inverse=np.transpose(v.T)
				change_of_basis_matrix=LA.inv(change_of_basis_matrix_inverse)
				new_jz=change_of_basis_matrix @ jz @ change_of_basis_matrix_inverse
				j=0
				magnetic_moments = LA.eigvalsh(new_jz[j:j+2,j:j+2])
				magnetic_moment_up = np.real(magnetic_moments[0])
				magnetic_moment_down = np.real(magnetic_moments[1])
				print("all close!")

			if (magnetic_moment_up < magnetic_moment_down):
				# switch moments
				temp = magnetic_moment_up
				magnetic_moment_up = magnetic_moment_down
				magnetic_moment_down = temp
				# switch energies
				temp = energy_up
				energy_up = energy_down
				energy_down = temp

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

energy_up_arr_full = np.array(res_energy_up)
energy_down_arr_full = np.array(res_energy_down)
magnetic_moment_up_arr_full = np.array(res_magnetic_moment_up)
magnetic_moment_down_arr_full = np.array(res_magnetic_moment_down)
print("check up and down magnetic moment arrays are the same: " + str(np.allclose((-1)*(np.flip(magnetic_moment_down_arr_full, 0)), magnetic_moment_up_arr_full, atol=1e-15)))
print("check up and down energy arrays are the same: " + str(np.allclose(np.flip(energy_down_arr_full, 0), energy_up_arr_full, atol=1e-15)))

energy_up_arr = energy_up_arr_full
energy_down_arr = energy_down_arr_full
magnetic_moment_up_arr = magnetic_moment_up_arr_full
magnetic_moment_down_arr = magnetic_moment_down_arr_full

# we should need only the magnetic moment and energy for one of the state and we choose to save only "up"
# to be sure, we check that the data for down exists in up by transposing Bz -> -Bz
# (and for the magnetic moment, also multiplying the results by -1)
if np.allclose((-1)*(np.flip(magnetic_moment_down_arr, 0)), magnetic_moment_up_arr, atol=1e-15):
	write_to_file('magnetic_moment_up_arr_%1.2f_0.014'%meanBx,magnetic_moment_up_arr, Bx, By, Bz)
	print('Wrote magnetic moment file: magnetic_moment_up_arr_%1.2f_0.014'%meanBx)
else:
	print('magnetic moment table not transposable!')
if np.allclose(np.flip(energy_down_arr, 0), energy_up_arr, atol=1e-15):
	write_to_file('energy_up_arr_%1.2f_0.014'%meanBx,energy_up_arr, Bx, By, Bz)
	print('Wrote energy file: energy_up_arr_%1.2f_0.014'%meanBx)
else:
	print('energy table not transposable!')


