import numpy as np
from numpy import linalg as LA
import math
import pandas as pd
import sys
import os
import config


def write_to_file(name, data, Bx, By, Bz):
	script_dir = os.path.dirname(os.path.abspath(__file__)) #<-- absolute dir the script is in
	rel_path = '/../LiHoF4/data/interactions/' + name + '.txt'
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
maxBx=3.0	# the real max is one less than this
maxBz=3.0	# the real max is one less than this
number=40	# the real number is twice this
#min_bz=0.0193
# next:
min_bz=np.geomspace(1,maxBz,num=number)[1] - np.geomspace(1,maxBz,num=number)[0]
Bx = np.geomspace(1,maxBx,num=number) - 1 + np.geomspace(1,maxBx,num=number)[1] - np.geomspace(1,maxBx,num=number)[0]
Bx = np.concatenate((np.flip(-1*Bx),Bx),axis=0)
Bx += meanBx
By = np.geomspace(1,maxBx,num=number) - 1 + np.geomspace(1,maxBx,num=number)[1] - np.geomspace(1,maxBx,num=number)[0]
By = np.concatenate((np.flip(-1*By),By),axis=0)
Bz = np.geomspace(1,maxBz,num=number) - 1 + min_bz
Bz = np.concatenate((np.flip(-1*Bz),Bz),axis=0)

H_zeeman = [[ -g_L * u_B * (Bx[i]*jx + By[j]*jy) for i in range(len(Bx)) ] for j in range(len(By))]

# full Hamiltonian
H = [ [( H_cf + H_zeeman_xy ) for H_zeeman_xy in H_zeeman_x ] for H_zeeman_x in H_zeeman ]

H=np.array(H)
assert np.allclose(H.transpose(0,1,3,2).conj(), H) # check Hermiticity

### diagonalization
res=[]
eigen_energies=[]
for hx in H:
	res_x=[]
	eigen_energies_x=[]
	for hxy in hx:
		w,v = LA.eigh(hxy)
		res_x.append(v.T)
		eigen_energies_x.append(w)
	res.append(res_x)
	eigen_energies.append(eigen_energies_x)

eigenstates1=np.array(res)
# 1st index is H, 2nd index is eigenvalue number (from lowest to highest)
change_of_basis_matrix_inverse=np.transpose(eigenstates1,axes=(0,1,3,2))
change_of_basis_matrix=LA.inv(change_of_basis_matrix_inverse)
new_jz=change_of_basis_matrix @ jz @ change_of_basis_matrix_inverse
new_jx=change_of_basis_matrix @ jx @ change_of_basis_matrix_inverse
new_jy=change_of_basis_matrix @ jy @ change_of_basis_matrix_inverse

# print(np.matrix(new_jz))

part_new_jz=np.array(((new_jz[:,:,0,0],new_jz[:,:,0,1]),(new_jz[:,:,1,0],new_jz[:,:,1,1])))
part_new_jx=np.array(((new_jx[:,:,0,0],new_jx[:,:,0,1]),(new_jx[:,:,1,0],new_jx[:,:,1,1])))
part_new_jy=np.array(((new_jy[:,:,0,0],new_jy[:,:,0,1]),(new_jy[:,:,1,0],new_jy[:,:,1,1])))

eigenvalues=[]
eigenstates=[]

j=0
for part_new_jz_i in part_new_jz.transpose(2,3,0,1):
	eigenvalues_i=[]
	eigenstates_i=[]
	for part_new_jz_ij in part_new_jz_i:
		#     print("part_new_jz_i")
		#     print(part_new_jz_i)
		w1,v1 = LA.eigh(part_new_jz_ij)
		eigenvalues_i.append(w1)
		#     print('Bx='+str(Bx[j]))
		#     j+=1
		for i in range(v1.T[:,0].size):
			v1.T[i,:]*=np.exp(-1j*np.angle(v1.T[i,0]))
		#     print(np.angle(v1.T))
		#     print(np.abs(v1.T))
		eigenstates_i.append(v1)
	eigenstates.append(eigenstates_i)
	eigenvalues.append(eigenvalues_i)


# # checks:
# alpha = eigenstates1[0, 0, :]
# beta = eigenstates1[0, 1, :]
# up_coord=v1[:,0]
# down_coord=v1[:,1]
# print('up_coord='+str(np.angle(up_coord)))
# print('down_coord='+str(np.angle(down_coord)))
# up=up_coord[0]*alpha + up_coord[1]*beta
# down=down_coord[0]*alpha + down_coord[1]*beta
# print('up='+str(up))
# print('down='+str(down))
# print("up@down="+str(np.conj(np.transpose(up)) @ down))
# # print("up@jz@down="+str(np.conj(np.transpose(up)) @ jz @ down))
# # print("down@jz@up="+str(np.conj(np.transpose(down)) @ jz @ up))
# # print("down@jz@down="+str(np.conj(np.transpose(down)) @ jz @ down))
# # print("up@jz@up="+str(np.conj(np.transpose(up)) @ jz @ up))

eigenvalues=np.array(eigenvalues)
eigenstates=np.array(eigenstates)

upjzup=[]
downjzdown=[]
upjxup=[]
downjxdown=[]
upjyup=[]
downjydown=[]
upjzdown=[]
downjzup=[]
upjxdown=[]
downjxup=[]
upjydown=[]
downjyup=[]

for x in range(eigenstates.shape[0]):
	upjzup_x=[]
	downjzdown_x=[]
	for y in range(eigenstates.shape[1]):
		upjzup_x.append(np.conj(eigenstates[x,y,:,0].T) @ part_new_jz.transpose(2,3,0,1)[x,y] @ eigenstates[x,y,:,0])
		downjzdown_x.append(np.conj(eigenstates[x,y,:,1].T) @ part_new_jz.transpose(2,3,0,1)[x,y] @ eigenstates[x,y,:,1])
	upjzup.append(upjzup_x)
	downjzdown.append(downjzdown_x)

upjzup=np.array(upjzup)
downjzdown=np.array(downjzdown)

data=np.broadcast_to(0.5*abs(upjzup-downjzdown),(2*number,2*number,2*number))

write_to_file('magnetic_moment_up_arr_%1.2f_chak'%meanBx, data, Bx, By, Bz)

# create zeeman term for different combination of Bx,By,Bz
res_energy_up=[]
res_energy_down=[]
for i, bz in enumerate(Bz):
	res_energy_y_up=[]
	res_energy_y_down=[]
	for by in By:
		res_energy_x_up=[]
		res_energy_x_down=[]
		for bx in Bx:
			H_zeeman = u_B*g_L*(bx*jx + by*jy + bz*jz)    # zeeman term
			H = H_cf - H_zeeman                 # full hamiltonian
			w,v = LA.eigh(H)
			# initially assume lower level is up and upper is down
			energy_up = w[0]
			energy_down = w[1]
			if bz<0:
				# switch energies
				temp = energy_up
				energy_up = energy_down
				energy_down = temp
			effective_bz = 1.1*math.sqrt(bx**2 + by**2)
			H_zeeman = u_B*g_L*(bx*jx + by*jy + effective_bz*jz)    # zeeman term
			H = H_cf - H_zeeman                 # full hamiltonian
			w,v = LA.eigh(H)
			magnetic_moment_up = np.real(np.diagonal(np.conj(v.T)@jz@v)[0])
			H_zeeman = u_B*g_L*(bx*jx + by*jy + -effective_bz*jz)    # zeeman term
			H = H_cf - H_zeeman                 # full hamiltonian
			w,v = LA.eigh(H)
			magnetic_moment_down = np.real(np.diagonal(np.conj(v.T)@jz@v)[0])

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

		res_energy_y_up.append(res_energy_x_up)
		res_energy_y_down.append(res_energy_x_down)
	res_energy_up.append(res_energy_y_up)
	res_energy_down.append(res_energy_y_down)
	print(str(100*i/Bz.size) + "%")

energy_up_arr = np.array(res_energy_up)
energy_down_arr = np.array(res_energy_down)
write_to_file('energy_up_arr_%1.2f_chak'%meanBx,energy_up_arr, Bx, By, Bz)