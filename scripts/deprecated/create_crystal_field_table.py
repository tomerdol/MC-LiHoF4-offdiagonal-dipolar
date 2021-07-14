"""
Create a magnetic moment and energy table for the Ho ion in LiHoF_4 including hyperfine interaction.
The 8 hf levels that constitute either the "up" or "down" states are assumed to be in thermal equilibrium
with the environment so the magnetic moment is taken as a Boltzmann average of these.
This is not necessarily a valid model and did not produce sensible results so it is deprecated for the time being.
"""
import numpy as np
from numpy import linalg as LA
import math
import sys

# constants
hbar=1
J=8
I=7/2
deg_J = 2 * J + 1
deg_I = 2 * I + 1
u_B=0.6717
g_L=5/4
A = 0.039 # hyperfine interaction


# J matrices
jplus = hbar * np.diag(np.array( [ math.sqrt(J*(J+1) - m*(m+1)) for m in np.arange(-J,J) ] ), 1)
jminus = hbar * np.diag(np.array( [ math.sqrt(J*(J+1) - m*(m-1)) for m in np.arange(-J+1,J+1) ] ),-1)
jx = (jplus + jminus) * 0.5
jy = (jplus - jminus) * (-0.5j)
jz = hbar * np.diag(np.arange(-J,J+1))
I_J = np.diag(np.ones(int(round(deg_J))))

# I matrices
iplus = hbar * np.diag(np.array( [ math.sqrt(I*(I+1) - m*(m+1)) for m in np.arange(-I,I) ] ), 1)
iminus = hbar * np.diag(np.array( [ math.sqrt(I*(I+1) - m*(m-1)) for m in np.arange(-I+1,I+1) ] ),-1)
ix = (iplus + iminus) * 0.5
iy = (iplus - iminus) * (-0.5j)
iz = hbar * np.diag(np.arange(-I,I+1))
I_I = np.diag(np.ones(int(round(deg_I))))

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
B06 = 4.64e-6
B44C = 0.0418
B46C = 8.12e-4
B46S = 1.137e-4

# crystal field Hamiltonian
H_cf = B02*O02 + B04*O04 + B06*O06 + B44C*O44C + B46C*O46C + B46S*O46S

# linear spacing for B:
meanBx=float(sys.argv[1])
Bx = meanBx + np.linspace(-1.,1.,num=21)
By = np.linspace(-1.,1.,num=21)
Bz = np.linspace(-1.,1.,num=1601)
# extend the ranges by -+1.0. the linear interpolation will give an intermediate
# value between that for Bx=-1.0 and Bx=-2.0 which is acceptable.
Bx = np.concatenate(([Bx[0]-1.],Bx,[Bx[-1]+1.]))
By = np.concatenate(([By[0]-1.],By,[By[-1]+1.]))
Bz = np.concatenate(([Bz[0]-1.],Bz,[Bz[-1]+1.]))

# with hf
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
            H = np.kron(H_cf - H_zeeman,I_I) + A * (np.kron(jx,ix) + np.kron(jy,iy) + np.kron(jz,iz)) # full hamiltonian
            w,v = LA.eigh(H)

            # get a new Jz matrix with respect to the new eigenenergy basis
            change_of_basis_matrix_inverse=np.transpose(v.T)
            change_of_basis_matrix=LA.inv(change_of_basis_matrix_inverse)
            new_jz=change_of_basis_matrix @ np.kron(jz,I_I) @ change_of_basis_matrix_inverse
            # for degenerate pairs of states, lift the degeneracy by diagonalizing the J_z block
            for j in range(16):
                if np.allclose(w[j],w[j:j+2],atol=10**-10):
                    new_jz[j:j+2,j:j+2]=np.diag(LA.eigvalsh(new_jz[j:j+2,j:j+2]))

            # ********* UP *********
            # with hf
            k=8
            arr = np.real(np.diag(new_jz)[:16])
            # get the 8 higher magnetic moments
            magnetic_moments_up = arr[np.sort(np.argpartition(arr, len(arr)-k)[-k:])]
            # Boltzmann weighted average
            energies_up = w[:16][np.sort(np.argpartition(arr, len(arr)-k)[-k:])]
            beta=1/1.53 # T=1.53 K
            magnetic_moment_up = np.average(magnetic_moments_up,weights=np.exp(-beta*energies_up))
            energy_up = np.average(energies_up,weights=np.exp(-beta*energies_up))
            # ******** DOWN ********
            magnetic_moments_down = arr[np.sort(np.argpartition(arr, len(arr)-k)[:k])]
            # Boltzmann weighted average
            energies_down = w[:16][np.sort(np.argpartition(arr, len(arr)-k)[:k])]
            beta=1/1.53
            magnetic_moment_down = np.average(magnetic_moments_down,weights=np.exp(-beta*energies_down))
            energy_down = np.average(energies_down,weights=np.exp(-beta*energies_down))

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
    sys.stdout.flush()

energy_up_arr = np.array(res_energy_up)
energy_down_arr = np.array(res_energy_down)
magnetic_moment_up_arr = np.array(res_magnetic_moment_up)
magnetic_moment_down_arr = np.array(res_magnetic_moment_down)
print("check up and down magnetic moment arrays are the same: " + str(np.allclose((-1)*(np.flip(magnetic_moment_down_arr, 0)), magnetic_moment_up_arr, atol=1e-15)))
print("check up and down energy arrays are the same: " + str(np.allclose(np.flip(energy_down_arr, 0), energy_up_arr, atol=1e-15)))

np.savez('cf_table_%1.2f.npy'%meanBx,energy_up_arr,energy_down_arr,magnetic_moment_up_arr,magnetic_moment_down_arr)
