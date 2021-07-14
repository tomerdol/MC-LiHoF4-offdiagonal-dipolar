"""
Find an initial guess for the critical temperature based on initial results (finite-size correlation length)
saved in /data/analysis. At least two different system sizes (L) must be given, as the script looks for the
value at which different system sizes cross.
Used by init_sim.sh following the initial temporary run.
"""
import numpy as np
import config
from phase_diagram_bin import find_initial_xc_from_arr, parse_arguments
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    """
    Find and print an initial guess for T_c for the given simulation with at least
    two different system sizes.
    :return: None
    """
    args = parse_arguments()
    
    all_L = args.L
    h_ex_list = args.h_ex
    mech_list = args.mech
    folderName_list = args.folder_list
    
    all_L.sort()
    
    arr_x=[]
    arr_y=[]
    
    fig, ax = plt.subplots()

    # read the correlation length data saved in /data/analysis/
    for L in all_L:
        fname = './' + config.system_name + '/data/analysis/sample_energy_'+str(L)+"_"+str(h_ex_list[0])+"_"+str(folderName_list[0])+"_"+str(mech_list[0])+".txt"
        array=np.genfromtxt(fname, skip_header=1)
        ax.plot(array[:,0],array[:,2],'o-',label='L='+str(L))
        arr_x.append(array[:,0])
        arr_y.append(array[:,2])

    # find where the different L curves approximately cross
    initial_tc=find_initial_xc_from_arr(arr_x,arr_y)
    # plot the curves and the crossing value (as a vertical red line)
    ax.set_yscale('log')
    ax.axvline(x=initial_tc, label=r'initial $T_c$', color='red')
    plt.ylabel(r'$\xi_{L} / L$')
    plt.xlabel('T')
    plt.legend(loc='best')
    fig.savefig('./' + config.system_name + '/figures/plot_%s_%s_%s.png'%(h_ex_list[0],folderName_list[0],mech_list[0]))

    # print out the result so it can be used by the calling script init_sim.sh
    print(initial_tc)

if __name__ == "__main__":
    main()