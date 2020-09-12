import numpy as np
from phase_diagram_bin import find_initial_xc_from_arr, parse_arguments
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import phase_diagram_bin.find_initial_xc_from_arr, phase_diagram_bin.parse_arguments

def main():
    args = parse_arguments()
    
    all_L = args.L
    h_ex_list = args.h_ex_list
    mech_list = args.mech
    folderName_list = args.folder_list
    
    all_L.sort()
    
    arr_x=[]
    arr_y=[]
    
    
    fig, ax = plt.subplots()
    
    for L in all_L:
        fname = "./R/sample_energy_"+str(L)+"_"+str(h_ex_list[0])+"_"+str(folderName_list[0])+"_"+str(mech_list[0])+".txt"
        array=np.genfromtxt(fname, skip_header=1)
        ax.plot(array[:,0],array[:,2],'o-',label='L='+str(L))
        arr_x.append(array[:,0])
        arr_y.append(array[:,2])
    
    initial_tc=find_initial_xc_from_arr(arr_x,arr_y)
    ax.set_yscale('log')
    ax.axvline(x=initial_tc, label=r'initial $T_c$', color='red')
    plt.ylabel(r'$\xi_{L} / L$')
    plt.xlabel('T')
    plt.legend(loc='best')
    fig.savefig('./R/graphs/plot_%s_%s_%s.png'%(h_ex_list[0],folderName_list[0],mech_list[0]))
    
    print(initial_tc)

if __name__ == "__main__":
    main()