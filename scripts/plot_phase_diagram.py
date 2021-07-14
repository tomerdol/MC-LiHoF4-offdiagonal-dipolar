import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
from deprecated import phase_diagram


def main():
    file_name_list = sys.argv[1:]
    
    all_L = file_name_list[0].split('_')[3:-1]
    
    # plot previous data
    fig, ax = plt.subplots()
    phase_diagram.plot_previous_data(ax)

    # plot w/o mechanism
    data = np.loadtxt(file_name_list[0],skiprows=1)
    if len(data.shape) > 1:
        x = data[:,0]
        y = data[:,2]
        x_err = data[:, 1]
        #x=(1.53-x[0])+x
    else:
        x = data[0]
        y=data[2]
        x_err=data[1]
    
    ax.errorbar(x,y,yerr=None, xerr=x_err, fmt='bs', label='w/o mechanism')
    
    # plot w/ mechanism
    data = np.loadtxt(file_name_list[1],skiprows=1)
    if len(data.shape) > 1:
        x = data[:,0]
        y = data[:,2]
        x_err = data[:, 1]
    #x=(1.53-x[0])+x
    else:
        x = data[0]
        y=data[2]
        x_err=data[1]
 
    ax.errorbar(x,y,yerr=None, xerr=x_err, fmt='gp', label='w/ mechanism')
    
    ax.set_xlabel(r'$T_c$')
    ax.set_ylabel(r'$B_x$')
    ax.legend()
    
    fig.savefig('../' + config.system_name + '/figures/phase_diagram_full_%s.png'%('_'.join(map(str,all_L))))
    
    plt.close()
    
if __name__ == "__main__":
    main()
