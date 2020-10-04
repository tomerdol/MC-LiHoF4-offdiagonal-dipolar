import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
import math
import csv
import os
import analysis_tools, bin_data
from itertools import cycle

def bootstrap_resample(X, n=None):
    """ Bootstrap resample an array_like
    Parameters
    ----------
    X : array_like
      data to resample
    n : int, optional
      length of resampled array, equal to len(X) if n==None
    Results
    -------
    returns X_resamples
    """
    if n == None:
        n = len(X)
        
    resample_i = np.floor(np.random.rand(n)*len(X)).astype(int)
    X_resample = X[resample_i]
    return X_resample

def plot_multiple(all_L, xdata, all_ydata, err_y, h_ex):
    markers=['o','s','^','D','v']
    toplot=[]
    fig = plt.figure()
    
    for i in range(0,len(all_L)):
        data_label = 'L=%d' %all_L[i]
        plt.errorbar(xdata,all_ydata[i],yerr=err_y[i], fmt=markers[i % len(markers)]+'-', capsize=3, label=data_label)
    ax = plt.gca()
    ax.set_yscale("log")
    plt.title('Correlation length vs. T')
    plt.xlabel('T')
    plt.ylabel(r'$\xi_{L} / L$')
    plt.legend(loc='best')
    axes = plt.gca()
    #axes.set_ylim([0.5,1])

    fig.savefig('../figures/plot_%s.png'%h_ex)

def get_binder(m2,m4):
    m2=np.mean(m2)
    m4=np.mean(m4)
    return 0.5*(3-(m4/(m2**2)))

def get_correlation_length(mk2,m2,L):
    m2=np.mean(m2)
    mk2=np.mean(mk2)
    try:
        correlation_length=math.sqrt((m2/mk2)-1)/(2*L*math.sin(math.pi/L))
    except ValueError:
        correlation_length=0
    return correlation_length

def calc_error_correction(arr, iter):
    binned_arr=np.copy(arr)
    error_est_arr = np.zeros(iter)
    for i in range(iter):
        error_est_arr[i]=np.std(binned_arr)
        binned_arr=get_binned_array(binned_arr)
    return error_est_arr
    
def get_binned_array(arr):
    if arr.size%2==0:
        return 0.5*(arr[::2]+arr[1::2])
    else:
        return 0.5*(arr[:-1:2]+arr[1::2])
    

def main_plot(simulations, boot_num, plot_options, to_plot=''):
    Nsigma=1.
    markers=['o','s','^','D','v']
    

    all_y_curves = []
    
    for i, sim in enumerate(simulations.itertuples()):
        y_to_plot=[]
        y_err_to_plot=[]
        
        
        y = bin_data.read_binned_data(sim, use_latest=False, use_bin=sim.eq_bin)

        single_ydata=[]
        for boot_index in range(boot_num):
            if to_plot=='':
                # this means a scaling function (Binder ratio or correlation length should be plotted, according to what is defined in plot_options)
                single_ydata.append(plot_options['func'](y['Magnetization^2'].sample(frac=1,replace=True),y['Magnetization^4'].sample(frac=1,replace=True),y['mk2'+plot_options['corr_length_axis']].sample(frac=1,replace=True),sim.L))  # current y value
            else:
                # this means 'to_plot' should be plotted
                single_ydata.append(plot_options['func'](y[to_plot].sample(frac=1,replace=True)))
            #print(boot_index)
        
        all_yboot=np.array(single_ydata)
        
        y_to_plot = pd.Series([np.mean(all_yboot,0), Nsigma * np.std(all_yboot,0)],index=['y_to_plot','y_to_plot_err'])
        #y_err_to_plot = Nsigma * np.std(all_yboot,0)
        
        y_to_plot=y_to_plot.append(pd.Series(sim._asdict()))
        #y_err_to_plot=y_err_to_plot.append(pd.Series(sim._asdict()))
        
        #print('Curve to be plotted (L=%s): '%sim.L)
        #print(y_to_plot)
        all_y_curves.append(y_to_plot)
        #all_y_curves_err.append(y_err_to_plot)
       
    all_y_curves = pd.DataFrame(all_y_curves)
    
    
    fig, ax = plt.subplots(figsize=(8,6))
    
    for (label, df), marker in zip(all_y_curves.groupby(['Bex','L','folderName','mech']), cycle(markers)):
        df.plot(x='T',y='y_to_plot', yerr='y_to_plot_err', ax=ax, label=str(label), capsize=3, marker=marker)
    
    ax.set_xlabel('T')
    
    ax.set_yscale(plot_options['axis_yscale'])
    plt.ylabel(plot_options['Name'])

    plt.legend(loc='best')
    plt.tight_layout()
    axes = plt.gca()

    fig.savefig('../figures/plot_%s_%s_%s_%s_%s.png'%(to_plot,'_'.join(map(str,simulations['Bex'].unique().tolist())),'_'.join(map(str,simulations['mech'].unique().tolist())),'_'.join(map(str,simulations['L'].unique().tolist())),'_'.join(map(str,simulations['folderName'].unique().tolist()))), dpi=300)
    fig.savefig('../figures/plot_%s_%s_%s_%s_%s.eps'%(to_plot,'_'.join(map(str,simulations['Bex'].unique().tolist())),'_'.join(map(str,simulations['mech'].unique().tolist())),'_'.join(map(str,simulations['L'].unique().tolist())),'_'.join(map(str,simulations['folderName'].unique().tolist()))))
    plt.close()
    #return ([xdata[l] for l in sorted(xdata.keys())],all_y_curves)
    return all_y_curves


def parse_arguments():  
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    
    parser = ArgumentParser(description="Analyzes Monte Carlo results and plots correlation length curves for LiHoF4", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument( "-L", nargs='+', type=int, required=True, help = "Linear system sizes. At least 2 required.")
    parser.add_argument( "-b", "--boot_num", type=int, default = 100, help = "Number of bootstrap samples.")
    parser.add_argument( "--h_ex", type=float, help = "External magnetic field value, Bex." , required=True)
    parser.add_argument( "-m", "--mech", nargs='+', choices=['true','false'], help = ("Whether internal fields are suppressed or not. \'false\' means "
    "that they aren't so the mechanism is on, and \'true\' means that they are and the mechanism is off." ), required=True)
    parser.add_argument( "-f", "--folder_list", nargs='+', type=str, help = "List of folders in \'data/results/\' in which results should be found. " , required=True)
    parser.add_argument( "--to_plot", type=str, nargs='?', default='corr_length', help = "Which observable should be plotted. Default is Correlation length / L.")
    
    args = parser.parse_args()
    
    return args



def main():
    args = parse_arguments()
    L = args.L                      
    boot_num = args.boot_num           
    h_ex = args.h_ex
    mech = args.mech
    folderName = args.folder_list 
    to_plot = args.to_plot
    
    simulations = analysis_tools.get_simulations(L, folderName, h_ex, mech)
    
    from fit6 import get_binder, get_correlation_length
    #plot_options = {'Name':'g', 'axis_yscale':'linear', 'func':get_binder}
    plot_options = {'Name':r'$\xi_{L} / L$', 'axis_yscale':'log', 'func':get_correlation_length}
    
    main_plot(simulations, boot_num, plot_options)
    #os.system("rsync -avzhe ssh ../figures/ tomerdol@newphysnet1:~/graphs/")

if __name__ == "__main__":
    main()
