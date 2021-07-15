"""
Plot various observables vs. T.
Uses results from a single long MC run.
Deprecated since independent runs are now used -- use plot_bin.py instead.
"""
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
import math
import csv
import os
import config

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

    fig.savefig('../'+config.system_name+'/figures/plot_%s.png'%h_ex)

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
    

def main_plot(all_L, L_equilibrated_min_value, tau_dict, boot_num, h_ex, mech, folderName, xdata, folder='../'+config.system_name+'/data/results', start_in_middle=True):
    all_yboot=[]
    Nsigma=1.
    markers=['o','s','^','D','v']
    toplot=[]
    fig = plt.figure()

    all_y_curves=[]
    for i, L in enumerate(all_L):
        y_to_plot=[]
        y_err_to_plot=[]
        data_label = 'L=%d' %L
        
        for T in xdata[L]:
            fname=folder+'/'+folderName[L]+'/table_'+str(L)+'_'+str(L)+'_'+str(h_ex)+'_'+str(T)+'_'+mech+'.txt'
            
            try:
                col_names = pd.read_csv(fname, delim_whitespace=True, error_bad_lines=False, comment='#', nrows=0).columns
                types_dict = {'index': int, 'swap': int}
                types_dict.update({col: np.float64 for col in col_names if col not in types_dict})
                y = pd.read_csv(fname, delim_whitespace=True, error_bad_lines=False, index_col='index', comment='#',dtype=types_dict)
                #y = pd.read_csv(fname, delim_whitespace=True, error_bad_lines=False, comment='#')
                #L_equilibrated_min_value[L]=int(y.shape[0]/2)
                if (not L in L_equilibrated_min_value) and start_in_middle==True:
                    L_equilibrated_min_value[L]=int(len(y.index)*0.5)
                mk2 = y['mk2x'].loc[L_equilibrated_min_value[L]::tau_dict[L]]
                energy = y['Energy'].loc[L_equilibrated_min_value[L]::tau_dict[L]]
                m = np.abs(y['Magnetization'].loc[L_equilibrated_min_value[L]::tau_dict[L]] )
                if 'Magnetization^2' in y.columns:
                    m2 = y['Magnetization^2'].loc[L_equilibrated_min_value[L]::tau_dict[L]]
                else: 
                    m2 = y['Magnetization'].loc[L_equilibrated_min_value[L]::tau_dict[L]]**2
            except Exception as e:
                print(fname + ' is empty')
                print(e)
                y=[]
                m2=[]
                mk2=[]
                m=[]
                energy=[]
            else:
                print('read %s : %s'%(fname,m2.size))
                pass
            single_ydata=[]
            for boot_index in range(boot_num):
                try:
                    single_ydata.append(get_correlation_length(mk2.sample(frac=1,replace=True),m2.sample(frac=1,replace=True),L))  # current y value
                    #single_ydata.append(get_binder(m2.sample(frac=1,replace=True),(m2**2).sample(frac=1,replace=True)))  # current y value
                    #single_ydata.append(np.mean(energy.sample(frac=1,replace=True)))  # current y value

                    #print('read ' + fname)
                except Exception as e:
                    print(fname + " not yet equilibrated")
                    raise e
                #print(boot_index)
            all_yboot=np.array(single_ydata)
            y_to_plot.append(np.mean(all_yboot,0))
            y_err_to_plot.append(Nsigma * np.std(all_yboot,0))
            
            del y
            del mk2
            del m2
            del m
            del energy
        print('Curve to be plotted (L=%s): '%L)
        print(y_to_plot)
        all_y_curves.append(y_to_plot)
        plt.errorbar(xdata[L],y_to_plot,yerr=y_err_to_plot, fmt=markers[i % len(markers)]+'-', capsize=3, label=data_label)
    ax = plt.gca()
    ax.set_yscale("log")
    #plt.title('Correlation length vs. T')
    #plt.title('Binder ratio vs. T')

    plt.xlabel('T')
    plt.ylabel(r'$\xi_{L} / L$')
    #plt.ylabel('g')
    #plt.ylabel('E')

    plt.legend(loc='best')
    plt.tight_layout()
    axes = plt.gca()
    fig.savefig('../'+config.system_name+'/figures/plot_%s_%s_%s_cl.png'%(h_ex,mech,'_'.join(map(str,all_L))))
    plt.close()
    #return ([xdata[l] for l in sorted(xdata.keys())],all_y_curves)
    return all_y_curves


def parse_arguments():  
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description="Analyzes Monte Carlo results and plots correlation length curves.", formatter_class=ArgumentDefaultsHelpFormatter, parents=[config.parse_arguments()], conflict_handler='resolve')
    parser.add_argument( "-L", nargs='+', type=int, required=True, help = "Linear system sizes. At least 2 required.")
    # parser.add_argument( "--h_ex", type=float, help = "External magnetic field value, Bex." , required=True)
    # parser.add_argument( "-m", "--mech", choices=['true','false'], help = ("Whether internal fields are suppressed or not. \'false\' means "
    # "that they aren't so the mechanism is on, and \'true\' means that they are and the mechanism is off." ), required=True)
    # parser.add_argument( "-f", "--folder_list", nargs='+', type=str, help = "List of folders in \'data/results\' in which results should be found. " , required=True)
    #parser.add_argument( "-o", "--overwrite_tmp", action='store_true', default=False, help = ("Overwrite parsed files in /tmp. If not given, existing files will "
    #"be used (which probably means older results)."))
    parser.add_argument( "--tmp", default=False, action="store_true", help = "Read from /tmp")
    args = parser.parse_args()
    if len(args.folder_list)>1 and len(args.folder_list)!=len(args.L): 
        parser.error("--folder_list and -L argument number mismatch.")

    config.system_name = args.system_name
    return args



def main():
    args = parse_arguments()
    #all_L = [5,6]
    all_L = args.L                      
    #L_equilibrated_min_value = {2:128,3:262144,4:262144, 5:94500, 6:12750, 7:2775, 8:1170} 
    if args.tmp:
        L_equilibrated_min_value = {l:0 for l in all_L}
        folder='/tmp'
    else:
        L_equilibrated_min_value = {}
        folder='../'+config.system_name+'/data/results'
    boot_num = args.boot_num           
    h_ex = args.h_ex         
    tau_dict={k:2 for k in all_L}
    mech = args.mech                   
    #boot_num = int(sys.argv[1])
    folderName_list = args.folder_list 
    #h_ex = float(sys.argv[2])
    #mech = sys.argv[3]
    #folderName = sys.argv[4]
   
    # create L-folder dict
    if len(folderName_list)==1:
        folderName_dict={k:folderName_list[0] for k in all_L}
    else: 
        folderName_dict={k:folderName_list[i] for i,k in enumerate(all_L)}
    xdata={}
    for L in all_L:
        with open('../temp_schedule_' + folderName_dict[L] + '_' + str(h_ex) + '_' + mech + '.txt','r') as temperature_schedule_file:
            reader = csv.reader(temperature_schedule_file)
            temp_list=list(reader) 
        xdata[L]=temp_list[0]
        xdata[L]=[float(i) for i in xdata[L]]
    
    main_plot(all_L, L_equilibrated_min_value, tau_dict, boot_num, h_ex, mech, folderName_dict, xdata, folder=folder)
    #os.system("rsync -avzhe ssh ../"+config.system_name+"/figures/ tomerdol@newphysnet1:~/graphs/")

if __name__ == "__main__":
    main()
