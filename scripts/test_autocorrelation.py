import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import pandas as pd
import math
import collections
import sys
import warnings
import os
import config

def calc_error_correction(arr, iter):
    binned_arr=np.copy(arr)
    error_est_arr = np.zeros((iter,2))
    for i in range(iter):
        M=binned_arr.size
        # error_est_arr[i,0]=math.sqrt(np.sum((binned_arr-np.mean(binned_arr))**2)/(M*(M-1)))
        # error_est_arr[i,1]=calc_error_correction_err(binned_arr)
        error_est_arr[i,0], error_est_arr[i,1] = calc_error_correction_err(binned_arr)
        binned_arr=get_binned_array(binned_arr)
    return error_est_arr
    
def calc_error_correction_err(arr):
    boot_results=[]
    M=arr.size
    for boot_index in range(500):
        resampled_arr=np.random.choice(arr,size=arr.size,replace=True)
        boot_results.append(math.sqrt(np.sum((resampled_arr-np.mean(resampled_arr))**2)/(M*(M-1)))) 
        #print(boot_index)
    boot_results=np.array(boot_results)
    return np.mean(boot_results), np.std(boot_results,0)
    
    

def get_binned_array(arr):
    if arr.size%2==0:
        return 0.5*(arr[::2]+arr[1::2])
    else:
        return 0.5*(arr[:-1:2]+arr[1::2])


def calc_asymptotic_err(arr):
    last_level=int(math.log2(arr.size))
    binned_arr=np.copy(arr)
    binning_level=0
    converged=False
    conv_queue_lower = []
    conv_queue_upper = []
    
    # first
    first_mean, first_err = calc_error_correction_err(binned_arr)
    conv_queue_lower.append(first_mean - first_err)
    conv_queue_upper.append(first_mean + first_err)
    
    while not converged and binning_level<=last_level and binned_arr.size >= 16:
        converged = test_convergence(conv_queue_lower,conv_queue_upper, 4)
        #print(conv_queue_lower)
        #print(conv_queue_upper)
        #print(converged)
        
        if not converged:
            binned_arr=get_binned_array(binned_arr)
            binning_level+=1
            
            curr_mean, curr_err = calc_error_correction_err(binned_arr)
            conv_queue_lower.append(curr_mean - curr_err)
            conv_queue_upper.append(curr_mean + curr_err)
    
    if converged:
        conv_value = max(conv_queue_lower[-4:])+0.5*(min(conv_queue_upper[-4:])-max(conv_queue_lower[-4:]))
        return (conv_value,first_mean)
    else:
        return (-1,first_mean)

def get_uncorrelated_binning_level(arr, min_len, min_consec_bins, error_arr=None):
    if min_consec_bins<=1:
        warnings.warn("Autocorrelation analysis did not converge at all. Using original timeseries... ",stacklevel=2)
        return arr

    if error_arr is not None:
        arr=np.zeros(2**error_arr.shape[0])    # dummy array
        print(arr.size)

    last_level=int(math.log2(arr.size))
    binned_arr=np.copy(arr)
    binning_level=0
    converged=False
    conv_queue_lower = []
    conv_queue_upper = []
    
    # first
    if error_arr is None:
        first_mean, first_err = calc_error_correction_err(binned_arr)
    else:
        first_mean, first_err = (error_arr[0,0], error_arr[0,1])

    conv_queue_lower.append(first_mean - first_err)
    conv_queue_upper.append(first_mean + first_err)
    
    while not converged and binning_level<=last_level and binned_arr.size >= min_len:
        converged = test_convergence(conv_queue_lower,conv_queue_upper, min_consec_bins)
        #print(conv_queue_lower)
        #print(conv_queue_upper)
        #print(converged)
        
        if not converged:
            binned_arr=get_binned_array(binned_arr)
            binning_level+=1
            
            if error_arr is None:
                curr_mean, curr_err = calc_error_correction_err(binned_arr)
            else:
                curr_mean, curr_err = (error_arr[binning_level,0], error_arr[binning_level,1])
            conv_queue_lower.append(curr_mean - curr_err)
            conv_queue_upper.append(curr_mean + curr_err)
    
    if converged:
        #binned_arr=np.copy(arr)
        #for i in range(binning_level-2):
        #    binned_arr=get_binned_array(binned_arr)
        #return binned_arr
        if error_arr is not None:
            return (max(conv_queue_lower[-4:])+0.5*(min(conv_queue_upper[-4:])-max(conv_queue_lower[-4:])), binning_level+1-min_consec_bins) 
        return binning_level+1-min_consec_bins
    else:
        warnings.warn("Autocorrelation analysis did not converge. Using less strict criterion for convergence (%s)."%(min_consec_bins-1),stacklevel=2)
        return get_uncorrelated_binning_level(arr,min_len,min_consec_bins-1)

def get_nth_level_binned_arr(arr, n):
    binned_arr=np.copy(arr)
    for i in range(n):
       binned_arr=get_binned_array(binned_arr)
    return binned_arr

# reads files from folder and performs binning analysis on the columns given in "cols_to_copy".
# the uncorrelated, processed, equilibrated timeseries of these columns are saved to /tmp/
# look out for warnings as they nofity a transition to less strict convergence criterion for the binning analysis
def save_uncorrelated_timeseries(T, cols_to_copy, L, L_equilibrated_min_value, Bex, folderName, mech, folder='../'+config.system_name+'/data/results'):
    for l in L:
        path='/tmp/'+folderName[l]
        try:
            os.makedirs(path,exist_ok=True)
        except OSError:
            print ("Creation of the directory %s failed" % path)
    

        for temperature in T[l]:
            fname=folder+'/'+folderName[l]+'/table_'+str(l)+'_'+str(l)+'_'+str(Bex)+'_'+str(temperature)+'_'+mech+'.txt'
            
            col_names = pd.read_csv(fname, delim_whitespace=True, error_bad_lines=False, comment='#', nrows=0).columns
            types_dict = {'index': int, 'swap': int}
            types_dict.update({col: np.float64 for col in col_names if col not in types_dict})                        
            dest=path+'/table_'+str(l)+'_'+str(l)+'_'+str(Bex)+'_'+str(temperature)+'_'+mech+'.txt'
            y = pd.read_csv(fname, delim_whitespace=True, error_bad_lines=False, index_col='index', comment='#',
            dtype=types_dict)
            
            print('checking file ' + fname + ' for autocorrelations... ', end = ' ')
            
            # find maximum converged binning level among all columns to test
            max_uncorrelated_binning_level=0
            for to_check_index, to_check_now in enumerate(cols_to_copy):
                if str(to_check_now)=="Energy":
                    a = y[to_check_now]/(4*int(l)**3)
                elif str(to_check_now)=="Magnetization":
                    a=y[to_check_now].abs()
                elif str(to_check_now)=="Magnetization^2":
                    a=y['Magnetization']**2
                else:
                    a=y[to_check_now]
                
                a=a.loc[L_equilibrated_min_value[l]:]
                curr_uncorrelated_binning_level = get_uncorrelated_binning_level(a, 16, 4)
                if curr_uncorrelated_binning_level>max_uncorrelated_binning_level:
                    max_uncorrelated_binning_level=curr_uncorrelated_binning_level
            
            # perform binning up to the maximum binning level for all columns and save
            uncorrelated_df = pd.DataFrame()
            for to_check_now in cols_to_copy:
                if str(to_check_now)=="Energy":
                    a = y[to_check_now]/(4*int(l)**3)
                elif str(to_check_now)=="Magnetization":
                    a=y[to_check_now].abs()
                elif str(to_check_now)=="Magnetization^2":
                    a=y['Magnetization']**2
                else:
                    a=y[to_check_now]
                a=a.loc[L_equilibrated_min_value[l]:]
                
                uncorrelated_df[to_check_now] = get_nth_level_binned_arr(a, max_uncorrelated_binning_level)

                #uncorrelated_df[to_check_now] = a

            print('copying file ' + fname + ' to /tmp/ with binning level %s'%max_uncorrelated_binning_level)
            uncorrelated_df.to_csv(dest, sep='\t', index_label='index')
            
            del y
    
    # when saving uncorrelated data, there is no need for skipping rows, so a dict with all values=1 is returned
    #return {k:1 for k in L}

def calc_tau(arr):
    asymptotic_err, simple_err = calc_asymptotic_err(arr)
    
    if (asymptotic_err<0):
        return -1
    else:        
        return 0.5*((asymptotic_err/simple_err)**2-1)
    
def calc_max_tau(T, to_check, L, L_equilibrated_min_value, Bex, folderName, mech, folder='../'+config.system_name+'/data/results'):
    max_tau_dict = {k:1 for k in L}
    for l in L:
        for temperature_index, temperature in enumerate(T):
            fname=folder+'/'+folderName+'/table_'+str(l)+'_'+str(l)+'_'+str(Bex)+'_'+str(temperature)+'_'+mech+'.txt'
            y = pd.read_csv(fname, delim_whitespace=True, error_bad_lines=False, index_col='index', comment='#')
            for to_check_index, to_check_now in enumerate(to_check):
                if str(to_check_now)=="Energy":
                    a = y[to_check_now]/(4*int(l)**3)
                elif str(to_check_now)=="Magnetization":
                    a=y[to_check_now].abs()
                elif str(to_check_now)=="Magnetization^2":
                    a=y['Magnetization']**2
                else:
                    a=y[to_check_now]
                
                a=a.loc[L_equilibrated_min_value[l]:]
                
                real_tau=calc_tau(a)
                curr_tau=math.ceil(real_tau)
                if (curr_tau<1):
                    print('not enough samples for autocorrelation test. L=%s, T=%s'%(l,temperature))
                if curr_tau > max_tau_dict[l]:
                    max_tau_dict[l]=curr_tau
                print("L=%s T=%s obs=%s : %s (%s)"%(l,temperature,to_check_now,curr_tau,real_tau))
            
            del y    
    print(max_tau_dict)
    return max_tau_dict

def test_convergence(conv_queue_lower, conv_queue_upper, min_consec_bins):
    # last 4 bins agree within error bars
    return len(conv_queue_lower)>=min_consec_bins and max(conv_queue_lower[-min_consec_bins:]) <= min(conv_queue_upper[-min_consec_bins:])

def plot_binning_analysis(arr,obs):
    from matplotlib.ticker import MaxNLocator
    last_level=int(math.log2(arr.size)-2)
    error_corrections=calc_error_correction(arr,last_level)
    error_and_value=get_uncorrelated_binning_level(arr,16,4,error_arr=error_corrections)
    plt.errorbar(range(last_level),error_corrections[:,0],yerr=error_corrections[:,1],fmt='.-',capsize=2,fillstyle='none',mew=1.0,linewidth=1.0)
    plt.hlines(error_and_value[0], plt.xlim()[0], plt.xlim()[1], colors='r')
    plt.axvline(x=error_and_value[1], linewidth=0.5, color='k')
    plt.title('%s Binned standard error vs. Binning level'%obs)
    plt.xlabel(r'Binning level $l$')
    plt.ylabel(r'%s $\Delta^{(l)}$'%obs)
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.gcf().savefig('binning_fig.png',dpi=400)

def main():
    obs = 'Energy'
    L = sys.argv[4]
    h_ex = sys.argv[1]
    folderName = sys.argv[2]
    mech = sys.argv[3]
    T=sys.argv[5]
    
    fname='../'+config.system_name+'/data/results/'+folderName+'/table_'+str(L)+'_'+str(L)+'_'+str(h_ex)+'_'+str(T)+'_'+mech+'.txt'
    y = pd.read_csv(fname, delim_whitespace=True, error_bad_lines=False, index_col='index', comment='#')
    print('done reading file!')
    
    arr=y[obs].values
    arr=arr[2**15-1:]
    
    plot_binning_analysis(arr,obs)
    
    #print(calc_tau(arr))
    


if __name__ == "__main__":
    main()
