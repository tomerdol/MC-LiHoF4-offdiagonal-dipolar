import numpy as np
import pandas as pd
import sys
import math
import csv
import scipy.stats as stats
import warnings

# return how many bins from the end are within the error bars of one another
# (>=3 is considered equilibrated)
# --------
# a: is the array in which we check equilibration
# min_bin: is the minimum bin from which to start checking (e.g. the first 3 bins being equilibrated does
# not means real equilibration)
# consec_bins_for_equilib: the minimum number of consecutive bins that agree within error bars for
# equilibration
import config


def check(a, min_bin, consec_bins_for_equilib):
    max_num_of_bins=0
    num_of_bins = int(math.log2(a.size+1))
    
    if num_of_bins-consec_bins_for_equilib<min_bin:
        raise Exception("minimum bin is larger than the number of bins minus the consec_bins_for_equilib so equilibration cannot occur. ")
    
    if num_of_bins>max_num_of_bins:
        max_num_of_bins=num_of_bins
    array = np.empty((num_of_bins,2), dtype=float)
    # here i, the bin index, runs from 0 to num_of_bins-1
    for i in range(0,num_of_bins):
        N = a.loc[2**(i)-1:2**(i+1)-2].size # number of samples in current bin
        if N > 0:
            # pandas .loc[] function is inclusive on both ends!
            array[i,0] = np.sum(a.loc[2**(i)-1:2**(i+1)-2]) / N
            array[i,1] = np.sum(a.loc[2**(i)-1:2**(i+1)-2]**2) / N
            if (array[i,1] - array[i,0]**2 > 0):
                array[i,1] = math.sqrt((array[i,1] - array[i,0]**2)/(N-1))
            else:
                array[i,1] = 0
        else:
            array[i,0] = None
            array[i,1] = None
        #print('finished bin %s in L=%s'%(i,l))
    
    curr_bin=min_bin
    eq=False
    
    while not eq and curr_bin<=num_of_bins-consec_bins_for_equilib:
        # cond1 is the standard condition, meaning all #consec_bins_for_equilib have overlapping errorbars
        # the condition itself checkes whether the maximum of the lower bound is below the minimum of the upper bounds
        cond1 = (np.max(array[curr_bin:curr_bin+consec_bins_for_equilib,0]-array[curr_bin:curr_bin+consec_bins_for_equilib,1]) 
        < np.min(array[curr_bin:curr_bin+consec_bins_for_equilib,0]+array[curr_bin:curr_bin+consec_bins_for_equilib,1]))
        
        # check if the #consec_bins_for_equilib have monotonically increasing averages
        monotonically_inc = np.all(array[curr_bin+1:curr_bin+consec_bins_for_equilib,0] >= array[curr_bin:curr_bin+consec_bins_for_equilib-1,0])
        # check if the #consec_bins_for_equilib have monotonically decreasing averages
        monotonically_dec = np.all(array[curr_bin+1:curr_bin+consec_bins_for_equilib,0] <= array[curr_bin:curr_bin+consec_bins_for_equilib-1,0])
        # cond2 is that the #consec_bins_for_equilib not have monotonically increasing or decreasing averages
        cond2 = not (monotonically_inc or monotonically_dec)
        
        if cond1 and cond2:
            eq=True
        curr_bin=curr_bin+1
    
    if eq:
        return curr_bin-1+consec_bins_for_equilib-1 # the index of the last bin in the equilibrated sequence
    else:
        #print(eq)
        #print(curr_bin)
        eq_bin_from_end=check_from_end(array)
        if eq_bin_from_end>=consec_bins_for_equilib:
            return num_of_bins-eq_bin_from_end+consec_bins_for_equilib-1
        else:
            return num_of_bins
            
        
def check_from_end(array):
    num_of_bins=np.size(array,0)
    #print(num_of_bins)
    n=1
    eq=True
    #common_range = [ array[-n,0]-array[-n,1], array[-n,0]+array[-n,1] ]
    while eq and n<num_of_bins:
        #extended_range = [ array[-n-1,0]-array[-n-1,1], array[-n-1,0]+array[-n-1,1] ]
        #common_range = [ max(extended_range[0],common_range[0]) , min(extended_range[1], common_range[1]) ]
        
        #if (common_range[0]>common_range[1] or array[-n-1:,0].max() > common_range[1] or array[-n-1:,0].min() < common_range[0]):
        if (np.max(array[-n-1:,0]-array[-n-1:,1]) > np.min(array[-n-1:,0]+array[-n-1:,1])):
            eq=False
        n=n+1
#    print('[%s,%s]'%(common_range[0],common_range[1]))
#    print(n-1)
    return (n-1)
    


def check_equilib(T, to_check, L, Bex, folderName, mech, folder='../' + config.system_name + '/data/results'):
    min_steps_dict={}
    all_equilibrated=True
    for l in L:
        print('L: ' + str(l))
        equilibrated_L = []
        num_of_bins = 0
        max_bin=0
        all_equilibrated_L=True
        for temperature_index, temperature in enumerate(T[l]):
            #print('temperature: ' + str(temperature))
            fname=folder+'/'+folderName[l]+'/table_'+str(l)+'_'+str(l)+'_'+str(Bex)+'_'+str(temperature)+'_'+mech+'.txt'
            col_names = pd.read_csv(fname, delim_whitespace=True, error_bad_lines=False, comment='#', nrows=0).columns
            types_dict = {'index': int, 'swap': int}
            types_dict.update({col: np.float64 for col in col_names if col not in types_dict})
            y = pd.read_csv(fname, delim_whitespace=True, error_bad_lines=False, index_col='index', comment='#',
            dtype=types_dict)
            #print('read csv %s' % fname)
            #to_check = y.drop('swap',axis=1).columns.values
            equilibrated_T = []
            
            for to_check_index, to_check_now in enumerate(to_check):
                if str(to_check_now)=="Energy":
                    a = y[to_check_now]/(4*int(l)**3)
                elif str(to_check_now)=="Magnetization":
                    a=y[to_check_now].abs()
                elif str(to_check_now)=="Magnetization^2":
                    a=y['Magnetization']**2
                else:
                    a=y[to_check_now]
                
                num_of_bins = max(math.floor(math.log2(a.size)), num_of_bins)
                bin=check(a, 6, 3)
                #print(bin)
                equilibrated = bin<num_of_bins-1   # min for equilibration is 3
                max_bin=max(max_bin,bin)
                all_equilibrated &= equilibrated
                all_equilibrated_L &= equilibrated
                equilibrated_T.append(equilibrated)
                #print("L=%s T=%s obs=%s : %s (%s)"%(l,temperature,to_check_now,equilibrated,num_of_bins-bin))
                
            equilibrated_L.append(equilibrated_T)
            del y
        t = pd.DataFrame(data=equilibrated_L,index=T[l],columns=to_check)
        print(t)
        print('min step: ' + (str(math.pow(2,max_bin+1)-1) if all_equilibrated_L else '-'))
        if not all_equilibrated_L:
            warnings.warn("Not all sizes equilibrated. Using only the last bin for those that weren't... ",stacklevel=2)
        min_steps_dict[l] = math.pow(2,max_bin+1)-1 if all_equilibrated_L else int(a.size/2)
    #print('All equilibrated: ' + str(all_equilibrated))
    return min_steps_dict

def main():
    to_check = ['Energy','Magnetization','Magnetization^2','mk2x']
    L = sys.argv[4:]
    Bex = sys.argv[1]
    folderName = sys.argv[2]
    mech = sys.argv[3]

    with open('../' + config.system_name + '/temp_schedule_' + folderName + '_' + str(Bex) + '_' + mech + '.txt','r') as temperature_schedule_file:
    #with open('../temperature_schedule_t.txt','r') as temperature_schedule_file:
        reader = csv.reader(temperature_schedule_file)
        temp_list=list(reader) 
    
    T=temp_list[0]
    T=[float(i) for i in T]
    
    T_dict=dict.fromkeys(L, T)
    folderName_dict=dict.fromkeys(L, folderName)
    check_equilib(T_dict, to_check, L, Bex, folderName_dict, mech)

if __name__ == "__main__":
    main()
