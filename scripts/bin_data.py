import numpy as np
import pandas as pd
import sys
import math
import csv
import os
import glob
import analysis_tools
import itertools
import subprocess

def bin_single_column_data(a, start_bin):
    num_of_bins = int(math.log2(a.size+1))
    
    array = np.empty((num_of_bins,2), dtype=float)
    # here i, the bin index, runs from 0 to num_of_bins-1
    for i in range(start_bin, num_of_bins):
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
    return array
    
def bin_data(a, l, start_bin):
    binned_data = pd.DataFrame(columns=a.drop(['Energy','Magnetization','swap'],axis=1).columns)
    
    # special columns:
    # Energy
    binned_data_with_error = bin_single_column_data(a['Energy']/(4*int(l)**3), start_bin)
    binned_data['Energy']=binned_data_with_error[:,0]
    binned_data['Energy_err']=binned_data_with_error[:,1]
    # Magnetization
    binned_data_with_error = bin_single_column_data(a['Magnetization'].abs(), start_bin)
    binned_data['|Magnetization|']=binned_data_with_error[:,0]
    binned_data['|Magnetization|_err']=binned_data_with_error[:,1]
    # Magnetization^2
    binned_data_with_error = bin_single_column_data(a['Magnetization']**2, start_bin)
    binned_data['Magnetization^2']=binned_data_with_error[:,0]
    binned_data['Magnetization^2_err']=binned_data_with_error[:,1]
    # Magnetization^4
    binned_data_with_error = bin_single_column_data(a['Magnetization']**4, start_bin)
    binned_data['Magnetization^4']=binned_data_with_error[:,0]
    binned_data['Magnetization^4_err']=binned_data_with_error[:,1]
    
    for col_name in a.drop(['Energy','Magnetization','swap'],axis=1).columns:
        binned_data_with_error = bin_single_column_data(a[col_name], start_bin)
        binned_data[col_name]=binned_data_with_error[:,0]
        binned_data[col_name+'_err']=binned_data_with_error[:,1]
    
    binned_data['bin'] = binned_data.index
    
    # reorder columns
    move_col_to_front(binned_data, 'Magnetization^4')
    move_col_to_front(binned_data, 'Magnetization^2')
    move_col_to_front(binned_data, '|Magnetization|')
    move_col_to_front(binned_data, 'Energy')
    move_col_to_front(binned_data, 'bin')
    
    return binned_data

    
def move_col_to_front(df, col_name):
    col = df[col_name]
    df.drop(labels=[col_name], axis=1,inplace = True)
    df.insert(0, col_name, col)
    
    return df

def mkdir(path):
    try:
        os.makedirs(path,exist_ok=True)
    except OSError:
        print ("Creation of the directory %s failed" % path)


def read_binned(sim, use_latest=True):
    path='../analysis/'+sim.folderName+'/binned_data/table_'+str(sim.L)+'_'+str(sim.L)+'_'+str(sim.Bex)+'_'+str(sim.T)+'_'+str(sim.mech)+'_'+'*'+'.txt'
    #print(glob.glob(path))
    arrays=[]
    seeds=[]
    max_bins=0
    min_bins=0
    first_iter=True
    for fname in glob.glob(path):
        curr_array=np.genfromtxt(fname,skip_header=1)
        max_bins = len(curr_array) if len(curr_array)>max_bins else max_bins
        min_bins = len(curr_array) if len(curr_array)<min_bins or first_iter else min_bins
        first_iter=False
        arrays.append(curr_array)
        seeds.append(fname.split("_")[-1].split(".")[0])    # extract and save seed from file name
    if abs(max_bins-min_bins)>1:
        print('one of the simulations might be lagging behind: max_bins=%s, min_bins=%s'%(max_bins,min_bins))
    if use_latest:
        # remove any runs that have length smaller than max_bins, so only the most advanced runs are used
        arrays = list(filter(lambda x: len(x)==max_bins, arrays))
        #print(list(map(len,arrays)))
    else:
        # remove last bins after min_bin so that all arrays have same dimensions
        arrays = [np.resize(a, (min_bins,a.shape[1])) for a in arrays]
        #print(list(map(len,arrays)))
    
    return (np.mean(arrays, axis=0), np.sqrt(np.std(arrays, axis=0)/(len(arrays)-1)))
    
def read_binned_data(sim, use_latest=False, use_bin=-1):
    """ Get the binned data as a pandas DataFrame 
    """
    path='../analysis/'+sim.folderName+'/binned_data/table_'+str(sim.L)+'_'+str(sim.L)+'_'+str(sim.Bex)+'_'+str(sim.T)+'_'+str(sim.mech)+'_'+'*'+'.txt'

    file_list = glob.glob(path) # list of all files that match the sim parameters
    arrays=[]
    
    # iterate over seeds (ind. runs):
    for fname in file_list:
        y = analysis_tools.get_table_data_by_fname(fname, print_prog=True)
        y['seed']=fname.split("_")[-1].split(".")[0]    # extract seed from file name
        arrays.append(y)
        
    all_tables = pd.concat(arrays)
    
    if use_bin==-1:
        # find last bin
        if use_latest:
            # this means use latest bins regardless of how many have finished
            last_bin=all_tables.index.max()
        else:
            count_rows = all_tables.groupby(all_tables.index).count()
            # last bin is the latest one where the number of independent runs it consists of is the same as the number bin 0 consists of
            last_bin = count_rows[count_rows<count_rows.iloc[0]].idxmax()-1
            
            if last_bin.isnull().all():
                last_bin = all_tables.index.max()
            else:
                
                if not np.all(last_bin==last_bin[0]):
                    raise Exception("Something wrong with finding the last bin. different observables seem to have different last bins: \n" + str(last_bin[0]))
                last_bin=last_bin[0]
        use_bin=last_bin
    return all_tables.loc[use_bin]
    

def bin_by_fname(fname, fname_bin, l, start_bin=0):
    y = analysis_tools.get_table_data_by_fname(fname)
#    pd.read_csv(fname, delim_whitespace=True, error_bad_lines=False, index_col='index', comment='#')
    
    binned_data=bin_data(y, l, start_bin)
    
    print(binned_data)
    #binned_data.to_csv(fname_bin,index_label='bin')
    with pd.option_context('display.float_format', '{:0.10f}'.format):
        with open(fname_bin,'w') as outfile:
            binned_data.to_string(outfile,index=False)


def last_index(fname):
    try:
        # get last line from the binned_data file (and suppress stderr)
        with open(os.devnull, 'w') as devnull:
            line = subprocess.check_output(['tail', '-1', fname], stderr=devnull)
        line = line.split()
        return int(line[0]) + 1
    except subprocess.CalledProcessError:
        # most probably the file just hasn't been created yet, but any other problem just means we bin the data from the beginning
        return 0
    except ValueError as e:
        # probably the file has been created but has no data yet (just headers)
        return -1 


def main_bin(simulations):
    

    for sim in simulations.itertuples(index=False):
        mkdir('../analysis/'+sim.folderName+'/binned_data')
        path='../analysis/'+sim.folderName+'/table_'+str(sim.L)+'_'+str(sim.L)+'_'+str(sim.Bex)+'_'+str(sim.T)+'_'+str(sim.mech)+'_'+'*'+'.txt'
        #print(glob.glob(path))
        for fname in glob.glob(path):
            tmp_fname_bin = fname.split('/')
            tmp_fname_bin.insert(3,'binned_data')
            fname_bin='/'.join(tmp_fname_bin)

            # first check if there are enough data for a new bin
            last_bin = last_index(fname_bin)
            last_sample = last_index(fname)
            if last_sample>=0 and (2*(2**(last_bin+1) - 1) <= last_sample): 
                print('rewriting bins for simulation: %s. seed: %s' % (str(sim),fname.split("_")[-1].split(".")[0]))
                bin_by_fname(fname, fname_bin, sim.L)
            else:
                print('not writing bins for simulation: %s. seed: %s' % (str(sim),fname.split("_")[-1].split(".")[0]))
                # not enough data for another bin
                pass
            
def main():
    L = sys.argv[4:]
    Bex = sys.argv[1]
    folderName = sys.argv[2]
    mech = sys.argv[3]
    
    simulations = analysis_tools.get_simulations(L, folderName, Bex, mech)
    main_bin(simulations) 
    
       
if __name__ == "__main__":
    main()
