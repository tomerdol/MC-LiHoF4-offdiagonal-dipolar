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
import config

def bin_single_column_data(a, start_bin):
    # samples that do no fill the last bin are practically discarded
    num_of_bins = int(math.log2(a.size+1))
    
    array = np.empty((num_of_bins,2), dtype=float)
    # here i, the bin index, runs from 0 to num_of_bins-1
    for i in range(start_bin, num_of_bins):
        N = a.loc[2**(i)-1:2**(i+1)-2].size # number of samples in current bin
        if N > 0:
            # pandas .loc[] function is inclusive on both ends!
            bin_arr = a.loc[2**(i)-1:2**(i+1)-2].to_numpy()
            if np.isnan(bin_arr).any():
                array[i,0] = np.nan
                array[i,1] = np.nan
            else:
                array[i,0] = np.sum(bin_arr) / N
                array[i,1] = np.sum(bin_arr**2) / N
            if (array[i,1] - array[i,0]**2 > 0):
                array[i,1] = math.sqrt((array[i,1] - array[i,0]**2)/(N-1))
            else:
                array[i,1] = 0
        else:
            array[i,0] = None
            array[i,1] = None
    return array
    
def bin_data(a, l, start_bin):
    binned_data = pd.DataFrame(columns=a.drop(['Energy','Magnetization'],axis=1).columns)
    
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

    for col_name in a.drop(['Energy','Magnetization'],axis=1).columns:
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


def read_binned_txt(sim, use_latest=True):
    path='../' + config.system_name + '/data/results/'+sim.folderName+'/binned_data/table_'+str(sim.L)+'_'+str(sim.L)+'_'+str(sim.Bex)+'_'+str(sim.T)+'_'+str(sim.mech)+'_'+'*'+'.txt'
    #print(glob.glob(path))
    arrays=[]
    seeds=[]
    max_bins=0
    min_bins=0
    first_iter=True
    files = glob.glob(path)
    if not files:
        # no files matching pattern found
        raise Exception("No binned files found matching the given pattern for " + str(sim))
    for fname in files:
        curr_array=np.genfromtxt(fname,skip_header=1)
        max_bins = len(curr_array) if len(curr_array)>max_bins else max_bins
        min_bins = len(curr_array) if len(curr_array)<min_bins or first_iter else min_bins
        first_iter=False
        arrays.append(curr_array)
        seeds.append(fname.split("_")[-1].split(".")[0])    # extract and save seed from file name
    if abs(max_bins-min_bins)>0:
        print('WARNING: one of the simulations might be lagging behind: max_bins=%s, min_bins=%s \n Simulation details: %s \n Seed: %s'%(max_bins,min_bins,sim,seeds[-1]), file=sys.stderr)
    if use_latest:
        # remove any runs that have length smaller than max_bins, so only the most advanced runs are used
        arrays = list(filter(lambda x: len(x)==max_bins, arrays))
        #print(list(map(len,arrays)))
    else:
        # remove last bins after min_bin so that all arrays have same dimensions
        arrays = [np.resize(a, (min_bins,a.shape[1])) for a in arrays]
        #print(list(map(len,arrays)))
    print(str(sim) + ": Using %s independent simulations with %s bins in each." % (len(arrays), len(arrays[0])))
    return (np.mean(arrays, axis=0), np.sqrt(np.std(arrays, axis=0)/(len(arrays)-1)))


def read_binned(sim, use_latest=True):
    fname='../' + config.system_name + '/data/results/'+sim.folderName+'/binned_data/table_'+str(sim.L)+'_'+str(sim.L)+'_'+str(sim.Bex)+'_'+str(sim.mech)+'.h5'
    if not os.path.isfile(fname):
        raise Exception("No binned data file found matching the given pattern for " + str(sim))

    arrays=[]
    seeds=[]
    max_bins=0
    min_bins=0
    first_iter=True
    with pd.HDFStore(fname, mode='r') as hdf_bin:
        for (path, subgroups, subkeys) in hdf_bin.walk():
            if path == '':    # we only need to iterate over the root level (seed groups)
                for subgroup in subgroups:
                    curr_array = hdf_bin.get(subgroup+"/T"+str(sim.T).replace('.','_'))
                    max_bins = len(curr_array) if len(curr_array)>max_bins else max_bins
                    min_bins = len(curr_array) if len(curr_array)<min_bins or first_iter else min_bins
                    first_iter=False
                    arrays.append(curr_array.to_numpy())
                    seeds.append(subgroup[1:])    # extract and save seed from the group name
                    # print(fname + '/' + subgroup + '/T'+str(sim.T).replace('.','_') + ' : ' + str(len(curr_array)))
    if abs(max_bins-min_bins)>0:
        print('WARNING: one of the simulations might be lagging behind: max_bins=%s, min_bins=%s \n Simulation details: %s \n Seed: %s'%(max_bins,min_bins,sim,seeds[-1]), file=sys.stderr)
    if use_latest:
        # remove any runs that have length smaller than max_bins, so only the most advanced runs are used
        arrays = list(filter(lambda x: len(x)==max_bins, arrays))
        #print(list(map(len,arrays)))
    else:
        # remove last bins after min_bin so that all arrays have same dimensions
        arrays = [np.resize(a, (min_bins,a.shape[1])) for a in arrays]
        #print(list(map(len,arrays)))
    print(str(sim) + ": Using %s independent simulations with %s bins in each." % (len(arrays), len(arrays[0])))
    return (np.mean(arrays, axis=0), np.sqrt(np.std(arrays, axis=0)/(len(arrays)-1)))


def read_binned_data_txt(sim, use_latest=False, use_bin=-1):
    """ Get the binned data as a pandas DataFrame
    """
    path='../' + config.system_name + '/data/results/'+sim.folderName+'/binned_data/table_'+str(sim.L)+'_'+str(sim.L)+'_'+str(sim.Bex)+'_'+str(sim.T)+'_'+str(sim.mech)+'_'+'*'+'.txt'

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
    if 'eq_bin' in sim._fields:
        # equilibration testing has occurred and there is an equilibrated bin.
        if int(sim.eq_bin) > use_bin:
            print('WARNING: using data before equilibration: Bin used: %s. Simulation details: %s'%(use_bin, sim))

    return all_tables.loc[use_bin]


def read_binned_data(sim, use_latest=False, use_bin=-1):
    """ Get the binned data as a pandas DataFrame 
    """
    fname='../' + config.system_name + '/data/results/'+sim.folderName+'/binned_data/table_'+str(sim.L)+'_'+str(sim.L)+'_'+str(sim.Bex)+'_'+str(sim.mech)+'.h5'
    if not os.path.isfile(fname):
        raise Exception("No binned data file found matching the given pattern for " + str(sim))

    arrays=[]

    with pd.HDFStore(fname, mode='r') as hdf_bin:
        # iterate over seeds (ind. runs):
        for (path, subgroups, subkeys) in hdf_bin.walk():
            if path=='':    # we only need to iterate over the root level (seed groups)
                for subgroup in subgroups:
                    try:
                        y = hdf_bin.get(subgroup+"/T"+str(sim.T).replace('.','_'))
                    except KeyError:
                        # this means that a seed group exists (otherwise it will not be in subgroups)
                        # but does not include the required temperature dataset, so we just skip it.
                        continue
                    print(fname + '/' + subgroup + '/T'+str(sim.T).replace('.','_') + ' : ' + str(len(y)))
                    y['seed']=subgroup[1:]    # extract seed from file name
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
    if 'eq_bin' in sim._fields:
        # equilibration testing has occurred and there is an equilibrated bin.
        if int(sim.eq_bin) > use_bin:
           print('WARNING: using data before equilibration: Bin used: %s. Simulation details: %s'%(use_bin, sim))

    return all_tables.loc[use_bin]


def bin_by_fname_txt(fname, fname_bin, l, start_bin=0):
    y = analysis_tools.get_table_data_by_fname(fname)
    #    pd.read_csv(fname, delim_whitespace=True, error_bad_lines=False, index_col='index', comment='#')

    binned_data=bin_data(y, l, start_bin)

    print(binned_data)
    #binned_data.to_csv(fname_bin,index_label='bin')
    with pd.option_context('display.float_format', '{:0.10f}'.format):
        with open(fname_bin,'w') as outfile:
            binned_data.to_string(outfile,index=False)


def bin_by_fname(fname, hdf_bin, T, seed, l, start_bin=0):
    y = analysis_tools.get_table_data_by_fname(fname)
#    pd.read_csv(fname, delim_whitespace=True, error_bad_lines=False, index_col='index', comment='#')
    
    binned_data = bin_data(y, l, start_bin)
    print(binned_data)

    hdf_bin.put(get_dataset_name(seed, T), binned_data, format='table', append=False)


def reverse_readline(filename, buf_size=8192):
    """A generator that returns the lines of a file in reverse order"""
    with open(filename) as fh:
        segment = None
        offset = 0
        fh.seek(0, os.SEEK_END)
        file_size = remaining_size = fh.tell()
        while remaining_size > 0:
            offset = min(file_size, offset + buf_size)
            fh.seek(file_size - offset)
            buffer = fh.read(min(remaining_size, buf_size))
            remaining_size -= buf_size
            lines = buffer.split('\n')
            # The first line of the buffer is probably not a complete line so
            # we'll save it and append it to the last line of the next buffer
            # we read
            if segment is not None:
                # If the previous chunk starts right from the beginning of line
                # do not concat the segment to the last line of new chunk.
                # Instead, yield the segment first
                if buffer[-1] != '\n':
                    lines[-1] += segment
                else:
                    yield segment
            segment = lines[0]
            for index in range(len(lines) - 1, 0, -1):
                if lines[index]:
                    yield lines[index]
        # Don't yield None if the file was empty
        if segment is not None:
            yield segment


def last_index(fname):
    try:
        # get last line from the binned_data file (and suppress stderr)
        for line in reverse_readline(fname):
            if line[0] != '#':
                return int(line.split()[0])
        # ret=0
        # while ret == 0:
        #     num_of_lines_from_end = 10
        #     with open(os.devnull, 'w') as devnull:
        #         lines = subprocess.check_output(['tail', '-'+str(num_of_lines_from_end), fname], stderr=devnull)
        #     lines = lines.splitlines(False)[-1].split() # Split to lines and then by whitespace
        #     ret = int(lines[0])
    except subprocess.CalledProcessError:
        # most probably the file just hasn't been created yet, but any other problem just means we bin the data from the beginning
        return 0
    except ValueError:
        # probably the file has been created but has no data yet (just headers)
        return -1
    except OSError:
        return -1


def main_bin_txt(simulations):
    for sim in simulations.itertuples(index=False):
        # if specific project name is given, the binned_data folder can be created just once, now
        if sim.folderName != '*': mkdir('../' + config.system_name + '/data/results/'+sim.folderName+'/binned_data')

        path='../' + config.system_name + '/data/results/'+sim.folderName+'/table_'+str(sim.L)+'_'+str(sim.L)+'_'+str(sim.Bex)+'_'+str(sim.T)+'_'+str(sim.mech)+'_'+'*'+'.txt'
        #print(glob.glob(path))
        for fname in glob.glob(path):
            # if no project name is given, the existence of a corresponding binned_data folder must be verified for each file
            if sim.folderName == '*': mkdir(os.path.dirname(fname) + '/binned_data')

            tmp_fname_bin = fname.split('/')
            tmp_fname_bin.insert(5,'binned_data')
            fname_bin='/'.join(tmp_fname_bin)

            # first check if there are enough data for a new bin
            last_bin = last_index(fname_bin)
            last_sample = last_index(fname)
            if last_sample<0:
                print('No data in simulation: %s. seed: %s. Skipping.' % (str(sim),fname.split("_")[-1].split(".")[0]))
            else:
                if (2*(2**(last_bin+1) - 1) <= last_sample):
                    print('rewriting bins for simulation: %s. seed: %s' % (str(sim),fname.split("_")[-1].split(".")[0]))
                    bin_by_fname_txt(fname, fname_bin, sim.L)
                else:
                    print('not writing bins for simulation: %s. last bin: %s. seed: %s' % (str(sim), last_bin, fname.split("_")[-1].split(".")[0]))
                    # not enough data for another bin
                    pass


def main_bin(simulations):
    for group, grouped_sim in simulations.groupby(['Bex','L','folderName','mech']):
    # for sim in simulations.itertuples(index=False):
        Bex=str(group[0])
        L=str(group[1])
        folderName=str(group[2])
        mech=str(group[3])
        mkdir('../' + config.system_name + '/data/results/'+folderName+'/binned_data')
        hdf_fname = '../' + config.system_name + '/data/results/'+folderName+'/binned_data/table_'+L+'_'+L+'_'+Bex+'_'+mech+'.h5'
        with pd.HDFStore(hdf_fname) as hdf_bin:
            # this basically loops over temperatures
            for sim in grouped_sim.itertuples(index=False):
                path='../' + config.system_name + '/data/results/'+sim.folderName+'/table_'+str(sim.L)+'_'+str(sim.L)+'_'+str(sim.Bex)+'_'+str(sim.T)+'_'+str(sim.mech)+'_'+'*'+'.txt'

                for fname in glob.glob(path):
                    seed = fname.split('/')[-1].split('_')[-1].split('.')[0]
                    # first check if there are enough data for a new bin
                    if get_dataset_name(seed, sim.T) not in hdf_bin:
                        last_bin=0
                    else:
                        last_bin = hdf_bin.get_storer(get_dataset_name(seed, sim.T)).nrows - 1  # last bin is the total number of rows minus 1
                    last_sample = last_index(fname)
                    if last_sample<0:
                        print('No data in simulation: %s. seed: %s. Skipping.' % (str(sim),fname.split("_")[-1].split(".")[0]))
                    else:
                        if (2*(2**(last_bin+1) - 1) <= last_sample):
                            print('rewriting bins for simulation: %s. seed: %s' % (str(sim),fname.split("_")[-1].split(".")[0]))
                            bin_by_fname(fname, hdf_bin, sim.T, seed, L)
                        else:
                            print('not writing bins for simulation: %s. last bin: %s. seed: %s' % (str(sim), last_bin, fname.split("_")[-1].split(".")[0]))
                            # not enough data for another bin
                            pass


def get_dataset_name(seed, T):
    return '/s'+seed+'/T'+str(T).replace('.','_')


def parse_arguments():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description="Analyzes Monte Carlo results and plots correlation length curves.", formatter_class=ArgumentDefaultsHelpFormatter, parents=[config.parse_arguments()], conflict_handler='resolve')
    args = parser.parse_args()
    config.system_name = args.system_name
    return args


def main():
    args = parse_arguments()
    L = args.L
    h_ex = args.h_ex
    mech = args.mech
    folderName = args.folder_list

    # if sys.argv[2] == '*':
    #     simulations=pd.DataFrame(columns=['L','folderName','Bex','mech','T'])
    #     simulations.loc[0] = '*'
    # else:
    #     L = sys.argv[5:]
    #     Bex = sys.argv[2]
    #     folderName = sys.argv[3]
    #     mech = sys.argv[4]
    simulations = analysis_tools.get_simulations(L, folderName, h_ex, mech)
    main_bin(simulations)
    # for testing:
    # print(read_binned_data(list(simulations.itertuples())[3]))
       
if __name__ == "__main__":
    main()
