"""An assortment of methods used for analyzing results from the main MC Java simulation.
"""
import csv
import os
import itertools
import numpy as np
import pandas as pd
import config
from collections import namedtuple


def is_list_type(x):
    """check whether the input is a type or list (list, tuple or set)."""
    return isinstance(x,tuple) or isinstance(x,list) or isinstance(x,set)


def listify(x):
    """make a variable into a list if it isn't already in one, to iterate over later."""
    if is_list_type(x):
        return x
    else:
        return [x]


def unlistify(x):
    """return first item from a list or just the input if it is not a list."""
    if is_list_type(x):
        return x[0]
    else:
        return x


def get_simulations(all_L, all_folderName, all_Bex, all_mech, T=None):
    """
    return a DataFrame with all simulations matching a combination of (L, folderName, Bex, mech).

    :param all_L: system sizes
    :param all_folderName: project name
    :param all_Bex: external magnetic field
    :param all_mech: mechanism flag
    :param T: temperature(s), optional -- if not given all temperatures will be returned.
                                        if given, the closest temperatures to those given will be returned.
    :return: pandas DataFrame with all simulations matching any combination of the given parameters
    """
    # listify all input variables so that we can use itertooles.product
    all_L=listify(all_L)
    all_folderName=listify(all_folderName)
    all_Bex=listify(all_Bex)
    all_mech=listify(all_mech)

    simulations=pd.DataFrame(columns=['L','folderName','Bex','mech','T'])
    # make a product of all possible combinations of the received parameters
    for x in itertools.product(all_L, all_folderName, all_Bex, all_mech):
        temp_schedule_file_name = '../' + config.system_name + '/temperature_schedules/temp_schedule_' + str(x[0]) + '_' + str(x[0]) + '_' + str(x[1]) + '_' + str(x[2]) + '_' + x[3] + '.txt'
        if os.path.exists(temp_schedule_file_name):
            # read the temperatures from the temp_schedule file
            with open(temp_schedule_file_name,'r') as temperature_schedule_file:
                reader = csv.reader(temperature_schedule_file)
                temp_list=list(reader)
                if not temp_list:
                    print(temp_schedule_file_name + " is empty.")
                    continue
            curr_simulations=pd.DataFrame(columns=simulations.columns)
            if T is None:   # no temperature given
                curr_simulations['T']=[float(i) for i in temp_list[0]]
            else:
                # make list of temperatures closest to the ones given in the list T[]
                temperatures_in_schedule=[float(i) for i in temp_list[0]]
                closest_temperatures = [min(temperatures_in_schedule, key=lambda x:abs(x-curr_T)) for curr_T in T]
                curr_simulations['T']=closest_temperatures
            curr_simulations[['L','folderName','Bex','mech']] = x
            simulations = simulations.append(curr_simulations)
        else:
            # combinations that do not correspond to an existing temperature schedule file are just ignored
            pass
    simulations.reset_index(inplace=True)
    print(simulations)
    return simulations


def get_simulation(L, folderName, Bex, mech, T):
    """
    return a namedtuple of a single simulation matching the given parameters.

    :param all_L: system sizes
    :param all_folderName: project name
    :param all_Bex: external magnetic field
    :param all_mech: mechanism flag
    :param T: temperature -- the closest temperature will be returned
    :return: a namedtuple of a single simulation matching the given parameters
    """
    # unlistify all input variables
    L=unlistify(L)
    folderName=unlistify(folderName)
    Bex=unlistify(Bex)
    mech=unlistify(mech)
    T=unlistify(T)
    
    simulation_tuple = namedtuple('Simulation', ['L', 'folderName','Bex','mech','T'])
    
    temp_schedule_file_name = '../' + config.system_name + '/temperature_schedules/temp_schedule_' + str(L) + '_' + str(L) + '_' + str(folderName) + '_' + str(Bex) + '_' + str(mech) + '.txt'
    if os.path.exists(temp_schedule_file_name):
        with open(temp_schedule_file_name,'r') as temperature_schedule_file:
            reader = csv.reader(temperature_schedule_file)
            temp_list=list(reader)
        
        temperatures=[float(i) for i in temp_list[0]]
        
        # get the closest temperature from the temperature schedule
        closest_temperature = min(temperatures, key=lambda x:abs(x-T))  # not the most efficient way, as temperatures should be sorted, but good enough for this
        
        # create the simulation namedtuple
        sim = simulation_tuple(L=L, folderName=folderName, Bex=Bex, mech=mech, T=closest_temperature)
    else:
        raise Exception("No temp_schedule found matching the given simulation parameters.")
    print(sim)
    return sim


def get_table_data_by_fname(fname, print_prog=True):
    """
    Get data from a table in a txt file with the given name.

    :param fname: name of the file to read
    :param print_prog: print which file was read and how many lines were read from it
    :return: pandas DataFrame with the data from the txt file
    """
    try:
        # first read the column names from the file
        col_names = pd.read_csv(fname, delim_whitespace=True, error_bad_lines=False, comment='#', nrows=0).columns
        # decide what type of data is expected in each column for more efficient reading
        types_dict = {'index': int, 'bin' : int, 'n' : int}
        types_dict.update({col: np.float64 for col in col_names if col not in types_dict})
        # the index column can be one of these three, depending which type of file is being read
        possible_index_column_names = ['index', 'bin', 'n']
        index_col = [name for name in possible_index_column_names if name in col_names]
        if len(index_col) != 1:
            raise Exception("Error finding index column. Possible matches: " + str(index_col))

        # read the data
        y = pd.read_csv(fname, delim_whitespace=True, error_bad_lines=False, index_col=index_col, comment='#',dtype=types_dict)

        if not y.index.is_unique:
            # non unique index read from file.
            # This is rare but could happen if the simulation is stopped exactly between the point at which it writes
            # the results to file and between the point at which it writes the checkpoint.
            # if this happens one might need to manually removed some lines from the results file
            print("Non unique indices found in file %s. Removing duplicates, but the file should also be inspected manually." % fname)
            y = y.reset_index().drop_duplicates(keep='first').set_index('index')
    except Exception as e:
        print(fname + ' is empty')
        print(e)
        y=[]
    else:
        if print_prog: print('read %s : %s'%(fname,y.shape[0]))

    return y


def get_table_data(T, L, folderName, Bex, mech, seed, folder):
    """Get data from a table in a txt file matching the given parameters."""
    fname=folder+'/'+folderName+'/table_'+str(L)+'_'+str(L)+'_'+str(Bex)+'_'+str(T)+'_'+mech+'_'+str(seed)+'.txt'
    return get_table_data_by_fname(fname)


def bootstrap(func, param, nboot, arr=None):
    """bootstrap the given function. if arr is given, it is used as the input.
    otherwise we assume func reads its input on its own and samples it as required."""
    res=[]
    for n in range(nboot):
        if arr is not None:
            ret = func(np.random.choice(arr, replace=True), *param)
        else:
            ret = func(*param)
        res.append(ret)
    return np.array(res)


def main():
    print(get_simulations([4,5],['parallel_test'],[0.0],['false']))


if __name__ == "__main__":
    main()
