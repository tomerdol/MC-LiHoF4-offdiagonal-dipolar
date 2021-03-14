import matplotlib

import config

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import analysis_tools
import sys
import glob


def plot_relax_func(simulations, to_plot):
    nonlinear_relax_funcs=[]
    
    for sim in simulations.itertuples():
        path='../' + config.system_name + '/data/results/'+sim.folderName+'/table_'+str(sim.L)+'_'+str(sim.L)+'_'+str(sim.Bex)+'_'+str(sim.T)+'_'+str(sim.mech)+'_'+'*'+'.txt'
        
        file_list = glob.glob(path)
        n=len(file_list)    # number of independent runs
        arrays=[]
        # iterate over seeds (ind. runs)
        
        
        for fname in file_list:
            y = analysis_tools.get_table_data_by_fname(fname)
            
            arrays.append(y)
            
        all_tables = pd.concat(arrays)
        
        # create new columns w/ abs and squared magnetization (otherwise averages are meaningless)
        all_tables['|m|']=all_tables['Magnetization'].abs()
        all_tables['m2']=all_tables['Magnetization']**2

        ind_runs_means = all_tables.groupby(all_tables.index).mean()
        print(ind_runs_means)
        
        nonlinear_relax_func = (ind_runs_means -  ind_runs_means.iloc[-1])/(ind_runs_means.iloc[0] - ind_runs_means.iloc[-1])
        #nonlinear_relax_func = (ind_runs_means -  ind_runs_means.iloc[-10000:].mean())/(ind_runs_means.iloc[0] - ind_runs_means.iloc[-10000:].mean())
        nonlinear_relax_func['Energy'].loc[:900].cumsum().plot()
        plt.gcf().savefig('../' + config.system_name + '/figures/plot_relax_time_cumsum.png',dpi=300)
        
        nonlinear_relax_func = pd.concat((nonlinear_relax_func,pd.DataFrame(sim._asdict(),index=nonlinear_relax_func.index)),axis=1)
        
        nonlinear_relax_funcs.append(nonlinear_relax_func)
    
    
    nonlinear_relax_funcs = pd.concat(nonlinear_relax_funcs)
    
    
    fig, ax = plt.subplots(figsize=(8,6))
    for label, df in nonlinear_relax_funcs.groupby(['Bex','L','folderName','mech','T']):
        df.head(1500).plot(y=to_plot, ax=ax, label=[lbl+' '+str(label) for lbl in to_plot])
    
    ax.set_xlabel('t [MCS]')
    
    try:
        plt.tight_layout()
    except:
        pass
    
    
    
    fig.savefig('../' + config.system_name + '/figures/plot_relax_func.png',dpi=300)


def plot_relax_times(simulations, to_plot):

    nonlinear_relax_times = []
    # for now, assuming this just means iterate by temperature:
    for i, sim in enumerate(simulations.itertuples()):
        path='../' + config.system_name + '/data/results/'+sim.folderName+'/table_'+str(sim.L)+'_'+str(sim.L)+'_'+str(sim.Bex)+'_'+str(sim.T)+'_'+str(sim.mech)+'_'+'*'+'.txt'
        
        file_list = glob.glob(path)
        n=len(file_list)    # number of independent runs
        arrays=[]
        # iterate over seeds (ind. runs)
        
        
        for fname in file_list:
            y = analysis_tools.get_table_data_by_fname(fname, print_prog=True)
            
            arrays.append(y)
            
        all_tables = pd.concat(arrays)
        
        # create new columns w/ abs and squared magnetization (otherwise averages are meaningless)
        all_tables['|m|']=all_tables['Magnetization'].abs()
        all_tables['m2']=all_tables['Magnetization']**2

        ind_runs_means = all_tables.groupby(all_tables.index).mean()
        
        nonlinear_relax_func = (ind_runs_means -  ind_runs_means.iloc[-1])/(ind_runs_means.iloc[0] - ind_runs_means.iloc[-1])
        #nonlinear_relax_func = (ind_runs_means -  ind_runs_means.iloc[-10000:].mean())/(ind_runs_means.iloc[0] - ind_runs_means.iloc[-10000:].mean())
        
        nonlinear_relax_time = nonlinear_relax_func.sum()
        print(nonlinear_relax_time)
        #nonlinear_relax_time['T']=sim.T
        nonlinear_relax_time=nonlinear_relax_time.append(pd.Series(sim._asdict()))
        print(nonlinear_relax_time)
        nonlinear_relax_times.append(nonlinear_relax_time)
    
    nonlinear_relax_times = pd.DataFrame(nonlinear_relax_times)
    
    #axes = nonlinear_relax_times.plot(x='T', y=to_plot, subplots=True, sharex=True)
    
    fig, ax = plt.subplots(figsize=(8,6))
    for label, df in nonlinear_relax_times.groupby(['Bex','L','folderName','mech']):
        df.plot(x='T',y=to_plot, ax=ax, label=[lbl+' '+str(label) for lbl in to_plot])
    
    ax.set_xlabel('T')
    
    try:
        plt.tight_layout()
    except:
        pass
    
    
    
    fig.savefig('../' + config.system_name + '/figures/plot_relax_times.png',dpi=300)

def main():
    to_plot = ['Energy','|m|','m2','mk2x']
    L = sys.argv[4:]
    Bex = sys.argv[1]
    folderName = sys.argv[2]
    mech = ['false','true']#sys.argv[3]

    # plot relaxation time as function of temperature
    simulations = analysis_tools.get_simulations(L, folderName, Bex, mech)
    plot_relax_times(simulations, to_plot)
    
    # plot relaxation function as function of t (MCS)
    #simulations = analysis_tools.get_simulations(L, folderName, Bex, mech, T=[1.61])
    #plot_relax_func(simulations, to_plot)
    
    


if __name__ == "__main__":
    main()