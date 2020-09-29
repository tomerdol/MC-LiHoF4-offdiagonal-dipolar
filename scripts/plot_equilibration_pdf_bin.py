import matplotlib
#matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import math
import csv
import analysis_tools
import bin_data

def check_equilibration(array, min_bin, consec_bins_for_equilib):
    #print(array)
    curr_bin=min_bin
    eq=False
    num_of_bins=array.shape[0] 
    #print(num_of_bins)
    while not eq and curr_bin<=num_of_bins-consec_bins_for_equilib:
        # cond1 is the standard condition, meaning all #consec_bins_for_equilib have overlapping errorbars
        # the condition itself checkes whether the maximum of the lower bound is below the minimum of the upper bounds
        cond = (np.max(array[curr_bin:curr_bin+consec_bins_for_equilib,0]-array[curr_bin:curr_bin+consec_bins_for_equilib,1]) 
        < np.min(array[curr_bin:curr_bin+consec_bins_for_equilib,0]+array[curr_bin:curr_bin+consec_bins_for_equilib,1]))
        
        if cond:
            eq=True
        curr_bin=curr_bin+1
        
    if eq:
        return curr_bin-1+consec_bins_for_equilib-1 # the index of the last bin in the equilibrated sequence
    else:
        return num_of_bins
     

def to_plot_col_index(to_plot):
    if to_plot=='index' or to_plot=='bin':
        return 0
    if to_plot=='Energy':
        return 1
    if to_plot=='|Magnetization|':
        return 2
    if to_plot=='Magnetization^2':
        return 3
    if to_plot=='mk2x':
        return 15

def save_equilibration_data(sim_name, equilibrated_bin):
    #sim_name[0]=Bex; sim_name[1]=L; sim_name[2]=folderName; sim_name[3]=mech
    path='../data/results/'+str(sim_name[2])+'/binned_data/equilib_data_'+str(sim_name[1])+'_'+str(sim_name[1])+'_'+str(sim_name[0])+'_'+str(sim_name[3])+'.txt'
    with open(path,'w') as outfile:
        outfile.write("%s\n" % equilibrated_bin)
        #for seed in seeds:
        #    outfile.write("%s\n" % seed)
    
def read_equilibration_data(sim_name):
    #sim_name[0]=Bex; sim_name[1]=L; sim_name[2]=folderName; sim_name[3]=mech
    path='../data/results/'+str(sim_name[2])+'/binned_data/equilib_data_'+str(sim_name[1])+'_'+str(sim_name[1])+'_'+str(sim_name[0])+'_'+str(sim_name[3])+'.txt'
    
    with open(path,'r') as file:
        eq_bin=file.readline()
    return int(eq_bin)
        
    

def main_check_equilibration(simulations, to_check):
    # a temperature ensemble is considered for equilibration (due to the use of parallel tempering)
    equilibrated_bin_dict={}
    for group_name, group_df in simulations.groupby(['Bex','L','folderName','mech']):
        if group_df['overwrite'].any():
            max_bin=0
            for i, sim in enumerate(group_df.itertuples()):
                data = bin_data.read_binned(sim, use_latest=False)
                for to_check_now in to_check:
                    a_index=data[0][:,0]
                    a=data[0][:,to_plot_col_index(to_check_now)]
                    a_err=data[1][:,to_plot_col_index(to_check_now)]
                    num_of_bins = len(a)
                    
                    equilibrated_bin = check_equilibration(np.stack([a,a_err],axis=-1),3,3)
                    if equilibrated_bin<num_of_bins:
                        max_bin = equilibrated_bin if equilibrated_bin>max_bin else max_bin
                    else:
                        max_bin = -1
            equilibrated_bin_dict[group_name]=max_bin+1
            save_equilibration_data(group_name, max_bin)
        else:
            equilibrated_bin_dict[group_name]=read_equilibration_data(group_name)
    # add column eq_bin with equilibrated bin
    simulations['eq_bin']=simulations.apply(lambda row: equilibrated_bin_dict[(row['Bex'],row['L'],row['folderName'],row['mech'])], axis=1)
    return simulations


def main_plot(simulations, to_plot, L, Bex, folderName, mech):
    markers=['o','s','^','D','v']	
    with PdfPages('../figures/plot_equilibration_%s_%s.pdf'%(Bex,mech)) as pdf:
        
        for temperature_index, temperature in enumerate(simulations['T'].unique()):
            fig, axes = plt.subplots(len(to_plot),1, figsize=(7,len(to_plot)*3))
            for i, sim in enumerate(simulations[simulations['T']==temperature].itertuples()):
                data=bin_data.read_binned(sim, use_latest=False)
                
                axes[0].set_title("Bex=%s , T=%1.5f"%(sim.Bex,temperature))
                for to_plot_index, to_plot_now in enumerate(to_plot):
                    a_index=data[0][:,0]
                    a=data[0][:,to_plot_col_index(to_plot_now)]
                    a_err=data[1][:,to_plot_col_index(to_plot_now)]
                    #a=a[12:]
                    #a_err=a_err[12:]
                    num_of_bins = len(a)
                    
                    
                    line = axes[to_plot_index].errorbar(a_index,a,yerr=a_err,fmt=markers[i % len(markers)]+'-',label=sim.L,capsize=2,fillstyle='none',mew=.7,linewidth=0.5)
                    equilibrated_bin = check_equilibration(np.stack([a,a_err],axis=-1),3,3)
                    axes[to_plot_index].annotate('',(equilibrated_bin,a[equilibrated_bin if equilibrated_bin<num_of_bins else equilibrated_bin-1]),xytext=(i,-4),textcoords='offset points',arrowprops=dict(arrowstyle='simple', color=line.lines[0].get_color()))
 
                    #axes[to_plot_index].annotate('',(equilibrated_bin,array[equilibrated_bin if equilibrated_bin<num_of_bins else equilibrated_bin-1,0]),xytext=(l_index,-4),textcoords='offset points',arrowprops=dict(arrowstyle='simple', color=line.lines[0].get_color()))
                    #if max_num_of_bins+1 > int(axes[to_plot_index].get_xticks(minor=False)[-1]):
                    #    axes[to_plot_index].set_xticks(np.arange(0,max_num_of_bins+1),minor=True)
                    
                    axes[to_plot_index].set_ylabel(to_plot_now)
                    axes[to_plot_index].set_xticks(range(math.floor(axes[to_plot_index].get_xlim()[0])+1,math.ceil(axes[to_plot_index].get_xlim()[1])))
            
            handles, labels = axes[0].get_legend_handles_labels()
            plt.figlegend(handles, labels, loc='upper right',prop={'size': 10})
            #fig.text(0.5, 0.04, r'MCS bin ($2^{i-1}$ to $2^{i}-1$)', ha='center')
            plt.xlabel(r'MCS bin ($2^{i}-1$ to $2(2^{i}-1)$)')
            plt.xticks(range(math.floor(plt.xlim()[0])+1,math.ceil(plt.xlim()[1])+1))
            
            try:
                fig.tight_layout()
            except:
                pass
            pdf.savefig()
            plt.close()

        #plt.title('%s vs. MCS bin (T=%s)'%(to_plot,T))
        #plt.xlabel(r'MCS bin ($2^{i-1}$ to $2^{i}-1$)')
        
        #plt.ylabel(to_plot)

        #plt.xscale('log')
        #plt.show()
        #fig.savefig('../figures/plot_equilibration_%s.pdf'%Bex)

def main():
    to_plot = ['Energy','|Magnetization|','Magnetization^2','mk2x']
    L = sys.argv[4:]
    Bex = sys.argv[1]
    folderName = sys.argv[2]
    mech = sys.argv[3]

    simulations = analysis_tools.get_simulations(L, folderName, Bex, mech)
    main_plot(simulations, to_plot, L, Bex, folderName, mech)
    print(main_check_equilibration(simulations, to_plot))


if __name__ == "__main__":
    main()
