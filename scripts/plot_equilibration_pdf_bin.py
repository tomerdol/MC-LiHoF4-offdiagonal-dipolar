"""
Plot and test the equilibration of a MC simulation
"""
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import sys
import math
import analysis_tools
import bin_data
import config


def check_equilibration(array, min_bin, consec_bins_for_equilib):
    """
    Check whether the given array which represents the binned time series of an observable from a MC simulation
    is equilibrated, i.e. have x consecutive bins that agree within error bars after min_bin.
    :param array: logarithmically binned time series of an observable from a MC simulation --
    first column is the means and the second column is the errors.
    :param min_bin: minimum bin in which to seek equilibration
    :param consec_bins_for_equilib: x -- the number of consecutive bins required for equilibration
    :return: the index of the last bin in the equilibrated sequence. if the number of bins is returned (last index+1),
    it is an indication that the series has not equilibrated.
    """
    # where to start looking
    curr_bin=min_bin
    eq=False
    num_of_bins=array.shape[0] 
    while not eq and curr_bin<=num_of_bins-consec_bins_for_equilib:
        # starting from the current bin, whether the next #consec_bins_for_equilib have overlapping error bars.
        # the condition itself checks whether the maximum of the lower bound is below the minimum of the upper bounds
        eq = (np.max(array[curr_bin:curr_bin+consec_bins_for_equilib,0]-array[curr_bin:curr_bin+consec_bins_for_equilib,1])
        < np.min(array[curr_bin:curr_bin+consec_bins_for_equilib,0]+array[curr_bin:curr_bin+consec_bins_for_equilib,1]))

        curr_bin=curr_bin+1
        
    if eq:
        return curr_bin-1+consec_bins_for_equilib-1 # the index of the last bin in the equilibrated sequence
    else:
        return num_of_bins
     

def to_plot_col_index(to_plot):
    """
    Translate column name to column index
    :param to_plot: column name
    :return: corresponding column index
    """
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
    """
    Save the equilibrated bin of the simulation in a text file in /binned_data/ to use later
    :param sim_name: list of the simulation parameters:
    sim_name[0]=Bex; sim_name[1]=L; sim_name[2]=folderName; sim_name[3]=mech
    :param equilibrated_bin: the equilibrated bin to save
    :return: None
    """
    #sim_name[0]=Bex; sim_name[1]=L; sim_name[2]=folderName; sim_name[3]=mech
    path='../' + config.system_name + '/data/results/'+str(sim_name[2])+'/binned_data/equilib_data_'+str(sim_name[1])+'_'+str(sim_name[1])+'_'+str(sim_name[0])+'_'+str(sim_name[3])+'.txt'
    with open(path,'w') as outfile:
        outfile.write("%s\n" % equilibrated_bin)

    
def read_equilibration_data(sim_name):
    """
    Read the equilibrated bin of the simulation from a text file in /binned_data/
    :param sim_name: list of the simulation parameters:
    sim_name[0]=Bex; sim_name[1]=L; sim_name[2]=folderName; sim_name[3]=mech
    :return: the index of the equilibrated bin as read from file
    """
    #sim_name[0]=Bex; sim_name[1]=L; sim_name[2]=folderName; sim_name[3]=mech
    path='../' + config.system_name + '/data/results/'+str(sim_name[2])+'/binned_data/equilib_data_'+str(sim_name[1])+'_'+str(sim_name[1])+'_'+str(sim_name[0])+'_'+str(sim_name[3])+'.txt'
    
    with open(path,'r') as file:
        eq_bin=file.readline()
    return int(eq_bin)


def main_check_equilibration(simulations, to_check):
    """
    Checks the equilibration of the given simulations and adds the equilibrated bin to a
    new 'eq_bin' column in the given DataFrame
    :param simulations: pandas DataFrame of all simulations to check for equilibration
    :param to_check: a list of observables to check for equilibration
    :return: the given simulations table with an additional column for the equilibrated bin
    and with simulations that could not be read removed
    """
    # a temperature ensemble is considered for equilibration (due to the use of parallel tempering)
    equilibrated_bin_dict={}
    simulations_to_remove=[]
    for group_name, group_df in simulations.groupby(['Bex','L','folderName','mech']):
        if group_df['overwrite'].any():
            max_bin=0
            for i, sim in enumerate(group_df.itertuples()):
                try:
                    data = bin_data.read_binned(sim, use_latest=True)
                except Exception as e:
                    print(e)
                    print("Dropping missing simulations from current analysis.")
                    simulations_to_remove.append(sim.Index)
                    continue
                # check for each observable whether it is equilibrated
                for to_check_now in to_check:
                    a=data[0][:,to_plot_col_index(to_check_now)]
                    a_err=data[1][:,to_plot_col_index(to_check_now)]
                    num_of_bins = len(a)
                    
                    equilibrated_bin = check_equilibration(np.stack([a,a_err],axis=-1),3,3)
                    if equilibrated_bin<num_of_bins:
                        # the equilibrated bin for a simulation group is the maximum
                        # among all temperatures and all observables checked
                        max_bin = equilibrated_bin if equilibrated_bin>max_bin else max_bin
                    else:
                        # not equilibrated
                        max_bin = -1
            equilibrated_bin_dict[group_name]=max_bin
            save_equilibration_data(group_name, max_bin)
        else:
            equilibrated_bin_dict[group_name]=read_equilibration_data(group_name)
    simulations.drop(simulations_to_remove,inplace=True)
    # add column eq_bin with equilibrated bin
    simulations['eq_bin']=simulations.apply(lambda row: equilibrated_bin_dict[(row['Bex'],row['L'],row['folderName'],row['mech'])], axis=1)
    print(simulations)
    return simulations


def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique))


def main_plot(simulations, to_plot):
    """
    Create a PDF with the equilibration plots of all observables and all temperature of the given simulations.
    The given simulations should all have the same mech, folderName and Bex, and should only differ by their L.
    :param simulations: pandas DataFrame containing all simulations whose equilibration process should be plotted
    :param to_plot: list of observables whose equilibration is to be plotted
    :return: None
    """
    markers = ['o','s','^','D','v']
    colors = ("red", "green", "blue", "yellow", "orange")
    # each different L has its own curve
    groups = simulations['L'].unique().tolist()
    fig, ax = plt.subplots()
    with PdfPages('../' + config.system_name + '/figures/plot_equilibration_%s_%s_%s_%s.pdf' %
                  ('_'.join(map(str,simulations['Bex'].unique().tolist())),
                   '_'.join(map(str,simulations['mech'].unique().tolist())),
                   '_'.join(map(str,simulations['L'].unique().tolist())),
                   '_'.join(map(str,simulations['folderName'].unique().tolist())))) as pdf:

        # each page shows a different T
        for temperature_index, temperature in enumerate(simulations['T'].unique()):
            pdf_fig, pdf_axes = plt.subplots(len(to_plot),1, figsize=(7,len(to_plot)*3))
            # find all simulations (possibly different L's) that have the currently plotted T
            for i, sim in enumerate(simulations[simulations['T']==temperature].itertuples()):
                try:
                    data = bin_data.read_binned(sim, use_latest=True)
                except Exception as e:
                    print(e)
                    print("Dropping missing simulations from current analysis.")
                    continue
                pdf_axes[0].set_title("T=%1.5f" % temperature)
                # iterate over the observables to plot
                for to_plot_index, to_plot_now in enumerate(to_plot):
                    a_index=data[0][:,0]
                    a=data[0][:,to_plot_col_index(to_plot_now)]
                    a_err=data[1][:,to_plot_col_index(to_plot_now)]
                    num_of_bins = len(a)
                    # plot the equilibration curve for this L
                    line = pdf_axes[to_plot_index].errorbar(a_index,a,yerr=a_err,fmt=markers[i % len(markers)]+'-',
                                                            color=colors[groups.index(sim.L) % len(colors)],label=sim.L,
                                                            capsize=2,fillstyle='none',mew=.7,linewidth=0.5)
                    # find the equilibrated bin for this curve
                    equilibrated_bin = check_equilibration(np.stack([a,a_err],axis=-1),3,3)
                    # if it is not equilibrated, print out its details
                    if equilibrated_bin >= num_of_bins:
                        print("Simulation " + str(sim) + " : " + str(to_plot_now) + " not equilibrated.")

                    # add this equilibrated bin to a different figure showing all equilibrated bins (ax)
                    ax.scatter(temperature, equilibrated_bin, label='L='+str(sim.L), c=colors[groups.index(sim.L) % len(colors)])
                    # add a small tick to indicate the equilibrated bin (offset its angle slightly so that multiple
                    # such ticks will not cover each other).
                    pdf_axes[to_plot_index].annotate('',(equilibrated_bin,a[equilibrated_bin if equilibrated_bin<num_of_bins else equilibrated_bin-1]),xytext=(i,-4),textcoords='offset points',arrowprops=dict(arrowstyle='simple', color=line.lines[0].get_color()))

                    pdf_axes[to_plot_index].set_ylabel(to_plot_now)
                    pdf_axes[to_plot_index].set_xticks(range(math.floor(pdf_axes[to_plot_index].get_xlim()[0])+1,math.ceil(pdf_axes[to_plot_index].get_xlim()[1]+1)))
            
            handles, labels = pdf_axes[0].get_legend_handles_labels()
            plt.figlegend(handles, labels, loc='upper right',prop={'size': 10})
            plt.xlabel(r'MCS bin ($2^{i}-1$ to $2(2^{i}-1)$)')
            plt.xticks(range(math.floor(plt.xlim()[0])+1,math.ceil(plt.xlim()[1])+1))
            
            try:
                pdf_fig.tight_layout()
            except:
                # not terrible
                pass

            pdf.savefig()
            plt.close()

    # a separate figure of all equilibrated bins
    legend_without_duplicate_labels(ax)
    ax.set_xlabel('T')
    ax.set_ylabel('Equilibrated bin')
    # add a horizontal line showing how many bins there are overall
    ax.axhline(num_of_bins-1)
    fig.savefig('../' + config.system_name + '/figures/plot_equilibration_%s_%s_%s_%s.png' %
                ('_'.join(map(str,simulations['Bex'].unique().tolist())),
                 '_'.join(map(str,simulations['mech'].unique().tolist())),
                 '_'.join(map(str,simulations['L'].unique().tolist())),
                 '_'.join(map(str,simulations['folderName'].unique().tolist()))))
    print(simulations)


def main():
    to_plot = ['Energy','|Magnetization|','Magnetization^2','mk2x']
    L = sys.argv[5:]
    Bex = sys.argv[2]
    config.system_name = sys.argv[1]
    folderName = sys.argv[3]
    mech = sys.argv[4]

    simulations = analysis_tools.get_simulations(L, folderName, Bex, mech)
    main_plot(simulations, to_plot)


if __name__ == "__main__":
    main()
