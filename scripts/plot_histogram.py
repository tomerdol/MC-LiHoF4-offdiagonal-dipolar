import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import glob
import analysis_tools
import pandas as pd

def parse_arguments():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

    parser = ArgumentParser(description="Plots histograms from lattice output files for LiHoF4", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument( "-L", nargs='+', type=int, required=True, help = "Linear system sizes.")
    parser.add_argument( "--h_ex", type=float, help = "External magnetic field value, Bex." , required=True)
    parser.add_argument( "-m", "--mech", nargs='+', choices=['true','false'], help = ("Whether internal fields are suppressed or not. \'false\' means "
                                                                                      "that they aren't so the mechanism is on, and \'true\' means that they are and the mechanism is off." ), required=True)
    parser.add_argument( "-f", "--folder_list", nargs='+', type=str, help = "List of folders in \'data/results/\' in which results should be found. " , required=True)
    parser.add_argument( "--to_plot", type=str, nargs='?', default='localBx', help = "Which observable should be plotted. Default is \'localBx\'")
    parser.add_argument( "-T", nargs='+', type=float, required=True, help = "Temperature. If not exact match, closest available temperature(s) will be used.")
    args = parser.parse_args()

    return args


def main_hist(simulations, to_plot):

    all_simulations=[]
    for i, sim in enumerate(simulations.itertuples()):
        path='../data/lattice_output/'+sim.folderName+'/table_'+str(sim.L)+'_'+str(sim.L)+'_'+str(sim.Bex)+'_'+str(sim.T)+'_'+str(sim.mech)+'_'+'*'+'.txt'

        file_list = glob.glob(path)
        data = pd.concat((analysis_tools.get_table_data_by_fname(f) for f in file_list))
        data['Simulation']=sim.Index
        all_simulations.append(data)
    all_simulations = pd.concat(all_simulations, axis=0, ignore_index=True)

    fig, ax = plt.subplots(figsize=(10,10))
    histograms_to_plot=[]
    histogram_labels=[]
    for (label, df) in all_simulations.groupby('Simulation'):
        histograms_to_plot.append(df[to_plot].to_numpy())
        histogram_labels.append(str(simulations.loc[int(label)].values.tolist()))
    plt.hist(histograms_to_plot, bins=30, alpha=0.5, label=histogram_labels)
    # plot = all_simulations.hist(column=to_plot, by='Simulation', grid=True, sharex=True)
    plt.legend()
    fig.savefig('../figures/foo.png')


def main():
    args = parse_arguments()
    L = args.L
    h_ex = args.h_ex
    mech = args.mech
    folderName = args.folder_list
    to_plot = args.to_plot
    T = args.T

    simulations = analysis_tools.get_simulations(L, folderName, h_ex, mech, T=T)
    main_hist(simulations, to_plot)


if __name__ == "__main__":
    main()