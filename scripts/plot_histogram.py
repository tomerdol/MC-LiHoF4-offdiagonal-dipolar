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
    parser.add_argument( "--to_plot", type=str, nargs='+', default='localBx', help = "Which observable should be plotted. Default is \'localBx\'")
    parser.add_argument( "-T", nargs='+', type=float, required=True, help = "Temperature. If not exact match, closest available temperature(s) will be used.")
    parser.add_argument("--flip", action="store_true", help="Also plot the flipped distribution. Useful for visualizing any asymmetry.")
    args = parser.parse_args()

    return args


def main_hist2(simulations, to_plot, flip=False):
    fig, ax = plt.subplots(figsize=(10,10))
    for to_plot_now in to_plot:
        for i, sim in enumerate(simulations.itertuples()):
            path='../data/lattice_output/'+sim.folderName+'/table_'+str(sim.L)+'_'+str(sim.L)+'_'+str(sim.Bex)+'_'+str(sim.T)+'_'+str(sim.mech)+'_'+'*'+'.txt'

            file_list = glob.glob(path)
            data=[]
            groups=[]
            num_independent_runs=len(file_list)
            for i, file in enumerate(file_list):
                shift = 0 if to_plot_now != 'localBx' else float(simulations['Bex'].iloc[0])
                temp_data=analysis_tools.get_table_data_by_fname(file)[to_plot_now].to_numpy() - shift
                data.append(temp_data)
                groups.append(np.full(temp_data.size,i))
            data=np.array(data)
            groups=np.array(groups)

            shared_bins = np.histogram_bin_edges(data, bins='auto')
            histograms = np.zeros((num_independent_runs, len(shared_bins)-1))
            for i in range(num_independent_runs):
                histograms[i], _ = np.histogram(data[groups == i], bins=shared_bins)

            np.mean(histograms,axis=0)
            plt.hist(np.mean(histograms,axis=0), bins=shared_bins, label=str(sim.tolist()))

            if flip:
                plt.hist(-np.mean(histograms,axis=0), bins=shared_bins, label="Flipped " + str(sim.tolist()), alpha=0.5)

    plt.legend()
    plt.title('Distribution of ' + to_plot_now)
    ax.set_ylabel('# of spins')
    ax.set_xlabel(to_plot_now + ('' if shift == 0 else ' - ' + str(shift)))
    fig.savefig('../figures/hist_%s_%s_%s_%s_%s.png'%('_'.join(map(str,simulations['Bex'].unique().tolist())),'_'.join(map(str,simulations['mech'].unique().tolist())),'_'.join(map(str,simulations['L'].unique().tolist())),'_'.join(map(str,simulations['folderName'].unique().tolist())),to_plot_now), dpi=300)
    plt.close(fig)


def main_hist(simulations, to_plot, flip=False):

    all_simulations=[]
    for i, sim in enumerate(simulations.itertuples()):
        path='../data/lattice_output/'+sim.folderName+'/table_'+str(sim.L)+'_'+str(sim.L)+'_'+str(sim.Bex)+'_'+str(sim.T)+'_'+str(sim.mech)+'_'+'*'+'.txt'

        file_list = glob.glob(path)
        data = pd.concat((analysis_tools.get_table_data_by_fname(f) for f in file_list))
        data['Simulation']=sim.Index
        all_simulations.append(data)
    all_simulations = pd.concat(all_simulations, axis=0, ignore_index=True)

    for to_plot_now in to_plot:
        fig, ax = plt.subplots(figsize=(10,10))
        histograms_to_plot=[]
        histogram_labels=[]
        shift = 0 if to_plot_now != 'localBx' else float(simulations['Bex'].iloc[0])
        for (label, df) in all_simulations.groupby('Simulation'):
            histograms_to_plot.append(df[to_plot_now].to_numpy() - shift)
            histogram_labels.append(str(simulations.loc[int(label)].values.tolist()))

        xmin = np.array(histograms_to_plot).min()
        xmax = np.array(histograms_to_plot).max()
        plt.hist(histograms_to_plot, bins=60, label=histogram_labels, range=(xmin,xmax))
        if flip:
            plt.hist([-hist for hist in histograms_to_plot], bins=60, label=["Flipped " + hist_label for hist_label in histogram_labels], alpha=0.5, range=(-xmax,-xmin))
        # plot = all_simulations.hist(column=to_plot, by='Simulation', grid=True, sharex=True)
        plt.legend()
        plt.title('Distribution of ' + to_plot_now)
        ax.set_ylabel('# of spins')
        ax.set_xlabel(to_plot_now + ('' if shift == 0 else ' - ' + str(shift)))
        fig.savefig('../figures/hist_%s_%s_%s_%s_%s.png'%('_'.join(map(str,simulations['Bex'].unique().tolist())),'_'.join(map(str,simulations['mech'].unique().tolist())),'_'.join(map(str,simulations['L'].unique().tolist())),'_'.join(map(str,simulations['folderName'].unique().tolist())),to_plot_now), dpi=300)
        plt.close(fig)


def main():
    from analysis_tools import listify
    args = parse_arguments()
    L = args.L
    h_ex = args.h_ex
    mech = args.mech
    folderName = args.folder_list
    to_plot = args.to_plot
    T = args.T
    to_plot=listify(to_plot)
    flip=args.flip

    simulations = analysis_tools.get_simulations(L, folderName, h_ex, mech, T=T)
    main_hist2(simulations, to_plot, flip=flip)


if __name__ == "__main__":
    main()