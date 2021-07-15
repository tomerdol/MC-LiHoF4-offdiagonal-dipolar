"""
Plot histograms of site-specific quantities read from final simulation states in /data/lattice_output/
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
import analysis_tools
from scipy.ndimage import gaussian_filter1d
import config


def parse_arguments():
    """Parse command line arguments. Uses the common parser from config.parse_arguments."""
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

    parser = ArgumentParser(description="Plots histograms from lattice output files for LiHoF4", formatter_class=ArgumentDefaultsHelpFormatter, parents=[config.parse_arguments()], conflict_handler='resolve')
    parser.add_argument( "-L", nargs='+', type=int, required=True, help = "Linear system sizes.")
    parser.add_argument( "--to_plot", type=str, nargs='+', default='localBx', help = "Which observable should be plotted. Default is \'localBx\'")
    parser.add_argument( "-T", nargs='+', type=float, required=True, help = "Temperature. If not exact match, closest available temperature(s) will be used.")
    parser.add_argument("--flip", action="store_true", help="Also plot the flipped distribution. Useful for visualizing any asymmetry.")
    parser.add_argument("--seed", type=str, default='*', help="Specify the seed of the simulation whose histogram is to be plotted.")
    args = parser.parse_args()
    config.system_name = args.system_name
    return args


def main_hist(simulations, to_plot, flip=False, seed='*'):
    """
    Plot histograms of site-specific quantities for the given simulations

    :param simulations: pandas DataFrame containing all simulations for which histograms should be plotted
    :param to_plot: which quantity(ies) to plot
    :param flip: whether to also plot the flipped histogram for each of the simulations (to easily
                    see if they are symmetric)
    :param seed: (optional) a specific seed for which to plot a histogram
    :return: None
    """
    from fit import str_with_err
    from plot_bin import format_label
    fig, ax = plt.subplots(figsize=(10,10))
    prop_iter = iter(plt.rcParams['axes.prop_cycle'])

    for to_plot_now in to_plot:
        num_of_simulations_to_plot=len(simulations)

        multiple_sim_data=[]
        multiple_sim_groups=[]
        for i, sim in enumerate(simulations.itertuples()):
            path='../' + config.system_name + '/data/lattice_output/'+sim.folderName+'/table_'+str(sim.L)+'_'+str(sim.L)+'_'+str(sim.Bex)+'_'+str(sim.T)+'_'+str(sim.mech)+'_'+str(seed)+'.txt'

            file_list = glob.glob(path)
            data=[]
            groups=[]
            num_independent_runs=len(file_list)
            for j, file in enumerate(file_list):
                shift = 0 if to_plot_now != 'localBx' else float(sim.Bex)
                temp_data=analysis_tools.get_table_data_by_fname(file)[to_plot_now].to_numpy() - shift
                data.append(temp_data)
                groups.append(np.full(temp_data.size,j))
            data=np.concatenate(data).ravel()
            groups=np.concatenate(groups).ravel()
            multiple_sim_data.append(data)
            multiple_sim_groups.append(groups)
            if flip:
                multiple_sim_data.append(-data)
                multiple_sim_groups.append(groups)
        shared_bins = np.histogram_bin_edges(np.concatenate(multiple_sim_data).ravel(), bins=120)

        for i, sim in enumerate(simulations.itertuples()):
            sim_index=i if not flip else i*2
            data=multiple_sim_data[sim_index]
            groups=multiple_sim_groups[sim_index]
            shared_bins_to_plot = shared_bins[:-1] + i*np.diff(shared_bins)/num_of_simulations_to_plot
            num_independent_runs=len(np.unique(groups))
            histograms = np.zeros((num_independent_runs, len(shared_bins)-1))
            #smooth_histogram_x = np.linspace(shared_bins[0],shared_bins[-1],num=200)
            #smooth_histogram = np.zeros((num_independent_runs, len(smooth_histogram_x)))
            simple_values=np.zeros(num_independent_runs)
            abs_values=np.zeros(num_independent_runs)
            standard_deviations=np.zeros(num_independent_runs)
            #window_size=5
            for j in range(num_independent_runs):
                simple_values[j] = np.mean(data[groups==j])
                abs_values[j] = np.mean(np.abs(data[groups==j]))
                standard_deviations[j] = np.sqrt(np.mean(data[groups==j]**2))
                histograms[j], _ = np.histogram(data[groups == j], bins=shared_bins)

                #for idx, x in enumerate(smooth_histogram_x):
                #    bin_width = shared_bins[1]-shared_bins[0]
                #    low_bound = x - bin_width*(window_size*0.5)
                #    upper_bound = x + bin_width*(window_size*0.5)
                #    smooth_histogram[i,idx] = np.sum((groups==i) & (data > low_bound) & (data<upper_bound))/(window_size)
            smooth_histogram_gauss = gaussian_filter1d(histograms,2,axis=1)
            simple_values_mean = simple_values.mean()
            simple_values_err = simple_values.std()/np.sqrt(simple_values.size-1)
            abs_values_mean = abs_values.mean()
            abs_values_err = abs_values.std()/np.sqrt(abs_values.size-1)
            sd_mean = standard_deviations.mean()
            sd_err = standard_deviations.std()/np.sqrt(standard_deviations.size-1)

            print(format_label(list(sim)[2:],format=['L','folderName','Bex','mech','T']))
            print("-"*len(format_label(list(sim)[2:],format=['L','folderName','Bex','mech','T'])))
            print("MEAN: " + str(simple_values_mean) + ' +- ' + str(simple_values_err))
            print("ABS: " + str(abs_values_mean) + ' +- ' + str(abs_values_err))
            print("SD: " + str(sd_mean) + ' +- ' + str(sd_err))
            plt.bar(shared_bins_to_plot, np.mean(histograms,axis=0), alpha=0.9, align='edge', width=np.diff(shared_bins)/num_of_simulations_to_plot, label="%s, MEAN=%s, ABS=%s, SD=%s" % (format_label(list(sim)[2:],format=['L','folderName','Bex','mech','T']), str_with_err(simple_values_mean,simple_values_err), str_with_err(abs_values_mean,abs_values_err), str_with_err(sd_mean,sd_err)),
                    yerr=np.std(histograms,axis=0)/np.sqrt(num_independent_runs-1), color=next(prop_iter)['color'])


            #plt.plot(smooth_histogram_x, np.mean(smooth_histogram,axis=0),label="Smooth %s 2"%format_label(list(sim)[2:],format=['L','folderName','Bex','mech','T']), color=next(prop_iter)['color'])
            plt.plot(shared_bins_to_plot + 0.5 * np.diff(shared_bins)/num_of_simulations_to_plot, np.mean(smooth_histogram_gauss, axis=0), label="Gaussian Smooth %s 2"%format_label(list(sim)[2:],format=['L','folderName','Bex','mech','T']), color=next(prop_iter)['color'])
            if flip:
                data=multiple_sim_data[sim_index+1]
                groups=multiple_sim_groups[sim_index+1]
                smooth_histogram_x = np.linspace(shared_bins[0],shared_bins[-1],num=200)
                smooth_histogram = np.zeros((num_independent_runs, len(smooth_histogram_x)))
                for i in range(num_independent_runs):
                    histograms[i], _ = np.histogram(data[groups == i], bins=shared_bins)
                    for idx, x in enumerate(smooth_histogram_x):
                        bin_width = shared_bins[1]-shared_bins[0]
                        low_bound = x - bin_width*(window_size*0.5)
                        upper_bound = x + bin_width*(window_size*0.5)
                        smooth_histogram[i,idx] = np.sum((groups==i) & (data > low_bound) & (data<upper_bound))/(window_size)
                plt.bar(shared_bins_to_plot, np.mean(histograms,axis=0), align='edge', width=np.diff(shared_bins)/num_of_simulations_to_plot, label='Flipped %s' % format_label(list(sim)[2:],format=['L','folderName','Bex','mech','T']), alpha=0.5,
                        yerr=np.std(histograms,axis=0)/np.sqrt(num_independent_runs-1), ecolor='b', color=next(prop_iter)['color'])
                plt.plot(smooth_histogram_x, np.mean(smooth_histogram,axis=0),label="Flipped smooth %s"%format_label(list(sim)[2:],format=['L','folderName','Bex','mech','T']), color=next(prop_iter)['color'])

    plt.legend()
    plt.grid()
    plt.title('Distribution of ' + to_plot_now)
    ax.set_ylabel('# of spins')
    ax.set_xlabel(to_plot_now + ('' if shift == 0 else ' - ' + str(shift)))
    plt.tight_layout()
    fig.savefig('../' + config.system_name + '/figures/hist_%s_%s_%s_%s_%s.png'%('_'.join(map(str,simulations['Bex'].unique().tolist())),'_'.join(map(str,simulations['mech'].unique().tolist())),'_'.join(map(str,simulations['L'].unique().tolist())),'_'.join(map(str,simulations['folderName'].unique().tolist())),to_plot_now), dpi=300)
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
    seed=args.seed

    simulations = analysis_tools.get_simulations(L, folderName, h_ex, mech, T=T)
    main_hist(simulations, to_plot, flip=flip, seed=seed)


if __name__ == "__main__":
    main()
