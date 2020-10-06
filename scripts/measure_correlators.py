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
    parser.add_argument( "--axis", type=str, nargs='+', choices=['x', 'y', 'z'], help = "Axis along which the correlator is computed. Options are \'x\', \'y\' and \'z\', and more than one option can be given.")
    args = parser.parse_args()

    return args


def main_corr_calc(simulations, axis):
    Nsigma=1.
    markers=['o','s','^','D','v']

    all_y_curves = []

    for i, sim in enumerate(simulations.itertuples()):
        y = bin_data.read_binned_data(sim, use_latest=False, use_bin=sim.eq_bin)

        single_ydata=[]
        for boot_index in range(boot_num):
            if to_plot=='':
                # this means a scaling function (Binder ratio or correlation length should be plotted, according to what is defined in plot_options)
                single_ydata.append(plot_options['func'](y['Magnetization^2'].sample(frac=1,replace=True),y['Magnetization^4'].sample(frac=1,replace=True),y['mk2'+plot_options['corr_length_axis']].sample(frac=1,replace=True),sim.L*plot_options['unit_cell_length']))  # current y value
            else:
                # this means 'to_plot' should be plotted
                single_ydata.append(plot_options['func'](y[to_plot].sample(frac=1,replace=True)))
            #print(boot_index)

        all_yboot=np.array(single_ydata)

        y_to_plot = pd.Series([np.mean(all_yboot,0), Nsigma * np.std(all_yboot,0)],index=['y_to_plot','y_to_plot_err'])
        y_to_plot=y_to_plot.append(pd.Series(sim._asdict()))
        #print('Curve to be plotted (L=%s): '%sim.L)
        #print(y_to_plot)
        all_y_curves.append(y_to_plot)
        #all_y_curves_err.append(y_err_to_plot)

    all_y_curves = pd.DataFrame(all_y_curves)

    fig, ax = plt.subplots(figsize=(8,6))

    for (label, df), marker in zip(all_y_curves.groupby(['Bex','L','folderName','mech']), cycle(markers)):
        df.plot(x='T',y='y_to_plot', yerr='y_to_plot_err', ax=ax, label=str(label), capsize=3, marker=marker)

    ax.set_xlabel('T')

    ax.set_yscale(plot_options['axis_yscale'])
    plt.ylabel(plot_options['Name'])

    plt.legend(loc='best')
    plt.tight_layout()
    axes = plt.gca()

    fig.savefig('../figures/plot_%s_%s_%s_%s_%s.png'%(to_plot,'_'.join(map(str,simulations['Bex'].unique().tolist())),'_'.join(map(str,simulations['mech'].unique().tolist())),'_'.join(map(str,simulations['L'].unique().tolist())),'_'.join(map(str,simulations['folderName'].unique().tolist()))), dpi=300)
    fig.savefig('../figures/plot_%s_%s_%s_%s_%s.eps'%(to_plot,'_'.join(map(str,simulations['Bex'].unique().tolist())),'_'.join(map(str,simulations['mech'].unique().tolist())),'_'.join(map(str,simulations['L'].unique().tolist())),'_'.join(map(str,simulations['folderName'].unique().tolist()))))
    plt.close()


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

    simulations = analysis_tools.get_simulations(L, folderName, h_ex, mech)
    main_corr_calc(simulations, to_plot, flip=flip)


if __name__ == "__main__":
    main()