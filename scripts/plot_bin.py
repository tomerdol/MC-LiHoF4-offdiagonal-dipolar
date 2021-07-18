"""
Plot various quantities vs. T
Replaces plot.py which used results from a single MC
simulation instead of many independent runs.
"""
import glob
import config
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import analysis_tools, bin_data
from itertools import cycle


def main_plot(simulations, boot_num, plot_options, to_plot='', shift_T=False):
    """
    Plot the MC simulation results for the simulations given in the 'simulations' table.

    :param simulations: pandas DataFrame containing all simulations to plot
    :param boot_num: number of bootstrap samples
    :param plot_options: dict of plot options (axis scales etc.)
    :param to_plot: (optional) which quantity to plot vs. T. if not given, then either the
        Binder ratio of the correlation length / L are plotted, according to what is defined in
        plot_options.
    :param shift_T: whether to shift the T axis according to hardcoded values so that different
        projects have the same Tc
    :return: all plotted curves in a padnas DataFrame
    """
    Nsigma = 1.
    markers = ['o','s','^','D','v']

    all_y_curves = []

    # iterate over the given simulations. each row is a data point in the final figure
    for i, sim in enumerate(simulations.itertuples()):
        y = bin_data.read_binned_data(sim, use_latest=True)

        single_ydata=[]
        # each point is the result of bootstrap sampling
        for boot_index in range(boot_num):
            if to_plot=='':
                # this means a scaling function (Binder ratio or correlation length should
                # be plotted, according to what is defined in plot_options)
                single_ydata.append(plot_options['func'](y['Magnetization^2'].sample(frac=1,replace=True),y['Magnetization^4'].sample(frac=1,replace=True),y['mk2'+plot_options['corr_length_axis']].sample(frac=1,replace=True),sim.L*plot_options['unit_cell_length']))  # current y value
            else:
                # this means 'to_plot' should be plotted
                single_ydata.append(plot_options['func'](y[to_plot].sample(frac=1,replace=True)))

        all_yboot=np.array(single_ydata)

        # make pandas Series with the bootstrap mean and error,
        # and add it to the simulations table
        y_to_plot = pd.Series([np.mean(all_yboot,0), Nsigma * np.std(all_yboot,0)],index=['y_to_plot','y_to_plot_err'])
        y_to_plot=y_to_plot.append(pd.Series(sim._asdict()))
        all_y_curves.append(y_to_plot)

    all_y_curves = pd.DataFrame(all_y_curves)
    # shift different T axes so that they have the same Tc. not in use
    if shift_T:
        all_y_curves.loc[(all_y_curves['mech']=='false') & (all_y_curves['Bex']==0.0),'T'] -= 1.5824758958651635
        all_y_curves.loc[(all_y_curves['mech']=='true') & (all_y_curves['Bex']==0.0),'T'] -= 1.7650756636778697
        all_y_curves.loc[(all_y_curves['mech']=='false') & (all_y_curves['Bex']==0.3),'T'] -= 1.591939036664854
        all_y_curves.loc[(all_y_curves['mech']=='true') & (all_y_curves['Bex']==0.3),'T'] -= 1.7626781407441243

    fig, ax = plt.subplots(figsize=(8,6))
    
    for (label, df), marker in zip(all_y_curves.groupby(['Bex','L','folderName','mech']), cycle(markers)):
        df.plot(x='T',y='y_to_plot', yerr='y_to_plot_err', ax=ax, label=format_label(label), capsize=3, marker=marker)
    
    if shift_T:
        ax.set_xlabel('T (shifted by $T_c$)')
    else:
        ax.set_xlabel('T')

    ax.set_yscale(plot_options['axis_yscale'])
    plt.ylabel(plot_options['Name'])

    plt.legend(loc='best')
    plt.tight_layout()

    # save both a png and an eps file
    fig.savefig('../' + config.system_name + '/figures/plot_%s_%s_%s_%s_%s.png'%(to_plot,'_'.join(map(str,simulations['Bex'].unique().tolist())),'_'.join(map(str,simulations['mech'].unique().tolist())),'_'.join(map(str,simulations['L'].unique().tolist())),'_'.join(map(str,simulations['folderName'].unique().tolist()))), dpi=300)
    fig.savefig('../' + config.system_name + '/figures/plot_%s_%s_%s_%s_%s.eps'%(to_plot,'_'.join(map(str,simulations['Bex'].unique().tolist())),'_'.join(map(str,simulations['mech'].unique().tolist())),'_'.join(map(str,simulations['L'].unique().tolist())),'_'.join(map(str,simulations['folderName'].unique().tolist()))))
    plt.close()

    return all_y_curves


def format_label(label_list, format=['Bex','L','folderName','mech']):
    """
    Create a string format that indicates the parameters of a simulation.

    :param label_list: list of labels that describe the parameters of a simulation, e.g. [0.0, 5, project1, true]
    :param format: the order of the labels in label_list -- what does each item
        describe, e.g. ['Bex','L','folderName','mech']
    :return: a string of the formatted data
    """
    mech_name = label_list[format.index('mech')]
    mech_name = 'incl.' if mech_name=='false' else 'excl.'
    return_string = '$B_{x}=%s$, L=%s, name=%s, ODD=%s' \
                    %(label_list[format.index('Bex')], label_list[format.index('L')], label_list[format.index('folderName')], mech_name)
    if 'T' in format:
        return_string += ', T=%s' % label_list[format.index('T')]
    return return_string


def plot_lattice_correlators(simulations, plot_options, axes, to_plot='spinSize', shift_T=True):
    """
    Plot lattice correlators based on the final state of the simulations saved in data/lattice_output/.

    :param simulations: pandas DataFrame containing all simulations whose lattice correlators should be plotted
    :param plot_options: plot options, e.g. axis scale
    :param axes: list of axes along which the correlators should be calculated and plotted
    :param to_plot: which spin-specific quantity to plot
    :param shift_T: whether to shift the T axis according to hardcoded values so that different
        projects have the same Tc
    :return: the plotted curves
    """
    import scipy.stats
    Nsigma=1.
    markers=['o','s','^','D','v']

    all_y_curves = []
    fig, ax = plt.subplots(figsize=(10,8))
    for i, sim in enumerate(simulations.itertuples()):
        Lx=int(sim.L)
        Lz=int(sim.L)
        # read the data from the 'lattice_output' file
        path='../' + config.system_name + '/data/lattice_output/'+sim.folderName+'/table_'+str(sim.L)+'_'+str(sim.L)+'_'+str(sim.Bex)+'_'+str(sim.T)+'_'+str(sim.mech)+'_'+'*'+'.txt'
        file_list = glob.glob(path)
        data=[]
        for i, file in enumerate(file_list):
            data.append(analysis_tools.get_table_data_by_fname(file)[to_plot].to_numpy())
        data=np.array(data)
        # np.roll shifts the array of spins, based on their order in the array,
        # this corresponds to multiplying s[x=0]*s[x=1] etc.
        x_correlator=np.mean(data*np.roll(data,Lx*Lz*4,axis=1)) - np.mean(data)**2
        x_correlator_err=scipy.stats.sem(np.mean(data*np.roll(data,Lx*Lz*4,axis=1),axis=0))

        y_correlator=np.mean(data*np.roll(data,Lz*4,axis=1)) - np.mean(data)**2
        y_correlator_err=scipy.stats.sem(np.mean(data*np.roll(data,Lz*4,axis=1),axis=0))

        z_correlator=np.mean(data*np.roll(data,4,axis=1)) - np.mean(data)**2
        z_correlator_err=scipy.stats.sem(np.mean(data*np.roll(data,4,axis=1),axis=0))

        y_to_plot = pd.Series([x_correlator, Nsigma * x_correlator_err, y_correlator, y_correlator_err, z_correlator, z_correlator_err],
                              index=['x_correlator','x_correlator_err','y_correlator','y_correlator_err','z_correlator','z_correlator_err'])
        y_to_plot=y_to_plot.append(pd.Series(sim._asdict()))
        all_y_curves.append(y_to_plot)

    all_y_curves = pd.DataFrame(all_y_curves)
    # shift different T axes so that they have the same Tc
    if shift_T:
        all_y_curves.loc[(all_y_curves['mech']=='false') & (all_y_curves['Bex']==0.0),'T'] -= 1.5559215348091224
        all_y_curves.loc[(all_y_curves['mech']=='true') & (all_y_curves['Bex']==0.0),'T'] -= 1.7650756636778697
        all_y_curves.loc[(all_y_curves['mech']=='false') & (all_y_curves['Bex']==0.3),'T'] -= 1.547928848272695
        all_y_curves.loc[(all_y_curves['mech']=='true') & (all_y_curves['Bex']==0.3),'T'] -= 1.7626781407441243
        all_y_curves.loc[(all_y_curves['mech']=='false') & (all_y_curves['Bex']==0.6),'T'] -= 1.534166216592495
        all_y_curves.loc[(all_y_curves['mech']=='false') & (all_y_curves['Bex']==1.0),'T'] -= 1.5051312414585143

    # plot the results
    for (label, df), marker in zip(all_y_curves.groupby(['Bex','L','folderName','mech']), cycle(markers)):
        for axis in axes:
            df.plot(x='T',y=axis+'_correlator', yerr=axis+'_correlator_err', ax=ax, label=format_label(label) + ' | ' + axis + ' correlator', capsize=3, marker=marker)

    if shift_T:
        ax.set_xlabel('T (shifted by $T_c$)')
    else:
        ax.set_xlabel('T')

    ax.set_yscale(plot_options['axis_yscale'])
    plt.ylabel(plot_options['Name'])

    plt.legend(loc='best')
    plt.tight_layout()

    # save both a png and an eps file
    fig.savefig('../' + config.system_name + '/figures/plot_corr_%s_%s_%s_%s_%s.png'%(to_plot,'_'.join(map(str,simulations['Bex'].unique().tolist())),'_'.join(map(str,simulations['mech'].unique().tolist())),'_'.join(map(str,simulations['L'].unique().tolist())),'_'.join(map(str,simulations['folderName'].unique().tolist()))), dpi=300)
    fig.savefig('../' + config.system_name + '/figures/plot_corr_%s_%s_%s_%s_%s.eps'%(to_plot,'_'.join(map(str,simulations['Bex'].unique().tolist())),'_'.join(map(str,simulations['mech'].unique().tolist())),'_'.join(map(str,simulations['L'].unique().tolist())),'_'.join(map(str,simulations['folderName'].unique().tolist()))))
    plt.close()

    return all_y_curves


def parse_arguments():
    """Parse command line arguments. Uses the common parser from config.parse_arguments."""
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    
    parser = ArgumentParser(description="Analyzes Monte Carlo results and plots correlation length curves.", formatter_class=ArgumentDefaultsHelpFormatter, parents=[config.parse_arguments()], conflict_handler='resolve')
    parser.add_argument( "--to_plot", type=str, nargs='?', default='corr_length', help = "Which observable should be plotted. E.g. \"|Magnetization|\", \"spin_correlator\", \"binder\", etc. Default is Correlation length / L.")
    
    args = parser.parse_args()
    config.system_name = args.system_name
    return args


def remove_max_temperature(simulations):
    """Remove the highest temperature from the simulations table for each simulation group"""
    ret=[]
    for _, sim in simulations.groupby(['Bex','L','folderName','mech']):
        ret.append(sim.loc[sim['T']!=sim['T'].max()])
    return pd.concat(ret).reset_index(drop=True)


def main():
    """
    Plot the results of a MC simulation vs. T.

    :return: pandas DataFrame with the plotted curves
    """
    args = parse_arguments()
    L = args.L                      
    boot_num = args.boot_num           
    h_ex = args.h_ex
    mech = args.mech
    folderName = args.folder_list 
    to_plot = args.to_plot

    simulations = analysis_tools.get_simulations(L, folderName, h_ex, mech)

    # if we plot the acceptance rate, we do not need the highest T
    # since it is never swapped for a higher T
    if to_plot == 'swap':
        simulations = remove_max_temperature(simulations)

    from fit import get_binder, get_correlation_length

    # plot packages for each quantity
    if to_plot == 'corr_length':
        if config.system_name == 'Fe8':
            corr_length_axis='z'
        elif config.system_name == 'LiHoF4':
            corr_length_axis='x'
        plot_options = {'Name':r'$\xi^{(%s)}_{L} / L$'%corr_length_axis, 'axis_yscale':'log', 'func':get_correlation_length, 'corr_length_axis':corr_length_axis, 'unit_cell_length':1.0}
        # we can also plot the z correlation length for LiHoF4 by uncommenting the next line:
        # plot_options = {'Name':r'$\xi^{(%s)}_{L} / L$'%corr_length_axis, 'axis_yscale':'log', 'func':get_correlation_length, 'corr_length_axis':corr_length_axis, 'unit_cell_length':2.077294686}
        to_plot=''  # the scaling functions are the default of main_plot, so nothing need to be given in to_plot
    elif to_plot == 'binder':
        plot_options = {'Name':'g', 'axis_yscale':'linear', 'func':get_binder, 'corr_length_axis':'x','unit_cell_length':1.0}
        to_plot=''  # the scaling functions are the default of main_plot, so nothing need to be given in to_plot
    elif to_plot.split('_')[-1] == 'correlator':
        # if the to_plot given is something like localBx_correlator, we conclude that we should plot a correlator
        # of the quantity localBx
        plot_options = {'Name':to_plot, 'axis_yscale':'linear'}
        return plot_lattice_correlators(simulations, plot_options, ['x','y','z'], to_plot=to_plot.split('_')[:-1], shift_T=False)
    else:
        # just plot the given quantity. 'corr_length_axis':'x', 'unit_cell_length':1.0 are not used in this case.
        plot_options = {'Name':to_plot, 'axis_yscale':'linear', 'func':lambda x: np.mean(x), 'corr_length_axis':'x', 'unit_cell_length':1.0}

    data = main_plot(simulations, boot_num, plot_options, to_plot=to_plot)
    return data


if __name__ == "__main__":
    data=main()
