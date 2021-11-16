"""
Performs a least-squares fit of multiple curves (for multiple linear system sizes)
to find the critical temperature and critical exponents.
Method names end with _bin since they work on data read from the binned_data directory.
A previous version of these methods that used the raw data from a single long MC simulation
can be found in revision 4009f266ea7f5fd18567f334b1181ece2c356947
"""
import matplotlib
matplotlib.use('Agg')
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
import sys
import pandas as pd
import math
import warnings
import analysis_tools
import bin_data
from itertools import cycle
import config


def f_bin(par, data):
    """
    Third-order polynomial function used to represent the scaling function for plotting and fitting:

    y(x) = p0 + p1*(x-x_c)*L^(1/v) + p2*(x-x_c)^2*L^(2/v) + p3*(x-x_c)^3*L^(3/v).

    :param par: 6-tuple parameters of the polynomial: (x_c, p0, p1, p2, p3, v)
    :param data: NumPy array of data:   first column is the linear system size of the point, L, and the
                                        second column is the x (T) value of the point
    :return: NumPy array of the values of the polynomial for each of the points defined by the received data
    """
    x_c, p0, p1, p2, p3, v = par
    L=data[:,0]
    x=data[:,1]
    try:
        return p0 + p1*((L**(1/v))*(x-x_c)) + p2*((L**(1/v))*(x-x_c))**2 + p3*((L**(1/v))*(x-x_c))**3
    except Exception as e:
        print(e)
        print(par, file=sys.stderr)


def fit_multiple_bin(data, initial_xc, bounds=False, input_p_global=[ -0.4, 0.2, -0.06, 0.6 ], max_nfev=None):
    """
    Fit curves for multiple system sizes that should all be described by the same universal function when rescaled.

    :param data: data to be fitted: formatted in a NumPy array of three columns: 0=L, 1=T, 2=y
    :param initial_xc: initial guess for the critical temperature (crossing of the different curves)
    :param bounds: whether to use bounds in the fitting process
    :param input_p_global: initial values for the fitting parameters
    :param max_nfev: maximum number of calls to the fitting function
    :return: best fit parameters
    """
    # format of data columns: 0-L, 1-T, 2-y
    xdata=data[:,0:2]
    ydata=data[:,2]
    # get min and max values that will be used to estimate the parameter p0
    max_val = ydata.max()
    min_val = ydata.min()

    p_global = [initial_xc, 0.5*(max_val+min_val)] + input_p_global

    if not bounds:
        p_best=optimize.least_squares(lambda p, x, y: f_bin(p, x) - y, p_global, args=(xdata,ydata), max_nfev=max_nfev, method='lm')
    else:
        p_best=optimize.least_squares(lambda p, x, y: f_bin(p, x) - y, p_global, args=(xdata,ydata), bounds=([0,min_val,-np.inf,-np.inf,-np.inf,0.0],[3,max_val,np.inf,np.inf,np.inf,np.inf]))

    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    # R^2 < 0.98
    if 1-p_best.cost/ss_tot < 0.98:
        p_best['success'] = False
    return p_best


def plot_multiple_bin(simulations, param, err_pfit, h_ex, mech, folderName, plot_options, index=''):
    """
    Plot multiple curves and their fitted universal function after rescaling.

    :param simulations: pandas DataFrame that list all simulations to plot, the values of their
        scaling functions and their errors (each T should appear separately)
    :param param: the found best fitting parameters
    :param err_pfit: (bootstrap) error of the fitted params
    :param h_ex: external Bx magnetic field
    :param mech:    whether internal transverse fields are suppressed (offdiagonal dipolar terms excluded) (="true")
                    or not (offdiagonal dipolar terms included) (="false")
    :param folderName: folder name
    :param plot_options: dict containing plot options
    :param index: (optional) index of this plot
    :return: None
    """
    markers=['o','s','^','D','v']

    # create 2 axes one next to the other:
    # the left one would have the raw curves (still bootstrapped)
    # the right one would have the rescaled curves that should all
    # fall on the same universal curve.
    fig, (ax1, ax2) = plt.subplots(1,2, sharey='row',figsize=(6.5,4))
    # create an additional, rescaled T axis
    simulations['rescaled_T']=(simulations['L']**(1/param[5]))*(simulations['T']-param[0])

    ax1.set_ylabel(plot_options['Name'],fontsize=12)

    # loop over simulations according to their L's
    for (L, df), marker in zip(simulations.groupby(['L']), cycle(markers)):
        data_label = 'L=%d' %L
        # plot crossing
        df.plot(x='T',y='scaling_func', yerr='scaling_func_err', ax=ax1, label=data_label, legend=False, capsize=3, marker=marker, linestyle='-',fillstyle='none')
        ax1.set_xlabel('T')
        
        # plot collapse (use same color as most recently plotted line)
        df.plot(x='rescaled_T',y='scaling_func', yerr='scaling_func_err', ax=ax2, label=data_label, capsize=3, marker=marker, color=ax1.get_lines()[-1].get_color(), linestyle='',fillstyle='none')
        ax2.set_xlabel(r'$L^{1/\nu} (T-T_c)$')

    
    # plot universal function approximation (this is just the polynomial w/ optimal parameters. f(x*L^(\nu)+x_c) is calculated)
    dummy_L=1
    universal_func_to_plot = f_bin(param, np.column_stack((int(dummy_L)*np.ones(50),param[0]+(dummy_L**param[5])*np.linspace(simulations['rescaled_T'].min(),simulations['rescaled_T'].max(), num=50) )))
    p_rescaled=ax2.plot(np.linspace(simulations['rescaled_T'].min(),simulations['rescaled_T'].max()), universal_func_to_plot, '-', color='plum')

    ax1.tick_params(direction='in',which='both')
    ax2.tick_params(direction='in',which='both')

    x = float(folderName.split("_")[-1])


    # add box with found T_c and nu values and the used Bx
    ob = offsetbox.AnchoredText(r'$T_c=%s$' '\n' r'$\nu=%s$' '\n' r'$B_x=%.2g$' '\n' r'$x=%.2g$'%(str_with_err(param[0],err_pfit[0]),str_with_err(param[-1],err_pfit[-1]),h_ex,x), loc=3, prop=dict(size=11))
    ax2.add_artist(ob)
    # set the y axis scale (linear or log)
    ax1.set_yscale(plot_options['axis_yscale'])
    ax2.set_yscale(plot_options['axis_yscale'])
    
    plt.tight_layout()

    # save both a png and an eps file
    fig.savefig('../' + config.system_name + '/figures/fit_%s_%s_%s_%s_%s.eps'%(h_ex,mech,'_'.join(map(str,simulations['L'].unique().tolist())), folderName, index),format='eps')
    fig.savefig('../' + config.system_name + '/figures/fit_%s_%s_%s_%s_%s.pdf'%(h_ex,mech,'_'.join(map(str,simulations['L'].unique().tolist())), folderName, index),format='pdf')
    fig.savefig('../' + config.system_name + '/figures/fit_%s_%s_%s_%s_%s.png'%(h_ex,mech,'_'.join(map(str,simulations['L'].unique().tolist())), folderName, index),dpi=300)


def str_with_err(value, error):
    """
    Format a string of a numerical value with its error in parentheses (only significant digits)
    :param value: numerical value, e.g. 1.5381
    :param error: error in value, e.g. 0.002
    :return: value with error in parentheses, e.g. 1.538(2)
    """
    if np.isnan(error):
        return value
    else:
        digits = -int(math.floor(math.log10(error)))
        return "{0:.{2}f}({1:.0f})".format(value, error*10**digits, digits)


def get_binder(m2,m4,mk2,L):
    """
    Calculate the Binder ratio -- a scaling function.

    :param m2: numpy array of independent values of m^2
    :param m4: numpy array of independent values of m^4
    :param mk2: |m(k_min)|^2: not used but accepted to be compatible with get_correlation_length
    :param L: linear system size: not used but accepted to be compatible with get_correlation_length
    :return: Binder ratio
    """
    m2=np.mean(m2)
    m4=np.mean(m4)
    return 0.5*(3-(m4/(m2**2)))


def get_correlation_length(m2,m4,mk2,L):
    """
    Calculate the finite-size correlation length, divided by the linear system size -- a scaling function.

    See K.-M. Tam and M. J. P. Gingras, Spin-Glass Transition at Nonzero Temperature in a Disordered Dipolar
    Ising System: The Case of LiHo_{x}Y_{1-x}F_{4}, Phys. Rev. Lett. 103, 087202 (2009).
    https://link.aps.org/doi/10.1103/PhysRevLett.103.087202
    and references therein.

    :param m2: m^2: not used but accepted to be compatible with get_binder
    :param m4: m^4: not used but accepted to be compatible with get_binder
    :param mk2: numpy array of independent values |m(k_min)|^2
    :param L: linear system size
    :return: the finite-size correlation length
    """
    m2=np.mean(m2)
    mk2=np.mean(mk2)
    if (m2/mk2)-1 <= 0:
        return 0
    return math.sqrt((m2/mk2)-1)/(2*L*math.sin(math.pi/L))


def parse_simulation_table(simulations):
    """ Validate that the given simulations are suitable for finite size scaling analysis.
        This means that each simulation is of a different L (a must), and that they all have the same mechanism
        toggle and the same B_ex value (not a must but doesn't make much sense otherwise)
    """
    if (simulations['L'] == simulations['L'].iloc[0]).all():
        raise Exception('For finite size scaling the given simulations should be of different linear system sizes (L). \n' + str(simulations))
    if not (simulations['Bex'] == simulations['Bex'].iloc[0]).all():
        raise Exception('For finite size scaling all simulations should have the same B_ex. \n' + str(simulations))
    if not (simulations['mech'] == simulations['mech'].iloc[0]).all():
        raise Exception('For finite size scaling all simulations should have either be with or without the mechanism. \n' + str(simulations))    
    
    xdata={}
    folderName={}
    all_L=[]
    for g in simulations.groupby(['Bex','L','folderName','mech']):
        L = g[1]['L'].iloc[0]
        temperature_list = g[1]['T'].tolist()
        all_L.append(L)
        xdata[L] = temperature_list
        folderName[L] = g[1]['folderName'].iloc[0]
    
    return (all_L, folderName, xdata)
    

def fit_bin(simulations, boot_num, min_x, max_x, initial_xc, fit_options):
    """
    Fit the given MC simulation data to a universal scaling function to find
    the critical temperature and critical exponents.
    Uses the bootstrap method to estimate errors.

    :param simulations: pandas DataFrame that includes all simulation (with a separate row for each T)
    :param boot_num: number of bootstrap samples to perform
    :param min_x: lower bound of T to use
    :param max_x: higher bound of T to use
    :param initial_xc: initial guess for the critical x (T)
    :param fit_options: options for the fitting and plotting (which scaling function to use etc.)
    :return: a 6-tuple of: (T_c, T_c error, nu, nu error, R^2, simulations)
    """
    parse_simulation_table(simulations)
    
    ps=[]
    # use just the defined range of temperatures
    simulations=simulations[(simulations['T']<max_x) & (simulations['T']>min_x)].reset_index(drop=True)

    # create an array for the bootstrap samples
    results=np.zeros((len(simulations),boot_num))
    
    for i, sim in enumerate(simulations.itertuples()):
        y = bin_data.read_binned_data(sim, use_latest=True)
        # for each simulation, calculate the given scaling function
        # using a randomly selected (w/ replacement) sample, and save
        # in a new column of the results array
        for boot_index in range(boot_num):
            results[i,boot_index] = fit_options['func'](y['Magnetization^2'].sample(frac=1,replace=True),y['Magnetization^4'].sample(frac=1,replace=True),y['mk2'+fit_options['corr_length_axis']].sample(frac=1,replace=True),sim.L*fit_options['unit_cell_length'])

    # some bootstrap samples could have extreme values and fitting them is hard.
    # therefore, we start by fitting the curves defined by the bootstrap mean which
    # should be easier. then we will use the values we found as better initial guesses
    # for all of the bootstrap samples.
    # find initial parameters (based on bootstrap mean of the results)
    data = np.column_stack((simulations[['L','T']].to_numpy(),np.mean(results,axis=1))).astype(np.float64)
    initial_fit_successful = False
    initial_fit_attempt = 0
    while not initial_fit_successful:
        initial_fit_param = fit_multiple_bin(data, initial_xc, bounds=False, max_nfev=40000)
        initial_fit_successful = initial_fit_param.success
        initial_fit_attempt += 1
        print('initial fit attempt #%s: R^2=%s, %s'%(initial_fit_attempt, 1-initial_fit_param.cost/(np.sum((data[:,2]-np.mean(data[:,2]))**2)),initial_fit_successful))
        if not initial_fit_successful:
            if initial_fit_attempt < boot_num:
                # try bootstrap samples, one of which might be easier to fit initially
                data = np.column_stack((simulations[['L','T']].to_numpy(),results[:,initial_fit_attempt])).astype(np.float64)
            else:
                # we ran out of bootstrap samples, so add noise to existing bootstrap mean
                print('Not enough bootstrap samples for initial fitting. Using random noise.')
                y_noise = 0.1 * np.random.normal(size=results.shape[0])
                data = np.column_stack((simulations[['L','T']].to_numpy(),y_noise+np.mean(results,axis=1))).astype(np.float64)
        if initial_fit_attempt > 20:
            raise Exception("Unable to find initial fitting parameters. Number of attempts > 20!")
    initial_fit_params = initial_fit_param.x[2:].tolist()   # ignore T_c and p0 for which we have good initial guesses anyway
    print('Initial fitting parameters found: ' + str(initial_fit_params))

    # start the actual bootstrap fitting
    for index, col in enumerate(results.T):
        # now iterating over the bootstrap datasets
        data = np.column_stack((simulations[['L','T']].to_numpy(),col)).astype(np.float64)
        curr_fit_param=fit_multiple_bin(data, initial_xc, bounds=False, input_p_global=initial_fit_params)

        if curr_fit_param.success:
            ps.append(curr_fit_param.x)
        else:
            # try fitting again with bounds on parameters
            # this doesn't work very well. if it happen a lot,
            # it might be necessary to investigate the process
            # and perhaps find better fitting parameters manually.
            # also, maybe manually define the range of temperatures to
            # a sensible value that has enough datapoints within it.
            default_err=np.seterr(all='raise')
            print('could not fit. trying again with bounds.',file=sys.stderr)
            print('info: ' + str(curr_fit_param),file=sys.stderr)
            curr_fit_param=fit_multiple_bin(data, initial_xc, bounds=True, input_p_global=initial_fit_params)
            if curr_fit_param.success:
                ps.append(curr_fit_param.x)
            else:
                warnings.warn("Fitting did not converge for bootstrap index: %s "%boot_index,stacklevel=2)
                simulations['scaling_func'] = data[:,2]
                markers=['o','s','^','D','v']

                fig, ax = plt.subplots(1,1,figsize=(5,5))
                for (L, df), marker in zip(simulations.groupby(['L']), cycle(markers)):
                    data_label = 'L=%d' %L
                    # plot crossing
                    df.plot(x='T',y='scaling_func', yerr='scaling_func_err', ax=ax, label=data_label, capsize=3, marker=marker, linestyle='-')
                    ax.set_xlabel('T')
                plt.legend(loc='best',framealpha=1)
                plt.tight_layout(pad=1.08)
                fig.savefig('../' + config.system_name + '/figures/unsuccessful_fit_%s_%s_%s_%s_%s.png'%(simulations['Bex'].iloc[0],simulations['mech'].iloc[0],'_'.join(map(str,simulations['L'].unique().tolist())), simulations['folderName'].iloc[0], boot_index),dpi=300)

            np.seterr(**default_err)

    # error is 1 sigma:
    Nsigma=1.

    # get bootstrap averages for y's
    mean_y = np.mean(results,axis=1)
    err_y = Nsigma*np.std(results,axis=1)
    
    # add results to simulations
    simulations['scaling_func'] = mean_y
    simulations['scaling_func_err'] = err_y
    
    # get bootstrap averages for parameters
    ps = np.array(ps)
    mean_pfit = np.mean(ps,0)
    err_pfit = Nsigma * np.std(ps,0)

    # plot the data and the fitted function
    plot_multiple_bin(simulations, mean_pfit, err_pfit, simulations['Bex'].iloc[0], simulations['mech'].iloc[0], simulations['folderName'].iloc[0], fit_options)
    pfit_bootstrap = mean_pfit
    perr_bootstrap = err_pfit
    print('parameters: '+str(pfit_bootstrap))
    print('x_c=' + str(pfit_bootstrap[0]) +'+-' + str(perr_bootstrap[0]))
    print('v=' + str(pfit_bootstrap[5]) +'+-' + str(perr_bootstrap[5]))
    
    # calculate R^2
    data_R2 = simulations[['L','T','scaling_func']].to_numpy()
    
    ss_res = np.sum((f_bin(pfit_bootstrap,data_R2[:,0:2]) - data_R2[:,2])**2)
    ss_tot = np.sum((mean_y-np.mean(mean_y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    print('R^2 = %s'%(str(r_squared)))
    
    plt.close()

    return pfit_bootstrap[0], perr_bootstrap[0], pfit_bootstrap[5], perr_bootstrap[5], r_squared, simulations


def parse_arguments():
    """Parse command line arguments. Uses the common parser from config.parse_arguments."""
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    
    parser = ArgumentParser(description="Analyzes Monte Carlo results and plots correlation length curves.", formatter_class=ArgumentDefaultsHelpFormatter, parents=[config.parse_arguments()], conflict_handler='resolve')
    parser.add_argument( "-L", nargs='+', type=int, required=True, help = "Linear system sizes. At least 2 required.")
    parser.add_argument( "-r", "--temperature_range", nargs=2, type=float, required=False, help = ("The temperature range to use when trying to fit."))
    parser.add_argument( "--tmp", default=False, action="store_true", help = "Read from /tmp")
    parser.add_argument( "-i", "--initial_Tc", type=float, help=("Initial guess for T_c. If not given than the middle value between x_min and x_max is used."))
    args = parser.parse_args()

    if len(args.folder_list)>1 and len(args.folder_list)!=len(args.L): 
        parser.error("--folder_list and L argument number mismatch.")
    
    if args.temperature_range is not None and args.temperature_range[0]>args.temperature_range[1]:
        parser.error("First argument of --temperature_range must be smaller than the second argument.")

    config.system_name = args.system_name
    return args


def main():
    args = parse_arguments()
    L = args.L
    boot_num = args.boot_num
    h_ex = args.h_ex
    mech = args.mech
    folderName = args.folder_list

    simulations = analysis_tools.get_simulations(L, folderName, h_ex, mech)
    
    # min & max x for fitting (should be around the apparent Tc)
    if args.temperature_range is not None:
        min_x=args.temperature_range[0]
        max_x=args.temperature_range[1]
    else:
        min_x=simulations['T'].min()
        max_x=simulations['T'].max()
    if min_x>max_x:
        raise Exception("temperature range mismatch. minimum T is larger than maximum T.")
    
    if args.initial_Tc is not None:
        initial_xc=args.initial_Tc
    else:
        initial_xc=0.5*(max_x+min_x)
    
    corr_length_axis='x'
    plot_options = {'Name':r'$\xi^{(%s)}_{L} / L$'%corr_length_axis, 'axis_yscale':'log', 'func':get_correlation_length, 'corr_length_axis':corr_length_axis, 'unit_cell_length':1.0}
    print('\n'.join(map(str,fit_bin(simulations, boot_num, min_x, max_x, initial_xc, plot_options))))


if __name__ == "__main__":
    main()
