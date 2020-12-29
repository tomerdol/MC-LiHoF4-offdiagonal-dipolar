import matplotlib
matplotlib.use('Agg')
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
import sys
import pandas as pd
import math
import csv
import os
import warnings
import analysis_tools
import bin_data
from itertools import cycle
#print('imports successful')

def f(x, L, par):
    x_c, p0, p1, p2, p3, v = par
    try:
        return p0 + p1*((L**(1/v))*(x-x_c)) + p2*((L**(1/v))*(x-x_c))**2 + p3*((L**(1/v))*(x-x_c))**3
    except Exception as e:
        print(e)
        print(par, file=sys.stderr)

def err(p, L, x, y):
    return f(x,L,p) - y

def err_global(p, all_L, all_x, all_y):
    err0=[]
    for i in range(0, len(all_L)):
        errc=err(p,all_L[i],all_x[all_L[i]],all_y[i])
        err1=np.concatenate((err0,errc))
        err0=err1
    return err0

def fit_multiple(all_xdata, all_ydata, all_L, initial_xc, bounds=False):
    all_xdata = {k:np.array(lst) for k,lst in all_xdata.items()}
    all_ydata = [np.array(lst) for lst in all_ydata]
    max_val = max([np.amax(arr) for arr in all_ydata])
    min_val = min([np.amin(arr) for arr in all_ydata])
    p_global=[ initial_xc, 0.5*(max_val+min_val), -0.4, 0.2, -0.06, 0.6]
    if not bounds:
        p_best=optimize.least_squares(err_global,p_global,args=(all_L,all_xdata,all_ydata))
    else:
        p_best=optimize.least_squares(err_global,p_global,args=(all_L,all_xdata,all_ydata),bounds=([0,min_val,-np.inf,-np.inf,-np.inf,0.0],[3,max_val,np.inf,np.inf,np.inf,np.inf]))
    #print(success)
    #err_toplot=err_global(p_best, all_L, all_xdata, all_ydata)
    return p_best

def f_bin(par, data):
    x_c, p0, p1, p2, p3, v = par
    L=data[:,0]
    x=data[:,1]
    try:
        return p0 + p1*((L**(1/v))*(x-x_c)) + p2*((L**(1/v))*(x-x_c))**2 + p3*((L**(1/v))*(x-x_c))**3
    except Exception as e:
        print(e)
        print(par, file=sys.stderr)

def fit_multiple_bin(data, initial_xc, bounds=False, input_p_global=[ -0.4, 0.2, -0.06, 0.6 ], max_nfev=None):
    # format of data columns: 0-L, 1-T, 2-y
    xdata=data[:,0:2]
    ydata=data[:,2]
    max_val = ydata.max()
    min_val = ydata.min()

    p_global = [initial_xc, 0.5*(max_val+min_val)] + input_p_global

    if not bounds:
        p_best=optimize.least_squares(lambda p, x, y: f_bin(p, x) - y, p_global, args=(xdata,ydata), max_nfev=max_nfev)
    else:
        p_best=optimize.least_squares(lambda p, x, y: f_bin(p, x) - y, p_global, args=(xdata,ydata), bounds=([0,min_val,-np.inf,-np.inf,-np.inf,0.0],[3,max_val,np.inf,np.inf,np.inf,np.inf]))
    #print(success)
    #err_toplot=err_global(p_best, all_L, all_xdata, all_ydata)
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    print('R^2 = ' + str(1-p_best.cost/ss_tot))
    # R^2 < 0.95
    if 1-p_best.cost/ss_tot < 0.95:
        p_best['success'] = False
    return p_best

def bootstrap_resample(X, n=None):
    """ Bootstrap resample an array_like
    Parameters
    ----------
    X : array_like
      data to resample
    n : int, optional
      length of resampled array, equal to len(X) if n==None
    Results
    -------
    returns X_resamples
    """
    if n == None:
        n = len(X)
        
    resample_i = np.floor(np.random.rand(n)*len(X)).astype(int)
    X_resample = X[resample_i]
    return X_resample

# DEPRECATED: A newer version is in the bin directory
def plot_multiple(f, all_L, param, err_pfit, xdata, all_ydata, err_y, h_ex, mech, index=''):
    markers=['o','s','^','D','v']
    toplot=[]
    toplot_rescaled=[]
    fig, (ax1, ax2) = plt.subplots(1,2, sharey='row',figsize=(10,7))
    rescaled_x={L:[L**(1/param[5])*(xi-param[0]) for xi in x] for L, x in xdata.items()}
    ax1.set_ylabel(r'$\xi_{L} / L$')

    for i in range(0,len(all_L)):
        L=all_L[i]
        toplot.append(f(np.linspace(xdata[L][0],xdata[L][-1]),all_L[i],param))
        toplot_rescaled.append(f(np.linspace(rescaled_x[L][0],rescaled_x[L][-1]),all_L[i],param))
    for i in range(0,len(all_L)):
        L=all_L[i]
        data_label = 'L=%d' %all_L[i]
        curve_label = 'L=%d fit' %all_L[i]
        # plot crossing
        # p=ax1.plot(np.linspace(xdata[L][0],xdata[L][-1]), toplot[i], '-')
        ax1.errorbar(xdata[L],all_ydata[L],yerr=err_y[L], fmt=markers[i % len(markers)], capsize=3, label=data_label, color=p[0].get_color() )
        ax1.set_xlabel('T')

        # plot collapse
        ax2.errorbar(rescaled_x[L],all_ydata[L],yerr=err_y[L], fmt=markers[i % len(markers)], capsize=3, label=data_label, color=p[0].get_color() )
        ax2.set_xlabel(r'$L^{1/\nu} (T-T_c)$')
    
    # plot universal function approximation (this is just the polynomial w/ optimal parameters. f(x*L^(\nu)+x_c) is calculated)
    dummy_L=1
    p_rescaled=ax2.plot(np.linspace(min([x[0] for x in rescaled_x.values()]),max([x[-1] for x in rescaled_x.values()])), f(param[0]+(dummy_L**param[5])*np.linspace(min([x[0] for x in rescaled_x.values()]),max([x[-1] for x in rescaled_x.values()])),dummy_L,param), '-', color='plum')

    plt.title(r'Correlation length vs. T ($h_{ex}=$'+str(h_ex)+')')
    #plt.title(r'Binder ratio vs. T ($h_{ex}=$'+str(h_ex)+')')

    #plt.ylabel('g')
    plt.legend(loc='best')
    #axes = plt.gca()
    #axes.set_ylim([0.5,1])
    #ax = fig.add_subplot(111)
    ob = offsetbox.AnchoredText(r'$T_c=%0.2f\pm%0.3f$' '\n' r'$\nu=%0.2f\pm%0.3f$' '\n' r'$B_x=%.2g$'%(round(param[0],2),err_pfit[0],param[-1],err_pfit[-1],h_ex), loc=3,
                                prop=dict(size=11))
    ax1.add_artist(ob)
    #ax2.add_artist(ob)
    ax1.set_yscale("log")
    ax2.set_yscale("log")
    plt.tight_layout(pad=1.08)
    #fig.savefig('../figures/fit_%s_%s_%s_%s.png'%(h_ex,mech,'_'.join(map(str,all_L)),index))
    fig.savefig('../figures/fit_%s_%s_%s_%s.eps'%(h_ex,mech,'_'.join(map(str,all_L)),index),format='eps')


def plot_multiple_bin(f, simulations, param, err_pfit, err_y, h_ex, mech, folderName, plot_options, index=''):
    
    markers=['o','s','^','D','v']
    
    toplot_rescaled=[]
    fig, (ax1, ax2) = plt.subplots(1,2, sharey='row',figsize=(10,7))
    simulations['rescaled_T']=(simulations['L']**(1/param[5]))*(simulations['T']-param[0])
    ax1.set_ylabel(plot_options['Name'],fontsize=16)

    for (L, df), marker in zip(simulations.groupby(['L']), cycle(markers)):
        curve_to_plot = f_bin(param, np.column_stack((int(L)*np.ones(50),np.linspace(df['T'].min(),df['T'].max(),num=50) )))
        data_label = 'L=%d' %L
        curve_label = 'L=%d fit' %L
        
        # plot crossing
        #p=ax1.plot(np.linspace(df['T'].min(),df['T'].max(),num=50), curve_to_plot, '-')
        #df.plot(x='T',y='scaling_func', yerr='scaling_func_err', ax=ax1, label=data_label, capsize=3, marker=marker, color=p[0].get_color(), linestyle='-')
        df.plot(x='T',y='scaling_func', yerr='scaling_func_err', ax=ax1, label=data_label, capsize=3, marker=marker, linestyle='-')
        ax1.set_xlabel('T')
        
        # plot collapse (use same color as most recently plotted line)
        df.plot(x='rescaled_T',y='scaling_func', yerr='scaling_func_err', ax=ax2, label=data_label, capsize=3, marker=marker, color=ax1.get_lines()[-1].get_color(), linestyle='')
        ax2.set_xlabel(r'$L^{1/\nu} (T-T_c)$')

    
    # plot universal function approximation (this is just the polynomial w/ optimal parameters. f(x*L^(\nu)+x_c) is calculated)
    dummy_L=1
    universal_func_to_plot = f_bin(param, np.column_stack((int(dummy_L)*np.ones(50),param[0]+(dummy_L**param[5])*np.linspace(simulations['rescaled_T'].min(),simulations['rescaled_T'].max(), num=50) )))
    p_rescaled=ax2.plot(np.linspace(simulations['rescaled_T'].min(),simulations['rescaled_T'].max()), universal_func_to_plot, '-', color='plum')
    
    plt.legend(loc='best',framealpha=1)
    
    #fig.suptitle(to_plot_disp + ' vs. T ($h_{ex}=$'+str(h_ex)+')')

    #plt.legend(loc='best')

    ob = offsetbox.AnchoredText(r'$T_c=%s$' '\n' r'$\nu=%s$' '\n' r'$B_x=%.2g$'%(str_with_err(param[0],err_pfit[0]),str_with_err(param[-1],err_pfit[-1]),h_ex), loc=3, prop=dict(size=11))
    ax2.add_artist(ob)
    ax1.set_yscale(plot_options['axis_yscale'])
    ax2.set_yscale(plot_options['axis_yscale'])
    
    #fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.tight_layout(pad=1.08)
    
    fig.savefig('../figures/fit_%s_%s_%s_%s_%s.eps'%(h_ex,mech,'_'.join(map(str,simulations['L'].unique().tolist())), folderName, index),format='eps')
    fig.savefig('../figures/fit_%s_%s_%s_%s_%s.png'%(h_ex,mech,'_'.join(map(str,simulations['L'].unique().tolist())), folderName, index),dpi=300)

def str_with_err(value, error):
    digits = -int(math.floor(math.log10(error)))
    return "{0:.{2}f}({1:.0f})".format(value, error*10**digits, digits)


def get_binder(m2,m4,mk2,L):
    m2=np.mean(m2)
    m4=np.mean(m4)
    return 0.5*(3-(m4/(m2**2)))


def get_correlation_length(m2,m4,mk2,L):
    m2=np.mean(m2)
    mk2=np.mean(mk2)
    return math.sqrt((m2/mk2)-1)/(2*L*math.sin(math.pi/L))

# DEPRECATED: A newer version is in the bin directory
def fit_main(all_L, L_equilibrated_min_value, tau_dict, boot_num, h_ex, mech, folderName, xdata, min_x, max_x, initial_xc, folder='../data/results', start_in_middle=True):
    all_yboot={k:[] for k in all_L}
    ps=[]
    xdata = {k:[x for x in v if x > min_x and x<max_x] for k,v in xdata.items()}
    #[x for x in xdata if x > min_x and x<max_x]
    for boot_index in range(boot_num):
        all_ydata=[]
        for L in all_L:
            singleL_ydata=[]
            for T in xdata[L]:
                fname=folder+'/'+folderName[L]+'/table_'+str(L)+'_'+str(L)+'_'+str(h_ex)+'_'+str(T)+'_'+mech+'.txt'
                col_names = pd.read_csv(fname, delim_whitespace=True, error_bad_lines=False, comment='#', nrows=0).columns
                types_dict = {'index': int, 'swap': int}
                types_dict.update({col: np.float64 for col in col_names if col not in types_dict})                                   
                y = pd.read_csv(fname, delim_whitespace=True, error_bad_lines=False, index_col='index', comment='#',dtype=types_dict)
                #y = pd.read_csv(fname, delim_whitespace=True, error_bad_lines=False, comment='#')
                #print('read %s'%fname)
                if (not L in L_equilibrated_min_value) and start_in_middle==True:
                    L_equilibrated_min_value[L]=int(len(y.index)*0.5)
                mk2 = y['mk2x'].loc[L_equilibrated_min_value[L]::tau_dict[L]]
                if 'Magnetization^2' in y.columns: 
                    m2 = y['Magnetization^2'].loc[L_equilibrated_min_value[L]::tau_dict[L]]
                else:
                    m2 = y['Magnetization'].loc[L_equilibrated_min_value[L]::tau_dict[L]]**2
                
                singleL_ydata.append(get_correlation_length(mk2.sample(frac=1,replace=True),m2.sample(frac=1,replace=True),L))  # current y value
                #singleL_ydata.append(get_binder(m2.sample(frac=1,replace=True),(m2**2).sample(frac=1,replace=True)))  # current y value

                del y
                del mk2
                del m2
            #print('ydata (L=%s): %s'%(L,singleL_ydata))
            all_ydata.append(singleL_ydata)
            all_yboot[L].append(singleL_ydata)
            default_err=np.seterr(all='raise')
        curr_fit_param=fit_multiple({k:[x for x in v if x > min_x and x<max_x] for k,v in xdata.items()},[[sub_list[ind] for ind, x in enumerate(xdata[all_L[L]]) if x>min_x and x<max_x] for L, sub_list in enumerate(all_ydata)], all_L, initial_xc, bounds=False)
#    plot_multiple(f, all_L, curr_fit_param, [0 for i in curr_fit_param], xdata, all_ydata, [0 for i in all_ydata], h, h_ex, mech, index=str(boot_index))
        if curr_fit_param.success:
            ps.append(curr_fit_param.x)
        else:
            # try fitting again with bounds on parameters
            print('could not fit. trying again with bounds. L=%s T=%s boot#=%s'%(L, T, boot_index),file=sys.stderr)
            curr_fit_param=fit_multiple({k:[x for x in v if x > min_x and x<max_x] for k,v in xdata.items()},[[sub_list[ind] for ind, x in enumerate(xdata[all_L[L]]) if x>min_x and x<max_x] for L, sub_list in enumerate(all_ydata)], all_L, initial_xc,bounds=True)
            if curr_fit_param.success:
                ps.append(curr_fit_param.x)
            else:
                warnings.warn("Fitting did not converge for bootstrap index: %s "%boot_index,stacklevel=2)
            np.seterr(**default_err)

        #print('boot index: ' + str(boot_index))
    # error is 1 sigma:
    Nsigma=1.
    # get bootstrap averages for y's
    #all_yboot = {L:np.array(all_yboot) for L in all_L} 
    mean_y = {L:np.mean(all_yboot[L],0) for L in all_L}
    err_y = {L:Nsigma*np.std(all_yboot[L],0) for L in all_L}
    #print(err_y)
    # get bootstrap averages for parameters
    ps = np.array(ps)
    mean_pfit = np.mean(ps,0)
    err_pfit = Nsigma * np.std(ps,0)
    plot_multiple(f, all_L, mean_pfit, err_pfit, xdata, mean_y, err_y, h_ex, mech)
    pfit_bootstrap = mean_pfit
    perr_bootstrap = err_pfit
    print('parameters: '+str(pfit_bootstrap))
    print('x_c=' + str(pfit_bootstrap[0]) +'+-' + str(perr_bootstrap[0]))
    print('v=' + str(pfit_bootstrap[5]) +'+-' + str(perr_bootstrap[5]))
    
    # calculate R^2
    ss_res = np.sum(err_global(pfit_bootstrap, all_L, xdata, np.array([mean_y[l] for l in sorted(xdata.keys())]))**2) 
    mean_y0=[]
    for L in all_L:
        mean_yc=mean_y[L]
        mean_y1=np.concatenate((mean_y0,mean_yc))
        mean_y0=mean_y1
    
    ss_tot = np.sum((mean_y0-np.mean(mean_y0))**2)
    r_squared = 1 - (ss_res / ss_tot)
    print('R^2 = %s'%(str(r_squared)))
    
    plt.close()
    
    # returns x_c, err(x_c), v, err(v), R^2
    return pfit_bootstrap[0], perr_bootstrap[0], pfit_bootstrap[5], perr_bootstrap[5], r_squared, [mean_y[l] for l in sorted(xdata.keys())]

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
    parse_simulation_table(simulations)
    
    #all_yboot={k:[] for k in all_L}
    ps=[]
    simulations=simulations[(simulations['T']<max_x) & (simulations['T']>min_x)].reset_index(drop=True)
    
    results=np.zeros((len(simulations),boot_num))
    
    for i, sim in enumerate(simulations.itertuples()):
        y = bin_data.read_binned_data(sim, use_latest=False)

        for boot_index in range(boot_num):
            results[i,boot_index] = fit_options['func'](y['Magnetization^2'].sample(frac=1,replace=True),y['Magnetization^4'].sample(frac=1,replace=True),y['mk2'+fit_options['corr_length_axis']].sample(frac=1,replace=True),sim.L*fit_options['unit_cell_length'])

    # find initial parameters (based on bootstrap mean of the results)
    data = np.column_stack((simulations[['L','T']].to_numpy(),np.mean(results,axis=1))).astype(np.float64)
    initial_fit_successful = False
    initial_fit_attempt=0
    while not initial_fit_successful:
        initial_fit_param=fit_multiple_bin(data, initial_xc, bounds=False, max_nfev=40000)
        initial_fit_successful = initial_fit_param.success
        initial_fit_attempt += 1
        data = np.column_stack((simulations[['L','T']].to_numpy(),results[:,initial_fit_attempt])).astype(np.float64)
        if initial_fit_attempt > 20:
            raise Exception("Unable to find initial fitting parameters. Number of attempts > 20!")
    initial_fit_params = initial_fit_param.x[2:].tolist()   # ignore T_c and p0 for which we have good initial guesses anyway

    for index, col in enumerate(results.T):
        # now iterating over the bootstrap datasets
        data = np.column_stack((simulations[['L','T']].to_numpy(),col)).astype(np.float64)
        curr_fit_param=fit_multiple_bin(data, initial_xc, bounds=False, input_p_global=initial_fit_params)

        if curr_fit_param.success:
            ps.append(curr_fit_param.x)
            print('successful fit: ' + str(curr_fit_param.x))
            #print(curr_fit_param)
            temp_simulations=simulations.copy(deep=True)
            temp_simulations['scaling_func']=col
            temp_simulations['scaling_func_err']=0.001
            plot_multiple_bin(f, temp_simulations, curr_fit_param.x, 0.001*np.ones(len(curr_fit_param.x)), 0.001, simulations['Bex'].iloc[0], simulations['mech'].iloc[0], simulations['folderName'].iloc[0], fit_options,index=index)
        else:
            # try fitting again with bounds on parameters
            default_err=np.seterr(all='raise')
            print('could not fit. trying again with bounds.',file=sys.stderr)
            print('info: ' + str(curr_fit_param),file=sys.stderr)
            curr_fit_param=fit_multiple_bin(data, initial_xc, bounds=True, input_p_global=initial_fit_params)
            if curr_fit_param.success:
                ps.append(curr_fit_param.x)
            else:
                warnings.warn("Fitting did not converge for bootstrap index: %s "%boot_index,stacklevel=2)
            np.seterr(**default_err)

        #print('boot index: ' + str(boot_index))
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
    
    plot_multiple_bin(f, simulations, mean_pfit, err_pfit, err_y, simulations['Bex'].iloc[0], simulations['mech'].iloc[0], simulations['folderName'].iloc[0], fit_options)
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
    
    # returns x_c, err(x_c), v, err(v), R^2
    return pfit_bootstrap[0], perr_bootstrap[0], pfit_bootstrap[5], perr_bootstrap[5], r_squared, simulations

def parse_arguments():  
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    
    parser = ArgumentParser(description="Analyzes Monte Carlo results and plots correlation length curves for LiHoF4", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument( "-L", nargs='+', type=int, required=True, help = "Linear system sizes. At least 2 required.")
    parser.add_argument( "-b", "--boot_num", type=int, default = 100, help = "Number of bootstrap samples.")
    parser.add_argument( "--h_ex", type=float, help = "External magnetic field value, Bex." , required=True)
    parser.add_argument( "-m", "--mech", choices=['true','false'], help = ("Whether internal fields are suppressed or not. \'false\' means "
    "that they aren't so the mechanism is on, and \'true\' means that they are and the mechanism is off." ), required=True)
    parser.add_argument( "-f", "--folder_list", nargs='+', type=str, help = "List of folders in \'data/results/\' in which results should be found. " , required=True)
    parser.add_argument( "-r", "--temperature_range", nargs=2, type=float, required=False, help = ("The temperature range to use when trying to fit."))
    #parser.add_argument( "-o", "--overwrite_tmp", action='store_true', default=False, help = ("Overwrite parsed files in /tmp. If not given, existing files will "
    #"be used (which probably means older results)."))
    parser.add_argument( "--tmp", default=False, action="store_true", help = "Read from /tmp")
    parser.add_argument( "-i", "--initial_Tc", type=float, help=("Initial guess for T_c. If not given than the middle value between x_min and x_max is used."))
    args = parser.parse_args()

    if len(args.folder_list)>1 and len(args.folder_list)!=len(args.L): 
        parser.error("--folder_list and L argument number mismatch.")
    
    if args.temperature_range is not None and args.temperature_range[0]>args.temperature_range[1]:
        parser.error("First argument of --temperature_range must be smaller than the second argument.")
        

    return args



def main():
    args = parse_arguments()
    L = args.L
    boot_num = args.boot_num
    h_ex = args.h_ex
    mech = args.mech
    folderName = args.folder_list
    
    simulations = analysis_tools.get_simulations(L, folderName, h_ex, mech)
    
    #min & max x for fitting (should be around the apparent Tc)
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
    
    #os.system("rsync -avzhe ssh ../figures/ tomerdol@newphysnet1:~/graphs/")
   
if __name__ == "__main__":
    main()
