"""
Create a phase diagram from the results of a MC simulation.
Replaces phase_diagram.py which used results from a single MC
simulation instead of many independent runs.
"""
import numpy as np


def parse_arguments():
    """Parse command line arguments. Uses the common parser from config.parse_arguments."""
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    import config

    parser = ArgumentParser(description="Analyzes Monte Carlo results to create a phase diagram for LiHoF4", formatter_class=ArgumentDefaultsHelpFormatter, parents=[config.parse_arguments()], conflict_handler='resolve')
    parser.add_argument( "-L", nargs='+', type=int, required=True, help = "Linear system sizes. At least 2 required.")
    parser.add_argument( "-o", "--overwrite_tmp", action='store_true', default=False, help = ("Overwrite parsed files in /tmp. If not given, existing files will "
    "be used (which probably means older results)."))
    parser.add_argument( "-s", "--scaling_func", choices=['binder','corr_length'], help = "Scaling function to use for finite-size scaling analysis. Possible choices are \'binder\' for the Binder ratio or \'corr_length\' for the finite-size correlation length"
                                                                                          " divided by the linear system size. Default is \'corr_length_x\'." , required=False, default='corr_length')
    parser.add_argument( "-a", "--corr_length_axis", nargs='?', choices=['x','y','z'], help = "The axis along which the correlation length should be measured for the finite-size correlation length.", required=False, default='x', const='x')
    parser.add_argument( "--delta", help = "The temperature range around the initial critical temperature that will be used for fitting.", required=False, default=0.04)

    args = parser.parse_args()
    if len(args.L)<2:
        parser.error("-L must have at least 2 different system sizes for finite size scaling analysis.")
    if len(args.folder_list)>1 and len(args.folder_list)!=len(args.L): 
        parser.error("--folder_list and -L argument number mismatch.")

    config.system_name = args.system_name
    return args


def find_initial_xc(all_y_curves):
    """Given a pandas DataFrame with scaling function data for various L's,
    find the approximate crossing of the different L curves.
    The data column can be named either 'y_to_plot' or 'scaling_func'.
    This is a wrapper for find_initial_xc_from_arr()."""
    arr_x=[]
    arr_y=[]
    y_column = 'y_to_plot' if 'y_to_plot' in all_y_curves.columns else 'scaling_func'
    # turn the data in the DataFrame into NumPy arrays for find_initial_xc_from_arr()
    for L, L_group in all_y_curves.groupby('L'):
        arr_x.append(L_group['T'].to_numpy())
        arr_y.append(L_group[y_column].to_numpy())

    return find_initial_xc_from_arr(arr_x,arr_y)


def find_initial_xc_from_arr(arr_x,arr_y):
    """
    find average x_c from crossing of all possible pairs of L's
    :param arr_x: list of numpy arrays for x values of the scaling function (each array is for a different L)
    :param arr_y: list of numpy arrays for y values of the scaling function (each array is for a different L)
    :return: approximate crossing of the different curves defined by the input
    """
    import itertools
    
    n=len(arr_y)
    num_of_pairs=int(0.5*n*(n-1))
    sum=0
    # inspect the given curves in pairs
    for pair in itertools.combinations(range(0,n), 2):
        # for each pair, find the approximate crossing point
        sum += find_initial_xc_from_pair(np.array([arr_x[int(pair[0])],arr_y[int(pair[0])]]),np.array([arr_x[int(pair[1])],arr_y[int(pair[1])]]))

    # average the results of the different curves
    initial_xc = sum/num_of_pairs
    return initial_xc


def find_initial_xc_from_pair(arr1,arr2):
    """
    find the crossing of two curves defined by the given numpy arrays.
    the format for both curves is arr[0,:]= xvalues and arr[1,:] = yvalues
    :param arr1: numpy array that defines the first curve.
    :param arr2: numpy array that defines the second curve.
    :return: approximate crossing point of the curves
    """
    # first, find the x value indices of the overlap area between the 2 curves
    idx_x_arr1=np.argwhere((arr1[0,:]>=arr2[0,:].min()) & (arr1[0,:]<=arr2[0,:].max())).flatten()
    idx_x_start = idx_x_arr1[0] # leftmost mutual x
    # interpolate to find values of the second curve for the x values of the first curve
    yinterp = np.interp(arr1[0,idx_x_arr1], arr2[0,:], arr2[1,:])
    # find the index (or indices if there is more than one crossing)
    # where the difference in y values of the two curves switches sign
    idx = np.argwhere(np.diff(np.sign(arr1[1,idx_x_arr1] - yinterp))).flatten() + idx_x_start

    if idx.size > 1:
        # more than one crossing. this is acceptable in some cases so return the intermediate value
        return arr1[0,idx[int(idx.size/2)]]
    elif idx.size==0:
        # no crossing
        if arr1[1,idx_x_start]>yinterp[0]:
            # if the first curve is always higher than the second,
            # then the crossing can happen to the left, so return
            # the leftmost x value of the overlap area
            return arr1[0,idx_x_start]
        else:
            # if the second curve is always higher than the first,
            # then return the rightmost x value of the overlap area
            return arr1[0,idx_x_arr1[-1]]
    else:
        # single clean crossing, return the average of the x values before and after the
        # curves switch places
        return 0.5*(arr1[0,idx[0]]+arr1[0,idx[0]+1])


def plot_previous_data(ax):
    """
    Plot previous data from experiments and numerical works.
    :param ax: axes object in which to plot
    :return: the axes object
    """
    # experimental data
    # P. Babkevich, N. Nikseresht, I. Kovacevic, J. O. Piatek, B. Dalla Piazza, C. Kraemer, K. W. Krämer,
    # K. Prokeš, S. Mat’aš, J. Jensen, and H. M. Rønnow, Phase Diagram of Diluted Ising Ferromagnet
    # LiHo_{x}Y_{1−x}F_4, Phys. Rev. B 94, 174443 (2016).
    # https://link.aps.org/doi/10.1103/PhysRevB.94.174443
    Tc_exp=np.array([1.529253677,1.531944444,1.530647491,1.52935218,1.525398936,1.525430129,1.524136459,1.521524494,1.503278601])
    Bx_exp=np.array([-0.002469136,0.091358025,0.190123457,0.29382716,0.402469136,0.496296296,0.604938272,0.748148148,0.994158712])
    
    # previous simulations
    # S. M. A. Tabei, M. J. P. Gingras, Y.-J. Kao, and T. Yavors’kii,
    # Perturbative Quantum Monte Carlo Study of LiHoF_4 in a Transverse Magnetic Field, Phys. Rev. B 78, 184408 (2008).
    # https://link.aps.org/doi/10.1103/PhysRevB.78.184408
    Tc_sim=np.array([1.531890332,1.495526696,1.470707071])
    Bx_sim=np.array([0.006624937,0.663602684,0.942186845])
    
    ax.plot(Tc_exp,Bx_exp,'m>',label='Experimental data')
    ax.plot(Tc_sim,Bx_sim,'y<',label='Previous simulations')
    
    return ax


def add_overwrite_col(overwrite_tmp, simulations):
    """
    Add an 'overwrite'
    :param overwrite_tmp:
    :param simulations:
    :return:
    """
    table = np.genfromtxt('simulation_plan',comments='$',dtype=str,encoding=None)
    overwrite_tmp_dict={}
    
    for group, _ in simulations.groupby(['Bex','L','folderName','mech']):
        #group[0]=Bex; group[1]=L; group[2]=folderName; group[3]=mech
        matching_row = table[(table[:,6] == group[2]) & (table[:,1] == str(group[1])) & (table[:,2] == str(group[0])) & (table[:,3] == str(group[3]))]
        
        path_to_saved_equilib_file='../' + config.system_name + '/data/results/'+str(group[2])+'/binned_data/equilib_data_'+str(group[1])+'_'+str(group[1])+'_'+str(group[0])+'_'+str(group[3])+'.txt'
        file_exists = os.path.exists(path_to_saved_equilib_file) and os.path.getsize(path_to_saved_equilib_file) > 0
        
        overwrite_tmp_dict[group] = not file_exists or (overwrite_tmp and not (len(matching_row)>0 and matching_row[0][0] == 'done'))
    
    simulations['overwrite']=simulations.apply(lambda row: overwrite_tmp_dict[(row['Bex'],row['L'],row['folderName'],row['mech'])], axis=1)
    return simulations

def validate_simulation_table(simulations):
    if simulations.empty:
        raise Exception('No simulations found matching the given arguments.')
    for Bex, df in simulations.groupby(['folderName','Bex']):
        if (df['L'] == df['L'].iloc[0]).all():
            raise Exception('For finite size scaling there should be different linear system sizes (L) for each given Bex. \n' \
                            + 'for Bex='+str(Bex)+' the given simulations were: \n' + str(df))
    return simulations

def main():
    args = parse_arguments()
    
    all_L = args.L
    boot_num = args.boot_num
    h_ex_list = args.h_ex
    mech_list = args.mech
    folderName_list = args.folder_list
    overwrite_tmp = args.overwrite_tmp
    unit_cell_length_by_axis = {'x':1.0, 'y':1.0, 'z':2.077294686}

    # option packages for plotting and fitting: either binder ratio or finite-size correlation length/L
    if args.scaling_func=='binder':
        plot_options = {'Name':'g', 'axis_yscale':'linear', 'func':fit.get_binder, 'corr_length_axis':args.corr_length_axis, 'unit_cell_length':unit_cell_length_by_axis[args.corr_length_axis]}
    elif args.scaling_func=='corr_length':
        plot_options = {'Name':r'$\xi^{(%s)}_{L} / L$'%args.corr_length_axis, 'axis_yscale':'log', 'func':fit.get_correlation_length, 'corr_length_axis':args.corr_length_axis, 'unit_cell_length':unit_cell_length_by_axis[args.corr_length_axis]}
    else:
        raise Exception('Invalid scaling function given! Either \'binder\' or \'corr_length\' are allowed.') 

    
    simulations = validate_simulation_table(analysis_tools.get_simulations(all_L, folderName_list, h_ex_list, mech_list))
    
    print('Adding overwrite column... ')
    simulations = add_overwrite_col(overwrite_tmp,simulations)
    print('done.')
    
    # re-bin data for those simulations that are marked 'overwrite'. Files should be changed only if there is enough data for new bins.
    print('Binning data... ')
    bin_data.main_bin(simulations[simulations['overwrite']])
    print('done.')
    
    print('Testing equilibration... ')
    to_check = ['Energy','|Magnetization|','Magnetization^2','mk2x']
    simulations = plot_equilibration_pdf_bin.main_check_equilibration(simulations, to_check)
    print('done.')
    # iterate over the 2 mech options. within the loop 'simulations_mech' is a DataFrame for just one of the options
    for mech, simulations_mech in simulations.groupby(['mech']):
        # file to write results to
        f = open("phase_diagram_%s_%s_%s_%s_res.txt"%(config.system_name, mech,'_'.join(map(str,all_L)), '_'.join(map(str,folderName_list))), "w")
        f.write('Tc\tTc_err\tBx\tSimulation_name\n')
        
        fig, ax = plt.subplots()
        ax=plot_previous_data(ax)
        
        # constant shift so that the result for Bex=0.0 is 1.53
        shift=None
        # iterate over the given folderNames's. within the loop 'simulations_mech_folderName' is a DataFrame for just one of the 'mech' options and just one of the folderName's
        # usually this loop (folderName) should have just one iteration. This is unless we want to see different phase diagram boundaries for different simulations (say, different Jex)
        for folderName, simulations_mech_folderName in simulations_mech.groupby(['folderName']):
            print('********************************************************')
            print('Starting work on folderName=%s'%folderName)
            print('********************************************************')
            # iterate over the given Bex's. within the loop 'simulations_mech_folderName_Bex' is a DataFrame for just one of the 'mech' options, just one of the folderName's and just one of the Bex's
            for Bex, simulations_mech_folderName_Bex in simulations_mech_folderName.groupby(['Bex']):
                print('********************************************************')
                print('Starting work on Bex=%s'%Bex)
                print('********************************************************')
                
                try:
                    #tau_dict = {k:1 for k in all_L}
                    print('Plotting finite size scaling graphs to find initial T_c guess...')
                    all_y_curves=plot_bin.main_plot(simulations_mech_folderName_Bex, boot_num, plot_options)
                    #print(all_y_curves)

                    initial_xc = find_initial_xc(all_y_curves)
                    
                    print('found initial T_c guess: %s'%initial_xc)
                except Exception as e:
                    print('error getting initial T_c guess. Skipping Bex=%s'%Bex)
                    print(e)
                    continue
                
                good_fit=False
                delta=min(float(args.delta),min(simulations_mech_folderName_Bex['T'].max()-initial_xc,initial_xc-simulations_mech_folderName_Bex['T'].min()))
                delta_step = delta*0.25
                #delta=0.04
                while not good_fit and delta<20*delta_step:
                    min_x = initial_xc-delta
                    max_x = initial_xc+delta
                    print('Starting fitting...')
                    x_c, x_c_err, v, v_err, r_squared, all_y_curves = fit.fit_bin(simulations_mech_folderName_Bex, boot_num, min_x, max_x, initial_xc, plot_options)

                    #print('x_c=%s, x_c_err=%s, v=%s, v_err=%s, r_squared=%s'%(x_c, x_c_err, v, v_err, r_squared))
                    #print('y_curves:')
                    #print(all_y_curves)
                    if r_squared<0.95 or v<=0 or x_c_err>0.01 or x_c<min_x or x_c>max_x:
                        print('x_c=%s, x_c_err=%s, v=%s, v_err=%s, r_squared=%s, max_x=%s, min_x=%s'%(x_c, x_c_err, v, v_err, r_squared, max_x, min_x))
                        print('Could not find good fit. Trying again.')
                        # retry
                        delta = delta+delta_step
                        try:
                            print(all_y_curves)
                            initial_xc = find_initial_xc(all_y_curves)    # this is the index in the smaller array

                        except Exception as e:
                            print('error finding next T_c guess')
                            print(e)
                            exit()
                        good_fit=False
                        print('new initial_xc=%s'%initial_xc)
                        print('new range=[%s,%s]'%(initial_xc-delta,initial_xc+delta))
                    else:
                        good_fit=True
                        print('Successful fit!')
                        print('x_c=%s, x_c_err=%s, v=%s, v_err=%s, r_squared=%s, max_x=%s, min_x=%s'%(x_c, x_c_err, v, v_err, r_squared, max_x, min_x))
                
                if good_fit:
                    print('Plotting found T_c for Bex=%s'%Bex)
                    if Bex==0:
                        shift=1.53-x_c
                    elif shift==None:
                        shift=0
                    ax.errorbar(x_c+shift,Bex,yerr=None, xerr=x_c_err, label=folderName, color='b')
                    print("%s\t%s\t%s\n"%(x_c,x_c_err,Bex))
                    # write to file
                    f.write("%s\t%s\t%s\t%s\n"%(x_c,x_c_err,Bex,folderName))
                else:
                    print('Could not find T_c for Bex=%s. Skipping.'%Bex)
                    continue
        f.close()
    
        #save fig
        fig.savefig('../' + config.system_name + '/figures/phase_diagram_%s_%s_%s.png'%(mech,'_'.join(map(str,all_L)),folderName_list[0]))
    #os.system("rsync -avzhe ssh ../"+config.system_name+"/figures/ tomerdol@newphysnet1:~/graphs/")
    
if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import os
    import plot_bin, plot_equilibration_pdf_bin, fit, bin_data, analysis_tools
    import config

    main()
