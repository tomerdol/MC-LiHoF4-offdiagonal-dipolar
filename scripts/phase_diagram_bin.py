import numpy as np

def check_exists_and_not_empty(T, L, Bex, folderName, mech):
    L_exists_dict={k:False for k in L}
    for l in L:
        exists=True
        for temperature in T[l]:
            path='/tmp/'+folderName[l]+'/table_'+str(l)+'_'+str(l)+'_'+str(Bex)+'_'+str(temperature)+'_'+mech+'.txt'
            exists = exists and (os.path.exists(path) and os.path.getsize(path) > 0)
                
        L_exists_dict[l] = exists
    return L_exists_dict

def parse_arguments():  
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    
    parser = ArgumentParser(description="Analyzes Monte Carlo results to create a phase diagram for LiHoF4", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument( "-L", nargs='+', type=int, required=True, help = "Linear system sizes. At least 2 required.")
    parser.add_argument( "-b", "--boot_num", type=int, default = 100, help = "Number of bootstrap samples.")
    parser.add_argument( "--h_ex_list", nargs='+', type=float, help = "List of external magnetic field values, Bex." , required=True)
    parser.add_argument( "-m", "--mech", nargs='+', choices=['true','false'], help = ("Whether internal fields are suppressed or not. \'false\' means "
    "that they aren't so the mechanism is on, and \'true\' means that they are and the mechanism is off." ), required=True)
    parser.add_argument( "-f", "--folder_list", nargs='+', type=str, help = "List of folders in \'data/results\' in which results should be found. " , required=True)
    parser.add_argument( "-o", "--overwrite_tmp", action='store_true', default=False, help = ("Overwrite parsed files in /tmp. If not given, existing files will "
    "be used (which probably means older results)."))
    parser.add_argument( "-s", "--scaling_func", choices=['binder','corr_length'], help = "Scaling function to use for finite-size scaling analysis. Possible choices are \'binder\' for the Binder ratio or \'corr_length\' for the finite-size correlation length"
                                                                                          " divided by the linear system size. Default is \'corr_length_x\'." , required=False, default='corr_length')
    parser.add_argument( "-a", "--corr_length_axis", nargs='?', choices=['x','y','z'], help = "The axis along which the correlation length should be measured for the finite-size correlation length.", required=False, default='x', const='x')

    args = parser.parse_args()
    if len(args.L)<2:
        parser.error("-L must have at least 2 different system sizes for finite size scaling analysis.")
    if len(args.folder_list)>1 and len(args.folder_list)!=len(args.L): 
        parser.error("--folder_list and -L argument number mismatch.")
        
    return args

def find_initial_xc(all_y_curves):
    arr_x=[]
    arr_y=[]
    y_column = 'y_to_plot' if 'y_to_plot' in all_y_curves.columns else 'scaling_func'
    for L, L_group in all_y_curves.groupby('L'):
        arr_x.append(L_group['T'].to_numpy())
        arr_y.append(L_group[y_column].to_numpy())

    return find_initial_xc_from_arr(arr_x,arr_y)

def find_initial_xc_from_arr(arr_x,arr_y):
    # find average x_c from crossing of all possible pairs of L's 
    import itertools
    
    n=len(arr_y)
    num_of_pairs=int(0.5*n*(n-1))
    sum=0
    for pair in itertools.combinations(range(0,n), 2):
        sum += find_initial_xc_from_pair(np.array([arr_x[int(pair[0])],arr_y[int(pair[0])]]),np.array([arr_x[int(pair[1])],arr_y[int(pair[1])]]))
    initial_xc = sum/num_of_pairs
    
    return initial_xc

def find_initial_xc_from_pair(arr1,arr2):
    # arr 1&2 are the curves for which a approximate crossing is found
    # the format is arr[0,:]= xvalues and arr[1,:] = yvalues
    # works only with 2 curves
    idx_x_arr1=np.argwhere((arr1[0,:]>=arr2[0,:].min()) & (arr1[0,:]<=arr2[0,:].max())).flatten()
    idx_x_start = idx_x_arr1[0] # leftmost mutual x
    yinterp = np.interp(arr1[0,idx_x_arr1], arr2[0,:], arr2[1,:])
    idx = np.argwhere(np.diff(np.sign(arr1[1,idx_x_arr1] - yinterp))).flatten() + idx_x_start

    #idx = np.argwhere(np.diff(np.sign(arr1 - arr2))).flatten()
    if idx.size > 1:
        return arr1[0,idx[int(idx.size/2)]]
        #raise Exception("did not find unique crossing point between the given curves")
    elif idx.size==0:
        if arr1[1,idx_x_start]>yinterp[0]:
            return arr1[0,idx_x_start]
        else:
            return arr1[0,idx_x_arr1[-1]]
    else:
        return 0.5*(arr1[0,idx[0]]+arr1[0,idx[0]+1])
    
def plot_previous_data(ax):
    # experimental data
    Tc_exp=np.array([1.529253677,1.531944444,1.530647491,1.52935218,1.525398936,1.525430129,1.524136459,1.521524494,1.503278601])
    Bx_exp=np.array([-0.002469136,0.091358025,0.190123457,0.29382716,0.402469136,0.496296296,0.604938272,0.748148148,0.994158712])
    
    #previous simulations
    Tc_sim=np.array([1.531890332,1.495526696,1.470707071])
    Bx_sim=np.array([0.006624937,0.663602684,0.942186845])
    
    ax.plot(Tc_exp,Bx_exp,'m>',label='Experimental data')
    ax.plot(Tc_sim,Bx_sim,'y<',label='Previous simulations')
    
    return ax
    
def copy_files_to_tmp(T, cols_to_copy, L, Bex, folderName, mech, folder='../data/results'):
    path='/tmp/'+folderName
    try:
        os.mkdir(path)
    except OSError:
        print ("Creation of the directory %s failed" % path)
    
    for l in L:
        for temperature in T:
            fname=folder+'/'+folderName+'/table_'+str(l)+'_'+str(l)+'_'+str(Bex)+'_'+str(temperature)+'_'+mech+'.txt'
            dest='/tmp/'+folderName+'/table_'+str(l)+'_'+str(l)+'_'+str(Bex)+'_'+str(temperature)+'_'+mech+'.txt'
            y = pd.read_csv(fname, delim_whitespace=True, error_bad_lines=False, index_col='index', comment='#')
            print('copying file ' + fname + ' to /tmp/')
            y[cols_to_copy].to_csv(dest, sep='\t')
            
            del y

# returns an alternate all_L list with only L's that do not already exist in /tmp.
def create_temp_all_L(all_L, L_exists_dict, overwrite_tmp_dict):
    new_all_L=[]
    #print(overwrite_tmp_dict)
    for l in all_L:
        if (not L_exists_dict[l]) or overwrite_tmp_dict[l]:
            new_all_L.append(l)
    return new_all_L

def add_overwrite_col(overwrite_tmp, simulations):
    table = np.genfromtxt('simulation_plan',comments='$',dtype=str,encoding=None)
    overwrite_tmp_dict={}
    
    for group, _ in simulations.groupby(['Bex','L','folderName','mech']):
        #group[0]=Bex; group[1]=L; group[2]=folderName; group[3]=mech
        matching_row = table[(table[:,6] == group[2]) & (table[:,1] == str(group[1])) & (table[:,2] == str(group[0])) & (table[:,3] == str(group[3]))]
        
        path_to_saved_equilib_file='../data/results/'+str(group[2])+'/binned_data/equilib_data_'+str(group[1])+'_'+str(group[1])+'_'+str(group[0])+'_'+str(group[3])+'.txt'
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
    h_ex_list = args.h_ex_list
    mech_list = args.mech
    folderName_list = args.folder_list
    overwrite_tmp = args.overwrite_tmp
    unit_cell_length_by_axis = {'x':1.0, 'y':1.0, 'z':2.077294686}

    # option packages for plotting and fitting: either binder ratio or finite-size correlation length/L
    if args.scaling_func=='binder':
        plot_options = {'Name':'g', 'axis_yscale':'linear', 'func':fit6.get_binder, 'corr_length_axis':args.corr_length_axis, 'unit_cell_length':unit_cell_length_by_axis[args.corr_length_axis]}
    elif args.scaling_func=='corr_length':
        plot_options = {'Name':r'$\xi^{(%s)}_{L} / L$'%args.corr_length_axis, 'axis_yscale':'log', 'func':fit6.get_correlation_length, 'corr_length_axis':args.corr_length_axis, 'unit_cell_length':unit_cell_length_by_axis[args.corr_length_axis]}
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
        f = open("phase_diagram_%s_%s_%s_res.txt"%(mech,'_'.join(map(str,all_L)), '_'.join(map(str,folderName_list))), "w")
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
                delta=min(0.04,min(simulations_mech_folderName_Bex['T'].max()-initial_xc,initial_xc-simulations_mech_folderName_Bex['T'].min()))
                #delta=0.04
                while not good_fit and delta<0.2:
                    min_x = initial_xc-delta
                    max_x = initial_xc+delta
                    print('Starting fitting...')
                    x_c, x_c_err, v, v_err, r_squared, all_y_curves = fit6.fit_bin(simulations_mech_folderName_Bex, boot_num, min_x, max_x, initial_xc, plot_options)

                    #print('x_c=%s, x_c_err=%s, v=%s, v_err=%s, r_squared=%s'%(x_c, x_c_err, v, v_err, r_squared))
                    #print('y_curves:')
                    #print(all_y_curves)
                    if r_squared<0.95 or v<=0 or x_c_err>0.01 or x_c<min_x or x_c>max_x:
                        print('x_c=%s, x_c_err=%s, v=%s, v_err=%s, r_squared=%s, max_x=%s, min_x=%s'%(x_c, x_c_err, v, v_err, r_squared, max_x, min_x))
                        print('Could not find good fit. Trying again.')
                        # retry
                        delta = delta+0.01
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
        fig.savefig('../figures/phase_diagram_%s_%s_%s.png'%(mech,'_'.join(map(str,all_L)),folderName_list[0]))
    #os.system("rsync -avzhe ssh ../figures/ tomerdol@newphysnet1:~/graphs/")
    
if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pandas as pd
    import csv, sys, os
    import plot, plot_bin, plot_equilibration_pdf_bin, check_equilibration, fit6, test_autocorrelation, bin_data, analysis_tools
    
    from shutil import copyfile
    main()
