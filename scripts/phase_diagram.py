import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import csv, sys, os
import plot, plot_equilibration_pdf, check_equilibration, fit6, test_autocorrelation
import itertools
from shutil import copyfile


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
    parser.add_argument( "-m", "--mech", choices=['true','false'], help = ("Whether internal fields are suppressed or not. \'false\' means "
    "that they aren't so the mechanism is on, and \'true\' means that they are and the mechanism is off." ), required=True)
    parser.add_argument( "-f", "--folder_list", nargs='+', type=str, help = "List of folders in \'/data/results\' in which results should be found. " , required=True)
    parser.add_argument( "-o", "--overwrite_tmp", action='store_true', default=False, help = ("Overwrite parsed files in /tmp. If not given, existing files will "
    "be used (which probably means older results)."))
    
    args = parser.parse_args()
    if len(args.L)<2:
        parser.error("-L must have at least 2 different system sizes for finite size scaling analysis.")
    if len(args.folder_list)>1 and len(args.folder_list)!=len(args.L): 
        parser.error("--folder_list and -L argument number mismatch.")


        
    return args

def find_initial_xc_idx(arr_x,arr_y):
    # find average x_c from crossing of all possible pairs of L's 
    n=len(arr_y)
    num_of_pairs=int(0.5*n*(n-1))
    sum=0
    for pair in itertools.combinations(range(0,n), 2):
        sum += find_initial_xc_from_pair(np.array([arr_x[int(pair[0])],arr_y[int(pair[0])]]),np.array([arr_x[int(pair[1])],arr_y[int(pair[1])]]))
    initial_xc = sum/num_of_pairs
    # return the index of the nearest value in the first given array (smallest L)
    return (np.abs(arr_x[0] - initial_xc)).argmin()

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

def overwrite_to_dict(overwrite_tmp,folderName_dict,all_L,h_ex):
    table = np.genfromtxt('simulation_plan',comments='$',dtype=str,encoding=None)
    overwrite_tmp_dict={}
    for l in all_L:
        overwrite_tmp_dict[l] = overwrite_tmp and table[(table[:,6] == folderName_dict[l]) & (table[:,1] == str(l)) & (table[:,2] == str(h_ex))][0][0] != 'done' 
    return overwrite_tmp_dict

def main():
    args = parse_arguments()
    
    all_L = args.L
    boot_num = args.boot_num
    h_ex_list = args.h_ex_list
    mech = args.mech
    folderName_list = args.folder_list
    overwrite_tmp = args.overwrite_tmp
    
    # create L-folder dict
    if len(folderName_list)==1:
        folderName_dict={k:folderName_list[0] for k in all_L}
    else: 
        folderName_dict={k:folderName_list[i] for i,k in enumerate(all_L)}
    
    to_check = ['Energy','Magnetization','Magnetization^2','mk2']
    
    # file to write results to
    f = open("phase_diagram_%s_%s_res.txt"%(mech,'_'.join(map(str,all_L))), "w")
    f.write('Tc\tTc_err\tBx\n')
    
    fig, ax = plt.subplots()
    ax=plot_previous_data(ax)
    
    # constant shift so that the result for Bex=0.0 is 1.53
    shift=None
    for h_ex in h_ex_list:
        print('********************************************************')
        print('Starting work on h_ex=%s'%h_ex)
        print('********************************************************')
        xdata={}
        for L in all_L:
            with open('../temp_schedule_' + folderName_dict[L] + '_' + str(h_ex) + '_' + mech + '.txt','r') as temperature_schedule_file:
                reader = csv.reader(temperature_schedule_file)
                temp_list=list(reader) 
            xdata[L]=temp_list[0]
            xdata[L]=[float(i) for i in xdata[L]]

        # check which L's need to be checked for equilibration, if overwrite_tmp==true all L's are returned
        overwrite_tmp_dict = overwrite_to_dict(overwrite_tmp,folderName_dict,all_L,h_ex)
        print('Starting equilibration tests...')
        temp_all_L=create_temp_all_L(all_L, check_exists_and_not_empty(xdata, all_L, h_ex, folderName_dict, mech), overwrite_tmp_dict)
        print('L=%s already exist in /tmp'%[item for item in all_L if item not in temp_all_L])
        L_equilibrated_min_value=check_equilibration.check_equilib(xdata, to_check, temp_all_L, h_ex, folderName_dict, mech, folder='../data/results')
        print('Finished equilibration tests.')
        print('Starting autocorrelation tests...')
        test_autocorrelation.save_uncorrelated_timeseries(xdata, to_check, temp_all_L, L_equilibrated_min_value, h_ex,
        folderName_dict, mech, folder='../data/results')
        print('Finished autocorrelation tests.')
        L_equilibrated_min_value = {k:0 for k in all_L}
        tau_dict = {k:1 for k in all_L}
        
        if -1 in L_equilibrated_min_value.values():
            print('%s not equilibrated. Skipping.'%h_ex)
            # skip this one
            continue
        try:
            #tau_dict = {k:1 for k in all_L}
            print('Plotting finite size scaling graphs to find initial T_c guess...')
            all_y_curves=plot.main_plot(all_L, L_equilibrated_min_value, tau_dict, boot_num, h_ex, mech, folderName_dict, xdata, folder='/tmp')
            #print(all_y_curves)
            initial_xc_idx = find_initial_xc_idx([xdata[l] for l in sorted(xdata.keys())],all_y_curves)
            initial_xc=0.5*(xdata[all_L[0]][initial_xc_idx] + xdata[all_L[0]][initial_xc_idx+1])
            print('found initial T_c guess: %s'%initial_xc)
        except Exception as e:
            print('error getting initial T_c guess. Skipping h_ex=%s'%h_ex)
            print(e)
            continue
        
        good_fit=False
        delta=0.06
        while not good_fit and delta<0.2:
            min_x = initial_xc-delta
            max_x = initial_xc+delta
            print('Starting fitting...')
            x_c, x_c_err, v, v_err, r_squared, all_y_curves = ( fit6.fit_main(all_L, L_equilibrated_min_value, tau_dict, boot_num, h_ex, 
                                                    mech, folderName_dict, xdata, min_x, max_x, initial_xc, folder='/tmp') )
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
                    initial_xc_idx = find_initial_xc_idx([[x for x in xdata[l] if x>min_x and x<max_x] for l in sorted(xdata.keys())],all_y_curves)    # this is the index in the smaller array
                    initial_xc=0.5*([x for x in xdata[all_L[0]] if x > min_x and x<max_x][initial_xc_idx] + [x for x in
                    xdata[all_L[0]] if x > min_x and x<max_x][initial_xc_idx+1])
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
            print('Plotting found T_c for h_ex=%s'%h_ex)
            if h_ex==0:
                shift=1.53-x_c
            elif shift==None:
                shift=0
            ax.errorbar(x_c+shift,h_ex,yerr=None, xerr=x_c_err, fmt='bs')
            print("%s\t%s\t%s\n"%(x_c,x_c_err,h_ex))
            # write to file
            f.write("%s\t%s\t%s\n"%(x_c,x_c_err,h_ex))
        else:
            print('Could not find T_c for h_ex=%s. Skipping.'%h_ex)
            continue
    f.close()
    
    #save fig
    fig.savefig('../figures/phase_diagram_%s_%s.png'%(mech,'_'.join(map(str,all_L))))
    #os.system("rsync -avzhe ssh ../figures/ tomerdol@newphysnet1:~/graphs/")
    
if __name__ == "__main__":
    main()
