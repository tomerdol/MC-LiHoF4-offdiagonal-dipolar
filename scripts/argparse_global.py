import argparse


parser = argparse.ArgumentParser(add_help=False)
# parser = argparse.ArgumentParser(description="Analyzes Monte Carlo results to create a phase diagram for LiHoF4", formatter_class=ArgumentDefaultsHelpFormatter)
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


all_L = args.L
boot_num = args.boot_num
h_ex_list = args.h_ex_list
mech_list = args.mech
folderName_list = args.folder_list
overwrite_tmp = args.overwrite_tmp