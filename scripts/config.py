"""
Configurations common to all Python scripts.

Includes:
- a system_name variable that can be set through other scripts as well.
- a command-line argument parser function with the arguments common to
most scripts (it extended in some of them)
"""
import argparse

system_name = 'LiHoF4'


def parse_arguments():
    """ Common argument parser """
    # parser = ArgumentParser(description="Analyzes Monte Carlo results to create a phase diagram for LiHoF4", formatter_class=ArgumentDefaultsHelpFormatter)
    parser = argparse.ArgumentParser(add_help=False, conflict_handler='resolve')
    parser.add_argument( "-sys", "--system_name", type=str, choices=['LiHoF4','Fe8'], default='LiHoF4', help = "Name of the system. Options are \'LiHoF4\' or \'Fe8\'.")
    parser.add_argument( "-L", nargs='+', type=int, required=True, help = "Linear system sizes.")
    parser.add_argument( "-L", nargs='+', type=int, required=True, help = "Linear system sizes.")
    parser.add_argument( "-b", "--boot_num", type=int, default = 100, help = "Number of bootstrap samples.")
    parser.add_argument( "-h_ex", nargs='+', type=float, help = "External magnetic field values, Bex." , required=True)
    parser.add_argument( "-m", "--mech", nargs='+', choices=['true','false'], help = ("Whether internal fields are suppressed or not. \'false\' means "
                                                                                      "that they aren't so the mechanism is on, and \'true\' means that they are and the mechanism is off." ), required=True)
    parser.add_argument( "-f", "--folder_list", nargs='+', type=str, help = "List of folders in \'data/results\' in which results should be found. " , required=True)

    return parser

