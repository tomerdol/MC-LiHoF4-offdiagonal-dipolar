"""
Calculates running time of the MC simulations based on the timestamps in the .o* files they create.

Usage: python3 get_avgs.py [maximum L to extrapolate to] [system name]
"""
import matplotlib
import config
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import os
from scipy.optimize import curve_fit
from datetime import datetime
from datetime import timezone

def pow_fit(x, a, b):
    """Fitting function: a*x^b"""
    return a * x**b

# earliest datetime to use
timestamp = datetime(2021, 7, 22, 14, 0, 0).replace(tzinfo=timezone.utc).timestamp()

max_L = float(sys.argv[1])  # for extrapolation
config.system_name = sys.argv[2]
output_directory = os.fsencode('../' + config.system_name + '/output/')

# hold the averages: keys are the file names
# values are lists of length 2 that hold the sum
# from all files (0) and the number of files used (1)
# which can be used to calculate the average.
averages_dict={}

# how many sweeps occur between timestamps. should correspond
# to the values in the appropriate parameter_{sys_name}.properties
obsPrintSweepNum={2:11000,3:180000,4:18000,5:3000, 6:1500,7:300,8:100,9:75,10:25}

for file in os.listdir(output_directory):
    fname = os.fsdecode(file)
    if ".o" in fname:
        if os.path.getmtime('../' + config.system_name + '/output/'+fname) > timestamp:
            try:
                # read the timestamps from the .o* file, non-time values will be discarded
                y=pd.read_csv('../' + config.system_name + '/output/' + fname,header=None,skiprows=0,parse_dates=True,infer_datetime_format=True, comment='#')
                ts = pd.Series(pd.to_datetime(y[0],errors='coerce'))
                ts = ts.dropna()

                # if we have just one timestamp or less we have nothing to do with this file
                if len(ts)>1:
                    avg = (ts-ts.shift(+1)).mean()
                    # parse the name of the file, different file names appear as different curves in the results
                    name = fname.partition(".o")[0]
                    if name.startswith('tr'):
                        name = name[2:]
                    elif name.startswith('r'):
                        name = name[1:]
                    name = '_'.join(name.split('_')[:2] + name.split('_')[3:])
                    print(name, avg)
                    # if it's not new, sum it. if it is, add it.
                    if name in averages_dict:
                        averages_dict[name][0] += avg
                        averages_dict[name][1] += 1
                    else:
                        averages_dict[name] = [avg, 1]
            except Exception as e:
                print(e)

# create a DataFrame with the file names and the average run times per 1000 MC sweeps
averages = pd.DataFrame([key.split('_') + [1000*value[0]/value[1]/obsPrintSweepNum[int(key.split('_')[0])]] for (key,value) in averages_dict.items()], columns = ['L','Bex','mech','runtime'])
averages['L'] = averages['L'].apply(lambda x: float(x))
print(averages)

# plot the run times vs. the linear system size, fit a power law and extrapolate up to the given L
fig, ax = plt.subplots(figsize=(8,6))

for (label, df) in averages.groupby(['Bex','mech']):
        runtimes=df['runtime'].to_numpy()/np.timedelta64(1,'h')
        Ls=df['L'].to_numpy(dtype=np.float64)
        line=ax.plot(Ls,runtimes,'o',label=label)
        try:
            if len(Ls)>2:
                fitting_parameters, covariance = curve_fit(pow_fit, Ls, runtimes,p0=[1,9])
                a, b = fitting_parameters
                new_x = np.linspace(Ls.min(),max_L)
                ax.plot(new_x,pow_fit(new_x,a,b),color=line[0].get_color())
        except Exception as e:
            print(e)
ax.set_xlabel('L')
ax.set_ylabel('Runtime for 1000 MCS [in hours]')
plt.legend()
fig.savefig('../' + config.system_name + '/figures/runtimes.png',dpi=360)
