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
    return a * x**b

timestamp = datetime(2021, 1, 16).replace(tzinfo=timezone.utc).timestamp() # datetime when timestamps printed to /output were changed. this script (in contrast to get_Avgs_old.py) only searches newer files

max_L = float(sys.argv[1])
config.system_name = sys.argv[2]
output_directory = os.fsencode('../' + config.system_name + '/output/')

averages_dict={}

obsPrintSweepNum={2:11000,3:180000,4:18000,5:3000, 6:1500,7:300,8:100,9:75,10:25}

for file in os.listdir(output_directory):
    fname = os.fsdecode(file)
    if ".o" in fname:
        if os.path.getmtime('../' + config.system_name + '/output/'+fname) > timestamp:
            try:
                #with open('../output/' + fname) as f:
                #    seed = f.readline().strip()
                #stream = os.popen("find ../data/results -name '*" + seed + ".txt' -not -path '*binned_data*' -exec head -n2 {} \\; -quit")
                #start_time = stream.read().strip()
                #start_time = start_time.split('#')[-1]
                y=pd.read_csv('../' + config.system_name + '/output/' + fname,header=None,skiprows=0,parse_dates=True,infer_datetime_format=True, comment='S')
                ts = pd.Series(pd.to_datetime(y[0],errors='coerce'))
                ts = ts.dropna()
                #ts = ts.iloc[:-1]
                #ts = pd.concat([pd.Series(pd.to_datetime(start_time)),ts])

                if len(ts)>1:
                    avg = (ts-ts.shift(+1)).mean()
                    name = fname.partition(".o")[0]
                    if name.startswith('tr'):
                        name = name[2:]
                    elif name.startswith('r'):
                        name = name[1:]
                    name = '_'.join(name.split('_')[:2] + name.split('_')[3:])
                    print(name, avg)
                    if name in averages_dict:
                        averages_dict[name][0] += avg
                        averages_dict[name][1] += 1
                    else:
                        averages_dict[name] = [avg, 1]
            except Exception as e:
                print(e)

averages = pd.DataFrame([key.split('_') + [1000*value[0]/value[1]/obsPrintSweepNum[int(key.split('_')[0])]] for (key,value) in averages_dict.items()], columns = ['L','Bex','mech','runtime'])
averages['L'] = averages['L'].apply(lambda x: float(x))
print(averages)

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
        #df.plot(x='L',y='runtime', ax=ax, label=label)
ax.set_xlabel('L')
ax.set_ylabel('Runtime for 1000 MCS [in hours]')
plt.legend()
fig.savefig('../' + config.system_name + '/figures/runtimes.png',dpi=360)
#for key in averages_dict:
#    avg_time = averages_dict[key][0]/averages_dict[key][1]
#    avg_time_per_thousand_sweeps = 1000*avg_time / obsPrintSweepNum[int(key.split('_')[0][1])]
#    print(key, '->', avg_time_per_thousand_sweeps)

