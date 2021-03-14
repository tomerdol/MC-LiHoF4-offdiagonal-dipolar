import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
from matplotlib.ticker import FormatStrFormatter
import config

all_L = [3,4]
h = 0.0
x=1.0
folderName='run_s'
mech='false'
#xdata = [1.0,1.02,1.04261426353,1.065,1.08704450252,1.11,1.13336810343,1.155,1.18166575047,1.2,1.23202156617,1.28452325787,1.33926227049,1.39633394583,1.45583768858,1.5178771395,1.58256035593]

all_Bex= [0.0,0.5]
#all_Bex=[0.0]
fig = plt.figure()
for Bex in all_Bex:
    with open('../' + config.system_name + '/temp_schedule_' + folderName + '_' + str(Bex) + '_' + mech + '.txt','r') as temperature_schedule_file:
        reader = csv.reader(temperature_schedule_file)
        temp_list=list(reader)
    xdata = temp_list[0][:-1]
    xdata=[float(i) for i in xdata]
    for L in all_L:
        singleL_ydata=[]
        for T in xdata:
            fname='../' + config.system_name + '/data/results/'+folderName+'/table_'+str(L)+'_'+str(L)+'_'+str(Bex)+'_'+str(T)+'_'+mech+'.txt'
            y = pd.read_csv(fname, delim_whitespace=True, error_bad_lines=False, comment='#')
            print(fname)
            try:
                if 1 in y.swap.value_counts(normalize=True):
                    rate = y.swap.value_counts(normalize=True)[1]
                else:
                    rate=0
                singleL_ydata.append(rate)  # current y value
            except:
                pass
            del y
        if len(singleL_ydata)>0:
            plt.plot(xdata,singleL_ydata,'*-',label='L=%s,Bex=%s'%(L,Bex))

plt.legend()
plt.title('PT acceptance rate vs. T')
plt.xlabel('T')
plt.ylabel('Acceptance rate')
plt.xticks(rotation=90)
fig.tight_layout()
fig.savefig('../' + config.system_name + '/figures/acceptance_rate.png')



