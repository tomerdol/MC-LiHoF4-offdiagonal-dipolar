import numpy as np
import pandas as pd
import sys

fname = sys.argv[1]

y=pd.read_csv(fname,header=None,parse_dates=True,infer_datetime_format=True, comment='S')

ts = pd.Series(pd.to_datetime(y[0],errors='coerce'))

avg = (ts-ts.shift(+1)).mean()

print(avg)

