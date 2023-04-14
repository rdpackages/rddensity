import rddensity
import pandas as pd
import numpy as np

# Load data base
data = pd.read_csv('../R/rddensity_senate.csv')

# rddensity: default options
print(repr(rddensity.rddensity(X=data)))

# rddensity: all statistics & default options
print(repr(rddensity.rddensity(X=data, useall=True)))

# rddensity: default statistic, restricted model & plugin standard errors
print(repr(rddensity.rddensity(X=data, fitselect="restricted", vce="plugin")))
 
### rdplotdensity: default options
rdd = rddensity.rddensity(X=data)
rddensity.rdplotdensity(rdd, data)

# rdbwdensity: default options
print(repr(rddensity.rdbwdensity(X=data)))

# rdbwdensity: compute bandwidth and then use them
bw = rddensity.rdbwdensity(X=data)
hr = bw.h.iloc[1,0]
print(repr(rddensity.rddensity(X=data, h=[10, hr])))

# Other examples
print(repr(rddensity.rddensity(X=data, kernel="uniform")))
print(repr(rddensity.rddensity(X=data, bwselect='diff')))
print(repr(rddensity.rddensity(X=data, h=[10, 15])))
print(repr(rddensity.rddensity(X=data, p=2, q=4)))
print(repr(rddensity.rddensity(X=data, c=5, useall=True)))

print(repr(rddensity.rdbwdensity(X=data, p=3, fitselect='restricted')))
print(repr(rddensity.rdbwdensity(X=data, kernel='uniform', vce='jackknife')))

