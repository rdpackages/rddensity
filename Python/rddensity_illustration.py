###########################################################################
## RDDENSITY Python Package
## Script for RDDENSITY Illustration
###########################################################################

from pathlib import Path

import pandas as pd
import plotnine as pn
import rddensity


def show(output):
    repr(output)
    print()


### Load rddensity package
### NOTE: depending on your system, you may need to install it first
# pip install rddensity

### Load data base
data = pd.read_csv(Path(__file__).resolve().with_name("rddensity_senate.csv"))
margin = data["margin"]

### Summary statistics
print(margin.describe())
print()

### rddensity: default options
show(rddensity.rddensity(X=margin))

### rddensity: all statistics & default options
show(rddensity.rddensity(X=margin, useall=True))

### rddensity: default statistic, restricted model & plugin standard errors
show(rddensity.rddensity(X=margin, fitselect="restricted", vce="plugin"))

### rdplotdensity: default options
rdd = rddensity.rddensity(X=margin)
rdplot = rddensity.rdplotdensity(rdd, margin)

### rdplotdensity: custom plot
rdplot_custom = (
    rddensity.rdplotdensity(
        rddensity.rddensity(X=margin),
        X=margin,
        plottype="line",
        xlabel="margin",
        plotRange=[-50, 50],
        plotN=[100],
    )
    + pn.ggtitle("Manipulation testing plot")
)

### rdbwdensity: default options
show(rddensity.rdbwdensity(X=margin))

### rdbwdensity: compute bandwidth and then use them
bw = rddensity.rdbwdensity(X=margin)
hr = bw.h.loc["r", "bw"]
show(rddensity.rddensity(X=margin, h=[10, hr]))

### Other examples
show(rddensity.rddensity(X=margin, kernel="uniform"))
show(rddensity.rddensity(X=margin, bwselect="diff"))
show(rddensity.rddensity(X=margin, h=[10, 15]))
show(rddensity.rddensity(X=margin, p=2, q=4))
show(rddensity.rddensity(X=margin, c=5, useall=True))

show(rddensity.rdbwdensity(X=margin, p=3, fitselect="restricted"))
show(rddensity.rdbwdensity(X=margin, kernel="uniform", vce="jackknife"))
