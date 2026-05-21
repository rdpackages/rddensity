###########################################################################
## RDDENSITY Python Package
## Script for RDDENSITY Plot Illustration
###########################################################################

from pathlib import Path
import math

import numpy as np
import pandas as pd
import plotnine as pn
from lpdensity import lpdensity
from scipy.stats import norm

import rddensity


### Load rddensity package
### NOTE: depending on your system, you may need to install it first
# pip install rddensity

### Load data base
data = pd.read_csv(Path(__file__).resolve().with_name("rddensity_senate.csv"))
margin = data["margin"]

### Plot and save results
rdd = rddensity.rddensity(X=margin)
rdplot = rddensity.rdplotdensity(rdd, X=margin)

### Construct the same plot
c = rdd.c
hl = rdd.h["left"]
hr = rdd.h["right"]
plot_range = [max(margin.min(), c - 3 * hl), min(margin.max(), c + 3 * hr)]
plot_n = [10, 10]

grid_left = np.linspace(plot_range[0], c, plot_n[0] + 1)
grid_left[plot_n[0]] = c
grid_right = np.linspace(c, plot_range[1], plot_n[1] + 1)
grid_right[0] = c

X = pd.DataFrame({"margin": margin}).dropna()
X_left = X[X["margin"] < c]
X_right = X[X["margin"] >= c]
scale_left = [(len(X_left) - 1) / (len(X) - 1)]
scale_right = [(len(X_right) - 1) / (len(X) - 1)]

Estl = lpdensity(
    data=X_left,
    grid=grid_left,
    bw=[hl],
    p=rdd.p,
    q=rdd.q,
    v=1,
    kernel=rdd.kernel,
    scale=scale_left,
    regularize=rdd.regularize,
    nLocalMin=rdd.nLocalMin,
    nUniqueMin=rdd.nUniqueMin,
    massPoints=rdd.massPoints,
)
Estr = lpdensity(
    data=X_right,
    grid=grid_right,
    bw=[hr],
    p=rdd.p,
    q=rdd.q,
    v=1,
    kernel=rdd.kernel,
    scale=scale_right,
    regularize=rdd.regularize,
    nLocalMin=rdd.nLocalMin,
    nUniqueMin=rdd.nUniqueMin,
    massPoints=rdd.massPoints,
)

### CI: below the cutoff
rdplot_left = Estl.Estimate.copy()
rdplot_left["cil"] = rdplot_left["f_q"] - norm.ppf(0.975) * rdplot_left["se_q"]
rdplot_left["ciu"] = rdplot_left["f_q"] + norm.ppf(0.975) * rdplot_left["se_q"]

### CI: above the cutoff
rdplot_right = Estr.Estimate.copy()
rdplot_right["cil"] = rdplot_right["f_q"] - norm.ppf(0.975) * rdplot_right["se_q"]
rdplot_right["ciu"] = rdplot_right["f_q"] + norm.ppf(0.975) * rdplot_right["se_q"]

### Histogram
N_left = int(((margin >= rdplot_left["grid"].min()) & (margin < c)).sum())
hist_n_left = math.ceil(min(math.sqrt(N_left), 10 * math.log(N_left) / math.log(10)))

N_right = int(((margin <= rdplot_right["grid"].max()) & (margin >= c)).sum())
hist_n_right = math.ceil(min(math.sqrt(N_right), 10 * math.log(N_right) / math.log(10)))

hist_breaks = np.concatenate(
    [
        np.linspace(rdplot_left["grid"].min(), c, hist_n_left + 1),
        np.linspace(c, rdplot_right["grid"].max(), hist_n_right + 1)[1:],
    ]
)
hist_scale = ((margin >= rdplot_left["grid"].min()) & (margin <= rdplot_right["grid"].max())).mean()
hist_data = pd.DataFrame({"margin": margin})

### Call plotnine
manual_plot = (
    pn.ggplot()
    + pn.theme_bw()
    + pn.geom_histogram(
        data=hist_data,
        mapping=pn.aes(x="margin", y=pn.after_stat(f"density * {hist_scale:.17g}")),
        breaks=hist_breaks,
        fill="#00BA38",
        color="white",
        alpha=0.2,
    )
    + pn.geom_ribbon(
        data=rdplot_left,
        mapping=pn.aes(x="grid", ymin="cil", ymax="ciu"),
        alpha=0.2,
        fill="black",
    )
    + pn.geom_ribbon(
        data=rdplot_right,
        mapping=pn.aes(x="grid", ymin="cil", ymax="ciu"),
        alpha=0.2,
        fill="red",
    )
    + pn.geom_line(data=rdplot_left, mapping=pn.aes(x="grid", y="f_p"), color="black")
    + pn.geom_line(data=rdplot_right, mapping=pn.aes(x="grid", y="f_p"), color="red")
    + pn.labs(x="margin", y="")
)
