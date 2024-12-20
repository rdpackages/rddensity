#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#ignore future warnings from pandas
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
from scipy.stats import norm
from scipy.stats import binomtest
import scipy.stats as spstat
import math
import statistics as stat
import pandas as pd
import plotnine as pn
from lpdensity import lpdensity

from . import funs
from . import rdbwdensity

def rdplotdensity(rdd, X, plotRange=None, plotN=[10], plotGrid=['es', 'qs'],
                  alpha=0.05, plottype='line', CItype='region', CIuniform=False,
                  CIsimul=2000, CIshade=0.2, bwselect=None,
                  hist=True, histBreaks=None, histfillshade=0.1, histlinecol='white',
                  title=None, xlabel=None, ylabel=None, legendtitle=None, legendgroups=None):

    """

	Parameters
	----------

    rdd : `rddensity` object
       returned by *rddensity*.
    X: Numeric vector or one dimentional matrix/dataframe
       the running variable.
    plotRange: Numeric.
        Specifies the lower and upper bound of the plotting region. Default is *[c-3hl, c+3hr]*
    plotN: Numeric.
        Specifies the bumber of grid points used for plotting on the two sides of the cutoff. Default is 10 on each side.
    plotGrid: String.
        Specifies position of grid points. *es* for evenly spaced, *qs* for quantile spaced.
    alpha: Numeric scalar between 0 and 1.
         The significance level for plotting confidence regions.
    plottype: String.
         *"line"*, *"points"* or *"both"* specifies how the estiamtes are plotted.
    CItype: String.
         *"region"* (default),  *"line"* or *"ebar"*, how the confidence region will be plotted.
    CIuniform: Boolean (default False)
         plotting pointwise confidence intervals or uniform confidence bands.
    CIsimul: Positive integer.
        Number of simulations used to construct confidence intervals (default 2000). Ignored if *CIuniform* is False.
    CIshade: Numeric, between 0 and 1.
        Opaquness of confidence region. Default is *0.2*.
    bwselect: String.
        Method for data-driven bandwidth selection. Default uses bandwidth from *rdd*. *"mse-dpi"*- mean squared error-optimal bandwidth selected for each grid points, *"imse-dpi"*- integrated MSE-optimal bandiwdth, common for all grid points, *"mse-rot"*- rule-of-thumb bandiwdth with Gaussian reference model, *"imse-rot"*-integrated rule-of-thumb bandiwdht with Gaussian reference model.
    hist: Boolean (default True).
        Adds histgram in background of plot.
    histBreaks: Numeric vector.
        breakpoints between histogram bars.
    histfillshade: Numeric between 0 and 1.
        Opaqueness of histrogram. Default is *0.1*.
    title: String.
        Title of the plot
    xlabel: String.
        Label for x-axis.
    ylabel: String.
        Label for y-axis.
    legendTitle: String.
        Title of legend.

    Returns
    -------
    plot: plotnine object.
        Can be customized further with plotnine options.

    See Also
    --------
    rddensity.rddensity
    rddensity.rdbwdensity
    """

    c = rdd.c
    p = rdd.p
    q = rdd.q
    hl = rdd.h['left']
    hr = rdd.h['right']
    kernel = rdd.kernel
    regularize = rdd.regularize
    nLocalMin = rdd.nLocalMin
    nUniqueMin = rdd.nUniqueMin
    massPoints = rdd.massPoints

    # missing value handling
    X = pd.DataFrame(X)
    if X.isnull().values.any():
        warnings.warn(f'{X.isnull().sum()} missing observation(s) are ignored.\n')

    X = X.dropna(axis=0)

    #check plot specifications
    if plotRange is None:
        plotRange = [max(min(X[X.columns[0]]), c-3*hl), min(max(X[X.columns[0]]), c + 3*hr)]
    elif len(plotRange) !=2:
        raise Exception("Plot range incorrectly specified.")
    elif (plotRange[0]>=c or plotRange[1]<= c):
        raise Exception("Plot range incorrectly specified.")

    if len(plotN) == 0:
        plotN = [10, 10]
    elif len(plotN) == 1:
        plotN = [plotN[0], plotN[0]]
    elif len(plotN) >2:
        raise Exception("Number of grid points incorrectly specified.")
    if (plotN[0]<=1 or plotN[1]<=1):
        raise Exception("Number of grid points incorrectly specified.")

    if len(plotGrid) == 0:
        plotGrid = 'es'
    else:
        plotGrid = plotGrid[0]

    if plotGrid not in ['es', 'qs']:
        raise Exception("Grid specification invalid.")

    temp_plot = pn.ggplot() +  pn.theme_bw()
    if hist == True and histBreaks is None:
        Xl = X[(X[X.columns[0]]>= plotRange[0]) & (X[X.columns[0]]<c)]
        Xr = X[(X[X.columns[0]]<= plotRange[1]) & (X[X.columns[0]]>=c)]
        histData = X[(X[X.columns[0]]>= plotRange[0]) & (X[X.columns[0]]<= plotRange[1])]
        temp_hist_n_l = len(Xl)
        temp_hist_n_l = math.ceil(min(np.sqrt(temp_hist_n_l), 10*math.log(temp_hist_n_l)/math.log(10)))
        temp_hist_n_r =  len(Xr)
        temp_hist_n_r = math.ceil(min(np.sqrt(temp_hist_n_r), 10*math.log(temp_hist_n_r)/math.log(10)))
        histBreaks_l = np.linspace(plotRange[0], c, temp_hist_n_l+1)
        histBreaks_r = np.linspace(c, plotRange[1], temp_hist_n_r+1)[1:(temp_hist_n_r+1)]
        histScale = X[(X[X.columns[0]]>= plotRange[0]) & (X[X.columns[0]]<= plotRange[1])].mean()

        temp_plot += pn.geom_histogram(data=histData, mapping=pn.aes(x=histData.columns[0], y=pn.mapping.after_stat('density')), bins=temp_hist_n_l+temp_hist_n_r, alpha=histfillshade)
#        temp_plot += pn.geom_histogram(data=Xr, mapping=pn.aes(x=Xr.columns[0], y=pn.mapping.after_stat('density'), fill=1, alpha=histfillshade/2), bins=temp_hist_n_r)


    scalel = (sum(X[X.columns[0]]<c)-1)/(len(X)-1)
    scaler = (sum(X[X.columns[0]]>=c)-1)/(len(X)-1)

    if plotGrid == 'es':
        gridl = np.linspace(plotRange[0], c, plotN[0]+1)
        gridl[plotN[0]] = c
        gridr = np.linspace(c, plotRange[1], plotN[1]+1)
        gridr[0] = c
    else:
        gridl = np.linspace(len(X[(X[X.columns[0]]<= plotRange[0])])/len(X), len(X[X[X.columns[0]]<c])/len(X), plotN[0]+1)
        gridl = np.quantile(X, gridl)
        gridl[plotN[0]] = c
        gridr = np.linspace(len(X[(X[X.columns[0]]< c)])/len(X), len(X[X[X.columns[0]]<=plotRange[1]])/len(X), plotN[1]+1)
        gridr = np.quantile(X, gridr)
        gridr[0] = c

    #calling lpdensity
    if bwselect is not None:
        if bwselect in ["mse-dpi", "imse-dpi", "mse-rot", "imse-rot"]:
            Estl = lpdensity(data = X[(X[X.columns[0]]< c)], grid=gridl, bwselect=bwselect,
                             p=p, q=q, v=1, kernel=kernel, scale=scalel, regularize=regularize,
                             nLocalMin=nLocalMin, nUniqueMin=nUniqueMin, massPoints=massPoints)
            Estr = lpdensity(data = X[(X[X.columns[0]]>= c)], grid=gridr, bwselect=bwselect,
                             p=p, q=q, v=1, kernel=kernel, scale=scaler, regularize=regularize,
                             nLocalMin=nLocalMin, nUniqueMin=nUniqueMin, massPoints=massPoints)
        else:
            raise Exception("Option bwselect incorrectly specified.")
    else:
        Estl = lpdensity(data = X[(X[X.columns[0]]< c)], grid=gridl, bw=[hl],
                             p=p, q=q, v=1, kernel=kernel, scale=[scalel], regularize=regularize,
                             nLocalMin=nLocalMin, nUniqueMin=nUniqueMin, massPoints=massPoints)
        Estr = lpdensity(data = X[(X[X.columns[0]]>= c)], grid=gridr, bw=[hr],
                             p=p, q=q, v=1, kernel=kernel, scale=[scaler], regularize=regularize,
                             nLocalMin=nLocalMin, nUniqueMin=nUniqueMin, massPoints=massPoints)

    #lpdensity plot
    if alpha<0 or alpha>1:
        raise Exception("Significance level incorrectly specified.")

    data_l = Estl.Estimate[['grid', 'f_p', 'f_q', 'se_p', 'se_q']]
    data_r = Estr.Estimate[['grid', 'f_p', 'f_q', 'se_p', 'se_q']]

    #critical value
    if CIuniform==True:
        if np.isnan(CIsimul):
            raise warnings.warn('Option CIsimul incorrectly specified. Will only plot pointwise confidence intervals.')
            z_val_l = spstat.norm.ppf(1-alpha/2)
            z_val_r = spstat.norm.ppf(1-alpha/2)
        elif math.ceil(CIsimul)<2:
            raise warnings.warn('Option CIsimul incorrectly specified. Will only plot pointwise confidence intervals.')
            z_val_l = spstat.norm.ppf(1-alpha/2)
            z_val_r = spstat.norm.ppf(1-alpha/2)
        else:
            CIsimul = math.ceil(CIsimul)
            corrMat_l = Estl.CovMat_q.mul(1/Estl.Estimate['se_q'], axis=0).mul(1/Estl.Estimate['se_q'], axis=1)
            corrMat_r = Estr.CovMat_q.mul(1/Estr.Estimate['se_q'], axis=0).mul(1/Estr.Estimate['se_q'], axis=1)
            #cehck PSD matrix
            if np.all(np.linalg.eigvals(corrMat_l)>=0):
                #simulate multivariate normal
                normalSimul_l = np.random.multivariate_normal(mean = np.repeat(0, len(Estl.CovMat_q)), cov = corrMat_l, size=CIsimul)
                #comput estimated quantile
                z_val_l = np.quantile(list(map(lambda x: max(abs(x)), normalSimul_l)), 1-alpha)

            else:
                #exception for non-PSD covariance matrix
                raise warnings.warn('Variance-Covariance is not positive semidefinite. Will only plot pointwise confidence intervals.')
                z_val_l = spstat.norm.ppf(1-alpha/2)
            if np.all(np.linalg.eigvals(corrMat_r)>=0):
                #simulate multivariate normal
                normalSimul_r = np.random.multivariate_normal(mean = np.repeat(0, len(Estr.CovMat_q)), cov = corrMat_r, size=CIsimul)
                #comput estimated quantile
                z_val_r = np.quantile(list(map(lambda x: max(abs(x)), normalSimul_r)), 1-alpha)

            else:
                #exception for non-PSD covariance matrix
                raise warnings.warn('Variance-Covariance is not positive semidefinite. Will only plot pointwise confidence intervals.')
                z_val_r = spstat.norm.ppf(1-alpha/2)
    else:
        z_val_l = spstat.norm.ppf(1-alpha/2)
        z_val_r = spstat.norm.ppf(1-alpha/2)

    data_l['CI_l'] = data_l['f_q'] - z_val_l * data_l['se_q']
    data_l['CI_r'] = data_l['f_q'] + z_val_l * data_l['se_q']
    data_r['CI_l'] = data_r['f_q'] - z_val_r * data_r['se_q']
    data_r['CI_r'] = data_r['f_q'] + z_val_r * data_r['se_q']

    if CItype in ['region', 'all']:
        temp_plot = temp_plot + pn.geom_ribbon(data=data_l, mapping=pn.aes(x='grid', ymin='CI_l', ymax='CI_r'), alpha=CIshade)
        temp_plot = temp_plot + pn.geom_ribbon(data=data_r, mapping=pn.aes(x='grid', ymin='CI_l', ymax='CI_r'), alpha=CIshade)

    if CItype in ['lines', 'all']:
        temp_plot = temp_plot + pn.geom_line(data=data_l, mapping=pn.aes(x='grid', y='CI_l')) + pn.geom_line(data=data_l, mapping=pn.aes(x='grid', y='CI_r'))
        temp_plot = temp_plot + pn.geom_line(data=data_r, mapping=pn.aes(x='grid', y='CI_l')) + pn.geom_line(data=data_r, mapping=pn.aes(x='grid', y='CI_r'))

    if CItype in ['ebar', 'all']:
        temp_plot = temp_plot + pn.geom_errorbar(data=data_l, mapping=pn.aes(x='grid', ymin='CI_l', ymax='CI_r'))
        temp_plot = temp_plot + pn.geom_errorbar(data=data_r, mapping=pn.aes(x='grid', ymin='CI_l', ymax='CI_r'))

    if plottype in ['line', 'both']:
        temp_plot = temp_plot + pn.geom_line(data=data_l, mapping=pn.aes(x='grid', y='f_p'))
        temp_plot = temp_plot + pn.geom_line(data=data_r, mapping=pn.aes(x='grid', y='f_p'))

    if plottype in ['points', 'both']:
        temp_plot = temp_plot + pn.geom_point(data=data_l, mapping=pn.aes(x='grid', y='f_p'))
        temp_plot = temp_plot + pn.geom_point(data=data_r, mapping=pn.aes(x='grid', y='f_p'))

    if xlabel is None:
        xlabel = ''

    if ylabel is None:
        ylabel = ''

    if title is None:
        title = ''

    temp_plot = temp_plot + pn.xlab(xlabel) + pn.ylab(ylabel)

    return(temp_plot)
