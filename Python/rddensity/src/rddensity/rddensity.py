#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#ignore future warnings from pandas
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
from scipy.stats import norm
from scipy.stats import binomtest
import math
import scipy.integrate as integrate
from numpy.linalg import inv
import scipy.optimize as optimize
import statistics as stat
import pandas as pd

from . import funs
from . import rdbwdensity


def rddensity(X, c=0, p=2, q=0,
              fitselect='unrestricted', kernel='triangular', vce='jackknife',
              h=[], bwselect='comb', useall=False,
              massPoints=True, regularize=True, nLocalMin=None, nUniqueMin=None,
              bino_flag=True, binoW=None, binoN=None, binoWStep=None, binoNStep=None, binoNW=10, binoP=[0.5]):
    r"""

	Parameters
	----------
    X: Numeric vector or one dimentional matrix/dataframe
        the running variable.
    c: Numeric
       Specifies the threshold or cutoff value in the support of *X*. Default is *0*.
    p: Nonnegative integer
        specifies the local polynomial order used to construct the density estimators. Default is *2* (local quadratic approximation).
    fitselect: String
        specifies the density estimation method. *unrestricted* (Default) for density estimation without any restrictions (two-sample, unrestricted inference). *restricted* for density estimation assuming equal distribution function and higher order dericatives.
    kernel: String
        specifies the kernel function used to construct the local polynomial estimators. Accepted kernels: *triangular* (Default), *epanechnikov* or *uniform*.
    vce: String
        specifies the procedure used to compute the variance-covariance matrix estimatior. *jackknife* (Default) for jackknife standard errors or *plugin* for asymptotic plug-in standard errors.
    massPoints: Boolean, Default *True*.
        Specifies wether to adjust for mass points in the data.
    useall: Boolean, Default *False*.
        If specified, will report two testing procedures: conventional test statistic (not valie when useing mse-optimal bandwidth) and robust bias-corrected statistic.
    h: Numeric
        Specifies the bandwidth used to construct the density estimators on the two sides of the cutoff. If not specified, the bandwidth h is computed using the companion function, *rdbwdensity*. If two bandwidths are specified, the first bandwidth is used for the data below the cutoff and the second bandwidth is used for the data above the cutoff.
    bwselect: String.
        Specified the bandwidth selection procedure to be used. *each*-based on MSE of each density estimator separately (two distinct bandwidths), *diff*- based on MSE of difference of two density estimators, gives one common bandwidth, *sum*-based on MSE of sum of two density estimators, gives one common bandwidth. *comb* (default)-bandiwdth is selected as a combination of the alternatives above. For *fitselect='unrestricted'*, it selects *median(each, diff, sum)*. For *fitselect='restricted'*, it selects *min(diff, sum)*
    regularize: Boolean, Default *True*.
        Specifies whether to conduct local sample size checking. When True, the bandwidth is chosen such that the local region includes at least *nLocalMin* observations and at least *nUniqueMin* unique observations.
    nLocalMin: Nonnegative integer
        Specifies the minimum number of observations in each local neighbourhood. This option will be ignored if set to *0* or if *regularize=False*. Default is *20+p+1*.
    nUniqueMin: Nonnegative integer
        Specifies the minimum number of unqieu observations in each local neighbourhood. This option will be ignored if set to *0* or if *regularize=False*. Default is *20+p+1*.
    bino: Boolean (Default True).
        Specifies whether to conduct binomial tests. By default the initial (smallest) window contains at least 20 observations on each side, and its length is also used as the increment for subsequent windows.
    binoW: Numeric.
        Specifies the half length(s) of the initial window. If two values are provided, they will be used for the data below and above the cutoff separately.
    binoN: Nonnegative integer.
        Specifies the minimum number of observations on each side of the cutoff used for the binomial test. This is ignored if *binoW* is provided.
    binoWStep: Numeric.
        Specifies the increment in half lengths.
    binoNStep: Nonnegative integer.
        Specifies the minimum increment in sample size (on each side of the cutoff). This is ignored if *binoWStep* is provided.
    binoNW: Nonnegative integer.
        Specifies the total number of windows. Default is *10*.
    binoP: Numeric.
        Specifies the null hypothesis of the binomial test. Default is *0.5*.

    Returns
    -------
    hat:
        left/right: density estimate to the left/right of the cutoff. diff: difference in estimated densities on the two sides of the cutoff.
    sd_asy:
        left/right: standard error for the estimated density to the left/right of the cutoff, diff: standard error for difference in estimated densities on the two sides of the cutoff. (based on asymptotic method)
    sd_jk:
        left/right: standard error for the estimated density to the left/right of the cutoff, diff: standard error for difference in estimated densities on the two sides of the cutoff. (based on jackknife method)
    test:
        t_asy/t_jk: t statistic for the density discontinuity test. p_asy/p_jk: p-value for the density discontinuity test.
    hat_p:
        Same as hat, without bias correction.
    bino:
        Binomial test results.
    h:
        bandwidth used to the left/right of the cutoff.
    n:
        full: full sample size, left/right: sample size to the left/right of the cutoff.
    X_min:
        Smallest observations to the left and right of the cutoff.
    X_max:
        Largest observations to the left and right of the cutoff.
    options:
        other options passed to the function are also stored within the object.

    See Also
    --------
    rddensity.rdbwdensity
    rddensity.rdplotdensity
    """
    if q==0:
        q = p+1
    #bandwidth
    if len(h) == 0:
        hl = 0
        hr = 0
    elif len(h) == 1:
        if h<=0:
            raise Exception("Bandwidth has to be positive.")
        else:
            hl = hr = h
    elif len(h) == 2:
        if min(h) <=0:
            raise Exception("Bandwidth has to be positive.")
        else:
            hl = h[0]
            hr = h[1]
    else:
        raise Exception("No more than two bandwidths are accepted.")

    # missing value handling
    X = pd.DataFrame(X)
    if X.isnull().values.any():
        warnings.warn(f'{X.isnull().sum()} missing observation(s) are ignored.\n')

    X = X.dropna(axis=0)

    #sample sizes
    X = X.sort_values(by=X.columns[0], ignore_index=True)
    n = len(X)
    nl = sum(X[X.columns[0]]<c)
    nr = sum(X[X.columns[0]]>=c)
    Xmin = X.min()[0]
    Xmax = X.max()[0]
    XUnique = funs.__rddensityUnique(X)
    freqUnique = XUnique["freq"]
    indexUnique = XUnique["indexLast"]
    XUnique = XUnique["unique"]
    nUnique = len(XUnique)
    nlUnique = sum(XUnique < c)[0]
    nrUnique = sum(XUnique >= c)[0]

    if (nUnique!=n & massPoints)==True:
        massPoints_flag = True
    else:
        massPoints_flag = False
    #Error Handling
    if c<Xmin or c>Xmax:
        raise Exception("The cutoff should be set within the range of the data.")
    if p<1 or p>7:
        raise Exception("p must be an integer between 1 and 7")
    if p>q:
        raise Exception("q cannot be smaller than p")
    if kernel not in ["triangular", "uniform", "epanechnikov"]:
        raise Exception("kernel incorrectly specified")
    if fitselect not in ['restricted', 'unrestricted']:
        raise Exception("fitselect incorrectly specified.")
    if vce not in ['plugin', 'jackknife']:
        raise Exception("vce incorrectly specified.")

    #Regularize
    if regularize is None:
        regularize = True
    elif type(regularize) != bool:
        raise Exception("Regularization parameter incorrectly specified.")

    # nLocalMin
    if nLocalMin is None:
        nLocalMin = 20 + p + 1
    if not np.isreal(nLocalMin) or np.isnan(nLocalMin):
        raise Exception("Option nLocalMin incorrectly specified.")
    elif math.ceil(nLocalMin)<0:
        raise Exception("Option nLocalMin incorrectly specified.")
    else:
        nLocalMin = math.ceil(nLocalMin)

    # nUniqueMin
    if nUniqueMin is None:
        nUniqueMin = 20 + p + 1
    if not np.isreal(nUniqueMin) or np.isnan(nUniqueMin):
        raise Exception("Option nUniqueMin incorrectly specified.")
    elif math.ceil(nUniqueMin)<0:
        raise Exception("Option nUniqueMin incorrectly specified.")
    else:
        nUniqueMin = math.ceil(nUniqueMin)

    # massPoints
    if massPoints is None:
        massPoints = True
    elif type(massPoints) != bool:
        raise Exception("Option massPoints incorrectly specified.")

    # bandwidth selection
    if (int(hl)>0 & int(hr>0)):
        bwselectl = "manual"
    else:
        bwselectl = "estimated"
        out = rdbwdensity.rdbwdensity(X=X, c=c, p=p, kernel=kernel, fitselect=fitselect,
                          vce=vce, regularize=regularize, nLocalMin=nLocalMin,
                          nUniqueMin=nUniqueMin, massPoints=massPoints).h
        if (fitselect=='unrestricted'):
            if bwselect=='each':
                if hl==0:
                    hl = out['bw'][0]
                if hr==0:
                    hr = out['bw'][1]
            elif bwselect=='diff':
                if hl==0:
                    hl = out['bw'][2]
                if hr==0:
                    hr = out['bw'][2]
            elif bwselect=='sum':
                if hl==0:
                    hl = out['bw'][3]
                if hr==0:
                    hr = out['bw'][3]
            elif bwselect=='comb':
                if hl==0:
                    hl = stat.median([out['bw'][0], out['bw'][2], out['bw'][3]])
                if hr==0:
                    hr = stat.median([out['bw'][1], out['bw'][2], out['bw'][3]])
        elif fitselect=='restricted':
            if bwselect=='diff':
                if hl==0:
                    hl = out['bw'][2]
                if hr==0:
                    hr = out['bw'][2]
            elif bwselect=='sum':
                if hl==0:
                    hl = out['bw'][3]
                if hr==0:
                    hr = out['bw'][3]
            elif bwselect=='comb':
                if hl==0:
                    hl = min([out['bw'][2], out['bw'][3]])
                if hr==0:
                    hr = min([out['bw'][2], out['bw'][3]])


    #data trimming
    X = X - c
    Y = np.array(range(n))/(n-1)
    if massPoints==True:
        Y = np.repeat(Y[indexUnique], freqUnique)

    Xh = X[(X>=-1*hl) & (X<=hr)].dropna()
    Yh = pd.DataFrame(Y[Xh.index])

    nlh = sum(Xh[Xh.columns[0]]<0)
    nrh = sum(Xh[Xh.columns[0]]>=0)
    nh = nrh + nlh


    #estimation
    fV_q = funs.__rddensity_fv(Y=Yh, X=Xh, nl=nl, nr=nr, nlh=nlh, nrh=nrh, hl=hl, hr=hr, p=q, s=1, kernel=kernel, fitselect=fitselect, vce=vce, massPoints=massPoints)
    T_asy = fV_q['hat'][2]/np.sqrt(fV_q['plugin'][2])
    T_jk = fV_q['hat'][2]/np.sqrt(fV_q['jackknife'][2])
    p_asy = 2*(1-norm.cdf(abs(T_asy)))
    p_jk = 2*(1-norm.cdf(abs(T_jk)))

    #binomial testing
    if bino_flag==True:
        XSort = abs(X).sort_values(X.columns[0]).reset_index(drop=True)
        XL = abs(X[X[X.columns[0]]<0]).reset_index(drop=True)
        XR = X[X[X.columns[0]]>=0].reset_index(drop=True)


        if isinstance(binoNW, (int, float)):
            if binoNW<=0:
                raise Exception("Option binoNW incorrectly specified.")
        else:
            binoNW = math.ceil(binoNW)


        binomTempLW = np.repeat(None, binoNW)
        binomTempRW = np.repeat(None, binoNW)

        #binoP check
        if len(binoP)>1:
            raise Exception("Option binoP incorrectly specified.")
        elif not isinstance(binoP[0], (int, float)):
            raise Exception("Option binoP incorrectly specified.")
        elif binoP[0] < 0 or binoP[0] > 1:
            raise Exception("Option binoP incorrectly specified.")

        #binoW and binoN check
        if binoW is None:
            if binoN is None:
                binoN = 20
                binomTempLW[0] = max(XL[XL.columns[0]][min(binoN, nl)-1], XR[XR.columns[0]][min(binoN, nr)-1])
                binomTempRW[0] = max(XL[XL.columns[0]][min(binoN, nl)-1], XR[XR.columns[0]][min(binoN, nr)-1])
            elif isinstance(binoN, (int, float)):
                if binoN <= 0 :
                    raise Exception("Option binoN incorrectly specified.")
                binoN = math.ceil(binoN)
                binomTempLW[0] = XSort[XSort.columns[0]][min(binoN, n)-1]
                binomTempRW[0] = XSort[XSort.columns[0]][min(binoN, n)-1]
            else:
                raise Exception("Option binoN incorrectly specified.")
        elif isinstance(binoW, (int, float)):
            if binoW <= 0:
                raise Exception("Option binoW incorrectly specified.")
            binomTempLW[0] = binoW
            binomTempRW[0] = binoW
            binoN = min(sum(XL[XL.columns[0]] <= binomTempLW[0]), sum(XR[XR.columns[0]] <= binomTempRW[0]))
        elif len(binoW) == 2:
            if min(binoW) <= 0:
                raise Exception("Option binoW incorrectly specified.")
            binomTempLW[0] = binoW[0]
            binomTempRW[0] = binoW[1]
            binoN = min(sum(XL[XL.columns[0]] <= binomTempLW[0]), sum(XR[XR.columns[0]] <= binomTempRW[0]))
        else:
            raise Exception("Option binoW incorrectly specified.")


        #binoWStep check
        if binoNW >1 :
            if binoWStep is None:
                if binoNStep is None:
                    if (binomTempLW[0] >= hl or binomTempRW[0] >= hr):
                        binomTempLW[:] = binomTempLW[0]
                        binomTempRW[:] = binomTempRW[0]
                        binoNW = 1
                    else:
                        if binomTempLW[0]*binoNW > hl:
                            binomTempLW[1:(binoNW)] = binomTempLW[0] + np.array(range(1, binoNW))*(hl-binomTempLW[0])/(binoNW-1)
                        else:
                            binomTempLW[1:(binoNW)] = binomTempLW[0] + np.array(range(1, binoNW))*binomTempLW[0]

                        if binomTempRW[0]*binoNW > hr:
                            binomTempRW[1:(binoNW)] = binomTempRW[0] + np.array(range(1, binoNW))*(hr-binomTempRW[0])/(binoNW-1)
                        else:
                            binomTempRW[1:(binoNW)] = binomTempRW[0] + np.array(range(1, binoNW))*binomTempRW[0]

                elif len(binoNStep) == 1:
                    if binoNStep[0] <= 0:
                        raise Exception("Option binoNStep incorrectly specified.")
                    binoNStep = math.ceil(binoNStep)
                    for jj in range(1, binoNW):
                        binomTempLW[jj] = binomTempLW[jj-1] + max(XL[XL.columns[0]][min(sum(XL<=binomTempLW[jj-1]) + binoNStep, nl)-1] - binomTempLW[jj-1],
                                                                  XR[XR.columns[0]][min(sum(XL<=binomTempRW[jj-1]) + binoNStep, nr)-1] - binomTempRW[jj-1])
                        binomTempRW[jj] = binomTempRW[jj-1] + max(XL[XL.columns[0]][min(sum(XL<=binomTempLW[jj-1]) + binoNStep, nl)-1] - binomTempLW[jj-1],
                                                                  XR[XR.columns[0]][min(sum(XL<=binomTempRW[jj-1]) + binoNStep, nr)-1] - binomTempRW[jj-1])
                else:
                    raise Exception("Option binoNStep incorrectly specified.")
            elif len(binoWStep) == 1:
                if binWStep[0] <= 0:
                    raise Exception("Option binoWStep incorrectly specified.")
                binomTempLW[1:(binoNW-1)] = binomTempLW[0] + np.array(range(binoNW))*binoWStep
                binomTempRW[1:(binoNW-1)] = binomTempRW[0] + np.array(range(binoNW))*binoWStep
            elif len(binoWStep) == 2:
                if min(binoWStep) <= 0:
                    raise Exception("Option binoWStep incorrectly specified.")
                binomTempLW[1:(binoNW-1)] = binomTempLW[0] + np.array(range(binoNW))*binoWStep[0]
                binomTempRW[1:(binoNW-1)] = binomTempRW[0] + np.array(range(binoNW))*binoWStep[1]
            else:
                raise Exception("Option binoWStep incorrectly specified.")

            binomTempLN = np.repeat(None, binoNW)
            binomTempRN = np.repeat(None, binoNW)
            binomTempPVal = np.repeat(None, binoNW)
            for jj in range(binoNW):
                binomTempLN[jj] = len(XL[XL[XL.columns[0]] <= binomTempLW[jj]])
                binomTempRN[jj] = len(XR[XR[XR.columns[0]] <= binomTempRW[jj]])
                binomTempPVal[jj] = binomtest(k=binomTempLN[jj], n =(binomTempLN[jj]+binomTempRN[jj]), p = binoP[0]).pvalue

        else:
            binomTempLN = binomTempRN = binomTempLW = binomTempRW = binomTempPVal = np.nan
    else:
        binomTempLN = binomTempRN = binomTempLW = binomTempRW = binomTempPVal = np.nan


    if useall == True:
        fV_p = funs.__rddensity_fv(Y=Yh, X=Xh, nl=nl, nr=nr, nlh=nlh, nrh=nrh, hl=hl, hr=hr, p=p, s=1, kernel=kernel, fitselect=fitselect, vce=vce, massPoints=massPoints)
        T_asy_p = fV_p['hat'][2]/np.sqrt(fV_p['plugin'][2])
        T_jk_p = fV_p['hat'][2]/np.sqrt(fV_p['jackknife'][2])
        p_asy_p = 2*(1-norm.cdf(abs(T_asy_p)))
        p_jk_p = 2*(1-norm.cdf(abs(T_jk_p)))

        result = CJMrddensity(
            hat = pd.Series(data={'left': fV_q['hat'][0], 'right': fV_q['hat'][1], 'diff':fV_q['hat'][2]}, index = ['left', 'right', 'diff']),
            sd_asy = pd.Series(data={'left': np.sqrt(fV_q['plugin'][0]), 'right': np.sqrt(fV_q['plugin'][1]), 'diff':np.sqrt(fV_q['plugin'][2])}, index = ['left', 'right', 'diff']),
            sd_jk = pd.Series(data={'left': np.sqrt(fV_q['jackknife'][0]), 'right': np.sqrt(fV_q['jackknife'][1]), 'diff':np.sqrt(fV_q['jackknife'][2])}, index = ['left', 'right', 'diff']),
            test = pd.Series(data={'t_asy': T_asy, 't_jk':T_jk, 'p_asy':p_asy, 'p_jk':p_jk}, index=['t_asy', 't_jk', 'p_asy', 'p_jk']),
            hat_p = pd.Series(data={'left': fV_p['hat'][0], 'right': fV_p['hat'][1], 'diff':fV_p['hat'][2]}, index = ['left', 'right', 'diff']),
            sd_asy_p = pd.Series(data={'left': np.sqrt(fV_p['plugin'][0]), 'right': np.sqrt(fV_p['plugin'][1]), 'diff':np.sqrt(fV_p['plugin'][2])}, index = ['left', 'right', 'diff']),
            sd_jk_p = pd.Series(data={'left': np.sqrt(fV_p['jackknife'][0]), 'right': np.sqrt(fV_p['jackknife'][1]), 'diff':np.sqrt(fV_p['jackknife'][2])}, index = ['left', 'right', 'diff']),
            test_p = pd.Series(data={'t_asy': T_asy_p, 't_jk':T_jk_p, 'p_asy':p_asy_p, 'p_jk':p_jk_p}, index=['t_asy', 't_jk', 'p_asy', 'p_jk']),
            n = pd.Series(data={'full':n, 'left':nl, 'right':nr, 'eff_left':nlh, 'eff_right':nrh}, index=['full', 'left', 'right', 'eff_left', 'eff_right']),
            h = pd.Series(data={'left':hl, 'right':hr}, index=['left', 'right']),
            fitselect=fitselect, kernel=kernel, vce=vce, c=c, p=p, q=q, useall=useall, bino_flag=bino_flag,
            regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin,
            massPoints=massPoints, massPoints_flag=massPoints_flag, bwselectl=bwselectl, bwselect=bwselect,
            binoN=binoN, binoW=binoW, binoNStep=binoNStep, binoWStep=binoWStep, binoNW= binoNW, binoP=binoP,
            X_min=pd.Series(data={'left': min(np.array(X[X<0].dropna()+c))[0],
                                  'right': min(np.array(X[X>=0].dropna()+c))[0]}, index=['left', 'right']),
            X_max=pd.Series(data={'left': max(np.array(X[X<0].dropna()+c))[0],
                                  'right': max(np.array(X[X>=0].dropna()+c))[0]}, index=['left', 'right']),
            bino = pd.Series(data={'leftN':binomTempLN, 'rightN':binomTempRN, 'leftWindow':binomTempLW, 'rightWindow':binomTempRW, 'pval':binomTempPVal}, index=['leftN', 'rightN', 'leftWindow', 'rightWindow', 'pval'])
        )
    else:
        result = CJMrddensity(
            hat = pd.Series(data={'left': fV_q['hat'][0], 'right': fV_q['hat'][1], 'diff':fV_q['hat'][2]}, index = ['left', 'right', 'diff']),
            sd_asy = pd.Series(data={'left': np.sqrt(fV_q['plugin'][0]), 'right': np.sqrt(fV_q['plugin'][1]), 'diff':np.sqrt(fV_q['plugin'][2])}, index = ['left', 'right', 'diff']),
            sd_jk = pd.Series(data={'left': np.sqrt(fV_q['jackknife'][0]), 'right': np.sqrt(fV_q['jackknife'][1]), 'diff':np.sqrt(fV_q['jackknife'][2])}, index = ['left', 'right', 'diff']),
            test = pd.Series(data={'t_asy': T_asy, 't_jk':T_jk, 'p_asy':p_asy, 'p_jk':p_jk}, index=['t_asy', 't_jk', 'p_asy', 'p_jk']),
            hat_p = pd.Series(data={'left': np.nan, 'right': np.nan, 'diff':np.nan}, index = ['left', 'right', 'diff']),
            sd_asy_p = pd.Series(data={'left': np.nan, 'right': np.nan, 'diff':np.nan}, index = ['left', 'right', 'diff']),
            sd_jk_p = pd.Series(data={'left': np.nan, 'right': np.nan, 'diff':np.nan}, index = ['left', 'right', 'diff']),
            test_p = pd.Series(data={'t_asy': np.nan, 't_jk':np.nan, 'p_asy':np.nan, 'p_jk':np.nan}, index=['t_asy', 't_jk', 'p_asy', 'p_jk']),
            n = pd.Series(data={'full':n, 'left':nl, 'right':nr, 'eff_left':nlh, 'eff_right':nrh}, index=['full', 'left', 'right', 'eff_left', 'eff_right']),
            h = pd.Series(data={'left':hl, 'right':hr}, index=['left', 'right']),
            fitselect=fitselect, kernel=kernel, vce=vce, c=c, p=p, q=q, useall=useall, bino_flag=bino_flag,
            regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin,
            massPoints=massPoints, massPoints_flag=massPoints_flag, bwselectl=bwselectl, bwselect=bwselect,
            binoN=binoN, binoW=binoW, binoNStep=binoNStep, binoWStep=binoWStep, binoNW= binoNW, binoP=binoP,
            X_min=pd.Series(data={'left': min(np.array(X[X<0].dropna()+c))[0],
                                  'right': min(np.array(X[X>=0].dropna()+c))[0]}, index=['left', 'right']),
            X_max=pd.Series(data={'left': max(np.array(X[X<0].dropna()+c))[0],
                                  'right': max(np.array(X[X>=0].dropna()+c))[0]}, index=['left', 'right']),
            bino = pd.Series(data={'leftN':binomTempLN, 'rightN':binomTempRN, 'leftWindow':binomTempLW, 'rightWindow':binomTempRW, 'pval':binomTempPVal}, index=['leftN', 'rightN', 'leftWindow', 'rightWindow', 'pval'])
        )
    return(result)



class CJMrddensity:
    """
    Class of rddensity function outputs.
    Object type returned by :py:meth:`~rddensity`.
    """

    def __init__(self, hat, sd_asy, sd_jk, test, hat_p, sd_asy_p, sd_jk_p, test_p, n, h,
                 fitselect, kernel, vce, c, p, q, regularize, nLocalMin, bino_flag,
                 nUniqueMin, massPoints, massPoints_flag, bwselectl, bwselect, binoN, binoW,
                 binoNStep, binoWStep, binoNW, binoP, useall,
                 X_min, X_max, bino):
        self.hat = hat
        self.sd_asy = sd_asy
        self.sd_jk = sd_jk
        self.test = test
        self.hat_p = hat_p
        self.sd_asy_p = sd_asy_p
        self.sd_jk_p = sd_jk_p
        self.test_p = test_p
        self.n = n
        self.h = h
        self.fitselect = fitselect
        self.kernel = kernel
        self.vce = vce
        self.c = c
        self.p = p
        self.q = q
        self.useall = useall
        self.bino_flag = bino_flag
        self.regularize = regularize
        self.nLocalMin = nLocalMin
        self.nUniqueMin = nUniqueMin
        self.massPoints = massPoints
        self.massPoints_flag = massPoints_flag
        self.bwselectl = bwselectl
        self.bwselect = bwselect
        self.binoN = binoN
        self.binoW = binoW
        self.binoNStep = binoNStep
        self.binoWStep = binoWStep
        self.binNW = binoNW
        self.binoP = binoP
        self.X_min = X_min
        self.X_max = X_max
        self.bino = bino

    def __str__(self):
        print('Call: rddensity')
        print('')
        fw = 30
        fw_r = 14
        print('Sample Size:'.ljust(fw), str(self.n['full']).rjust(25))
        print('Cutoff:'.ljust(fw), str(self.c).rjust(25))
        print('Model:'.ljust(fw), str(self.fitselect).rjust(25))
        print('Kernel:'.ljust(fw), str(self.kernel).rjust(25))
        print('VCE:'.ljust(fw), str(self.vce).rjust(25))
        return('')

    def __repr__(self):
        print('Manipulation testing using local polynomial density estimation')

        fw = 22
        n_dec = 4
        print('Number of obs:'.ljust(fw), str(self.n['full']).rjust(25))
        print('Model:'.ljust(fw), str(self.fitselect).rjust(25))
        print('Kernel:'.ljust(fw), str(self.kernel).rjust(25))

        if self.bwselectl != 'manual':
            print('BW method:'.ljust(fw), str(self.bwselectl).rjust(25))
        else:
            print('BW method:'.ljust(fw), str(self.bwselectl).rjust(25))

        print('VCE:'.ljust(fw), str(self.vce).rjust(25))

        print('')

        print(str('c = %s'%self.c).ljust(fw), 'Left of c'.rjust(fw), 'Right of c'.rjust(fw))
        print('Number of obs:'.ljust(fw), str(self.n['left']).rjust(fw), str(self.n['right']).rjust(fw))
        print('Eff. number of obs:'.ljust(fw), str(self.n['eff_left']).rjust(fw), str(self.n['eff_right']).rjust(fw))
        print('Order est. (p):'.ljust(fw), str(self.p).rjust(fw), str(self.p).rjust(fw))
        print('Order bias. (q):'.ljust(fw), str(self.q).rjust(fw), str(self.q).rjust(fw))
        print('BW est.'.ljust(fw), str(round(self.h['left'], n_dec)).rjust(fw), str(round(self.h['right'], n_dec)).rjust(fw))

        print('')

        print('Method:'.ljust(fw), 'T'.rjust(fw), 'P > |T|'.rjust(fw))
        if self.useall == True:
            if self.vce == 'plugin':
                print('Conventional'.ljust(fw), str(round(self.test_p['t_asy'], n_dec)).rjust(fw), str(round(self.test_p['p_asy'], n_dec)).rjust(fw))
            else:
                print('Conventional'.ljust(fw), str(round(self.test_p['t_jk'], n_dec)).rjust(fw), str(round(self.test_p['p_jk'], n_dec)).rjust(fw))
        if self.vce == 'plugin':
            print('Robust'.ljust(fw), str(round(self.test['t_asy'], n_dec)).rjust(fw), str(round(self.test['p_asy'], n_dec)).rjust(fw))
        else:
            print('Robust'.ljust(fw), str(round(self.test['t_jk'], n_dec)).rjust(fw), str(round(self.test['p_jk'], n_dec)).rjust(fw))


        print('')

        if self.h['left'] > (self.c - self.X_min['left']):
            warnings.warn('Bandwidth hl greater than the range of the data.')

        if self.h['right'] > (self.X_max['right'] - self.c):
            warnings.warn('Bandwidth hr greater than the range of the data.')

        if self.n['eff_left'] <20 or self.n['eff_right'] <20:
            warnings.warn('Bandwidth may be too small.')

        if self.massPoints_flag:
            warnings.warn("There are repeated observations. Point estimates and standard errors have been adjusted. Use option massPoints=FALSE to suppress this feature.")

        if self.bino_flag == True:
            print('P-values of binomial tests (H0: p = ', self.binoP, ').')
            print('')

            if np.array_equal(self.bino['leftWindow'], self.bino['rightWindow']) == True:
                print('Window Length/2'.ljust(fw), '< c'. rjust(fw), '>= c'.rjust(fw+2), 'P>|T|'.rjust(fw+5))

                for jj in range(len(self.bino['leftWindow'])):
                    print_flag = False
                    if jj ==0:
                        print_flag = True
                    else:
                        if self.bino['leftWindow'][jj] != self.bino['leftWindow'][jj-1] or self.bino['rightWindow'][jj] != self.bino['rightWindow'][jj-1]:
                            print_flag = True
                    if print_flag==True:
                        print()
                        print(str(round(self.bino['leftWindow'][jj], n_dec)).ljust(fw),
                              str(self.bino['leftN'][jj]).rjust(fw), ' ', str(self.bino['rightN'][jj]).rjust(fw), '   ',
                              str(round(self.bino['pval'][jj], n_dec)).rjust(fw))


        return('')
