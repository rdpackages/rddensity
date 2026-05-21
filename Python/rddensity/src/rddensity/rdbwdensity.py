#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#ignore future warnings from pandas
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
from scipy.stats import norm
import math
import scipy.integrate as integrate
from numpy.linalg import inv
import scipy.optimize as optimize
from scipy.stats import norm
import statistics as stat
import pandas as pd

from . import funs

def rdbwdensity(X, c=0, p=2,
                fitselect='unrestricted', kernel='triangular', vce='jackknife',
                massPoints=True, regularize=True, nLocalMin=None, nUniqueMin=None):

    """

	Parameters
	----------
    X: Numeric vector or one dimentional matrix/dataframe
       the running variable.
    c: Numeric
       specifies the threshold or cutoff value in the support of *X*. Default is *0*.
    p: Nonnegative integer,
        specifies the local polynomial order used to construct the density estimators. Default is *2* (local quadratic approximation).
    fitselect: String
        specifies the density estimation method. *unrestricted* (Default) for density estimation without any restrictions (two-sample, unrestricted inference). *restricted* for density estimation assuming equal distribution function and higher order dericatives.
    kernel: String
        specifies the kernel function used to construct the local polynomial estimators. Accepted kernels: *triangular* (Default), *epanechnikov* or *uniform*.
    vce: String
        specifies the procedure used to compute the variance-covariance matrix estimatior. *jackknife* (Default) for jackknife standard errors or *plugin* for asymptotic plug-in standard errors.
    massPoints: Boolean, Default *True*.
        Specifies wether to adjust for mass points in the data.
    regularize: Boolean, Default *True*.
        Specifies whether to conduct local sample size checking. When True, the bandwidth is chosen such that the local region includes at least *nLocalMin* observations and at least *nUniqueMin* unique observations.
    nLocalMin: Nonnegative integer
        specifies the minimum number of observations in each local neighbourhood. This option will be ignored if set to *0* or if *regularize=False*. Default is *20+p+1*.
    nUniqueMin: Nonnegative integer
        specifies the minimum number of unqieu observations in each local neighbourhood. This option will be ignored if set to *0* or if *regularize=False*. Default is *20+p+1*.

    Returns
    -------

    h
        Bandwidths for density discontinuity test, left and right of the cutoff, asymptotic variance and bias.
    n
        full-full sample size, *left*/*right*: sample size to the left/right of the cutoff.
    X_min
        Smallest observations to the left and right of the cutoff.
    X_max
        Largest observations to the left and right of the cutoff.
    options
        other options passed to the function are also stored within the object.

    See Also
    --------
    rddensity.rddensity
    rddensity.rdplotdensity

    """

    # missing values
    X = funs._as_numeric_vector(X)
    missing = np.isnan(X)
    if np.any(missing):
        warnings.warn(f'{int(np.sum(missing))} missing observation(s) are ignored.\n')
        X = X[~missing]

    if len(X) == 0:
        raise Exception("X cannot be empty after removing missing values.")

    # sample sizes
    X = np.sort(X)
    n = len(X)
    nl = int(np.searchsorted(X, c, side="left"))
    nr = n - nl
    Xmin = X[0]
    Xmax = X[-1]
    XUnique = funs.__rddensityUnique(X)
    freqUnique = np.asarray(XUnique["freq"], dtype=int)
    indexUnique = np.asarray(XUnique["indexLast"], dtype=int)
    XUnique = np.asarray(XUnique["unique"], dtype=float)
    nUnique = len(XUnique)
    nlUnique = int(np.searchsorted(XUnique, c, side="left"))
    nrUnique = nUnique - nlUnique

    X = X-c
    XUnique = XUnique - c
    Xmu = float(np.mean(X))
    Xsd = float(np.std(X, ddof=1))

    if nUnique != n and massPoints:
        massPoints_flag = True
    else:
        massPoints_flag = False
    #Error Handling
    if c<Xmin or c>Xmax:
        raise Exception("The cutoff should be set within the range of the data.")
    if p<1 or p>7:
        raise Exception("p must be an integer between 1 and 7")
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

    # selecting preliminary bandwidth
    fhatb = 1/(np.power(funs.__rddensity_H(Xmu/Xsd, p+2), 2) *norm.pdf(Xmu/Xsd))
    fhatc = 1/(np.power(funs.__rddensity_H(Xmu/Xsd, p), 2) *norm.pdf(Xmu/Xsd))

    Cb = [25884.444444494150957,3430865.4551236177795,845007948.04262602329,330631733667.03808594,187774809656037.3125,145729502641999264,146013502974449876992]
    Cc = [4.8000000000000246914,548.57142857155463389,100800.00000020420703,29558225.458100609481,12896196859.612621307,7890871468221.609375,6467911284037581]
    # bn is for higher-order derivative estimation
    bn = np.power((2*p+1)/4 * fhatb*Cb[p-1]/n, 1/(2*p+5))
    bn = bn * Xsd
    # cn is for density estimation
    cn = np.power(1/(2*p)*fhatc * Cc[p-1]/n, 1/(2*p+1))
    cn = cn*Xsd

    if regularize==True:
        #bw should not exceed range of data
        bn = min(bn, np.max(np.abs(XUnique)))
        cn = min(cn, np.max(np.abs(XUnique)))

        #nLocalMin check
        if nLocalMin is not None:
            bn = max(bn,
                     np.sort(np.abs(X[X<0]))[min(20+p+2, nl)-1],
                     X[X>= 0][min(20 + p+2, nr)-1])
            cn = max(cn,
                     np.sort(np.abs(X[X<0]))[min(20+p, nl)-1],
                     X[X>= 0][min(20 + p, nr)-1])

        if nUniqueMin is not None:
            bn = max(bn,
                     np.sort(np.abs(XUnique[XUnique<0]))[min(20+p+2, nlUnique)-1],
                     XUnique[XUnique >= 0][min(20 + p+2, nrUnique)-1])
            cn = max(cn,
                     np.sort(np.abs(XUnique[XUnique<0]))[min(20+p, nlUnique)-1],
                     XUnique[XUnique >= 0][min(20 + p, nrUnique)-1])

    #estimate bandwidth
    Y = np.array(range(n))/(n-1)
    Y0 = Y
    if massPoints==True:
        Y = np.repeat(Y[indexUnique], freqUnique)

    mask_b = np.abs(X) <= bn
    Yb = Y[mask_b]
    Xb = X[mask_b]
    nlb = int(np.sum(Xb < 0))
    nrb = len(Xb) - nlb
    mask_c = np.abs(X) <= cn
    Yc = Y[mask_c]
    Xc = X[mask_c]
    nlc = int(np.sum(Xc < 0))
    nrc = len(Xc) - nlc

    fV_b = funs.__rddensity_fv(Y=Yb, X=Xb, nl=nl, nr=nr, nlh=nlb, nrh=nrb, hl=bn, hr=bn, p=p+2, s=p+1, kernel=kernel, fitselect=fitselect, vce=vce, massPoints=massPoints)
    fV_c = funs.__rddensity_fv(Y=Yc, X=Xc, nl=nl, nr=nr, nlh=nlc, nrh=nrc, hl=cn, hr=cn, p=p,   s=1,   kernel=kernel, fitselect=fitselect, vce=vce, massPoints=massPoints)
    fV_b_values = fV_b.to_numpy(dtype=float)
    fV_c_values = fV_c.to_numpy(dtype=float)

    hn = np.zeros((4,3))
    if vce=='plugin':
        hn[:, 1] = n*cn*fV_c_values[:, 2]
    else:
        hn[:, 1] = n*cn*fV_c_values[:, 1]

    if fitselect=='unrestricted':
        S = funs.__Sgenerate(p=p, low=0, up=1, kernel=kernel)
        C = funs.__Cgenerate(k=p+1, p=p, low=0, up=1, kernel=kernel)
        bias_constant = np.matmul(inv(S), C)[1, 0]
        hn[0, 2] = fV_b_values[0, 3] * bias_constant * np.power(-1, p)
        hn[1, 2] = fV_b_values[1, 3] * bias_constant
        hn[2, 2] = hn[1, 2] - hn[0, 2]
        hn[3, 2] = hn[1, 2] + hn[0, 2]
    else:
        Splus = funs.__Splusgenerate(p=p, kernel=kernel)
        Cplus = funs.__Cplusgenerate(k=p+1, p=p, kernel=kernel)
        Psi = funs.__Psigenerate(p=p)
        Sinv = inv(fV_c_values[1, 0]*Splus+fV_c_values[0,0]*np.matmul(Psi, np.matmul(Splus, Psi)))
        C = fV_b_values[0, 3]*(fV_c_values[1, 0]*Cplus+ np.power(-1, p+1)*fV_c_values[0,0]*np.matmul(Psi, Cplus))
        temp = np.matmul(Sinv, C)
        hn[0,2] = temp[1, 0]
        hn[1, 2] = temp[2, 0]
        hn[2, 2] = hn[1, 2] - hn[0, 2]
        hn[3, 2] = hn[1, 2] + hn[0, 2]

    hn[:, 2] = np.power(hn[:, 2], 2)
    hn[:, 0] = np.power(1/(2*p)*hn[:, 1]/hn[:, 2]/n, 1/(2*p+1))

    for i in range(4):
        if hn[i, 1]<0:
            hn[i, 0] = 0
            hn[i, 1] = np.nan
        if np.isnan(hn[i, 0]):
            hn[i, 0] = 0

    #bw regularization
    if regularize:
        unique_left_range = abs(XUnique[0])
        unique_right_range = XUnique[nUnique-1]
        unique_max_range = max(unique_left_range, unique_right_range)
        hn[0,0] = min(hn[0,0], unique_left_range)
        hn[1,0] = min(hn[1,0], unique_right_range)
        hn[2,0] = min(hn[2,0], unique_max_range)
        hn[3,0] = min(hn[3,0], unique_max_range)

        if nLocalMin > 0:
            hlmin = np.sort(np.abs(X[X<0]))[::-1][min(nl, nLocalMin)-1]
            hrmin = X[X>= 0][min(nr, nLocalMin)-1]
            hn[0,0] = max(hn[0,0], hlmin)
            hn[1,0] = max(hn[1,0], hlmin)
            hn[2,0] = max(hn[2,0], hlmin, hrmin)
            hn[3,0] = max(hn[3,0], hlmin, hrmin)

        if nUniqueMin>0:
            hlmin = np.sort(np.abs(XUnique[XUnique<0]))[::-1][min(nlUnique, nUniqueMin)-1]
            hrmin = XUnique[XUnique >= 0][min(nrUnique, nUniqueMin)-1]
            hn[0,0] = max(hn[0,0], hlmin)
            hn[1,0] = max(hn[1,0], hlmin)
            hn[2,0] = max(hn[2,0], hlmin, hrmin)
            hn[3,0] = max(hn[3,0], hlmin, hrmin)

    X_left = X[X < 0] + c
    X_right = X[X >= 0] + c
    hn = pd.DataFrame(hn, columns=['bw', 'variance', 'biassq'], index=['l', 'r', 'diff', 'sum'])
    result = bw_output(hn, n= pd.Series(data={'full':n, 'left':nl, 'right':nr}, index=['full', 'left', 'right']),
                       fitselect=fitselect, kernel=kernel, vce=vce, c=c, p=p,
                          regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin,
                          massPoints=massPoints, massPoints_flag=massPoints_flag,
                       X_min=pd.Series(data={'left': np.min(X_left),
                                             'right': np.min(X_right)}, index=['left', 'right']),
                       X_max=pd.Series(data={'left': np.max(X_left),
                                             'right': np.max(X_right)}, index=['left', 'right'])
                       )


    return(result)

class bw_output:
    """
    Class of rdbwdensity function outputs.
    Object type returned by :py:meth:`~rdbwdensity`.
    """
    def __init__(self, h, n, fitselect, kernel, vce, c, p,
                 regularize, nLocalMin, nUniqueMin, massPoints, massPoints_flag, X_min, X_max):
        self.h = h
        self.n = n
        self.fitselect = fitselect
        self.kernel = kernel
        self.vce = vce
        self.c = c
        self.p = p
        self.regularize = regularize
        self.nLocalMin = nLocalMin
        self.nUniqueMin = nUniqueMin
        self.massPoints = massPoints
        self.massPoints_flag = massPoints_flag
        self.X_min = X_min
        self.X_max = X_max

    def __str__(self):
        print('Call: rdbwdensity')
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
        print('Bandwidth selection for manupulation testing')

        if self.massPoints_flag:
            warnings.warn("There are repeated observations. Point estimates and standard errors have been adjusted. Use option massPoints=FALSE to suppress this feature.")

        fw = 20
        n_dec = 4
        print('Number of obs:'.ljust(fw), str(self.n['full']).rjust(25))
        print('Model:'.ljust(fw), str(self.fitselect).rjust(25))
        print('Kernel:'.ljust(fw), str(self.kernel).rjust(25))
        print('VCE:'.ljust(fw), str(self.vce).rjust(25))

        print('')
        print('')

        print('Cutoff:'.ljust(fw), 'Left of c'.rjust(fw), 'Right of c'.rjust(fw))
        print('Number of obs:'.ljust(fw), str(self.n['left']).rjust(fw), str(self.n['right']).rjust(fw))
        print('Min Running var:'.ljust(fw), str(round(self.X_min['left'], n_dec)).rjust(fw), str(round(self.X_min['right'], n_dec)).rjust(fw))
        print('Max Running var:'.ljust(fw), str(round(self.X_max['left'], n_dec)).rjust(fw), str(round(self.X_max['right'], n_dec)).rjust(fw))
        print('Order est. (p):'.ljust(fw), str(self.p).rjust(fw), str(self.p).rjust(fw))

        print('')
        print('')

        print('Target'.ljust(fw), 'Bandwidth'.rjust(fw), 'Variance'.rjust(fw), 'Bias^2'.rjust(fw))
        print('left density'.ljust(fw),
              str(round(self.h.iloc[0,0], n_dec)).rjust(fw),
              str(round(self.h.iloc[0,1], n_dec)).rjust(fw),
              str(round(self.h.iloc[0,2], n_dec)).rjust(fw),)
        print('right density'.ljust(fw),
              str(round(self.h.iloc[1,0], n_dec)).rjust(fw),
              str(round(self.h.iloc[1,1], n_dec)).rjust(fw),
              str(round(self.h.iloc[1,2], n_dec)).rjust(fw),)
        print('diff. densities'.ljust(fw),
              str(round(self.h.iloc[2,0], n_dec)).rjust(fw),
              str(round(self.h.iloc[2,1], n_dec)).rjust(fw),
              str(round(self.h.iloc[2,2], n_dec)).rjust(fw),)
        print('sum densities'.ljust(fw),
              str(round(self.h.iloc[3,0], n_dec)).rjust(fw),
              str(round(self.h.iloc[3,1], n_dec)).rjust(fw),
              str(round(self.h.iloc[3,2], n_dec)).rjust(fw),)
        return('')
