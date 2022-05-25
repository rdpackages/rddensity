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
    X = pd.DataFrame(X)
    if X.isnull().values.any():
        warnings.warn(X.isnull().sum() + 'missing observation(s) are ignored.\n')

    X = X.dropna(axis=0)
    # sample sizes
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

    X = X-c
    XUnique = XUnique - c
    XUnique = pd.DataFrame(XUnique)
    Xmu = X.mean()[0]
    Xsd = X.std()[0]

    if (nUnique!=n & massPoints)==True:
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

    X = X.iloc[:, 0]
    if regularize==True:
        #bw should not exceed range of data
        bn = min(bn, max(abs(XUnique)[0]))
        cn = min(cn, max(abs(XUnique)[0]))

        #nLocalMin check
        if nLocalMin is not None:
            bn = max(bn,
                     np.sort(abs(X[X<0].reset_index(drop=True)), axis=0)[min(20+p+2, nl)-1],
                     (X[X>= 0].reset_index(drop=True))[min(20 + p+2, nr)-1])
            cn = max(cn,
                     np.sort(abs(X[X<0].reset_index(drop=True)), axis=0)[min(20+p, nl)-1],
                     (X[X>= 0].reset_index(drop=True))[min(20 + p, nr)-1])

        if nUniqueMin is not None:
            bn = max(bn,
                     np.sort(abs(XUnique[XUnique<0].dropna().reset_index(drop=True)), axis=0)[min(20+p+2, nlUnique)-1],
                     (XUnique[XUnique >= 0].dropna().reset_index(drop=True)).loc[min(20 + p+2, nrUnique)-1, 0])
            cn = max(cn,
                     np.sort(abs(XUnique[XUnique<0].dropna().reset_index(drop=True)), axis=0)[min(20+p, nlUnique)-1],
                     (XUnique[XUnique >= 0].dropna().reset_index(drop=True)).loc[min(20 + p, nrUnique)-1, 0])

    #estimate bandwidth
    Y = np.array(range(n))/(n-1)
    Y0 = Y
    if massPoints==True:
        Y = np.repeat(Y[indexUnique], freqUnique)

    Y = pd.DataFrame(Y)
    Yb = pd.DataFrame(Y[abs(X)<=bn])
    Xb = X[abs(X)<=bn].reset_index(drop=True)
    nlb = sum(Xb<0)
    nrb = sum(Xb>=0)
    Yc = pd.DataFrame(Y[abs(X)<=cn])
    Xc = X[abs(X)<=cn].reset_index(drop=True)
    nlc = sum(Xc<0)
    nrc = sum(Xc>=0)

    X = pd.DataFrame(X)
    Xb = pd.DataFrame(Xb)
    Xc = pd.DataFrame(Xc)

    hn = np.zeros((4,3))
    hn = pd.DataFrame(hn, columns=['bw', 'variance', 'biassq'], index=['l', 'r', 'diff', 'sum'])
    fV_b = funs.__rddensity_fv(Y=Yb, X=Xb, nl=nl, nr=nr, nlh=nlb, nrh=nrb, hl=bn, hr=bn, p=p+2, s=p+1, kernel=kernel, fitselect=fitselect, vce=vce, massPoints=massPoints)
    fV_c = funs.__rddensity_fv(Y=Yc, X=Xc, nl=nl, nr=nr, nlh=nlc, nrh=nrc, hl=cn, hr=cn, p=p,   s=1,   kernel=kernel, fitselect=fitselect, vce=vce, massPoints=massPoints)

    if vce=='plugin':
        hn.iloc[:, 1] = n*cn*fV_c.iloc[:, 2]
    else:
        hn.iloc[:, 1] = n*cn*fV_c.iloc[:, 1]

    if fitselect=='unrestricted':
        S = funs.__Sgenerate(p=p, low=0, up=1, kernel=kernel)
        C = funs.__Cgenerate(k=p+1, p=p, low=0, up=1, kernel=kernel)
        hn.iloc[0, 2] = fV_b.iloc[0, 3] * np.matmul(inv(S), C)[1] * np.power(-1, p)
        hn.iloc[1, 2] = fV_b.iloc[1, 3] * np.matmul(inv(S), C)[1]
        hn.iloc[2, 2] = hn.iloc[1, 2] - hn.iloc[0, 2]
        hn.iloc[3, 2] = hn.iloc[1, 2] + hn.iloc[0, 2]
    else:
        Splus = funs.__Splusgenerate(p=p, kernel=kernel)
        Cplus = funs.__Cplusgenerate(k=p+1, p=p, kernel=kernel)
        Psi = funs.__Psigenerate(p=p)
        Sinv = inv(fV_c.iloc[1, 0]*Splus+fV_c.iloc[0,0]*np.matmul(Psi, np.matmul(Splus, Psi)))
        C = fV_b.iloc[0, 3]*(fV_c.iloc[1, 0]*Cplus+ np.power(-1, p+1)*fV_c.iloc[0,0]*np.matmul(Psi, Cplus))
        temp = np.matmul(Sinv, C)
        hn.iloc[0,2] = temp[1]
        hn.iloc[1, 2] = temp[2]
        hn.iloc[2, 2] = hn.iloc[1, 2] - hn.iloc[0, 2]
        hn.iloc[3, 2] = hn.iloc[1, 2] + hn.iloc[0, 2]

    hn.iloc[:, 2] = np.power(hn.iloc[:, 2], 2)
    hn.iloc[:, 0] = np.power(1/(2*p)*hn.iloc[:, 1]/hn.iloc[:, 2]/n, 1/(2*p+1))

    for i in range(4):
        if hn.iloc[i, 1]<0:
            hn.iloc[i, 0] = 0
            hn.iloc[i, 1] = np.nan
        if np.isnan(hn.iloc[i, 0]):
            hn.iloc[i, 0] = 0

    #bw regularization
    if regularize:
        hn.iloc[0,0] = min(hn.iloc[0,0], abs(XUnique.iloc[0, 0]))
        hn.iloc[1,0] = min(hn.iloc[1,0], XUnique.iloc[nUnique-1, 0])
        hn.iloc[2,0] = min(hn.iloc[2,0], max(abs(XUnique.iloc[0, 0]), XUnique.iloc[nUnique-1, 0]))
        hn.iloc[3,0] = min(hn.iloc[3,0], max(abs(XUnique.iloc[0, 0]), XUnique.iloc[nUnique-1, 0]))

        if nLocalMin > 0:
            hlmin = np.sort(abs(X[X<0]))[::-1][min(nl, nLocalMin)-1]
            hrmin = (X[X>= 0].dropna().reset_index(drop=True)).iloc[min(nr, nLocalMin)-1, 0]
            hn.iloc[0,0] = max(hn.iloc[0,0], hlmin)
            hn.iloc[1,0] = max(hn.iloc[1,0], hlmin)
            hn.iloc[2,0] = max(hn.iloc[2,0], hlmin, hrmin)
            hn.iloc[3,0] = max(hn.iloc[3,0], hlmin, hrmin)

        if nUniqueMin>0:
            hlmin = np.sort(abs(XUnique[XUnique<0].dropna().reset_index(drop=True)))[::-1][min(nlUnique, nUniqueMin)-1]
            hrmin = (XUnique[XUnique >= 0].dropna().reset_index(drop=True)).iloc[min(nrUnique, nUniqueMin)-1, 0]
            hn.iloc[0,0] = max(hn.iloc[0,0], hlmin)
            hn.iloc[1,0] = max(hn.iloc[1,0], hlmin)
            hn.iloc[2,0] = max(hn.iloc[2,0], hlmin, hrmin)
            hn.iloc[3,0] = max(hn.iloc[3,0], hlmin, hrmin)

    result = bw_output(hn, n= pd.Series(data={'full':n, 'left':nl, 'right':nr}, index=['full', 'left', 'right']),
                       fitselect=fitselect, kernel=kernel, vce=vce, c=c, p=p,
                          regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin,
                          massPoints=massPoints, massPoints_flag=massPoints_flag,
                       X_min=pd.Series(data={'left': min(np.array(X[X<0].dropna()+c))[0],
                                             'right': min(np.array(X[X>=0].dropna()+c))[0]}, index=['left', 'right']),
                       X_max=pd.Series(data={'left': max(np.array(X[X<0].dropna()+c))[0],
                                             'right': max(np.array(X[X>=0].dropna()+c))[0]}, index=['left', 'right'])
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
