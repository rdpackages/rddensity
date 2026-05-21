#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#ignore future warnings from pandas
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

from functools import lru_cache
import numpy as np
from scipy.stats import norm
import sympy
from sympy.abc import x, y
import math
import scipy.integrate as integrate
from numpy.linalg import inv
import scipy.optimize as optimize
import statistics as stat
import pandas as pd


def _as_numeric_vector(value, name="X"):
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        raise Exception(f"{name} should be numeric.") from None

    if array.ndim == 0:
        array = array.reshape(1)
    elif array.ndim == 1:
        array = array.reshape(-1)
    elif array.ndim == 2 and 1 in array.shape:
        array = array.reshape(-1)
    else:
        raise Exception(f"{name} should be one-dimensional.")

    if len(array) == 0:
        raise Exception(f"{name} cannot be empty.")
    return array


def __rddensityUnique(x):
# Find unique elements and their frequencies in a numeric vector. This function
#  has a similar performance as the built-in R function \code{unique}.
    values = _as_numeric_vector(x)
    n = len(values)
    #if x has no elements
    if n==0:
        return({'unique':np.array([], dtype=float), 'freq':np.array([], dtype=int), 'indexFirst':np.array([], dtype=int), 'indexLast':np.array([], dtype=int)})
    if n==1:
        return({'unique':values.copy(), 'freq': np.array([1]), 'indexFirst':np.array([0]), 'indexLast':np.array([0])})

    unique_mask = np.r_[values[1:] != values[:-1], True]
    unique = values[unique_mask]
    index_last = np.flatnonzero(unique_mask)
    freq = np.diff(np.r_[-1, index_last])
    index_first = index_last - freq + 1

    return({'unique':unique, 'freq' : freq, 'indexFirst':index_first, 'indexLast': index_last})


# Generate S matrix
# Output matches R
def __Sgenerate(p, low=-1, up=1, kernel='triangular'):
    return _Sgenerate_cached(p, float(low), float(up), kernel).copy()


@lru_cache(maxsize=None)
def _Sgenerate_cached(p, low, up, kernel):
    S = np.zeros((p+1, p+1))
    for i in range(1,p+2):
        for j in range(1,p+2):
            if kernel=='uniform':
                integrand = lambda l: 0.5*np.power(l, i+j-2)
            elif kernel=='epanechnikov':
                integrand = lambda l: 0.75*np.power(l, i+j-2)*(1-l**2)
            else:
                integrand = lambda l: (1-abs(l))*np.power(l, i+j-2)

            S[i-1, j-1] = integrate.quad(integrand, low, up)[0]

    return(S)


def __Sminusgenerate(p, kernel='triangular'):
    S = __Sgenerate(p, low=-1, up=0, kernel=kernel)
    temp = np.zeros((p+2, p+2))
    temp[0:1, 0:1] = S[0:1, 0:1]
    if p>1:
        # insert row of zeros
        temp = np.insert(S, 2, np.zeros((1, p+1)), 0)
        # insert column of zeros
        temp = np.insert(temp, 2, np.zeros((1, p+2)), 1)
    return(temp)

def __Splusgenerate(p, kernel='triangular'):
    S = __Sgenerate(p, low=0, up=1, kernel=kernel)
    # insert row of zeros
    temp = np.insert(S, 1, np.zeros((1, p+1)), 0)
    # insert column of zeros
    temp = np.insert(temp, 1, np.zeros((1, p+2)), 1)
    return(temp)

# Generate C matrix
def __Cgenerate(k, p, low=-1, up=1, kernel='triangular'):
    return _Cgenerate_cached(k, p, float(low), float(up), kernel).copy()


@lru_cache(maxsize=None)
def _Cgenerate_cached(k, p, low, up, kernel):
    C = np.zeros((p+1, 1))
    for i in range(1, p+2):
        if kernel=='uniform':
            integrand = lambda l: 0.5*np.power(l, i+k-1)
        elif kernel=='epanechnikov':
            integrand = lambda l: 0.75*np.power(l, i+k-1)*(1-l**2)
        else:
            integrand = lambda l: (1-abs(l))*np.power(l, i+k-1)

        C[i-1, 0] = integrate.quad(integrand, low, up)[0]

    return(C)


def __Cminusgenerate(k, p, kernel='triangular'):
    C = __Cgenerate(k=k, p=p, kernel=kernel, low=-1, up=0)
    temp = np.zeros((p+2, 1))
    temp[0:1] = C[0:1]
    if p>1:
        # insert row of zeros
        temp = np.insert(C, 2, np.zeros((1, 1)), 0)
    return(temp)

def __Cplusgenerate(k, p, kernel='triangular'):
    C = __Cgenerate(k=k, p=p, kernel=kernel, low=0, up=1)
    temp = np.zeros((p+2, 1))
    # insert row of zeros
    temp = np.insert(C, 1, np.zeros((1, 1)), 0)

    return(temp)

# Generate G matrix
def __Ggenerate(p, low=-1, up=1, kernel='triangular'):
    return _Ggenerate_cached(p, float(low), float(up), kernel).copy()


@lru_cache(maxsize=None)
def _Ggenerate_cached(p, low, up, kernel):
    G = np.zeros((p+1, p+1))
    for i in range(1, p+2):
        for j in range(1, p+2):
            if kernel=='uniform':
                def integrand_1(x,y):
                    return(np.power(x,i) * np.power(y, (j-1))*0.25)

                def integrand_2(x,y):
                    return(np.power(x,i-1) * np.power(y, j)*0.25)

                def x_integral_1(y):
                    return (integrate.quad(integrand_1, low, y, args=(y))[0])

                def x_integral_2(y):
                    return (integrate.quad(integrand_2, y, up, args=(y))[0])

                G[i-1, j-1] = integrate.quad(x_integral_1, low, up)[0] + integrate.quad(x_integral_2, low, up)[0]

            elif kernel == 'epanechnikov':
                def integrand_1(x,y):
                    return(np.power(x,i) * np.power(y, (j-1))*np.power(0.75, 2)*(1-x**2)*(1-y**2))

                def integrand_2(x,y):
                    return(np.power(x,i-1) * np.power(y, j)*np.power(0.75, 2)*(1-x**2)*(1-y**2))

                def x_integral_1(y):
                    return (integrate.quad(integrand_1, low, y, args=(y))[0])

                def x_integral_2(y):
                    return (integrate.quad(integrand_2, y, up, args=(y))[0])

                G[i-1, j-1] = integrate.quad(x_integral_1, low, up)[0] + integrate.quad(x_integral_2, low, up)[0]

            elif kernel=='triangular':
                def integrand_1(x,y):
                    return(np.power(x,i) * np.power(y, (j-1))*(1-abs(x))*(1-abs(y)))

                def integrand_2(x,y):
                    return(np.power(x,i-1) * np.power(y, j)*(1-abs(x))*(1-abs(y)))

                def x_integral_1(y):
                    return (integrate.quad(integrand_1, low, y, args=(y))[0])

                def x_integral_2(y):
                    return (integrate.quad(integrand_2, y, up, args=(y))[0])

                G[i-1, j-1] = integrate.quad(x_integral_1, low, up)[0] + integrate.quad(x_integral_2, low, up)[0]

    return(G)


def __Gminusgenerate(p, kernel='triangular'):
    G = __Ggenerate(p, low=-1, up=0, kernel=kernel)
    temp = np.zeros((p+2, p+2))
    temp[0:1, 0:1] = G[0:1, 0:1]
    if p>1:
        # insert row of zeros
        temp = np.insert(G, 2, np.zeros((1, p+1)), 0)
        # insert column of zeros
        temp = np.insert(temp, 2, np.zeros((1, p+2)), 1)
    return(temp)

def __Gplusgenerate(p, kernel='triangular'):
    G = __Ggenerate(p, low=0, up=1, kernel=kernel)
    # insert row of zeros
    temp = np.insert(G, 1, np.zeros((1, p+1)), 0)
    # insert column of zeros
    temp = np.insert(temp, 1, np.zeros((1, p+2)), 1)
    return(temp)

def __Psigenerate(p):
    if p>1:
        temp = np.append(np.array([1, 0, 0]), np.power(-1, range(2, p+1)))
    else:
        temp = np.array([1,0,0])

    temp = np.diag(temp)
    temp[1,2] = temp[2,1] = -1
    return(temp)


def __rddensity_fv(Y, X, nl, nr, nlh, nrh, hl, hr, p, s, kernel, fitselect, vce, massPoints):
    n = nl + nr
    nh = nlh +nrh
    Y = np.asarray(Y, dtype=float).reshape(-1)
    X = np.asarray(X, dtype=float).reshape(-1)


    # Construct kernel weights
    W = np.repeat(np.nan, nh)
    if kernel=='uniform':
        W[0:nlh] = 1/(2*hl)
        W[nlh:nh] = 1/(2*hr)
    elif kernel=='triangular':
        W[0:nlh] = (1 + X[0:nlh]/hl)/hl
        W[nlh:nh] = (1-X[nlh:nh]/hr)/hr
    else:
        W[0:nlh] = 0.75*((1 + X[0:nlh]/hl)**2)/hl
        W[nlh:nh] = 0.75*((1 + X[nlh:nh]/hr)**2)/hr


    # Construct the design matrix and the bandwidth matrix
    if fitselect=='restricted':
        if nh ==0:
            Xp=None
        else:
            Xp = np.zeros((nh, p+2))
            Xp[:, 0] = 1
            Xp[0:nlh, 1] = X[0:nlh]/hl
            Xp[nlh:nh, 1] = 0
            Xp[0:nlh, 2] = 0
            Xp[nlh:nh, 2] = X[nlh:nh]/hr

            if p>1:
                for j in range(3, p+2):
                    Xp[0:nlh, j] = np.power(X[0:nlh]/hl, j-1)
                    Xp[nlh:nh, j] = np.power(X[nlh:nh]/hr, j-1)
                    v = np.append(np.array([0, 1, 1]), range(2, p+1))
            else:
                v = np.array([0, 1, 1])
            HpInv = np.diag(1 / np.power(hl, v))
    else:
        if nh ==0:
            Xp=None
        else:
            Xp = np.zeros((nh, 2*p+2))
            Hp = np.repeat(0.0, 2*p+2)
            for power in range(p+1):
                Xp[0:nlh, 2*power] = np.power(X[0:nlh]/hl, power)
                Xp[nlh:nh, 2*power] = 0
                Xp[0:nlh, 2*power+1] = 0
                Xp[nlh:nh, 2*power+1] = np.power(X[nlh:nh]/hr, power)
                Hp[2*power] = np.power(hl, power)
                Hp[2*power+1] = np.power(hr, power)
            HpInv = np.diag(1 / Hp)


    out = np.full((4,4), np.nan)
    out_index = ('l', 'r', 'diff', 'sum')
    out_columns = ('hat', 'jackknife', 'plugin', 's')

    #X'WX inverse matrix
    nrow, ncol = Xp.shape
    XpW = Xp * W[:, None]

    try:
        Sinv = inv(XpW.T @ Xp)
    except:
        return pd.DataFrame(out, columns=out_columns, index=out_index)

    # point estimates
    b = (XpW.T @ Y) @ Sinv @ HpInv

    if fitselect=='restricted':
        out[0, 0] = b[1]
        out[1, 0] = b[2]
        out[2, 0] = b[2] - b[1]
        out[3, 0] = b[2] + b[1]
        out[0, 3] = b[s+1]
        out[1, 3] = b[s+1]
        out[2, 3] = 0
        out[3, 3] = 2 * out[0,3]
    else:
        out[0, 0] = b[2]
        out[1, 0] = b[3]
        out[2, 0] = b[3] - b[2]
        out[3, 0] = b[3] + b[2]
        out[0, 3] = b[2*s]
        out[1, 3] = b[2*s+1]
        out[2, 3] = out[1, 3] - out[0, 3]
        out[3, 3] = out[1, 3] + out[0, 3]

    # Jackknife
    if vce=='jackknife':
        L = np.zeros((nrow, ncol))
        if massPoints==True:
            _, indexunique, frequnique = np.unique(X, return_index=True, return_counts=True)
            for jj in range(ncol):
                xpw_col = np.r_[XpW[:, jj], 0.0]
                L[:, jj] = np.repeat((np.cumsum(xpw_col[::-1])/(n-1))[indexunique], frequnique)[::-1]
        else:
            L[0, :] = np.sum(XpW[1:, :], axis=0)/(n-1)
            for i in range(1, nh):
                L[i, :] = L[i-1, :] - XpW[i, :]/(n-1)

        V = HpInv @ Sinv @ (L.T @ L) @ Sinv @ HpInv

        if fitselect=='restricted':
            out[0, 1] = V[1, 1]
            out[1, 1] = V[2, 2]
            out[2, 1] = V[1, 1] + V[2, 2] - 2*V[1, 2]
            out[3, 1] = V[1, 1] + V[2, 2] + 2*V[1, 2]
        else:
            out[0, 1] = V[2, 2]
            out[1, 1] = V[3, 3]
            out[2, 1] = V[2, 2] + V[3, 3] - 2*V[2, 3]
            out[3, 1] = V[2, 2] + V[3, 3] + 2*V[2, 3]

    # plug in
    if vce=='plugin':
        if fitselect=='unrestricted':
            S = __Sgenerate(p, low=0, up=1, kernel=kernel)
            G = __Ggenerate(p, low=0, up=1, kernel=kernel)
            Sinv_kernel = inv(S)
            V = Sinv_kernel @ G @ Sinv_kernel
            out[0, 2] = out[0, 0] * V[1, 1]/(n*hl)
            out[1, 2] = out[1, 0] * V[1, 1]/(n*hr)
            out[2, 2] = out[0, 2] + out[1, 2]
            out[3, 2] = out[0, 2] + out[1, 2]
        else:
            S = __Splusgenerate(p=p, kernel=kernel)
            G = __Gplusgenerate(p=p, kernel=kernel)
            Psi = __Psigenerate(p=p)
            Sm = np.matmul(Psi, np.matmul(S, Psi))
            Gm = np.matmul(Psi, np.matmul(G, Psi))

            Ainv = inv(out[0, 0] * Sm + out[1, 0]*S)
            V = Ainv @ (np.power(out[0, 0], 3)*Gm + np.power(out[1, 0], 3)*G) @ Ainv
            out[0, 2] = V[1, 1]/(n*hl)
            out[1, 2] = V[2, 2]/(n*hl)
            out[2, 2] = (V[1, 1] + V[2, 2] - 2*V[1, 2])/(n*hl)
            out[3, 2] = (V[1, 1] + V[2, 2] + 2*V[1, 2])/(n*hl)

    for i in range(4):
        for j in range(1, 3):
            if not np.isnan(out[i, j]):
                if out[i, j]<0:
                    out[i, j] = np.nan

    out = pd.DataFrame(out, columns=out_columns, index=out_index)
    return(out)


def __h_opt_density(x, p, n, dgp_f1, dgp_fp1, f_low, f_up, kernel='triangular'):
    if kernel != 'triangular' & kernel!='epanechnikov' & kernel != 'uniform':
        warnings.warn("invalid kernel input given. Using triangular kernel instead.")
        kernel = 'triangular'

    if x==f_low | x==f_up:
        c_low = 0
    else:
        c_low = -1

    c_up = 1

    e_mat = np.zeros((p+1, 1))
    e_mat[1] = 1

    S1 = __Sgenerate(p=p, low=c_low, up=c_up, kernel=kernel)
    Cp1 = __Cgenerate(k=p+1, p=p, low = c_low, up=c_up, kernel=kernel)
    G =__Ggenerate(p=p, low=c_low, up=c_up, kernel=kernel)

    kappa = np.power(n, -1/(2*p+1)) * np.power(np.matmul(np.transpose(e_mat), np.matmul(inv(S1), np.matmul(G, np.matmul(inv(S1), e_mat)))), 1/(2*p+1)) * np.power(np.matmul(abs(np.transpose(e_mat)), np.matmul(inv(S1), Cp1)), -2/(2*p+1)) * np.power(math.factorial(p+1), 2/(2*p+1)) * np.power(2*p, -1/(2*p+1))

    biassq = np.power(abs(dgp_fp1), -2/(2*p+1))
    var1 = np.power(dgp_f1, 1/(2*p+1))

    h_opt = kappa * biassq * var1

    return(h_opt)



def __h_opt_density_res(p, n, dgp_f1_l, dgp_f1_r, dgp_fp1_l, dgp_fp1_r, kernel='triangular'):
    if kernel != 'triangular' & kernel!='epanechnikov' & kernel != 'uniform':
        warnings.warn("invalid kernel input given. Using triangular kernel instead.")
        kernel = 'triangular'

    e_mat = np.zeros((p+2, 1))
    if p%2:
        e[1] = -1
        e[2] = -1
    else:
        e[1] = 1
        e[2] = 1

    Sminus = __Sminusgenerate(p=p, kernel=kernel)
    Splus = __Splusgenerate(p=p, kernel=kernel)
    Sinv = inv(Sminus*dgp_f1_l+ Splus*dgp_f1_r)
    C = dgp_f1_l * dgp_f1_l * __Cminusgenerate(k=p+1, p=p, kernel=kernel) + dgp_f1_r * dgp_f1_r * __Cplusgenerate(k=p+1, p=p, kernel=kernel)
    G = np.power(dgp_f1_l, 3) * __Gminusgenerate(p=p, kernel=kernel) + np.power(dgp_f1_r, 3) * __Gplusgenerate(p=p, kernel=kernel) + np.power(dgp_f1_l, 2) * dgp_f1_r* np.matmul(Sminus[1, :], Splus[0, :]) + np.power(dgp_f1_l, 2) * dgp_f1_r* np.matmul(Splus[0, :], Sminus[1, :])

    h_opt = np.power(n, -1/(2*p+1)) * np.power(np.matmul(e_mat, Sinv) * np.matmul(G, np.matmul(Sinv, e_mat)), 1/(2*p+1))* np.power(np.matmul(abs(e_mat), np.matmul(Sinv, C)), -2/(2*p+1)) * np.power(math.factorial(p+1), 2/(2*p+1)) * np.power(2*p, -1/(2*p+1))

    return(h_opt)

def __rddensity_H(x, p):
    if p==0:
        out = 1
    if p==1:
        out = x
    if p==2:
        out = x**2 -1
    if p==3:
        out = np.power(x, 3) - 3*x
    if p==4:
        out = np.power(x, 4) - 6 * np.power(x, 2) + 3
    if p==5:
        out = np.power(x, 5) - 10 * np.power(x, 3) + 15*x
    if p==6:
        out = np.power(x, 6) - 15 * np.power(x, 4) + 45*np.power(x, 2) - 15
    if p==7:
        out = np.power(x, 7) - 21 * np.power(x, 5) + 105*np.power(x, 3) - 105*x
    if p==8:
        out = np.power(x, 8) - 28*np.power(x, 6)  + 210*np.power(x, 4) - 420*np.power(x, 2) + 105
    if p==9:
        out = np.power(x, 9) - 36*np.power(x, 7) + 378*np.power(x, 5) - 1260*np.power(x, 3) + 945*x
    if p==10:
        out = np.power(x, 10) - 45*np.power(x, 8) + 630*np.power(x, 6) - 3150*np.power(x, 4) + 4725*np.power(x, 2) - 945
    return(out)
