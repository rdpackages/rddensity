********************************************************************************
** RDDENSITY Stata Package 
** Do-file for Empirical Illustration
** Authors: Matias D. Cattaneo, Michael Jansson and Xinwei Ma
********************************************************************************
** net install rddensity, from(https://raw.githubusercontent.com/rdpackages/rddensity/master/stata) replace
** net install lpdensity, from(https://raw.githubusercontent.com/nppackages/lpdensity/master/stata) replace
********************************************************************************
clear all
set more off 
 
********************************************************************************
** Load data 
********************************************************************************
use "rddensity_senate.dta", clear
sum margin

********************************************************************************
** rddensity: default options
********************************************************************************
rddensity margin

********************************************************************************
** rddensity: with plot
********************************************************************************
rddensity margin, plot 
rddensity margin, plot plot_range(-50 50) hist_range(-50 50)

********************************************************************************
** rddensity: all statistics & default options
********************************************************************************
rddensity margin, all

********************************************************************************
** rddensity: default statistic, restricted model & plugin standard errors
********************************************************************************
rddensity margin, fitselect(restricted) vce(plugin)

********************************************************************************
** rdbwdensity: default options
********************************************************************************
rdbwdensity margin

********************************************************************************
** rdbwdensity: compute bandwidth and then use for rddensity
********************************************************************************
qui rdbwdensity margin
mat h = e(h)
local hr = h[2,1]
rddensity margin, h(10 `hr')

********************************************************************************
** Other examples
********************************************************************************
rddensity margin, kernel(uniform)
rddensity margin, bwselect(diff)
rddensity margin, h(10 15)
rddensity margin, p(2) q(4)
rddensity margin, c(5) all

rdbwdensity margin, p(3) fitselect(restricted)
rdbwdensity margin, kernel(uniform) vce(jackknife)
