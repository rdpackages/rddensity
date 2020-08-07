###########################################################################
## RDDENSITY R Package
## Do-file for RDDENSITY Illustration
## Authors: Matias D. Cattaneo, Michael Jansson and Xinwei Ma
###########################################################################
### Clear R environment
rm(list=ls(all=TRUE))
setwd("...")

### Install R library
### NOTE: depending on your system, you may need to do it as root
#install.packages('rddensity')

### Load rddensity package
library(rddensity)

### Load data base
data("rddensity_senate")

### Summary stats
summary(margin)

### rddensity: default options
summary(rddensity(X = margin))

### rddensity: all statistics & default options
summary(rddensity(X = margin, all=TRUE))

### rddensity: default statistic, restricted model & plugin standard errors
summary(rddensity(X = margin, fitselect="restricted", vce="plugin"))

### rdplotdensity: default options
rdd <- rddensity(X = margin)
rdplotdensity(rdd, margin)

### rdplotdensity: change color and title, restrict to [-50, 50], and add more evaluation points
rdplotdensity(rdd, margin, lcol = c("black", "black"), xlabel = "margin",
              plotRange = c(-50, 50), plotN = 100)

### rdbwdensity: default options
summary(rdbwdensity(X = margin))

### rdbwdensity: compute bandwidth and then use them
tmp <- rdbwdensity(X = margin)
hr <- tmp$h[2,1]
summary(rddensity(X = margin, h=c(10,hr)))

### Other examples
summary(rddensity(X = margin, kernel = "uniform"))
summary(rddensity(X = margin, bwselect = "diff"))
summary(rddensity(X = margin, h = c(10,15)))
summary(rddensity(X = margin, p = 2, q = 4))
summary(rddensity(X = margin, c = 5, all = TRUE))

summary(rdbwdensity(X = margin, p = 3, fitselect = "restricted"))
summary(rdbwdensity(X = margin, kernel = "uniform", vce = "jackknife"))


