###########################################################################
## RDDENSITY R Package
## Do-file for RDDENSITY Plot Illustration
## Authors: Matias D. Cattaneo, Michael Jansson and Xinwei Ma  
###########################################################################
### Clear R environment
rm(list=ls(all=TRUE))

### Install R library
### NOTE: depending on your system, you may need to do it as root
#install.packages('rddensity') 

### Load rddensity package
library(rddensity)
library(ggplot2)

### Load data base
data("rddensity_senate")

### Plot and save results
rdplot <- rdplotdensity(rddensity(X = margin), X = margin)

### Construct the same plot
# CI: below the cutoff
rdplotLeft <- as.data.frame(rdplot$Estl$Estimate)
rdplotLeft$cil <- rdplotLeft$f_q - qnorm(0.975) * rdplotLeft$se_q
rdplotLeft$ciu <- rdplotLeft$f_q + qnorm(0.975) * rdplotLeft$se_q
# CI: above the cutoff
rdplotRight <- as.data.frame(rdplot$Estr$Estimate)
rdplotRight$cil <- rdplotRight$f_q - qnorm(0.975) * rdplotRight$se_q
rdplotRight$ciu <- rdplotRight$f_q + qnorm(0.975) * rdplotRight$se_q
# histogram
NLeft <- sum(margin>=min(rdplotLeft$grid) & margin<0)
histNLeft <- ceiling(min(sqrt(NLeft), 10 * log(NLeft)/log(10)))

NRight <- sum(margin<=max(rdplotRight$grid) & margin>=0)
histNRight <- ceiling(min(sqrt(NRight), 10 * log(NRight)/log(10)))

histBreaks <- c(seq(min(rdplotLeft$grid), 0, length.out = histNRight+1), seq(0, max(rdplotRight$grid), length.out = histNRight+1)[2:(histNRight+1)])
histScale <- mean(margin>=min(rdplotLeft$grid) & margin<=max(rdplotRight$grid))

### Call ggplot()
ggplot() + theme_bw() + 
  geom_histogram(data=as.data.frame(margin), aes(x=margin, y=..density..*histScale), breaks=histBreaks, fill=3, col="white", alpha=0.2) +
  geom_ribbon(data=rdplotLeft , aes(x=grid, ymin=cil, ymax=ciu), alpha=0.2, fill="black") + 
  geom_ribbon(data=rdplotRight, aes(x=grid, ymin=cil, ymax=ciu), alpha=0.2, fill="red") + 
  geom_line(data=rdplotLeft , aes(x=grid, y=f_p), col="black") +
  geom_line(data=rdplotRight, aes(x=grid, y=f_p), col="red") + 
  labs(x="margin", y="")
  