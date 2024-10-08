% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rddensity.R
\name{rddensity}
\alias{rddensity}
\title{Manipulation Testing Using Local Polynomial Density Estimation}
\usage{
rddensity(
  X,
  c = 0,
  p = 2,
  q = 0,
  fitselect = "",
  kernel = "",
  vce = "",
  massPoints = TRUE,
  h = c(),
  bwselect = "",
  all = FALSE,
  regularize = TRUE,
  nLocalMin = NULL,
  nUniqueMin = NULL,
  bino = TRUE,
  binoW = NULL,
  binoN = NULL,
  binoWStep = NULL,
  binoNStep = NULL,
  binoNW = 10,
  binoP = 0.5
)
}
\arguments{
\item{X}{Numeric vector or one dimensional matrix/data frame, the running variable.}

\item{c}{Numeric, specifies the threshold or cutoff value in the support of \code{X},
which determines the two samples (e.g., control and treatment units in RD settings).
Default is \code{0}.}

\item{p}{Nonnegative integer, specifies the local polynomial order used to construct
the density estimators. Default is \code{2} (local quadratic approximation).}

\item{q}{Nonnegative integer, specifies the local polynomial order used to construct
the bias-corrected density estimators. Default is \code{p+1} (local cubic approximation
for default \code{p=2}).}

\item{fitselect}{String, specifies the density estimation method.
\code{"unrestricted"} for density estimation without any restrictions (two-sample,
unrestricted inference). This is the default option.
\code{"restricted"} for density estimation assuming equal distribution function and
higher-order derivatives.}

\item{kernel}{String, specifies the kernel function used to construct the local
polynomial estimators.
\code{"triangular"}: \code{K(u)=(1-|u|)*(|u|<=1)}. This is the
default option.
\code{"epanechnikov"}: \code{K(u) = 0.75*(1-u^2)*(|u|<=1)}.
\code{"uniform"}: \code{K(u) = 0.5 * (|u|<=1)}.}

\item{vce}{String, specifies the procedure used to compute the variance-covariance matrix estimator.
\code{"plugin"} for asymptotic plug-in standard errors.
\code{"jackknife"} for jackknife standard errors. This is the default option.}

\item{massPoints}{\code{TRUE} (default) or \code{FALSE}, specifies whether to adjust for
mass points in the data.}

\item{h}{Numeric, specifies the bandwidth used to construct the density estimators on the two
sides of the cutoff.  If not specified, the bandwidth h is computed by the companion command
\code{\link{rdbwdensity}}  If two bandwidths are specified, the first bandwidth is used
for the data below the cutoff and the second bandwidth is used for the data above the cutoff.}

\item{bwselect}{String, specifies the bandwidth selection procedure to be used.
\code{"each"} based on MSE of each density estimator separately (two distinct bandwidths, \code{hl} and \code{hr}).
\code{"diff"} based on MSE of difference of two density estimators (one common bandwidth, \code{hl=hr}).
\code{"sum"} based on MSE of sum of two density estimators (one common bandwidth, \code{hl=hr}).
\code{"comb"} bandwidth is selected as a combination of the alternatives above. This is the default option.
For \code{fitselect="unrestricted"}, it selects \code{median(each,diff,sum)}. For \code{fitselect = "restricted"},
it selects \code{min(diff,sum)}.}

\item{all}{\code{TRUE} or \code{FALSE} (default), if specified, will report two
testing procedures: conventional test statistic (not valid when using MSE-optimal bandwidth choice)
and robust bias-corrected statistic.}

\item{regularize}{\code{TRUE} (default) or \code{FALSE}, specifies whether to conduct local sample size checking.
When set to \code{TRUE}, the bandwidth is chosen such that the local region includes
at least \code{nLocalMin} observations and at least \code{nUniqueMin} unique observations.}

\item{nLocalMin}{Nonnegative integer, specifies the minimum number of observations in each local neighborhood.
This option will be ignored if set to \code{0}, or if \code{regularize=FALSE} is used. Default is \code{20+p+1}.}

\item{nUniqueMin}{Nonnegative integer, specifies the minimum number of unique observations in
each local neighborhood. This option will be ignored if set to \code{0}, or if \code{regularize=FALSE} is used.
Default is \code{20+p+1}.}

\item{bino}{\code{TRUE} (default) or \code{FALSE}, specifies whether to conduct binomial tests. By default,
the initial (smallest) window contains at least 20 observations on each side, and its length is also used as the increment
for subsequent windows. This feature is based on the \code{\link{binom.test}} function.}

\item{binoW}{Numeric, specifies the half length(s) of the initial window. If two values are provided, they will
be used for the data below and above the cutoff separately.}

\item{binoN}{Nonnegative integer, specifies the minimum number of observations on each side of the cutoff used for
the binomial test. This option will be ignored if \code{binoW} is provided.}

\item{binoWStep}{Numeric, specifies the increment in half length(s).}

\item{binoNStep}{Nonnegative integer, specifies the minimum increment in sample size (on each side of the cutoff).
This option will be ignored if \code{binoWStep} is provided.}

\item{binoNW}{Nonnegative integer, specifies the total number of windows. Default is \code{10}.}

\item{binoP}{Numeric, specifies the null hypothesis of the binomial test. Default is \code{0.5}.}
}
\value{
\item{hat}{\code{left}/\code{right}: density estimate to the left/right of cutoff; \code{diff}: difference in
  estimated densities on the two sides of cutoff.}
\item{sd_asy}{\code{left}/\code{right}: standard error for the estimated density to the left/right of the
  cutoff; \code{diff}: standard error for the difference in estimated densities. (Based on
  asymptotic formula.)}
\item{sd_jk}{\code{left}/\code{right}: standard error for the estimated density to the left/right of the
  cutoff; \code{diff}: standard error for the difference in estimated densities. (Based on the
  jackknife method.)}
\item{test}{\code{t_asy}/\code{t_jk}: t-statistic for the density discontinuity test, with standard error
  based on asymptotic formula or the jackknife; \code{p_asy}/\code{p_jk}: p-value for the density
  discontinuity test, with standard error based on asymptotic formula or the jackknife.}
\item{hat_p}{Same as \code{hat}, without bias correction (only available when \code{all=TRUE}).}
\item{sd_asy_p}{Same as \code{sd_asy}, without bias correction (only available when \code{all=TRUE}).}
\item{sd_jk_p}{Same as \code{sd_jk}, without bias correction (only available when \code{all=TRUE}).}
\item{test_p}{Same as \code{test}, without bias correction (only available when \code{all=TRUE}).}
\item{N}{\code{full}: full sample size; \code{left}/\code{right}: sample size to the left/right of the cutoff;
  \code{eff_left}/\code{eff_right}: effective sample size to the left/right of the cutoff (this depends
  on the bandwidth).}
\item{h}{\code{left}/\code{right}: bandwidth used to the left/right of the cutoff.}
\item{opt}{Options passed to the function.}
\item{bino}{Binomial test results. \code{leftWindow}/\code{rightWindow}: window lengths.
   \code{leftN}/\code{rightN}: number of observations. \code{pval}: p-values.}
\item{X_min}{\code{left}/\code{right}: the samllest observation to the left/right of the cutoff.}
\item{X_max}{\code{left}/\code{right}: the largest observation to the left/right of the cutoff.}
}
\description{
\code{rddensity} implements manipulation testing procedures using the
  local polynomial density estimators proposed in Cattaneo, Jansson and Ma (2020),
  and implements graphical procedures with valid confidence bands using the results
  in Cattaneo, Jansson and Ma (2022, 2023).  In addition, the command provides complementary
  manipulation testing based on finite sample exact binomial testing following the
  esults in Cattaneo, Frandsen and Titiunik (2015) and Cattaneo, Frandsen and
  Vazquez-Bare (2017). For an introduction to manipulation testing see McCrary (2008).

 A companion \code{Stata} package is described in Cattaneo, Jansson and Ma (2018).

Companion commands: \code{\link{rdbwdensity}} for data-driven bandwidth selection, and
  \code{\link{rdplotdensity}} for density plots.

Related Stata and R packages useful for inference in regression discontinuity (RD)
  designs are described in the website: \url{https://rdpackages.github.io/}.
}
\examples{
### Continuous Density
set.seed(42)
x <- rnorm(2000, mean = -0.5)
rdd <- rddensity(X = x, vce = "jackknife")
summary(rdd)

### Bandwidth selection using rdbwdensity()
rddbw <- rdbwdensity(X = x, vce = "jackknife")
summary(rddbw)

### Plotting using rdplotdensity()
# 1. From -2 to 2 with 25 evaluation points at each side
plot1 <- rdplotdensity(rdd, x, plotRange = c(-2, 2), plotN = 25)

# 2. Plotting a uniform confidence band
set.seed(42) # fix the seed for simulating critical values
plot2 <- rdplotdensity(rdd, x, plotRange = c(-2, 2), plotN = 25, CIuniform = TRUE)

### Density discontinuity at 0
x[x > 0] <- x[x > 0] * 2
rdd2 <- rddensity(X = x, vce = "jackknife")
summary(rdd2)
plot3 <- rdplotdensity(rdd2, x, plotRange = c(-2, 2), plotN = 25)

}
\references{
Cattaneo, M. D., B. Frandsen, and R. Titiunik. 2015. Randomization Inference in the Regression Discontinuity Design: An Application to the Study of Party Advantages in the U.S. Senate. \emph{Journal of Causal Inference} 3(1): 1-24. \doi{10.1515/jci-2013-0010}

Cattaneo, M. D., M. Jansson, and X. Ma. 2018. Manipulation Testing based on Density Discontinuity. \emph{Stata Journal} 18(1): 234-261. \doi{10.1177/1536867X1801800115}

Cattaneo, M. D., M. Jansson, and X. Ma. 2020. Simple Local Polynomial Density Estimators. \emph{Journal of the American Statistical Association}, 115(531): 1449-1455. \doi{10.1080/01621459.2019.1635480}

Cattaneo, M. D., M. Jansson, and X. Ma. 2022. lpdensity: Local Polynomial Density Estimation and Inference. \emph{Journal of Statistical Software}, 101(2): 1–25. \doi{10.18637/jss.v101.i02}

Cattaneo, M. D., M. Jansson, and X. Ma. 2023. Local Regression Distribution Estimators. \emph{Journal of Econometrics}, 240(2): 105074. \doi{10.1016/j.jeconom.2021.01.006}

Cattaneo, M. D., R. Titiunik and G. Vazquez-Bare. 2017. Comparing Inference Approaches for RD Designs: A Reexamination of the Effect of Head Start on Child Mortality. \emph{Journal of Policy Analysis and Management} 36(3): 643-681. \doi{10.1002/pam.21985}

McCrary, J. 2008. Manipulation of the Running Variable in the Regression Discontinuity Design: A Density Test. \emph{Journal of Econometrics} 142(2): 698-714. \doi{10.1016/j.jeconom.2007.05.005}
}
\seealso{
\code{\link{rdbwdensity}}, \code{\link{rdplotdensity}}
}
\author{
Matias D. Cattaneo, Princeton University  \email{cattaneo@princeton.edu}.

Michael Jansson, University of California Berkeley.  \email{mjansson@econ.berkeley.edu}.

Xinwei Ma (maintainer), University of California San Diego. \email{x1ma@ucsd.edu}.
}
