################################################################################
#' @title Manipulation Testing Using Local Polynomial Density Estimation
#'
#' @description \code{rddensity} implements manipulation testing procedures using the
#'   local polynomial density estimators proposed in Cattaneo, Jansson and Ma (2020),
#'   and implements graphical procedures with valid confidence bands using the results
#'   in Cattaneo, Jansson and Ma (2022, 2023).  In addition, the command provides complementary
#'   manipulation testing based on finite sample exact binomial testing following the
#'   esults in Cattaneo, Frandsen and Titiunik (2015) and Cattaneo, Frandsen and
#'   Vazquez-Bare (2017). For an introduction to manipulation testing see McCrary (2008).
#'
#'  A companion \code{Stata} package is described in Cattaneo, Jansson and Ma (2018).
#'
#' Companion commands: \code{\link{rdbwdensity}} for data-driven bandwidth selection, and
#'   \code{\link{rdplotdensity}} for density plots.
#'
#' Related Stata and R packages useful for inference in regression discontinuity (RD)
#'   designs are described in the website: \url{https://rdpackages.github.io/}.
#'
#' @param X Numeric vector or one dimensional matrix/data frame, the running variable.
#' @param c Numeric, specifies the threshold or cutoff value in the support of \code{X},
#'   which determines the two samples (e.g., control and treatment units in RD settings).
#'   Default is \code{0}.
#' @param p Nonnegative integer, specifies the local polynomial order used to construct
#'   the density estimators. Default is \code{2} (local quadratic approximation).
#' @param q Nonnegative integer, specifies the local polynomial order used to construct
#'   the bias-corrected density estimators. Default is \code{p+1} (local cubic approximation
#'   for default \code{p=2}).
#' @param fitselect String, specifies the density estimation method.
#'   \code{"unrestricted"} for density estimation without any restrictions (two-sample,
#'   unrestricted inference). This is the default option.
#'   \code{"restricted"} for density estimation assuming equal distribution function and
#'   higher-order derivatives.
#' @param kernel String, specifies the kernel function used to construct the local
#'   polynomial estimators.
#'   \code{"triangular"}: \code{K(u)=(1-|u|)*(|u|<=1)}. This is the
#'   default option.
#'   \code{"epanechnikov"}: \code{K(u) = 0.75*(1-u^2)*(|u|<=1)}.
#'   \code{"uniform"}: \code{K(u) = 0.5 * (|u|<=1)}.
#' @param vce String, specifies the procedure used to compute the variance-covariance matrix estimator.
#'   \code{"plugin"} for asymptotic plug-in standard errors.
#'   \code{"jackknife"} for jackknife standard errors. This is the default option.
#' @param massPoints \code{TRUE} (default) or \code{FALSE}, specifies whether to adjust for
#'   mass points in the data.
#' @param all \code{TRUE} or \code{FALSE} (default), if specified, will report two
#'   testing procedures: conventional test statistic (not valid when using MSE-optimal bandwidth choice)
#'   and robust bias-corrected statistic.
#' @param h Numeric, specifies the bandwidth used to construct the density estimators on the two
#'   sides of the cutoff.  If not specified, the bandwidth h is computed by the companion command
#'   \code{\link{rdbwdensity}}  If two bandwidths are specified, the first bandwidth is used
#'   for the data below the cutoff and the second bandwidth is used for the data above the cutoff.
#' @param bwselect String, specifies the bandwidth selection procedure to be used.
#'   \code{"each"} based on MSE of each density estimator separately (two distinct bandwidths, \code{hl} and \code{hr}).
#'   \code{"diff"} based on MSE of difference of two density estimators (one common bandwidth, \code{hl=hr}).
#'   \code{"sum"} based on MSE of sum of two density estimators (one common bandwidth, \code{hl=hr}).
#'   \code{"comb"} bandwidth is selected as a combination of the alternatives above. This is the default option.
#'   For \code{fitselect="unrestricted"}, it selects \code{median(each,diff,sum)}. For \code{fitselect = "restricted"},
#'   it selects \code{min(diff,sum)}.
#' @param regularize \code{TRUE} (default) or \code{FALSE}, specifies whether to conduct local sample size checking.
#'   When set to \code{TRUE}, the bandwidth is chosen such that the local region includes
#'   at least \code{nLocalMin} observations and at least \code{nUniqueMin} unique observations.
#' @param nLocalMin Nonnegative integer, specifies the minimum number of observations in each local neighborhood.
#'   This option will be ignored if set to \code{0}, or if \code{regularize=FALSE} is used. Default is \code{20+p+1}.
#' @param nUniqueMin Nonnegative integer, specifies the minimum number of unique observations in
#'   each local neighborhood. This option will be ignored if set to \code{0}, or if \code{regularize=FALSE} is used.
#'   Default is \code{20+p+1}.
#' @param bino \code{TRUE} (default) or \code{FALSE}, specifies whether to conduct binomial tests. By default,
#'   the initial (smallest) window contains at least 20 observations on each side, and its length is also used as the increment
#'   for subsequent windows. This feature is based on the \code{\link{binom.test}} function.
#' @param binoW Numeric, specifies the half length(s) of the initial window. If two values are provided, they will
#'   be used for the data below and above the cutoff separately.
#' @param binoN Nonnegative integer, specifies the minimum number of observations on each side of the cutoff used for
#'    the binomial test. This option will be ignored if \code{binoW} is provided.
#' @param binoWStep Numeric, specifies the increment in half length(s).
#' @param binoNStep Nonnegative integer, specifies the minimum increment in sample size (on each side of the cutoff).
#'   This option will be ignored if \code{binoWStep} is provided.
#' @param binoNW Nonnegative integer, specifies the total number of windows. Default is \code{10}.
#' @param binoP Numeric, specifies the null hypothesis of the binomial test. Default is \code{0.5}.
#'
#' @return
#' \item{hat}{\code{left}/\code{right}: density estimate to the left/right of cutoff; \code{diff}: difference in
#'   estimated densities on the two sides of cutoff.}
#' \item{sd_asy}{\code{left}/\code{right}: standard error for the estimated density to the left/right of the
#'   cutoff; \code{diff}: standard error for the difference in estimated densities. (Based on
#'   asymptotic formula.)}
#' \item{sd_jk}{\code{left}/\code{right}: standard error for the estimated density to the left/right of the
#'   cutoff; \code{diff}: standard error for the difference in estimated densities. (Based on the
#'   jackknife method.)}
#' \item{test}{\code{t_asy}/\code{t_jk}: t-statistic for the density discontinuity test, with standard error
#'   based on asymptotic formula or the jackknife; \code{p_asy}/\code{p_jk}: p-value for the density
#'   discontinuity test, with standard error based on asymptotic formula or the jackknife.}
#' \item{hat_p}{Same as \code{hat}, without bias correction (only available when \code{all=TRUE}).}
#' \item{sd_asy_p}{Same as \code{sd_asy}, without bias correction (only available when \code{all=TRUE}).}
#' \item{sd_jk_p}{Same as \code{sd_jk}, without bias correction (only available when \code{all=TRUE}).}
#' \item{test_p}{Same as \code{test}, without bias correction (only available when \code{all=TRUE}).}
#' \item{N}{\code{full}: full sample size; \code{left}/\code{right}: sample size to the left/right of the cutoff;
#'   \code{eff_left}/\code{eff_right}: effective sample size to the left/right of the cutoff (this depends
#'   on the bandwidth).}
#' \item{h}{\code{left}/\code{right}: bandwidth used to the left/right of the cutoff.}
#' \item{opt}{Options passed to the function.}
#' \item{bino}{Binomial test results. \code{leftWindow}/\code{rightWindow}: window lengths.
#'    \code{leftN}/\code{rightN}: number of observations. \code{pval}: p-values.}
#' \item{X_min}{\code{left}/\code{right}: the samllest observation to the left/right of the cutoff.}
#' \item{X_max}{\code{left}/\code{right}: the largest observation to the left/right of the cutoff.}
#'
#' @author
#' Matias D. Cattaneo, Princeton University  \email{cattaneo@princeton.edu}.
#'
#' Michael Jansson, University of California Berkeley.  \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @references
#' Cattaneo, M. D., B. Frandsen, and R. Titiunik. 2015. Randomization Inference in the Regression Discontinuity Design: An Application to the Study of Party Advantages in the U.S. Senate. \emph{Journal of Causal Inference} 3(1): 1-24. \doi{10.1515/jci-2013-0010}
#'
#' Cattaneo, M. D., M. Jansson, and X. Ma. 2018. Manipulation Testing based on Density Discontinuity. \emph{Stata Journal} 18(1): 234-261. \doi{10.1177/1536867X1801800115}
#'
#' Cattaneo, M. D., M. Jansson, and X. Ma. 2020. Simple Local Polynomial Density Estimators. \emph{Journal of the American Statistical Association}, 115(531): 1449-1455. \doi{10.1080/01621459.2019.1635480}
#'
#' Cattaneo, M. D., M. Jansson, and X. Ma. 2022. lpdensity: Local Polynomial Density Estimation and Inference. \emph{Journal of Statistical Software}, 101(2): 1–25. \doi{10.18637/jss.v101.i02}
#'
#' Cattaneo, M. D., M. Jansson, and X. Ma. 2023. Local Regression Distribution Estimators. \emph{Journal of Econometrics}, 240(2): 105074. \doi{10.1016/j.jeconom.2021.01.006}
#'
#' Cattaneo, M. D., R. Titiunik and G. Vazquez-Bare. 2017. Comparing Inference Approaches for RD Designs: A Reexamination of the Effect of Head Start on Child Mortality. \emph{Journal of Policy Analysis and Management} 36(3): 643-681. \doi{10.1002/pam.21985}
#'
#' McCrary, J. 2008. Manipulation of the Running Variable in the Regression Discontinuity Design: A Density Test. \emph{Journal of Econometrics} 142(2): 698-714. \doi{10.1016/j.jeconom.2007.05.005}
#'
#' @seealso \code{\link{rdbwdensity}}, \code{\link{rdplotdensity}}
#'
#' @examples
#' ### Continuous Density
#' set.seed(42)
#' x <- rnorm(2000, mean = -0.5)
#' rdd <- rddensity(X = x, vce = "jackknife")
#' summary(rdd)
#'
#' ### Bandwidth selection using rdbwdensity()
#' rddbw <- rdbwdensity(X = x, vce = "jackknife")
#' summary(rddbw)
#'
#' ### Plotting using rdplotdensity()
#' # 1. From -2 to 2 with 25 evaluation points at each side
#' plot1 <- rdplotdensity(rdd, x, plotRange = c(-2, 2), plotN = 25)
#'
#' # 2. Plotting a uniform confidence band
#' set.seed(42) # fix the seed for simulating critical values
#' plot2 <- rdplotdensity(rdd, x, plotRange = c(-2, 2), plotN = 25, CIuniform = TRUE)
#'
#' ### Density discontinuity at 0
#' x[x > 0] <- x[x > 0] * 2
#' rdd2 <- rddensity(X = x, vce = "jackknife")
#' summary(rdd2)
#' plot3 <- rdplotdensity(rdd2, x, plotRange = c(-2, 2), plotN = 25)
#'
#' @export
rddensity <- function(X, c=0, p=2, q=0, fitselect="", kernel="", vce="", massPoints=TRUE,
                      h=c(), bwselect="", all=FALSE,
                      regularize=TRUE, nLocalMin=NULL, nUniqueMin=NULL,
                      bino=TRUE, binoW=NULL, binoN=NULL, binoWStep=NULL, binoNStep=NULL, binoNW=10, binoP=0.5) {

  ################################################################################
  # default values
  ################################################################################

  if (q==0) { q <- p+1 }
  if (kernel == "") { kernel <- "triangular" }
  kernel <- tolower(kernel)
  if (fitselect == "") { fitselect <- "unrestricted" }
  fitselect <- tolower(fitselect)
  if (bwselect == "") { bwselect <- "comb" }
  bwselect <- tolower(bwselect)
  if (vce == "") { vce <- "jackknife" }
  vce <- tolower(vce)
  # end of default values

  ################################################################################
  # bandwidth values
  ################################################################################
  if (length(h) == 0) {
    hl <- hr <- 0
  }
  if (length(h) == 1) {
    hl <- hr <- h
    if (h <= 0) {
      stop("Bandwidth has to be positive.")
    }
  }
  if (length(h) == 2) {
    hl <- h[1]; hr <- h[2]
    if (min(h) <= 0) {
      stop("Bandwidth has to be positive.")
    }
  }
  if (length(h) > 2) {
    stop("No more than two bandwidths are accepted.")
  }

  ################################################################################
  # missing value handling
  ################################################################################
  X <- as.vector(X)
  if (any(is.na(X))) {
    warning(paste(sum(is.na(X)), " missing ", switch((sum(is.na(X))>1)+1, "observation is", "observations are"), " ignored.\n", sep=""))
    X <- X[!is.na(X)]
  }

  ################################################################################
  # sample sizes
  ################################################################################
  X <- sort(X)
  N <- length(X); Nl <- sum(X<c); Nr <- sum(X>=c); Xmin <- min(X); Xmax <- max(X)
  XUnique     <- rddensityUnique(X)
  freqUnique  <- XUnique$freq
  indexUnique <- XUnique$indexLast
  XUnique     <- XUnique$unique
  NUnique     <- length(XUnique)
  NlUnique    <- sum(XUnique <  c)
  NrUnique    <- sum(XUnique >= c)

  if (sum(freqUnique != 1) > 0 & massPoints) {
    masspoints_flag <- 1
  } else {
    masspoints_flag <- 0
  }
  # end of sample sizes

  ################################################################################
  # error handling
  ################################################################################
  if (c <= Xmin | c >= Xmax) { stop("The cutoff should be set within the range of the data.") }
#  if (Nl <= 10 | Nr <= 10) { stop("Not enough observations to perform calculations.") }
  if (p!=1 & p!=2 & p!=3 & p!=4 & p!=5 & p!=6 & p!= 7) { stop("p must be an integer between 1 and 7.") }
  if (p > q) { stop("q cannot be smaller than p.") }
  if (kernel!="uniform" & kernel!= "triangular" & kernel!="epanechnikov") { stop("kernel incorrectly specified.") }
  if (fitselect!="unrestricted" & fitselect!="restricted") { stop("fitselect incorrectly specified.") }
  if (fitselect=="restricted" & hl!=hr) { stop("Bandwidths must be equal in restricted model.") }
  if (bwselect!="each" & bwselect!="diff" & bwselect!="sum" & bwselect!="comb") { stop("bwselect incorrectly specified.") }
  if (fitselect=="restricted" & bwselect=="each") { stop("bwselect=each is not available in the restricted model.") }
  if (vce!="plugin" & vce!="jackknife") { stop("vce incorrectly specified.") }

  # regularize
  if (length(regularize) == 0) {
    regularize <- TRUE
  } else if (length(regularize) > 1 | !regularize[1]%in%c(TRUE, FALSE)) {
    stop("Regularization parameter incorrectly specified.\n")
  }

  # nLocalMin
  if (length(nLocalMin) == 0) { nLocalMin <- 20 + p + 1 }
  if (!is.numeric(nLocalMin) | is.na(nLocalMin)) {
    stop("Option nLocalMin incorrectly specified.\n")
  } else if (ceiling(nLocalMin) < 0) {
    stop("Option nLocalMin incorrectly specified.\n")
  } else {
    nLocalMin <- ceiling(nLocalMin)
  }

  # nUniqueMin
  if (length(nUniqueMin) == 0) { nUniqueMin <- 20 + p + 1 }
  if (!is.numeric(nUniqueMin) | is.na(nUniqueMin)) {
    stop("Option nUniqueMin incorrectly specified.\n")
  } else if (ceiling(nUniqueMin) < 0) {
    stop("Option nUniqueMin incorrectly specified.\n")
  } else {
    nUniqueMin <- ceiling(nUniqueMin)
  }

  # massPoints
  if (length(massPoints) == 0) {
    massPoints <- TRUE
  } else if (length(massPoints) > 1 | !massPoints[1]%in%c(TRUE, FALSE)) {
    stop("Option massPoints incorrectly specified.\n")
  }
  # end of error handling

  ################################################################################
  # bandwidth selection
  ################################################################################
  if (hl > 0 & hr > 0) { bwselectl <- "mannual"} else { bwselectl <- "estimated"
    out <- rdbwdensity(X=X, c=c, p=p, kernel=kernel, fitselect=fitselect, vce=vce, regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin, massPoints=massPoints)$h
    if (fitselect=="unrestricted" & bwselect=="each" & hl==0)  hl = out[1,1]
  	if (fitselect=="unrestricted" & bwselect=="each" & hr==0)  hr = out[2,1]
		if (fitselect=="unrestricted" & bwselect=="diff" & hl==0)  hl = out[3,1]
		if (fitselect=="unrestricted" & bwselect=="diff" & hr==0)  hr = out[3,1]
		if (fitselect=="unrestricted" & bwselect=="sum"  & hl==0)  hl = out[4,1]
		if (fitselect=="unrestricted" & bwselect=="sum"  & hr==0)  hr = out[4,1]
		if (fitselect=="unrestricted" & bwselect=="comb" & hl==0)  hl = median(c(out[1,1], out[3,1], out[4,1]))
		if (fitselect=="unrestricted" & bwselect=="comb" & hr==0)  hr = median(c(out[2,1], out[3,1], out[4,1]))

		if (fitselect=="restricted" & bwselect=="diff" & hl==0)  hl = out[3,1]
		if (fitselect=="restricted" & bwselect=="diff" & hr==0)  hr = out[3,1]
		if (fitselect=="restricted" & bwselect=="sum"  & hl==0)  hl = out[4,1]
		if (fitselect=="restricted" & bwselect=="sum"  & hr==0)  hr = out[4,1]
		if (fitselect=="restricted" & bwselect=="comb" & hl==0)  hl = min(c(out[3,1], out[4,1]))
		if (fitselect=="restricted" & bwselect=="comb" & hr==0)  hr = min(c(out[3,1], out[4,1]))
  }
  # end of bandwidth selection

  ################################################################################
  # data trimming
  ################################################################################
  X <- X - c
  Y <- (0:(N-1)) / (N-1)
  if (massPoints) {
    Y <- rep(Y[indexUnique], times=freqUnique)
  }
  Xh <- X[(X>=-1*hl) & (X<=hr)]; Yh <- Y[(X>=-1*hl) & (X<=hr)]
  Nlh <- sum(Xh < 0); Nrh <- sum(Xh >= 0)
  #if (Nlh < 5 | Nrh < 5) { stop("Not enough observations to perform calculation.") }
  Nh <- Nlh + Nrh
  # end of data trimming

  ################################################################################
  # estimation
  ################################################################################
  fV_q <- rddensity_fV(Y=Yh, X=Xh, Nl=Nl, Nr=Nr, Nlh=Nlh, Nrh=Nrh, hl=hl, hr=hr, p=q, s=1, kernel=kernel, fitselect=fitselect, vce=vce, massPoints)
  T_asy <- fV_q[3, 1] / sqrt(fV_q[3, 3]); T_jk <- fV_q[3, 1] / sqrt(fV_q[3, 2])
  p_asy <- pnorm(abs(T_asy), lower.tail=FALSE) * 2; p_jk <- pnorm(abs(T_jk), lower.tail=FALSE) * 2

  ################################################################################
  # binomial test
  ################################################################################
  if (bino) {
    XSort <- sort(abs(X))
    XL <- sort(abs(X[X < 0]))
    XR <- sort(X[X >= 0])

    # binoNW
    if (length(binoNW) == 0) {
      binoNW <- 10
    } else if (length(binoNW) > 1) {
      stop("Option binoNW incorrectly specified.\n")
    } else if (binoNW<=0) {
      stop("Option binoNW incorrectly specified.\n")
    } else {
      binoNW <- ceiling(binoNW)
    }
    binomTempLW <- rep(NA, binoNW)
    binomTempRW <- rep(NA, binoNW)

    # binoP check
    if (length(binoP) > 1 | !all(is.numeric(binoP))){
      stop("Option binoP incorrectly specified.\n")
    }
    if (binoP < 0 | binoP > 1) {
      stop("Option binoP incorrectly specified.\n")
    }

    # binoW check
    if (length(binoW) == 0) {
      # binoN check
      if (length(binoN) == 0) {
        binoN <- 20
        # default, minimum 20 obs on each side
        binomTempLW[1] <- binomTempRW[1] <- max(XL[min(binoN, Nl)], XR[min(binoN, Nr)])
      } else if (length(binoN) == 1) {
        if (binoN <= 0) {
          stop("Option binoN incorrectly specified.\n")
        }
        binoN <- ceiling(binoN)
        binomTempLW[1] <- binomTempRW[1] <- max(XL[min(binoN, Nl)], XR[min(binoN, Nr)])
      } else {
        stop("Option binoN incorrectly specified.\n")
      }
    } else if (length(binoW) == 1) {
      if (binoW <= 0) {
        stop("Option binoW incorrectly specified.\n")
      }
      binomTempLW[1] <- binomTempRW[1] <- binoW
      binoN <- min(sum(XL <= binomTempLW[1]), sum(XR <= binomTempRW[1]))
    } else if (length(binoW) == 2) {
      if (binoW[1] <= 0 | binoW[2] <= 0) {
        stop("Option binoW incorrectly specified.\n")
      }
      binomTempLW[1] <- binoW[1]
      binomTempRW[1] <- binoW[2]
      binoN <- min(sum(XL <= binomTempLW[1]), sum(XR <= binomTempRW[1]))
    } else {
      stop("Option binoW incorrectly specified.\n")
    }

    # binoWStep check
    if (binoNW > 1) {
      if (length(binoWStep) == 0) {
        # binoNStep check
        if (length(binoNStep) == 0) {
          # default
          if (binomTempLW[1] >= hl | binomTempRW[1] >= hr) {
            binomTempLW <- binomTempLW[1]
            binomTempRW <- binomTempRW[1]
            binoNW <- 1
          } else {
            if (binomTempLW[1] * binoNW > hl) {
              binomTempLW[2:binoNW] <- binomTempLW[1] + (1:(binoNW-1)) * (hl - binomTempLW[1]) / (binoNW-1)
            } else {
              binomTempLW[2:binoNW] <- binomTempLW[1] + (1:(binoNW-1)) * binomTempLW[1]
            }

            if (binomTempRW[1] * binoNW > hr) {
              binomTempRW[2:binoNW] <- binomTempRW[1] + (1:(binoNW-1)) * (hr - binomTempRW[1]) / (binoNW-1)
            } else {
              binomTempRW[2:binoNW] <- binomTempRW[1] + (1:(binoNW-1)) * binomTempRW[1]
            }
          }
        } else if (length(binoNStep) == 1) {
          if (binoNStep <= 0) {
            stop("Option binoNStep incorrectly specified.\n")
          }
          binoNStep <- ceiling(binoNStep)
          for (jj in 2:binoNW) {
            binomTempLW[jj] <- binomTempLW[jj-1] + max(XL[min(sum(XL<=binomTempLW[jj-1]) + binoNStep, Nl)] - binomTempLW[jj-1],
                                                       XR[min(sum(XR<=binomTempRW[jj-1]) + binoNStep, Nr)] - binomTempRW[jj-1])

            binomTempRW[jj] <- binomTempRW[jj-1] + max(XL[min(sum(XL<=binomTempLW[jj-1]) + binoNStep, Nl)] - binomTempLW[jj-1],
                                                       XR[min(sum(XR<=binomTempRW[jj-1]) + binoNStep, Nr)] - binomTempRW[jj-1])
          }
        } else {
          stop("Option binoNStep incorrectly specified.\n")
        }
      } else if (length(binoWStep) == 1) {
        if (binoWStep <= 0) {
          stop("Option binoWStep incorrectly specified.\n")
        }
        binomTempLW[2:binoNW] <- binomTempLW[1] + (1:(binoNW-1)) * binoWStep
        binomTempRW[2:binoNW] <- binomTempRW[1] + (1:(binoNW-1)) * binoWStep
      } else if (length(binoWStep) == 2) {
        if (binoWStep[1] <= 0 | binoWStep[2] <= 0) {
          stop("Option binoWStep incorrectly specified.\n")
        }
        binomTempLW[2:binoNW] <- binomTempLW[1] + (1:(binoNW-1)) * binoWStep[1]
        binomTempRW[2:binoNW] <- binomTempRW[1] + (1:(binoNW-1)) * binoWStep[2]
      } else {
        stop("Option binoWStep incorrectly specified.\n")
      }
    }

    binomTempLN <- binomTempRN <- binomTempPVal <- rep(NA, binoNW)
    for (jj in 1:binoNW) {
      binomTempLN[jj] <- sum(XL <= binomTempLW[jj])
      binomTempRN[jj] <- sum(XR <= binomTempRW[jj])
      binomTempPVal[jj] <- binom.test(x = c(binomTempLN[jj], binomTempRN[jj]), p = binoP)$p.value
    }
  } else {
    binomTempLN <- binomTempRN <- binomTempLW <- binomTempRW <- binomTempPVal <- NA
  }

  if (all) {
    fV_p <- rddensity_fV(Y=Yh, X=Xh, Nl=Nl, Nr=Nr, Nlh=Nlh, Nrh=Nrh, hl=hl, hr=hr, p=p, s=1, kernel=kernel, fitselect=fitselect, vce=vce, massPoints)
    T_asy_p <- fV_p[3, 1] / sqrt(fV_p[3, 3]); T_jk_p <- fV_p[3, 1] / sqrt(fV_p[3, 2])
    p_asy_p <- pnorm(abs(T_asy_p), lower.tail=FALSE) * 2; p_jk_p <- pnorm(abs(T_jk_p), lower.tail=FALSE) * 2

    result <- list( hat=   list(left=fV_q[1,1], right=fV_q[2,1], diff=fV_q[3,1]),
                    sd_asy=list(left=sqrt(fV_q[1,3]), right=sqrt(fV_q[2,3]), diff=sqrt(fV_q[3,3])),
                    sd_jk= list(left=sqrt(fV_q[1,2]), right=sqrt(fV_q[2,2]), diff=sqrt(fV_q[3,2])),
                    test=  list(t_asy=T_asy, t_jk=T_jk, p_asy=p_asy, p_jk=p_jk),

                    hat_p=   list(left=fV_p[1,1], right=fV_p[2,1], diff=fV_p[3,1]),
                    sd_asy_p=list(left=sqrt(fV_p[1,3]), right=sqrt(fV_p[2,3]), diff=sqrt(fV_p[3,3])),
                    sd_jk_p= list(left=sqrt(fV_p[1,2]), right=sqrt(fV_p[2,2]), diff=sqrt(fV_p[3,2])),
                    test_p=  list(t_asy=T_asy_p, t_jk=T_jk_p, p_asy=p_asy_p, p_jk=p_jk_p),

                    N=     list(full=N, left=Nl, right=Nr, eff_left=Nlh, eff_right=Nrh),
                    h=     list(left=hl, right=hr),
                    opt=   list(fitselect=fitselect, kernel=kernel, bwselectl=bwselectl,
                                vce=vce, c=c, p=p, q=q, all=all,
                                regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin,
                                massPoints=massPoints, masspoints_flag=masspoints_flag,
                                bino=bino, binoN=binoN, binoW=binoW, binoNStep=binoNStep, binoWStep=binoWStep, binoNW=binoNW, binoP=binoP),
                    X_min      =list(left=min(X[X<0])+c, right=min(X[X>=0])+c),
                    X_max      =list(left=max(X[X<0])+c, right=max(X[X>=0])+c),
                    bino = list(LeftN=binomTempLN, RightN=binomTempRN, LeftWindow=binomTempLW, RightWindow=binomTempRW, pval=binomTempPVal))
  } else {
    result <- list( hat=   list(left=fV_q[1,1], right=fV_q[2,1], diff=fV_q[3,1]),
                    sd_asy=list(left=sqrt(fV_q[1,3]), right=sqrt(fV_q[2,3]), diff=sqrt(fV_q[3,3])),
                    sd_jk= list(left=sqrt(fV_q[1,2]), right=sqrt(fV_q[2,2]), diff=sqrt(fV_q[3,2])),
                    test=  list(t_asy=T_asy, t_jk=T_jk, p_asy=p_asy, p_jk=p_jk),

                    hat_p=   list(left=NA, right=NA, diff=NA),
                    sd_asy_p=list(left=NA, right=NA, diff=NA),
                    sd_jk_p= list(left=NA, right=NA, diff=NA),
                    test_p=  list(t_asy=NA, t_jk=NA, p_asy=NA, p_jk=NA),

                    N=     list(full=N, left=Nl, right=Nr, eff_left=Nlh, eff_right=Nrh),
                    h=     list(left=hl, right=hr),
                    opt=   list(fitselect=fitselect, kernel=kernel, bwselectl=bwselectl,
                                vce=vce, c=c, p=p, q=q, all=all,
                                regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin,
                                massPoints=massPoints, masspoints_flag=masspoints_flag,
                                bino=bino, binoN=binoN, binoW=binoW, binoNStep=binoNStep, binoWStep=binoWStep, binoNW=binoNW, binoP=binoP),
                    X_min      =list(left=min(X[X<0])+c, right=min(X[X>=0])+c),
                    X_max      =list(left=max(X[X<0])+c, right=max(X[X>=0])+c),
                    bino = list(LeftN=binomTempLN, RightN=binomTempRN, LeftWindow=binomTempLW, RightWindow=binomTempRW, pval=binomTempPVal))
  }


  class(result) <- "CJMrddensity"
  return(result)
}

################################################################################
#' Internal function.
#'
#' @param object Class \code{CJMrddensity} objects.
#'
#' @keywords internal
#' @export
summary.CJMrddensity <- function(object, ...) {
  x <- object
  cat("\nManipulation testing using local polynomial density estimation.\n")
  cat("\n")

  cat(paste(format("Number of obs =", width=22), toString(x$N$full), sep="")); cat("\n")
  cat(paste(format("Model =", width=22), x$opt$fitselect, sep="")); cat("\n")
  cat(paste(format("Kernel =", width=22), x$opt$kernel, sep="")); cat("\n")
  if (x$opt$bwselectl!="mannual") {
    cat(paste(format("BW method =", width=22), x$opt$bwselect, sep="")); cat("\n")
  } else {
    cat(paste(format("BW method =", width=22), x$opt$bwselectl, sep="")); cat("\n")
  }
  cat(paste(format("VCE method =", width=22), x$opt$vce, sep="")); cat("\n")
  cat("\n")

  cat(paste(format(paste("c = ", toString(round(x$opt$c, 3)), sep=""), width=22), format("Left of c", width=20), format("Right of c", width=20), sep="")); cat("\n")
  cat(paste(format("Number of obs", width=22), format(toString(x$N$left), width=20), format(toString(x$N$right), width=20), sep="")); cat("\n")
  cat(paste(format("Eff. Number of obs", width=22), format(toString(x$N$eff_left), width=20), format(toString(x$N$eff_right), width=20), sep="")); cat("\n")
  #cat(paste(format("Min Running var.", width=22), format(toString(round(x$X_min$left, 3)), width=20), format(toString(round(x$X_min$right, 3)), width=20), sep="")); cat("\n")
  #cat(paste(format("Max Running var.", width=22), format(toString(round(x$X_max$left, 3)), width=20), format(toString(round(x$X_max$right, 3)), width=20), sep="")); cat("\n")
  cat(paste(format("Order est. (p)", width=22), format(toString(x$opt$p), width=20), format(toString(x$opt$p), width=20), sep="")); cat("\n")
  cat(paste(format("Order bias (q)", width=22), format(toString(x$opt$q), width=20), format(toString(x$opt$q), width=20), sep="")); cat("\n")
  cat(paste(format("BW est. (h)", width=22), format(toString(round(x$h$left, 3)), width=20), format(toString(round(x$h$right, 3)), width=20), sep="")); cat("\n")
  cat("\n")

  cat(paste(format("Method", width=22), format("T", width=20), format("P > |T|", width=20), sep="")); cat("\n")
  if (x$opt$all) {
    if (x$opt$vce == "plugin") {
      cat(paste(format("Conventional", width=22), format(toString(round(x$test_p$t_asy, 4)), width=20), format(toString(round(x$test_p$p_asy, 4)), width=20), sep="")); cat("\n")
    } else {
      cat(paste(format("Conventional", width=22), format(toString(round(x$test_p$t_jk, 4)), width=20), format(toString(round(x$test_p$p_jk, 4)), width=20), sep="")); cat("\n")
    }
  }
  if (x$opt$vce == "plugin") {
    cat(paste(format("Robust", width=22), format(toString(round(x$test$t_asy, 4)), width=20), format(toString(round(x$test$p_asy, 4)), width=20), sep="")); cat("\n")
  } else {
    cat(paste(format("Robust", width=22), format(toString(round(x$test$t_jk, 4)), width=20), format(toString(round(x$test$p_jk, 4)), width=20), sep="")); cat("\n")
  }
  cat("\n")

  if (x$h$left > x$opt$c - x$X_min$left) {
    warning("Bandwidth hl greater than the range of the data.")
  }
  if (x$h$right>x$X_max$right - x$opt$c) {
    warning("Bandwidth hr greater than the range of the data.")
  }
  if (x$N$eff_left < 20 | x$N$eff_right < 20) {
    warning("Bandwidth h may be too small.")
  }
	if (x$opt$masspoints_flag) {
	  warning("There are repeated observations. Point estimates and standard errors have been adjusted. Use option massPoints=FALSE to suppress this feature.")
	}

  if (x$opt$bino) {
    cat(paste("\nP-values of binomial tests (H0: p=", x$opt$binoP, ").\n", sep=""))
    cat("\n")

    if (all(x$bino$LeftWindow==x$bino$RightWindow)) {
      cat(paste(format("Window Length / 2", width=21, justify="left"), " ",  "     <c", " ", "    >=c", " ", "   P>|T|", sep="")); cat("\n")
      for (jj in 1:length(x$bino$LeftWindow)) {
        print_flag <- FALSE
        if (jj == 1) {
          print_flag <- TRUE
        } else {
          if (x$bino$LeftWindow[jj] != x$bino$LeftWindow[jj-1] | x$bino$RightWindow[jj] != x$bino$RightWindow[jj-1]){
            print_flag <- TRUE
          }
        }
        if (print_flag) {
          cat(format(sprintf("%5.3f", x$bino$LeftWindow[jj]), width=21, justify="left"), " ",
              format(x$bino$LeftN[jj], width=7), " ", format(x$bino$RightN[jj], width=7), "    ",
              format(sprintf("%1.4f", x$bino$pval[jj]), width=6, justify="left"),
              sep=""); cat("\n")
        }
      }
    } else {
      cat(paste(format("Window Length", width=21, justify="left"), " ",  "     <c", " ", "    >=c", " ", "   P>|T|", sep="")); cat("\n")
      for (jj in 1:length(x$bino$LeftWindow)) {
        print_flag <- FALSE
        if (jj == 1) {
          print_flag <- TRUE
        } else {
          if (x$bino$LeftWindow[jj] != x$bino$LeftWindow[jj-1] | x$bino$RightWindow[jj] != x$bino$RightWindow[jj-1]){
            print_flag <- TRUE
          }
        }
        if (print_flag) {
          cat(format(sprintf("%5.3f", x$bino$LeftWindow[jj]), width=9, justify="left"), " + ",
              format(sprintf("%5.3f", x$bino$RightWindow[jj]), width=9, justify="left"), " ",
              format(x$bino$LeftN[jj], width=7), " ", format(x$bino$RightN[jj], width=7), "    ",
              format(sprintf("%1.4f", x$bino$pval[jj]), width=6, justify="left"),
              sep=""); cat("\n")
        }
      }
    }
  }
}

################################################################################
#' Internal function.
#'
#' @param x Class \code{CJMrddensity} objects.
#'
#' @keywords internal
#' @export
print.CJMrddensity <- function(x, ...) {
  cat("Call:\n")
  cat("rddensity.\n")
  cat("Sample size:\ ", x$N$full, ". ", "Cutoff: ", x$opt$c, ".\n", sep="")
  cat("Model:\ ", x$opt$fitselect, ". ", "Kernel: ", x$opt$kernel, ". ", "VCE: ", x$opt$vce, sep="")
}
