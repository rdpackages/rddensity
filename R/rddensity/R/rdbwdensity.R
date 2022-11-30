################################################################################
#' @title Bandwidth Selection for Manipulation Testing
#'
#' @description \code{rdbwdensity} implements several data-driven bandwidth selection
#'   methods useful to construct manipulation testing procedures using the local
#'   polynomial density estimators proposed in Cattaneo, Jansson and Ma (2020).
#'
#'  A companion \code{Stata} package is described in Cattaneo, Jansson and Ma (2018).
#'
#' Companion command: \code{\link{rddensity}} for manipulation (density discontinuity)
#'   testing.
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
#' @param regularize \code{TRUE} (default) or \code{FALSE}, specifies whether to conduct local sample size checking.
#'   When set to \code{TRUE}, the bandwidth is chosen such that the local region includes
#'   at least \code{nLocalMin} observations and at least \code{nUniqueMin} unique observations.
#' @param nLocalMin Nonnegative integer, specifies the minimum number of observations in each local neighborhood.
#'   This option will be ignored if set to \code{0}, or if \code{regularize=FALSE} is used. Default is \code{20+p+1}.
#' @param nUniqueMin Nonnegative integer, specifies the minimum number of unique observations in
#'   each local neighborhood. This option will be ignored if set to \code{0}, or if \code{regularize=FALSE} is used.
#'   Default is \code{20+p+1}.
#'
#' @return
#' \item{h}{Bandwidths for density discontinuity test, left and right to the cutoff, and asymptotic variance and bias.}
#' \item{N}{\code{full}: full sample size; \code{left}/\code{right}: sample size to the left/right of the cutoff.}
#' \item{opt}{Options passed to the function.}
#' \item{X_min}{Smallest observations to the left and right of the cutoff.}
#' \item{X_max}{Largest observations to the left and right of the cutoff.}
#'
#' @author
#' Matias D. Cattaneo, Princeton University  \email{cattaneo@princeton.edu}.
#'
#' Michael Jansson, University of California Berkeley.  \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @references
#' Cattaneo, M. D., M. Jansson, and X. Ma. 2018. \href{https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2018_Stata.pdf}{Manipulation Testing based on Density Discontinuity}. \emph{Stata Journal} 18(1): 234-261. \doi{10.1177/1536867X1801800115}
#'
#' Cattaneo, M. D., M. Jansson, and X. Ma. 2020. \href{https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA.pdf}{Simple Local Polynomial Density Estimators}. \emph{Journal of the American Statistical Association}, 115(531): 1449-1455. \doi{10.1080/01621459.2019.1635480}
#'
#' @seealso \code{\link{rddensity}}
#'
#' @examples
#' # Generate a random sample
#' set.seed(42)
#' x <- rnorm(2000, mean = -0.5)
#'
#' # Bandwidth selection
#' summary(rdbwdensity(X = x, vce="jackknife"))
#'
#' @export
rdbwdensity <- function(X, c=0, p=2, fitselect="", kernel="", vce="", massPoints=TRUE, regularize=TRUE, nLocalMin=NULL, nUniqueMin=NULL) {

  ################################################################################
  # default values
  ################################################################################
  if (kernel == "") { kernel <- "triangular" }
  kernel <- tolower(kernel)
  if (fitselect == "") { fitselect <- "unrestricted" }
  fitselect <- tolower(fitselect)
  if (vce == "") { vce <- "jackknife" }
  vce <- tolower(vce)
  # end of default values

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
  X <- sort(X, decreasing=FALSE)
  N <- length(X); Nl <- sum(X<c); Nr <- sum(X>=c); Xmin <- min(X); Xmax <- max(X)
  XUnique     <- rddensityUnique(X)
  freqUnique  <- XUnique$freq
  indexUnique <- XUnique$indexLast
  XUnique     <- XUnique$unique
  NUnique     <- length(XUnique)
  NlUnique    <- sum(XUnique <  c)
  NrUnique    <- sum(XUnique >= c)

  X <- X - c; Xmu <- mean(X); Xsd <- sd(X)
  XUnique <- XUnique - c

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
  if (kernel!="uniform" & kernel!= "triangular" & kernel!="epanechnikov") { stop("kernel incorrectly specified.") }
  if (fitselect!="unrestricted" & fitselect!="restricted") { stop("fitselect incorrectly specified.") }
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
  # select preliminary bandwidth
  ################################################################################
  fhatb <- 1 / (rddensity_H(Xmu / Xsd, p+2)^2 * dnorm(Xmu / Xsd))
  fhatc <- 1 / (rddensity_H(Xmu / Xsd, p)^2 * dnorm(Xmu / Xsd))
  # these constants are for uniform kernel
  Cb <- c(25884.444444494150957,3430865.4551236177795,845007948.04262602329,330631733667.03808594,187774809656037.3125,145729502641999264,146013502974449876992)
  Cc <- c(4.8000000000000246914,548.57142857155463389,100800.00000020420703,29558225.458100609481,12896196859.612621307,7890871468221.609375,6467911284037581)
  bn <- ((2*p+1)/4 * fhatb * Cb[p] / N)^(1/(2*p+5))
  cn <- (1/(2*p) * fhatc * Cc[p] / N)^(1/(2*p+1))
  bn <- bn * Xsd; cn <- cn * Xsd

  # bn is for higher-order derivative estimation
  # cn is for density estimation
  if (regularize) {
    # bandwidth should not exceed the range of data
    bn <- min(bn, max(abs(XUnique)))
    cn <- min(cn, max(abs(XUnique)))

    # nLocalMin check
    if (nLocalMin > 0) {
      bn <- max(bn, sort(abs(X[X < 0]))[min(20 + p+2 + 1, Nl)], (X[X >= 0])[min(20 + p+2 + 1, Nr)])
      cn <- max(cn, sort(abs(X[X < 0]))[min(20 + p   + 1, Nl)], (X[X >= 0])[min(20 + p   + 1, Nr)])
    }

    # nUniqueMin check
    if (nUniqueMin > 0) {
      bn <- max(bn, sort(abs(XUnique[XUnique < 0]))[min(20 + p+2 + 1, NlUnique)], (XUnique[XUnique >= 0])[min(20 + p+2 + 1, NrUnique)])
      cn <- max(cn, sort(abs(XUnique[XUnique < 0]))[min(20 + p   + 1, NlUnique)], (XUnique[XUnique >= 0])[min(20 + p   + 1, NrUnique)])
    }
  }
  # end of selection of preliminary bandwidth
#cat(bn)
#cat("\n")
#cat(cn)
  ################################################################################
  # estimate main bandwidth
  ################################################################################
  # mass points correction for the empirical distribution function
  Y <- (0:(N-1)) / (N-1)
Y0 <- Y
  if (massPoints) {
    Y <- rep(Y[indexUnique], times=freqUnique)
  }

  Yb <- Y[abs(X) <= bn]; Xb <- X[abs(X) <= bn]; Yc <- Y[abs(X) <= cn]; Xc <- X[abs(X) <= cn]
  Nlb <- sum(Xb < 0); Nrb <- sum(Xb >= 0); Nlc <- sum(Xc < 0); Nrc <- sum(Xc >= 0)

  hn <- matrix(NA, ncol=3, nrow=4)
  colnames(hn) <- c("bw", "variance", "biassq"); rownames(hn) <- c("l", "r", "diff", "sum")
  fV_b <- rddensity_fV(Y=Yb, X=Xb, Nl=Nl, Nr=Nr, Nlh=Nlb, Nrh=Nrb, hl=bn, hr=bn, p=p+2, s=p+1, kernel=kernel, fitselect=fitselect, vce=vce, massPoints)
  fV_c <- rddensity_fV(Y=Yc, X=Xc, Nl=Nl, Nr=Nr, Nlh=Nlc, Nrh=Nrc, hl=cn, hr=cn, p=p,   s=1,   kernel=kernel, fitselect=fitselect, vce=vce, massPoints)

  if (vce == "plugin") { hn[, 2] <- N * cn * fV_c[, 3]  } else { hn[, 2] <- N * cn * fV_c[, 2] }
  if (fitselect == "unrestricted") {
    S <- Sgenerate(p=p, low=0, up=1, kernel=kernel); C <- Cgenerate(k=p+1, p=p, low=0, up=1, kernel=kernel)
    hn[1, 3] <- fV_b[1, 4] * (solve(S) %*% C)[2] * (-1)^p
    hn[2, 3] <- fV_b[2, 4] * (solve(S) %*% C)[2]
    hn[3, 3] <- hn[2, 3] - hn[1, 3]; hn[4, 3] <- hn[2, 3] + hn[1, 3]
  } else {
    Splus <- Splusgenerate(p=p, kernel=kernel); Cplus <- Cplusgenerate(k=p+1, p=p, kernel=kernel); Psi <- Psigenerate(p=p)
    Sinv <- solve(fV_c[2, 1] * Splus + fV_c[1, 1] * Psi%*%Splus%*%Psi)
    C <- fV_b[1, 4] * (fV_c[2, 1] * Cplus + (-1)^(p+1) * fV_c[1, 1] * Psi%*%Cplus)
    temp <- Sinv%*%C
    hn[1, 3] <- temp[2]; hn[2, 3] <- temp[3]; hn[3, 3] <- hn[2, 3] - hn[1, 3]; hn[4, 3] <- hn[2, 3] + hn[1, 3]
  }
#print(hn)
  hn[, 3] <- hn[, 3]^2
  hn[, 1] <- (1/(2*p) * hn[, 2] / hn[, 3] / N)^(1/(2*p+1))
  # end of estimating main bandwidth

  for (i in 1:4) {
    if (hn[i, 2] < 0) { hn[i, 1] <- 0; hn[i, 2] <- NA }
    if (is.na(hn[i, 1])) { hn[i, 1] <- 0 }
  }
#print(hn)
  # bandwidth regularization
  if (regularize) {
    # bandwidth should not exceed the range of data
    hn[1,1] <- min(hn[1,1], abs(XUnique[1]))
    hn[2,1] <- min(hn[2,1], XUnique[NUnique])
    hn[3,1] <- min(hn[3,1], max(abs(XUnique[1]), XUnique[NUnique]))
    hn[4,1] <- min(hn[4,1], max(abs(XUnique[1]), XUnique[NUnique]))

    # nLocalMin check
    if (nLocalMin > 0) {
      hlMin <- sort(abs(X[X < 0]),decreasing=FALSE)[min(Nl, nLocalMin)]
      hrMin <- (X[X >= 0])[min(Nr, nLocalMin)]
      hn[1,1] <- max(hn[1,1], hlMin)
      hn[2,1] <- max(hn[2,1], hrMin)
      hn[3,1] <- max(hn[3,1], hlMin, hrMin)
      hn[4,1] <- max(hn[4,1], hlMin, hrMin)
    }

    # nUniqueMin check
    if (nUniqueMin > 0) {
      hlMin <- sort(abs(XUnique[XUnique < 0]),decreasing=FALSE)[min(NlUnique, nUniqueMin)]
      hrMin <- (XUnique[XUnique >= 0])[min(NrUnique, nUniqueMin)]
      hn[1,1] <- max(hn[1,1], hlMin)
      hn[2,1] <- max(hn[2,1], hrMin)
      hn[3,1] <- max(hn[3,1], hlMin, hrMin)
      hn[4,1] <- max(hn[4,1], hlMin, hrMin)
    }
  }

  result <- list(h=hn, N=list(full=N, left=Nl, right=Nr),
                 opt=list(fitselect=fitselect, kernel=kernel, vce=vce, c=c, p=p,
                          regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin,
                          massPoints=massPoints, masspoints_flag=masspoints_flag),
                 X_min      =list(left=min(X[X<0])+c, right=min(X[X>=0])+c),
                 X_max      =list(left=max(X[X<0])+c, right=max(X[X>=0])+c))

  class(result) <- "CJMrdbwdensity"
  return(result)
}

################################################################################
#' Internal function.
#'
#' @param object Class \code{CJMrdbwdensity} objects.
#'
#' @keywords internal
#' @export
summary.CJMrdbwdensity <- function(object, ...) {
  x <- object
  cat("\nBandwidth selection for manipulation testing.\n")
  cat("\n")

  cat(paste(format("Number of obs =", width=20), toString(x$N$full), sep="")); cat("\n")
  cat(paste(format("Model =", width=20), x$opt$fitselect, sep="")); cat("\n")
  cat(paste(format("Kernel =", width=20), x$opt$kernel, sep="")); cat("\n")
  cat(paste(format("VCE method =", width=20), x$opt$vce, sep="")); cat("\n")
  cat("\n")

  cat(paste(format(paste("Cutoff c = ", toString(round(x$opt$c, 3)), sep=""), width=20), format("Left of c", width=20), format("Right of c", width=20), sep="")); cat("\n")
  cat(paste(format("Number of obs", width=20), format(toString(x$N$left), width=20), format(toString(x$N$right), width=20), sep="")); cat("\n")
  cat(paste(format("Min Running var.", width=20), format(toString(round(x$X_min$left, 3)), width=20), format(toString(round(x$X_min$right, 3)), width=20), sep="")); cat("\n")
  cat(paste(format("Max Running var.", width=20), format(toString(round(x$X_max$left, 3)), width=20), format(toString(round(x$X_max$right, 3)), width=20), sep="")); cat("\n")
  cat(paste(format("Order est. (p)", width=20), format(toString(x$opt$p), width=20), format(toString(x$opt$p), width=20), sep="")); cat("\n")
  cat("\n")

  cat(paste(format("Target", width=20), format("Bandwidth", width=20), format("Variance", width=20), format("Bias^2", width=20), sep="")); cat("\n")
  cat(paste(format("left density", width=20), format(toString(round(x$h[1,1], 4)), width=20), format(toString(round(x$h[1,2], 4)), width=20), format(toString(round(x$h[1,3], 4)), width=20), sep="")); cat("\n")
  cat(paste(format("right density", width=20), format(toString(round(x$h[2,1], 4)), width=20), format(toString(round(x$h[2,2], 4)), width=20), format(toString(round(x$h[2,3], 4)), width=20), sep="")); cat("\n")
  cat(paste(format("diff. densities", width=20), format(toString(round(x$h[3,1], 4)), width=20), format(toString(round(x$h[3,2], 4)), width=20), format(toString(round(x$h[3,3], 4)), width=20), sep="")); cat("\n")
  cat(paste(format("sum densities", width=20), format(toString(round(x$h[4,1], 4)), width=20), format(toString(round(x$h[4,2], 4)), width=20), format(toString(round(x$h[4,3], 4)), width=20), sep="")); cat("\n")
  cat("\n")

  if (x$opt$masspoints_flag) {
    warning("There are repeated observations. Point estimates and standard errors have been adjusted. Use option massPoints=FALSE to suppress this feature.")
  }
}

################################################################################
#' Internal function.
#'
#' @param x Class \code{CJMrdbwdensity} objects.
#'
#' @keywords internal
#' @export
print.CJMrdbwdensity <- function(x, ...) {
  cat("Call:\n")
  cat("rdbwdensity\n")
  cat("Sample size:\ ", x$N$full, ". ", "Cutoff: ", x$opt$c, ".\n", sep="")
  cat("Model:\ ", x$opt$fitselect, ". ", "Kernel: ", x$opt$kernel, ". ", "VCE: ", x$opt$vce, sep="")
}

