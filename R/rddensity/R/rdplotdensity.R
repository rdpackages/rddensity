################################################################################
#' @title Density Plotting for Manipulation Testing
#'
#' @description \code{rdplotdensity} constructs density plots. It is based on the
#'   local polynomial density estimator proposed in Cattaneo, Jansson and Ma (2020, 2023).
#'   A companion \code{Stata} package is described in Cattaneo, Jansson and Ma (2018).
#'
#' Companion command: \code{\link{rddensity}} for manipulation (density discontinuity) testing.
#'
#' Related Stata and R packages useful for inference in regression discontinuity (RD)
#'   designs are described in the website: \url{https://rdpackages.github.io/}.
#'
#' @details
#' Bias correction is only used for the construction of confidence intervals/bands, but not for point
#' estimation. The point estimates, denoted by \code{f_p}, are constructed using local polynomial estimates of order
#' \code{p}, while the centering of the confidence intervals/bands, denoted by \code{f_q}, are constructed using local
#' polynomial estimates of order \code{q}. The confidence intervals/bands take the form:
#' \code{[f_q - cv * SE(f_q) , f_q + cv * SE(f_q)]}, where \code{cv} denotes the appropriate critical value and
#' \code{SE(f_q)} denotes a standard error estimate
#' for the centering of the confidence interval/band. As a result, the confidence intervals/bands may not be
#' centered at the point estimates because they have been bias-corrected. Setting \code{q} and \code{p} to be equal
#' results on centered at the point estimate confidence intervals/bands, but requires undersmoothing for valid
#' inference (i.e., (I)MSE-optimal bandwdith for the density point estimator cannot be used).  Hence the bandwidth
#' would need to be specified manually when \code{q=p}, and the point estimates will not be (I)MSE optimal. See
#' Cattaneo, Jansson and Ma (2022, 2023) for details, and also Calonico, Cattaneo, and Farrell (2018, 2022) for
#' robust bias correction methods.
#'
#' Sometimes the density point estimates may lie outside of the confidence intervals/bands, which can happen if
#' the underlying distribution exhibits high curvature at some evaluation point(s).  One possible solution in this
#' case is to increase the polynomial order \code{p} or to employ a smaller bandwidth.
#'
#' @param rdd Object returned by \code{\link{rddensity}}
#' @param X Numeric vector or one dimensional matrix/data frame, the running variable.
#' @param plotRange Numeric, specifies the lower and upper bound of the plotting region.  Default is
#'   \code{[c-3*hl,c+3*hr]} (three bandwidths around the cutoff).
#' @param plotN Numeric, specifies the number of grid points used for plotting on the two sides of the cutoff.
#'   Default is \code{c(10,10)} (i.e., 10 points are used on each side).
#' @param plotGrid String, specifies how the grid points are positioned.  Options are \code{es} (evenly spaced)
#'   and \code{qs} (quantile spaced).
#' @param alpha Numeric scalar between 0 and 1, the significance level for plotting
#'   confidence regions. If more than one is provided, they will be applied to the two sides
#'   accordingly.
#' @param type String, one of \code{"line"} (default), \code{"points"} or \code{"both"}, how
#'   the point estimates are plotted. If more than one is provided, they will be applied to the two sides
#'   accordingly.
#' @param lty Line type for point estimates, only effective if \code{type} is \code{"line"} or
#'   \code{"both"}. \code{1} for solid line, \code{2} for dashed line, \code{3} for dotted line.
#'   For other options, see the instructions for \code{ggplot2} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to the two sides accordingly.
#' @param lwd Line width for point estimates, only effective if \code{type} is \code{"line"} or
#'   \code{"both"}. Should be strictly positive. For other options, see the instructions for
#'   \code{ggplot2} or \code{\link{par}}. If more than one is provided, they will be applied
#'   to the two sides accordingly.
#' @param lcol Line color for point estimates, only effective if \code{type} is \code{"line"} or
#'   \code{"both"}. \code{1} for black, \code{2} for red, \code{3} for green, \code{4} for blue.
#'   For other options, see the instructions for \code{ggplot2} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to the two sides
#'   accordingly.
#' @param pty Scatter plot type for point estimates, only effective if \code{type} is \code{"points"} or
#'   \code{"both"}. For options, see the instructions for \code{ggplot2} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to the two sides
#'   accordingly.
#' @param pwd Scatter plot size for point estimates, only effective if \code{type} is \code{"points"} or
#'   \code{"both"}. Should be strictly positive. If more than one is provided, they will be applied to the two sides
#'   accordingly.
#' @param pcol Scatter plot color for point estimates, only effective if \code{type} is \code{"points"} or
#'   \code{"both"}. \code{1} for black, \code{2} for red, \code{3}
#'   for green, \code{4} for blue.
#'   For other options, see the instructions for \code{ggplot2} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to the two sides
#'   accordingly.
#' @param CItype String, one of \code{"region"} (shaded region, default), \code{"line"} (dashed lines),
#'   \code{"ebar"} (error bars), \code{"all"} (all of the previous) or \code{"none"} (no confidence region),
#'   how the confidence region should be plotted. If more than one is provided, they will be applied to the two sides
#'   accordingly.
#' @param CIuniform \code{TRUE} or \code{FALSE} (default), plotting either pointwise confidence intervals (\code{FALSE}) or
#'   uniform confidence bands (\code{TRUE}).
#' @param CIsimul Positive integer, the number of simulations used to construct critical values (default is 2000). This
#'   option is ignored if \code{CIuniform=FALSE}.
#' @param CIshade Numeric, opaqueness of the confidence region, should be between 0 (transparent) and
#'   1. Default is 0.2. If more than one is provided, they will be applied to the two sides
#'   accordingly.
#' @param CIcol Color of the confidence region. \code{1} for black, \code{2} for red, \code{3}
#'   for green, \code{4} for blue.
#'   For other options, see the instructions for \code{ggplot2} or \code{\link{par}}. If
#'   more than one is provided, they will be applied to the two sides
#'   accordingly.
#' @param bwselect String, the method for data-driven bandwidth selection. Available options
#'   are (1) \code{"mse-dpi"} (mean squared error-optimal bandwidth selected for each grid point);
#'   (2) \code{"imse-dpi"} (integrated MSE-optimal bandwidth, common for all grid points);
#'   (3) \code{"mse-rot"} (rule-of-thumb bandwidth with Gaussian reference model);
#'   and (4) \code{"imse-rot"} (integrated rule-of-thumb bandwidth with Gaussian reference model).
#'   If omitted, bandwidths returned by \code{rddensity} will be used.
#' @param hist \code{TRUE} (default) or \code{FALSE}, whether adding a histogram to the background.
#' @param histBreaks Numeric vector, giving the breakpoints between histogram cells.
#' @param histFillCol Color of the histogram cells.
#' @param histFillShade Opaqueness of the histogram cells, should be between 0 (transparent) and
#'   1. Default is 0.2.
#' @param histLineCol Color of the histogram lines.
#' @param title,xlabel,ylabel Strings, title of the plot and labels for x- and y-axis.
#' @param legendTitle String, title of legend.
#' @param legendGroups String Vector, group names used in legend.
#' @param noPlot No density plot will be generated if set to \code{TRUE}.
#'
#' @return
#' \item{Estl, Estr}{Matrices containing estimation results:
#'   (1) \code{grid} (grid points),
#'   (2) \code{bw} (bandwidths),
#'   (3) \code{nh} (number of observations in each local neighborhood),
#'   (4) \code{nhu} (number of unique observations in each local neighborhood),
#'   (5) \code{f_p} (point estimates with p-th order local polynomial),
#'   (6) \code{f_q} (point estimates with q-th order local polynomial, only if option \code{q} is nonzero),
#'   (7) \code{se_p} (standard error corresponding to \code{f_p}), and (8) \code{se_q} (standard error
#'   corresponding to \code{f_q}).
#'   Variance-covariance matrix corresponding to \code{f_p}.
#'   Variance-covariance matrix corresponding to \code{f_q}.
#'   A list containing options passed to the function.}
#' \item{Estplot}{A stadnard \code{ggplot} object is returned, hence can be used for further customization.}
#'
#' @author
#' Matias D. Cattaneo, Princeton University  \email{cattaneo@princeton.edu}.
#'
#' Michael Jansson, University of California Berkeley.  \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @references
#' Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2018. On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference. \emph{Journal of the American Statistical Association} 113(522): 767-779. \doi{10.1080/01621459.2017.1285776}
#'
#' Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2022. Coverage Error Optimal Confidence Intervals for Local Polynomial Regression. \emph{Bernoulli}, 28(4): 2998-3022. \doi{10.3150/21-BEJ1445}
#'
#' Cattaneo, M. D., M. Jansson, and X. Ma. 2018. Manipulation Testing based on Density Discontinuity. \emph{Stata Journal} 18(1): 234-261. \doi{10.1177/1536867X1801800115}
#'
#' Cattaneo, M. D., M. Jansson, and X. Ma. 2020. Simple Local Polynomial Density Estimators. \emph{Journal of the American Statistical Association}, 115(531): 1449-1455. \doi{10.1080/01621459.2019.1635480}
#'
#' Cattaneo, M. D., M. Jansson, and X. Ma. 2022. lpdensity: Local Polynomial Density Estimation and Inference. \emph{Journal of Statistical Software}, 101(2): 1–25. \doi{10.18637/jss.v101.i02}
#'
#' Cattaneo, M. D., M. Jansson, and X. Ma. 2023. Local Regression Distribution Estimators. \emph{Journal of Econometrics}, 240(2): 105074. \doi{10.1016/j.jeconom.2021.01.006}
#'
#' @seealso \code{\link{rddensity}}
#'
#' @examples
#' # Generate a random sample with a density discontinuity at 0
#' set.seed(42)
#' x <- rnorm(2000, mean = -0.5)
#' x[x > 0] <- x[x > 0] * 2
#'
#' # Estimation
#' rdd <- rddensity(X = x)
#' summary(rdd)
#'
#' # Density plot (from -2 to 2 with 25 evaluation points at each side)
#' plot1 <- rdplotdensity(rdd, x, plotRange = c(-2, 2), plotN = 25)
#'
#' # Plotting a uniform confidence band
#' set.seed(42) # fix the seed for simulating critical values
#' plot3 <- rdplotdensity(rdd, x, plotRange = c(-2, 2), plotN = 25, CIuniform = TRUE)
#'
#' @export
rdplotdensity <- function(rdd, X, plotRange = NULL, plotN = 10, plotGrid = c("es", "qs"),
                          alpha = 0.05,
                          type = NULL,
                          lty = NULL, lwd = NULL, lcol = NULL,
                          pty = NULL, pwd = NULL, pcol = NULL,
                          CItype = NULL,
                          CIuniform=FALSE, CIsimul=2000, CIshade = NULL, CIcol = NULL,
                          bwselect = NULL,
                          hist=TRUE, histBreaks=NULL, histFillCol=3, histFillShade=0.2, histLineCol="white",
                          title = "", xlabel = "", ylabel = "", legendTitle = NULL, legendGroups = NULL,
                          noPlot = FALSE){

  # obtain options from rddensity result
  c       <- rdd$opt$c
  p       <- rdd$opt$p
  q       <- rdd$opt$q
  hl      <- rdd$h$left
  hr      <- rdd$h$right
  kernel  <- rdd$opt$kernel
  regularize <- rdd$opt$regularize
  nLocalMin  <- rdd$opt$nLocalMin
  nUniqueMin <- rdd$opt$nUniqueMin
  massPoints <- rdd$opt$massPoints

  ################################################################################
  # missing value handling
  ################################################################################
  X <- as.vector(X)
  if (any(is.na(X))) {
    warning(paste(sum(is.na(X)), " missing ", switch((sum(is.na(X))>1)+1, "observation is", "observations are"), " ignored.\n", sep=""))
    X <- X[!is.na(X)]
  }

  # check grid specifications
  if (length(plotRange) == 0) {
    plotRange <- c( max(min(X), c - 3*hl), min(max(X), c + 3 * hr) )
  } else if (length(plotRange) != 2) {
    stop("Plot range incorrectly specified.\n")
  } else if (plotRange[1] >= c | plotRange[2] <= c) {
    stop("Plot range incorrectly specified.\n")
  }

  if (length(plotN) == 0) {
    plotN <- c(10, 10)
  } else if (length(plotN) == 1) {
    plotN <- c(plotN, plotN)
  } else if (length(plotN) > 2) {
    stop("Number of grid points incorrectly specified.\n")
  }
  if (plotN[1] <=1 | plotN[2] <=1) {
    stop("Number of grid points incorrectly specified.\n")
  }

  if (length(plotGrid) == 0) {
    plotGrid <- "es"
  } else {
    plotGrid <- plotGrid[1]
  }
  if (!plotGrid%in%c("es", "qs")) {
    stop("Grid specification invalid.\n")
  }

  if (hist & is.null(histBreaks)) {
    temp_hist_n_l <- sum(X >= plotRange[1] & X < c)
    temp_hist_n_l <- ceiling(min(sqrt(temp_hist_n_l), 10 * log(temp_hist_n_l)/log(10)))
    temp_hist_n_r <- sum(X <= plotRange[2] & X >= c)
    temp_hist_n_r <- ceiling(min(sqrt(temp_hist_n_r), 10 * log(temp_hist_n_r)/log(10)))

    # implemented in version 2.6, so that obs  (mass point) on the cutoff are included in the bin right of c
    histBreaks <- (plotRange[1] - plotRange[2]) / (1e8) +
      c(seq(plotRange[1], c, length.out = temp_hist_n_l+1), seq(c, plotRange[2], length.out = temp_hist_n_r+1)[2:(temp_hist_n_r+1)])
  }

  # some preparation
  #scalel <- (sum(X <= c) - 1) / (length(X) - 1)
  #scaler <- (sum(X >= c) - 1) / (length(X) - 1)

  # implemented in version 2.6
  scalel <- (sum(X <  c) ) / (length(X) - 1 )
  scaler <- (sum(X >= c) ) / (length(X) - 1 )

  if (plotGrid == "es") {
    gridl <- seq(plotRange[1], c, length.out=plotN[1])
    gridl[plotN[1]] <- c
    gridr <- seq(c, plotRange[2], length.out=plotN[2])
    gridr[1] <- c
  } else {
    gridl <- seq(mean(X <= plotRange[1]), mean(X < c), length.out=plotN[1])
    gridl <- quantile(X, gridl)
    gridr <- seq(mean(X < c), mean(X <= plotRange[2]), length.out=plotN[2])
    gridr <- quantile(X, gridr)
    gridl[plotN[1]] <- c
    gridr[1] <- c
  }

  # call lpdensity
  if (!is.null(bwselect)) {
    if (bwselect%in%c("mse-dpi", "imse-dpi", "mse-rot", "imse-rot")) {
      Estl <- lpdensity(data=X[X< c], grid=gridl, bwselect=bwselect, p=p, q=q, v=1, kernel=kernel, scale=scalel,
                        regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin, massPoints=massPoints)
      Estr <- lpdensity(data=X[X>=c], grid=gridr, bwselect=bwselect, p=p, q=q, v=1, kernel=kernel, scale=scaler,
                        regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin, massPoints=massPoints)
    } else {
      stop("Option bwselect incorrectly specified.\n")
    }
  } else {
    Estl <- lpdensity(data=X[X< c], grid=gridl, bw=hl, p=p, q=q, v=1, kernel=kernel, scale=scalel,
                      regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin, massPoints=massPoints)
    Estr <- lpdensity(data=X[X>=c], grid=gridr, bw=hr, p=p, q=q, v=1, kernel=kernel, scale=scaler,
                      regularize=regularize, nLocalMin=nLocalMin, nUniqueMin=nUniqueMin, massPoints=massPoints)
  }


  # call lpdensity.plot
  Estplot <- lpdensity.plot(Estl, Estr, alpha = alpha,
                            type = type,
                            lty = lty, lwd = lwd,
                            lcol = lcol, pty = pty, pwd = pwd, pcol = pcol,
                            CItype = CItype, CIuniform=CIuniform, CIsimul=CIsimul, CIshade = CIshade, CIcol = CIcol,
                            hist=hist, histData=X, histBreaks=histBreaks, histFillCol=histFillCol, histFillShade=histFillShade, histLineCol=histLineCol,
                            title = title, xlabel = xlabel, ylabel = ylabel, legendTitle = legendTitle, legendGroups = legendGroups) +
    theme(legend.position = "none")

  if (!noPlot) { print(Estplot) }

  return(list(Estl=Estl, Estr=Estr, Estplot=Estplot))
}
