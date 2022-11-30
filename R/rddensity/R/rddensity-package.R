################################################################################
#' @title rddensity: Manipulation Testing Based on Density Discontinuity
#'
#' @description Density discontinuity testing (a.k.a. manipulation testing) is commonly
#'   employed in regression discontinuity designs and other program evaluation settings
#'   to detect perfect self-selection (manipulation) around a cutoff where
#'   treatment/policy assignment changes.
#'
#' This package implements manipulation testing procedures using the
#'   local polynomial density estimators proposed in Cattaneo, Jansson and Ma (2020),
#'   and implements graphical procedures with valid confidence bands using the results
#'   in Cattaneo, Jansson and Ma (2022a,b).  In addition, this package provides complementary
#'   manipulation testing based on finite sample exact binomial testing following the
#'   esults in Cattaneo, Frandsen and Titiunik (2015) and Cattaneo, Frandsen and
#'   Vazquez-Bare (2017).
#'
#' A companion \code{Stata} package is described in Cattaneo, Jansson and Ma (2018).
#'
#' Commands: \code{\link{rddensity}} for manipulation (density discontinuity) testing.
#'   \code{\link{rdbwdensity}} for data-driven bandwidth selection, and
#'   \code{\link{rdplotdensity}} for density plots.
#'
#' Related Stata and R packages useful for inference in regression discontinuity (RD)
#'   designs are described in the website: \url{https://rdpackages.github.io/}.
#'
#' @author
#' Matias D. Cattaneo, Princeton University  \email{cattaneo@princeton.edu}.
#'
#' Michael Jansson, University of California Berkeley.  \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of California San Diego. \email{x1ma@ucsd.edu}.
#'
#' @references
#' Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2018. \href{https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2018_JASA.pdf}{On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference}. \emph{Journal of the American Statistical Association} 113(522): 767-779. \doi{10.1080/01621459.2017.1285776}
#'
#' Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2022. \href{https://cattaneo.princeton.edu/papers/Calonico-Cattaneo-Farrell_2022_Bernoulli.pdf}{Coverage Error Optimal Confidence Intervals for Local Polynomial Regression}. \emph{Bernoulli}, 28(4): 2998-3022. \doi{10.3150/21-BEJ1445}
#'
#' Cattaneo, M. D., B. Frandsen, and R. Titiunik. 2015. \href{https://rdpackages.github.io/references/Cattaneo-Frandsen-Titiunik_2015_JCI.pdf}{Randomization Inference in the Regression Discontinuity Design: An Application to the Study of Party Advantages in the U.S. Senate.} \emph{Journal of Causal Inference} 3(1): 1-24. \doi{10.1515/jci-2013-0010}
#'
#' Cattaneo, M. D., M. Jansson, and X. Ma. 2018. \href{https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2018_Stata.pdf}{Manipulation Testing based on Density Discontinuity}. \emph{Stata Journal} 18(1): 234-261. \doi{10.1177/1536867X1801800115}
#'
#' Cattaneo, M. D., M. Jansson, and X. Ma. 2020. \href{https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA.pdf}{Simple Local Polynomial Density Estimators}. \emph{Journal of the American Statistical Association}, 115(531): 1449-1455. \doi{10.1080/01621459.2019.1635480}
#'
#' Cattaneo, M. D., M. Jansson, and X. Ma. 2022a. \href{https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2022_JoE.pdf}{Local Regression Distribution Estimators}. \emph{Journal of Econometrics}, forthcoming. \doi{10.1016/j.jeconom.2021.01.006}
#'
#' Cattaneo, M. D., M. Jansson, and X. Ma. 2022b. \href{https://nppackages.github.io/references/Cattaneo-Jansson-Ma_2022_JSS.pdf}{lpdensity: Local Polynomial Density Estimation and Inference}. \emph{Journal of Statistical Software}, 101(2), 1–25. \doi{10.18637/jss.v101.i02}
#'
#' Cattaneo, M. D., R. Titiunik and G. Vazquez-Bare. 2017. \href{https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2017_JPAM.pdf}{Comparing Inference Approaches for RD Designs: A Reexamination of the Effect of Head Start on Child Mortality}. \emph{Journal of Policy Analysis and Management} 36(3): 643-681. \doi{10.1002/pam.21985}
#'
#' McCrary, J. 2008. Manipulation of the Running Variable in the Regression Discontinuity Design: A Density Test. \emph{Journal of Econometrics} 142(2): 698-714. \doi{10.1016/j.jeconom.2007.05.005}
#'
#' @importFrom graphics text
#' @importFrom stats dnorm
#' @importFrom stats integrate
#' @importFrom stats median
#' @importFrom stats pnorm
#' @importFrom stats sd
#' @importFrom stats quantile
#' @importFrom stats binom.test
#' @import lpdensity
#' @import ggplot2
#'
#' @aliases rddensity-package
"_PACKAGE"
