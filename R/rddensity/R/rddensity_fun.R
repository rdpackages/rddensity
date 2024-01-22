################################################################################
#' Internal function.
#'
#' Find unique elements and their frequencies in a numeric vector. This function
#'   has a similar performance as the built-in R function \code{unique}.
#'
#' @param x Numeric vector, already sorted in ascending order.
#'
#' @return
#' \item{unique}{A vector containing unique elements in \code{x}.}
#' \item{freq}{The frequency of each element in \code{x}.}
#' \item{index}{The last occurrence of each element in \code{x}.}
#'
#' @keywords internal
rddensityUnique <- function(x) {
  n <- length(x)
  # if x has one or no element
  if (n == 0) return(list(unique=NULL, freq=c(), indexFirst=c(), indexLast=c()))
  if (n == 1) return(list(unique=x, freq=1, indexFirst=1, indexLast=1))

  # else
  uniqueIndex <- c(x[2:n] != x[1:(n-1)], TRUE)
  unique <- x[uniqueIndex]
  nUnique <- length(unique)

  # all are distinct
  if (nUnique == n) return(list(unique=unique, freq=rep(1,length(x)), indexFirst=1:n, indexLast=1:n))
  # all are the same
  if (nUnique == 1) return(list(unique=unique, freq=n, indexFirst=1, indexLast=n))

  # otherwise
  freq <- (cumsum(!uniqueIndex))[uniqueIndex]
  freq <- freq - c(0, freq[1:(nUnique-1)]) + 1

  return(list(unique=unique, freq=freq, indexFirst=c(1, ((1:n)[uniqueIndex]+1)[1:(nUnique-1)]), indexLast=(1:n)[uniqueIndex]))
}

#rddensityUnique(1:10)
#rddensityUnique(c(1,1,2,3,3,3, 4,4,5,5,5,5,5,5))

################################################################################
#' Internal function, generate matrices.
#'
#' \code{Sgenerate} generates a matrix.
#'
#' This is an internal function, and should not be called by users.
#'
#' @param p Integer, polynomial order.
#' @param low,up Numeric, lower and upper limits of integration, -1 and 1 by
#'   default.
#' @param kernel String, the kernel function, can be \code{triangular} (default),
#'   \code{uniform} or \code{epanechnikov}.
#'
#' @return Returns a matrix.
#'
#' @keywords internal
Sgenerate <- function(p, low=-1, up=1, kernel="triangular") {
  popwarning <- FALSE
  S <- matrix(rep(0, (p+1)^2), ncol=(p+1))
  for (i in 1:(p+1)) {
    for (j in 1:(p+1)) {
      if (kernel == "uniform") {
        integrand <- function(x) { x^(i+j-2)*0.5 }
      } else if (kernel == "epanechnikov") {
        integrand <- function(x) { x^(i+j-2)*0.75*(1-x^2) }
      } else {
        integrand <- function(x) { x^(i+j-2)*(1-abs(x)) }
      }
      S[i,j] <- (integrate(integrand, lower=low, upper=up)$value)
    }
  }
  if (popwarning) {warning(text)}
  return(S)
}

################################################################################
#' Internal function, generate matrices.
#'
#' \code{Sminusgenerate} generates a matrix.
#'
#' This is an internal function, and should not be called by users.
#'
#' @param p Integer, polynomial order.
#' @param kernel String, the kernel function, can be \code{triangular} (default),
#'   \code{uniform} or \code{epanechnikov}.
#'
#' @return Returns a matrix.
#'
#' @keywords internal
Sminusgenerate <- function(p, kernel="triangular") {
  S <- Sgenerate(p, low=-1, up=0, kernel=kernel)
  temp <- matrix(0, ncol=p+2, nrow=p+2)
  temp[1:2, 1:2] <- S[1:2, 1:2]
  if (p>1) {
    temp[1:2, 4:(p+2)] <- S[1:2, 3:(p+1)]
    temp[4:(p+2), 1:2] <- S[3:(p+1), 1:2]
    temp[4:(p+2), 4:(p+2)] <- S[3:(p+1), 3:(p+1)]
  }
  return(temp)
}

################################################################################
#' Internal function, generate matrices.
#'
#' \code{Splusgenerate} generates a matrix.
#'
#' This is an internal function, and should not be called by users.
#'
#' @param p Integer, polynomial order.
#' @param kernel String, the kernel function, can be \code{triangular} (default),
#'   \code{uniform} or \code{epanechnikov}.
#'
#' @return Returns a matrix.
#'
#' @keywords internal
Splusgenerate <- function(p, kernel="triangular") {
  S <- Sgenerate(p, low=0, up=1, kernel=kernel)
  temp <- matrix(0, ncol=p+2, nrow=p+2)
  temp[1, 1] <- S[1, 1]
  temp[1, 3:(p+2)] <- S[1, 2:(p+1)]
  temp[3:(p+2), 1] <- S[2:(p+1), 1]
  temp[3:(p+2), 3:(p+2)] <- S[2:(p+1), 2:(p+1)]
  return(temp)
}

################################################################################
#' Internal function, generate matrices.
#'
#' \code{Cgenerate} generates a matrix.
#'
#' This is an internal function, and should not be called by users.
#'
#' @param k,p Integer, polynomial order.
#' @param kernel String, the kernel function, can be \code{triangular} (default),
#'   \code{uniform} or \code{epanechnikov}.
#'
#' @return Returns a matrix.
#'
#' @keywords internal
Cgenerate <- function(k, p, low=-1, up=1, kernel="triangular") {
  popwarning <- FALSE
  C <- matrix(rep(0, (p+1)), ncol=1)
  for (i in 1:(p+1)) {
    if (kernel == "uniform") {
      integrand <- function(x) { x^(i+k-1)*0.5 }
    } else if (kernel == "epanechnikov") {
      integrand <- function(x) { x^(i+k-1)*0.75*(1-x^2) }
    }
    else {
      integrand <- function(x) { x^(i+k-1)*(1-abs(x)) }
    }
    C[i,1] <- (integrate(integrand, lower=low, upper=up)$value)
  }
  if (popwarning) {warning(text)}
  return(C)
}

################################################################################
#' Internal function, generate matrices.
#'
#' \code{Cminusgenerate} generates a matrix.
#'
#' This is an internal function, and should not be called by users.
#'
#' @param k,p Integer, polynomial order.
#' @param kernel String, the kernel function, can be \code{triangular} (default),
#'   \code{uniform} or \code{epanechnikov}.
#'
#' @return Returns a matrix.
#'
#' @keywords internal
Cminusgenerate <- function(k, p, kernel="triangular") {
  C <- Cgenerate(k=k, p=p, kernel=kernel, low=-1, up=0)
  temp <- matrix(0, ncol=1, nrow=p+2)
  temp[1:2] <- C[1:2]
  if(p>1) {
    temp[4:(p+2)] <- C[3:(p+1)]
  }
  return (temp)
}

################################################################################
#' Internal function, generate matrices.
#'
#' \code{Cminusgenerate} generates a matrix.
#'
#' This is an internal function, and should not be called by users.
#'
#' @param k,p Integer, polynomial order.
#' @param kernel String, the kernel function, can be \code{triangular} (default),
#'   \code{uniform} or \code{epanechnikov}.
#'
#' @return Returns a matrix.
#'
#' @keywords internal
Cplusgenerate <- function(k, p, kernel="triangular") {
  C <- Cgenerate(k=k, p=p, kernel=kernel, low=0, up=1)
  temp <- matrix(0, ncol=1, nrow=p+2)
  temp[1] <- C[1]
  temp[3:(p+2)] <- C[2:(p+1)]
  return (temp)
}

################################################################################
#' Internal function, generate matrices.
#'
#' \code{Ggenerate} generates a matrix.
#'
#' This is an internal function, and should not be called by users.
#'
#' @param p Integer, polynomial order.
#' @param low,up Numeric, lower and upper limits of integration, -1 and 1 by
#'   default.
#' @param kernel String, the kernel function, can be \code{triangular} (default),
#'   \code{uniform} or \code{epanechnikov}.
#'
#' @return Returns a matrix.
#'
#' @keywords internal
Ggenerate <- function(p, low=-1, up=1, kernel="triangular") {
  popwarning <- FALSE
  G <- matrix(rep(0, (p+1)^2), ncol=(p+1))
  for (i in 1:(p+1)) {
    for (j in 1:(p+1)) {
      if (kernel == "uniform") {
        G[i,j] <- integrate(function(y) {
          sapply(y, function(y) {
            integrate(function(x) x^i * y^(j-1)*0.25, low, y)$value
          })
        }, low, up)$value +
          integrate(function(y) {
            sapply(y, function(y) {
              integrate(function(x) x^(i-1) * y^j*0.25, y, up)$value
            })
          }, low, up)$value
      } else if (kernel == "epanechnikov") {
        G[i,j] <- integrate(function(y) {
          sapply(y, function(y) {
            integrate(function(x) x^i * y^(j-1) * 0.75^2 *
                        (1-x^2) * (1-y^2), low, y)$value
          })
        }, low, up)$value +
          integrate(function(y) {
            sapply(y, function(y) {
              integrate(function(x) x^(i-1) * y^j * 0.75^2 *
                          (1-x^2) * (1-y^2), y, up)$value
            })
          }, low, up)$value
      } else {
        G[i,j] <- integrate(function(y) {
          sapply(y, function(y) {
            integrate(function(x) x^i * y^(j-1) *
                        (1-abs(x)) * (1-abs(y)), low, y)$value
          })
        }, low, up)$value +
          integrate(function(y) {
            sapply(y, function(y) {
              integrate(function(x) x^(i-1) * y^j *
                          (1-abs(x)) * (1-abs(y)), y, up)$value
            })
          }, low, up)$value
      }
    }
  }
  if (popwarning) {warning(text)}
  return(G)
}

################################################################################
#' Internal function, generate matrices.
#'
#' \code{Gminusgenerate} generates a matrix.
#'
#' This is an internal function, and should not be called by users.
#'
#' @param p Integer, polynomial order.
#' @param kernel String, the kernel function, can be \code{triangular} (default),
#'   \code{uniform} or \code{epanechnikov}.
#'
#' @return Returns a matrix.
#'
#' @keywords internal
Gminusgenerate <- function(p, kernel="triangular") {
  G <- Ggenerate(p, low=-1, up=0, kernel=kernel)
  temp <- matrix(0, ncol=p+2, nrow=p+2)
  temp[1:2, 1:2] <- G[1:2, 1:2]
  if (p>1) {
    temp[1:2, 4:(p+2)] <- G[1:2, 3:(p+1)]
    temp[4:(p+2), 1:2] <- G[3:(p+1), 1:2]
    temp[4:(p+2), 4:(p+2)] <- G[3:(p+1), 3:(p+1)]
  }
  return(temp)
}

################################################################################
#' Internal function, generate matrices.
#'
#' \code{Gplusgenerate} generates a matrix.
#'
#' This is an internal function, and should not be called by users.
#'
#' @param p Integer, polynomial order.
#' @param kernel String, the kernel function, can be \code{triangular} (default),
#'   \code{uniform} or \code{epanechnikov}.
#'
#' @return Returns a matrix.
#'
#' @keywords internal
Gplusgenerate <- function(p, kernel="triangular") {
  G <- Ggenerate(p, low=0, up=1, kernel=kernel)
  temp <- matrix(0, ncol=p+2, nrow=p+2)
  temp[1, 1] <- G[1, 1]
  temp[1, 3:(p+2)] <- G[1, 2:(p+1)]
  temp[3:(p+2), 1] <- G[2:(p+1), 1]
  temp[3:(p+2), 3:(p+2)] <- G[2:(p+1), 2:(p+1)]
  return(temp)
}

################################################################################
#' Internal function, generate matrices.
#'
#' \code{Psigenerate} generates a matrix.
#'
#' This is an internal function, and should not be called by users.
#'
#' @param p Integer, polynomial order.
#'
#' @return Returns a matrix.
#'
#' @keywords internal
Psigenerate <- function(p) {
  if (p > 1) {
    temp <- c(1, 0, 0, (-1)^(2:p))
  } else {
    temp <- c(1, 0, 0)
  }
  temp <- diag(temp)
  temp[2, 3] <- temp[3, 2] <- -1
  return(temp)
}

################################################################################
#' Internal function, implements density test.
#'
#' \code{rddensity_fV} generates primitive results for the density test.
#'
#' This is an internal function, and should not be called by users.
#'   NOTE: data is assumed to be on ascending order.
#'
#' @param Y Numeric vector or data matrix, the estimated c.d.f.
#' @param X Numeric vector or data matrix, the running variable
#' @param Nl,Nr Integers, sample sizes to the left and right of cutoff.
#' @param Nlh,Nrh Integers, sample sizes to the left and right of cutoff,
#'   within bandwidth
#' @param hl,hr Numeric, bandwidth to the left and right of cutoff.
#' @param p Integer, polynomial order
#' @param s Integer, higher order derivative estimate
#' @param kernel String, the kernel function, can be \code{triangular} (default),
#'   \code{uniform} or \code{epanechnikov}.
#' @param fitselect String, the model, can be \code{restricted} or \code{unrestricted}
#' @param vce String, specifies the procedure used to compute the variance-covariance matrix estimator. Options are:
#'   \code{"plugin"} for asymptotic plug-in standard errors. \code{"jackknife"} for jackknife standard errors.
#' @param massPoints Boolean, whether whether point estimates and standard errors
#'
#' @return Returns a data frame for further use.
#'
#' @keywords internal
rddensity_fV <- function(Y, X, Nl, Nr, Nlh, Nrh, hl, hr, p, s,
                         kernel, fitselect, vce, massPoints) {

  N <- Nl + Nr
  Nh <- Nlh + Nrh
  Y <- matrix(Y, ncol=1)
  X <- matrix(X, ncol=1)

  # Construct the kernel weightings
  W <- rep(NA, Nh)
  if (kernel == "uniform") {
    W[1:Nlh] <- 1 / (2 * hl); W[(Nlh+1):Nh] <- 1 / (2 * hr)
  } else if (kernel == "triangular") {
    W[1:Nlh] <- (1 + X[1:Nlh]/hl) / hl; W[(Nlh+1):Nh] <- (1 - X[(Nlh+1):Nh]/hr) / hr
  } else {
    W[1:Nlh] <- 0.75 * (1 - (X[1:Nlh]/hl)^2) / hl; W[(Nlh+1):Nh] <- 0.75 * (1 - (X[(Nlh+1):Nh]/hr)^2) / hr
  }

  # Construct the design matrix and the bandwidth matrix
  if (fitselect == "restricted") {
    Xp <- matrix(NA, ncol=p+2, nrow=Nh)
    Xp[, 1] <- 1
    Xp[1:Nlh, 2] <- X[1:Nlh] / hl; Xp[(Nlh+1):Nh, 2] <- 0
    Xp[1:Nlh, 3] <- 0; Xp[(Nlh+1):Nh, 3] <- X[(Nlh+1):Nh] / hr
    if (p>1) {
      for (j in 4:(p+2)) {
        Xp[1:Nlh, j] <- (X[1:Nlh] / hl)^(j-2); Xp[(Nlh+1):Nh, j] <- (X[(Nlh+1):Nh] / hr)^(j-2)
      }
      v <- c(0, 1, 1, 2:p)
    } else {
      v <- c(0, 1, 1)
    }
    HpInv <- diag(1/(hl^v))
    Hp <- diag(hl^v)
  } else {
    Xp <- matrix(NA, ncol=2*p+2, nrow=Nh)
    Hp <- rep(NA, 2*p+2)
    for (j in 1:(2*p+2)) {
      if (j %% 2) {
        Xp[1:Nlh, j] <- (X[1:Nlh] / hl)^((j-1)/2); Xp[(Nlh+1):Nh, j] <- 0
        Hp[j] <- hl^((j-1)/2)
      } else {
        Xp[1:Nlh, j] <- 0; Xp[(Nlh+1):Nh, j] <- (X[(Nlh+1):Nh] / hr)^((j-2)/2)
        Hp[j] <- hr^((j-2)/2)
      }
    }
    HpInv <- diag(1/Hp)
    Hp <- diag(Hp)
  }

  out <- matrix(NA, ncol=4, nrow=4)
  colnames(out) <- c("hat", "jackknife", "plugin", "s")
  rownames(out) <- c("l", "r", "diff", "sum")

  # X'WX inverse matrix
  XpW <- sweep(Xp, MARGIN=1, STATS=W, FUN="*")
  Sinv <- try(solve(crossprod(XpW, Xp), tol=0), silent=TRUE)

  if (typeof(Sinv) == "character") {
    return(data.frame(out))
  }
  # point estimates
  b <- HpInv %*% Sinv %*% crossprod(XpW, Y)

  if (fitselect == "restricted") {
    out[1, 1] <- b[2]
    out[2, 1] <- b[3]
    out[3, 1] <- b[3] - b[2]
    out[4, 1] <- b[3] + b[2]
    out[1, 4] <- out[2, 4] <- b[s+2]
    out[3, 4] <- 0
    out[4, 4] <- 2 * out[1, 4]
  } else {
    out[1, 1] <- b[3]
    out[2, 1] <- b[4]
    out[3, 1] <- b[4] - b[3]
    out[4, 1] <- b[4] + b[3]
    out[1, 4] <- b[2*s+1]
    out[2, 4] <- b[2*s+2]
    out[3, 4] <- out[2, 4] - out[1, 4]
    out[4, 4] <- out[2, 4] + out[1, 4]
  }

  # Jackknife
  if (vce == "jackknife") {
    L <- matrix(0, nrow=dim(Xp)[1], ncol=dim(Xp)[2])
    # mass points correction
    if (massPoints) {
      XUnique     <- rddensityUnique(X)
      freqUnique  <- XUnique$freq
      indexUnique <- XUnique$indexFirst
      for (jj in 1:ncol(L)) {
        L[, jj] <- rep(((cumsum(c(0, XpW[Nh:1, jj])) / (N - 1))[Nh:1])[indexUnique], times=freqUnique)
      }
    } else {
      L[1, ] <- colSums(XpW[2:Nh, ]) / (N - 1)
      for (i in 2:Nh) {
        L[i, ] <- L[i-1, ] - XpW[i, ] / (N - 1)
      }
    }

    V = HpInv %*% Sinv %*% (t(L) %*% L) %*% Sinv %*% HpInv
    if (fitselect == "restricted") {
      out[1, 2] <- V[2, 2]
      out[2, 2] <- V[3, 3]
      out[3, 2] <- V[2, 2] + V[3, 3] - 2 * V[2,3]
      out[4, 2] <- V[2, 2] + V[3, 3] + 2 * V[2,3]
    } else {
      out[1, 2] <- V[3, 3]
      out[2, 2] <- V[4, 4]
      out[3, 2] <- V[3, 3] + V[4, 4] - 2 * V[3,4]
      out[4, 2] <- V[3, 3] + V[4, 4] + 2 * V[3,4]
    }
  }


  # plugin
  if (vce == "plugin") {
    if (fitselect=="unrestricted") {
      S <- Sgenerate(p, low=0, up=1, kernel=kernel)
      G <- Ggenerate(p, low=0, up=1, kernel=kernel)
      V <- solve(S) %*% G %*% solve(S)
      out[1, 3] <- out[1, 1] * V[2, 2] / (N * hl)
      out[2, 3] <- out[2, 1] * V[2, 2] / (N * hr)
      out[3, 3] <- out[4, 3] <- out[1, 3] + out[2, 3]
    } else {
      S <- Splusgenerate(p=p, kernel=kernel)
      G <- Gplusgenerate(p=p, kernel=kernel)
      Psi <- Psigenerate(p=p)
      Sm <- Psi %*% S %*% Psi; Gm <- Psi %*% G %*% Psi
      V <- solve(out[1, 1] * Sm + out[2, 1] * S) %*% (out[1, 1]^3 * Gm + out[2, 1]^3 * G) %*% solve(out[1, 1] * Sm + out[2, 1] * S)
      out[1, 3] <- V[2, 2] / (N * hl)
      out[2, 3] <- V[3, 3] / (N * hl)
      out[3, 3] <- (V[2, 2] + V[3, 3] - 2 * V[2,3]) / (N * hl)
      out[4, 3] <- (V[2, 2] + V[3, 3] + 2 * V[2,3]) / (N * hl)
    }
  }


  for (i in 1:4) {
    for (j in 2:3) {
      if (!is.na(out[i, j])) if (out[i, j] < 0) { out[i,j] <- NA }
    }
  }

  return(data.frame(out))
}

################################################################################
#' Internal function, calculates theoretical optimal bandwidth.
#'
#' \code{h_opt_density} calculates theoretical optimal bandwidth.
#'
#' This is an internal function, and should not be called by users.
#'
#' @param x Numeric vector or data matrix, the running variable.
#' @param p Integer, polynomial order.
#' @param N Integer, sample size.
#' @param dgp_F1,dgp_Fp1 Numeric, theoretical d.g.p.
#' @param f_low,f_up Numeric, lower and upper boundaries.
#' @param kernel String, the kernel function, can be \code{triangular} (default),
#'   \code{uniform} or \code{epanechnikov}.
#'
#' @return Returns MSE optimal bandwidth.
#'
#' @keywords internal
h_opt_density <- function(x, p, N, dgp_F1, dgp_Fp1, f_low, f_up, kernel="triangular") {

  # warning message
  if ((kernel != "triangular") & (kernel != "epanechnikov") & (kernel != "uniform")) {
    text <- paste("No kernel as ", toString(kernel),
                  ", triangular kernel (default) is used.", sep="")
    message(text)
    kernel <- "triangular"
  }
  # end of warning message

  if (x==f_low | x==f_up) {
    c_low <- 0
  } else {
    c_low <- -1
  }
  c_up <- 1

  e <- matrix(rep(0, p+1), ncol=1)
  e[2] <- 1

  S1    <- Sgenerate(p=p, low=c_low, up=c_up, kernel=kernel)
  Cp1   <- Cgenerate(k=p+1, p=p, low=c_low, up=c_up, kernel=kernel)
  G     <- Ggenerate(p=p, low=c_low, up=c_up, kernel=kernel)

  kappa <- N^(-1/(2*p+1)) * (t(e) %*% solve(S1) %*% G %*% solve(S1) %*% e)^(1/(2*p+1)) *
    (abs(t(e) %*% solve(S1) %*% Cp1))^(-2/(2*p+1)) * factorial(p+1)^(2/(2*p+1)) * (2*p)^(-1/(2*p+1))

  biassq <- (abs(dgp_Fp1))^(-2/(2*p+1))
  var1 <-  (dgp_F1)^(1/(2*p+1))

  h.opt <- kappa * biassq * var1
  return(as.numeric(h.opt))
}

################################################################################
#' Internal function, calculates theoretical optimal bandwidth.
#'
#' \code{h_opt_density_res} calculates theoretical optimal bandwidth for restricted
#'   model.
#'
#' This is an internal function, and should not be called by users.
#'
#' @param p Integer, polynomial order.
#' @param N Integer, sample size.
#' @param dgp_F1_l,dgp_F1_r,dgp_Fp1_l,dgp_Fp1_r Numeric, theoretical d.g.p.
#' @param kernel String, the kernel function, can be \code{triangular} (default),
#'   \code{uniform} or \code{epanechnikov}.
#'
#' @return Returns MSE optimal bandwidth.
#'
#' @keywords internal
h_opt_density_res <- function(p, N, dgp_F1_l, dgp_F1_r, dgp_Fp1_l, dgp_Fp1_r, kernel="triangular") {

  # warning message
  if ((kernel != "triangular") & (kernel != "epanechnikov") & (kernel != "uniform")) {
    text <- paste("No kernel as ", toString(kernel),
                  ", triangular kernel (default) is used.", sep="")
    message(text)
    kernel <- "triangular"
  }
  # end of warning message

  e <- matrix(rep(0, p+2), ncol=1)
  if (p %% 2) {
    e[2] <- -1; e[3] <- 1
  } else {
    e[2] <- 1; e[3] <- 1
  }

  S.minus <- Sminusgenerate(p=p, kernel=kernel); S.plus <- Splusgenerate(p=p, kernel=kernel)
  S.inv    <- solve(S.minus * dgp_F1_l + S.plus * dgp_F1_r)
  C   <- dgp_F1_l * dgp_Fp1_l * Cminusgenerate(k=p+1, p=p, kernel=kernel) +
    dgp_F1_r * dgp_Fp1_r * Cplusgenerate(k=p+1, p=p, kernel=kernel)
  G       <- dgp_F1_l^3 * Gminusgenerate(p=p, kernel=kernel) + dgp_F1_r^3 * Gplusgenerate(p=p, kernel=kernel) +
    dgp_F1_l^2 * dgp_F1_r * tcrossprod(S.minus[2, ], S.plus[1, ]) +
    dgp_F1_l^2 * dgp_F1_r * tcrossprod(S.plus[1, ], S.minus[2, ])

  h.opt <- N^(-1/(2*p+1)) * (t(e) %*% S.inv %*% G %*% S.inv %*% e)^(1/(2*p+1)) *
    (abs(t(e) %*% S.inv %*% C))^(-2/(2*p+1)) * factorial(p+1)^(2/(2*p+1)) * (2*p)^(-1/(2*p+1))

  return(as.numeric(h.opt))
}

################################################################################
#' Internal function, normal distribution related quantities.
#'
#' \code{rddensity_H} calculates normal densities and derivatives.
#'
#' This is an internal function, and should not be called by users.
#'
#' @param x Numeric, point of evaluation.
#' @param p Integer, polynomial order
#'
#' @return Returns density or derivatives.
#'
#' @keywords internal
rddensity_H <- function(x, p){
  if (p==0)  out = 1
  if (p==1)  out = x
  if (p==2)  out = x^2 - 1
  if (p==3)  out = x^3 - 3*x
  if (p==4)  out = x^4 - 6*x^2 + 3
  if (p==5)  out = x^5 - 10*x^3 + 15*x
  if (p==6)  out = x^6 - 15*x^4 + 45*x^2 - 15
  if (p==7)  out = x^7 - 21*x^5 + 105*x^3 - 105*x
  if (p==8)  out = x^8 - 28*x^6 + 210*x^4 - 420*x^2 + 105
  if (p==9)  out = x^9 - 36*x^7 + 378*x^5 - 1260*x^3 + 945*x
  if (p==10) out = x^10 - 45*x^8 + 630*x^6 - 3150*x^4 + 4725*x^2 - 945
  return(out)
}
