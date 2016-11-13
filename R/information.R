# -*- coding: utf-8 -*-
#
# information.R
# copyright (c) by Alexander Lehmann <afwlehmann@googlemail.com>
#


#' Bayesian Information Criterion
#'
#' Calculates the BIC with respect to the given \code{model}. In general, the
#' BIC is defined as \eqn{-2 \cdot \ln L + k \cdot \ln N} where \eqn{k} denotes
#' the number of free model parameters, \eqn{k} denotes the number of free model
#' parameters, and \eqn{L} denotes the maximized value of the likelihood for the
#' given model.
#' @param model the model, i.e. an instance of \code{gmm}
#' @param N the number of samples of the dataset used for training
#' @return the BIC
#' @export
calcBIC <- function(model, N) {
  numGaussians <- length(model$lambda)
  dim <- length(model$mu[[1]])
  llh <- tail(model$llh, 1)
  #    lambda  mu   Sigma (symmetric incl. diagonal)
  k <- ( 1 +  dim + (dim+1)*dim/2 ) * numGaussians
  -2 * llh + k * log(N)
}


#' Akaike Information Criterion
#'
#' Calculates the AIC with respect to the given \code{model}. In general, the
#' AIC is defined as \eqn{2 \cdot k - 2 \ln L} where \eqn{k} denotes the
#' number of free model parameters and \eqn{L} denotes the maximized value of
#' the likelihood for the given model.
#' @param model the model, i.e. an instance of \code{gmm}
#' @return the AIC
#' @export
calcAIC <- function(model) {
  numGaussians <- length(model$lambda)
  dim <- length(model$mu[[1]])
  llh <- tail(model$llh, 1)
  #    lambda  mu   Sigma (symmetric incl. diagonal)
  k <- ( 1 +  dim + (dim+1)*dim/2 ) * numGaussians
  2 * k - 2 * llh
}


#' Akaike Information Criterion (corrected)
#'
#' Calculates the AICc with respect to the given \code{model}. In general, the
#' AICc is defined as \deqn{2 \cdot k - 2 \ln L + \frac{2 \cdot k(k+1)}{N-k-1} =
#' \textit{AIC} + \frac{2 \cdot k(k+1)}{N-k-1}} where \eqn{k} denotes the number
#' of free model parameters and \eqn{L} denotes the maximized value of the
#' likelihood for the given model.
#' @param model the model, i.e. an instance of \code{gmm}
#' @param N the number of samples of the dataset used for training
#' @return the AICc
#' @export
calcAICc <- function(model, N) {
  numGaussians <- length(model$lambda)
  dim <- length(model$mu[[1]])
  llh <- tail(model$llh, 1)
  #    lambda  mu   Sigma
  k <- ( 1 +  dim + (dim+1)*dim/2 ) * numGaussians
  2 * k - 2 * llh + 2 * k * (k + 1) / (N - k - 1)
}


#' Differential entropy of a multivariate normal distribution.
#' @param Sigma the covariance matrix (or a vector of the matrix diagonal)
#' @return the differential entropy
#' @export
demvnorm <- function(Sigma) {
  # Make sure that Sigma is a matrix.
  if (is.vector(Sigma))
    Sigma <- diag(Sigma, ncol=length(Sigma))

  # Note that the result does not depend on the mean of the distribution, since
  #
  #   H(X) = -\int p(x) \log p(x) dx
  #
  # where w.l.o.g. for a zero-mean distribution
  #         p(x) = 1/sqrt((2pi)^D |Σ|) * exp(-1/2 x^T Σ^-1 x)  \\
  # => \log p(x) = -1/2 ( D \log(2pi) + \log |Σ| + x^T Σ^-1 x)
  #
  # from which it follows that
  #
  #   H(X) = -\int p(x) \log p(x) dx
  #        = -\int p(x) (-1/2) ( D \log(2pi) + \log |Σ| + x^T Σ^-1 x) dx
  #        = 1/2 (D \log(2pi) + \log |Σ|) \int p(x) dx + 1/2 \int p(x) x^T Σ^-1 x dx
  #        = 1/2 (D \log(2pi) + \log |Σ| + D)
  #
  D <- ncol(Sigma)
  aux <- qr(Sigma)
  invSigma <- qr.solve(aux)
  logDet <- sum(log(abs(diag(aux$qr))))
  0.5 * (D * log(2*pi) + logDet + D)
}


#' Differential entropy of the given (Semi-Wrapped) Gaussian Mixture Model.
#' @param model an instance of \code{gmm} or \code{swgmm}
#' @return the differential entropy
#' @importFrom cubature adaptIntegrate
#' @export
degmm <- function(model) {
  stopifnot(inherits(model, "gmm") || inherits(model, "swgmm"))

  # Determine the dimension of the data.
  D <- length(model$mu[[1]])

  # Precompute what can be precomputed.
  tmp <- lapply(model$Sigma, function(S) qr(S))
  {
    model$invSigma  <- lapply(tmp, function(q) qr.solve(q))
    model$logDet    <- lapply(tmp, function(q) sum(log(abs(diag(q$qr)))))
    model$logLambda <- log(model$lambda)
  }
  rm(tmp); gc()

  # Auxiliary function: Logarithm of the probability density of a GMM.
  logPDFGMM <- function(x) {
    # Compute the partial probabilities.
    tmp <- mapply(function(ll,m,invS,ld) {
      xPrime <- x - m
      md <- sum((xPrime %*% invS) * xPrime)
      ll - 0.5 * (D * 1.8378770664093453 + ld + md)
    }, model$logLambda, model$mu, model$invSigma, model$logDet)
    # Accumulate the result.
    Reduce(logSum, tmp)
  }

  # Auxiliary function: Logarithm of the probability density of a SW-GMM.
  logPDFSWGMM <- function(x) {
    # Compute the partial probabilities
    tmp <- mapply(function(ll,m,invS,ld) {
      XPrime <- fastRowPlus(model$wGrid, x - m)
      md <- Reduce(logSum, -0.5 * rowSums((XPrime %*% invS) * XPrime))
      ll - 0.5 * (D * 1.8378770664093453 + ld) + md
    }, model$logLambda, model$mu, model$invSigma, model$logDet)
    # Accumulate the result.
    Reduce(logSum, tmp)
  }

  if (inherits(model, "swgmm")) {
    integrand <- function(x) {
      stopifnot(is.vector(x))
      lp <- logPDFSWGMM(x)
      if (is.infinite(lp)) 0 else exp(lp) * lp
    }
    periodics <- which(model$wraps > 0)
    -adaptIntegrate(mapInfPeriodic, rep(-1, D), rep(+1, D), periodics=periodics, fun=integrand)$integral
  } else if (inherits(model, "gmm")) {
    integrand <- function(x) {
      stopifnot(is.vector(x))
      lp <- logPDFGMM(x)
      if (is.infinite(lp)) 0 else exp(lp) * lp
    }
    -adaptIntegrate(mapInf, rep(-1, D), rep(+1, D), fun=integrand)$integral
  } else stop("This is a bug!")
}


#' Discrete entropy of a (Semi-Wrapped) Gaussian Mixture Model.
#' @param X a matrix whose rows represent the multivariate samples
#' @param model an instance of \code{gmm} or \code{swgmm}
#' @return the discrete entropy (in nats)
#' @export
entropy <- function(X, model) {
  p <- dgmm(X, model=model)
  -sum(p * log(p))
}
