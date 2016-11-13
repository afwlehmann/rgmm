# -*- coding: utf-8 -*-
#
# pdf.R
# copyright (c) by Alexander Lehmann <afwlehmann@googlemail.com>
#


#' Probability density of a given ((Semi-)Wrapped) Gaussian Mixture Model.
#' @param X a vector or a matrix where each row holds one sample and the columns
#'          denote the variables
#' @param model the model
#' @return a vector holding the probabilities
#' @export
dgmm <- function(X, model) {
  if (inherits(model, "gmm")) {
    if (inherits(model, "swgmm")) {
      # Semi-Wrapped Gaussian Mixture Model
      tmp <- exp(mapply(
          function(ll,m,S) ll + logdswmvnorm(X, mu=m, Sigma=S, wGrid=model$wGrid),
          log(model$lambda), model$mu, model$Sigma
      ))
      # `tmp` is either a matrix or a vector.
      if (is.matrix(tmp)) rowSums(tmp) else sum(tmp)
    } else {
      # Gaussian Mixture Model
      tmp <- exp(mapply(
        function(ll,m,S) ll + logdmvnorm(X, mu=m, Sigma=S),
        log(model$lambda), model$mu, model$Sigma
      ))
      # `tmp` is either a matrix or a vector.
      if (is.matrix(tmp)) rowSums(tmp) else sum(tmp)
    }
  } else
    stop("Invalid model with class ", class(model))
}


#' Probability density of a (semi-)wrapped multivariate Gaussian.
#'
#' The probability density for a (semi-)wrapped multivariate Gaussian for
#' periodic data.
#'
#' The probability density for a (semi-)wrapped multivariate Gaussian for
#' \eqn{2\pi}-periodic data is defined as
#'
#' \deqn{ p(x) = \frac{1}{\sqrt{(2\pi)^D \cdot |\Sigma|}} \cdot
#' \sum_{w \in \mathcal{Z}^D} e^{(x + 2\pi w - \mu)^T
#' \Sigma^{-1} (x + 2\pi w - \mu)} }
#'
#' where \eqn{D}, \eqn{\mu} and \eqn{\Sigma} denote the data's dimension, mean
#' and covariance matrix, and \eqn{\mathcal{Z} = -\infty \ldots
#' \infty} represents the set of integers.
#'
#' \eqn{w} represents the displacement (aka the wrapping) along the periodic
#' dimensions as vertices of the \eqn{D}-dimensional hypercube. For example, if
#' the first 2 out of 5 dimensions of a semi-wrapped multivariate Gaussian were
#' \eqn{2\pi}-periodic, then \eqn{w \in \{(2\pi k_1, 2\pi k_2, 0, 0, 0)^T | -1
#' \le k_1, k_2 \le 1\}}, assuming one wrap in every corresponding direction.
#'
#' Finally, note that as the computational costs are substantial, depending on
#' the variance of the data, 1 or 2 wraps are usually considered
#' sufficient. Also, the number of wraps need not be the same per direction.
#' @param X a vector or a matrix where each row holds one sample and the columns
#' denote the variables. Note that, when given as a vector, the shape of
#' \code{mu} and \code{Sigma} determine the shape of \code{X}. Also note that,
#' \emph{as opposed to \code{mixtools:dmvnorm}}, if neither neither \code{mu}
#' nor \code{Sigma} are given, \code{X} is assumed to be univariate and the
#' result is equal to that of \code{dnorm}.
#' @param mu the distribution's mean vector
#' @param Sigma the distribution's covariance matrix
#' @param wGrid a matrix whose rows represent the displacements along the
#' periodic dimensions as vertices of the \eqn{D}-dimensional hypercube where
#' \eqn{D} denotes the dimension of the data (see Details for further
#' explanations). If \code{is.null(wGrid)}, a default \eqn{2\pi}-periodic grid
#' will be chosen for all dimensions
#' @return a vector holding the probabilities
#' @export
dswmvnorm <- function(X, mu=NULL, Sigma=NULL, wGrid=NULL) {
  exp(logdswmvnorm(X, mu, Sigma, wGrid))
}


#' Log probability density of a (semi-wrapped) multivariate Gaussian
#' @seealso \code{\link{dswmvnorm}}
#' @param X a vector or a matrix where each row holds one sample and the columns
#' denote the variables. Note that, when given as a vector, the shape of
#' \code{mu} and \code{Sigma} determine the shape of \code{X}. Also note that,
#' \emph{as opposed to \code{mixtools:dmvnorm}}, if neither neither \code{mu}
#' nor \code{Sigma} are given, \code{X} is assumed to be univariate and the
#' result is equal to that of \code{dnorm}.
#' @param mu the distribution's mean vector
#' @param Sigma the distribution's covariance matrix
#' @param wGrid a matrix whose rows represent the displacements along the
#' periodic dimensions as vertices of the \eqn{D}-dimensional hypercube where
#' \eqn{D} denotes the dimension of the data (see Details for further
#' explanations). If \code{is.null(wGrid)}, a default \eqn{2\pi}-periodic grid
#' will be chosen for all dimensions
#' @return a vector holding the logarithm of the probabilities
#' @export
logdswmvnorm <- function(X, mu=NULL, Sigma=NULL, wGrid=NULL) {
  # Make sure that Sigma is a matrix.
  if (is.vector(Sigma))
    Sigma <- diag(Sigma, ncol=length(Sigma)) # see ?diag

  # Determine the dimension of the samples.
  D <- if (is.matrix(X)) ncol(X)
       else if (!is.null(mu)) length(mu)
       else if (!is.null(Sigma)) ncol(Sigma)
       else 1
  
  # Make sure that X is a matrix.
  if (is.vector(X))
    X <- matrix(X, ncol=D)

  # If no mean was given, assume zero mean.
  if (is.null(mu))
      mu <- rep(0, D)
  
  # If no covariance matrix was given, assume identity.
  if (is.null(Sigma))
    Sigma <- diag(D)

  # Compute a displacement grid if none was given, defaulting to +/- one
  # 2pi-periodic wraps along all axes.
  if (is.null(wGrid))
    wGrid <- (2*pi) * as.matrix(do.call(expand.grid, lapply(1:D, function(i) seq(-1, 1))))

  # For numerical stability (we're dealing with a lot of very small
  # probabilities here) the following calculations ought to be done in
  # log-space. Hence for
  #
  #    p(x) = 1 / √sqrt((2Π)^D * |Σ|) * exp(-1/2 * x^T Σ^-1 x)
  #
  # where without loss of generality x is centered around the distribution's
  # mean:
  #
  #    ln p(x) = -0.5 (D * ln(2Π) + ln |Σ| + x^T Σ^-1 x)
  #
  # When accumulating the Mahalanobis distance, however, ATTENTION must be paid
  # to a) the correct summation in log-space and b) the correct scaling:
  #
  #    ln p(x) = -0.5 (D * ln(2Π) + ln |Σ|) + \lnsum_w -0.5 * (x-w)^T Σ^-1 (x-w)
  aux <- qr(Sigma)
  logDet <- sum(log(abs(diag(aux$qr))))
  invSigma <- qr.solve(aux)

  .X <- cbind(X, +1)
  .I <- diag(ncol(X))
  rm(X)
  md <- Reduce(logSum, lapply(1:nrow(wGrid), function(w) {
    Xw <- .X %*% rbind(.I, wGrid[w,] - mu)
    -0.5 * rowSums((Xw %*% invSigma) * Xw)
  }))
  
  -0.5 * (D * 1.8378770664093453 + logDet) + md
}


#' Probability density of a multivariate Gaussian.
#'
#' This is a replacement for \code{mixtools::dmvnorm} that tries to be more
#' reasonable when \code{X} is given as a vector.
#' @examples
#' # Univariate Gaussian with mean 3 and stddev 2
#' dmvnorm(1, mu=3, Sigma=2)
#'
#' # Univariate Gaussian with mean 3 and stddev 2 (NOTE that the shape of
#' # the data, when given as a vector, depends on the shape of mu and/or Sigma)
#' dmvnorm(1:50, mu=3, Sigma=2)
#'
#' # Multivariate Gaussian with mean c(3,4) and covariance matrix ((3,-2),(-2,4))
#' dmvnorm(matrix(1:10, ncol=2), mu=c(3,4), Sigma=matrix(c(3,-2,-2,4), ncol=2))
#'
#' # In the multivariate case, one or both of mu and Sigma determine the shape
#' # of the samples when given as a vector.
#' dmvnorm(1:10, mu=c(3,4), Sigma=matrix(c(3,-2,-2,4), ncol=2))
#' @param X a vector or a matrix where each row holds one sample and the columns
#' denote the variables. Note that, when given as a vector, the shape of
#' \code{mu} and \code{Sigma} determine the shape of \code{X}. Also note that,
#' \emph{as opposed to \code{mixtools:dmvnorm}}, if neither neither \code{mu}
#' nor \code{Sigma} are given, \code{X} is assumed to be univariate and the
#' result is equal to that of \code{dnorm}.
#' @param mu the distribution's mean vector
#' @param Sigma the distribution's covariance matrix
#' @return a vector holding the probabilities
#' @export
dmvnorm <- function(X, mu=NULL, Sigma=NULL) {
  exp(logdmvnorm(X, mu, Sigma))
}


#' Log probability density of a multivariate Gaussian.
#' @seealso \code{\link{dmvnorm}}
#' @param X a vector or a matrix where each row holds one sample and the columns
#' denote the variables. Note that, when given as a vector, the shape of
#' \code{mu} and \code{Sigma} determine the shape of \code{X}. Also note that,
#' \emph{as opposed to \code{mixtools:dmvnorm}}, if neither neither \code{mu}
#' nor \code{Sigma} are given, \code{X} is assumed to be univariate and the
#' result is equal to that of \code{dnorm}.
#' @param mu the distribution's mean vector
#' @param Sigma the distribution's covariance matrix
#' @return a vector holding the logarithm of the probabilities
#' @export
logdmvnorm <- function(X, mu=NULL, Sigma=NULL) {
  # Make sure that Sigma is a matrix.
  if (is.vector(Sigma))
    Sigma <- diag(Sigma, ncol=length(Sigma)) # see ?diag

  # Determine the dimension of the samples.
  D <- if (is.matrix(X)) ncol(X)
       else if (!is.null(mu)) length(mu)
       else if (!is.null(Sigma)) ncol(Sigma)
       else 1
  
  # Make sure that X is a matrix.
  if (is.vector(X))
    X <- matrix(X, ncol=D)

  # Center the given data only if necessary.
  if (!is.null(mu))
    X <- cbind(X, -1) %*% rbind(diag(D), mu) #fastRowMinus(X, mu)

  # If no covariance matrix was given, assume identity.
  if (is.null(Sigma))
    Sigma <- diag(D)

  # For numerical stability (we're dealing with a lot of very small
  # probabilities here), the following calculations ought to be done in
  # log-space. Hence for p(x) = 1/√sqrt((2Π)^D * |Σ|) * exp(-1/2 * x^T Σ^-1 x),
  # where without loss of generality x is centered around the distribution's
  # mean:
  #
  #   ln p(x) = -0.5 (D * ln(2Π) + ln |Σ| + x^T Σ^-1 x)
  aux <- qr(Sigma)
  logDet <- sum(log(abs(diag(aux$qr))))
  invSigma <- qr.solve(aux)
  md <- rowSums((X %*% invSigma) * X)
  -0.5 * (D * 1.8378770664093453 + logDet + md)
}
