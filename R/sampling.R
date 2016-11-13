# -*- coding: utf-8 -*-
#
# sampling.R
# copyright (c) by Alexander Lehmann <afwlehmann@googlemail.com>
#


#' Draw samples from a given ((Semi-)Wrapped) Gaussian Mixture Model.
#' @param N the number of samples to be drawn
#' @param model an instance of \code{gmm}
#' @return a matrix whose rows contain the samples drawn from the given
#'  \code{model}
rgmm <- function(N, model) {
  stopifnot(inherits(model, "gmm"))

  breaks <- Reduce(`+`, model$lambda, accumulate=T)
  r <- cut(runif(N, 0, 1), breaks=c(0, breaks), labels=F)
  do.call(rbind,
          mapply(function(i, m, S) {
            K <- sum(r == i)
            if (K > 0) rmvnorm(K, m, S) else NULL
          }, seq(length(model$mu)), model$mu, model$Sigma))
}


#' Draw random samples from a multivariate Gaussian.
#' @param N the number of samples
#' @param mu the mean. If \code{mean} is \code{NULL}, zero mean is assumed for
#' every variable
#' @param Sigma the covariance matrix or a vector of variances. If \code{Sigma}
#' is \code{NULL}, unit variance is assumed for every variable
#' @return a matrix with one sample per row
#' @export
rmvnorm <- function(N, mu=NULL, Sigma=NULL) {
  if (is.vector(Sigma))
    Sigma <- diag(Sigma, nrow=length(Sigma))

  D <-
    if (is.null(mu)) {
      if (is.null(Sigma))
        1
      else
        ncol(Sigma)
    } else
      length(mu)

  A <- if (is.null(Sigma)) diag(D) else chol(Sigma)
  mu <- if (is.null(mu)) rep(0, D) else mu

  samples <- matrix(rnorm(D * N), ncol=D)
  fastRowPlus(samples %*% A, mu)
}
