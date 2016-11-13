# -*- coding: utf-8 -*-
#
# gmm.R
# copyright (c) by Alexander Lehmann <afwlehmann@googlemail.com>
#


#' EM estimation of a Gaussian Mixture Model.
#'
#' EM estimation of a Gaussian Mixture Model that fits the given dataset
#' \code{X}.
#' @param X a N x M matrix with N samples and M variables
#' @param K an integer specifying the number of Gaussians
#' @param lambda a vector of real mixing coefficients
#' @param mu a list of mean vectors
#' @param Sigma a list of covariance matrices
#' @param epsilon the covergence criterion (0 < epsilon < 1), iteration will
#' stop once last and current likelihood differ by less than \code{epsilon *
#' 100} percent
#' @param maxIter the maximum number of iterations
#' @param verb be verbose
#' @return an instance of \code{gmm} whose parameters represent the model
#'  with the best log-likelihood, i.e.
#'  \item{lambda}{a vector of real mixing coefficients}
#'  \item{mu}{a list of mean vectors}
#'  \item{Sigma}{a list of covariance matrices}
#'  \item{llh}{a vector of the log-likelihood for every iteration step}
#' @export
computeGMM <- function(X, K=2, lambda=NULL, mu=NULL, Sigma=NULL,
                       epsilon=0.05, maxIter=250, verb=TRUE)
{
  stopifnot(0 < K)
  stopifnot(0 < epsilon && epsilon < 1)
  stopifnot(0 < maxIter)

  # Make sure that X is a matrix.
  X <- as.matrix(X)
  N <- nrow(X)

  model <- expectationMaximization(X, K, lambda, mu, Sigma, epsilon, maxIter, verb)
  structure(model, class=c("gmm"))
}


#' EM estimation of a (Semi-)Wrapped Gaussian Mixture Model.
#'
#' EM estimation of a (Semi-)Wrapped Gaussian Mixture Model that fits the given
#' dataset \code{X}.  Means, covariance matrices and mixing coefficients can be
#' given e.g. as an initial guess.\cr
#' Periodic variables will be considered \eqn{2\pi}-periodic. The number of wraps
#' can be specified per variable through the use of \code{wraps}.\cr
#' \code{trafo} can be given to apply a linear transformation to the computed
#' displacement grid, e.g. for whitened datasets whose axes are not aligned with
#' the original axes anymore.
#' @seealso For further details, in particular concerning the displacement grid,
#' see \code{\link{dswmvnorm}}.
#' @param X a N x M matrix with N samples and M variables
#' @param wraps an integer vector of length \code{ncol(X)}, describing the
#' number of wraps per variable, e.g. \code{c(2,0,1)} for a data set with two
#' periodic variables of which the first should be wrapped twice whereas the
#' last should be wrapped only once
#' @param K an integer specifying the number of Gaussians
#' @param lambda a vector of real mixing coefficients
#' @param mu a list of mean vectors
#' @param Sigma a list of covariance matrices
#' @param trafo linear transformation of the displacement grid, defaults to
#' identity
#' @param invTrafo inverse of \code{trafo}. If \code{invTrafo} is \code{NULL}
#' it will be computed via QR decomposition
#' @param epsilon the covergence criterion (0 < epsilon < 1), iteration will
#' stop once last and current likelihood differ by less than \code{epsilon *
#' 100} percent
#' @param maxIter the maximum number of iterations
#' @param verb be verbose
#' @return an instance of \code{swgmm} (also: \code{gmm}) whose parameters
#' represent the model with the best log-likelihood, i.e.
#'  \item{lambda}{a vector of mixing coefficients}
#'  \item{mu}{a list of mean vectors}
#'  \item{Sigma}{a list of covariance matrices}
#'  \item{wraps}{a copy of the given \code{wraps}}
#'  \item{wGrid}{for internal use}
#'  \item{llh}{a vector of the log-likelihood for every iteration step}
#' @export
computeSWGMM <- function(X, wraps, K=2, lambda=NULL, mu=NULL, Sigma=NULL,
                         trafo=diag(NCOL(X)), invTrafo=NULL,
                         epsilon=0.05, maxIter=250, verb=TRUE)
{
  wraps <- as.integer(wraps)
  if (all(wraps == 0))
    warning("No periodic variables. Consider using computeGMM instead.")

  stopifnot(all(0 <= wraps) && length(wraps) == NCOL(X))
  stopifnot(0 < K)
  stopifnot(0 < epsilon && epsilon < 1)
  stopifnot(0 < maxIter)

  # If a transformation was given, compute its inverse. The inverse is needed
  # because we need to wrap the means of the Gaussians during each iteration
  # step and for the computation of the "safe" grid at the end of the process.
  invTrafo <-
    if (is.null(invTrafo)) {
      tryCatch(qr.solve(qr(trafo)),
               error=function(e) stop("`trafo` is singular!"))
    } else invTrafo
  
  # Make sure that X is a matrix.
  X <- as.matrix(X)

  # Initialize the 2pi-periodic displacement grid and apply the possibly given
  # transformation.
  .wGrid <- (2*pi) * as.matrix(
      do.call(expand.grid, lapply(wraps, function(nw) seq(-nw, nw)))
    ) %*% trafo

  # Y contains a copy of the data points in X for every possible displacement.
  Y <- do.call(rbind, lapply(split(.wGrid, row(.wGrid)), function(w) fastRowPlus(X, w)))

  # And clean up a bit.
  rm(X)
  
  model <- expectationMaximization(Y, K, lambda, mu, Sigma, epsilon, maxIter, verb)
  structure(c(model, list(wGrid=.wGrid, wraps=wraps)),
            class=c("swgmm", "gmm"))
}


expectationMaximization <- function(X, K, lambda, mu, Sigma, epsilon, maxIter, verb) {
  # Finalize the model parameters.
  ip <- estimateGMMParameters(X, as.integer(K), lambda, mu, Sigma)
  K      <- ip$K
  lambda <- ip$lambda
  mu     <- ip$mu
  Sigma  <- ip$Sigma

  # Auxiliary function that computes responsibilities and log-likelihood.
  computeResponsibilities <- function() {
    gamma <- exp(mapply(function(ll,m,S) ll + logdmvnorm(X, m, S),
                        log(lambda), mu, Sigma))
    # The sum of the n-th row of `gamma` represents the total probability
    # of the n-th row of `X` for this model.
    rs <- rowSums(gamma)
    # NaNs might occur along with very small `rs`, but they are dealt with
    # during the singularity check, so there's no need to handle them here.
    list(llh=sum(log(rs)), gamma=gamma/rs)
  }

  # Compute the responsibilities once during initialization.
  resp <- computeResponsibilities()
  llh <- resp$llh
  gamma <- resp$gamma
  if (verb)
    cat(sprintf("%.18e\n", resp$llh))

  .X <- cbind(X, -1)      # auxiliary copy for faster column-wise subtraction
  .I <- diag(ncol(.X)-1)  # likewise
  bestModel <- NULL
  diff <- 1
  for (iter in 1:maxIter) {
    # Check for singularities.
    if (!is.finite(llh) || any(is.nan(gamma)) || any(sapply(Sigma, det) < 1e-100)) {
      # Oops.
      warning("Parameter reestimation due to singularity!")
      ip <- estimateGMMParameters(X, K=K, nstart=2)
      lambda <- ip$lambda
      mu     <- ip$mu
      Sigma  <- ip$Sigma
    } else {
      # Everything seems fine.
      Nk     <- colSums(gamma)
      lambda <- Nk / nrow(X)
      mu     <- lapply(1:K, function(i) colSums(gamma[,i] * X) / Nk[i])
      Sigma  <- lapply(1:K, function(i) {
        Xprime <- .X %*% rbind(.I, mu[[i]]) #fastRowMinus(X, mu[[i]])
        crossprod(sqrt(gamma[,i]) * Xprime) / Nk[i]
      })
    }

    # Recompute the responsibilities.
    resp <- computeResponsibilities()
    gamma <- resp$gamma
    if (verb)
      cat(sprintf("%.18e (%d)\n", resp$llh, iter))

    # Explicit call to the garbage collector.
    gc()

    # Update `bestModel` if applicable.
    # Note that the EM algorithm is guaranteed to converge, i.e. llh_{t+1) >=
    # llh_{t}, where t denotes the iteration. This mean that there's actually no
    # need to save the `bestModel`. However, due to its liability to
    # singularities the algorithm may start over with a new set of parameters,
    # which is why it makes sense to store the `bestModel` so far.
    if (is.null(bestModel) || resp$llh > max(llh))
      bestModel <- list(lambda = lambda, mu = mu, Sigma = Sigma)

    # See if we have converged.
    diff <- if (tail(llh,1) > resp$llh) tail(llh,1) - resp$llh
            else resp$llh - tail(llh,1)
    llh <- c(llh, resp$llh)
    if (is.finite(diff) && diff <= log1p(epsilon))
      break
  }

  if (iter >= maxIter && (!is.finite(diff) || diff > log1p(epsilon)))
    warning("EM did not converge in ", iter, " iterations.")

  c(bestModel, list(llh=llh))
}


#' Estimate GMM parameters for a given dataset.
#'
#' Estimation of any missing parameters is done via clustering of the given
#' dataset and subsequent analysis of the found clusters.
#' @param X a N x M matrix with N samples and M variables
#' @param K the number of Gaussians
#' @param lambda a vector of mixing coefficients
#' @param mu a list of means
#' @param Sigma a list of covariance matrices
#' @param ... parameters passed on to \code{kmeans}
#' @return a list consisting of
#'  \item{K}{the number of Gaussians}
#'  \item{lambda}{a vector of mixing coefficients}
#'  \item{mu}{a list of means}
#'  \item{Sigma}{a list of covariance matrices}
#' @export
estimateGMMParameters <- function(X, K=NULL, lambda=NULL, mu=NULL, Sigma=NULL, ...) {
  if (is.null(K) && is.null(lambda) && is.null(mu) && is.null(Sigma))
    stop("Cannot guess K without lambda, mu or Sigma!")

  tmp <- c(K, length(lambda), length(mu), length(Sigma))
  K <- max(tmp)
  if ((!is.null(lambda) && tmp[2] != K) ||
      (!is.null(mu)     && tmp[3] != K) ||
      (!is.null(Sigma)  && tmp[4] != K))
  {
    stop("All of lambda, mu and Sigma must have equal length.")
  }

  # Make sure that X is a matrix.
  if (!is.matrix(X))
    X <- as.matrix(X)

  if (nrow(X) < K*2)
      stop("Not enough data! Provide at least K*2 samples.")
  
  dim <- ncol(X)

  # Check that all means and covariance matrices have matching dimensions.
  if (!all(mapply(function(m, S)
                    is.vector(m) && length(m) == dim &&
                    is.matrix(S) && length(S) == dim^2,
                  if (!is.null(mu)) mu else list(seq(dim)),          # the latter
                  if (!is.null(Sigma)) Sigma else list(diag(dim))))) # are dummies
    stop("Both mu and Sigma have to match the shape of X.")

  # Use K-Means in order to determine any missing parameters.
  if (is.null(lambda) || is.null(mu) || is.null(Sigma)) {
    km <- kmeans(X,
                 if (is.null(mu)) K else matrix(unlist(mu), ncol=dim, byrow=T),
                 ...)

    if (is.null(mu))
      mu <- lapply(split(km$centers, seq(K)), as.vector)
    
    if (is.null(Sigma))
      Sigma <- lapply(split(as.data.frame(X), km$cluster), cov)
    
    if (is.null(lambda))
      lambda <- km$size / sum(km$size)
  }

  # If any of the covariance matrices is singular, perform a random split
  # of the whole dataset, recompute and double-check.
  if (any(sapply(Sigma, Negate(is.finite))) ||
      any(sapply(Sigma, det) < 1e-100)) {
    # Random split, see https://gist.github.com/dsparks/3695362.
    splitIndices <- local({
      .draw <- rnorm(nrow(X))
      cut(.draw, quantile(.draw, 0:K/K), include.lowest=T, labels=F)
    })
    splitX <- split(as.data.frame(X), splitIndices)
    mu <- lapply(splitX, colMeans)
    Sigma <- lapply(splitX, cov)
    lambda <- sapply(splitX, nrow) / nrow(X)
    rm(splitX)
    if (any(sapply(Sigma, Negate(is.finite))) ||
        any(sapply(Sigma, det) < 1e-100)) {
      stop("Parameter estimation failed. The dataset appears to have outliers.")
    }
  }

  list(K=K, lambda=lambda, mu=mu, Sigma=Sigma)
}


#' Inverse transformation of a given Gaussian Mixture Model.
#'
#' If the data that were used during the EM-estimation of the given model were
#' subject to some linear transformation `H`, for example due to whitening of
#' the dataset prior to the EM-estimation, then consequently the model's
#' parameters reflect that transformation as well.\cr This function exists
#' merely for convenience. It applies the given inverse transformation
#' \code{invH} to all of the model's parameters plus it adds the given
#' \code{delta} to the model's means.
#' @param model an instance of \code{gmm}
#' @param invH the matrix inverse of the original linear transformation
#' @param delta a vector to be added to each of the model's means
#' @return an "inversely transformed" model. Note that the log-likelihood
#'  history \code{llh} of the resulting model will not be removed, however it
#'  will be marked as invalid by adding a boolean attribute \code{invalid} with
#'  value \code{TRUE}.
#' @export
invTransform <- function(model, invH, delta) {
  stopifnot(inherits(model, "gmm"))

  # Means (transform and "uncenter"):
  model$mu <- lapply(model$mu, function(m) as.vector(m %*% invH) + delta)
  
  # Covariance matrices (transform only):
  # Let ∑be the covariance matrix of the dataset X, and let X' = X H for a
  # linear transformation H. Then X = X' H^{-1} and
  #
  #   ∑ = X^T X
  #     = (X' H^{-1})^T X' H^{-1}
  #     = H^{-T} X'^T X' H^{-1}
  #     = H^{-T} ∑' H^{-1}
  #
  # where ∑' is the covariance matrix of X'.
  model$Sigma <- lapply(model$Sigma, function(s) t(invH) %*% s %*% invH)

  if (inherits(model, "swgmm"))
      model$wGrid <- model$wGrid %*% invH

  attr(model$llh, "invalid") <- TRUE

  model
}


#' Marginal distribution of a (Semi-Wrapped) Gaussian Mixture Model.
#'
#' Marginal distribution of the given SW-GMM \code{model}. Only those variables
#' with corresponding indices in \code{margin} are kept.
#' @examples
#' # Create example model.
#' model <- structure(
#'   list(mu=list(c(-14, -14), c(0, 0), c(15, -15)),
#'        Sigma=list(35*diag(2), 25*diag(2), 18*diag(2)),
#'        lambda=c(0.3, 0.5, 0.2)),
#'   class="gmm"
#' )
#' integrate(function(x) dgmm(x, marginal(model, 1)), -Inf, +Inf)
#' @param model an instance of \code{gmm}
#' @param margin an integer vector of the indices of the variables to keep
#' @return a marginalized copy of the given \code{model}. Note that the
#'  log-likelihood history \code{llh} of the resulting model will not be
#'  removed, however it will be marked as invalid by adding a boolean attribute
#'  \code{invalid} with value \code{TRUE}.\cr If the given \code{model} is an
#'  instance of \code{swgmm}, an instance of \code{swgmm} will be returned if
#'  and only if the marginalized model has any \code{wraps > 0}.
#' @export
marginal <- function(model, margin) {
  stopifnot(inherits(model, "gmm"))
  if (any(margin < 1) || any(margin > length(model$mu[[1]])))
    stop("Invalid margin: ", margin)

  # Copy the given model and change only what needs to be changed.
  result       <- model

  result$mu    <- lapply(model$mu, function(m) m[margin])
  result$Sigma <- lapply(model$Sigma, function(S) S[margin,margin,drop=F])

  if (sum(result$wraps[margin]) > 0) {
    result$wraps <- result$wraps[margin]
    result$wGrid <- as.matrix(do.call(expand.grid, lapply(margin, function(m) {
      values <- result$wGrid[,margin]
      unique(replace(values, abs(values) < 1e-15, 0))
    })))
  } else {
    result$wraps <- NULL
    result$wGrid <- NULL
    class(result) <- setdiff(class(result), c("swgmm"))
  }

  if (!is.null(result$llh))
      attr(result$llh, "invalid") <- TRUE

  result
}
