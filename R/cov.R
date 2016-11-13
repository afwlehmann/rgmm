# -*- coding: utf-8 -*-
#
# cov.R
# copyright (c) by Alexander Lehmann <afwlehmann@googlemail.com>
#


#' Map finite to infinite respective 2pi-periodic coordinates and apply function.
#'
#' Performs a non-linear change of variables for each of the given
#' /non-periodic/ coordinates from \eqn{[-1,+1]} to \eqn{[-Inf,
#' +Inf]}{[-\infty,+\infty]}, while each of the /periodic/ coordinates is
#' linearly mapped from \eqn{[-1,+1]} to \eqn{[-pi, +pi]}{[-\pi,+\pi]}.
#' @seealso \code{\link{mapInf}}
#' @param x a vector representing the point where to evaluate \code{fun}
#' @param periodics the indices of the 2pi-periodic variables in \code{x}
#' @param fun the function to be applied
#' @param ... passed on to \code{fun}
#' @export
mapInfPeriodic <- function(x, periodics, fun, ...) {
  stopifnot(is.vector(x))

  # Non-linear mapping from [-1,+1] to [-Inf,+Inf] except for those variables
  # given by the indices in `periodics`, which will be linearly mapped to
  # [-pi,+pi].
  y <- x^2
  xPrime <- x / (1 - y)
  auxDetJ <- (1 + y) / (1 - y)^2

  if (sum(periodics) > 0) {
    xPrime[periodics] <- x[periodics] * pi
    auxDetJ[periodics] <- pi
  }

  # See `mapInf` for why the Jacobian's determinant is equal to the product of
  # its diagonal in this case.
  detJ <- prod(replace(auxDetJ, is.infinite(auxDetJ), 1))
  fun( xPrime, ... ) * detJ
}


#' Map finite to infinite coordinates and apply function.
#'
#' Performs a non-linear change of variables for each of the given coordinates
#' from \eqn{[-1,+1]} to \eqn{[-Inf, +Inf]}{[-\infty,+\infty]}.
#' 
#' Rule for integration by substitution (univariate):
#' \deqn{
#'   \int_{\phi(a)}^{\phi(b)} f(x) dx = \int_a^b f(\phi(t)) \phi'(t) dt
#' }
#' 
#' And in the multivariate case:
#' \deqn{
#'   \int_{\phi(a)}^{\phi(b)} f(x) dx = \int_a^b f(\phi(t)) | J_{\phi}(t) | dt
#' }
#' 
#' where \eqn{J_{\phi}} denotes the Jacobian of \eqn{\phi} and consequently
#' \eqn{| J_{\phi}(t) |} the determinant of the Jacobian of \eqn{\phi} evaluated
#' at \eqn{t}.
#' @param x a vector representing the point where to evaluate \code{fun}
#' @param fun the function to be applied
#' @param ... passed on to \code{fun}
#' @export
mapInf <- function(x, fun, ...) {
  stopifnot(is.vector(x))

  y <- x^2
  xPrime <- x / (1 - y)
  auxDetJ <- (1 + y) / (1 - y)^2
  
  # If \eqn{\phi} substitutes all variables independently from each other, the
  # Jacobian matrix is in fact a diagonal matrix. Hence its determinant is equal
  # to the product of its diagonal elements:
  #
  #  \phi( (x, y, z) )     = ( x/(1-x^2), y/(1-y^2), z/(1-z^2) )
  #
  #                          [ (1+x^2)/(1-x^2)^2  0                   0                 ]
  #  J_{\phi}( (x, y, z) ) = [ 0                  (1+y^2)/(1-y^2)^2   0                 ]
  #                          [ 0                  0                   (1+z^2)/(1-z^2)^2 ]
  detJ <- prod(replace(auxDetJ, is.infinite(auxDetJ), 1))
  
  fun( xPrime, ... ) * detJ
}


