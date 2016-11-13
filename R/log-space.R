# -*- coding: utf-8 -*-
#
# log-space.R
# copyright (c) by Alexander Lehmann <afwlehmann@googlemail.com>
#


#' The logarithm of the sum of `a` and `b` where both are given in log-space.
#'
#' Proof:
#' \deqn{
#' \begin{array}{rl}
#' \ln(a + b) &= \ln a + \ln(a+b) - \ln a \\
#' &= \ln a + \ln\frac{a+b}{a} \\
#' &= \ln a + \ln (1 + \frac{b}{a}) \\
#' &= \ln a + \ln (1 + e^{\ln b - \ln a})
#' \end{array}
#' }
#' @param a a vector of values in log-space
#' @param b a vector of values in log-space
#' @return the sum of `a` and `b` in log-space
#' @examples
#' all.equal( log(7), logSum(log(3), log(4)) )
#' all.equal( log(7), logSum(log(4), log(3)) )
#' all.equal( log(7), logSum(log(0), log(7)) )
#' all.equal( log(7), logSum(log(7), log(0)) )
#' all.equal( log(0), logSum(log(0), log(0)) )
#' all.equal( Inf   , logSum(Inf, Inf)       )
#' all.equal( Inf   , logSum(Inf, -Inf)      )
#' all.equal( Inf   , logSum(-Inf, Inf)      )
#' all.equal( Inf   , logSum(Inf, log(3))    )
#' all.equal( Inf   , logSum(log(3), Inf)    )
#' @export
logSum <- function(a, b) {
  upper <- pmax(a, b)
  lower <- pmin(a, b)
  infInd <- is.infinite(lower)
  r <- upper + log1p(exp(lower - upper))
  replace(r, infInd, upper[infInd])
}
