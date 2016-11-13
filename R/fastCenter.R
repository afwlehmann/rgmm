#
# fastRowPlusMinus.R
# copyright (c) by Alexander Lehmann <afwlehmann@googlemail.com>
#


#' Subtract a given vector from every row of a given matrix.
#' @param X a matrix
#' @param mu a vector to be subtracted from every row of \code{X}
#' @return a centered matrix
#' @export
fastRowMinus <- function(X, mu) {
  N <- nrow(X)
  M <- ncol(X)

  if (length(mu) != M)
    stop("The number of columns of `X` must match the length of `mu`.")
  
  matrix(.C("fastRowMinus",
            d=as.double(X),
            as.integer(N),
            as.integer(M),
            as.double(mu),
            PACKAGE="rgmm")$d,
         nrow=N)
}


#' Add a given vector to every row of a given matrix.
#' @param X a matrix
#' @param mu a vector to be subtracted from every row of \code{X}
#' @return a centered matrix
#' @export
fastRowPlus <- function(X, mu) {
  N <- nrow(X)
  M <- ncol(X)

  if (length(mu) != M)
    stop("The number of columns of `X` must match the length of `mu`.")
  
  matrix(.C("fastRowPlus",
            d=as.double(X),
            as.integer(N),
            as.integer(M),
            as.double(mu),
            PACKAGE="rgmm")$d,
         nrow=N)
}
