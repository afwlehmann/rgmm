# -*- coding: utf-8 -*-
#
# plot.R
# copyright (c) by Alexander Lehmann <afwlehmann@googlemail.com>
#


#' Multidimensional contour plot of a given Gaussian Mixture Model.
#'
#' All possible combinations of the model's axes are layed out in a grid where
#' each cell is then filled with a two-dimensional contour plot of the given
#' GMM. A corresponding two-dimensional histogram will be plotted if a dataset
#' is provided.
#' @param x the model (an instance of \code{gmm})
#' @param dataset a matrix or data-frame whose rows each contain one
#'  d-dimensional observation where d corresponds to the dimensionality of the
#'  given \code{x}. When given as a vector, the vector will be converted to
#'  a one-column matrix. May be \code{NULL} as long as \code{limits} is given.
#' @param limits a list of 2-tuples, each describing the extents of the
#'  corresponding axis. May be \code{NULL} as long as \code{dataset} is given. In
#'  the latter case, the axes' extents are derived from the columns of
#'  \code{dataset}.
#' @param resolution a vector describing the resolution per axis
#' @param col the color to use for the contour plot
#' @param oma the outer margin, see \code{par}
#' @param ... passed on to the \code{contour} function
#' @importFrom gplots hist2d
#' @method plot gmm
#' @export plot.gmm
plot.gmm <- function(x, dataset=NULL, limits=NULL,
                     resolution=100, col="magenta", oma=c(2,1,3.5,1), ...)
{
  stopifnot(inherits(x, "gmm"))

  D <- length(x$mu[[1]])

  # Make sure that `dataset` is a (suitable) matrix.
  if (is.vector(dataset))
    dataset <- matrix(dataset, ncol=D)

  # Make sure we have proper `limits`.
  if (is.null(limits)) {
    if (is.null(dataset))
      stop("Either `limits`, `dataset`, or both must be given.")
    limits <- lapply(1:D, function(i) range(dataset[,i]))
  } else if (length(limits) != D)
    stop("The length of `limits` does not correspond to the model.")

  # First, make sure that `resolution` is a vector of length `D` as we're going
  # to use it as the `dim` parameter in the next call to `array`.
  resolution <- as.integer(resolution)
  if (length(resolution) != D)
    resolution <- rep(resolution, length.out=D)

  # Plot all vs. all.
  parBackup <- par(no.readonly=TRUE)
  par(oma=oma, mar=rep(0.2, 4), mgp=c(3, 0.75, 0), las=1)
  layout(matrix(seq(1, D^2), ncol=D), respect=TRUE)
  for (i in 1:D) {
    for (j in 1:D) {
      xmm <- range(limits[[i]])
      ymm <- range(limits[[j]])

      plot.new()
      plot.window(xlim=xmm, ylim=ymm)
      box()
      
      if (i == j) {
        text(x=sum(xmm)/2, y=sum(ymm)/2, labels=colnames(dataset)[i])
      } else {
        if (!is.null(dataset)) {
          par(new=TRUE)
          suppressWarnings(
            hist2d(dataset[,c(i,j)], xaxt="n", yaxt="n", xaxs="i", yaxs="i",
                   xlim=xmm, ylim=ymm, col=c("black", topo.colors(31)))
            )
        }
        par(new=TRUE)
        m <- marginal(x, c(i,j))
        dx <- seq(xmm[1], xmm[2], length.out=resolution[i])
        dy <- seq(ymm[1], ymm[2], length.out=resolution[j])
        dz <- local({
          grid <- as.matrix(expand.grid(dx, dy))
          matrix(dgmm(grid, model=marginal(x, c(i,j))), ncol=resolution[i])
        })
        contour(dx, dy, dz,
                xaxt="n", yaxt="n", xaxs="i", yaxs="i", col=col,
                levels=sapply(1:10, function(i) 5*10^(-i)),
                xlim=xmm, ylim=ymm, add=TRUE, ...)
      }

      # Top/bottom axis labels
      if (i %% 2 == 1 && j == D)
        axis(1, at=pretty(limits[[i]], n=5), ...)  # bottom
      else if (i %% 2 == 0 && j == 1)
        axis(3, at=pretty(limits[[i]], n=5), ...)  # top
      # Left/right axis labels
      if (j %% 2 == 1 && i == D)
        axis(4, at=pretty(limits[[j]], n=5), ...)  # right
      else if (j %% 2 == 0 && i == 1)
        axis(2, at=pretty(limits[[j]], n=5), ...)  # left
    }
  }
  par(parBackup)
}
