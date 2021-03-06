% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gmm.R
\name{estimateGMMParameters}
\alias{estimateGMMParameters}
\title{Estimate GMM parameters for a given dataset.}
\usage{
estimateGMMParameters(X, K = NULL, lambda = NULL, mu = NULL,
  Sigma = NULL, ...)
}
\arguments{
\item{X}{a N x M matrix with N samples and M variables}

\item{K}{the number of Gaussians}

\item{lambda}{a vector of mixing coefficients}

\item{mu}{a list of means}

\item{Sigma}{a list of covariance matrices}

\item{...}{parameters passed on to \code{kmeans}}
}
\value{
a list consisting of
 \item{K}{the number of Gaussians}
 \item{lambda}{a vector of mixing coefficients}
 \item{mu}{a list of means}
 \item{Sigma}{a list of covariance matrices}
}
\description{
Estimation of any missing parameters is done via clustering of the given
dataset and subsequent analysis of the found clusters.
}

