#
# test-information.R
# copyright (c) by Alexander Lehmann <afwlehmann@googlemail.com>
#


suppressMessages(library(cubature))


context("information")


test_that("demvnorm is correct", {
  integralSolution <- function(Sigma) {
    f <- function(x) {
      stopifnot(is.vector(x))
      lp <- logdmvnorm(x, Sigma=Sigma)
      if (is.infinite(lp)) 0 else exp(lp) * lp
    }
    D <- NCOL(Sigma)
    -adaptIntegrate(mapInf, rep(-1, D), rep(1, D), fun=f, maxEval=50000)$integral
  }

  # unit variance
  for (i in 1:3) {
    expect_equal(demvnorm(diag(i)), integralSolution(diag(i)),
                 tolerance=1e-4, info=paste("dim ==", i))
  }

  # "random" variance
  mapply(function(S, i) {
    expect_equal(demvnorm(diag(S, ncol=i)), integralSolution(diag(S, ncol=i)),
                 tolerance=1e-4, info=paste("dim ==", i, "non-unit"))
  }, list(c(53.28), c(89.34, 55.59), c(1.20, 7.72, 4.88)), 1:3)
})


test_that("degmm is correct for GMMs", {
  integralSolution <- function(model) {
    f <- function(x) {
      stopifnot(is.vector(x))
      p  <- dgmm(x, model)
      lp <- log(p)
      if (is.infinite(lp)) 0 else p * lp
    }
    D <- length(model$mu[[1]])
    -adaptIntegrate(mapInf, rep(-1, D), rep(1, D), fun=f, maxEval=50000)$integral
  }

  # "random" variance
  mus     <- list(c(0),     c(-7,23),        c(42,-49,3))
  Sigmas  <- list(c(53.28), c(89.34, 55.59), c(1.20, 7.72, 4.88))
  lambdas <- list(c(1),     c(0.4, 0.6),     c(0.4, 0.1, 0.5))
  mapply(function(l, m, S) {
    model <- structure(list(lambda=l, mu=m, Sigma=S), class=c("gmm"))
    expect_equal(degmm(model), integralSolution(model))
  }, lambdas, mus, Sigmas)
})


test_that("degmm is correct for SW-GMMs", {
  integralSolution <- function(model) {
    f <- function(x) {
      stopifnot(is.vector(x))
      p  <- dgmm(x, model)
      lp <- log(p)
      if (is.infinite(lp)) 0 else p * lp
    }
    D <- length(model$mu[[1]])
    periodics <- which(model$wraps > 0)
    -adaptIntegrate(mapInfPeriodic, rep(-1, D), rep(+1, D),
                    periodics=periodics, fun=f, maxEval=50000)$integral
  }

  lambda <- c(0.4299832, 0.5700168)
  mu <- list(c(5.005915, 1132.952), c(2.52771, 929.7838))
  Sigma <- list(
    matrix(c(7.953732, -26.48475, -26.484753, 186416.44314), nrow=2, byrow=T),
    matrix(c(13.98248, -752.3661, -752.36607, 471346.9741), nrow=2, byrow=T))
  wraps <- c(1,0)
  wGrid <- matrix(c(
    -18.849556,  8.131516e-20,
    -12.566371,  2.710505e-20,
     -6.283185,  1.355253e-20,
      0.000000,  0.000000e+00,
      6.283185, -1.355253e-20,
     12.566371, -2.710505e-20,
     18.849556, -8.131516e-20), ncol=2, byrow=T)
  model <- structure(list(lambda=lambda, mu=mu, Sigma=Sigma, wraps=wraps, wGrid=wGrid),
                     class=c("swgmm", "gmm"))
  expect_equal(degmm(model), integralSolution(model))
})
