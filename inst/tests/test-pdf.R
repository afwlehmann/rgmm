#
# test-pdf.R
# copyright (c) by Alexander Lehmann <afwlehmann@googlemail.com>
#


suppressMessages(library(cubature))


context("pdf")


test_that("dmvnorm handles univariate case", {
  set.seed(1087456098)
  expect_that(length(dmvnorm(1:100)), equals(100),
              info="neither mu nor Sigma")
  expect_that(length(dmvnorm(1:100, mu=3)), equals(100),
              info="mu")
  expect_that(length(dmvnorm(1:100, Sigma=2)), equals(100),
              info="Sigma (as vector)")
  expect_that(length(dmvnorm(1:100, mu=3, Sigma=2)), equals(100),
              info="mu and Sigma (as vector)")
  expect_that(length(dmvnorm(1:100, Sigma=matrix(2))), equals(100),
              info="Sigma (as matrix)")
  expect_that(length(dmvnorm(1:100, mu=3, Sigma=matrix(2))), equals(100),
              info="mu and Sigma (as matrix)")
  expect_equal(dmvnorm(-50:50), dnorm(-50:50),
               info="correct results")
  expect_equal(integrate(dmvnorm, -Inf, +Inf)$value, 1, label="integral")
})


test_that("dmvnorm handles multivariate case (sample points as vector)", {
  set.seed(1087456098)
  tmp <- rnorm(5 * 100)
  # This is only relevant once either or both of mu and Sigma
  # are given since that determines the shape of `tmp`.
  expect_that(length(dmvnorm(tmp, mu=3:7)), equals(100),
              info="mu")
  expect_that(length(dmvnorm(tmp, Sigma=2:6)), equals(100),
              info="Sigma (as vector)")
  expect_that(length(dmvnorm(tmp, mu=3:7, Sigma=2:6)), equals(100),
              info="mu and Sigma (as vector)")
  expect_that(length(dmvnorm(tmp, Sigma=diag(2:6))), equals(100),
              info="Sigma (as matrix)")
  expect_that(length(dmvnorm(tmp, mu=3:7, Sigma=diag(2:6))), equals(100),
              info="mu and Sigma (as matrix)")
  x <- seq(-50, 50, length.out=500)
  y <- seq(-30, 30, length.out=500)
  expect_equal(dmvnorm(as.matrix(expand.grid(x, y))),
               as.vector(outer(x, y, function(a, b) dnorm(a) * dnorm(b))),
               info="correct results")
  if (require(cubature)) {
    # Integration by substitution (phi = x / (1-x^2)).
    subst <- function(x) {
      y <- x^2
      detJ <- prod( (1+y) / (1-y)^2 )
      dmvnorm(x / (1-y), mu=rep(0,2)) * detJ
    }
    expect_equal(adaptIntegrate(subst, rep(-1, 2), rep(1, 2))$integral, 1,
                 label="integral")
  }
})


test_that("dmvnorm handles multivariate case (sample points as matrix)", {
  set.seed(1087456098)
  tmp <- matrix(rnorm(5 * 100), ncol=5)
  expect_that(length(dmvnorm(tmp)), equals(100),
              info="neither mu nor Sigma")
  expect_that(length(dmvnorm(tmp, mu=3:7)), equals(100),
              info="mu")
  expect_that(length(dmvnorm(tmp, Sigma=2:6)), equals(100),
              info="Sigma (as vector)")
  expect_that(length(dmvnorm(tmp, mu=3:7, Sigma=2:6)), equals(100),
              info="mu and Sigma (as vector)")
  expect_that(length(dmvnorm(tmp, Sigma=diag(2:6))), equals(100),
              info="Sigma (as matrix)")
  expect_that(length(dmvnorm(tmp, mu=3:7, Sigma=diag(2:6))), equals(100),
              info="mu and Sigma (as matrix)")
  x <- seq(-50, 50, length.out=500)
  y <- seq(-30, 30, length.out=500)
  expect_equal(dmvnorm(as.matrix(expand.grid(x, y))),
               as.vector(outer(x, y, function(a, b) dnorm(a) * dnorm(b))),
               info="correct results")
  if (require(cubature)) {
    # Integration by substitution (phi = x / (1-x^2)).
    subst <- function(x) {
      if (is.vector(x))
        x <- matrix(x, nrow=1)
      y <- x^2
      detJ <- prod( (1+y) / (1-y)^2 )
      dmvnorm(x / (1-y)) * detJ
    }
    expect_equal(adaptIntegrate(subst, rep(-1, 2), rep(1, 2))$integral, 1,
                 label="integral")
  }
})


test_that("dswmvnorm handles univariate case", {
  set.seed(1087456098)
  expect_that(length(dswmvnorm(1:100)), equals(100),
              info="neither mu nor Sigma")
  expect_that(length(dswmvnorm(1:100, mu=3)), equals(100),
              info="mu")
  expect_that(length(dswmvnorm(1:100, Sigma=2)), equals(100),
              info="Sigma (as vector)")
  expect_that(length(dswmvnorm(1:100, mu=3, Sigma=2)), equals(100),
              info="mu and Sigma (as vector)")
  expect_that(length(dswmvnorm(1:100, Sigma=matrix(2))), equals(100),
              info="Sigma (as matrix)")
  expect_that(length(dswmvnorm(1:100, mu=3, Sigma=matrix(2))), equals(100),
              info="mu and Sigma (as matrix)")
  expect_equal(integrate(dswmvnorm, -pi, +pi)$value, 1, label="integral")
})


test_that("dswmvnorm handles multivariate case (sample points as vector)", {
  set.seed(1087456098)
  tmp <- rnorm(5 * 100)
  # This is only relevant once either or both of mu and Sigma
  # are given since that determines the shape of `tmp`.
  expect_that(length(dswmvnorm(tmp, mu=3:7)), equals(100),
              info="mu")
  expect_that(length(dswmvnorm(tmp, Sigma=2:6)), equals(100),
              info="Sigma (as vector)")
  expect_that(length(dswmvnorm(tmp, mu=3:7, Sigma=2:6)), equals(100),
              info="mu and Sigma (as vector)")
  expect_that(length(dswmvnorm(tmp, Sigma=diag(2:6))), equals(100),
              info="Sigma (as matrix)")
  expect_that(length(dswmvnorm(tmp, mu=3:7, Sigma=diag(2:6))), equals(100),
              info="mu and Sigma (as matrix)")
  if (require(cubature)) {
    aux <- function(x) dswmvnorm(x, mu=rep(0,2))
    expect_equal(adaptIntegrate(aux, rep(-pi, 2), rep(pi, 2))$integral, 1,
                 tol=1e-6, label="integral")
  }
})


test_that("dswmvnorm handles multivariate case (sample points as matrix)", {
  set.seed(1087456098)
  tmp <- matrix(rnorm(5 * 100), ncol=5)
  expect_that(length(dswmvnorm(tmp)), equals(100),
              info="neither mu nor Sigma")
  expect_that(length(dswmvnorm(tmp, mu=3:7)), equals(100),
              info="mu")
  expect_that(length(dswmvnorm(tmp, Sigma=2:6)), equals(100),
              info="Sigma (as vector)")
  expect_that(length(dswmvnorm(tmp, mu=3:7, Sigma=2:6)), equals(100),
              info="mu and Sigma (as vector)")
  expect_that(length(dswmvnorm(tmp, Sigma=diag(2:6))), equals(100),
              info="Sigma (as matrix)")
  expect_that(length(dswmvnorm(tmp, mu=3:7, Sigma=diag(2:6))), equals(100),
              info="mu and Sigma (as matrix)")
  if (require(cubature)) {
    aux <- function(x) {
      if (is.vector(x))
        x <- matrix(x, nrow=1)
      dswmvnorm(x)
    }
    expect_equal(adaptIntegrate(aux, rep(-pi, 2), rep(pi, 2))$integral, 1,
                 tol=1e-6, label="integral")
  }
})


test_that("dgmm yields correct results", {
  model <- structure(
    list(mu=list(c(-14, -14), c(0, 0), c(15, -15)),
         Sigma=list(35*diag(2), 25*diag(2), 18*diag(2)),
         lambda=c(0.3, 0.5, 0.2)),
    class=c("gmm")
    )
  expect_equal(adaptIntegrate(mapInf, rep(-1, 2), rep(+1, 2), maxEval=10000,
                              fun=dgmm, model=model)$integral,
               1, tolerance=1e-5)

  model$mu <- lapply(model$mu, function(m) m %% (2*pi))
  model$Sigma <- lapply(model$Sigma, function(S) S * 1e-2)
  class(model) <- c("swgmm", "gmm")
  expect_equal(adaptIntegrate(dgmm, rep(0, 2), rep(2*pi, 2), maxEval=50000,
                              model=model)$integral,
               1, tolerance=1e-5)
})


test_that("integral of marginalized (SW-)GMM equals 1", {
  model <- structure(
    list(mu=list(c(-14, -14), c(0, 0), c(15, -15)),
         Sigma=list(35*diag(2), 25*diag(2), 18*diag(2)),
         lambda=c(0.3, 0.5, 0.2),
         llh=seq(10)), # fake llh
    class=c("gmm")
    )
  m <- marginal(model, 1)
  expect_equal(integrate(dgmm, -Inf, +Inf, model=m)$value, 1, info="GMM")

  model$mu <- lapply(model$mu, function(m) m %% (2*pi))
  model$Sigma <- lapply(model$Sigma, function(S) S * 1e-2)
  m <- marginal(model, 1)
  class(m) <- c("swgmm", "gmm")
  expect_equal(integrate(dgmm, 0, 2*pi, model=m)$value, 1, info="SWGMM")
})
