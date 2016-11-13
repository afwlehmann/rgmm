#
# test-pdf.R
# copyright (c) by Alexander Lehmann <afwlehmann@googlemail.com>
#


suppressMessages(library(cubature))


context("sampling")


test_that("rmvnorm draws the desired number of samples", {
  expect_equal(nrow(rmvnorm(10)), 10)
})


test_that("rmvnorm handles mu", {
  set.seed(1087456098)
  expect_equal(mean(rmvnorm(500000, mu=5)), 5, tol=1e-4,
               info="univariate")
  expect_equal(colMeans(rmvnorm(500000, mu=c(-3,4))), c(-3,4), tol=1e-3,
               info="multivariate")
})


test_that("rmvnorm handles Sigma", {
  set.seed(1087456098)
  # Univariate sigma (this is basically the same as given as a vector).
  expect_equal(var(rmvnorm(500000, Sigma=15)), as.matrix(15), tol=1e-3,
               info="univariate")
  # Sigma given as a vector.
  tmp <- rmvnorm(2500000, Sigma=c(2, 5))
  expect_equal(cov(tmp), diag(c(2, 5)), tol=1e-2,
               info="multivariate w/ Sigma given as a vector")
  # Sigma given as a matrix.
  tmp <- rmvnorm(2500000, Sigma=diag(c(2, 5)))
  expect_equal(cov(tmp), diag(c(2, 5)), tol=1e-2,
               info="multivariate w/ Sigma given as a matrix")
})


test_that("rgmm works as expected", {
  set.seed(1087456098)

  lambda <- c(0.3, 0.7)
  mu     <- list(c(-10, 4), c(3, -5))
  Sigma  <- list(diag(c(3, 4)), diag(2))

  model  <- structure(list(lambda=lambda, mu=mu, Sigma=Sigma), class=c("gmm"))
  X <- rgmm(100000, model)
  modelPrime <- computeGMM(X, K=2, lambda=lambda, mu=mu, Sigma=Sigma, verb=F)

  createIntegrand <- function(modelP, modelQ) {
    function(x) {
      #stopifnot(length(x) == 2L)
      #x <- matrix(x, nrow=1)
      p   <- dgmm(x, modelP)
      q   <- dgmm(x, modelQ)
      lpq <- log(p) - log(q)
      # Note that is.finite is FALSE for both +/- Inf _AND_ NaN.
      # This is important because -Inf - (-Inf) returns NaN!
      if (is.finite(lpq)) lpq * p else 0
    }
  }

  klDivPQ <- adaptIntegrate(createIntegrand(model, modelPrime),
                            rep(-100, 2), rep(+100, 2), maxEval=5000)$integral
  klDivQP <- adaptIntegrate(createIntegrand(modelPrime, model),
                            rep(-100, 2), rep(+100, 2), maxEval=5000)$integral
  expect_equal(klDivPQ, klDivQP, tolerance=1e-4, info="KL(P||Q) == KL(Q||P)")
})
