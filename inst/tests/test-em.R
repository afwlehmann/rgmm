#
# test-em.R
# copyright (c) by Alexander Lehmann <afwlehmann@googlemail.com>
#


context("em")


test_that("computeGMM works on univariate data", {
  set.seed(1231236713)
  
  # K=1
  model <- computeGMM(rnorm(100000, mean=2, sd=1.5), K=1, verb=F)
  expect_equal(length(model$lambda), 1, info="K=1")
  expect_equal(length(model$mu), 1, info="K=1")
  expect_equal(length(model$Sigma), 1, info="K=1")
  expect_equal(model$lambda, 1.0, tolerance=0, info="K=1")
  expect_equal(model$mu[[1]], 2.0, tolerance=1e-3, info="K=1")
  expect_equal(model$Sigma[[1]], matrix(2.25), tolerance=1e-2, info="K=1")
  
  # K=2
  X <- c(rnorm(100000, mean=0), rnorm(50000, mean=10, sd=2))
  model <- computeGMM(X, K=2, epsilon=1e-10, verb=F)
  expect_equal(length(model$lambda), 2, info="K=2")
  expect_equal(length(model$mu), 2, info="K=2")
  expect_equal(length(model$Sigma), 2, info="K=2")
  expect_equal(model$lambda, c(2/3, 1/3), tolerance=1e-2, info="K=2")
  expect_equal(model$mu, list(0.0, 10.0), tolerance=1e-2, info="K=2")
  expect_equal(model$Sigma,
               list(matrix(1), matrix(4)), tolerance=1e-1, info="K=2")
})


test_that("computeGMM works on multivariate data", {
  set.seed(1231236713)

  # K = 1
  model <- computeGMM(rmvnorm(100000, mu=c(5, 3), Sigma=diag(c(2, 0.5))),
                      K=1, verb=F)
  expect_equal(length(model$lambda), 1, info="K=1")
  expect_equal(length(model$mu), 1, info="K=1")
  expect_equal(length(model$Sigma), 1, info="K=1")
  expect_equal(model$lambda, 1.0, tolerance=0, info="K=1")
  expect_equal(model$mu, list(c(5.0, 3.0)), tolerance=1e-3, info="K=2")
  expect_equal(model$Sigma,
               list(diag(c(2.0, 0.5))), tolerance=1e-2, info="K=2")

  # K = 2
  X <- rbind(rmvnorm(100000, mu=c(0,0)),
             rmvnorm(50000, mu=c(5, 3), Sigma=diag(c(2, 0.5))))
  model <- computeGMM(X, K=2, epsilon=1e-10, verb=F)
  expect_equal(length(model$lambda), 2, info="K=2")
  expect_equal(length(model$mu), 2, info="K=2")
  expect_equal(length(model$Sigma), 2, info="K=2")
  expect_equal(model$lambda, c(2/3, 1/3), tolerance=1e-2, info="K=2")
  expect_equal(model$mu,
               list(c(0,0), c(5,3)), tolerance=1e-2, info="K=2")
  expect_equal(model$Sigma,
               list(diag(2), diag(c(2, 0.5))), tolerance=1e-1, info="K=2")
})


# FIXME: Temporarily disabled
# test_that("EM issues warning on singularities", {
#   X <- c(rnorm(10, mean=0), 10)
#   expect_that(computeGMM(X, K=2,
#                          lambda=c(.8, .2),
#                          mu=list(0, 10),
#                          Sigma=list(matrix(1), matrix(1)),
#                          verb=F),
#               gives_warning())
# })
