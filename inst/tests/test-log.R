#
# test-log.R
# copyright (c) by Alexander Lehmann <afwlehmann@googlemail.com>
#


context("log")


test_that("logSum yields correct results in every scenario", {
  expect_equal(logSum(log(3), log(4)), log(7))
  expect_equal(logSum(log(4), log(3)), log(7))
  expect_equal(logSum(log(0), log(7)), log(7))
  expect_equal(logSum(log(7), log(0)), log(7))
  expect_equal(logSum(log(0), log(0)), log(0))
  expect_equal(logSum(Inf, Inf)      , Inf   )
  expect_equal(logSum(Inf, -Inf)     , Inf   )
  expect_equal(logSum(-Inf, Inf)     , Inf   )
  expect_equal(logSum(Inf, log(3))   , Inf   )
  expect_equal(logSum(log(3), Inf)   , Inf   )
})
