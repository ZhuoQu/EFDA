# tests/test-outlying_depth.R
library(testthat)
context("Testing the outlying_depth function")

test_that("The result of the outlying_depth function is a vector", {
  # Test case 1
  library(spam);n = 200;p = 4;Sigmainv <- .25^abs(outer(1:p,1:p,"-"))
  Sigmainv <- as.spam( Sigmainv, eps=1e-4)
  Sigma <- solve(Sigmainv);data <- rmvnorm(n, mu = rep(0, p), Sigma = Sigma)
  result1 <- outlying_depth(data, "RP")
  expect_true(is.vector(result1)) ##
  expect_length(result1, n)
})


