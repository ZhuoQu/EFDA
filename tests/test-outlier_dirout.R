# tests/test-outlier_dirout.R
library(testthat)
context("Testing the outlier_dirout function")

test_that("The result of the outlier_dirout is a list", {
  ## Test case 1
  library(spam);n = 200; timegrid = 6;Sigmainv1 <- .25^abs(outer(1:timegrid,1:timegrid,"-"))
  Sigmainv1 <- as.spam(Sigmainv1, eps=1e-4)
  Sigma1 <- solve(Sigmainv1); data1 <- rmvnorm(n, mu = rep(0, timegrid), Sigma = Sigma1)
  Sigmainv2 <- .4^abs(outer(1:timegrid,1:timegrid,"-"))
  Sigmainv2 <- as.spam(Sigmainv2, eps=1e-4)
  Sigma2 <- solve(Sigmainv2); data2 <- rmvnorm(n, mu = rep(0, timegrid), Sigma = Sigma2)
  data <- list(data1, data2)

  result1 <- outlier_dirout(data, sq_vo = FALSE, depth.dir = "RP")
  expect_true(is.list(result1))
  expect_true(is.vector(result1[[1]]),
              info = "The element should be a vector.") ### outlier
  expect_true(is.vector(result1[[2]]),
              info = "The element should be a vector.") ### median
})
