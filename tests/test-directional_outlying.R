library(testthat)# tests/test-directional_outlying.R

context("Testing the directional outlyingness function")

test_that("The result of the directional outlyingness function is a list", {
  ## Test case 1
  library(spam);n = 200; timegrid = 6;Sigmainv1 <- .25^abs(outer(1:timegrid,1:timegrid,"-"))
  Sigmainv1 <- as.spam(Sigmainv1, eps=1e-4)
  Sigma1 <- solve(Sigmainv1); data1 <- rmvnorm(n, mu = rep(0, timegrid), Sigma = Sigma1)
  Sigmainv2 <- .4^abs(outer(1:timegrid,1:timegrid,"-"))
  Sigmainv2 <- as.spam(Sigmainv2, eps=1e-4)
  Sigma2 <- solve(Sigmainv2); data2 <- mvrnorm(n, mu = rep(0, timegrid), Sigma = Sigma2)
  data <- list(data1, data2)
  result1 <- directional_outlying(data, diroutmatrix = FALSE, depth.dir = "RP", d_value = TRUE, sq_vo = FALSE)

  expect_true(is.list(result1))
  expect_true(is.matrix(result1[[1]])) ## mean outlyingness
  expect_length(result1[[2]], n) ## variational outlyingness
})


