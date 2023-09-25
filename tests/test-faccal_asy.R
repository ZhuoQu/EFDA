# tests/test-faccal_asy.R
library(testthat)
context("Testing the factor and the cutoff value of the approximate F distribution function")

test_that("the factor and cutoff value of the approximate F distribution function", {
  # Test case 1
  result1 <- faccal_asy(1000, 5)
  expect_length(result1, 2)

  # Test case 2
  result2 <- faccal_asy(100, 3)
  expect_length(result2, 2)
})
