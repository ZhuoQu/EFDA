# tests/test-faccal_num.R
library(testthat)
context("Testing the factor and the cutoff value given the population size and its components")

test_that("the factor and cutoff value given n and p", {
  ## Test case 1
  result1 <- faccal_num(1000, 5)
  expect_length(result1, 2)

  ## Test case 2
  result2 <- faccal_num(100, 3)
  expect_length(result2, 2)
})
