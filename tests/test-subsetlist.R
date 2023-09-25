# tests/test-subsetlist.R
library(testthat)
context("Testing the subsetlist function")

test_that("extract the subsetlist", {
  # Test case 1
  result1 <-  subsetlist(list(l1 = c(1, 2, 3),
                              l2 = c("d", "efd", "g"),
                              l3 = c(12,345), l4 = 1), c(1, 3))
  expect_length(result1, 2)

  # Test case 2

  result2 <- subsetlist(list(list(center = 1, element = c(1, 2))), 1)
  expect_length(result2, 1)
})

