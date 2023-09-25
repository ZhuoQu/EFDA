# tests/test-msplot_modified.R
library(testthat)
context("Testing the msplot_modified function")

test_that("The result of the msplot_modified function is a list", {
  library(fda)
  library(roahd)
  canada_temp <- t(CanadianWeather$monthlyTemp)
  canada_precip <- t(CanadianWeather$monthlyPrecip)
  ## 1. Univariate Functional Data
  result1 <- msplot_modified(data = canada_temp,
                plot = TRUE, sq_vo = FALSE)
  expect_true(is.list(result1))
  expect_true(is.vector(result1[[1]]),
              info = "The element should be a vector.") ### mo
  expect_true(is.vector(result1[[2]]),
              info = "The element should be a vector.") ### vo

  ### 2. Bivariate Functional Data
  result2 <- msplot_modified(data = list(temp = canada_temp, precip = canada_precip),
                plot = TRUE, sq_vo = FALSE)
  expect_true(is.list(result2))
  expect_true(is.matrix(result2[[1]]),
              info = "The element should be a matrix.") ### mo
  expect_true(is.vector(result2[[2]]),
              info = "The element should be a vector.") ### vo
})
