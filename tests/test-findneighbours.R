# tests/test-findneighbours.R
library(testthat)
context("Testing findneighbours function based on the distance matrix and the threshold")

test_that("Testing findneighbours function based on the distance matrix and the threshold", {
  # Test case 1
  library(fda)
  canada_temp <- t(CanadianWeather$monthlyTemp) ## a matrix of 35 stations and 12 months
  n <- nrow(canada_temp); nt <- ncol(canada_temp); time <- 1:12
  complete_canada.temp_list <- lapply(1:n, function(k) {
    res <- list(argvals = time,
                subj = rep(k, 12), y = canada_temp[k, ])
    return(res)
  })
  dist1 <- elastic_time_distance(complete_canada.temp_list)
  result1 <- findneighbours(dist1)
  expect_true(is.matrix(result1),
              info = "The result should be a matrix.")
  expect_identical(mode(unlist(result1)), "logical",
                   info = "All elements should have a mode of 'logical'.")


  # Test case 2
  canada_temp <- t(CanadianWeather$monthlyTemp) ## a matrix of 35 stations and 12 months
  canada_precip <- t(CanadianWeather$monthlyPrecip) ## a matrix of 35 stations and 12 months
  sparse_canada_list <- lapply(1:n, function(k) {
    observed_loc <- sort(sample.int(12, size = 9, replace = FALSE))
    res <- list(argvals = time[observed_loc],
                subj = rep(k, 9),
                y = cbind(canada_temp[k, observed_loc], canada_precip[k, observed_loc]))
    return(res)
  })
  dist2 <- elastic_time_distance(sparse_canada_list)

  result2 <- findneighbours(dist2)
  expect_true(is.matrix(result2),
              info = "The result should be a matrix.")
  expect_identical(mode(unlist(result2)), "logical",
                   info = "All elements should have a mode of 'logical'.")
})
