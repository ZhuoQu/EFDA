# tests/test-bootstrap_implementation.R
library(testthat)
context("Testing bootstrap_implementation.R")

test_that("Testing the result of the bootstrap_implementation function is a list", {
  ## 1. Univariate sparse functional data
   library(refund); list("cd4"); argval <- -18:42;
   bootstrap_cd4 <- bootstrap_implementation(bootstrap_times = 100, argval = -18:42, true_data = NULL, sparse_data = cd4, confid_alpha = 0.05)
   expect_true(is.list(bootstrap_cd4))
   ## 2. Bivariate sparse functional data
   library(fda)
   canada_temp <- t(CanadianWeather$monthlyTemp) ## a matrix of 35 stations and 12 months
   canada_precip <- t(CanadianWeather$monthlyPrecip) ## a matrix of 35 stations and 12 months
   true_data <- list(temp = canada_temp, precip = canada_precip);
   n <- nrow(canada_temp); nt <- ncol(canada_temp); missing_prob <- 0.2; argval <- 1:12
   random_matrix <- matrix(runif(n * nt), nrow = n, ncol = nt)
   # Create a mask based on the missingness probability
   missing_mask <- random_matrix <= missing_prob
   sparse_data <- lapply(true_data, function(k) {
    k[missing_mask] <- NA
    return(k)
   })
  # Replace values in the full matrix with NA where the mask is TRUE
   bootstrap_canada <- bootstrap_implementation(bootstrap_times = 10, argval, true_data, sparse_data, confid_alpha = 0.05)
  expect_true(is.list(bootstrap_canada))
})
