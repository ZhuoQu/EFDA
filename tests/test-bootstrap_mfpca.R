# tests/test-bootstrap_mfpca.R
library(testthat)
context("Testing bootstrap mfpca.R")

test_that("Testing the bootstrap mfpca method is a list", {
  ## 1. Univariate sparse functional data
  library(refund); list("cd4"); argval <- -18:42;
  result1 <- bootstrap_mfpca(argval, true_data = NULL, sparse_data = cd4, isbootstrap = TRUE)
  expect_true(is.list(result1))
  ## 2. Bivariate sparse functional data
   library(fda)
   canada_temp <- t(CanadianWeather$monthlyTemp) ## a matrix of 35 stations and 12 months
   canada_precip <- t(CanadianWeather$monthlyPrecip) ## a matrix of 35 stations and 12 months
   true_data <- list(temp = canada_temp, precip = canada_precip);
   n <- nrow(canada_temp); nt <- ncol(canada_temp); missing_prob <- 0.3; argval <- 1:12
   random_matrix <- matrix(runif(n * nt), nrow = n, ncol = nt)
  ## Create a mask based on the missingness probability
   missing_mask <- random_matrix <= missing_prob
   sparse_data <- lapply(true_data, function(k) {
    k[missing_mask] <- NA
    return(k)
   })
  ## Replace values in the full matrix with NA where the mask is TRUE
   result2 <- bootstrap_mfpca(argval, true_data,sparse_data, isbootstrap = TRUE)
   expect_true(is.list(result2))
})
