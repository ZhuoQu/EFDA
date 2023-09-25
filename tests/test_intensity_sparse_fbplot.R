# tests/test-intensity_sparse_fbplot.R
library(testthat)
context("Testing the intensity sparse functional boxplot function")

test_that("The result of the intensity sparse functional boxplot is a list", {
  ##  1. Univariate complete functional data
  library(fda)
  library(roahd)
  canada_temp <- t(CanadianWeather$monthlyTemp)
  result1 <- intensity_sparse_fbplot(fit = canada_temp,
                           sparse = canada_temp,
                           time_index = 1:12,
                           depth = MBD(canada_temp))
  expect_true(is.list(result1))

  ##  2. Univariate sparse functional data
  ##library(refund); list("cd4"); argval <- -18:42;
  ##bootstrap_cd4 <- bootstrap_implementation(bootstrap_times = 50, argval = -18:42, true_data = NULL, sparse_data = cd4, confid_alpha = 0.05)
  ##fit_cd4 <- bootstrap_cd4$bootstrap_fit[[1]]
  ##depth_cd4 <- MBD(fit_cd4) ## modified band depth
  ##result2 <- intensity_sparse_fbplot(fit_cd4, sparse = cd4, time_index = NULL, depth = depth_cd4)
  ##expect_true(is.list(result2))
  ###  3. Bivariate complete functional data
  library(fda); library(mrfDepth)
  canada_temp <- t(CanadianWeather$monthlyTemp) ## a matrix of 35 stations and 12 months
  canada_precip <- t(CanadianWeather$monthlyPrecip) ## a matrix of 35 stations and 12 months
  true_canada <- list(temp = canada_temp, precip = canada_precip);
  n <- nrow(canada_temp); nt <- ncol(canada_temp); missing_prob <- 0.3; argval <- 1:12
  canada_array <- array(NA, dim = c(nt, n, 2)); canada_array[, , 1] <- t(canada_temp); canada_array[, , 2] <- t(canada_precip);
  mfd_canada <- mfd(x = canada_array, z = canada_array, type = "sdepth")
  result3 <- intensity_sparse_fbplot(fit = true_canada,
                                     sparse = true_canada,
                                     time_index = NULL,
                                     depth = as.numeric(mfd_canada$MFDdepthZ))
  expect_true(is.list(result3))

  ###  4. Bivariate sparse functional data
  ##canada_temp <- t(CanadianWeather$monthlyTemp) ## a matrix of 35 stations and 12 months
  ##canada_precip <- t(CanadianWeather$monthlyPrecip) ## a matrix of 35 stations and 12 months
  ##true_canada <- list(temp = canada_temp, precip = canada_precip);
  ##n <- nrow(canada_temp); nt <- ncol(canada_temp); missing_prob <- 0.3; argval <- 1:12
  ##random_matrix <- matrix(runif(n * nt), nrow = n, ncol = nt)
  ## Create a mask based on the missingness probability
  ##missing_mask <- random_matrix <= missing_prob
  ##sparse_canada <- lapply(true_canada, function(k) {
  ##  k[missing_mask] <- NA
  ##  return(k)
  ##})
  ###   Replace values in the full matrix with NA where the mask is TRUE
  ##bt_canada <- bootstrap_implementation(bootstrap_times = 100, argval, true_data, sparse_data, confid_alpha = 0.05)
  ##fit_canada <- bt_canada$bootstrap_fit
  ##fit_canada_array <- array(NA, dim = c(nt, n, 2));
  ##fit_canada_array[, , 1] <- t(fit_canada[[1]]); fit_canada_array[, , 2] <- t(fit_canada[[2]]);
  ##mfd_fit_canada <- mfd(x = fit_canada_array, z = fit_canada_array, type = "sdepth")
  ##result4 <- intensity_sparse_fbplot(fit = fit_canada, sparse = sparse_canada,
  ##                         time_index = NULL, depth = as.numeric(mfd_fit_canada$MFDdepthZ))
  ##expect_true(is.list(result4))

})
