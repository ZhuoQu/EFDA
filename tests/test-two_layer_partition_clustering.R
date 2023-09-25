# tests/test-two_layer_partition_clustering.R
library(testthat)
context("Testing the two layer partition clustering method")

test_that("Testing the result of the two layer partition clustering method based on the logical neighbour matrix", {
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
  neighbour_matrix1 <- findneighbours(dist1)

  result1 <- two_layer_partition_clustering(neighbour_matrix1)

  expect_is(result1, "list",
            info = "The result should be a list.")

  expect_length(result1, 2)

  expect_gt(length(result1[[1]]), length(result1[[2]]))

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

  neighbour_matrix2 <- findneighbours(dist2)
  result2 <- two_layer_partition_clustering(neighbour_matrix2)

  expect_is(result2, "list",
            info = "The result should be a list.")

  expect_length(result2, 2)
  expect_gt(length(result2[[1]]), length(result2[[2]]))
})
