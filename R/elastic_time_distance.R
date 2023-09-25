### We assume that pd is a list. Each sublist is a list with three components: argvals, subj and y ( a matrix of length|t| * variable numbers)
#' Description of elastic_time_distance.R.
#' This function returns the distance matrix for the multivariate functional data.
#' #' See Robust two-layer partition clustering of sparse multivariate functional data,
#' Zhuo Qu, Wenlin Dai, and Marc G. Genton, https://www.sciencedirect.com/science/article/pii/S2452306223000138
#' This function implements algorithm 1 in the paper.
#'
#' @param pd A list. Each sublist is a list with three components: argvals, subj and y. If the data is univariate, y is a vector.
#' If the data is multivariate, then y is a matrix with |argvals| rows and |variables| columns.
#' @return a non-negative, symmetric distance matrix
#' @export
#' @examples
## 1. Univariate complete functional data
#' library(fda); canada_temp <- t(CanadianWeather$monthlyTemp) ## a matrix of 35 stations and 12 months
#'  n <- nrow(canada_temp); nt <- ncol(canada_temp); time <- 1:12
#'  complete_canada.temp_list <- lapply(1:n, function(k) {
#'  res <- list(argvals = time,
#'              subj = rep(k, 12), y = canada_temp[k, ])
#'   return(res)
#'   })
#'   dist_complete_canada.temp <- elastic_time_distance(complete_canada.temp_list)
#'
#'
##  2. Univariate sparse functional data
#'library(fda); canada_temp <- t(CanadianWeather$monthlyTemp) ## a matrix of 35 stations and 12 months
#'  n <- nrow(canada_temp); nt <- ncol(canada_temp); missing_prob <- 0.3; time <- 1:12
#'  random_matrix <- matrix(runif(n * nt), nrow = n, ncol = nt)
##  Create a mask based on the missingness probability
#'  missing_mask <- random_matrix <= missing_prob
#'  sparse_canada.temp_list <- lapply(1:n, function(k) {
#'  observed <- which(missing_mask[k, ] != TRUE)
#'  res <- list(argvals = time[observed],
#'              subj = rep(k, length(observed)), y = canada_temp[k, observed])
#'   return(res)
#'   })
#'   dist_sparse_canada.temp <- elastic_time_distance(sparse_canada.temp_list)
#'
##  3. Bivariate complete functional data
#' canada_temp <- t(CanadianWeather$monthlyTemp) ## a matrix of 35 stations and 12 months
#'  canada_precip <- t(CanadianWeather$monthlyPrecip) ## a matrix of 35 stations and 12 months
#'  complete_canada_list <- lapply(1:n, function(k) {
#'  res <- list(argvals = time,
#'              subj = rep(k, 12), y = cbind(canada_temp[k, ], canada_precip[k, ]))
#'   return(res)
#'   })
#'   dist_complete_canada <- elastic_time_distance(complete_canada_list)

##  4. Bivariate sparse functional data
#' canada_temp <- t(CanadianWeather$monthlyTemp) ## a matrix of 35 stations and 12 months
#'  canada_precip <- t(CanadianWeather$monthlyPrecip) ## a matrix of 35 stations and 12 months
#'  sparse_canada_list <- lapply(1:n, function(k) {
#'  observed_loc <- sort(sample.int(12, size = 9, replace = FALSE))
#'  res <- list(argvals = time[observed_loc],
#'              subj = rep(k, 9), y = cbind(canada_temp[k, observed_loc], canada_precip[k, observed_loc]))
#'   return(res)
#'   })
#'   dist_sparse_canada <- elastic_time_distance(sparse_canada_list)


elastic_time_distance <- function(pd) {

  DO <- lapply(pd, function(k) {k$y})
  n <- length(pd)
  if ( "numeric" %in% class(pd[[1]]$y) ) {
    nbvar <- 1
  } else {
    nbvar <- ncol(pd[[1]]$y)
  }
  ## define the standard time grid
  max_time <- max(unlist(sapply(pd, function(k) {k$argvals})))
  max_grid_length <- max(sapply(pd, function(k) {length(k$argvals)}))
  sample_time <- seq(1, max_time, length.out = max_grid_length) / max_time

  ## define the time grid for each subject
  time_warping_index <- function(index) {
    time_orig <- pd[[index]]$argvals / max_time
    time_warping <- sapply(1:length(sample_time), function(k) {
      ind <- which(
        abs(time_orig-sample_time[k]) == min(abs(time_orig-sample_time[k]))
      )
      if (length(ind) > 1) {
        ind <- ind[1]
      }
      return (ind)
    })
    return (time_warping)
  }
  ## define the distance between the rowindex-th curve and the colindex-th curve [rowindex, colindex]
  distance_warping <- function(rowindex, colindex) {
    time_r_index <- time_warping_index(rowindex)
    time_c_index <- time_warping_index(colindex)
    if (nbvar > 1) {
      difference <- DO[[rowindex]][time_r_index, ] - DO[[colindex]][time_c_index, ]
      if ("matrix" %in% class(difference)) {
        val <- max(apply(difference, 1, function(k) { sqrt(sum(k^2)) } ), na.rm = TRUE)
      } else {
        val <- sqrt(sum(difference^2))
      }
    } else {
      difference <- DO[[rowindex]][time_r_index] - DO[[colindex]][time_c_index]
      val <- max(abs(difference))
    }
    return (val)
  }

  ############################# calculate the distance matrix
  upt <- lapply(1:(n-1), function (nr) {
    sapply((nr+1):n, function (nc) {
      distance_warping(nr, nc)
    })
  })

  numer <- unlist(upt)

  distance_matrix <- matrix(0, nrow = n, ncol = n)

  distance_matrix[lower.tri(distance_matrix, diag = FALSE)] <- numer

  distance_matrix <- t(distance_matrix)

  lowerTriangle(distance_matrix) <- upperTriangle(distance_matrix, byrow = TRUE)

  return (distance_matrix)
}
