#' Description of bootstrap_implementation.R.
#' A function to estimate the univariate functional data/multivariate functional data
#' with bootstrap functional principal component analysis.
#' @param bootstrap_times Bootstrap times. By default, it is 100.
#' @param argval A list of numeric vectors or a single numeric vector, giving the sampling points in the domains. It has the same
#' definition as in the funData function in the package "funData".
#' @param true_data If data is a simulated functional data, true_data is a list of p (the number of variables), and each sublist is a matrix
#' with n (the number of observations) rows and length(argval) columns. If data comes from the real data with missing values, then true_data = NULL.
#' By default, true_data = NULL.
#' @param sparse_data sparse_data is a list of p (the number of variables), and each sublist is a matrix
#' with n (the number of observations) rows and length(argval) columns. NA is allowed in each sublist.
#' @param confid_alpha represents the level of significance or the probability of making a Type I error. By default, confid_alpha = 0.05.
#' @return The result is a list of 3 elements.
#' 1. fit. This is the same data structure as the sparse_data in the argument.
#' 2. sample_index. This is the sampling index from the population. If isbootstrap = False, sample_index is 1:n. Otherwise, it is the random sampling from
#' the population with replacement.
#' 3. mrse: is the summation of square of residuals / summation of square of the original data.
#' @export
#' @examples
## 1. Univariate sparse functional data
#' library(refund); list("cd4"); argval <- -18:42;
#' bootstrap_cd4 <- bootstrap_implementation(bootstrap_times = 1, argval = -18:42, true_data = NULL, sparse_data = cd4, confid_alpha = 0.05)
## 2. Bivariate sparse functional data
#' library(fda)
#' canada_temp <- t(CanadianWeather$monthlyTemp) ## a matrix of 35 stations and 12 months
#' canada_precip <- t(CanadianWeather$monthlyPrecip) ## a matrix of 35 stations and 12 months
#' true_data <- list(temp = canada_temp, precip = canada_precip);
#' n <- nrow(canada_temp); nt <- ncol(canada_temp); missing_prob <- 0.3; argval <- 1:12
#' random_matrix <- matrix(runif(n * nt), nrow = n, ncol = nt)
## Create a mask based on the missingness probability
#' missing_mask <- random_matrix <= missing_prob
#' sparse_data <- lapply(true_data, function(k) {
#'  k[missing_mask] <- NA
#'  return(k)
#' })
## Replace values in the full matrix with NA where the mask is TRUE
#' bt_canada <- bootstrap_implementation(bootstrap_times = 1, argval, true_data, sparse_data, confid_alpha = 0.05)


bootstrap_implementation <- function(bootstrap_times = 100, argval,
                                     true_data,
                                     sparse_data, confid_alpha = 0.05) {
  ##source("bootstrap_mfpca.R")
  if ("matrix" %in% class(sparse_data)) {
    p <- 1
    n <- nrow(sparse_data)
  } else if ("list" %in% class(sparse_data)) {
    p <- length(sparse_data)
    n <- nrow(sparse_data[[1]])
  } else {
    stop ("sparse_data needs to be either a matrix (displaying univariate functional data) or
          a list (displaying multivariate functional data)!!!")
  }
  mfpca_full <- bootstrap_mfpca(argval, true_data,
                                #observed_data,
                                sparse_data, isbootstrap = FALSE)
  if (length(true_data) == 0) {
    true_data <- mfpca_full$fit
  }

  train_bootstrap <- lapply(1:bootstrap_times, function(bt) {
    cat(bt, "\n")
    mfpca <- bootstrap_mfpca(argval, true_data, sparse_data, isbootstrap = TRUE)
    mfpca_fit <- lapply(1:p, function(k) {
    temp_mat <- matrix(NA, nrow = n, ncol = length(argval))
    temp_mat[mfpca$sample_index, ] <- mfpca$fit[[k]]
    return (temp_mat)
    })
  curve_specific_bias <- function(obs_index) {
    lapply(1:p, function(k) {
      bias <- mfpca_fit[[k]][obs_index, ] - mfpca_full$fit[[k]][obs_index, ]
      return(bias)
      })
  }

  curve_specific_covariance <- function(obs_index) {
    lapply(1:p, function(k) {
      measurement_var <- var(mfpca_fit[[k]][obs_index, ] - mfpca_full$fit[[k]][obs_index, ])
    })
  }
  exp_res <- lapply(1:n, curve_specific_bias)
    return(list(model_dif_cov = mfpca$model_error_variance,
                fit = mfpca_fit,
                exp_res = exp_res))
  })

  bootstrap_fit <- lapply(1:p, function(l) {
    fit <- matrix(0, nrow = n, ncol = length(argval))
    for (obs in 1:n) {
      temp_mat <- t(sapply(train_bootstrap, function(k) {
        k$fit[[l]][obs, ]}))
      fit_value <- apply(temp_mat, 2, function(w) {
        mean(w, na.rm = TRUE)})
      fit[obs, ] <- fit_value
    }
    return(fit)
  })

  bootstrap_naive_lower <- lapply(1:p, function(l) {
    fit <- matrix(0, nrow = n, ncol = length(argval))
    for (obs in 1:n) {
      temp_mat <- t(sapply(train_bootstrap, function(k) {
      k$fit[[l]][obs, ]}))
      fit_value <- apply(temp_mat, 2, function(w) {
       quantile(w, prob = confid_alpha / 2, na.rm = TRUE)})
      fit[obs, ] <- fit_value
    }
    return(fit)
  })

  bootstrap_naive_upper <- lapply(1:p, function(l) {
    fit <- matrix(0, nrow = n, ncol = length(argval))
    for (obs in 1:n) {
      temp_mat <- t(sapply(train_bootstrap, function(k) {
        k$fit[[l]][obs, ]}))
      fit_value <- apply(temp_mat, 2, function(w) {
        quantile(w, prob = 1 - confid_alpha / 2, na.rm = TRUE)})
      fit[obs, ] <- fit_value
    }
    return(fit)
  })

  mrse_numerator <- lapply(1:p, function(k) {apply(bootstrap_fit[[k]] - true_data[[k]], 1, function(m) {sum(m^2)})})
  mrse_denominator <- lapply(1:p, function(k) {apply(true_data[[k]], 1, function(m) {sum(m^2)})})
  numer <- mrse_numerator[[1]]
  denom <- mrse_denominator[[1]]
  if (p > 1) {
    for (nb in 2:p) {
      numer <- numer + mrse_numerator[[nb]]
      denom <- denom + mrse_denominator[[nb]]
    }
  }
  mrse <- mean(numer / denom)

  return(list(bootstrap_fit = bootstrap_fit,
              full_fit = mfpca_full$fit,
              bootstrap_mrse = mrse, full_mrse = mfpca_full$mrse,
              CI_lower = bootstrap_naive_lower, CI_upper = bootstrap_naive_upper))
}


