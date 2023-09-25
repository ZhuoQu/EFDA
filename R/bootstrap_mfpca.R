#' Description of bootstrap_mfpca.R.
#' A function to estimate the univariate functional data/multivariate functional data with the time grid and the original data.
#' @param argval A list of numeric vectors or a single numeric vector, giving the sampling points in the domains. It has the same
#' definition as in the funData function in the package "funData".
#' @param true_data If data is a simulated functional data, true_data is a list of p (the number of variables), and each sublist is a matrix
#' with n (the number of observations) rows and length(argval) columns. If data comes from the real data with missing values, then true_data = NULL.
#' By default, true_data = NULL.
#' @param sparse_data sparse_data is a list of p (the number of variables), and each sublist is a matrix
#' with n (the number of observations) rows and length(argval) columns. NA is allowed in each sublist.
#' @param isbootstrap is a logical variable of whether to apply bootstrap in the functional principal component analysis.
#' @return The result is a list of 3 elements.
#' 1. fit. This is the same data structure as the sparse_data in the argument.
#' 2. sample_index. This is the sampling index from the population. If isbootstrap = False, sample_index is 1:n. Otherwise, it is the random sampling from
#' the population with replacement.
#' 3. mrse: is the summation of square of residuals / summation of square of the original data.
#' @export
#' @examples
## 1. Univariate sparse functional data
#' library(refund); list("cd4"); argval <- -18:42;
#' bootstrap_mfpca(argval, true_data = NULL, sparse_data = cd4, isbootstrap = TRUE)
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
#' bootstrap_mfpca(argval, true_data,sparse_data, isbootstrap = TRUE)

##################################################################################
bootstrap_mfpca <- function(argval, true_data = NULL,
                            sparse_data, isbootstrap) {
  if ("matrix" %in% class(sparse_data)) {
  sparse_data <- list(sparse_data)
  }
  nc <- ncol(sparse_data[[1]])
  p <- length(sparse_data)
  n <- nrow(sparse_data[[1]])

  obj <- lapply(1:p, function(i) {
    Ly <- sparse_data[[i]]
    res <- funData::funData(list(argval), Ly)
    return(res) }  )
  expres <- list(type = "uFPCA")
  mFData <- multiFunData(obj)
  newmFData <- mFData

  sample_index <- 1:n
  estimated_observed_p <- sapply(1:p, function(i) {
    apply(sparse_data[[i]], 1, function(j) {
      1 - sum(is.na(j)) / nc})})

  min_observed_points <- nc * min(estimated_observed_p)

  M <- ifelse(min_observed_points <= 3, 4, 9)
  if (length(true_data) == 0) {
    fit_mfpca <- MFPCA::MFPCA(newmFData, M,
                              uniExpansions = lapply(1:p, function(i) {expres}),
                              fit = TRUE)
    true_data <- lapply(fit_mfpca$fit, function(k) {k@X})
  }
  new_true_data <- true_data
  newsparse_data <- sparse_data

  if (isbootstrap == TRUE) {
    sample_nb <- base::sample(1:n, n, replace = TRUE)
    # sample_nb <- base::sample(which(apply(sparse_data[[1]], 1, function(k){sum(!is.na(k))}) > 5),
    #                           round(0.9 * n), replace = TRUE)
    # sample_nb <- c(sample_nb, sample(which(apply(sparse_data[[1]], 1, function(k){sum(!is.na(k))}) < 5),
    #                                   n - round(0.9 * n), replace = TRUE))
    sample_index <- sort(sample_nb)

    for (i in 1:p) {
      newmFData[[i]]@X <- mFData[[i]]@X[sample_index, ]
      newsparse_data[[i]] <- sparse_data[[i]][sample_index, ]
      new_true_data[[i]] <- true_data[[i]][sample_index, ]
      #new_observed_data[[i]] <- observed_data[[i]][sample_index, ]
    }
  }

  estimated_observed_p <- sapply(1:p, function(i) {
    apply(sparse_data[[i]], 1, function(j) {
      1 - sum(is.na(j)) / nc})})

  fit_mfpca <- MFPCA::MFPCA(newmFData, M,
                    uniExpansions = lapply(1:p, function(i) {expres}),
                    fit = TRUE)

  eigenfun <- fit_mfpca$functions
  gamma <- diag(fit_mfpca$values)
  fitY <- lapply(1:p, function(i) {fit_mfpca$fit[[i]]@X})
  #residual <- lapply(1:p, function(i) {fitY[[i]] - new_observed_data[[i]]}) ## fit minus observation with epsilon
  estimation_error <- lapply(1:p, function(i) {fitY[[i]] - new_true_data[[i]]}) ## fit minus true data with out epsilon

  mrse_numerator <- lapply(estimation_error, function(k) {apply(k, 1, function(m) {sum(m^2)})})
  if ("list" %in% class(new_true_data)) {
    mrse_denominator <- lapply(new_true_data, function(k) {apply(k, 1, function(m) {sum(m^2)})})
  } else {
    mrse_denominator <- apply(new_true_data, 1, function(m) {sum(m^2)})
    mrse_denominator <- list(mrse_denominator)
  }
  numer <- mrse_numerator[[1]]
  denom <- mrse_denominator[[1]]
  if (p >= 2) {
    for (nb in 2:p) {
      numer <- numer + mrse_numerator[[nb]]
      denom <- denom + mrse_denominator[[nb]]
    }
  }
  mrse <- mean(numer / denom)

  result <- list(fit = fitY, ### p lists of n * length(sett) matrix
                 #residual = residual,
                 sample_index = sample_index,
                 mrse = mrse)
  ### fit = score %*% Phi + mean
  return (result)
}

