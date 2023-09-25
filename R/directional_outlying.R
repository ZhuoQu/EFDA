#' Description of directional_outlying.R.
#' Obtain the directional outlyingness for univariate functional data or multivariate functional data
#'
#' @param data is either a list of length = number of variables, and each sublist is a nobs * length(timegrid) matrix,
#' or data is a matrix of nobs rows and length(timegrid) columns.
#' @param diroutmatrix is a logical variable of whether to obtain the directional outlyingness. By default,
#' diroutmatrix = FALSE
#' @param depth.dir is a chracter variable in c("RP","MhD","SD","HS"), where "RP" represents the random
#' projection depth, "MhD" represents the Mahalanobis depth, "SD" represents the simplicial depth, and
#' "HS" represents the random halfspace depth. By default, depth.dir = "RP".
#' @param d_value is a logical variable of whether to return the squared Mahalanobis distance
#' @param sq_vo is a logical variable of whether to use the variational outlyingness or the
#' square root of the variational outlyingness. By default, sq_vo is FALSE.
#' @return a list including the index of the mean outlyingness and variational outlyingness
#' @export
#' @examples
#' library(spam);n = 200;p = 4;Sigmainv <- .25^abs(outer(1:n,1:n,"-"))
#' Sigmainv <- as.spam( Sigmainv, eps=1e-4)
#' Sigma <- solve( Sigmainv);data <- rmvnorm(n, mu = rep(0, p), Sigma = Sigma)
#' directional_outlying(data, diroutmatrix = FALSE, depth.dir = "RP", d_value = TRUE, sq_vo = FALSE)

directional_outlying <- function(data, diroutmatrix = FALSE,
                   depth.dir = "RP", d_value = TRUE, sq_vo = FALSE) {
  #####################
  ##source("faccal_num.R")
  ##source("outlying_depth.R")
  #Univariate cases   #
  #####################
  if ("list" %in% class(data)) {
    mat <- data[[1]]
  } else {
    mat <- data
  }
  time_nb <- ncol(mat)
  obs_nb <- nrow(mat)

  if ("matrix" %in% class(data) ) {
    medvec <- apply(mat, 1, median)
    madvec <- apply(mat, 1, mad)
    outmat <- abs((mat - medvec) / (madvec))
    signmat <- sign((mat - medvec))
    dirout <- outmat * signmat
    out_avr <- apply(dirout, 1, FUN = function(y) mean(y, na.rm = TRUE))
    out_var <- apply(dirout, 1, FUN = function(y) var(y, na.rm = TRUE))
  }
  #####################
  #Multivariate cases #
  #####################
  if ("list" %in% class(data)) {
    dirout <- array(0, dim = c(obs_nb, time_nb, length(data)))
    for (j in 1:time_nb) {
      data_mat <- sapply(data, function(mm) {mm[, j]})
      out <- outlying_depth(data_mat, depth.dir)
      me <- data_mat[order(out)[1], ]
      dir <- t(apply(data_mat, 1, function(jm) {
        (jm - me) / (sum((jm - me)^2))^(1 / 2)
      }))
      dir[order(out)[1], ] <- 0
      for (i in 1:obs_nb) {
        dirout[i, j, ] <- out[i] * dir[i, ]
      }
    }
    out_avr <- apply(dirout, c(1, 3), FUN = function(y) mean(y, na.rm = TRUE))
    norm_dif <- sapply(1:time_nb, function(kk) {
      apply(dirout[, kk, ] - out_avr, 1, function(w) {
        sum(w^2)})
    })
    out_var <- apply(norm_dif, 1, mean)
  }

  ##########################################
  ################ Logical Conditions######
  if (sq_vo == TRUE) {
    out_var <- sqrt(out_var)
  }
  result <- list(out_avr = out_avr, out_var = out_var)
  if (diroutmatrix) {
    result <- list.append(result, dirout = dirout) ## function from the package rlist
  }
  if (d_value) {
    matr <- cbind(out_avr, out_var)
    ans <- cov.rob(matr, method = "mcd", nsamp = "best")
    ## Compute a multivariate location and scale estimate
    ## with a high breakdown point – this can be thought
    ## of as estimating the mean and covariance of the
    ## good part of the data.
    cov <- ans$cov ## the final estimate of scatter
    me <- ans$center ## the final estimate of location.
    mhd <- mahalanobis(matr, me, cov)
    ##Returns the squared Mahalanobis distance of all rows in x and the vector \muμ = center with respect to \SigmaΣ = cov. This is (for vector x) defined as
    ##D^2 = (x - \mu)' \Sigma^{-1} (x - \mu)D
    result <- list.append(result, d_value = mhd, cov = cov)
  }
  return(result)
}



