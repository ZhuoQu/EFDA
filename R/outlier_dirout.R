#' Description of outlier_dirout.R.
#' Obtain the index of outlier, index of median curve and the measure of outlyingness
#' with directional outlyingness method
#' for univariate functional data or multivariate functional data.
#'
#' @param data is either a list of length p, and each sublist is a nobs * length(timegrid) matrix,
#' or data is a matrix of nobs rows and length(timegrid) columns.
#' @param sq_vo is a logical variable of whether to use the variational outlyingness or the
#' square root of the variational outlyingness. By default, sq_vo is FALSE.
#' @param depth.dir is a chracter variable in c("RP","MhD","SD","HS"), where "RP" represents the random
#' projection depth, "MhD" represents the Mahalanobis depth, "SD" represents the simplicial depth, and
#' "HS" represents the random halfspace depth. By default, depth.dir = "RP".
#' @return a list including the index of the mean outlyingness and variational outlyingness
#' @export
#' @examples
#' library(spam);n = 200;p = 4;Sigmainv <- .25^abs(outer(1:p,1:p,"-"))
#' Sigmainv <- as.spam( Sigmainv, eps=1e-4)
#'Sigma <- solve( Sigmainv);data <- rmvnorm(n, mu = rep(0, p), Sigma = Sigma)
#' outlier_dirout(data, sq_vo = FALSE)


outlier_dirout <- function(data, sq_vo = FALSE, depth.dir = "RP") {
  ##source("directional_outlying.R")
  ##source("faccal_num.R")
  result <- directional_outlying(data, diroutmatrix = TRUE,
                                 depth.dir = depth.dir, d_value = TRUE, sq_vo)
  mo <- result$out_avr
  vo <- result$out_var
  mah_dist <- result$d_value
  if ("list" %in% class(data)) {
    d <- length(data) ## number of variables
    n <- dim(data[[1]])[1] ## number of observations
  } else if ("matrix" %in% class(data)){
    d <- 1 ## number of variables
    n <- nrow(data) ## number of observations
  } else {
    stop("data needs to be either a matrix or a list!!!")
  }
  ## n is the observation number
  ## dim is the dimension of the t(MO, vo).
  fact <- faccal_num(n, d + 1) ## result includes S (covariance of the (MO, Vo) (p + 1) * (p + 1) )
  fac <- fact$fac1
  cutoff <- fact$fac2
  cutoff <- cutoff / fac #cut off value for testing/outlier detection#
  out.dir <- which(mah_dist > cutoff)
  medcurve <- which.min(mah_dist)
  return(list(outlier = out.dir, median = medcurve, mo = mo, vo = vo))
}
