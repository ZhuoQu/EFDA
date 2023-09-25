#' Description of outlying_depth.R.
#' Obtain the directional outlyingness for multivariate data
#'
#' @param data_mat  a matrix with no. rows = number of observations, and no. columns = number of variables.
#' @param depth.dir chracter variable in "RP" (Random Projection depth), "MhD" (Mahalanobis depth),
#'  "SD" (Simplicial depth), "HS" (Random Halfspace depth).
#' @return a vector with length = number of observations. It measures the outlyingness of the data which is an inverse of the depth.
#' @export
#' @examples
#' library(spam); n = 200; p <- 4; Sigmainv <- .25^abs(outer(1:p,1:p,"-"))
#' Sigmainv <- as.spam( Sigmainv, eps=1e-4)
#' Sigma <- solve(Sigmainv)
#' data <- rmvnorm(n, mean = rep(0, p), Sigma = Sigma)
#' outlying_depth(data, "RP")
outlying_depth <- function(data_mat, depth.dir) {
  if (depth.dir == "RP") {
    out <- 1 / mdepth.RP(data_mat, proj = 200)$dep - 1
  } else if (depth.dir == "MhD") {
    out <- 1 / mdepth.MhD(data_mat)$dep - 1
  } else if (depth.dir == "SD") {
    out <- 1 / mdepth.SD(data_mat)$dep - 1
  } else if (depth.dir == "HS") {
    out <- 1 / mdepth.HS(data_mat)$dep - 1
  }
  return (out)
}
