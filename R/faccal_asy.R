#' Description of faccal_asy.R.
#' to estimate the the factor and the cutoff value of the approximate F distribution of factor times the distance matrix
#'  between any observation and minimum determinant covariance mean.
#' See the Appendix in The Distribution of Robust Distances (2005), Johanna Hardin & David M Rocke,
#' Journal of Computational and Graphical Statistics,14(4), Pages 928--946
#'  https://www.tandfonline.com/doi/epdf/10.1198/106186005X77685?needAccess=true
#' @param n is the number of observations in the multivariate data
#' @param dim is the dimension of the data
#' @return a list showing the factor and the cutoff value
#' @export
#' @examples faccal_asy(100, 3)

faccal_asy <- function(n, dim) {
  h <- floor((n + dim + 1) / 2) ##  the minimum number of points which must not be outlying
  alpha <- (n -h) /n ## (A.1) in the paper
  q_alpha <- qchisq(p = 1 - alpha, df = dim) ## (A.2) in the paper
  c_alpha <- (1 - alpha) / pchisq(q = q_alpha, df = dim + 2) ## (A.3) in the paper
  c2 <- -pchisq(q = q_alpha, df = dim + 2) / 2 ## (A.4) in the paper
  c3 <- -pchisq(q = q_alpha, df = dim + 4) / 2 ## (A.5) in the paper
  c4 <- 3 * c3 ## (A.6) in the paper
  b1 <- c_alpha * (c3 - c4) / (1 - alpha) ## (A.7) in the paper
  b2 <- 0.5 + c_alpha / (1 - alpha) * (c3 - q_alpha * (c2 + (1 - alpha) / 2) / dim) ## (A.8) in the paper
  v1 <- (1 - alpha) * b1^2 * (alpha * (c_alpha * q_alpha / dim - 1)^2 - 1) -
    2 * c3 * c_alpha^2 * (3 * (b1 - dim * b2)^2 +
                            (dim + 2) * b2 * (2 * b1 - dim * b2)) ## (A.9) in the paper
  v2 <- n * (b1 * (b1 - dim * b2) * (1 - alpha))^2 * c_alpha^2 ## (A.10) in the paper
  v <- v1 / v2  ## (A.11) in the paper
  m_asy <- 2 / (c_alpha^2 * v) ## (A.12) in the paper

  m <- m_asy * exp(0.725 - 0.00663 * dim - 0.078 * log(n)) ## from Equation (3.4) in the paper


  c <- pchisq(q = q_alpha, df = dim + 2) / (1 - alpha) ## Equation before (3.4) in the paper

  fac <- c * (m - dim + 1) / (dim * m) ## c(m - dim + 1)/(dim * m) times d^2_{S*}(X_i, bar(X)^*) ~ F_{dim, m - dim + 1}
  cutoff <- qf(p = 0.993, df1 = dim, df2 = m - dim + 1)
  return (list(fac1 = fac, fac2 = cutoff))
}
