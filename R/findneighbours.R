#' Description of findneighbours.R.
#' This function returns the logical matrix. Each row shows whether each subject is the neighbor.
#' @param distance_matrix is the matrix of distances for the population.
#' @param theta is the probability (0 < theta < 1) for the cut-off quantile. By default, theta = 0.05.
#' @return a logical matrix
#' @export
#' @examples
#'library(fda); canada_temp <- t(CanadianWeather$monthlyTemp) ## a matrix of 35 stations and 12 months
#' n <- nrow(canada_temp); nt <- ncol(canada_temp); time <- 1:12
#' complete_canada.temp_list <- lapply(1:n, function(k) {
#' res <- list(argvals = time,
#'             subj = rep(k, 12), y = canada_temp[k, ])
#'  return(res)
#'  })
#'  dist_complete_canada.temp <- elastic_time_distance(complete_canada.temp_list)
#'  neighbors_canada.temp <- findneighbours(dist_complete_canada.temp)
#' no.of.neighbors <- apply(neighbors_canada.temp, 1, sum)

findneighbours <- function(distance_matrix, theta = 0.05) {
  numer <- distance_matrix[lower.tri(distance_matrix, diag = FALSE)]
  cutoffvalue <- quantile(numer[numer != 0], probs = theta)
  neighbour_matrix <- distance_matrix < cutoffvalue
  return (neighbour_matrix)
}
