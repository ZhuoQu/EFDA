% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/two_layer_partition_clustering.R
\name{two_layer_partition_clustering}
\alias{two_layer_partition_clustering}
\title{
This function implements the two-layer partition clustering step in the robust two-layer partition clustering method.}
\usage{
two_layer_partition_clustering(neighbour_matrix)
}
\arguments{
\item{neighbour_matrix}{is the logical matrix of whether each subject is the neighbour or not.}
}
\value{
A list of two elements:
1. group is the first-layer partition clustering. group is a list, where each sublist is a subgroup including center and element.
2. clust is the second-layer partition clustering, where each sublist is a subclust including center and element.
}
\description{
Description of two_layer_partition_clustering.R.
This function implements the two-layer partition clustering step.
See Robust two-layer partition clustering of sparse multivariate functional data,
Zhuo Qu, Wenlin Dai, and Marc G. Genton, https://www.sciencedirect.com/science/article/pii/S2452306223000138
This function implements algorithm 2 in the paper.
}
\examples{
library(fda); canada_temp <- t(CanadianWeather$monthlyTemp) ## a matrix of 35 stations and 12 months
n <- nrow(canada_temp); nt <- ncol(canada_temp); time <- 1:12
complete_canada.temp_list <- lapply(1:n, function(k) {
res <- list(argvals = time,
            subj = rep(k, 12), y = canada_temp[k, ])
 return(res)
 })
 dist_complete_canada.temp <- elastic_time_distance(complete_canada.temp_list)
 neighbors_canada.temp <- findneighbours(dist_complete_canada.temp)
no.of.neighbors <- apply(neighbors_canada.temp, 1, sum)
two_layer_partition_clustering(neighbors_canada.temp)
}
