#' Description of two_layer_partition_clustering.R.
#' This function implements the two-layer partition clustering.
#' See Robust two-layer partition clustering of sparse multivariate functional data,
#' Zhuo Qu, Wenlin Dai, and Marc G. Genton, https://www.sciencedirect.com/science/article/pii/S2452306223000138
#' This function implements algorithm 2 in the paper.
#' @param neighbour_matrix is the logical matrix of whether each subject is the neighbour or not.
#' @return A list of two elements:
#' 1. group is the first-layer partition clustering. group is a list, where each sublist is a subgroup including center and element.
#' 2. clust is the second-layer partition clustering, where each sublist is a subclust including center and element.
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
#' two_layer_partition_clustering(neighbors_canada.temp)


two_layer_partition_clustering <- function(neighbour_matrix) {
  ##source("subsetlist.R")
  ################### group partition ########################
  M <- 0
  S <- 1:nrow(neighbour_matrix)
  group <- list()
  while (length(S) > 0) {
    M <- M + 1
    if (length(S) == 1) {
      element <- S
      group[[M]] <- list(center = S, element = S)
    } else {
      center <- S[which(apply(neighbour_matrix[S, S], 1, sum)
                      == max(apply(neighbour_matrix[S, S], 1, sum)))]
      if (length(center) > 1) {
        center <- center[1]
      }
      element <- intersect(which(neighbour_matrix[center, ] == TRUE), S)

      group[[M]] <- list(center = center, element = element)
    }
    S <- setdiff(S, element)
  }
##################### group partition  ###################################

##################### Cluster Partition #################################
  P <- group
  I <- 0
  cluster_iteration <- list()
  while (length(P) > 0) {
    I <- I + 1
    cluster_iteration[[I]] <- P[[1]]$element
    remove_index <- 1
    neighbour_bound <- c()
    for (l in cluster_iteration[[I]]) {
      neighbour_bound <- union(neighbour_bound, which(neighbour_matrix[l, ] == TRUE))
    }
    if (length(P) > 1) {
      for (l in 2:length(P)) {
        for (cen in P[[l]]$center) {
          if (cen %in% neighbour_bound) {
            cluster_iteration[[I]] <- union(cluster_iteration[[I]], P[[l]]$element)
            remove_index <- c(remove_index, l)
            break
          }
        }
      }
    }
    P <- subsetlist(P, setdiff(1:length(P), remove_index))
  }
##################### Cluster Partition #############################
  clust <- vector(mode = "list", length = length(cluster_iteration))
  clust <- lapply(1:length(cluster_iteration),
                  function(l) {
                    element = cluster_iteration[[l]]
                    if (length(element) > 1) {
                      neighbour_within_clust <- apply(neighbour_matrix[element, element], 1, sum)
                      center <- element[which(neighbour_within_clust == max(neighbour_within_clust))]
                    } else {
                      center = element
                    }
                    return (list(center = center, element = element))
                  })
  return (list(group, clust))
}
