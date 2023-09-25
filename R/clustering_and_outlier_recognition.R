#' Description of two_layer_partition_clustering.R.
#' This function implements the two-layer partition clustering.
#' See Robust two-layer partition clustering of sparse multivariate functional data,
#' Zhuo Qu, Wenlin Dai, and Marc G. Genton, https://www.sciencedirect.com/science/article/pii/S2452306223000138
#' This function implements algorithm 3 in the paper.
#' @param cluster_result is the list of the clustering function. Each list includes center, and element
#' @param minprop minimal required number of elements in a cluster / the number of samples in the population.
#' By default, minprop = 0.05.
#' @param distance_matrix is the matrix of distances for the population.
#' @param alpha_c is used to defines the alpha_c quantile of the distances. By default, alpha_c = 0.85.
#' between curves inside a cluster and the cluster center.
#' @return A list of three elements:
#' 1. clust is the result of the robust two-layer partition clustering, where each sublist is a subclust including center and element.
#' 2. isolate is the result of outliers in the robust two-layer partition clustering.
#' 3. criteria_matrix is a matrix of the distance between each potential outlier and the center of all clust. The cut-off values are the distance
#' threshold determined by alpha_c.
#' @export
#' @examples
#' library(fda); canada_temp <- t(CanadianWeather$monthlyTemp) ## a matrix of 35 stations and 12 months
#'  n <- nrow(canada_temp); nt <- ncol(canada_temp); time <- 1:12
#'  complete_canada.temp_list <- lapply(1:n, function(k) {
#'  res <- list(argvals = time,
#'              subj = rep(k, 12), y = canada_temp[k, ])
#'   return(res)
#'   })
#'   dist_complete_canada.temp <- elastic_time_distance(complete_canada.temp_list)
#'   neighbors_canada.temp <- findneighbours(dist_complete_canada.temp)
#'  no.of.neighbors <- apply(neighbors_canada.temp, 1, sum)
#'  two_layer_partition <- two_layer_partition_clustering(neighbors_canada.temp)
#'  cluster_canada.temp <- two_layer_partition[[2]]
#'  clust_outlier_canada.temp <- clustering_and_outlier_recognition(cluster_canada.temp, minprop = 0.1, distance_matrix = dist_complete_canada.temp)

clustering_and_outlier_recognition <- function(cluster_result,
                                               minprop = 0.05,
                                               distance_matrix,
                                               alpha_c = 0.85) {
  ##source("subsetlist.R")
  ###################################################
  basis_cluster_index <- which(lapply(cluster_result, function(k) {
    length(k$element) > minprop * nrow(distance_matrix)}) == T)

  isolatecenter <- subsetlist(cluster_result, setdiff(1:length(cluster_result), basis_cluster_index))

  maincenter <- subsetlist(cluster_result, basis_cluster_index)

  if (length(isolatecenter) == 0) {
    return (list(clust = maincenter))
  } else if (length(maincenter) == 0) {
    return (list(isolate = isolatecenter))
  } else {
    Distance_to_gravity <- lapply(maincenter, function(l) {
    distance_matrix[l$center[1], l$element] })

    distance_threshold <- sapply(Distance_to_gravity, function(l) {
      sort(l)[ceiling(alpha_c * length(l))]
    })
    isolate <- c()
    ###############################################
    if (length(maincenter) > 1) {
      D_g <- t(sapply(isolatecenter, function(nr) {
        sapply(maincenter, function(nc) {
          return (distance_matrix[nr$center[1], nc$center[1]])
        })
      }))
      D_g_bool <- t(apply(D_g, 1, function(k) {
        temp <- c()
        for (l in 1:length(k)) {
          temp <- c(temp, k[l] <= distance_threshold[l])
        }
        return (temp)
      }))
    } else {
      D_g <- sapply(isolatecenter, function(nr) {
      return (distance_matrix[nr$center[1], maincenter[[1]]$center[1]])
     })
      D_g <- as.matrix(D_g)
      D_g_bool <- D_g <= distance_threshold
    }
    rownames(D_g) <- lapply(isolatecenter, function(k) {
      return (k$center[1])
    })
    for (l in 1:nrow(D_g)) {
      if (length(unique(D_g_bool[l, ])) == 1 &&
          unique(D_g_bool[l, ]) == FALSE) {
        isolate <- c(isolate, isolatecenter[[l]]$element)
      } else {
        cand_index <- which(D_g_bool[l, ] == TRUE)
        main_quantile <- sapply(cand_index, function(num) {
          return (rank(c(D_g[l, num], Distance_to_gravity[[num]])[1] /
                         (length(Distance_to_gravity[[num]]) + 1) ) )
        })
        index <- which(main_quantile == min(main_quantile))
        if (length(index) > 1) {
          index <- which(D_g[l, ] == min(D_g[l, ]))
        }
        maincenter[[index]]$element <- c(maincenter[[index]]$element,
                                         isolatecenter[[index]]$element)
      }
    }
    s <- rbind(D_g, distance_threshold)
    return (list(clust = maincenter,
                 isolate = sort(isolate),
                 criteria_matrix = s))
  }
}



#############################################
# plot(Distance_to_gravity[[1]], rank(Distance_to_gravity[[1]]) / length(Distance_to_gravity[[1]]),
#      xlim = range(c(unlist(Distance_to_gravity), D_g)), ylim = c(0, 1), type = "n", xlab ="Distance to the Center",
#      ylab ="Empirical cdf")
# if (length(maincenter) > 1) {
#   lapply(1:length(maincenter), function(l) {
#   candidate_curve <- c(Distance_to_gravity[[l]], D_g[, l])
#   points(candidate_curve, rank(candidate_curve) / length(1 + candidate_curve), col = l)
#   #points(D_g[,l],rank(D_g[,l])/length(1+candidate_curve),col=l, pch=4)
# })
# } else {
#   candidate_curve <- c(Distance_to_gravity[[1]], D_g)
#   points(candidate_curve, rank(candidate_curve) / length(1 + candidate_curve))
#
# }
#############################################
