#' Description of sparse_fbplot.R.
#' This function visualizes the functional data with sparse functional boxplot.
#' When the original data is complete, the sparse functional boxplot becomes functional boxplot.
#' When the original data is sparse, the sparse functional boxplot shows the sparseness proportion information
#' in the central region.
#'
#' @param fit depicts the fitted data. fit is a list with p (the number of variables). And each list is a n*t matrix, where n is the number of
#' observations, and t is the number of the time grid.
#' @param sparse depicts the original data. sparse is a list with p (the number of variables). And each list is a n*t matrix, where n is the number of
#' observations, and t is the number of the time grid.
#' @param time_index is the coordinate of time or index. If not given, 1,...length(x) is provided.
#' By default, time_index = NULL.
#' @param depth is the depth value for n observations.
#' @param two_stage is a logical variable of whether to obtain the two-stage sparse functional boxplot or
#' the sparse functional boxplot. By default, two-stage = TRUE.
#' @param sq_vo is a condition of whether to apply the square root procedure for variational outlyingness.
#' By default, sq_vo = FALSE.
#' @param plot is a logical variable of whether to obtain a figure. By default, plot = TRUE.
#' @param xlab is the title of the x-axis.
#' @param ylab is the title of the y-axis.
#' @param title is the title of the figure.
#' @param yrange is the range of y axis.
#' @param cex.main is the magnification to be used for the title.
#' @param cex.lab is the magnification to be used for the label.
#' @param cex.axis is the magnification to be used for the axis.
#' @param medlabel is a logical variable of whether to show the index of the median curve.
#' By default, medlabel = TRUE.
#' @param outlabel is a logical variable of whether to show the index of the median curve.
#' By default, outlabel = TRUE.
#' @param prob shows the central region of some probability. By default, prob = 0.5.
#' @param factor is a threshold factor to detect outliers. By default, factor = 1.5.
#' @param outliercolor.fb is the color for functional outliers. By default, outliercolor.fb = 2.
#' @param barcol is the colour of the mincurve and maxcurve. By default, barcol = 4.
#' @param outliercolor.dir is the color for functional directional outliers.
#' By default, outliercolor.dir = 3.
#' @param fullout is a logical variable of whether to show the full outliers in the plot.
#' @return plot the sparse functional boxplot.
#' The result is a list of four elements.
#' 1. sparse_density_ct is a list of p.
#' Each list is matrix of t rows and 2 columns, where the first column shows the observation proportion,
#' and the second column shows the sparse proportion at this time point.
#' 2. fb_outlier_index is a list of p. In each list, it shows the outlier from functional boxplot if any
#' 3. dir_outlier_index is a vector showing possible outliers from the directional outlyingness
#' 4. med_index is a index showing the index of the median, which is the observation
#' with the highest depth excluding the directional outlier.
##################     Value     #################################

#' @export
#' @examples
## 1. Univariate complete functional data
#' library(fda)
#' library(roahd)
#' canada_temp <- t(CanadianWeather$monthlyTemp)
#'sparse_fbplot(fit = canada_temp,
#'              sparse = canada_temp,
#'              time_index = 1:12, depth = MBD(canada_temp))
## 2. Univariate sparse functional data
#'  library(refund); list("cd4"); argval <- -18:42;
#'  bootstrap_cd4 <- bootstrap_implementation(bootstrap_times = 100, argval = -18:42, true_data = NULL, sparse_data = cd4, confid_alpha = 0.05)
#' fit_cd4 <- bootstrap_cd4$bootstrap_fit[[1]]
#'  depth_cd4 <- MBD(fit_cd4) ## modified band depth
#'  sparse_fbplot(fit_cd4, sparse = cd4, time_index = NULL, depth = depth_cd4)

## 3. Bivariate complete functional data
#'  library(fda); library(mrfDepth)
#'  canada_temp <- t(CanadianWeather$monthlyTemp) ## a matrix of 35 stations and 12 months
#'  canada_precip <- t(CanadianWeather$monthlyPrecip) ## a matrix of 35 stations and 12 months
#'  true_canada <- list(temp = canada_temp, precip = canada_precip);
#'  n <- nrow(canada_temp); nt <- ncol(canada_temp); missing_prob <- 0.3; argval <- 1:12
#'  canada_array <- array(NA, dim = c(nt, n, 2)); canada_array[, , 1] <- t(canada_temp); canada_array[, , 2] <- t(canada_precip);
#'  mfd_canada <- mfd(x = canada_array, z = canada_array, type = "sdepth")
#'  sparse_fbplot(fit = true_canada, sparse = true_canada, time_index = NULL, depth = as.numeric(mfd_canada$MFDdepthZ))

## 4. Bivariate sparse functional data
#'  library(fda);
#'  canada_temp <- t(CanadianWeather$monthlyTemp) ## a matrix of 35 stations and 12 months
#'  canada_precip <- t(CanadianWeather$monthlyPrecip) ## a matrix of 35 stations and 12 months
#'  true_canada <- list(temp = canada_temp, precip = canada_precip);
#'  n <- nrow(canada_temp); nt <- ncol(canada_temp); missing_prob <- 0.3; argval <- 1:12
#'  random_matrix <- matrix(runif(n * nt), nrow = n, ncol = nt)
##  Create a mask based on the missingness probability
#'  missing_mask <- random_matrix <= missing_prob
#'  sparse_canada <- lapply(true_canada, function(k) {
#'   k[missing_mask] <- NA
#'   return(k)
#'   })
##  Replace values in the full matrix with NA where the mask is TRUE
#'   bt_canada <- bootstrap_implementation(bootstrap_times = 100, argval,
#'   true_data = true_canada, sparse_data = sparse_canada, confid_alpha = 0.05)
#'  fit_canada <- bt_canada$bootstrap_fit
#'  fit_canada_array <- array(NA, dim = c(nt, n, 2));
#'  fit_canada_array[, , 1] <- t(fit_canada[[1]]); fit_canada_array[, , 2] <- t(fit_canada[[2]]);
#'  mfd_fit_canada <- mfd(x = fit_canada_array, z = fit_canada_array, type = "sdepth")
#'  sparse_fbplot(fit = fit_canada, sparse = sparse_canada,
#'                time_index = NULL, depth = as.numeric(mfd_fit_canada$MFDdepthZ))




###############################################################################################
sparse_fbplot <- function (fit, sparse, time_index = NULL, depth = NULL,
                         two_stage = TRUE, sq_vo = FALSE, plot = TRUE,
                         xlab = NULL, ylab = NULL, title = NULL,
                         yrange = NULL,
                         cex.main = 0.9, cex.lab = 0.85, cex.axis = 0.85,
                         medlabel = TRUE, outlabel = TRUE,
                         prob = 0.5, factor = 1.5,
                         outliercolor.fb = 2, barcol = 4,
                         outliercolor.dir = 3, fullout = FALSE) {
  #source("directional_outlying.R")
  #source("fbplot_min_max_curve.R")
  #source("fbplot_polygon.R")
  #source("fbplot_outlier.R")

  if ("list"%in% class(fit)) {
    n <- nrow(fit[[1]])
    tp <- ncol(fit[[1]])
    p <- length(fit)
  } else {
    n <- nrow(fit)
    tp <- ncol(fit)
    p <- 1
  }

  if (length(time_index) == 0) {
    time_index <- 1:tp
  }

  if (class(depth) != "numeric" || length(depth) != n) {
    stop("Depth is not given yet!")
  }
  index <- order(depth, decreasing = TRUE)
  ampli_factor <- ifelse(medlabel == TRUE || outlabel == TRUE, 1 / 40,
                         1 / 120)
  xrange <- c(range(time_index)[1] - ampli_factor * (range(time_index)[2] - range(time_index)[1]),
              range(time_index)[2] + ampli_factor * (range(time_index)[2] - range(time_index)[1]))
  if (plot) {
     colnum <- min(p, 3)
     m <- matrix(1:colnum, nrow = 1, ncol = colnum, byrow = TRUE)
     layout(mat = m,
           widths = rep(1 / colnum, colnum))
    par(mai = c(0.7, 0.8, 0.4, 0.1), mar = c(4.5, 4.2, 2, 0.5),
        mgp = c(3, 1.5, 0))
  }
  if (two_stage == TRUE) {
     outl_detect <- outlier_dirout(fit, sq_vo)
     dir_outlier_index <- outl_detect$outlier
     #med_index <- outl_detect$median
  } else {
    #med_index <- which(depth == max(depth))
    dir_outlier_index <- NULL
  }
  normal_index <- setdiff(index, dir_outlier_index)
  med_index <- which(depth == max(depth[normal_index]))

  sparse_ds <- function(p, variable_index) {
    if (p > 1) {
    subfit <- t(fit[[variable_index]])
    subsparse <- t(sparse[[variable_index]])
    } else {
      subfit <- t(fit)
      subsparse <- t(sparse)
    }
    medavg <- matrix(subfit[, med_index], ncol = length(med_index), nrow = tp)
    med_value <- apply(medavg, 1, mean)

    m <- floor(length(normal_index) * prob)
    center <- subfit[, normal_index[1:m]]
    sp_center <- subsparse[, normal_index[1:m]]

    inf <- apply(center, 1, min) ### shown in figure
    sup <- apply(center, 1, max) ### shown in figure
    dist <- factor * (sup - inf)
    upper <- sup + dist
    lower <- inf - dist
    outly <- (subfit <= lower) + (subfit >= upper)
    outcol <- colSums(outly)
    removefb <- which(outcol > 0)
    fb_outlier_index <- setdiff(removefb, dir_outlier_index)
    if (length(fb_outlier_index) == 0) {
      good <- subfit[, normal_index]
    } else {
      good <- subfit[, setdiff(normal_index,fb_outlier_index)]
    }

    maxcurve <- apply(good, 1, max) ### shown in figure
    mincurve <- apply(good, 1, min) ### shown in figure

    sparse_density_center <- t(sapply(1:nrow(sp_center),
                                      function(nnr) {
      resf <- data.frame(obs = length(which(!is.na(sp_center[nnr, ]))) / m,
                       spa = length(which(is.na(sp_center[nnr, ]))) / m)
      }))

    if (plot) {
      fbplot_min_max_curve(time_index, xrange, yrange, inf, sup,
                           subfit, mincurve, maxcurve,
                           barcol = 4,
                           xlab[variable_index],
                          ylab[variable_index],
                          title[variable_index],
                          cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis)
      if (fullout == FALSE) {
        fbplot_outlier(two_stage, outlabel,
                       outliercolor.fb = 2,
                       outliercolor.dir = 3,
                       fb_outlier_index, dir_outlier_index,
                       time_index,
                       subfit, subsparse)
        fbplot_polygon(time_index, xrange,
                       inf, sup,
                       subfit, subsparse,
                       sparse_density_center,
                       med_index, med_value,
                       medlabel, barcol = 4)
      } else {
        fbplot_polygon(time_index, xrange,
                       inf, sup,
                       subfit, subsparse,
                       sparse_density_center,
                       med_index, med_value,
                       medlabel, barcol = 4)
        fbplot_outlier(two_stage, outlabel,
                       outliercolor.fb = 2,
                       outliercolor.dir = 3,
                       fb_outlier_index, dir_outlier_index,
                       time_index,
                       subfit, subsparse)
      }
    }
    return(list(fb_outlier_index = fb_outlier_index,
                sparse_density_ct = sparse_density_center))
  }
  #####################################################################################
  execute<- lapply(1:p, function(k) {sparse_ds(p, k)})
  sparse_density <- lapply(execute, function(k) {k$sparse_density_ct})
  fb_outlier_index <- lapply(execute, function(k) {k$fb_outlier_index})

  return(list(sparse_density_ct = sparse_density,
              fb_outlier_index = fb_outlier_index,
              dir_outlier_index = dir_outlier_index,
              med_index = med_index))
}

