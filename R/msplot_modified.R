#' Description of msplot_modified.R.
#' Obtain the Magnitude-Shape plot
#'
#' @param data is a list with d (the number of variables). And each list is a n*time_nb matrix, where n is the number of
#' observations, and time_nb is the number of the time grid.
#' @param depth.dir is a chracter variable in c("RP","MhD","SD","HS"), where "RP" represents the random
#' projection depth, "MhD" represents the Mahalanobis depth, "SD" represents the simplicial depth, and
#' "HS" represents the random halfspace depth. By default, depth.dir = "RP".
#' @param plot is a logical variable of whether to show a figure.
#' @param col.out is the color of outliers. By default, col.out = "green".
#' @param col.med is the color of the median. By default, col.med = "purple".
#' @param col.normal is the color of the non-outlyinge curves. By default, col.normal = "black".
#' @param plot.type is either "parallel" or "scatter". By default, plot.type = "scatter".
#' @param sq_vo is a logical variable of whether to use the variational outlyingness or the
#' square root of the variational outlyingness. By default, sq_vo is FALSE.
#' @param datatopic is the description of your data. This will be shown in the title of the figure.
#' @return a list of 4 elements:
#' 1. mo. If data is a list, then mo is a matrix of n * p, where n is the number of samples,
#' and p is the number of variables.
#' 2. vo. vo is a vector of n, and n is the number of samples in the population.
#' 3. the index of outliers from the directional outlyingness
#' 4. the index of the median curve from the directional outlyingness.
#' @export
#' @examples
#' library(fda)
#' library(roahd)
#' canada_temp <- t(CanadianWeather$monthlyTemp)
#' canada_precip <- t(CanadianWeather$monthlyPrecip)
## 1. Univariate Functional Data
#' msplot_modified(data = canada_temp,
#'              plot = TRUE, sq_vo = FALSE)
### 2. Bivariate Functional Data
#' msplot_modified(data = list(temp = canada_temp, precip = canada_precip),
#'              plot = TRUE, sq_vo = FALSE)
#'
msplot_modified <- function(data, depth.dir = "RP", plot = FALSE,
                            col.out = "green", col.med = "purple",
                            col.normal = "black", plot.type = "scatter",
                            sq_vo = FALSE, datatopic = NULL) {
  ###pairwise plots of variation of outlyingness (VO) against mean outlyingness (MO)###
  ##source("outlier_dirout.R")
  if ("list"%in% class(data)) {
    n <- nrow(data[[1]])
    time_nb <- ncol(data[[1]])
    d <- length(data)
  } else {
    n <- nrow(data)
    time_nb <- ncol(data)
    d <- 1
    }
#####################
  outl_detect <- outlier_dirout(data, sq_vo = sq_vo, depth.dir = depth.dir)
  out.dir <- outl_detect$outlier
  medcurve <- outl_detect$median
  mo <- outl_detect$mo
  vo <- outl_detect$vo

  if (plot) {
    M <- cbind(mo, vo)
    col.point <- rep(col.normal, n)
    col.point[medcurve] <- col.med
    col.point[out.dir] <- col.out
    pch <- rep(20, n)
    pch[out.dir] <- 19
    pch[medcurve] <- 17
    outcol <- rep("white", n)
    outcol[out.dir] <- 1:length(out.dir)
    outcol[medcurve] <- "black"

    if (plot.type == "parallel") {
      data.ms <- data.frame(MO = mo, VO = vo, col = col.point)
      paral <- ggparcoord(data.ms, columns = 1:ncol(data.ms), groupColumn = (d + 2), showPoints = TRUE) ### from library GGally
      return(r1 = list(mo = mo, vo = vo, out.dir = out.dir, medcurve = medcurve, paral))
      }

    if (plot.type == "scatter") {
      if (d > 1) {
        MO <- (apply(mo ^ 2, 1, sum))^(1 / 2)
      } else {
        MO <- abs(mo)
        }
      outlabel <- rep('', n)
      outlabel[out.dir] <- out.dir
      outlabel[medcurve] <- medcurve
      ms.data <- data.frame(x = MO, y = vo, col = col.point, outlabel = outlabel)
      plot(MO, vo, type = "n", ylab = "VO", xlab = "||MO||",
          main = paste("MS plot: ", datatopic))
      for (i in setdiff(1:length(outlabel), c(out.dir,medcurve))) {
        points(MO[i], vo[i], pch = pch[i],
               col = col.normal, cex = 0.71)}
      for (i in out.dir) {
        points(MO[i], vo[i], pch = pch[i],
               col = col.out, cex = 0.7)
        }
      points(MO[medcurve], vo[medcurve], col = col.med, pch = pch[medcurve])
      text(MO + (range(MO)[2] - range(MO)[1]) / 70,
      vo + runif(length(vo), (range(vo)[2] - range(vo)[1]) / 45,
      (range(vo)[2] - range(vo)[1]) / 45), labels = outlabel,
      cex = 0.7)
      }
  }
  return (list(mo = mo, vo = vo, out.dir = out.dir, medcurve = medcurve))
}






