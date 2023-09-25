#' Description of subsetlist.R.
#' More information about what subsetlist.R function does.
#'
#' @param myList Input parameter is a list
#' @param elementNames Input parameter is the element index that of interests.
#' @return a subset list
#' @export
#' @examples subsetlist(list(l1 = c(1, 2, 3), l2 = c("d", "efd", "g"), l3 =c(12,345), l4 = 1), c(1, 3))

subsetlist <- function(myList, elementNames) {
  if (length(elementNames) == 0) {
    return (NULL)
  } else if (length(elementNames) > 1 && all(elementNames > 0)) {
    lapply(elementNames, FUN = function(x) { myList[[x]] })
  } else if (length(elementNames) == 1 && all(elementNames > 0)) {
    P <- list()
    P[[1]] <- myList[[elementNames]]
    return (P)
  } else {
    stop("The second argument should be a positive integer or a vector of positive integers")
  }
}
