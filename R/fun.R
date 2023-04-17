#' Standardize a variable
#'
#' Standardizes a numeric vector to mean 0 and variance 1. Note that missing
#' values are ignored.
#'
#' @param x A numeric vector.
#' @returns A standardized variable
stdVar <- function (x) {
  (x - mean(x, na.rm=T))/sd(x, na.rm=T)
}
