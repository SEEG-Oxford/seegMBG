# Miscellaneous utility functions

#' @name badRows
#' @rdname badRows
#'
#' @title Find rows with NA values in a dataframe or matrix
#'
#' @description For each row, if any element is NA TRUE is returned.
#'
#' @param x a dataframe or matrix
#'
#' @return a logical vector
#'
#' @family utility
#'
#' @examples
#' # fake data matrix with a missing row
#' x <- matrix(rnorm(3 * 4), nrow = 3)
#' x[2, 3] <- NA
#' badRows(x)
#'
badRows <- function (x) apply(x, 1, function(row) any(is.na(row)))

#' @name cellIdx
#' @rdname cellIdx
#'
#' @title Find non-NA Raster or RasterBrick Cell numbers
#'
#' @description Get the cell numbers matching all non-NA cells.
#'
#' @param x a Raster or RasterBrick object
#'
#' @return a numeric vector foi cell numbers
#'
#' @family utility
cellIdx <- function (x) which(!is.na(getValues(x[[1]])))
