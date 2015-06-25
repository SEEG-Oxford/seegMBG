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
#' @export
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
#' @title Find non-NA \code{Raster} or \code{RasterBrick} cell numbers
#'
#' @description Get the cell numbers matching all non-NA cells.
#'
#' @param x a \code{Raster} or \code{RasterBrick} object
#'
#' @return a numeric vector of cell numbers
#'
#' @family utility
#'
#' @export
#' @import raster
#'
cellIdx <- function (x) which(!is.na(getValues(x[[1]])))

# get formula left-hand-side
getLHS <- function(f) {

  # check input
  stopifnot(inherits(f, 'formula'))

  # get element
  if(length(f) == 2) {
    stop ('this formula has no left-hand side')
  } else {
    ans <- f[[2]]
  }

  # convert to string
  ans <- deparse(ans)

  return (ans)
}

# get formula right-hand-side
getRHS <- function(f, split = TRUE) {

  # check input
  stopifnot(inherits(f, 'formula'))

  # get element
  if(length(f) == 2) {
    ans <- f[[2]]
  } else if (length(f) == 3) {
    ans <- f[[3]]
  } else {
    stop ('I have no idea what this is')
  }

  # convert to string
  ans <- deparse(ans)


  # optionally split up predictors
  if (split) {
    ans <- strsplit(ans, '\\+')[[1]]
  }

  return (ans)

}

#' @name merge.formula
#' @rdname merge.formula
#'
#' @title Merge Two Formulae
#'
#' @description Merge two formulae to combine predictors (right hand side).
#'
#' @param x a \code{formula} object with both response and predictor terms
#'
#' @param y a \code{formula} object with predictor terms (response ignored)
#'
#' @return a \code{formula} with the response form \code{x} and the predictors
#'  from both \code{x} and code{y}
#'
#' @family utility
#'
#' @export
#'
#' @examples
#' f <- y ~ x
#' g <- ~ a * b + c
#'
#' fg <- merge(f, g)
#' fg
#'
merge.formula <- function(x, y, ...){

  # modified from:
  # https://stevencarlislewalker.wordpress.com/2012/08/06/merging-combining-adding-together-two-formula-objects-in-r/

  # get character strings of the names for the responses
  # (i.e. left hand sides, lhs)
  lhs <- getLHS(x)

  # get character strings of the right hand sides
  rhs <- c(getRHS(x), getRHS(y))

  # put the two sides together
  ans <- reformulate(rhs, lhs)

  # set the environment of the result
  environment(ans) <- parent.frame()

  return (ans)
}

#' @rdname merge.formula
#'
#' @export
#' @examples
#'
#' \\
#'
#' fg <- f + g
#' fg
#'
`+.formula` <- function (x, y) {
  ans <- merge(e1, e2)
  environment(ans) <- parent.frame()
  return (ans)
}
