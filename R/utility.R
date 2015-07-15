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
  ans <- paste0(deparse(ans), collapse = '')

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
  ans <- paste0(deparse(ans), collapse = '')

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
#' @param \dots further arguments passed to or from other methods.
#'
#' @return a \code{formula} with the response form \code{x} and the predictors
#'  from both \code{x} and \code{y}
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
#'
#' fg <- f + g
#' fg
#'
`+.formula` <- function (x, y) {
  ans <- merge(x, y)
  environment(ans) <- parent.frame()
  return (ans)
}

#' @name expand
#' @rdname expand
#'
#' @title Expand a Vector into a Multi-column Matrix
#' @description Create a matrix by combining multiple copies of a vector as
#'  columns
#'
#' @param x a vector (data to repeat into \code{n} columns)
#' @param n an integer (number of times to repeat \code{x})
#'
#' @export
#'
#' @return a matrix with dimensions \code{length(x), n}, with \code{x} repeated
#'  in each column
#'
#' @family utility
#'
#' @examples
#'
#' x <- rnorm(5)
#' expand(x, 3)
#'
expand <- function(x, n) {

  # check input
  stopifnot(length(n) == 1)
  stopifnot(is.vector(x))

  # replicate x as an n-column matrix
  ans <- matrix(rep(x, n), nrow = length(x))

  return (ans)

}

#' @name rowProds
#' @rdname rowProds
#'
#' @title Form Row Products
#' @description Form row products for numeric arrays
#'
#' @param x an array of two or more dimensions, containing numeric, complex,
#'  integer or logical values, or a numeric data frame.
#' @param na.rm logical. Should missing values (including NaN) be omitted from
#'  the calculations?
#' @param dims integer: Which dimensions are regarded as rows to
#'  sum over. The product is over dimensions dims+1.
#'
#' @export
#'
#' @return A numeric or complex array of suitable size, or a vector if the
#'  result is one-dimensional.
#'
#' @family utility
#'
#' @examples
#'
#' x <- expand(rnorm(5), 3)
#' rowProds(x)
#'
rowProds <- function(x, na.rm = FALSE, dims = 1) {
  apply(x, dims, prod, na.rm = na.rm)
}

#' @name aggMatrix
#' @rdname aggMatrix
#'
#' @title Aggregate Values Accross Rows of a Matrix
#' @description Apply a function columnwise to groups of rows in a matrix,
#'  according to a grouping index. Essentially applies the same
#'  \code{\link{tapply}} to multiple colums of a matrix.
#'
#' @param x a matrix containing numeric, complex, integer or logical values,
#'  or a  data frame.
#' @param index list of one or more factors, each of same length as
#'  \code{dim(x)[1]} (as in \code{tapply}). The elements are coerced
#'  to factors by as.factor.
#' @param fun the function to be applied. In the case of functions
#'  like \code{+}, \code{\%*\%}, etc., the function name must be backquoted
#'  or quoted.
#' @param margin over which margin should the operation take place
#'  (i.e. which margin matches \code{index}).
#' @param \dots other arguments to be passed to \code{fun}.
#'
#' @export
#'
#' @return a matrix with dimensions \code{dim(x)[1], length(index)}
#'
#' @family utility
#'
#' @examples
#'
#' x <- expand(rnorm(5), 3)
#' aggMatrix(x, c(1, 1, 1, 2, 2))
#'
aggMatrix <- function(x, index, fun = sum, margin = 1, ...) {
  # aggregate rows of a matrix for dataframe using function and according to
  # index
  agg <- function(x) tapply(x, INDEX = index, FUN = fun, ...)
  margin <- ifelse(margin == 1, 2, 1)
  apply(x, margin, agg)
}


