# functions to transform covariates

#' @name pcaTrans
#' @rdname pcaTrans
#'
#' @title PCA transform a \code{Raster*} object
#'
#' @description Given a matrix of point coordinates and a Raster* object of
#'  covariates, return a \code{Raster*} object of the same size with the layers
#'  giving principal components from a PCA rotation based on the data values
#'  at the coordinates.
#'
#' @param coords a two-column matrix or dataframe giving the location about
#' which to carry out the PCA analysis
#'
#' @param covs a \code{Raster*} object containing covariates to rotate
#'
#' @return a \code{Raster*}
#'
#' @family transform
#'
#' @export
#' @import raster
#'
pcaTrans <- function(coords, covs) {


  # get covariate data at point locations
  vals <- extract(covs, coords)
  vals <- na.omit(vals)

  # do pca analysis
  pca <- prcomp(vals, retx = FALSE, center = TRUE, scale. = TRUE)

  # find non-missing cells
  cell_idx <- getCellIdx(covs)

  # extract covariate values
  vals <- raster::extract(covs, cell_idx)

  # convert to a data.frame
  vals <- data.frame(vals)

  # get PCA predictions
  vals_trans <- predict(pca, vals)

  # set new raster values
  trans_ras <- insertRaster(covs = covs,
                            new_vals = vals_trans,
                            idx = cell_idx)
  return (trans_ras)

}
