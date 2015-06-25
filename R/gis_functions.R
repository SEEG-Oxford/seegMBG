# GIS manipulation functions

#' @name bufferMask
#' @rdname bufferMask
#'
#' @title Create a mask by distance from points
#'
#' @description Create a mask RasterLayer with NA vaules in all cells more
#'  than a specified distance from a set of points.
#'
#' @param mask a Raster* object, used as the starting mask. All cells which
#'  are NA in \code{mask} will be NA in the resulting RasterLayer.
#'
#' @param points a set of points (given as a two-column matrix of coordinates)
#'  around which top calculate the buffer.
#'
#' @param buffer the distance from points beyond which to set cells to NA.
#'  If the coordinate system of \code{mask} is latitude/longitude, this will be
#'  in metres, otherwise in the units of the \code{mask}'s coordinate system.
#'
#' @return A RasterLayer with the same extent, resolution and coordinate system
#' as \code{mask} and with non-NA cells being 0.
#'
#' @family GIS
#'
bufferMask <- function (mask, points, buffer = 0) {

  # set mask to 0
  mask <- mask * 0

  # get cell numbers for all points
  points_cells <- extract(mask,
                          points,
                          cellnumbers = TRUE)[, 1]

  # find spatially unique records (on grid)
  unique_idx <- !duplicated(points_cells)

  # find cvalues (and numbers) of cells falling under buffered, spatially
  # unique points
  mask_vals <- extract(mask,
                       points[unique_idx, ],
                       cellnumbers = TRUE,
                       buffer = buffer)

  # get as a vector of unique cells falling under the mask
  mask_cells <- unique(unlist(lapply(mask_vals, '[', , 1)))

  # set these cells to NA
  mask[mask_cells] <- NA

  # and return the new mask
  return (mask)

}

#' @name insertRaster
#' @rdname insertRaster
#'
#' @title Insert new values into a Raster* template
#'
#' @description Given a template raster, new set of values for non-missing
#'  cells and optional index for non-missing cells, create a Raster* object
#'  containing the new values. The raster object will be replicate to have
#'  as many layers as the new dataframe has columns.
#'
#' @param raster a Raster* object with non-NA cells in which to insert new_vals
#' @param new_vals a dataframe of new values to put in the non-NA cells of
#'  \code{raster}, with values from each row in a separate layer
#' @param idx an optional index giving the cell numbers corresponding to the
#'  rows of \code{new_vals}.
#'
#' @return a RasterBrick object with the same extent and resolution as
#'  \code{raster} and as many layers as the number of columns in
#'  \code{new_vals}.
#'
#' @family GIS
#'
insertRaster <- function (raster, new_vals, idx = NULL) {


  # calculate cell index if not provided
  if (is.null(idx)) idx <- getCellIdx(raster)

  # check the index makes superficial sense
  stopifnot(length(idx) == nrow(new_vals))
  stopifnot(max(idx) <= ncell(raster))

  # create results raster
  n <- ncol(new_vals)
  raster_new <- brick(replicate(n,
                              raster[[1]],
                              simplify = FALSE))
  names(raster_new) <- colnames(new_vals)

  # update the values
  for(i in 1:n) {
    raster_new[[i]][idx] <- new_vals[, i]
  }

  return (raster_new)

}








