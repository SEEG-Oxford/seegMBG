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
