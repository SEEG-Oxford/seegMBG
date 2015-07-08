# GIS manipulation functions

#' @name bufferMask
#' @rdname bufferMask
#'
#' @title Create a mask by distance from points
#'
#' @description Create a mask \code{RasterLayer} with NA values in all cells more
#'  than a specified distance from a set of points.
#'
#' @param mask a \code{Raster*} object, used as the starting mask. All cells which
#'  are NA in \code{mask} will be NA in the resulting \code{RasterLayer}.
#'
#' @param points a set of points (given as a two-column matrix of coordinates)
#'  around which top calculate the buffer.
#'
#' @param buffer the distance from points beyond which to set cells to NA.
#'  If the coordinate system of \code{mask} is latitude/longitude, this will be
#'  in metres, otherwise in the units of the \code{mask}'s coordinate system.
#'
#' @return A \code{RasterLayer} with the same extent, resolution and coordinate system
#' as \code{mask} and with non-NA cells being 0.
#'
#' @family GIS
#'
#' @export
#' @import raster
#'
bufferMask <- function (mask, points, buffer = 0) {

  # set mask to 0
  mask <- mask * 0

  # get cell numbers for all points
  points_cells <- raster::extract(mask,
                          points,
                          cellnumbers = TRUE)[, 1]

  # find spatially unique records (on grid)
  unique_idx <- !duplicated(points_cells)

  # find values (and numbers) of cells falling under buffered, spatially
  # unique points
  mask_vals <- raster::extract(mask,
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
#'  rows of \code{new_vals}. If not provided, it will be calculated from
#'  \code{raster} using \code{\link{cellIdx}}.
#'
#' @return a RasterBrick object with the same extent and resolution as
#'  \code{raster} and as many layers as the number of columns in
#'  \code{new_vals}.
#'
#' @family GIS
#'
#' @export
#' @import raster
#'
insertRaster <- function (raster, new_vals, idx = NULL) {


  # calculate cell index if not provided
  if (is.null(idx)) idx <- cellIdx(raster)

  # check the index makes superficial sense
  stopifnot(length(idx) == nrow(new_vals))
  stopifnot(max(idx) <= ncell(raster))

  # create results raster
  n <- ncol(new_vals)
  raster_new <- raster::brick(replicate(n,
                              raster[[1]],
                              simplify = FALSE))
  names(raster_new) <- colnames(new_vals)

  # update the values
  for(i in 1:n) {
    raster_new[[i]][idx] <- new_vals[, i]
  }

  return (raster_new)

}


#' @name ll2cart
#' @rdname ll2cart
#'
#' @title Convert from Latitude/Longitude to Cartesian Coordinates
#'
#' @description Convert coordinates given as latitude and longitude to
#'  three-dimensional coordinates on a sphere. This can be used to account
#'  for great-circle distances when fitting spatial models.
#'
#' @param longlat a matrix or dataframe giving the longitudes (first column)
#'  and latitudes (second column) of a set of locations
#' @param radius radius of the sphere in on which the Cartesian coodinates are
#'  defined. By default the approximate radius of the earth in kilometres.
#'
#' @export
#'
#' @return ll2cart: a three-column dataframe giving the x, y and z positions
#'  for each set of coordinates
#'
ll2cart <- function (longlat, radius = 6371) {

  # check inputs
  stopifnot(is.finite(radius) & radius > 0)
  stopifnot(inherits(longlat, 'matrix') | inherits(longlat, 'data.frame'))
  stopifnot(ncol(longlat) == 2)

  # extract required columns
  longitude <- longlat[, 1]
  latitude <- longlat[, 2]

  # convert
  ans <- data.frame(x = radius * cos(latitude) * cos(longitude),
                    y = radius * cos(latitude) * sin(longitude),
                    z = radius * sin(latitude))

  return (ans)

}


#' @name cart2ll
#' @rdname ll2cart
#'
#' @param xyz a matrix or dataframe giving the cartesian coordinates
#'  (in order x, y, z) of a set of locations
#'
#' @export
#'
#' @return cart2ll: a two-column dataframe giving the longitudes and latitudes
#'  for each set of coordinates
#'
cart2ll <- function (xyz, radius = 6371) {

  # check inputs
  stopifnot(is.finite(radius) & radius > 0)
  stopifnot(inherits(longlat, 'matrix') | inherits(longlat, 'data.frame'))
  stopifnot(ncol(xyz) == 3)

  # extract columns
  x <- xyz[, 1]
  y <- xyz[, 2]
  z <- xyz[, 3]

  # convert
  ans <- data.frame(longitude = atan2(y, x),
                    latitude = asin(z / radius))

  return (ans)

}
