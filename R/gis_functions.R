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
#' @family GIS
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
#' @family GIS
#' @export
#'
#' @return cart2ll: a two-column dataframe giving the longitudes and latitudes
#'  for each set of coordinates
#'
cart2ll <- function (xyz, radius = 6371) {

  # check inputs
  stopifnot(is.finite(radius) & radius > 0)
  stopifnot(inherits(xyz, 'matrix') | inherits(xyz, 'data.frame'))
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

#' @rdname getArea
#' @name getArea
#'
#' @title Return the Area of an sp Object
#'
#' @description Determine the area of an sp polygon object,
#'  summing across sub-polygons if necesssary
#'
#' @param sp an object of class \code{sp} represnting one or
#'  more polygons
#'
#' @family GIS
#' @export
#'
#' @return a scalar numeric giving the total area
#'
getArea <- function (sp) {
  # get the area of an sp polygon object
  areas <- sapply(slot(sp, "polygons"), slot, "area")
  area <- sum(areas)
  return (area)
}


#' @title safeMask
#' @rdname safeMask
#'
#' @title Mask a RasterLayer by an sp Object Safely
#'
#' @description Mask a \code{RasterLayer} object using a
#'  \code{sp} object representing one or more polygons in a
#'  such that cells falling under small polygons are masked out.
#'  Note this function is slower than \code{raster::mask}
#'  and has different behaviour.
#'
#' @family GIS
#' @export
#' @import raster
#'
#' @param raster a \code{RasterLayer} object to mask.
#' @param sp an \code{sp} object by which to mask \code{raster}.
#'
#' @return a \code{RasterLayer} object, the same as \code{raster} but with
#'  any cell falling under \code{sp} set to \code{NA}
safeMask <- function(raster, sp) {

  # extract by small polygons
  tmp <- raster::extract(raster, sp, cellnumbers = TRUE, small = TRUE)

  # get all cell numbers under polygons
  cells <- na.omit(unlist(lapply(tmp, function(x) x[, 1])))

  # get those not under polygons
  cells_chuck <- (1:ncell(raster))[-cells]

  # set to NA
  raster[cells_chuck] <- NA

  # return the masked raster
  return (raster)

}


#' @title getPoints
#' @rdname getPoints
#'
#' @title Generate Spatial Integration Points for a Shapefile
#' @description Use random sampling and K-means clustering to generate
#'  a set of coordinates and corresponding weights that can be used
#'  to carry out spatial integration over an area represented by an
#'  \code{sp} object. The weights can be purely spatial, or can
#'  determined by the values of a raster representing the process of
#'  interest - for example population density. The raster should be
#'  such that random points can be sampled from it according to cell
#'  values using \code{seegSDM::bgSample}.
#'
#' @param shape an \code{sp} object containing one or more polygons
#'  representing the region for which integration points are required
#' @param raster a \code{RasterLayer} object, optionally containing values
#'  determining weights for integration
#' @param n the number of integration points required. Rounded up if not an
#'  integer.
#' @param perpixel whether \code{n} gives the expected number of points per
#'  valid (non-NA and non-zero) pixel, or else the total number of points
#' @param prob whether to weight the integration points by the values of
#'  \code{raster}. Pixels with value 0 will never be sampled from and
#'  negative pixels will cause an error. If all cells are 0 or missing,
#'  prob will be set to FALSE.
#'
#' @import seegSDM
#'
#' @family GIS
#' @export
#'
#' @return a three-column matrix giving the coordinates and corresponding
#'  weights for the spatial integration points over \code{sp}. Note that
#'  if there are fewer unique points found than \code{n}, only the unique
#'  points will be returned. If there are no non-missing cells, a dataframe
#'  with 0 rows will be returned and a warning issued.
#'
getPoints <- function (shape,
                       raster,
                       n = 10,
                       perpixel = FALSE,
                       prob = FALSE) {

  # if prob is FALSE, a shit-ton of points are sampled
  # uniformly at random within sp; if prob is TRUE they are biased
  # by population within the polygon.
  # The resulting points are then k-means clustered to yield
  # n representative points

  # check raster
  stopifnot(inherits(raster, 'RasterLayer'))

  # resize it
  raster <- raster::crop(raster, shape)
  raster <- safeMask(raster, shape)


  # get cell values to check them
  vals <- getValues(raster)

  # if there are no valid cells, return no integration points
  # as a dataframe with no rows & issue a warning
  if (all(is.na(vals))) {
    warning ('no non-NA cells found in raster for this polygon, no integration points returned')
    ans <- data.frame(x = NA, y = NA, weights = NA)[0, ]
    return (ans)
  }

  # if prob = TRUE and all valid cells are 0, set prob to FALSE
  if (prob) {

    if (all(na.omit(vals) == 0)) {
      prob <- FALSE
    } else {
      # otherwise set them to NAs
      raster[raster == 0] <- NA
    }
  }

  if (perpixel) {
    # get number of valid pixels
    n_valid <- length(seegSDM:::notMissingIdx(raster))
  } else {
    # otherwise get one times this many
    n_valid <- 1
  }

  # correct total number to integer
  n <- ceiling(n * n_valid)

  # get representative points in an sp object by uniform random
  # sampling then kmeans clustering

  # sample points
  x <- seegSDM::bgSample(raster,
                         10000,
                         prob = prob,
                         spatial = FALSE,
                         replace = TRUE)

  # coerce x to be a dataframe
  x <- as.data.frame(x)

  # make sure there aren't more centres than unique datapoints
  n_unique <- nrow(unique(x))
  n <- pmin(n, n_unique)

  # k-means cluster the data
  kmn <- kmeans(x, n)

  # get the cluster centres
  u <- kmn$centers

  # get the weights (proportion of points falling in that area)
  weights <- summary(factor(kmn$cluster)) / length(kmn$cluster)

  # remove any rownames from the inducing points
  rownames(u) <- NULL

  # make sure they have their column names
  colnames(u) <- colnames(x)

  # add the weights on
  u <- cbind(u, weights = weights)

  # and return them
  return (u)

}
