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
#' @return a \code{Raster*} object with the same extent, resolutiona and number
#'  of layers as \code{covs} but with each layer giving the location on a
#'  different principal component axis.
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
  cell_idx <- cellIdx(covs)

  # extract covariate values
  vals <- raster::extract(covs, cell_idx)

  # convert to a data.frame
  vals <- data.frame(vals)

  # get PCA predictions
  vals_trans <- predict(pca, vals)

  # set new raster values
  trans_ras <- insertRaster(raster = covs,
                            new_vals = vals_trans,
                            idx = cell_idx)
  return (trans_ras)

}



#' @name gamTrans
#' @rdname gamTrans
#'
#' @title Carry out a model-based covariate transformation using a GAM
#'
#' @description Define an optimal set of univariate covariate
#'  transformations of a set of model covariates by fitting a generalised
#'  additive model with univariate smoothers to data, and then using the
#'  smoothers to spline-transform the covariates.
#'  This makes use of the\code{type = 'terms'} argument in
#'  \code{\link{predict.gam}}.
#'  This function also makes use of
#'
#' @param coords a two-column matrix of coordinates of records
#'
#' @param response an object acting as thge response object in the GAM
#'  model (e.g. a vector of counts, or a matrix for binomial data)
#'
#' @param covs a \code{Raster*} object giving the spatial covariates
#'  for the main part of the model
#'
#' @param family the distribution family for the gam
#'
#' @param condition an optional vector of 1s and 0s of the same length as
#'  the number of records in \code{coords} and \code{response} and stating
#'  whether the record should also be modelled using covariates in
#'  \code{condition_covs} (1 if so and 0 if not). This enables the construction
#'  of slightly more complex models, such as those with an explicitly modelled
#'  observation process. This is achieved by passing \code{condition} to the
#'  \code{by} argument in \code{mgcv::s} when fitting smooths for
#'  the condition covariates, as well as adding the condition as an intercept.
#'
#' @param condition_covs an optional \code{Raster*} object giving the spatial covariates
#'  for the conditional part of the model
#'
#' @param extra_terms an optional formula object (of the form \code{~ s(x, k = 2)}
#'  or similar which can be concatenated onto the model formula)
#'  specifying further model components (in \code{extra_data}) not provided in
#'  the spatial covariates.
#'
#' @param extra_data an optional dataframe giving the covariates referred to in
#'  \code{extra_terms}
#'
#' @param bam whether to fit the model using \code{mgcv::bam} (the default),
#'  otherwise \code{mgcv::gam} is used instead
#'
#' @param bs a two letter character string indicating the (penalized) smoothing
#'   basis to use. (eg "tp" for thin plate regression spline, "cr" for cubic
#'   regression spline). See \code{\link[mgcv]{smooth.terms}} for an over view
#'   of what is available.
#'
#' @param predict whether to transform the rasters after fitting the model.
#'  If set to \code{FALSE} this can enable model tweaking before the final
#'  transformations are applied, without the computaitonal cost of prediction
#'
#' @param \dots other arguments to be passed to \code{mgcv::bam} or
#'  \code{mgcv::gam}
#'
#' @return a three-element named list containing:
#'  \itemize{
#'    \item{model}{the fitted \code{bam} or \code{gam} model object}
#'    \item{trans}{if \code{predict = TRUE} a \code{Raster*} object of the
#'     same extent, resolution and number of layers as \code{covs}, but with
#'     the values of each layer having been optimally spline-transformed.
#'     Otherwise \code{NULL}}
#'    \item{trans_cond}{if \code{predict = TRUE} and \code{condition} is not
#'     \code{NULL} a \code{Raster*} object of the same extent, resolution and
#'     number of layers as \code{condition_covs}, but with the values of each layer
#'     having been optimally spline-transformed. Otherwise \code{NULL}}
#'  }
#'
#' @export
#' @import mgcv
#' @import raster
#'
gamTrans <- function(coords,
                     response,
                     covs,
                     family = gaussian,
                     condition = NULL,
                     condition_covs = NULL,
                     extra_terms = NULL,
                     extra_data = NULL,
                     bam = TRUE,
                     bs = 'tp',
                     predict = TRUE,
                     ...) {

  # whether there's a conditional bit
  cond <- !is.null(condition)

  stopifnot(inherits(extra_terms, 'formula'))

  # check inputs
  stopifnot(inherits(covs, 'Raster'))
  if (cond)
    stopifnot(inherits(condition_covs, 'Raster'))

  # add 'cond_' onto the conditional covariate names to prevent naming conflicts
  # with the disease model
  if (cond)
    names(condition_covs) <- paste0('cond_', names(condition_covs))

  # get covariate names
  cov_names <- names(covs)
  if (cond)
    cond_names <- names(condition_covs)


  # ~~~~~~~~~~~~~
  # build formula

  cov_terms_string <- paste(sprintf('s(%s, bs = "%s")',
                                    cov_names,
                                    bs),
                            collapse = ' + ')

  cov_terms <- reformulate(cov_terms_string)

  f <- response ~ 1
  f <- f + cov_terms

  # if required, add conditional terms
  if (cond) {
    cov_terms_string <- paste(sprintf('s(%s, bs = "%s", by = condition)',
                                      cond_names,
                                      bs),
                              collapse = ' + ')

    f <- f + cond_terms + ~ condition_intercept
  }

  # if required, add extra terms
  if (!is.null(extra_terms))
    f <- f + extra_terms

  # ~~~~~~~~~~~~~
  # get training data

  # extract covariates
  data <- data.frame(extract(covs, coords))

  # optionally combine this with the conditional and extra data
  if (cond)
    data <- cbind(data,
                  data.frame(extract(condition_covs,
                                     coords)),
                  condition_intercept = condition,
                  condition)
  if (!is.null(extra_data))
    data <- cbind(data,
                  extra_data)

  # find any missing values and remove corresponding rows
  rem_idx <- badRows(data)

  data <- data[!rem_idx, ]

  if (is.vector(response)) {
    response <- response[!rem_idx]
  } else {
    response <- response[!rem_idx, ]
  }

  # fit the model
  if (bam) {
    m <- mgcv::bam(f, data = data, family = family, ...)
  } else {
    m <- mgcv::gam(f, data = data, family = family, ...)
  }

  # ~~~~~~~~~~
  # optionally apply transformations

  if (predict) {

    # find index for non-missing cells
    cell_idx <- cellIdx(covs)

    # transform the main covariates

    # extract covariate values
    vals <- raster::extract(covs, cell_idx)

    # convert to a data.frame
    vals <- data.frame(vals)

    # add on the extra data
    vals <- cbind(vals, extra_data[rep(1, nrow(vals)), , drop = FALSE])

    # optionally add on dummy data for the condition intercept and variables
    if (cond) {
      cond_data <- data.frame(matrix(0,
                                     nrow = 1,
                                     ncol = length(cond_names)))
      names(cond_data) <- cond_names
      vals <- cbind(vals,
                    condition = 0,
                    condition_intercept = 0,
                    cond_data)
    }

    # get the transformations of these values
    vals_trans <- predict(m,
                          newdata = vals,
                          type = 'terms')

    # find any condition terms
    cond_terms_idx <- grep('):condition$', colnames(vals_trans))

    # remove the condition ones and the condition index
    if (length(cond_terms_idx) > 0) {
      vals_trans <- vals_trans[, -cond_terms_idx]
    }
    vals_trans <- vals_trans[, !(colnames(vals_trans) %in% c('condition', 'condition_intercept'))]

    # format the names
    colnames(vals_trans) <- gsub(')', '', colnames(vals_trans))
    colnames(vals_trans) <- gsub('^s\\(', '', colnames(vals_trans))

    # keep only the names that are in the raster
    vals_trans <- vals_trans[, colnames(vals_trans) %in% cov_names]

    # set new raster values
    trans_ras <- insertRaster(raster = covs,
                              new_vals = vals_trans,
                              idx = cell_idx)

    # optionally apply the transformation to the conditional terms

    if (cond) {

      # remove the previous variables
      rm(vals, vals_trans)

      # extract covariate values
      vals <- raster::extract(condition_covs, cell_idx)

      # convert to a data.frame
      vals <- data.frame(vals)

      # add on the extra data
      vals <- cbind(vals, extra_data[rep(1, nrow(vals)), ])

      # add on dummy data for the main variables
      main_data <- data.frame(matrix(0,
                                     nrow = 1,
                                     ncol = length(cov_names)))

      names(main_data) <- cov_names

      vals <- cbind(vals,
                    condition = 1,
                    condition_intercept = 0,
                    main_data)

      # get the transformations of these values
      vals_trans <- predict(m,
                            newdata = vals,
                            type = 'terms')

      # keep only the condition ones
      vals_trans <- vals_trans[, grep('):condition$', colnames(vals_trans))]

      # format the names
      colnames(vals_trans) <- gsub('):condition$', '', colnames(vals_trans))
      colnames(vals_trans) <- gsub('^s\\(', '', colnames(vals_trans))

      # keep only the names that are in the raster
      vals_trans <- vals_trans[, colnames(vals_trans) %in% cond_names]

      # remove the `cond_` bit
      colnames(vals_trans) <- gsub('^cond_', '', colnames(vals_trans))

      # set new raster values
      trans_cond_ras <- insertRaster(raster = condition_covs,
                                     new_vals = vals_trans,
                                     idx = cell_idx)

    } else {
      trans_cond_ras <- NULL
    }

  } else {
    # if not poredicting, set these to NULL
    trans_ras <- trans_cond_ras <- NULL

  }

  # return the three components
  return (list(model = m,
               trans = trans_ras,
               trans_cond = trans_cond_ras))

}
