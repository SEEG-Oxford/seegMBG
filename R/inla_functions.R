# functions for prediction from INLA MBG models

#' @name predictINLA
#' @rdname predictINLA
#'
#' @title Prediction from INLA MBG models
#'
#' @description Given a fitted \code{\link{inla}} object for a spatial or
#'  spatio-temporal geostatistical model, make predictions to a new dataset
#'  either from the maximum \emph{a posterior} parameter set, or as samples
#'  from the predictive posterior.
#'
#' @param inla a fitted \code{\link{inla}} object with fixed effects terms
#'  and a spatial or spatio-temporal random effect model. If
#'  \code{method = "sample"} then \code{inla} must have been fitted with the
#'  argument \code{control.compute = list(config=TRUE)}.
#' @param data a dataframe giving all of the covariates for the fixed effects
#'  part of the model and the temporal term is applicable.
#' @param mesh the \code{\link{inla.mesh}} object used to construct the spatial
#'  random effect.
#' @param coords a character vector of length two giving the names of the
#'  columns in \code{data} which contain the coordinates of locations at which
#'  to make predictions.
#' @param type the type of prediction required. The default is on the scale of
#'  the linear predictors; the alternative \code{"response"} is on the scale of
#'  the response variable. Thus for a default binomial model the default
#'  predictions are of log-odds (probabilities on logit scale) and
#'  \code{type = "response"} gives the predicted probabilities.
#' @param method whether to draw \code{n} samples from the predict posterior
#'  (if \code{method = "sample"}, the default) or to make predictions at the
#'  maximum \emph{a posteriori} estimates of the model parameters (if
#'  \code{method = "MAP"}).  If \code{method = "sample"} then \code{inla} must
#'  have been fitted with the argument
#'  \code{control.compute = list(config=TRUE)}.
#' @param n if \code{method = "sample"}, the number of samples to draw from the
#'  predictive posterior.
#' @param ncpu the number of cores to use to make predictions. If
#'  \code{ncpu = 1} then predicitons will be made sequentially, otherwise for
#'  \code{n > 1} a snowfall cluster will be initiated and predictions run in
#'  parallel.
#'
#' @return a matrix with the same number of rows as \code{data} giving the
#'  predicted values at these locations
#'
predictINLA <- function(inla,
                        data,
                        mesh,
                        coords = c('lat', 'long'),
                        type = c('link', 'response'),
                        method = c('sample', 'MAP'),
                        n = 1,
                        ncpu = 1) {

  # match multiple-choice arguments
  type <- match.arg(type)
  method <- match.arg(method)

  # check incoming data types
  stopifnot(inherits(inla, 'inla'))
  stopifnot(inla$model.random == 'SPDE2 model')
  stopifnot(inherits(data, 'data.frame'))
  stopifnot(inherits(mesh, 'inla.mesh'))

  # check incoming data shapes
  stopifnot(length(coords) == 2 & all(coords %in% names(data)))
  stopifnot(all(inla$names.fixed %in% names(data)))

  # ~~~~~~~~~~~~~~~~~~
  # get required data

  # ~~~
  # get link functions
  link <- inla$misc$linkfunctions$names

  # throw an error if there is more than one
  if (length(link) > 1) stop ('multiple link functions not yet supported')

  # otherwise, get the function
  link <- get(paste0('inla.link.', link))

              # ~~~
              # extract coordinates
              coords <- data[, coords]

              # ~~~
              # find spatial term
              spatial_term <- names(inla$summary.random)
              stopifnot (length(spatial_term) == 1)

              # ~~~
              # get parameters

              params <- getINLAParameters(inla = inla,
                                          method = method,
                                          n = n)


              # make predictions
              # loop, optionally in parallel

              # optionally switch to response scale
              if (type == 'response')
                results <- link(results, inverse = TRUE)

              # return results
              return (results)

}

#' @name getINLAParameters
#' @rdname getINLAParameters
#'
#' @title Extracting parameters form INLA MBG models
#'
#' @description Given a fitted \code{\link{inla}} object for a spatial or
#'  spatio-temporal geostatistical model, extract model parameters
#'  either from the maximum \emph{a posterior} parameter set, or as samples
#'  from the predictive posterior.
#'
#' @param inla a fitted \code{\link{inla}} object with fixed effects terms
#'  and a spatial or spatio-temporal random effect model. If
#'  \code{method = "sample"} then \code{inla} must have been fitted with the
#'  argument \code{control.compute = list(config=TRUE)}.
#' @param method whether to draw \code{n} samples from the predict posterior
#'  (if \code{method = "sample"}, the default) or to make predictions at the
#'  maximum \emph{a posteriori} estimates of the model parameters (if
#'  \code{method = "MAP"}).  If \code{method = "sample"} then \code{inla} must
#'  have been fitted with the argument
#'  \code{control.compute = list(config=TRUE)}.
#' @param n if \code{method = "sample"}, the number of samples to draw from the
#'  predictive posterior.
#'
#' @return a list with three elements:
#'  \itemize{
#'    \item params parameters
#'    \item names parameter names
#'    \item indices
#'      \itemize{
#'        \item fixed_idx index to fixed effect parameters in params and names
#'        \item spatial_idx index to spatial random effect parameters in params and names
#'        \item hyper_idx index to hyperparameters in params and names
#'      }
#'  }
#'
#'
getINLAParameters <- function (inla, method = c('sample', 'MAP'), n = 1) {

  # get the method
  method <- match.arg(method)

  if (method == 'MAP') {
    # if MAP estimate, get a single estimate of the model parameters

    # ~~~
    # hyper parameters
    hyper_names <- rownames(inla$summary.hyperpar)
    params_hyper <- inla$summary.hyperpar$mode

    # ~~~
    # fixed parameters
    fixed_names <- rownames(inla$summary.fixed)
    params_fixed <- inla$summary.fixed$mode

    # ~~~
    # spatial parameters
    spatial_names <- rownames(inla$summary.random[[1]])
    params_spatial <- inla$summary.random[[1]]$mode

    # ~~~
    # combine the parameters
    fixed_idx <- seq_along(params_fixed)
    spatial_idx <- seq_along(params_spatial) + length(params_fixed)
    hyper_idx <- seq_along(params_hyper) + length(params_fixed) + length(params_spatial)

    param_names <- c(fixed_names, spatial_names, hyper_names)
    params <- c(params_fixed, params_spatial, params_hyper)
    params <- list(params)


  } else if (method == 'sample') {

    # if samples, get multiple draws of the model parameters
    suppressWarnings(params <- inla.posterior.sample(inla,
                                                     n = n))

    # ~~~
    # hyper parameters

    # extract the hypers
    hypers <- lapply(params,
                     '[[',
                     'hyperpar')

    # keep their names
    hyper_names <- names(hypers[[1]])

    # ~~~
    # latent parameters
    params <- lapply(params,
                     '[[',
                     'latent')

    # keep parameter names
    param_names <- rownames(params[[1]])

    # find training data predictors
    predictor_idx <- grep('^Predictor.',
                          param_names)

    # remove them from parameter sets and names
    param_names <- param_names[-predictor_idx]
    params <- lapply(params,
                     '[',
                     -predictor_idx,
                     drop = FALSE)

    # find the spatial indices
    spatial_prefix <- paste0('^', spatial_term, '.')
    spatial_idx <- grep(spatial_prefix,
                        param_names)

    # and the fixed indices should be those that remain
    fixed_idx <- (1:length(param_names))[-spatial_idx]

    # now tidy up the parameter names
    param_names[spatial_idx] <- gsub(spatial_prefix, '', param_names[spatial_idx])
    param_names[fixed_idx] <- gsub('.1$', '', param_names[fixed_idx])

    # convert spatial names to condensed numeric names (i.e. remove leading 0s)
    param_names[spatial_idx] <- as.character(as.numeric(param_names[spatial_idx]))

    # ~~~
    # combine them
    hyper_idx <- max(fixed_idx) + seq_along(hyper_names)
    param_names[hyper_idx] <- hyper_names
    for (i in 1:length(params)) params[[i]] <- c(params[[i]], hypers[[i]])

  } else {

    stop ('invalid method')

  }

  # return the results as a list
  ans <- list(params,
              names = param_names,
              indices = list(fixed_idx,
                             spatial_idx,
                             hyper_idx))

  return (ans)

}
