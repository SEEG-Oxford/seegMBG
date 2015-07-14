# functions for demographic analysis

#' @name periodMortality
#' @rdname periodMortality
#'
#' @title Synthetic Cohort Life Table Child Mortality Estimation
#' @description Compute child mortality estimates using the DHS three-cohort
#'  life table approach.
#'
#' @param age_death numeric vector giving the age of death in months of each
#'  child. Note that if the child is still alive at the time of interview,
#'  corresponding eleents of \code{age_death} should be given a number greater
#'  than \code{max(windows_upper)}
#' @param birth_int numeric vector giving the time in month between the
#'  child's birth and the interview
#' @param cluster_id vector giving the cluster identifier (can be numeric,
#'  cahracter or factor)
#' @param windows_lower,windows_upper numeric vectors giving the
#'  non-overlapping lower and upper ages in months of each survival window
#' @param ages_lower,ages_upper numeric vectors giving the
#'  non-overlapping lower and upper ages in months of the age windows for
#'  which to estimate mortality rates
#' @param period the length of time in months for which mortality rates
#'  should be estimated, either vector or scalar
#' @param delay the length of time in months prior to the interview date
#'  to end the period, either vector or scalar.
#'  I.e. the period runs from \code{period + delay}
#'  months to \code{delay} months before the interview month.
#' @param glm whether to infer the window-specific survival probabilities
#'  using a binomial random effects model across cluster, window and cohort.
#'  If \code{FALSE} probabilities are calculated as the raw ratio of the
#'  number that survived to the number exposed and may therefore contain
#'  zeros and indeterminate values.
#' @param verbose whether to regularly report the stage of the analysis
#' @param \dots other arguments to pass to \code{INLA::inla}
#'
#' @export
#'
#' @import INLA
#'
#' @return a list with the same number of elements as \code{ages}, each
#'  element being a dataframe containing the estimated mortality rates
#'  and effective sample size (weighted number exposed in each suvival
#' window).
#'
periodMortality <- function (age_death,
                             birth_int,
                             cluster_id,
                             windows_lower = c(0, 1, 3, 6, 12, 24, 36, 48),
                             windows_upper = c(0, 2, 5, 11, 23, 35, 47, 59),
                             ages_lower = c(0, 0),
                             ages_upper = c(12, 60),
                             period = 60,
                             delay = max(windows_upper - windows_lower),
                             glm = FALSE,
                             verbose = TRUE,
                             ...) {

  if (verbose) message('formatting data')

  n <- length(age_death)
  nw <- length(windows_lower)
  na <- length(ages_lower)

  # expand period and delay if needed
  if(length(period) == 1) period <- rep(period, n)
  if(length(delay) == 1) delay <- rep(delay, n)

  # check all argument are the right size
  stopifnot(length(birth_int) == n)
  stopifnot(length(cluster_id) == n)
  stopifnot(length(period) == n)
  stopifnot(length(delay) == n)
  stopifnot(length(windows_upper) == nw)
  stopifnot(length(ages_upper) == na)

  # get unique clusters
  clusters <- unique(cluster_id)

  if (any(windows_upper[-1] <= windows_lower[-nw])) {
    stop ('windows_upper and windows_lower appear to overlap')
  }

  if (any(ages_upper[-1] <= ages_lower[-na])) {
    stop ('ages_upper and ages_lower appear to overlap')
  }

  # window age limit matrices
  upper_mat <- t(expand(windows_upper, n))
  lower_mat <- t(expand(windows_lower, n))

  # delay and period matrices
  delay_mat <- expand(delay, nw)
  period_mat <- expand(period, nw)

  # data matrices
  age_death <- expand(age_death, nw)
  birth_int <- expand(birth_int, nw)

  # empty objects to store total deaths and exposures for each window
  deaths <- exposed <- 0

  # loop through cohorts
  for(cohort in c('A', 'B', 'C')) {

    if (verbose) message(paste('processing cohort', cohort))

    # cohort end date matrix
    start_mat <- switch(cohort,
                        B = upper_mat,
                        A = period_mat + lower_mat,
                        C = lower_mat )

    # cohort start date matrix
    end_mat <- switch(cohort,
                      B = period_mat + lower_mat,
                      A = period_mat + upper_mat,
                      C = upper_mat)

    # and add gap between end of period and the interview date
    start_mat <- start_mat + delay_mat
    end_mat <- end_mat + delay_mat

    # add effective number exposed
    exposed_cohort <- (birth_int <= end_mat &  # entered cohort before cohort end date
                         birth_int >= start_mat &  # entered cohort after cohort start date
                         age_death >= lower_mat) # alive at start of cohort

    # and effective number that died
    deaths_cohort <- (exposed_cohort > 0 &  # actually exposed this time
                        age_death <= upper_mat &  # died before end of window
                        age_death >= lower_mat)  # died after start of window

    # get the cohort weight
    weight <- ifelse(cohort == 'B', 1, 0.5)

    # otherwise accumulate these raw numbers with weights
    exposed <- exposed + exposed_cohort * weight
    deaths <- deaths + deaths_cohort * weight

  } # cohort loop

  # check they're all sane
  stopifnot(all(exposed <= 2))
  stopifnot(all(deaths <= 2))

  # aggregate by cluster
  exposed_agg <- aggMatrix(exposed, cluster_id)
  deaths_agg <- aggMatrix(deaths, cluster_id)

  if (glm) {

    # set up the glm data frame
    if (verbose) message('formatting data for glm')

    # round up the deaths and exposures
    exposed_agg <- round(exposed_agg)
    deaths_agg <- round(deaths_agg)

    # create long format dataframe by stacking on top of each other
    # if using a glm, set up dataframe
    df <- data.frame(died = as.vector(deaths_agg),
                     exposed = as.vector(exposed_agg),
                     window = rep(1:nw, each = length(clusters)),
                     cluster = rep(clusters, nw))

    if (verbose) message('running glm')

    # define the formula
    f <- died ~ 1 + f(window, model = 'iid') + f(cluster, model = 'iid')

    # fit the model
    m <- inla(f,
              data = df,
              family = 'binomial',
              Ntrials = df$exposed,
              control.predictor = list(compute = TRUE),
              ...)

    # get fitted mortality probabilities
    p <- m$summary.fitted.values$mode

    # reformat to a matrix
    p <- matrix(p, ncol = nw)

    # get the fitted values
    survival_rates <- 1 - p  # matrix of cluster by window survival probabilities

  } else {

    # get raw survival rates for each window, for each cluster
    survival_rates <- 1 - deaths_agg / exposed_agg

  }

  if (verbose) message('computing survival probabilities')

  # get results
  ans <- list()

  # loop through required ages
  for (i in 1:na) {

    index <- windows_lower < ages_upper[i] &
      windows_upper > ages_lower[i]

    mortality <- 1 - rowProds(survival_rates[, index, drop = FALSE])
    sample_size <- rowSums(exposed_agg[, index, drop = FALSE])

    ans[[i]] <- data.frame(mortality, sample_size)

  }

  # add names to list
  names(ans) <- paste0('ages_', ages_lower, '_', ages_upper)

  # return this
  return (ans)

}


# function to count number exposed/died for each cluster in each
# age group
#' @name periodTabulate
#' @rdname periodTabulate
#'
#' @title Tabulate Exposures/Deaths from CBH records
#' @description Count the number exposed and number dying
#'  in a given set of age ranges and time periods
#'
#' @param age_death numeric vector giving the age of death in months of each
#'  child. Note that if the child is still alive at the time of interview,
#'  corresponding eleents of \code{age_death} should be given a number greater
#'  than \code{max(windows_upper)}
#' @param birth_int numeric vector giving the time in month between the
#'  child's birth and the interview
#' @param cluster_id vector giving the cluster identifier (can be numeric,
#'  cahracter or factor)
#' @param windows_lower,windows_upper numeric vectors giving the
#'  non-overlapping lower and upper ages in months of each survival window
#' @param period the length of time in months for which mortality rates
#'  should be estimated, either vector or scalar
#' @param nperiod the number of consecutive periods to calculate for. I.e.
#'  if \code{period = 12} and \code{nperiod = 3}, numbers will be returned for
#'  periods 0-12, 13-24 and 25-36 months prior to the interview (plus delay).
#' @param delay the length of time in months prior to the interview date
#'  to end the period, either vector or scalar.
#'  I.e. the period runs from \code{period + delay}
#'  months to \code{delay} months before the interview month.
#' @param verbose whether to regularly report the stage of the analysis
#'
#' @export
#'
#' @return a dataframe with number of rows equal to all combinations of
#'  clusters, age bins and periods; i.e.
#'  \code{length(unique(cluster_id)) * length(windows_lower) * nperiod},
#'  and four columns:
#'  \itemize{
#'    \item \code{cluster_id} identifier for the cluster
#'    \item \code{exposed} number exposed in a given period and age bin
#'    \item \code{died} number exposed and died in a given period and age bin
#'    \item \code{period} the corresponding period number
#'    \item \code{age_bin} the corresponding age bin number
#'  }
#'
periodTabulate <- function (age_death,
                            birth_int,
                            cluster_id,
                            windows_lower = c(0, 1, 3, 6, 12, 24, 36, 48),
                            windows_upper = c(0, 2, 5, 11, 23, 35, 47, 59),
                            period = 60,
                            nperiod = 1,
                            delay = max(windows_upper - windows_lower),
                            verbose = TRUE) {

  # get unique clusters
  clusters <- unique(cluster_id)

  # sizes of things
  n <- length(age_death)
  nw <- length(windows_lower)
  np <- nperiod
  ncl <- length(clusters)

  # expand period and delay if needed
  if(length(period) == 1) period <- rep(period, n)
  if(length(delay) == 1) delay <- rep(delay, n)

  # check all argument are the right size
  stopifnot(length(birth_int) == n)
  stopifnot(length(cluster_id) == n)
  stopifnot(length(period) == n)
  stopifnot(length(delay) == n)
  stopifnot(length(windows_upper) == nw)

  if (any(windows_upper[-1] <= windows_lower[-nw])) {
    stop ('windows_upper and windows_lower appear to overlap')
  }

  # window age limit matrices
  upper_mat <- t(expand(windows_upper, n))
  lower_mat <- t(expand(windows_lower, n))

  # data matrices
  age_death <- expand(age_death, nw)
  birth_int <- expand(birth_int, nw)

  # delay and period matrices
  delay_mat <- expand(delay, nw)
  period_mat <- expand(period, nw)

  # set up dataframe to store the results
  ans <- data.frame(cluster_id = rep(clusters, each = nw * np),
                    exposed = rep(NA, ncl * nw * np),
                    died = NA,
                    period = rep(rep(1:np, each = nw), ncl),
                    age_bin = rep(1:nw, ncl * np))

  for (p in 1:np) {

    # notify the user
    if (verbose & p > 1) message(paste('processing period', p))


    # get the extra delay matrix
    extra_delay_mat <- period_mat * (p - 1)

    # empty objects to store total deaths and exposures for each window
    deaths <- exposed <- 0

    # loop through cohorts
    for (cohort in c('A', 'B', 'C')) {

      # notify the user
      if (verbose) message(paste('processing cohort', cohort))

      # cohort end date matrix
      start_mat <- switch(cohort,
                          B = upper_mat,
                          A = period_mat + lower_mat,
                          C = lower_mat )

      # cohort start date matrix
      end_mat <- switch(cohort,
                        B = period_mat + lower_mat,
                        A = period_mat + upper_mat,
                        C = upper_mat)

      # add the delays - gap between end of period and interview date
      # and extra delays for different periods
      start_mat <- start_mat + delay_mat + extra_delay_mat
      end_mat <- end_mat + delay_mat + extra_delay_mat

      # add effective number exposed
      exposed_cohort <- (birth_int <= end_mat &  # entered cohort before cohort end date
                           birth_int >= start_mat &  # entered cohort after cohort start date
                           age_death >= lower_mat) # alive at start of cohort

      # and effective number that died
      deaths_cohort <- (exposed_cohort > 0 &  # actually exposed this time
                          age_death <= upper_mat)  # died before end of window

      # get the cohort weight
      weight <- ifelse(cohort == 'B', 1, 0.5)

      # otherwise accumulate these raw numbers with weights
      exposed <- exposed + exposed_cohort * weight
      deaths <- deaths + deaths_cohort * weight

    } # cohort loop

    # check they're all sane
    stopifnot(all(exposed <= 2))
    stopifnot(all(deaths <= 2))

    # aggregate by cluster
    exposed_agg <- aggMatrix(exposed, cluster_id)
    deaths_agg <- aggMatrix(deaths, cluster_id)

    # add to results
    ans$exposed[ans$period == p] <- as.vector(t(exposed_agg))
    ans$died[ans$period == p] <- as.vector(t(deaths_agg))

  }

  # round the columns up
  ans$exposed <- ceiling(ans$exposed)
  ans$died <- ceiling(ans$died)

  # return the combined data
  return (ans)

}
