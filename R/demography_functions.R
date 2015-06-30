# functions for demographic analysis

# function to calculate mortality rates
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
                      glm = FALSE) {

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

  if (glm) {

    # if using a glm, set up dataframe
    df <- data.frame(died = NA,
                     exposed = NA,
                     window = factor(NA),
                     cohort = factor(NA),
                     cluster = factor(NA),
                     weight = NA)

    # remove dummy row
    df <- df[0, ]

  }

  # loop through cohorts
  for(cohort in c('A', 'B', 'C')) {

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
    deaths_cohort <- (exposed > 0 &  # actually exposed this time
                        age_death <= upper_mat &  # died before end of window
                        age_death >= lower_mat)  # died after start of window

    # get the cohort weight
    weight <- ifelse(cohort == 'B', 0.5, 0.25)

    if (glm) {

      # if glm probabilities wanted, format training data

      # first aggregate exposures/deaths by cluster
      exposed_agg <- aggMatrix(exposed_cohort, cluster_id)
      deaths_agg <- aggMatrix(deaths_cohort, cluster_id)

      # convert to long format by stacking on top of each other
      # if using a glm, set up dataframe
      df_tmp <- data.frame(died = as.vector(exposed_agg),
                       exposed = as.vector(deaths_agg),
                       window = rep(1:nw, each = n),
                       cluster = rep(clusters, nw),
                       weight = weight)

      # add to the master matrix
      df <- rbind(df, df_tmp)

    } else {

      # otherwise accumulate these raw numbers with weights
      exposed <- exposed + exposed_cohort * weight
      deaths <- deaths + deaths_cohort * weight

    }

  } # cohort loop

  if (glm) {

    # fit the glm

    # define the formula
    f <- died ~ 1 + f(window, model = 'iid') + f(cluster, model = 'iid')

    # enable weights in inla
    inla.setOption("enable.inla.argument.weights", TRUE)

    # fit the model
    m <- inla(f,
              data = df,
              family = 'binomial',
              weights = df$weight,
              Ntrials = df$exposed,
              control.compute = list(config = TRUE))

    # switch weights off again
    inla.setOption("enable.inla.argument.weights", FALSE)

    # get fitted mortality probabilities
    p <- m$summary.fitted.values$mode

    # keep only on cohort's amount (should all be the same)
    p <- p[1:(length(p) / 3)]

    # reformat to a matrix
    p <- matrix(p, ncol = nw)

    # get the fitted values
    survival_rates <- 1 - p  # matrix of cluster by window survival probabilities

  } else {

    # check they're all sane
    stopifnot(all(exposed <= 1))
    stopifnot(all(deaths <= 1))

    # aggregate by cluster
    exposed_agg <- aggMatrix(exposed, cluster_id)
    deaths_agg <- aggMatrix(deaths, cluster_id)

    # get mortality rates for each window, for each cluster
    survival_rates <- 1 - deaths_agg / exposed_agg

  }

  # get results
  ans <- list()
  # loop through required ages
  for (i in 1:na) {

    index <- windows_lower < ages_upper[i] &
      windows_upper > ages_lower[i]

    mortality <- 1 - rowProds(survival_rates[, index])
    sample_size <- rowSums(exposed_agg[, index])

    ans[[i]] <- data.frame(mortality, sample_size)

  }

  # add names to list
  names(ans) <- paste0('ages_', ages_lower, '_', ages_upper)

  # return this
  return (ans)

}

