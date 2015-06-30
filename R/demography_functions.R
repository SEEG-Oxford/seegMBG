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
#'  should be estimated
#' @param delay the length of time in months prior to the interview date
#'  to end the period. I.e. the period runs from \code{period + delay}
#'  months to \code{delay} months before the interview month.
#'
#' @export
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
                      delay = max(windows_upper - windows_lower)) {

  n <- length(age_death)
  nw <- length(windows_lower)
  na <- length(ages_lower)

  # check all argument are the right size
  stopifnot(length(birth_int) == n)
  stopifnot(length(cluster_id) == n)
  stopifnot(length(windows_upper) == nw)
  stopifnot(length(ages_upper) == na)

  if (any(windows_upper[-1] <= windows_lower[-nw])) {
    stop ('windows_upper and windows_lower appear to overlap')
  }

  if (any(ages_upper[-1] <= ages_lower[-na])) {
    stop ('ages_upper and ages_lower appear to overlap')
  }

  # window age limit matrices
  upper_mat <- t(expand(windows_upper, n))
  lower_mat <- t(expand(windows_lower, n))

  # data matrices
  age_death <- expand(age_death, nw)
  birth_int <- expand(birth_int, nw)

  # empty objects to store total deaths and exposures for each window
  deaths <- exposed <- 0

  # loop through cohorts
  for(cohort in c('A', 'B', 'C')) {

    # cohort end date
    start <- switch(cohort,
                    B = windows_upper,
                    A = period + windows_lower,
                    C = windows_lower )

    # cohort start date
    end <- switch(cohort,
                  B = period + windows_lower,
                  A = period + windows_upper,
                  C = windows_upper)

    # turn into matrices, and add gap between end of period and the interview
    # date
    start_mat <- t(expand(start, n)) + delay
    end_mat <- t(expand(end, n)) + delay

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

    # accumulate these, with weights
    exposed <- exposed + exposed_cohort * weight
    deaths <- deaths + deaths_cohort * weight

  }

  # check they're all sane
  stopifnot(all(exposed <= 1))
  stopifnot(all(deaths <= 1))

  # aggregate by cluster
  exposed_agg <- aggMatrix(exposed, cluster_id)
  deaths_agg <- aggMatrix(deaths, cluster_id)

  # get mortality rates for each window, for each cluster
  survival_rates <- 1 - deaths_agg / exposed_agg

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

