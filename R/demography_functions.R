# functions for demographic analysis

#' @name periodMortality
#' @rdname periodMortality
#'
#' @template period_args
#'
#' @param ages_lower,ages_upper numeric vectors giving the
#'  non-overlapping lower and upper ages in months of the age windows for
#'  which to estimate mortality rates
#'
#' @param glm whether to infer the window-specific survival probabilities
#'  using a binomial random effects model across cluster, window and cohort.
#'  If \code{FALSE} probabilities are calculated as the raw ratio of the
#'  number that survived to the number exposed and may therefore contain
#'  zeros and indeterminate values.
#'
#' @param \dots other arguments to pass to \code{INLA::inla}
#'
#' @return \code{periodMortality}: a list with one element for each period,
#'  each element being a dataframe with the same number of columns as elements
#'  in \code{ages_lower}, each column giving the estimated mortality rates for
#'  the corresponding age bin in each cluster.
#'
#' @export
#' @import INLA
#'
periodMortality <- function (age_death,
                             birth_int,
                             cluster_id,
                             windows_lower = c(0, 1, 3, 6, 12, 24, 36, 48),
                             windows_upper = c(0, 2, 5, 11, 23, 35, 47, 59),
                             ages_lower = c(0, 0),
                             ages_upper = c(12, 60),
                             period = 60,
                             method = c('monthly', 'direct'),
                             cohorts = c('one', 'three'),
                             inclusion = c('enter', 'exit', 'both', 'either'),
                             mortality = c('bin', 'monthly'),
                             nperiod = 1,
                             delay = max(windows_upper - windows_lower),
                             glm = FALSE,
                             verbose = TRUE,
                             n_cores = 1,
                             ...) {

  # check the number of cores is valid


  # throw a warning is monthly rates are wanted, but monthly mortalities
  # not used to calculate them
  if (mortality =='monthly' && method != 'monthly') {
    stop ("The argument mortality = 'monthly' can only be used if method = 'monthly'")
  }

  # check incoming data
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

  if (any(windows_upper[-1] <= windows_lower[-nw])) {
    stop ('windows_upper and windows_lower appear to overlap')
  }

  if (any(ages_upper[-1] <= ages_lower[-na])) {
    stop ('ages_upper and ages_lower appear to overlap')
  }

  # tabulate data for analysis
  df_all <- periodTabulate(age_death = age_death,
                           birth_int = birth_int,
                           cluster_id = cluster_id,
                           windows_lower = windows_lower,
                           windows_upper = windows_upper,
                           period = period,
                           method = method,
                           cohorts = cohorts,
                           inclusion = inclusion,
                           mortality = mortality,
                           nperiod = nperiod,
                           delay = delay,
                           verbose = verbose,
                           n_cores = n_cores)


  # loop through periods

  # results for all periods
  res <- list()

  for (j in 1:nperiod) {

    df <- df_all[df_all$period == j, ]

    # fit the required model
    if (glm) {

      if (n_cores == 1) {
        inla_threads <- NULL
      } else {
        inla_threads <- n_cores
      }

      # round up the deaths and exposures
      df$exposed <- ceiling(df$exposed)
      df$died <- ceiling(df$died)

      # make sure there aren't more deaths than exposures, from rounding errors
      df$exposed <- pmax(df$exposed, df$died)

      # define the formula
      f <- died ~ 1 + f(age_bin, model = 'iid') + f(cluster_id, model = 'iid')

      # notify the user the model's running
      if (verbose) message('running glm')

      # fit the model
      m <- inla(f,
                data = df,
                family = 'binomial',
                Ntrials = df$exposed,
                control.predictor = list(compute = TRUE),
                num.threads = inla_threads,
                ...)

      # get fitted mortality probabilities
      p <- m$summary.fitted.values$mode

      # reformat to a matrix
      p <- matrix(p, ncol = nw, byrow = TRUE)

      # get the fitted values
      survival_mat <- 1 - p  # matrix of cluster by window survival probabilities

    } else {

      # get raw survival rates for each window, for each cluster
      survival_rates <- 1 - df$died / df$exposed

      # reformat as a matrix
      survival_mat <- matrix(survival_rates,
                             ncol = nw,
                             byrow = TRUE)

    }

    if (verbose) message('computing survival probabilities')

    # get results
    ans <- list()

    # loop through age bins for which mortality estimates are needed
    for (i in 1:na) {

      # if it's a one-month window, match exactly
      if (ages_lower[i] == ages_upper[i]) {
        index <- windows_lower == ages_upper[i] &
          windows_upper == ages_lower[i]
      } else{
        index <- windows_lower <= ages_upper[i] &
          windows_upper > ages_lower[i]
      }

      # get bin mortality estimates
      if (mortality == 'bin') {

        # if the whole bin is wanted, it's the product of all probabilities
        mort <- 1 - rowProds(survival_mat[, index, drop = FALSE])

      } else if (mortality == 'monthly') {

        # if the monthly mortality estimate is needed, it's the
        # mean rate weighted by the number of months

        # get number of months in each bin in the required age bin
        months <- 1 + windows_upper[index] - windows_lower[index]

        # get weighted mean rate
        mort <- 1 - apply(survival_mat[, index, drop = FALSE],
                          1,
                          weighted.mean,
                          w = months)

      }

      ans[[i]] <- mort

    }

    # convert to a dataframe
    ans <- data.frame(ans)

    # add names
    names(ans) <- paste0('ages_', ages_lower, '_', ages_upper)
    rownames(ans) <- unique(df$cluster_id)

    res[[j]] <- ans

  }

  # name and return the results
  names(res) <- paste0('period_', 1:nperiod)

  return (res)

}


# function to count number exposed/died for each cluster in each
# age group
#' @name periodTabulate
#' @rdname periodMortality
#'
#' @export
#'
#' @return \code{periodTabulate}: a dataframe with number of rows equal to all combinations of
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
                            method = c('monthly', 'direct'),
                            cohorts = c('one', 'three'),
                            inclusion = c('enter', 'exit', 'both', 'either'),
                            mortality = c('bin', 'monthly'),
                            nperiod = 1,
                            delay = NULL,
                            verbose = TRUE,
                            n_cores = 1) {

  # parallelism by recursion
  if (n_cores > 1) {

    # let the user know
    message(sprintf('running periodTabulate on %s cores', n_cores))

    # check data sizes
    stopifnot(length(birth_int) == length(age_death))
    stopifnot(length(cluster_id) == length(age_death))

    # get all of the unique clusters
    clusters <- unique(cluster_id)

    # split up data by finding indices for clusters
    n <- length(age_death)
    indices <- parallel::splitIndices(n, n_cores)

    # convert these into the actual cluster ids
    cluster_groups <- lapply(indices, function(i, x) {x[i]}, clusters)

    # define function to act on groups of clusters
    parfun <- function (cluster_group,
                        age_death,
                        birth_int,
                        cluster_id,
                        ...) {

      idx <- which(cluster_id %in% cluster_group)

      periodTabulate(age_death = age_death[idx],
                     birth_int = birth_int[idx],
                     cluster_id = cluster_id[idx],
                     ...)

    }

    # set up cluster
    sfInit(parallel = TRUE, cpus = n_cores)

    # run chunks in parallel
    ans_list <- sfLapply(cluster_groups,
                         parfun,
                         age_death,
                         birth_int,
                         cluster_id,
                         windows_lower = windows_lower,
                         windows_upper = windows_upper,
                         period = period,
                         method = method,
                         cohorts = cohorts,
                         inclusion = inclusion,
                         mortality = mortality,
                         nperiod = nperiod,
                         delay = delay,
                         verbose = verbose,
                         n_cores = 1)

    # stop the cluster
    sfStop()

    # recombine results into ans
    ans <- do.call(rbind, ans_list)

    return (ans)
  }

  # get the tabulation method
  method <- match.arg(method)
  cohorts <- match.arg(cohorts)
  inclusion <- match.arg(inclusion)
  mortality <- match.arg(mortality)

  # if no delay is specified, get the correct one for the method
  if (is.null(delay)) {
    delay <- switch(cohorts,
                    one = 0,
                    three = max(windows_upper - windows_lower))
  }

  # match up the cohorts
  cohort_names <- switch(cohorts,
                         one = 'B',
                         three = c('A', 'B', 'C'))

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

  # set up dataframe to store the results
  ans <- data.frame(cluster_id = rep(clusters, each = nw * np),
                    exposed = rep(NA, ncl * nw * np),
                    died = NA,
                    period = rep(rep(1:np, each = nw), ncl),
                    age_bin = rep(1:nw, ncl * np))

  for (p in 1:np) {

    # notify the user
    if (verbose & np > 1) message(paste('\nprocessing period', p))

    if (method == 'monthly') {

      # loop through the age windows
      for (w in 1:nw) {

        # component monthly age windows
        new_windows <- seq(windows_lower[w],
                           windows_upper[w],
                           by = 1)

        # account for the 0th month not being made into a sequence
        if (length(new_windows) == 1 && new_windows == 0) {
          new_windows <- c(0, 0)
        }

        # number of new windows
        n_nw <- length(new_windows) - 1

        # recursively call this function with method = direct, calculating
        # numbers for monthly bins
        res_tmp <- periodTabulate(age_death = age_death,
                                  birth_int = birth_int,
                                  cluster_id = cluster_id,
                                  windows_lower = new_windows[-n_nw],
                                  windows_upper = new_windows[-1],
                                  period = period,
                                  method = "direct",
                                  cohorts = cohorts,
                                  inclusion = inclusion,
                                  mortality = mortality,
                                  nperiod = 1,
                                  delay = delay + period * (p - 1),
                                  verbose = verbose)

        # this returns cluster-level deaths & exposures for each age bin
        # next, combine them to get expected numbers of period exposures
        # and deaths

        # aggregate monthly deaths and exposures
        exposed_mnth <- tapply(res_tmp$exposed, res_tmp$cluster_id, sum)
        died_mnth <- tapply(res_tmp$died, res_tmp$cluster_id, sum)

        if (mortality == 'bin') {

          # get expected period exposures
          exposed_per <- exposed_mnth / n_nw

          # get expected period deaths
          rate <- 1 - (1 - (died_mnth / (exposed_mnth))) ^ n_nw

          died_per <- rate * exposed_per

          # set any 0 total monthly exposures to 0 in the periods
          exposed_per[exposed_mnth == 0] <- 0
          died_per[exposed_mnth == 0] <- 0

          # set exposed and died resultsto these period figures
          exposed_res <- exposed_per
          died_res <- died_per

        } else {
          # otherwise if monthly numbers needed

          # set exposed and died results to the monthly figures
          exposed_res <- exposed_mnth
          died_res <- died_mnth

        }

        # insert these into the results dataframe
        idx_insert <- which(ans$period == p & ans$age_bin == w)

        ans$exposed[idx_insert] <- exposed_res
        ans$died[idx_insert] <- died_res


      }

    } else if (method == 'direct') {

      # get matrices

      # window age limit matrices
      upper_mat <- t(expand(windows_upper, n))
      lower_mat <- t(expand(windows_lower, n))

      # age bin range matrix
      age_range_mat <- upper_mat - lower_mat + 1

      # data matrices
      age_death <- expand(age_death, nw)
      birth_int <- expand(birth_int, nw)

      # delay and period matrices
      delay_mat <- expand(delay, nw)
      period_mat <- expand(period, nw)

      # get the extra delay matrix, for when multiple periods are needed
      extra_delay_mat <- period_mat * (p - 1)

      # empty objects to store total deaths and exposures for each window
      deaths <- exposed <- 0

      # loop through cohorts
      for (cohort in cohort_names) {

        # notify the user
        if (verbose & length(cohort_names) > 1) {
          message(paste('processing cohort', cohort))
        }

        # define inclusion offsets
        # i.e. if they only need to enter the age bin in the period,
        # can reduce the upper limit of birth-to-interview time


        # offset for upper limit
        start_offset <- switch(inclusion,
                               enter = -age_range_mat,
                               exit = 0,
                               both = 0,
                               either = -age_range_mat)

        # offset for lower limit
        end_offset <- switch(inclusion,
                             enter = 0,
                             exit = age_range_mat,
                             both = 0,
                             either = age_range_mat)

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

        # add in the inclusion offsets
        start_mat <- start_mat + start_offset
        end_mat <- end_mat + end_offset

        # add effective number exposed
        exposed_cohort <- (birth_int <= end_mat &  # entered cohort before cohort end date
                             birth_int >= start_mat &  # entered cohort after cohort start date
                             age_death >= lower_mat) # alive at start of cohort

        # and effective number that died
        deaths_cohort <- (exposed_cohort > 0 &  # actually exposed this time
                            age_death <= upper_mat)  # died before end of window

        # get the cohort weight
        weight <- ifelse(cohort == 'B', 1, 0.5)

        # accumulate these raw numbers with weights
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


  }

  # return the combined data
  return (ans)

}
