# functions in this file have the same arguments structure
# Args:
#   mu: untransformed linear predictor matrix
#   draws: A named list returned by extract_draws containing 
#          all required data and samples
# Returns:
#   transformed linear predictor representing the mean
#   of the response distribution
fitted_gaussian <- function(draws) {
  fitted_default(draws)
}

fitted_student <- function(draws) {
  fitted_default(draws)
}

fitted_cauchy <- function(draws) {
  fitted_default(draws)
}

fitted_skew_normal <- function(draws) {
  fitted_default(draws)
}

fitted_lognormal <- function(draws) {
  draws$sigma <- get_sigma(
    draws$sigma, data = draws$data, dim = dim_mu(draws)
  )
  draws$mu <- ilink(draws$mu, draws$f$link)
  if (!is_trunc(draws$data)) {
    draws$mu <- with(draws, exp(mu + sigma^2 / 2))
  } else {
    draws$mu <- fitted_trunc(draws)
  }
  draws$mu
}

fitted_binomial <- function(draws) {
  trials <- as_draws_matrix(draws$data$trials, dim_mu(draws))
  draws$mu <- ilink(draws$mu, draws$f$link) 
  if (!is_trunc(draws$data)) {
    # scale mu from [0,1] to [0,trials]
    draws$mu <- draws$mu * trials 
  } else {
    draws$mu <- fitted_trunc(draws)
  }
  draws$mu
}

fitted_bernoulli <- function(draws) {
  fitted_default(draws)
}

fitted_poisson <- function(draws) {
  fitted_default(draws)
}

fitted_negbinomial <- function(draws) {
  fitted_default(draws)
}

fitted_geometric <- function(draws) {
  fitted_default(draws)
}

fitted_exponential <- function(draws) {
  fitted_default(draws)
}

fitted_gamma <- function(draws) {
  fitted_default(draws)
}

fitted_weibull <- function(draws) {
  draws$shape <- get_shape(
    draws$shape, data = draws$data, dim = dim_mu(draws)
  )
  draws$mu <- ilink(draws$mu / draws$shape, draws$f$link)
  if (!is_trunc(draws$data)) {
    draws$mu <- with(draws, mu * gamma(1 + 1 / shape))
  } else {
    draws$mu <- fitted_trunc(draws)
  }
  draws$mu
}

fitted_frechet <- function(draws) {
  fitted_default(draws)
}

fitted_gen_extreme_value <- function(draws) {
  draws$sigma <- get_sigma(
    draws$sigma, data = draws$data, dim = dim_mu(draws)
  )
  draws$xi <- get_dpar(draws$xi)
  draws$mu <- ilink(draws$mu, draws$f$link)
  if (!is_trunc(draws$data)) {
    draws$mu <- with(draws, mu + sigma * (gamma(1 - xi) - 1) / xi)
  } else {
    draws$mu <- fitted_trunc(draws)
  }
  draws$mu
}

fitted_inverse.gaussian <- function(draws) {
  fitted_default(draws)
}

fitted_exgaussian <- function(draws) {
  draws$mu <- ilink(draws$mu, draws$f$link) + get_dpar(draws$beta)
  fitted_trunc(draws)
}

fitted_wiener <- function(draws) {
  # mu is the drift rate
  draws$mu <- ilink(draws$mu, draws$f$link)
  draws$bs <- get_dpar(draws$bs)
  draws$ndt <- get_dpar(draws$ndt)
  draws$bias <- get_dpar(draws$bias)
  with(draws,
   ndt - bias / mu + bs / mu * 
     (exp(- 2 * mu * bias) - 1) / (exp(-2 * mu * bs) - 1)
  )
}

fitted_beta <- function(draws) {
  fitted_default(draws)
}

fitted_von_mises <- function(draws) {
  fitted_default(draws)
}

fitted_asym_laplace <- function(draws) {
  draws$quantile <- get_dpar(draws$quantile)
  draws$sigma <- get_sigma(
    draws$sigma, data = draws$data, dim = dim_mu(draws)
  )
  draws$mu <- ilink(draws$mu, draws$f$link)
  with(draws, 
    mu + sigma * (1 - 2 * quantile) / (quantile * (1 - quantile))
  )
}

fitted_hurdle_poisson <- function(draws) {
  draws$hu <- get_zi_hu(draws, par = "hu")
  draws$mu <- adjust_old_forked(draws$mu, draws$hu)
  draws$mu <- ilink(draws$mu, draws$f$link)
  with(draws, mu / (1 - exp(-mu)) * (1 - hu))
}

fitted_hurdle_negbinomial <- function(draws) {
  draws$shape <- get_shape(
    draws$shape, data = draws$data, dim = dim_mu(draws)
  )
  draws$hu <- get_zi_hu(draws, par = "hu")
  draws$mu <- adjust_old_forked(draws$mu, draws$hu)
  draws$mu <- ilink(draws$mu, draws$f$link)
  with(draws, mu / (1 - (shape / (mu + shape))^shape) * (1 - hu))
}

fitted_hurdle_gamma <- function(draws) {
  draws$hu <- get_zi_hu(draws, par = "hu")
  draws$mu <- adjust_old_forked(draws$mu, draws$hu)
  draws$mu <- ilink(draws$mu, draws$f$link)
  with(draws, mu * (1 - hu))
}

fitted_hurdle_lognormal <- function(draws) {
  sigma <- get_sigma(
    draws$sigma, data = draws$data, dim = dim_mu(draws)
  )
  draws$hu <- get_zi_hu(draws, par = "hu")
  draws$mu <- ilink(draws$mu, draws$f$link)
  with(draws, exp(mu + sigma^2 / 2) * (1 - hu))
}

fitted_zero_inflated_poisson <- function(draws) {
  draws$zi <- get_zi_hu(draws, par = "zi")
  draws$mu <- adjust_old_forked(draws$mu, draws$zi)
  ilink(draws$mu, draws$f$link) * (1 - draws$zi) 
}

fitted_zero_inflated_negbinomial <- function(draws) {
  draws$zi <- get_zi_hu(draws, par = "zi")
  draws$mu <- adjust_old_forked(draws$mu, draws$zi)
  ilink(draws$mu, draws$f$link) * (1 - draws$zi) 
}

fitted_zero_inflated_binomial <- function(draws) {
  draws$zi <- get_zi_hu(draws, par = "zi")
  draws$mu <- adjust_old_forked(draws$mu, draws$zi)
  draws$mu <- ilink(draws$mu, draws$f$link) * (1 - draws$zi)
  trials <- draws$data[["trials"]]
  if (!is.null(draws$data$N_trait)) {
    # deprecated as of brms 1.0.0
    J <- seq_len(ceiling(length(trials) / 2))
    trials <- trials[J]
  }
  trials <- as_draws_matrix(trials, dim_mu(draws))
  draws$mu * trials
}

fitted_zero_inflated_beta <- function(draws) {
  draws$zi <- get_zi_hu(draws, par = "zi")
  draws$mu <- adjust_old_forked(draws$mu, draws$zi)
  ilink(draws$mu, draws$f$link) * (1 - draws$zi)
}

fitted_zero_one_inflated_beta <- function(draws) {
  draws$zoi <- get_dpar(draws$zoi)
  draws$coi <- get_dpar(draws$coi)
  draws$zoi * draws$coi + 
    ilink(draws$mu, draws$f$link) * (1 - draws$zoi)
}

fitted_categorical <- function(draws) {
  fitted_catordinal(draws)
}

fitted_cumulative <- function(draws) {
  draws$disc <- get_disc(draws, ncat = draws$data[["ncat"]])
  draws$mu <- draws$disc * draws$mu
  fitted_catordinal(draws)
}

fitted_sratio <- function(draws) {
  draws$disc <- get_disc(draws, ncat = draws$data[["ncat"]])
  draws$mu <- draws$disc * draws$mu 
  fitted_catordinal(draws)
}

fitted_cratio <- function(draws) {
  draws$disc <- get_disc(draws, ncat = draws$data[["ncat"]])
  draws$mu <- draws$disc * draws$mu 
  fitted_catordinal(draws)
}

fitted_acat <- function(draws) {
  draws$disc <- get_disc(draws, ncat = draws$data[["ncat"]])
  draws$mu <- draws$disc * draws$mu 
  fitted_catordinal(draws)
}

fitted_mixture <- function(draws) {
  families <- family_names(draws$f)
  draws$theta <- get_theta(draws)
  out <- 0
  for (j in seq_along(families)) {
    fitted_fun <- paste0("fitted_", families[j])
    fitted_fun <- get(fitted_fun, asNamespace("brms"))
    dpars <- valid_dpars(families[j])
    tmp_draws <- list(
      f = draws$f$mix[[j]],
      nsamples = draws[["nsamples"]],
      data = draws[["data"]]
    )
    for (ap in dpars) {
      tmp_draws[[ap]] <- draws[[paste0(ap, j)]]
    }
    if (length(dim(draws$theta)) == 3L) {
      theta <- draws$theta[, , j]
    } else {
      theta <- draws$theta[, j]
    }
    out <- out + theta * fitted_fun(tmp_draws)
  }
  out
}

# ------ fitted helper functions ------

fitted_default <- function(draws) {
  # default fitted values
  draws$mu <- ilink(draws$mu, draws$f$link)
  draws$mu <- fitted_lagsar(draws)
  fitted_trunc(draws)
}

fitted_catordinal <- function(draws) {
  # fitted values for categorical and ordinal models
  get_density <- function(s) {
    # get probabilities of each category
    do.call(dens, c(args, list(eta = draws$mu[, s, ])))
  }
  ncat <- draws$data[["ncat"]]
  args <- list(seq_len(ncat), ncat = ncat, link = draws$f$link)
  dens <- paste0("d", draws$f$family)
  draws$mu <- abind(lapply(seq_len(ncol(draws$mu)), get_density), along = 3)
  aperm(draws$mu, perm = c(1, 3, 2))
}

fitted_lagsar <- function(draws) {
  if (!is.null(draws[["lagsar"]])) {
    stopifnot(draws$f$family %in% c("gaussian", "student"))
    .fitted_lagsar <- function(s) {
      W_new <- with(draws, diag(data$N) - lagsar[s, ] * data$W)
      as.numeric(solve(W_new) %*% draws$mu[s, ])
    }
    draws$mu <- do.call(rbind, lapply(1:draws$nsamples, .fitted_lagsar))
  }
  draws$mu
}

adjust_old_forked <- function(mu, par) {
  # for compatibility with zi / hu models 
  # using old multivariate syntax
  # Args:
  #   par: samples of zi or hu parameter
  if (isTRUE(ncol(mu) == 2 * ncol(par))) {
    mu <- mu[, seq_len(ncol(mu) / 2)]
  }
  mu
}

as_draws_matrix <- function(x, dim) {
  # expand data to dimension appropriate for
  # vectorized multiplication with posterior samples
  stopifnot(length(dim) == 2L, length(x) %in% c(1, dim[2]))
  matrix(x, nrow = dim[1], ncol = dim[2], byrow = TRUE)
}

dim_mu <- function(draws) {
  c(nrow(draws$mu), draws$data$N)
}

is_trunc <- function(data) {
  any(data[["lb"]] > - Inf) || any(data[["ub"]] < Inf)
}

fitted_trunc <- function(draws) {
  # prepares data required for truncation and calles the 
  # family specific truncation function for fitted values
  if (is_trunc(draws$data)) {
    lb <- as_draws_matrix(draws$data[["lb"]], dim_mu(draws))
    ub <- as_draws_matrix(draws$data[["ub"]], dim_mu(draws))
    fitted_trunc_fun <- paste0("fitted_trunc_", draws$f$family)
    fitted_trunc_fun <- try(
      get(fitted_trunc_fun, asNamespace("brms")), 
      silent = TRUE
    )
    if (is(fitted_trunc_fun, "try-error")) {
      stop2("Fitted values on the respone scale not yet implemented ",
            "for truncated '", draws$f$family, "' models.")
    } else {
      trunc_args <- nlist(draws, lb, ub)
      draws$mu <- do.call(fitted_trunc_fun, trunc_args)
    }
  }
  draws$mu
}

# ----- family specific truncation functions -----
# Args:
#   mu: draws: output of extract_draws
#   lb: lower truncation bound
#   ub: upper truncation bound
# Returns:
#   samples of the truncated mean parameter
fitted_trunc_gaussian <- function(draws, lb, ub) {
  draws$sigma <- get_sigma(
    draws$sigma, data = draws$data, dim = dim_mu(draws)
  )
  zlb <- (lb - draws$mu) / draws$sigma
  zub <- (ub - draws$mu) / draws$sigma
  # truncated mean of standard normal; see Wikipedia
  trunc_zmean <- (dnorm(zlb) - dnorm(zub)) / (pnorm(zub) - pnorm(zlb))  
  draws$mu + trunc_zmean * draws$sigma  
}

fitted_trunc_student <- function(draws, lb, ub) {
  draws$sigma <- get_sigma(
    draws$sigma, data = draws$data, dim = dim_mu(draws)
  )
  draws$nu <- get_dpar(draws$nu)
  zlb <- (lb - draws$mu) / draws$sigma
  zub <- (ub - draws$mu) / draws$sigma
  # see Kim 2008: Moments of truncated Student-t distribution
  G1 <- with(draws,
    gamma((nu - 1) / 2) * nu^(nu / 2) / 
     (2 * (pt(zub, df = nu) - pt(zlb, df = nu))
      * gamma(nu / 2) * gamma(0.5))
  )
  A <- with(draws, (nu + zlb^2) ^ (-(nu - 1) / 2))
  B <- with(draws, (nu + zub^2) ^ (-(nu - 1) / 2))
  trunc_zmean <- G1 * (A - B)
  draws$mu + trunc_zmean * draws$sigma 
}

fitted_trunc_lognormal <- function(draws, lb, ub) {
  # mu has to be on the linear scale
  draws$sigma <- get_sigma(
    draws$sigma, data = draws$data, dim = dim_mu(draws)
  )
  m1 <- with(draws, 
    exp(mu + sigma^2 / 2) * 
      (pnorm((log(ub) - mu) / sigma - sigma) - 
       pnorm((log(lb) - mu) / sigma - sigma))
  )
  with(draws, 
    m1 / (plnorm(ub, meanlog = mu, sdlog = sigma) - 
          plnorm(lb, meanlog = mu, sdlog = sigma))
  )
}

fitted_trunc_gamma <- function(draws, lb, ub) {
  draws$shape <- get_shape(
    draws$shape, data = draws$data, dim = dim_mu(draws)
  )
  # mu becomes the scale parameter
  draws$mu <- draws$mu / draws$shape
  # see Jawitz 2004: Moments of truncated continuous univariate distributions
  m1 <- with(draws, 
    mu / gamma(shape) * 
      (incgamma(ub / mu, 1 + shape) - incgamma(lb / mu, 1 + shape))
  )
  with(draws, 
    m1 / (pgamma(ub, shape, scale = mu) - pgamma(lb, shape, scale = mu))
  )
}

fitted_trunc_exponential <- function(draws, lb, ub) {
  # see Jawitz 2004: Moments of truncated continuous univariate distributions
  # mu is already the scale parameter
  inv_mu <- 1 / draws$mu
  m1 <- with(draws, mu * (incgamma(ub / mu, 2) - incgamma(lb / mu, 2)))
  with(draws, m1 / (pexp(ub, rate = inv_mu) - pexp(lb, rate = inv_mu)))
}

fitted_trunc_weibull <- function(draws, lb, ub) {
  # see Jawitz 2004: Moments of truncated continuous univariate distributions
  # mu is already the scale parameter
  draws$shape <- get_shape(
    draws$shape, data = draws$data, dim = dim_mu(draws)
  )
  a <- 1 + 1 / draws$shape
  m1 <- with(draws,
    mu * (incgamma((ub / mu)^shape, a) - incgamma((lb / mu)^shape, a))
  )
  with(draws,
    m1 / (pweibull(ub, shape, scale = mu) - pweibull(lb, shape, scale = mu))
  )
}

fitted_trunc_binomial <- function(draws, lb, ub) {
  lb <- ifelse(lb < -1, -1, lb)
  max_value <- max(draws$data$trials)
  ub <- ifelse(ub > max_value, max_value, ub)
  trials <- draws$data$trials
  if (length(trials) > 1) {
    trials <- as_draws_matrix(trials, dim_mu(draws))
  }
  args <- list(size = trials, prob = draws$mu)
  fitted_trunc_discrete(dist = "binom", args = args, lb = lb, ub = ub)
}

fitted_trunc_poisson <- function(draws, lb, ub) {
  lb <- ifelse(lb < -1, -1, lb)
  max_value <- 3 * max(draws$mu)
  ub <- ifelse(ub > max_value, max_value, ub)
  args <- list(lambda = draws$mu)
  fitted_trunc_discrete(dist = "pois", args = args, lb = lb, ub = ub)
}

fitted_trunc_negbinomial <- function(draws, lb, ub) {
  lb <- ifelse(lb < -1, -1, lb)
  max_value <- 3 * max(draws$mu)
  ub <- ifelse(ub > max_value, max_value, ub)
  draws$shape <- get_shape(
    draws$shape, data = draws$data, dim = dim_mu(draws)
  )
  args <- list(mu = draws$mu, size = draws$shape)
  fitted_trunc_discrete(dist = "nbinom", args = args, lb = lb, ub = ub)
}

fitted_trunc_geometric <- function(draws, lb, ub) {
  lb <- ifelse(lb < -1, -1, lb)
  max_value <- 3 * max(draws$mu)
  ub <- ifelse(ub > max_value, max_value, ub)
  args <- list(mu = draws$mu, size = 1)
  fitted_trunc_discrete(dist = "nbinom", args = args, lb = lb, ub = ub)
}

fitted_trunc_discrete <- function(dist, args, lb, ub) {
  stopifnot(is.matrix(lb), is.matrix(ub))
  message(
    "Computing fitted values for truncated ", 
    "discrete models may take a while."
  )
  pdf <- get(paste0("d", dist), mode = "function")
  cdf <- get(paste0("p", dist), mode = "function")
  mean_kernel <- function(x, args) {
    # just x * density(x)
    x * do.call(pdf, c(x, args))
  }
  if (any(is.infinite(c(lb, ub)))) {
    stop("lb and ub must be finite")
  }
  # simplify lb and ub back to vector format 
  vec_lb <- lb[1, ]
  vec_ub <- ub[1, ]
  min_lb <- min(vec_lb)
  # array of dimension S x N x length((lb+1):ub)
  mk <- lapply((min_lb + 1):max(vec_ub), mean_kernel, args = args)
  mk <- do.call(abind, c(mk, along = 3))
  m1 <- vector("list", ncol(mk))
  for (n in seq_along(m1)) {
    # summarize only over non-truncated values for this observation
    J <- (vec_lb[n] - min_lb + 1):(vec_ub[n] - min_lb)
    m1[[n]] <- rowSums(mk[, n, ][, J, drop = FALSE])
  }
  rm(mk)
  m1 <- do.call(cbind, m1)
  m1 / (do.call(cdf, c(list(ub), args)) - do.call(cdf, c(list(lb), args)))
}
