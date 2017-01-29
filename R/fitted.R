# functions in this file have the same arguments structure
# Args:
#   mu: untransformed linear predictor matrix
#   draws: A named list returned by extract_draws containing 
#          all required data and samples
# Returns:
#   transformed linear predictor representing the mean
#   of the response distribution
fitted_gaussian <- function(mu, draws) {
  fitted_default(mu, draws)
}

fitted_student <- function(mu, draws) {
  fitted_default(mu, draws)
}

fitted_cauchy <- function(mu, draws) {
  fitted_default(mu, draws)
}

fitted_lognormal <- function(mu, draws) {
  dim <- dim_mu(mu, draws)
  sigma <- get_sigma(draws$sigma, data = draws$data, dim = dim)
  mu <- ilink(mu, draws$f$link)
  if (!is_trunc(draws$data)) {
    mu <- exp(mu + sigma^2 / 2)  
  } else {
    mu <- fitted_trunc(mu, draws)
  }
  mu
}

fitted_binomial <- function(mu, draws) {
  dim <- dim_mu(mu, draws)
  trials <- matrix(draws$data$trials, nrow = dim[1], 
                   ncol = dim[2], byrow = TRUE)
  mu <- ilink(mu, draws$f$link) 
  if (!is_trunc(draws$data)) {
    # scale mu from [0,1] to [0,trials]
    mu <- mu * trials 
  } else {
    mu <- fitted_trunc(mu, draws)
  }
  mu
}

fitted_bernoulli <- function(mu, draws) {
  fitted_default(mu, draws)
}

fitted_poisson <- function(mu, draws) {
  fitted_default(mu, draws)
}

fitted_negbinomial <- function(mu, draws) {
  fitted_default(mu, draws)
}

fitted_geometric <- function(mu, draws) {
  fitted_default(mu, draws)
}

fitted_exponential <- function(mu, draws) {
  fitted_default(mu, draws)
}

fitted_gamma <- function(mu, draws) {
  fitted_default(mu, draws)
}

fitted_weibull <- function(mu, draws) {
  dim <- dim_mu(mu, draws)
  shape <- get_shape(draws$shape, data = draws$data, dim = dim)
  mu <- ilink(mu / shape, draws$f$link)
  if (!is_trunc(draws$data)) {
    mu <- mu * gamma(1 + 1 / shape) 
  } else {
    mu <- fitted_trunc(mu, draws)
  }
  mu
}

fitted_frechet <- function(mu, draws) {
  fitted_default(mu, draws)
}

fitted_gen_extreme_value <- function(mu, draws) {
  dim <- dim_mu(mu, draws)
  sigma <- get_sigma(draws$sigma, data = draws$data, dim = dim)
  xi <- get_auxpar(draws$xi)
  mu <- ilink(mu, draws$f$link)
  if (!is_trunc(draws$data)) {
    mu <- mu + sigma * (gamma(1 - xi) - 1) / xi
  } else {
    mu <- fitted_trunc(mu, draws)
  }
  mu
}

fitted_inverse.gaussian <- function(mu, draws) {
  fitted_default(mu, draws)
}

fitted_exgaussian <- function(mu, draws) {
  mu <- ilink(mu, draws$f$link) + get_auxpar(draws$beta)
  fitted_trunc(mu, draws)
}

fitted_wiener <- function(mu, draws) {
  delta <- ilink(mu, draws$f$link)
  bs <- get_auxpar(draws$bs)
  ndt <- get_auxpar(draws$ndt)
  bias <- get_auxpar(draws$bias)
  ndt - bias / delta + bs / delta * 
    (exp(- 2 * delta * bias) - 1) / (exp(-2 * delta * bs) - 1)
}

fitted_beta <- function(mu, draws) {
  fitted_default(mu, draws)
}

fitted_von_mises <- function(mu, draws) {
  fitted_default(mu, draws)
}

fitted_asym_laplace <- function(mu, draws) {
  dim <- dim_mu(mu, draws)
  mu <- ilink(mu, draws$f$link)
  quantile <- get_auxpar(draws$quantile)
  sigma <- get_sigma(draws$sigma, data = draws$data, dim = dim)
  mu + sigma * (1 - 2 * quantile) / (quantile * (1 - quantile))
}

fitted_hurdle_poisson <- function(mu, draws) {
  hu <- get_theta(draws, par = "hu")
  mu <- adjust_old_forked(mu, hu)
  mu <- ilink(mu, draws$f$link)
  mu / (1 - exp(-mu)) * (1 - hu)
}

fitted_hurdle_negbinomial <- function(mu, draws) {
  dim <- dim_mu(mu, draws)
  shape <- get_shape(draws$shape, data = draws$data, dim = dim)
  hu <- get_theta(draws, par = "hu")
  mu <- adjust_old_forked(mu, hu)
  mu <- ilink(mu, draws$f$link)
  mu / (1 - (shape / (mu + shape))^shape) * (1 - hu)
}

fitted_hurdle_gamma <- function(mu, draws) {
  hu <- get_theta(draws, par = "hu")
  mu <- adjust_old_forked(mu, hu)
  mu <- ilink(mu, draws$f$link)
  mu * (1 - hu)
}

fitted_hurdle_lognormal <- function(mu, draws) {
  dim <- dim_mu(mu, draws)
  sigma <- get_sigma(draws$sigma, data = draws$data, dim = dim)
  hu <- get_theta(draws, par = "hu")
  mu <- ilink(mu, draws$f$link)
  exp(mu + sigma^2 / 2) * (1 - hu)
}

fitted_zero_inflated_poisson <- function(mu, draws) {
  zi <- get_theta(draws, par = "zi")
  mu <- adjust_old_forked(mu, zi)
  ilink(mu, draws$f$link) * (1 - zi) 
}

fitted_zero_inflated_negbinomial <- function(mu, draws) {
  zi <- get_theta(draws, par = "zi")
  mu <- adjust_old_forked(mu, zi)
  ilink(mu, draws$f$link) * (1 - zi) 
}

fitted_zero_inflated_binomial <- function(mu, draws) {
  zi <- get_theta(draws, par = "zi")
  mu <- adjust_old_forked(mu, zi)
  mu <- ilink(mu, draws$f$link) * (1 - zi)
  trials <- draws$data[["trials"]]
  if (!is.null(draws$data$N_trait)) {
    # deprecated as of brms 1.0.0
    J <- seq_len(ceiling(length(trials) / 2))
    trials <- trials[J]
  }
  dim <- dim_mu(mu, draws)
  trials <- matrix(trials, nrow = dim[1], ncol = dim[2], byrow = TRUE)
  mu * trials
}

fitted_zero_inflated_beta <- function(mu, draws) {
  zi <- get_theta(draws, par = "zi")
  mu <- adjust_old_forked(mu, zi)
  ilink(mu, draws$f$link) * (1 - zi) 
}

fitted_categorical <- function(mu, draws) {
  fitted_catordinal(mu, draws)
}

fitted_cumulative <- function(mu, draws) {
  disc <- get_disc(draws, ncat = draws$data[["ncat"]])
  mu <- disc * mu 
  fitted_catordinal(mu, draws)
}

fitted_sratio <- function(mu, draws) {
  disc <- get_disc(draws, ncat = draws$data[["ncat"]])
  mu <- disc * mu 
  fitted_catordinal(mu, draws)
}

fitted_cratio <- function(mu, draws) {
  disc <- get_disc(draws, ncat = draws$data[["ncat"]])
  mu <- disc * mu 
  fitted_catordinal(mu, draws)
}

fitted_acat <- function(mu, draws) {
  disc <- get_disc(draws, ncat = draws$data[["ncat"]])
  mu <- disc * mu 
  fitted_catordinal(mu, draws)
}

# ------ fitted helper functions ------

fitted_default <- function(mu, draws) {
  # default fitted values
  mu <- ilink(mu, draws$f$link)
  fitted_trunc(mu, draws)
}

fitted_catordinal <- function(mu, draws) {
  # fitted values for categorical and ordinal models
  get_density <- function(s) {
    # get probabilities of each category
    do.call(dens, c(args, list(eta = mu[, s, ])))
  }
  ncat <- draws$data[["ncat"]]
  args <- list(seq_len(ncat), ncat = ncat, link = draws$f$link)
  dens <- paste0("d", draws$f$family)
  out <- abind(lapply(seq_len(ncol(mu)), get_density), along = 3)
  aperm(out, perm = c(1, 3, 2))
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

dim_mu <- function(mu, draws) {
  c(nrow(mu), draws$data$N)
}

is_trunc <- function(data) {
  any(data[["lb"]] > - Inf) || any(data[["ub"]] < Inf)
}

fitted_trunc <- function(mu, draws) {
  # prepares data required for truncation and calles the 
  # family specific truncation function for fitted values
  if (is_trunc(draws$data)) {
    lb <- matrix(draws$data[["lb"]], nrow = nrow(mu), 
                 ncol = ncol(mu), byrow = TRUE)
    ub <- matrix(draws$data[["ub"]], nrow = nrow(mu), 
                 ncol = ncol(mu), byrow = TRUE)
    fitted_trunc_fun <- paste0("fitted_trunc_", draws$f$family)
    fitted_trunc_fun <- try(get(fitted_trunc_fun, asNamespace("brms")),
                            silent = TRUE)
    if (is(fitted_trunc_fun, "try-error")) {
      stop2("Fitted values on the respone scale not yet implemented ",
            "for truncated '", draws$f$family, "' models.")
    } else {
      dim <- dim_mu(mu, draws)
      trunc_args <- nlist(mu, lb, ub, draws, dim)
      mu <- do.call(fitted_trunc_fun, trunc_args)
    }
  }
  mu
}

# ----- family specific truncation functions -----
# Args:
#   mu: (usually) samples of the untruncated mean parameter
#   lb: lower truncation bound
#   ub: upper truncation bound
#   draws: output of extract_draws
#   dim: vector of length 2: c(Nsamples, Nobs)
#   ...: ignored arguments
# Returns:
#   samples of the truncated mean parameter
fitted_trunc_gaussian <- function(mu, lb, ub, draws, dim) {
  sigma <- get_sigma(draws$sigma, data = draws$data, dim = dim)
  zlb <- (lb - mu) / sigma
  zub <- (ub - mu) / sigma
  # truncated mean of standard normal; see Wikipedia
  trunc_zmean <- (dnorm(zlb) - dnorm(zub)) / (pnorm(zub) - pnorm(zlb))  
  mu + trunc_zmean * sigma  
}

fitted_trunc_student <- function(mu, lb, ub, draws, dim) {
  sigma <- get_sigma(draws$sigma, data = draws$data, dim = dim)
  nu <- get_auxpar(draws$nu)
  zlb <- (lb - mu) / sigma
  zub <- (ub - mu) / sigma
  # see Kim 2008: Moments of truncated Student-t distribution
  G1 <- gamma((nu - 1) / 2) * nu^(nu / 2) / 
    (2 * (pt(zub, df = nu) - pt(zlb, df = nu)) * gamma(nu / 2) * gamma(0.5))
  A <- (nu + zlb^2) ^ (-(nu - 1) / 2)
  B <- (nu + zub^2) ^ (-(nu - 1) / 2)
  trunc_zmean <- G1 * (A - B)
  mu + trunc_zmean * sigma 
}

fitted_trunc_lognormal <- function(mu, lb, ub, draws, dim) {
  # mu has to be on the linear scale
  sigma <- get_sigma(draws$sigma, data = draws$data, dim = dim)
  m1 <- exp(mu + sigma^2 / 2) * (pnorm((log(ub) - mu) / sigma - sigma) - 
                                   pnorm((log(lb) - mu) / sigma - sigma))
  m1 / (plnorm(ub, meanlog = mu, sdlog = sigma) - 
          plnorm(lb, meanlog = mu, sdlog = sigma))
}

fitted_trunc_gamma <- function(mu, lb, ub, draws, dim) {
  shape <- get_shape(draws$shape, data = draws$data, dim = dim)
  scale <- mu / shape
  # see Jawitz 2004: Moments of truncated continuous univariate distributions
  m1 <- scale / gamma(shape) * 
    (incgamma(ub / scale, 1 + shape) - incgamma(lb / scale, 1 + shape))
  m1 / (pgamma(ub, shape, scale = scale) - pgamma(lb, shape, scale = scale))
}

fitted_trunc_exponential <- function(mu, lb, ub, ...) {
  # see Jawitz 2004: Moments of truncated continuous univariate distributions
  # mu is already the scale parameter
  inv_mu <- 1 / mu
  m1 <- mu * (incgamma(ub / mu, 2) - incgamma(lb / mu, 2))
  m1 / (pexp(ub, rate = inv_mu) - pexp(lb, rate = inv_mu))
}

fitted_trunc_weibull <- function(mu, lb, ub, draws, dim) {
  # see Jawitz 2004: Moments of truncated continuous univariate distributions
  # mu is already the scale parameter
  shape <- get_shape(draws$shape, data = draws$data, dim = dim)
  a <- 1 + 1 / shape
  m1 <- mu * (incgamma((ub / mu)^shape, a) - incgamma((lb / mu)^shape, a))
  m1 / (pweibull(ub, shape, scale = mu) - pweibull(lb, shape, scale = mu))
}

fitted_trunc_binomial <- function(mu, lb, ub, draws, dim) {
  lb <- ifelse(lb < -1, -1, lb)
  max_value <- max(draws$data$trials)
  ub <- ifelse(ub > max_value, max_value, ub)
  trials <- draws$data$trials
  if (length(trials) > 1) {
    trials <- matrix(trials, nrow = dim[1], ncol = dim[2], byrow = TRUE)
  }
  args <- list(size = trials, prob = mu)
  message(paste("Computing fitted values for a truncated binomial model.",
                "This may take a while."))
  fitted_trunc_discrete(dist = "binom", args = args, lb = lb, ub = ub)
}

fitted_trunc_poisson <- function(mu, lb, ub, draws, ...) {
  lb <- ifelse(lb < -1, -1, lb)
  max_value <- 3 * max(draws$data$Y)
  ub <- ifelse(ub > max_value, max_value, ub)
  args <- list(lambda = mu)
  message(paste("Computing fitted values for a truncated poisson model.",
                "This may take a while."))
  fitted_trunc_discrete(dist = "pois", args = args, lb = lb, ub = ub)
}

fitted_trunc_negbinomial <- function(mu, lb, ub, draws, dim) {
  lb <- ifelse(lb < -1, -1, lb)
  max_value <- 3 * max(draws$data$Y)
  ub <- ifelse(ub > max_value, max_value, ub)
  shape <- get_shape(draws$shape, data = draws$data, dim = dim)
  args <- list(mu = mu, size = shape)
  message(paste("Computing fitted values for a truncated negbinomial model.",
                "This may take a while."))
  fitted_trunc_discrete(dist = "nbinom", args = args, lb = lb, ub = ub)
}

fitted_trunc_geometric <- function(mu, lb, ub, draws, ...) {
  lb <- ifelse(lb < -1, -1, lb)
  max_value <- 3 * max(draws$data$Y)
  ub <- ifelse(ub > max_value, max_value, ub)
  args <- list(mu = mu, size = 1)
  message(paste("Computing fitted values for a truncated geometric model.",
                "This may take a while."))
  fitted_trunc_discrete(dist = "nbinom", args = args, lb = lb, ub = ub)
}

fitted_trunc_discrete <- function(dist, args, lb, ub) {
  stopifnot(is.matrix(lb), is.matrix(ub))
  pdf <- get(paste0("d", dist), mode = "function")
  cdf <- get(paste0("p", dist), mode = "function")
  get_mean_kernel <- function(x, args) {
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
  mk <- lapply((min_lb + 1):max(vec_ub), get_mean_kernel, args = args)
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
