fitted_response <- function(draws, mu) {
  # comnpute fitted values on the response scale
  # Args:
  #   draws: output of extract_draws
  #   mu: linear predictor matrix
  # Returns: 
  #   (usually) an S x N matrix containing samples of 
  #   the response distribution's mean 
  data <- draws$data
  if (is.2PL(draws$f)) {
    # the second part of mu is the log discriminability
    mu <- mu[, 1:data$N_trait] * exp(mu[, (data$N_trait + 1):data$N])
  }
  is_trunc <- !(is.null(data$lb) && is.null(data$ub))
  
  # compute (mean) fitted values
  if (draws$f$family == "binomial") {
    if (length(data$max_obs) > 1L) {
      trials <- matrix(rep(data$max_obs, nrow(mu)), 
                       nrow = nrow(mu), byrow = TRUE)
    } else {
      trials <- data$max_obs
    }
    mu <- ilink(mu, draws$f$link) 
    if (!is_trunc) {
      # scale mu from [0,1] to [0,max_obs]
      mu <- mu * trials 
    }
  } else if (draws$f$family == "lognormal") {
    sigma <- get_sigma(draws$sigma, data = draws$data, i = nrow(mu))
    mu <- ilink(mu, draws$f$link)
    if (!is_trunc) {
      # compute untruncated lognormal mean
      mu <- exp(mu + sigma^2 / 2)  
    }
  } else if (draws$f$family == "weibull") {
    shape <- get_shape(draws$shape, data = draws$data)
    mu <- ilink(mu / shape, draws$f$link)
    if (!is_trunc) {
      # compute untruncated weibull mean
      mu <- mu * gamma(1 + 1 / shape) 
    }
  } else if (is.ordinal(draws$f) || is.categorical(draws$f)) {
    mu <- fitted_catordinal(mu, max_obs = data$max_obs, family = draws$f)
  } else if (is.hurdle(draws$f)) {
    mu <- fitted_hurdle(mu, shape = draws$shape, N_trait = data$N_trait,
                        family = draws$f)
  } else if (is.zero_inflated(draws$f)) {
    mu <- fitted_zero_inflated(mu, N_trait = data$N_trait, family = draws$f)
    if (draws$f$family == "zero_inflated_binomial") {
      if (length(data$max_obs) > 1L) {
        trials <- data$max_obs[1:ceiling(length(data$max_obs) / 2)]
        trials <- matrix(rep(trials, nrow(mu)), nrow = nrow(mu), 
                         byrow = TRUE)
      } else {
        trials <- data$max_obs
      }
      mu <- mu * trials
    }
  } else {
    mu <- ilink(mu, draws$f$link)
  }
  # fitted values for truncated models
  if (is_trunc) {
    lb <- ifelse(is.null(data$lb), -Inf, data$lb)
    ub <- ifelse(is.null(data$ub), Inf, data$ub)
    fitted_trunc_fun <- try(get(paste0("fitted_trunc_", draws$f$family), 
                                mode = "function"))
    if (is(fitted_trunc_fun, "try-error")) {
      stop(paste("fitted values on the respone scale not implemented",
                 "for truncated", family, "models"))
    } else {
      mu <- fitted_trunc_fun(mu = mu, lb = lb, ub = ub, draws = draws)
    }
  } 
  mu
}

fitted_catordinal <- function(mu, max_obs, family) {
  # compute fitted values for categorical and ordinal families
  ncat <- max(max_obs)
  # get probabilities of each category
  get_density <- function(s) {
    do.call(paste0("d", family$family), 
            list(1:ncat, eta = mu[, s, ], ncat = ncat, link = family$link))
  }
  aperm(abind(lapply(1:ncol(mu), get_density), along = 3), perm = c(1, 3, 2))
}

fitted_hurdle <- function(mu, shape, N_trait, family) {
  # Args:
  #   mu: linear predictor matrix
  #   shape: shape parameter samples
  #   N_trait: Number of observations of the main response variable
  n_base <- 1:N_trait
  n_hu <- n_base + N_trait
  pre_mu <- ilink(mu[, n_base], family$link)
  # adjust pre_mu as it is no longer the mean of the truncated distributions
  if (family$family == "hurdle_poisson") {
    adjusted_mu <- pre_mu / (1 - exp(-pre_mu))
  } else if (family$family == "hurdle_negbinomial") {
    adjusted_mu <- pre_mu / (1 - (shape / (pre_mu + shape))^shape)
  } else {
    adjusted_mu <- pre_mu
  }
  # incorporate hurdle part
  pre_mu * (1 - ilink(mu[, n_hu], "logit")) 
}

fitted_zero_inflated <- function(mu, N_trait, family) {
  # Args:
  #   mu: linear predictor matrix
  #   N_trait: Number of observations of the main response variable
  n_base <- 1:N_trait
  n_zi <- n_base + N_trait
  # incorporate zero-inflation part
  ilink(mu[, n_base], family$link) * (1 - ilink(mu[, n_zi], "logit")) 
}

#------------- helper functions for truncated models ----------------
# Args:
#   mu: (usually) samples of the untruncated mean parameter
#   lb: lower truncation bound
#   ub: upper truncation bound
#   draws: output of extract_draws
#   ...: ignored arguments
# Returns:
#   samples of the truncated mean parameter
fitted_trunc_gaussian <- function(mu, lb, ub, draws) {
  sigma <- get_sigma(draws$sigma, data = draws$data, i = nrow(mu))
  zlb <- (lb - mu) / sigma
  zub <- (ub - mu) / sigma
  # truncated mean of standard normal; see Wikipedia
  trunc_zmean <- (dnorm(zlb) - dnorm(zub)) / (pnorm(zub) - pnorm(zlb))  
  mu + trunc_zmean * sigma  
}

fitted_trunc_student <- function(mu, lb, ub, draws) {
  sigma <- get_sigma(draws$sigma, data = draws$data, i = nrow(mu))
  nu <- draws$nu
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

fitted_trunc_lognormal <- function(mu, lb, ub, draws) {
  # mu has to be on the linear scale
  sigma <- get_sigma(draws$sigma, data = draws$data, i = nrow(mu))
  m1 <- exp(mu + sigma^2 / 2) * (pnorm((log(ub) - mu) / sigma - sigma) - 
                                   pnorm((log(lb) - mu) / sigma - sigma))
  m1 / (plnorm(ub, meanlog = mu, sdlog = sigma) - 
          plnorm(lb, meanlog = mu, sdlog = sigma))
}

fitted_trunc_gamma <- function(mu, lb, ub, draws) {
  shape <- get_shape(draws$shape, data = draws$data)
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

fitted_trunc_weibull <- function(mu, lb, ub, draws) {
  # see Jawitz 2004: Moments of truncated continuous univariate distributions
  # mu is already the scale parameter
  shape <- get_shape(draws$shape, data = draws$data)
  a <- 1 + 1 / shape
  m1 <- mu * (incgamma((ub / mu)^shape, a) - incgamma((lb / mu)^shape, a))
  m1 / (pweibull(ub, shape, scale = mu) - pweibull(lb, shape, scale = mu))
}

fitted_trunc_binomial <- function(mu, lb, ub, draws) {
  lb <- max(lb, -1)
  ub <- min(ub, max(draws$data$trials))
  trials <- draws$data$trials
  if (length(trials) > 1) {
    trials <- matrix(rep(trials, nrow(mu)), ncol = draws$data$N, byrow = TRUE)
  }
  args <- list(size = trials, prob = mu)
  message(paste("Computing fitted values for a truncated binomial model.",
                "This may take a while."))
  fitted_trunc_discrete(dist = "binom", args = args, lb = lb, ub = ub)
}

fitted_trunc_poisson <- function(mu, lb, ub, draws) {
  lb <- max(lb, -1)
  ub <- min(ub, 3 * max(draws$data$Y))
  args <- list(lambda = mu)
  message(paste("Computing fitted values for a truncated poisson model.",
                "This may take a while."))
  fitted_trunc_discrete(dist = "pois", args = args, lb = lb, ub = ub)
}

fitted_trunc_negbinomial <- function(mu, lb, ub, draws) {
  lb <- max(lb, -1)
  ub <- min(ub, 3 * max(draws$data$Y))
  shape <- get_shape(draws$shape, data = draws$data)
  args <- list(mu = mu, size = shape)
  message(paste("Computing fitted values for a truncated negbinomial model.",
                "This may take a while."))
  fitted_trunc_discrete(dist = "nbinom", args = args, lb = lb, ub = ub)
}

fitted_trunc_geometric <- function(mu, lb, ub, draws) {
  lb <- max(lb, -1)
  ub <- min(ub, 3 * max(draws$data$Y))
  args <- list(mu = mu, size = 1)
  message(paste("Computing fitted values for a truncated geometric model.",
                "This may take a while."))
  fitted_trunc_discrete(dist = "nbinom", args = args, lb = lb, ub = ub)
}

fitted_trunc_discrete <- function(dist, args, lb, ub) {
  pdf <- get(paste0("d", dist), mode = "function")
  cdf <- get(paste0("p", dist), mode = "function")
  get_mean_kernel <- function(x, args) {
    # just x * density(x)
    x * do.call(pdf, c(x, args))
  }
  if (any(is.infinite(c(lb, ub)))) {
    stop("lb and ub must be finite")
  }
  # array of dimension S x N x length(lg:ub)
  mean_kernel <- lapply((lb + 1):ub, get_mean_kernel, args = args)
  mean_kernel <- do.call(abind, c(mean_kernel, along = 3))
  m1 <- apply(mean_kernel, MARGIN = 2, FUN = rowSums)
  m1 / (do.call(cdf, c(ub, args)) - do.call(cdf, c(lb, args)))
}