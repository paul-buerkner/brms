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
    sigma <- get_sigma(draws$sigma, data = draws$data, 
                       i = nrow(mu), sobs = FALSE)
    mu <- ilink(mu, draws$f$link)
    if (!is_trunc) {
      # compute untruncated lognormal mean
      mu <- exp(mu + sigma^2 / 2)  
    }
  } else if (draws$f$family == "weibull") {
    shape <- get_shape(draws$shape, data = draws$data,
                       i = nrow(mu), sobs = FALSE)
    mu <- ilink(mu / shape, draws$f$link)
    if (!is_trunc) {
      # compute untruncated weibull mean
      mu <- mu * gamma(1 + 1 / shape) 
    }
  } else if (is.ordinal(draws$f) || is.categorical(draws$f)) {
    mu <- fitted_catordinal(mu, max_obs = data$max_obs, family = draws$f)
  } else if (is.hurdle(draws$f)) {
    shape <- get_shape(draws$shape, data = draws$data,
                       i = nrow(mu), sobs = FALSE)
    mu <- fitted_hurdle(mu, shape = shape, N_trait = data$N_trait,
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
    lb <- matrix(data$lb, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
    ub <- matrix(data$ub, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
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
  sigma <- get_sigma(draws$sigma, data = draws$data, 
                     i = nrow(mu), sobs = FALSE)
  zlb <- (lb - mu) / sigma
  zub <- (ub - mu) / sigma
  # truncated mean of standard normal; see Wikipedia
  trunc_zmean <- (dnorm(zlb) - dnorm(zub)) / (pnorm(zub) - pnorm(zlb))  
  mu + trunc_zmean * sigma  
}

fitted_trunc_student <- function(mu, lb, ub, draws) {
  sigma <- get_sigma(draws$sigma, data = draws$data, 
                     i = nrow(mu), sobs = FALSE)
  nu <- get_auxpar(draws$nu, i = nrow(mu), sobs = FALSE)
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
  sigma <- get_sigma(draws$sigma, data = draws$data, 
                     i = nrow(mu), sobs = FALSE)
  m1 <- exp(mu + sigma^2 / 2) * (pnorm((log(ub) - mu) / sigma - sigma) - 
                                   pnorm((log(lb) - mu) / sigma - sigma))
  m1 / (plnorm(ub, meanlog = mu, sdlog = sigma) - 
          plnorm(lb, meanlog = mu, sdlog = sigma))
}

fitted_trunc_gamma <- function(mu, lb, ub, draws) {
  shape <- get_shape(draws$shape, data = draws$data, 
                     i = nrow(mu), sobs = FALSE)
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
  shape <- get_shape(draws$shape, data = draws$data,
                     i = nrow(mu), sobs = FALSE)
  a <- 1 + 1 / shape
  m1 <- mu * (incgamma((ub / mu)^shape, a) - incgamma((lb / mu)^shape, a))
  m1 / (pweibull(ub, shape, scale = mu) - pweibull(lb, shape, scale = mu))
}

fitted_trunc_binomial <- function(mu, lb, ub, draws) {
  lb <- ifelse(lb < -1, -1, lb)
  max_value <- max(draws$data$trials)
  ub <- ifelse(ub > max_value, max_value, ub)
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
  lb <- ifelse(lb < -1, -1, lb)
  max_value <- 3 * max(draws$data$Y)
  ub <- ifelse(ub > max_value, max_value, ub)
  args <- list(lambda = mu)
  message(paste("Computing fitted values for a truncated poisson model.",
                "This may take a while."))
  fitted_trunc_discrete(dist = "pois", args = args, lb = lb, ub = ub)
}

fitted_trunc_negbinomial <- function(mu, lb, ub, draws) {
  lb <- ifelse(lb < -1, -1, lb)
  max_value <- 3 * max(draws$data$Y)
  ub <- ifelse(ub > max_value, max_value, ub)
  shape <- get_shape(draws$shape, data = draws$data,
                     i = nrow(mu), sobs = FALSE)
  args <- list(mu = mu, size = shape)
  message(paste("Computing fitted values for a truncated negbinomial model.",
                "This may take a while."))
  fitted_trunc_discrete(dist = "nbinom", args = args, lb = lb, ub = ub)
}

fitted_trunc_geometric <- function(mu, lb, ub, draws) {
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