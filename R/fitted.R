fitted_response <- function(draws, mu) {
  # comnpute fitted values on the response scale
  # Args:
  #   draws: output of extract_draws
  #   mu: linear predictor matrix
  # Returns: 
  #   (usually) an S x N matrix containing samples of 
  #   the response distribution's mean 
  data <- draws$data
  is_trunc <- !(is.null(data$lb) && is.null(data$ub))
  # compute (mean) fitted values
  dim <- c(nrow(mu), draws$data$N)
  if (draws$f$family == "binomial") {
    if (length(data$max_obs) > 1L) {
      trials <- matrix(data$max_obs, nrow = dim[1], 
                       ncol = dim[2], byrow = TRUE)
    } else {
      trials <- data$max_obs
    }
    mu <- ilink(mu, draws$f$link) 
    if (!is_trunc) {
      # scale mu from [0,1] to [0,max_obs]
      mu <- mu * trials 
    }
  } else if (draws$f$family == "lognormal") {
    sigma <- get_sigma(draws$sigma, data = draws$data, dim = dim)
    mu <- ilink(mu, draws$f$link)
    if (!is_trunc) {
      # compute untruncated lognormal mean
      mu <- exp(mu + sigma^2 / 2)  
    }
  } else if (draws$f$family == "weibull") {
    shape <- get_shape(draws$shape, data = draws$data, dim = dim)
    mu <- ilink(mu / shape, draws$f$link)
    if (!is_trunc) {
      # compute untruncated weibull mean
      mu <- mu * gamma(1 + 1 / shape) 
    }
  } else if (is.ordinal(draws$f) || is.categorical(draws$f)) {
    mu <- fitted_catordinal(mu, max_obs = data$max_obs, family = draws$f)
  } else if (is.hurdle(draws$f)) {
    shape <- get_shape(draws$shape, data = draws$data, dim = dim)
    hu <- get_theta(draws, par = "hu")
    mu <- fitted_hurdle(mu, hu = hu, shape = shape, family = draws$f)
  } else if (is.zero_inflated(draws$f)) {
    zi <- get_theta(draws, par = "zi")
    mu <- fitted_zero_inflated(mu, zi = zi, family = draws$f)
    if (draws$f$family == "zero_inflated_binomial") {
      if (length(data$max_obs) > 1L) {
        if (!is.null(draws$data$N_trait)) {
          # deprecated as of brms 1.0.0
          J <- seq_len(ceiling(length(data$max_obs) / 2))
          trials <- data$max_obs[J]
        } else {
          trials <- data$max_obs
        }
        trials <- matrix(trials, nrow = dim[1], ncol = dim[2], byrow = TRUE)
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
      trunc_args <- nlist(mu, lb, ub, draws, dim)
      mu <- do.call(fitted_trunc_fun, trunc_args)
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

fitted_hurdle <- function(mu, hu, shape, family) {
  # Args:
  #   mu: linear predictor matrix
  #   hu: hurdle probability samples
  #   shape: shape parameter samples
  #   family: the model family
  if (isTRUE(ncol(mu) == 2 * ncol(hu))) {
    # old multivariate syntax model
    mu <- mu[, seq_len(ncol(mu) / 2)]
  }
  pre_mu <- ilink(mu, family$link)
  # adjust pre_mu as it is no longer the mean of the truncated distributions
  if (family$family == "hurdle_poisson") {
    adjusted_mu <- pre_mu / (1 - exp(-pre_mu))
  } else if (family$family == "hurdle_negbinomial") {
    adjusted_mu <- pre_mu / (1 - (shape / (pre_mu + shape))^shape)
  } else {
    adjusted_mu <- pre_mu
  }
  # incorporate hurdle part
  pre_mu * (1 - hu)
}

fitted_zero_inflated <- function(mu, zi, family) {
  # Args:
  #   mu: linear predictor matrix
  #   zi: zero-inflation probability samples
  #   family: the model family
  if (isTRUE(ncol(mu) == 2 * ncol(zi))) {
    # old multivariate syntax model
    mu <- mu[, seq_len(ncol(mu) / 2)]
  }
  ilink(mu, family$link) * (1 - zi) 
}

#------------- helper functions for truncated models ----------------
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