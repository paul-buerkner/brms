# All functions in this file have the same arguments structure
# Args:
#  i: the column of draws to use i.e. the ith obervation 
#     in the initial data.frame 
#  draws: A named list returned by extract_draws containing 
#         all required data and samples
#  ntrys: number of trys in rejection sampling for truncated discrete models
#  ...: ignored arguments
# Returns:
#   A vector of length draws$nsamples containing samples
#   from the posterior predictive distribution
predict_gaussian <- function(i, draws, ...) {
  sigma <- get_sigma(draws$sigma, data = draws$data, method = "predict", i = i)
  args <- list(mean = ilink(get_eta(i, draws), draws$f$link), sd = sigma)
  rng_continuous(nrng = draws$nsamples, dist = "norm",
                 args = args, data = draws$data)
}

predict_student <- function(i, draws, ...) {
  sigma <- get_sigma(draws$sigma, data = draws$data, method = "predict", i = i)
  args <- list(df = draws$nu, mu = ilink(get_eta(i, draws), draws$f$link), sigma = sigma)
  rng_continuous(nrng = draws$nsamples, dist = "student",
                 args = args, data = draws$data)
}

predict_cauchy <- function(i, draws, ...) {
  sigma <- get_sigma(draws$sigma, data = draws$data, method = "predict", i = i)
  args <- list(df = 1, mu = ilink(get_eta(i, draws), draws$f$link), sigma = sigma)
  rng_continuous(nrng = draws$nsamples, dist = "student",
                 args = args, data = draws$data)
}

predict_lognormal <- function(i, draws, ...) {
  sigma <- get_sigma(draws$sigma, data = draws$data, method = "predict", i = i)
  args <- list(meanlog = ilink(get_eta(i, draws), draws$f$link), sdlog = sigma)
  rng_continuous(nrng = draws$nsamples, dist = "lnorm",
                 args = args, data = draws$data)
}

predict_gaussian_multi <- function(i, draws, ...) {
  # currently no truncation available
  obs <- seq(i, draws$data$N, draws$data$N_trait)
  eta <- get_eta(obs, draws)
  .fun <- function(s) {
    rmulti_normal(1, Sigma = draws$Sigma[s, , ],
                  mu = ilink(eta[s, ], draws$f$link))
  }
  do.call(rbind, lapply(1:draws$nsamples, .fun))
}

predict_student_multi <- function(i, draws, ...) {
  # currently no truncation available
  obs <- seq(i, draws$data$N, draws$data$N_trait)
  eta <- get_eta(obs, draws)
  .fun <- function(s) {
    rmulti_student(1, df = draws$nu[s, ], 
                  mu = ilink(eta[s, ], draws$f$link),
                  Sigma = draws$Sigma[s, , ])
  }
  do.call(rbind, lapply(1:draws$nsamples, .fun))
}

predict_cauchy_multi <- function(i, draws, ...) {
  # currently no truncation available
  obs <- seq(i, draws$data$N, draws$data$N_trait)
  eta <- get_eta(obs, draws)
  .fun <- function(s) {
    rmulti_student(1, df = 1, mu = ilink(eta[s, ], draws$f$link),
                  Sigma = draws$Sigma[s, , ])
  }
  do.call(rbind, lapply(1:draws$nsamples, .fun))
}

predict_gaussian_cov <- function(i, draws, ...) {
  # currently, only ARMA1 processes are implemented
  obs <- with(draws$data, begin_tg[i]:(begin_tg[i] + nobs_tg[i] - 1))
  eta <- get_eta(obs, draws)
  args <- list(sigma = draws$sigma, se2 = draws$data$se2[obs], 
               nrows = length(obs))
  if (!is.null(draws$ar) && is.null(draws$ma)) {
    # AR1 process
    args$ar <- draws$ar
    Sigma <- do.call(get_cov_matrix_ar1, args)
  } else if (is.null(draws$ar) && !is.null(draws$ma)) {
    # MA1 process
    args$ma <- draws$ma
    Sigma <- do.call(get_cov_matrix_ma1, args)
  } else {
    # ARMA1 process
    args[c("ar", "ma")] <- draws[c("ar", "ma")]
    Sigma <- do.call(get_cov_matrix_arma1, args)
  }
  .fun <- function(s) {
    rmulti_normal(1, mu = ilink(eta[s, ], draws$f$link), 
                  Sigma = Sigma[s, , ])
  }
  do.call(rbind, lapply(1:draws$nsamples, .fun))
}

predict_student_cov <- function(i, draws, ...) {
  # currently, only ARMA1 processes are implemented
  obs <- with(draws$data, begin_tg[i]:(begin_tg[i] + nobs_tg[i] - 1))
  eta <- get_eta(obs, draws)
  args <- list(sigma = draws$sigma, se2 = draws$data$se2[obs], 
               nrows = length(obs))
  if (!is.null(draws$ar) && is.null(draws$ma)) {
    # AR1 process
    args$ar <- draws$ar
    Sigma <- do.call(get_cov_matrix_ar1, args)
  } else if (is.null(draws$ar) && !is.null(draws$ma)) {
    # MA1 process
    args$ma <- draws$ma
    Sigma <- do.call(get_cov_matrix_ma1, args)
  } else {
    # ARMA1 process
    args[c("ar", "ma")] <- draws[c("ar", "ma")]
    Sigma <- do.call(get_cov_matrix_arma1, args)
  }
  .fun <- function(s) {
    rmulti_student(1, df = draws$nu[s, ], 
                   mu = ilink(eta[s, ], draws$f$link), 
                   Sigma = Sigma[s, , ])
  }
  do.call(rbind, lapply(1:draws$nsamples, .fun))
}

predict_cauchy_cov <- function(i, draws, ...) {
  draws$nu <- matrix(rep(1, draws$nsamples))
  predict_student_cov(i = i, draws = draws, ...) 
}

predict_gaussian_fixed <- function(i, draws, ...) {
  stopifnot(i == 1)
  eta <- get_eta(1:nrow(draws$data$V), draws)
  .fun <- function(s) {
    rmulti_normal(1, mu = ilink(eta[s, ], draws$f$link), 
                  Sigma = draws$data$V)
  }
  do.call(rbind, lapply(1:draws$nsamples, .fun))
}

predict_student_fixed <- function(i, draws, ...) {
  stopifnot(i == 1)
  eta <- get_eta(1:nrow(draws$data$V), draws)
  .fun <- function(s) {
    rmulti_student(1, df = draws$nu[s, ], Sigma = draws$data$V,
                   mu = ilink(eta[s, ], draws$f$link))
  }
  do.call(rbind, lapply(1:draws$nsamples, .fun))
}

predict_cauchy_fixed <- function(i, draws, ...) {
  stopifnot(i == 1)
  draws$nu <- matrix(rep(1, draws$nsamples))
  predict_student_fixed(i, draws = draws, ...)
}

predict_binomial <- function(i, draws, ntrys = 5, ...) {
  trials <- ifelse(length(draws$data$max_obs) > 1, 
                   draws$data$max_obs[i], draws$data$max_obs) 
  args <- list(size = trials, prob = ilink(get_eta(i, draws), draws$f$link))
  rng_discrete(nrng = draws$nsamples, dist = "binom",
               args = args, data = draws$data, ntrys = ntrys)
}

predict_bernoulli <- function(i, draws, ...) {
  # truncation not useful
  if (!is.null(draws$data$N_trait)) {  # 2PL model
    eta <- get_eta(i, draws) * exp(get_eta(i + draws$data$N_trait, draws))
  } else {
    eta <- get_eta(i, draws)
  }
  rbinom(length(eta), size = 1, prob = ilink(eta, draws$f$link))
}

predict_poisson <- function(i, draws, ntrys = 5, ...) {
  args <- list(lambda = ilink(get_eta(i, draws), draws$f$link))
  rng_discrete(nrng = draws$nsamples, dist = "pois",
               args = args, data = draws$data, ntrys = ntrys)
}

predict_negbinomial <- function(i, draws, ntrys = 5, ...) {
  shape <- get_shape(draws$shape, data = draws$data, method = "predict", i = i)
  args <- list(mu = ilink(get_eta(i, draws), draws$f$link), size = shape)
  rng_discrete(nrng = draws$nsamples, dist = "nbinom",
               args = args, data = draws$data, ntrys = ntrys)
}

predict_geometric <- function(i, draws, ntrys = 5, ...) {
  args <- list(mu = ilink(get_eta(i, draws), draws$f$link), size = 1)
  rng_discrete(nrng = draws$nsamples, dist = "nbinom",
               args = args, data = draws$data, ntrys = ntrys)
}

predict_exponential <- function(i, draws, ...) {
  args <- list(rate = 1 / ilink(get_eta(i, draws), draws$f$link))
  rng_continuous(nrng = draws$nsamples, dist = "exp",
                 args = args, data = draws$data)
}

predict_gamma <- function(i, draws, ...) {
  shape <- get_shape(draws$shape, data = draws$data, method = "predict", i = i)
  args <- list(shape = shape, 
               scale = ilink(get_eta(i, draws), draws$f$link) / shape)
  rng_continuous(nrng = draws$nsamples, dist = "gamma",
                 args = args, data = draws$data)
}

predict_weibull <- function(i, draws, ...) {
  shape <- get_shape(draws$shape, data = draws$data, method = "predict", i = i)
  args <- list(shape = shape, 
               scale = ilink(get_eta(i, draws) / shape, draws$f$link))
  rng_continuous(nrng = draws$nsamples, dist = "weibull",
                 args = args, data = draws$data)
}

predict_inverse.gaussian <- function(i, draws, ...) {
  args <- list(mean = ilink(get_eta(i, draws), draws$f$link), 
               shape = draws$shape)
  rng_continuous(nrng = draws$nsamples, dist = "invgauss",
                 args = args, data = draws$data)
}

predict_beta <- function(i, draws, ...) {
  mu <- ilink(get_eta(i, draws), draws$f$link)
  args <- list(shape1 = mu * draws$phi, shape2 = (1 - mu) * draws$phi)
  rng_continuous(nrng = draws$nsamples, dist = "beta",
                 args = args, data = draws$data)
}

predict_hurdle_poisson <- function(i, draws, ...) {
  # theta is the bernoulii hurdle parameter
  theta <- ilink(get_eta(i + draws$data$N_trait, draws), "logit")
  lambda <- ilink(get_eta(i, draws), draws$f$link)
  ndraws <- draws$nsamples
  # compare with theta to incorporate the hurdle process
  hu <- runif(ndraws, 0, 1)
  # sample from a truncated poisson distribution
  # by adjusting lambda and adding 1
  t = -log(1 - runif(ndraws) * (1 - exp(-lambda)))
  ifelse(hu < theta, 0, rpois(ndraws, lambda = lambda - t) + 1)
}

predict_hurdle_negbinomial <- function(i, draws, ...) {
  # theta is the bernoulii hurdle parameter
  theta <- ilink(get_eta(i + draws$data$N_trait, draws), "logit")
  mu <- ilink(get_eta(i, draws), draws$f$link)
  ndraws <- draws$nsamples
  # compare with theta to incorporate the hurdle process
  hu <- runif(ndraws, 0, 1)
  # sample from an approximative(!) truncated negbinomial distribution
  # by adjusting mu and adding 1
  t = -log(1 - runif(ndraws) * (1 - exp(-mu)))
  ifelse(hu < theta, 0, rnbinom(ndraws, mu = mu - t, size = draws$shape) + 1)
}

predict_hurdle_gamma <- function(i, draws, ...) {
  # theta is the bernoulii hurdle parameter
  theta <- ilink(get_eta(i + draws$data$N_trait, draws), "logit")
  scale <- ilink(get_eta(i, draws), draws$f$link) / draws$shape
  ndraws <- draws$nsamples
  # compare with theta to incorporate the hurdle process
  hu <- runif(ndraws, 0, 1)
  ifelse(hu < theta, 0, rgamma(ndraws, shape = draws$shape, scale = scale))
}

predict_zero_inflated_beta <- function(i, draws, ...) {
  # theta is the bernoulii hurdle parameter
  theta <- ilink(get_eta(i + draws$data$N_trait, draws), "logit")
  mu <- ilink(get_eta(i, draws), draws$f$link)
  shape1 <- mu * draws$phi
  shape2 <- (1 - mu) * draws$phi
  ndraws <- draws$nsamples
  # compare with theta to incorporate the hurdle process
  hu <- runif(ndraws, 0, 1)
  ifelse(hu < theta, 0, rbeta(ndraws, shape1 = shape1, shape2 = shape2))
}

predict_zero_inflated_poisson <- function(i, draws, ...) {
  # theta is the bernoulii zero-inflation parameter
  theta <- ilink(get_eta(i + draws$data$N_trait, draws), "logit")
  lambda <- ilink(get_eta(i, draws), draws$f$link)
  ndraws <- draws$nsamples
  # compare with theta to incorporate the zero-inflation process
  zi <- runif(ndraws, 0, 1)
  ifelse(zi < theta, 0, rpois(ndraws, lambda = lambda))
}

predict_zero_inflated_negbinomial <- function(i, draws, ...) {
  # theta is the bernoulii zero-inflation parameter
  theta <- ilink(get_eta(i + draws$data$N_trait, draws), "logit")
  mu <- ilink(get_eta(i, draws), draws$f$link)
  ndraws <- draws$nsamples
  # compare with theta to incorporate the zero-inflation process
  zi <- runif(ndraws, 0, 1)
  ifelse(zi < theta, 0, rnbinom(ndraws, mu = mu, size = draws$shape))
}

predict_zero_inflated_binomial <- function(i, draws, ...) {
  # theta is the bernoulii zero-inflation parameter
  theta <- ilink(get_eta(i + draws$data$N_trait, draws), "logit")
  trials <- ifelse(length(draws$data$max_obs) > 1, 
                   draws$data$max_obs[i], draws$data$max_obs)
  prob <- ilink(get_eta(i, draws), draws$f$link)
  ndraws <- draws$nsamples
  # compare with theta to incorporate the zero-inflation process
  zi <- runif(ndraws, 0, 1)
  ifelse(zi < theta, 0, rbinom(ndraws, size = trials, prob = prob))
}

predict_categorical <- function(i, draws, ...) {
  ncat <- ifelse(length(draws$data$max_obs) > 1, 
                 draws$data$max_obs[i], draws$data$max_obs) 
  p <- pcategorical(1:ncat, eta = get_eta(i, draws, ordinal = TRUE)[, 1, ], 
                    ncat = ncat, link = draws$f$link)
  first_greater(p, target = runif(draws$nsamples, min = 0, max = 1))
}

predict_cumulative <- function(i, draws, ...) {
  predict_ordinal(i = i, draws = draws, family = "cumulative")
}

predict_sratio <- function(i, draws, ...) {
  predict_ordinal(i = i, draws = draws, family = "sratio")
}

predict_cratio <- function(i, draws, ...) {
  predict_ordinal(i = i, draws = draws, family = "cratio")
}

predict_acat <- function(i, draws, ...) {
  predict_ordinal(i = i, draws = draws, family = "acat")
}  

predict_ordinal <- function(i, draws, family, ...) {
  ncat <- ifelse(length(draws$data$max_obs) > 1, 
                 draws$data$max_obs[i], draws$data$max_obs)
  p <- pordinal(1:ncat, eta = get_eta(i, draws, ordinal = TRUE)[, 1, ], 
                ncat = ncat, family = family, link = draws$f$link)
  first_greater(p, target = runif(draws$nsamples, min = 0, max = 1))
}

#---------------predict helper-functions----------------------------

rng_continuous <- function(nrng, dist, args, data) {
  # random numbers from (possibly truncated) continuous distributions
  # Args:
  #   nrng: number of random values to generate
  #   dist: name of a distribution for which the functions
  #         p<dist>, q<dist>, and r<dist> are available
  #   args: dditional arguments passed to the distribution functions
  #   data: data initially passed to Stan
  # Returns:
  #   a vector of random values draws from the distribution
  if (is.null(data$lb) && is.null(data$ub)) {
    # sample as usual
    rdist <- paste0("r",dist)
    out <- do.call(rdist, c(nrng, args))
  } else {
    # sample from truncated distribution
    lb <- ifelse(is.null(data$lb), -Inf, data$lb)
    ub <- ifelse(is.null(data$ub), Inf, data$ub)
    pdist <- paste0("p",dist)
    qdist <- paste0("q",dist)
    plb <- do.call(pdist, c(lb, args))
    pub <- do.call(pdist, c(ub, args))
    rng <- list(runif(nrng, min = plb, max = pub))
    out <- do.call(qdist, c(rng, args))
    # remove infinte values caused by numerical imprecision
    out[out %in% c(-Inf, Inf)] <- NA
  }
  out
}

rng_discrete <- function(nrng, dist, args, data, ntrys) {
  # random numbers from (possibly truncated) discrete distributions
  # currently rejection sampling is used for truncated distributions
  # Args:
  #   nrng: number of random values to generate
  #   dist: name of a distribution for which the functions
  #         p<dist>, q<dist>, and r<dist> are available
  #   args: dditional arguments passed to the distribution functions
  #   data: data initially passed to Stan
  #   number of trys in rejection sampling for truncated models
  # Returns:
  #   a vector of random values draws from the distribution
  rdist <- get(paste0("r",dist), mode = "function")
  if (is.null(data$lb) && is.null(data$ub)) {
    # sample as usual
    do.call(rdist, c(nrng, args))
  } else {
    # sample from truncated distribution via rejection sampling
    lb <- ifelse(is.null(data$lb), -Inf, data$lb)
    ub <- ifelse(is.null(data$ub), Inf, data$ub)
    rng <- matrix(do.call(rdist, c(nrng * ntrys, args)), ncol = ntrys)
    apply(rng, 1, extract_valid_sample, lb = lb, ub = ub)
  }
}

extract_valid_sample <- function(rng, lb, ub) {
  # extract the first valid predicted value 
  # per Stan sample per observation 
  # Args:
  #   rng: draws to be check against truncation boundaries
  #   lb: lower bound
  #   ub: upper bound
  # Returns:
  #   a valid truncated sample or else the closest boundary
  valid_rng <- match(TRUE, rng > lb & rng <= ub)
  if (is.na(valid_rng)) {
    # no valid truncated value found
    # set sample to lb or ub
    # 1e-10 is only to identify the invalid draws later on
    ifelse(max(rng) <= lb, lb + 1 - 1e-10, ub + 1e-10)
  } else {
    rng[valid_rng]
  }
}

get_pct_invalid <- function(x, data) {
  # percentage of invalid predictions of truncated discrete models
  # Args:
  #   x: matrix of predicted values
  #   data: data initially passed to Stan
  if (!(is.null(data$lb) && is.null(data$ub))) {
    lb <- ifelse(is.null(data$lb), -Inf, data$lb)
    ub <- ifelse(is.null(data$ub), Inf, data$ub)
    x <- c(x)[!is.na(c(x))]
    sum(x < (lb + 1) | x > ub) / length(x) 
  } else {
    0
  }
}