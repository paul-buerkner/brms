# All functions in this file have the same arguments structure
#
# Args:
#  n: the column of samples to use i.e. the nth obervation in the initial data.frame 
#  data: the data as passed to Stan
#  samples: samples obtained through Stan. Must at least contain variable eta
#  link: the link function
#  ntrys: number of trys in rejection sampling for truncated discrete models
#  ...: ignored arguments
#
# Returns:
#   A vector of length nrow(samples) containing samples from the posterior predictive distribution
predict_gaussian <- function(n, data, samples, link, ...) {
  sigma <- get_sigma(samples$sigma, data = data, method = "predict", n = n)
  args <- list(mean = ilink(samples$eta[, n], link), sd = sigma)
  rng_continuous(nrng = nrow(samples$eta), dist = "norm",
                 args = args, data = data)
}

predict_student <- function(n, data, samples, link, ...) {
  sigma <- get_sigma(samples$sigma, data = data, method = "predict", n = n)
  args <- list(df = samples$nu, mu = ilink(samples$eta[, n], link), sigma = sigma)
  rng_continuous(nrng = nrow(samples$eta), dist = "student",
                 args = args, data = data)
}

predict_cauchy <- function(n, data, samples, link, ...) {
  sigma <- get_sigma(samples$sigma, data = data, method = "predict", n = n)
  args <- list(df = 1, mu = ilink(samples$eta[, n], link), sigma = sigma)
  rng_continuous(nrng = nrow(samples$eta), dist = "student",
                 args = args, data = data)
}

predict_lognormal <- function(n, data, samples, link, ...) {
  # link is currently ignored for lognormal models
  # as 'identity' is the only valid link
  sigma <- get_sigma(samples$sigma, data = data, method = "predict", n = n)
  args <- list(meanlog = samples$eta[, n], sdlog = sigma)
  rng_continuous(nrng = nrow(samples$eta), dist = "lnorm",
                 args = args, data = data)
}

predict_multinormal <- function(n, data, samples, link, ...) {
  # currently no truncation available
  obs <- seq(n, data$N, data$N_trait)
  .fun <- function(i) {
    rmultinormal(1, Sigma = samples$Sigma[i, , ],
                 mu = ilink(samples$eta[i, obs], link))
  }
  do.call(rbind, lapply(1:nrow(samples$eta), .fun))
}

predict_gaussian_ar1 <- function(n, data, samples, link, ...) {
  # weights, truncation and censoring not allowed
  rows <- with(data, begin_tg[n]:(begin_tg[n] + nrows_tg[n] - 1))
  eta_part <- samples$eta[, rows]
  squared_se_part <- data$squared_se[rows]
  # both sigma and SEs are present!
  Sigma <- get_cov_matrix_ar1(ar = samples$ar, sigma = samples$sigma, 
                              sq_se = squared_se_part, nrows = length(rows)) 

  .fun <- function(i) {
    rmultinormal(1, mu = ilink(eta_part[i, ], link), 
                 Sigma = Sigma[i, , ])
  }
  do.call(rbind, lapply(1:nrow(samples$eta), .fun))
}

predict_binomial <- function(n, data, samples, link, ntrys, ...) {
  max_obs <- ifelse(length(data$max_obs) > 1, data$max_obs[n], data$max_obs) 
  args <- list(size = max_obs, prob = ilink(samples$eta[, n], link))
  rng_discrete(nrng = nrow(samples$eta), dist = "binom",
               args = args, data = data, ntrys = ntrys)
}  

predict_bernoulli <- function(n, data, samples, link, ...) {
  # truncation not useful
  rbinom(nrow(samples$eta), size = 1, 
         prob = ilink(samples$eta[, n], link))
}

predict_poisson <- function(n, data, samples, link, ntrys, ...) {
  args <- list(lambda = ilink(samples$eta[, n], link))
  rng_discrete(nrng = nrow(samples$eta), dist = "pois",
               args = args, data = data, ntrys = ntrys)
}

predict_negbinomial <- function(n, data, samples, link, ntrys, ...) {
  args <- list(mu = ilink(samples$eta[, n], link), size = samples$shape)
  rng_discrete(nrng = nrow(samples$eta), dist = "nbinom",
               args = args, data = data, ntrys = ntrys)
}

predict_geometric <- function(n, data, samples, link, ntrys, ...) {
  args <- list(mu = ilink(samples$eta[, n], link), size = 1)
  rng_discrete(nrng = nrow(samples$eta), dist = "nbinom",
               args = args, data = data, ntrys = ntrys)
}

predict_exponential <-  function(n, data, samples, link, ...) {
  args <- list(rate = 1 / ilink(samples$eta[, n], link))
  rng_continuous(nrng = nrow(samples$eta), dist = "exp",
                 args = args, data = data)
}

predict_gamma <- function(n, data, samples, link, ...) {
  args <- list(shape = samples$shape,
               scale = ilink(samples$eta[, n], link) / samples$shape)
  rng_continuous(nrng = nrow(samples$eta), dist = "gamma",
                 args = args, data = data)
}

predict_weibull <- function(n, data, samples, link, ...) {
  args <- list(shape = samples$shape,
               scale = ilink(samples$eta[, n] / samples$shape, link))
  rng_continuous(nrng = nrow(samples$eta), dist = "weibull",
                 args = args, data = data)
}

predict_inverse.gaussian <- function(n, data, samples, link, ...) {
  args <- list(mean = ilink(samples$eta[, n], link), shape = samples$shape)
  rng_continuous(nrng = nrow(samples$eta), dist = "invgauss",
                 args = args, data = data)
}

predict_hurdle_poisson <- function(n, data, samples, link, ...) {
  # theta is the bernoulii hurdle parameter
  theta <- ilink(samples$eta[, n + data$N_trait], "logit")
  lambda <- ilink(samples$eta[, n], link)
  nsamples <- nrow(samples$eta)
  # compare with theta to incorporate the hurdle process
  hu <- runif(nsamples, 0, 1)
  # sample from a truncated poisson distribution
  # by adjusting lambda and adding 1
  t = -log(1 - runif(nsamples) * (1 - exp(-lambda)))
  ifelse(hu < theta, 0, rpois(nsamples, lambda = lambda - t) + 1)
}

predict_hurdle_negbinomial <- function(n, data, samples, link, ...) {
  # theta is the bernoulii hurdle parameter
  theta <- ilink(samples$eta[, n + data$N_trait], "logit")
  mu <- ilink(samples$eta[, n], link)
  nsamples <- nrow(samples$eta)
  # compare with theta to incorporate the hurdle process
  hu <- runif(nsamples, 0, 1)
  # sample from an approximative(!) truncated negbinomial distribution
  # by adjusting mu and adding 1
  t = -log(1 - runif(nsamples) * (1 - exp(-mu)))
  ifelse(hu < theta, 0, rnbinom(nsamples, mu = mu - t, size = samples$shape) + 1)
}

predict_hurdle_gamma <- function(n, data, samples, link, ...) {
  # theta is the bernoulii hurdle parameter
  theta <- ilink(samples$eta[, n + data$N_trait], "logit")
  scale <- ilink(samples$eta[, n], link) / samples$shape
  nsamples <- nrow(samples$eta)
  # compare with theta to incorporate the hurdle process
  hu <- runif(nsamples, 0, 1)
  ifelse(hu < theta, 0, rgamma(nsamples, shape = samples$shape, scale = scale))
}

predict_zero_inflated_poisson <- function(n, data, samples, link, ...) {
  # theta is the bernoulii zero-inflation parameter
  theta <- ilink(samples$eta[, n + data$N_trait], "logit")
  lambda <- ilink(samples$eta[, n], link)
  nsamples <- nrow(samples$eta)
  # compare with theta to incorporate the zero-inflation process
  zi <- runif(nsamples, 0, 1)
  ifelse(zi < theta, 0, rpois(nsamples, lambda = lambda))
}

predict_zero_inflated_negbinomial <- function(n, data, samples, link, ...) {
  # theta is the bernoulii zero-inflation parameter
  theta <- ilink(samples$eta[, n + data$N_trait], "logit")
  mu <- ilink(samples$eta[, n], link)
  nsamples <- nrow(samples$eta)
  # compare with theta to incorporate the zero-inflation process
  zi <- runif(nsamples, 0, 1)
  ifelse(zi < theta, 0, rnbinom(nsamples, mu = mu, size = samples$shape))
}

predict_categorical <- function(n, data, samples, link, ...) {
  ncat <- ifelse(length(data$max_obs) > 1, data$max_obs[n], data$max_obs) 
  p <- pcategorical(1:ncat, eta = samples$eta[, n, ], 
                    ncat = ncat, link = link)
  first_greater(p, target = runif(nrow(samples$eta), min = 0, max = 1))
}

predict_cumulative <- function(n, data, samples, link, ...) {
  predict_ordinal(n = n, data = data, samples = samples, link = link, 
                  family = "cumulative")
}

predict_sratio <- function(n, data, samples, link, ...) {
  predict_ordinal(n = n, data = data, samples = samples, link = link, 
                  family = "sratio")
}

predict_cratio <- function(n, data, samples, link, ...) {
  predict_ordinal(n = n, data = data, samples = samples, link = link, 
                  family = "cratio")
}

predict_acat <- function(n, data, samples, link, ...) {
  predict_ordinal(n = n, data = data, samples = samples, link = link, 
                  family = "acat")
}  

predict_ordinal <- function(n, data, samples, family, link, ...) {
  ncat <- ifelse(length(data$max_obs) > 1, data$max_obs[n], data$max_obs)
  p <- pordinal(1:ncat, eta = samples$eta[, n, ], ncat = ncat, 
                family = family, link = link)
  first_greater(p, target = runif(nrow(samples$eta), min = 0, max = 1))
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
  #   a vector of random values samples from the distribution
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
  #   a vector of random values samples from the distribution
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
  #   rng: samples to be check against truncation boundaries
  #   lb: lower bound
  #   ub: upper bound
  # Returns:
  #   a valid truncated sample or else the closest boundary
  valid_rng <- match(TRUE, rng > lb & rng <= ub)
  if (is.na(valid_rng)) {
    # no valid truncated value found
    # set sample to lb or ub
    # 1e-10 is only to identify the invalid samples later on
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