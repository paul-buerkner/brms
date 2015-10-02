# All functions in this file have the same arguments structure
#
# Args:
#  n: the column of samples to use i.e. the nth obervation in the initial data.frame 
#  data: the data as passed to Stan
#  samples: samples obtained through Stan. Must at least contain variable eta
#  link: the link function
#
# Returns:
#   A vector of length nrow(samples) containing samples from the posterior predictive distribution
predict_gaussian <- function(n, data, samples, link) {
  sigma <- if (!is.null(samples$sigma)) samples$sigma
           else data$sigma
  rnorm(nrow(samples$eta), mean = ilink(samples$eta[, n], link), sd = sigma)
}

predict_student <- function(n, data, samples, link) {
  sigma <- if (!is.null(samples$sigma)) samples$sigma
           else data$sigma
  rstudent(nrow(samples$eta), df = samples$nu, mu = ilink(samples$eta[, n], link), 
           sigma = sigma)
}

predict_cauchy <- function(n, data, samples, link) {
  sigma <- if (!is.null(samples$sigma)) samples$sigma
           else data$sigma
  rstudent(nrow(samples$eta), df = 1, mu = ilink(samples$eta[, n], link), sigma = sigma)
}

predict_lognormal <- function(n, data, samples, link) {
  # link is currently ignored for lognormal models
  # as 'identity' is the only valid link
  sigma <- if (!is.null(samples$sigma)) samples$sigma
           else data$sigma
  rlnorm(nrow(samples$eta), meanlog = samples$eta[, n], sdlog = sigma)
}

predict_multinormal <- function(n, data, samples, link) {
  .fun <- function(i) {
    rmultinormal(1, Sigma = samples$Sigma[i, , ],
                 mu = samples$eta[i, seq(n, data$N, data$N_trait)])
  }
  do.call(rbind, lapply(1:nrow(samples$eta), .fun))
}

predict_binomial <- function(n, data, samples, link) {
  max_obs <- ifelse(length(data$max_obs) > 1, data$max_obs[n], data$max_obs) 
  rbinom(nrow(samples$eta), size = max_obs, 
         prob = ilink(samples$eta[, n], link))
}  

predict_bernoulli <- function(n, data, samples, link) {
  rbinom(nrow(samples$eta), size = 1, 
         prob = ilink(samples$eta[, n], link))
}

predict_poisson <- function(n, data, samples, link) {
  rpois(nrow(samples$eta), lambda = ilink(samples$eta[, n], link))
}

predict_negbinomial <- function(n, data, samples, link) {
  rnbinom(nrow(samples$eta), mu = ilink(samples$eta[, n], link), 
          size = samples$shape)
}

predict_geometric <- function(n, data, samples, link) {
  rnbinom(nrow(samples$eta), mu = ilink(samples$eta[, n], link), size = 1)
}

predict_exponential <-  function(n, data, samples, link) {
  rexp(nrow(samples$eta), rate = ilink(-samples$eta[, n], link))
}

predict_gamma <- function(n, data, samples, link) {
  rgamma(nrow(samples$eta), shape = samples$shape,
         scale = ilink(samples$eta[, n], link) / samples$shape)
}

predict_weibull <- function(n, data, samples, link) {
  rweibull(nrow(samples$eta), shape = samples$shape,
           scale = 1 / (ilink(-samples$eta[, n] / samples$shape, link)))
}

predict_inverse.gaussian <- function(n, data, samples, link) {
  rinv_gaussian(nrow(samples$eta), mu = ilink(samples$eta[, n], link), 
                lambda = samples$shape)
}

predict_hurdle_poisson <- function(n, data, samples, link) {
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

predict_hurdle_negbinomial <- function(n, data, samples, link) {
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

predict_hurdle_gamma <- function(n, data, samples, link) {
  # theta is the bernoulii hurdle parameter
  theta <- ilink(samples$eta[, n + data$N_trait], "logit")
  scale <- ilink(samples$eta[, n], link) / samples$shape
  nsamples <- nrow(samples$eta)
  # compare with theta to incorporate the hurdle process
  hu <- runif(nsamples, 0, 1)
  ifelse(hu < theta, 0, rgamma(nsamples, shape = samples$shape, scale = scale))
}

predict_zero_inflated_poisson <- function(n, data, samples, link) {
  # theta is the bernoulii zero-inflation parameter
  theta <- ilink(samples$eta[, n + data$N_trait], "logit")
  lambda <- ilink(samples$eta[, n], link)
  nsamples <- nrow(samples$eta)
  # compare with theta to incorporate the zero-inflation process
  zi <- runif(nsamples, 0, 1)
  ifelse(zi < theta, 0, rpois(nsamples, lambda = lambda))
}

predict_zero_inflated_negbinomial <- function(n, data, samples, link) {
  # theta is the bernoulii zero-inflation parameter
  theta <- ilink(samples$eta[, n + data$N_trait], "logit")
  mu <- ilink(samples$eta[, n], link)
  nsamples <- nrow(samples$eta)
  # compare with theta to incorporate the zero-inflation process
  zi <- runif(nsamples, 0, 1)
  ifelse(zi < theta, 0, rnbinom(nsamples, mu = mu, size = samples$shape))
}

predict_categorical <- function(n, data, samples, link) {
  ncat <- ifelse(length(data$max_obs) > 1, data$max_obs[n], data$max_obs) 
  p <- pcategorical(1:ncat, eta = samples$eta[, n, ], 
                    ncat = ncat, link = link)
  first_greater(p, target = runif(nrow(samples$eta), min = 0, max = 1))
}

predict_cumulative <- function(n, data, samples, link) {
  predict_ordinal(n = n, data = data, samples = samples, link = link, 
                  family = "cumulative")
}

predict_sratio <- function(n, data, samples, link) {
  predict_ordinal(n = n, data = data, samples = samples, link = link, 
                  family = "sratio")
}

predict_cratio <- function(n, data, samples, link) {
  predict_ordinal(n = n, data = data, samples = samples, link = link, 
                  family = "cratio")
}

predict_acat <- function(n, data, samples, link) {
  predict_ordinal(n = n, data = data, samples = samples, link = link, 
                  family = "acat")
}  

predict_ordinal <- function(n, data, samples, family, link) {
  ncat <- ifelse(length(data$max_obs) > 1, data$max_obs[n], data$max_obs)
  p <- pordinal(1:ncat, eta = samples$eta[, n, ], ncat = ncat, 
                family = family, link = link)
  first_greater(p, target = runif(nrow(samples$eta), min = 0, max = 1))
}