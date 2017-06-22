# All functions in this file have the same arguments structure
# Args:
#  i: the column of draws to use i.e. the ith obervation 
#     in the initial data.frame 
#  draws: A named list returned by extract_draws containing 
#         all required data and samples
#  data: currently ignored; included for compatibility 
#        with loo::loo.function      
# Returns:
#   A vector of length draws$nsamples containing the pointwise 
#   log-likelihood fo the ith observation 
loglik_gaussian <- function(i, draws, data = data.frame()) {
  args <- list(mean = ilink(get_eta(draws$mu, i), draws$f$link), 
               sd = get_sigma(draws$sigma, data = draws$data, i = i))
  # loglik_censor computes the conventional loglik in case of no censoring 
  out <- loglik_censor(dist = "norm", args = args, i = i, data = draws$data)
  out <- loglik_truncate(out, cdf = pnorm, args = args, 
                         i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_student <- function(i, draws, data = data.frame()) {
  args <- list(df = get_auxpar(draws$nu, i = i), 
               mu = ilink(get_eta(draws$mu, i), draws$f$link), 
               sigma = get_sigma(draws$sigma, data = draws$data, i = i))
  out <- loglik_censor(dist = "student_t", args = args, 
                       i = i, data = draws$data)
  out <- loglik_truncate(out, cdf = pstudent_t, args = args, 
                         i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_cauchy <- function(i, draws, data = data.frame()) {
  args <- list(df = 1, mu = ilink(get_eta(draws$mu, i), draws$f$link), 
               sigma = get_sigma(draws$sigma, data = draws$data, i = i))
  out <- loglik_censor(dist = "student_t", args = args, 
                       i = i, data = draws$data)
  out <- loglik_truncate(out, cdf = pstudent_t, args = args, 
                         i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_lognormal <- function(i, draws, data = data.frame()) {
  sigma <- get_sigma(draws$sigma, data = draws$data, i = i)
  args <- list(meanlog = ilink(get_eta(draws$mu, i), draws$f$link), sdlog = sigma)
  out <- loglik_censor(dist = "lnorm", args = args, i = i, data = draws$data)
  out <- loglik_truncate(out, cdf = plnorm, args = args, 
                         i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_skew_normal <- function(i, draws, data = data.frame()) {
  sigma <- get_sigma(draws$sigma, data = draws$data, i = i)
  alpha <- get_auxpar(draws$alpha, i = i)
  mu <- ilink(get_eta(draws$mu, i), draws$f$link)
  args <- nlist(mu, sigma, alpha)
  out <- loglik_censor(
    dist = "skew_normal", args = args, i = i, data = draws$data
  )
  out <- loglik_truncate(
    out, cdf = pskew_normal, args = args, i = i, data = draws$data
  )
  loglik_weight(out, i = i, data = draws$data)
}

loglik_gaussian_mv <- function(i, draws, data = data.frame()) {
  if (!is.null(draws$data$N_trait)) {
    obs <- seq(i, draws$data$N, draws$data$N_trait)
    mu <- ilink(get_eta(draws$mu, obs), draws$f$link)
  } else {
    mu <- ilink(get_eta(draws$mu, i), draws$f$link)[, 1, ]
  }
  out <- sapply(1:draws$nsamples, function(s) 
    dmulti_normal(draws$data$Y[i, ], Sigma = draws$Sigma[s, , ], 
                  mu = mu[s, ], log = TRUE))
  # no truncation allowed
  loglik_weight(out, i = i, data = draws$data)
}

loglik_student_mv <- function(i, draws, data = data.frame()) {
  if (!is.null(draws$data$N_trait)) {
    obs <- seq(i, draws$data$N, draws$data$N_trait)
    mu <- ilink(get_eta(draws$mu, obs), draws$f$link)
  } else {
    mu <- ilink(get_eta(draws$mu, i), draws$f$link)[, 1, ]
  }
  nu <- get_auxpar(draws$nu, obs)
  out <- sapply(1:draws$nsamples, function(s) 
    dmulti_student_t(draws$data$Y[i, ], df = nu[s, ], mu = mu[s, ],
                   Sigma = draws$Sigma[s, , ], log = TRUE))
  # no truncation allowed
  loglik_weight(out, i = i, data = draws$data)
}

loglik_cauchy_mv <- function(i, draws, data = data.frame()) {
  if (!is.null(draws$data$N_trait)) {
    obs <- seq(i, draws$data$N, draws$data$N_trait)
    mu <- ilink(get_eta(draws$mu, obs), draws$f$link)
  } else {
    mu <- ilink(get_eta(draws$mu, i), draws$f$link)[, 1, ]
  }
  out <- sapply(1:draws$nsamples, function(s) 
    dmulti_student_t(draws$data$Y[i, ], df = 1, mu = mu[s, ],
                   Sigma = draws$Sigma[s, , ], log = TRUE))
  # no truncation allowed
  loglik_weight(out, i = i, data = draws$data)
}

loglik_gaussian_cov <- function(i, draws, data = data.frame()) {
  # currently, only ARMA1 processes are implemented
  obs <- with(draws$data, begin_tg[i]:(begin_tg[i] + nobs_tg[i] - 1))
  args <- list(sigma = get_sigma(draws$sigma, data = draws$data, i = obs), 
               se2 = draws$data$se2[obs], nrows = length(obs))
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
  mu <- ilink(get_eta(draws$mu, obs), draws$f$link)
  out <- sapply(1:draws$nsamples, function(s)
    dmulti_normal(
      draws$data$Y[obs], mu = mu[s, ], 
      Sigma = Sigma[s, , ], log = TRUE
    )
  )
  # weights, truncation and censoring not allowed
  out
}

loglik_student_cov <- function(i, draws, data = data.frame()) {
  # currently, only ARMA1 processes are implemented
  obs <- with(draws$data, begin_tg[i]:(begin_tg[i] + nobs_tg[i] - 1))
  args <- list(sigma = get_sigma(draws$sigma, data = draws$data, i = obs), 
               se2 = draws$data$se2[obs], nrows = length(obs))
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
  mu <- ilink(get_eta(draws$mu, obs), draws$f$link)
  nu <- get_auxpar(draws$nu, obs)
  out <- sapply(1:draws$nsamples, function(s)
    dmulti_student_t(
      draws$data$Y[obs], df = nu[s, ], mu = mu[s, ], 
      Sigma = Sigma[s, , ], log = TRUE
    )
  )
  # weights, truncation and censoring not yet allowed
  out
}

loglik_cauchy_cov <- function(i, draws, data = data.frame()) {
  draws$nu <- matrix(rep(1, draws$nsamples))
  loglik_student_cov(i = i, draws = draws, data = data)
}

loglik_gaussian_lagsar <- function(i, draws, data = data.frame()) {
  stopifnot(i == 1)
  .loglik_gaussian_lagsar <- function(s) {
    W_new <- with(draws, diag(data$N) - lagsar[s, ] * data$W)
    mu <- as.numeric(solve(W_new) %*% mu[s, ])
    Sigma <- solve(crossprod(W_new)) * sigma[s]^2
    dmulti_normal(draws$data$Y, mu = mu, Sigma = Sigma, log = TRUE)
  }
  mu <- ilink(get_eta(draws$mu), draws$f$link)
  sigma <- get_sigma(draws$sigma, data = draws$data)
  # weights, truncation and censoring not yet allowed
  sapply(1:draws$nsamples, .loglik_gaussian_lagsar)
}

loglik_student_lagsar <- function(i, draws, data = data.frame()) {
  stopifnot(i == 1)
  .loglik_student_lagsar <- function(s) {
    W_new <- with(draws, diag(data$N) - lagsar[s, ] * data$W)
    mu <- as.numeric(solve(W_new) %*% mu[s, ])
    Sigma <- solve(crossprod(W_new)) * sigma[s]^2
    dmulti_student_t(
      draws$data$Y, df = nu[s], mu = mu, 
      Sigma = Sigma, log = TRUE
    )
  }
  N <- draws$data$N
  mu <- ilink(get_eta(draws$mu), draws$f$link)
  sigma <- get_sigma(draws$sigma, data = draws$data)
  nu <- get_auxpar(draws$nu)
  # weights, truncation and censoring not yet allowed
  sapply(1:draws$nsamples, .loglik_student_lagsar)
}

loglik_gaussian_errorsar <- function(i, draws, data = data.frame()) {
  stopifnot(i == 1)
  .loglik_gaussian_errorsar <- function(s) {
    W_new <- with(draws, diag(data$N) - errorsar[s, ] * data$W)
    Sigma <- solve(crossprod(W_new)) * sigma[s]^2
    dmulti_normal(draws$data$Y, mu = mu[s, ], Sigma = Sigma, log = TRUE)
  }
  mu <- ilink(get_eta(draws$mu), draws$f$link)
  sigma <- get_sigma(draws$sigma, data = draws$data)
  # weights, truncation and censoring not yet allowed
  sapply(1:draws$nsamples, .loglik_gaussian_errorsar)
}

loglik_student_errorsar <- function(i, draws, data = data.frame()) {
  stopifnot(i == 1)
  .loglik_student_errorsar <- function(s) {
    W_new <- with(draws, diag(data$N) - errorsar[s, ] * data$W)
    Sigma <- solve(crossprod(W_new)) * sigma[s]^2
    dmulti_student_t(
      draws$data$Y, df = nu[s], mu = mu[s, ], 
      Sigma = Sigma, log = TRUE
    )
  }
  mu <- ilink(get_eta(draws$mu), draws$f$link)
  sigma <- get_sigma(draws$sigma, data = draws$data)
  nu <- get_auxpar(draws$nu)
  # weights, truncation and censoring not yet allowed
  sapply(1:draws$nsamples, .loglik_student_errorsar)
}

loglik_gaussian_fixed <- function(i, draws, data = data.frame()) {
  stopifnot(i == 1)
  mu <- ilink(get_eta(draws$mu, 1:nrow(draws$data$V)), draws$f$link)
  ulapply(1:draws$nsamples, function(s) 
    dmulti_normal(draws$data$Y, mu = mu[s, ], Sigma = draws$data$V, 
                  log = TRUE))
}

loglik_student_fixed <- function(i, draws, data = data.frame()) {
  stopifnot(i == 1)
  mu <- ilink(get_eta(draws$mu, 1:nrow(draws$data$V)), draws$f$link)
  nu <- get_auxpar(draws$nu, 1:nrow(draws$data$V))
  sapply(1:draws$nsamples, function(s) 
    dmulti_student_t(draws$data$Y, df = nu[s, ], mu = mu[s, ],
                   Sigma = draws$data$V, log = TRUE))
}
  
loglik_cauchy_fixed <- function(i, draws, data = data.frame()) {
  stopifnot(i == 1)
  mu <- ilink(get_eta(draws$mu, 1:nrow(draws$data$V)), draws$f$link)
  sapply(1:draws$nsamples, function(s) 
    dmulti_student_t(draws$data$Y, df = 1, mu = mu[s, ],
                   Sigma = draws$data$V, log = TRUE))
}

loglik_binomial <- function(i, draws, data = data.frame()) {
  trials <- draws$data$trials[i]
  args <- list(size = trials, prob = ilink(get_eta(draws$mu, i), draws$f$link))
  out <- loglik_censor(dist = "binom", args = args, i = i, data = draws$data)
  out <- loglik_truncate(out, cdf = pbinom, args = args, 
                         i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}  

loglik_bernoulli <- function(i, draws, data = data.frame()) {
  eta <-  get_eta(draws$mu, i)
  args <- list(size = 1, prob = ilink(eta, draws$f$link))
  out <- loglik_censor(dist = "binom", args = args, i = i, 
                       data = draws$data)
  # no truncation allowed
  loglik_weight(out, i = i, data = draws$data)
}

loglik_poisson <- function(i, draws, data = data.frame()) {
  args <- list(lambda = ilink(get_eta(draws$mu, i), draws$f$link))
  out <- loglik_censor(dist = "pois", args = args, i = i, 
                       data = draws$data)
  out <- loglik_truncate(out, cdf = ppois, args = args, 
                         i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_negbinomial <- function(i, draws, data = data.frame()) {
  shape <- get_shape(draws$shape, data = draws$data, i = i)
  args <- list(mu = ilink(get_eta(draws$mu, i), draws$f$link), size = shape)
  out <- loglik_censor(dist = "nbinom", args = args, i = i, data = draws$data)
  out <- loglik_truncate(out, cdf = pnbinom, args = args, 
                         i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_geometric <- function(i, draws, data = data.frame()) {
  args <- list(mu = ilink(get_eta(draws$mu, i), draws$f$link), size = 1)
  out <- loglik_censor(dist = "nbinom", args = args, i = i, data = draws$data)
  out <- loglik_truncate(out, cdf = pnbinom, args = args,
                         i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_exponential <- function(i, draws, data = data.frame()) {
  args <- list(rate = 1 / ilink(get_eta(draws$mu, i), draws$f$link))
  out <- loglik_censor(dist = "exp", args = args, i = i, data = draws$data)
  out <- loglik_truncate(out, cdf = pexp, args = args,
                         i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_gamma <- function(i, draws, data = data.frame()) {
  shape <- get_shape(draws$shape, data = draws$data, i = i)
  args <- list(shape = shape, 
               scale = ilink(get_eta(draws$mu, i), draws$f$link) / shape)
  out <- loglik_censor(dist = "gamma", args = args, i = i, data = draws$data)
  out <- loglik_truncate(out, cdf = pgamma, args = args,
                         i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_weibull <- function(i, draws, data = data.frame()) {
  shape <- get_shape(draws$shape, data = draws$data, i = i)
  scale <- ilink(get_eta(draws$mu, i) / shape, draws$f$link)
  args <- list(shape = shape, scale = scale)
  out <- loglik_censor(dist = "weibull", args = args, 
                       i = i, data = draws$data)
  out <- loglik_truncate(out, cdf = pweibull, args = args,
                         i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_frechet <- function(i, draws, data = data.frame()) {
  nu <- get_auxpar(draws$nu, i = i)
  scale <- ilink(get_eta(draws$mu, i), draws$f$link) / gamma(1 - 1 / nu)
  args <- list(scale = scale, shape = nu)
  out <- loglik_censor(dist = "frechet", args = args, 
                       i = i, data = draws$data)
  out <- loglik_truncate(out, cdf = pfrechet, args = args,
                         i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_gen_extreme_value <- function(i, draws, data = data.frame()) {
  sigma <- get_sigma(draws$sigma, data = draws$data, i = i)
  xi <- get_auxpar(draws$xi, i = i)
  mu <- ilink(get_eta(draws$mu, i), draws$f$link)
  args <- nlist(mu, sigma, xi)
  out <- loglik_censor(dist = "gen_extreme_value", args = args, 
                       i = i, data = draws$data)
  out <- loglik_truncate(out, cdf = pgen_extreme_value, 
                         args = args, i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_inverse.gaussian <- function(i, draws, data = data.frame()) {
  args <- list(mu = ilink(get_eta(draws$mu, i), draws$f$link), 
               shape = get_shape(draws$shape, data = draws$data, i = i))
  out <- loglik_censor(dist = "inv_gaussian", args = args, 
                       i = i, data = draws$data)
  out <- loglik_truncate(out, cdf = pinv_gaussian, args = args,
                         i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_exgaussian <- function(i, draws, data = data.frame()) {
  args <- list(mu = ilink(get_eta(draws$mu, i), draws$f$link), 
               sigma = get_sigma(draws$sigma, data = draws$data, i = i),
               beta = get_auxpar(draws$beta, i = i))
  out <- loglik_censor(dist = "exgaussian", args = args, 
                       i = i, data = draws$data)
  out <- loglik_truncate(out, cdf = pexgaussian, args = args,
                         i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_wiener <- function(i, draws, data = data.frame()) {
  args <- list(delta = ilink(get_eta(draws$mu, i), draws$f$link), 
               alpha = get_auxpar(draws$bs, i = i),
               tau = get_auxpar(draws$ndt, i = i),
               beta = get_auxpar(draws$bias, i = i),
               resp = draws$data[["dec"]][i])
  out <- do.call("dwiener", c(draws$data$Y[i], args, log = TRUE))
  loglik_weight(out, i = i, data = draws$data)
}

loglik_beta <- function(i, draws, data = data.frame()) {
  mu <- ilink(get_eta(draws$mu, i), draws$f$link)
  phi <- get_auxpar(draws$phi, i)
  args <- list(shape1 = mu * phi, shape2 = (1 - mu) * phi)
  out <- loglik_censor(dist = "beta", args = args, i = i, data = draws$data)
  out <- loglik_truncate(out, cdf = pbeta, args = args,
                         i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_von_mises <- function(i, draws, data = data.frame()) {
  args <- list(mu = ilink(get_eta(draws$mu, i), draws$f$link), 
               kappa = get_auxpar(draws$kappa, i = i))
  out <- loglik_censor(dist = "von_mises", args = args, 
                       i = i, data = draws$data)
  out <- loglik_truncate(out, cdf = pvon_mises, args = args,
                         i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_asym_laplace <- function(i, draws, ...) {
  args <- list(mu = ilink(get_eta(draws$mu, i), draws$f$link), 
               sigma = get_sigma(draws$sigma, data = draws$data, i = i),
               quantile = get_auxpar(draws$quantile, i = i))
  out <- loglik_censor(dist = "asym_laplace", args = args, 
                       i = i, data = draws$data)
  out <- loglik_truncate(out, cdf = pvon_mises, args = args,
                         i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_hurdle_poisson <- function(i, draws, data = data.frame()) {
  theta <- get_zi_hu(draws, i, par = "hu")
  args <- list(lambda = ilink(get_eta(draws$mu, i), draws$f$link))
  out <- loglik_hurdle_discrete(pdf = dpois, theta = theta, 
                                args = args, i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_hurdle_negbinomial <- function(i, draws, data = data.frame()) {
  theta <- get_zi_hu(draws, i, par = "hu")
  args <- list(mu = ilink(get_eta(draws$mu, i), draws$f$link), 
               size = get_shape(draws$shape, data = draws$data, i = i))
  out <- loglik_hurdle_discrete(pdf = dnbinom, theta = theta, 
                                args = args, i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_hurdle_gamma <- function(i, draws, data = data.frame()) {
  theta <- get_zi_hu(draws, i, par = "hu")
  shape <- get_shape(draws$shape, data = draws$data, i = i)
  args <- list(shape = shape, 
               scale = ilink(get_eta(draws$mu, i), draws$f$link) / shape)
  out <- loglik_hurdle_continuous(pdf = dgamma, theta = theta, 
                                  args = args, i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_hurdle_lognormal <- function(i, draws, data = data.frame()) {
  theta <- get_zi_hu(draws, i, par = "hu")
  sigma <- get_sigma(draws$sigma, data = draws$data, i = i)
  args <- list(meanlog = ilink(get_eta(draws$mu, i), draws$f$link), sdlog = sigma)
  out <- loglik_hurdle_continuous(pdf = dlnorm, theta = theta, 
                                  args = args, i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_zero_inflated_poisson <- function(i, draws, data = data.frame()) {
  theta <- get_zi_hu(draws, i, par = "zi")
  args <- list(lambda = ilink(get_eta(draws$mu, i), draws$f$link))
  out <- loglik_zero_inflated(pdf = dpois, theta = theta, 
                              args = args, i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_zero_inflated_negbinomial <- function(i, draws, data = data.frame()) {
  theta <- get_zi_hu(draws, i, par = "zi")
  args <- list(mu = ilink(get_eta(draws$mu, i), draws$f$link), 
               size = get_shape(draws$shape, data = draws$data, i = i))
  out <- loglik_zero_inflated(pdf = dnbinom, theta = theta, 
                              args = args, i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_zero_inflated_binomial <- function(i, draws, data = data.frame()) {
  trials <- draws$data$trials[i] 
  theta <- get_zi_hu(draws, i, par = "zi")
  args <- list(size = trials, prob = ilink(get_eta(draws$mu, i), draws$f$link))
  out <- loglik_zero_inflated(pdf = dbinom, theta = theta, 
                              args = args, i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_zero_inflated_beta <- function(i, draws, data = data.frame()) {
  theta <- get_zi_hu(draws, i, par = "zi")
  mu <- ilink(get_eta(draws$mu, i), draws$f$link)
  phi <- get_auxpar(draws$phi, i)
  args <- list(shape1 = mu * phi, shape2 = (1 - mu) * phi)
  # zi_beta is technically a hurdle model
  out <- loglik_hurdle_continuous(pdf = dbeta, theta = theta, 
                                  args = args, i = i, data = draws$data)
  loglik_weight(out, i = i, data = draws$data)
}

loglik_zero_one_inflated_beta <- function(i, draws, data = data.frame()) {
  zoi <- get_auxpar(draws$zoi, i)
  coi <- get_auxpar(draws$coi, i)
  if (draws$data$Y[i] %in% c(0, 1)) {
    out <- dbinom(1, size = 1, prob = zoi, log = TRUE) + 
      dbinom(draws$data$Y[i], size = 1, prob = coi, log = TRUE)
  } else {
    phi <- get_auxpar(draws$phi, i)
    mu <- ilink(get_eta(draws$mu, i), draws$f$link)
    args <- list(shape1 = mu * phi, shape2 = (1 - mu) * phi)
    out <- dbinom(0, size = 1, prob = zoi, log = TRUE) + 
      do.call(dbeta, c(draws$data$Y[i], args, log = TRUE))
  }
  loglik_weight(out, i = i, data = draws$data)
}

loglik_categorical <- function(i, draws, data = data.frame()) {
  if (draws$f$link == "logit") {
    p <- cbind(rep(0, draws$nsamples), get_eta(draws$mu, i)[, 1, ])
    out <- p[, draws$data$Y[i]] - log(rowSums(exp(p)))
  } else {
    stop(paste("Link", draws$f$link, "not supported"))
  }
  loglik_weight(out, i = i, data = draws$data)
}

loglik_cumulative <- function(i, draws, data = data.frame()) {
  ncat <- draws$data$ncat
  disc <- get_disc(draws, i, ncat)
  eta <- disc * get_eta(draws$mu, i)
  y <- draws$data$Y[i]
  if (y == 1) { 
    out <- log(ilink(eta[, 1, 1], draws$f$link))
  } else if (y == ncat) {
    out <- log(1 - ilink(eta[, 1, y - 1], draws$f$link)) 
  } else {
    out <- log(ilink(eta[, 1, y], draws$f$link) - 
                  ilink(eta[, 1, y - 1], draws$f$link))
  }
  loglik_weight(out, i = i, data = draws$data)
}

loglik_sratio <- function(i, draws, data = data.frame()) {
  ncat <- draws$data$ncat
  disc <- get_disc(draws, i, ncat)
  eta <- disc * get_eta(draws$mu, i)
  y <- draws$data$Y[i]
  q <- sapply(1:min(y, ncat - 1), function(k) 
    1 - ilink(eta[, 1, k], draws$f$link))
  if (y == 1) {
    out <- log(1 - q[, 1]) 
  } else if (y == 2) {
    out <- log(1 - q[, 2]) + log(q[, 1])
  } else if (y == ncat) {
    out <- rowSums(log(q))
  } else {
    out <- log(1 - q[, y]) + rowSums(log(q[, 1:(y - 1)]))
  }
  loglik_weight(out, i = i, data = draws$data)
}

loglik_cratio <- function(i, draws, data = data.frame()) {
  ncat <- draws$data$ncat
  disc <- get_disc(draws, i, ncat)
  eta <- disc * get_eta(draws$mu, i)
  y <- draws$data$Y[i]
  q <- sapply(1:min(y, ncat-1), function(k) 
    ilink(eta[, 1, k], draws$f$link))
  if (y == 1) {
    out <- log(1 - q[, 1])
  }  else if (y == 2) {
    out <- log(1 - q[, 2]) + log(q[, 1])
  } else if (y == ncat) {
    out <- rowSums(log(q))
  } else {
    out <- log(1 - q[, y]) + rowSums(log(q[, 1:(y - 1)]))
  }
  loglik_weight(out, i = i, data = draws$data)
}

loglik_acat <- function(i, draws, data = data.frame()) {
  ncat <- draws$data$ncat
  disc <- get_disc(draws, i, ncat)
  eta <- disc * get_eta(draws$mu, i)
  y <- draws$data$Y[i]
  if (draws$f$link == "logit") { # more efficient calculation 
    q <- sapply(1:(ncat - 1), function(k) eta[, 1, k])
    p <- cbind(rep(0, nrow(eta)), q[, 1], 
               matrix(0, nrow = nrow(eta), ncol = ncat - 2))
    if (ncat > 2) {
      p[, 3:ncat] <- sapply(3:ncat, function(k) rowSums(q[, 1:(k - 1)]))
    }
    out <- p[, y] - log(rowSums(exp(p)))
  } else {
    q <- sapply(1:(ncat - 1), function(k) 
      ilink(eta[, 1, k], draws$f$link))
    p <- cbind(apply(1 - q[, 1:(ncat - 1)], 1, prod), 
               matrix(0, nrow = nrow(eta), ncol = ncat - 1))
    if (ncat > 2) {
      p[, 2:(ncat - 1)] <- sapply(2:(ncat - 1), function(k) 
        apply(as.matrix(q[, 1:(k - 1)]), 1, prod) * 
          apply(as.matrix(1 - q[, k:(ncat - 1)]), 1, prod))
    }
    p[, ncat] <- apply(q[, 1:(ncat - 1)], 1, prod)
    out <- log(p[, y]) - log(apply(p, 1, sum))
  }
  loglik_weight(out, i = i, data = draws$data)
}

loglik_mixture <- function(i, draws, data = data.frame()) {
  families <- family_names(draws$f)
  theta <- get_theta(draws, i = i)
  out <- array(NA, dim = dim(theta))
  for (j in seq_along(families)) {
    loglik_fun <- paste0("loglik_", families[j])
    loglik_fun <- get(loglik_fun, asNamespace("brms"))
    auxpars <- valid_auxpars(families[j])
    tmp_draws <- list(
      f = draws$f$mix[[j]],
      nsamples = draws[["nsamples"]],
      data = draws[["data"]]
    )
    for (ap in auxpars) {
      tmp_draws[[ap]] <- draws[[paste0(ap, j)]]
    }
    out[, j] <- exp(log(theta[, j]) + loglik_fun(i, tmp_draws))
  }
  if (isTRUE(draws[["pp_mixture"]])) {
    out <- log(out) - log(rowSums(out))
  } else {
    out <- log(rowSums(out))
  }
  out
}

# ----------- loglik helper-functions -----------

loglik_censor <- function(dist, args, i, data) {
  # compute (possibly censored) loglik values
  # Args:
  #   dist: name of a distribution for which the functions
  #         d<dist> (pdf) and p<dist> (cdf) are available
  #   args: additional arguments passed to pdf and cdf
  #   data: data initially passed to Stan
  # Returns:
  #   vector of loglik values
  pdf <- get(paste0("d", dist), mode = "function")
  cdf <- get(paste0("p", dist), mode = "function")
  if (is.null(data$cens) || data$cens[i] == 0) {
    do.call(pdf, c(data$Y[i], args, log = TRUE))
  } else if (data$cens[i] == 1) {
    do.call(cdf, c(data$Y[i], args, lower.tail = FALSE, log.p = TRUE))
  } else if (data$cens[i] == -1) {
    do.call(cdf, c(data$Y[i], args, log.p = TRUE))
  } else if (data$cens[i] == 2) {
    log(do.call(cdf, c(data$rcens[i], args)) - 
          do.call(cdf, c(data$Y[i], args)))
  }
}

loglik_truncate <- function(x, cdf, args, i, data) {
  # adjust logLik in truncated models
  # Args:
  #   x: vector of loglik values
  #   cdf: a cumulative distribution function 
  #   args: arguments passed to cdf
  #   i: observation number
  #   data: data initally passed to Stan
  # Returns:
  #   vector of loglik values
  if (!(is.null(data$lb) && is.null(data$ub))) {
    lb <- if (is.null(data$lb)) -Inf else data$lb[i]
    ub <- if (is.null(data$ub)) -Inf else data$ub[i]
    x - log(do.call(cdf, c(ub, args)) - do.call(cdf, c(lb, args)))
  } else {
    x
  }
}

loglik_weight <- function(x, i, data) {
  # weight loglik values according to defined weights
  # Args:
  #   x: vector of loglik values
  #   i: observation number
  #   data: data initially passed to Stan
  # Returns:
  #   vector of loglik values
  if ("weights" %in% names(data)) {
    x * data$weights[i]
  } else {
    x
  }
}

loglik_hurdle_discrete <- function(pdf, theta, args, i, data) {
  # loglik values for discrete hurdle models
  # Args:
  #  pdf: a probability density function 
  #  theta: bernoulli hurdle parameter
  #  args: arguments passed to pdf
  #  data: data initially passed to Stan
  # Returns:
  #   vector of loglik values
  if (data$Y[i] == 0) {
    dbinom(1, size = 1, prob = theta, log = TRUE)
  } else {
    dbinom(0, size = 1, prob = theta, log = TRUE) + 
      do.call(pdf, c(data$Y[i], args, log = TRUE)) -
      log(1 - do.call(pdf, c(0, args)))
  }
}

loglik_hurdle_continuous <- function(pdf, theta, args, i, data) {
  # loglik values for continuous hurdle models
  # does not call log(1 - do.call(pdf, c(0, args)))
  # Args:
  #   same as loglik_hurdle_discrete
  if (data$Y[i] == 0) {
    dbinom(1, size = 1, prob = theta, log = TRUE)
  } else {
    dbinom(0, size = 1, prob = theta, log = TRUE) + 
      do.call(pdf, c(data$Y[i], args, log = TRUE))
  }
}

loglik_zero_inflated <- function(pdf, theta, args, i, data) {
  # loglik values for zero-inflated models
  # Args:
  #  pdf: a probability density function 
  #  theta: bernoulli zero-inflation parameter
  #  args: arguments passed to pdf
  #  data: data initially passed to Stan
  # Returns:
  #   vector of loglik values
  if (data$Y[i] == 0) {
    log(dbinom(1, size = 1, prob = theta) + 
        dbinom(0, size = 1, prob = theta) *
          do.call(pdf, c(0, args)))
  } else {
    dbinom(0, size = 1, prob = theta, log = TRUE) +
      do.call(pdf, c(data$Y[i], args, log = TRUE))
  }
}
