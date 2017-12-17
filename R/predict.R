predict_internal <- function(draws, ...) {
  UseMethod("predict_internal")
}

#' @export
predict_internal.mvbrmsdraws <- function(draws, ...) {
  if (length(draws$mvpars$rescor)) {
    draws$mvpars$Mu <- get_Mu(draws)
    draws$mvpars$Sigma <- get_Sigma(draws)
    out <- predict_internal.brmsdraws(draws, ...)
  } else {
    out <- lapply(draws$resps, predict_internal, ...)
    along <- ifelse(length(out) > 1L, 3, 2)
    out <- do.call(abind, c(out, along = along))
  }
  out
}

#' @export
predict_internal.brmsdraws <- function(draws, summary = TRUE, transform = NULL,
                                       sort = FALSE, robust = FALSE, 
                                       probs = c(0.025, 0.975), ...) {
  for (dp in names(draws$dpars)) {
    draws$dpars[[dp]] <- get_dpar(draws, dpar = dp)
  }
  predict_fun <- paste0("predict_", draws$f$fun)
  predict_fun <- get(predict_fun, asNamespace("brms"))
  N <- choose_N(draws)
  out <- lapply(seq_len(N), predict_fun, draws = draws, ...)
  if (grepl("_mv$", draws$f$fun)) {
    out <- do.call(abind, c(out, along = 3))
    out <- aperm(out, perm = c(1, 3, 2))
    dimnames(out)[[3]] <- names(draws$resps)
  } else {
    out <- do.call(cbind, out) 
  }
  # percentage of invalid samples for truncated discrete models
  # should always be zero for all other models
  pct_invalid <- get_pct_invalid(out, lb = draws$data$lb, ub = draws$data$ub) 
  if (pct_invalid >= 0.01) {
    warning2(
      round(pct_invalid * 100), "% of all predicted values ", 
      "were invalid. Increasing argument 'ntrys' may help."
    )
  }
  out <- reorder_obs(out, draws$old_order, sort = sort)
  # transform predicted response samples before summarizing them 
  if (!is.null(transform)) {
    out <- do.call(transform, list(out))
  }
  if (summary) {
    if (is_ordinal(draws$f) || is_categorical(draws$f)) {
      # compute frequencies of categories 
      out <- posterior_table(out, levels = seq_len(max(draws$data$ncat)))
    } else {
      out <- posterior_summary(out, probs = probs, robust = robust)
    }
  }
  out
}

# All predict_<family> functions have the same arguments structure
# Args:
#  i: the column of draws to use i.e. the ith obervation 
#     in the initial data.frame 
#  draws: A named list returned by extract_draws containing 
#         all required data and samples
#  ...: ignored arguments
# Returns:
#   A vector of length draws$nsamples containing samples
#   from the posterior predictive distribution
predict_gaussian <- function(i, draws, ...) {
  args <- list(
    mean = get_dpar(draws, "mu", i = i), 
    sd = get_dpar(draws, "sigma", i = i)
  )
  rng_continuous(
    nrng = draws$nsamples, dist = "norm", args = args, 
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_student <- function(i, draws, ...) {
  args <- list(
    df = get_dpar(draws, "nu", i = i), 
    mu = get_dpar(draws, "mu", i = i), 
    sigma = get_dpar(draws, "sigma", i = i)
  )
  rng_continuous(
    nrng = draws$nsamples, dist = "student_t", args = args, 
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_lognormal <- function(i, draws, ...) {
  args <- list(
    meanlog = get_dpar(draws, "mu", i = i), 
    sdlog = get_dpar(draws, "sigma", i = i)
  )
  rng_continuous(
    nrng = draws$nsamples, dist = "lnorm", args = args, 
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_skew_normal <- function(i, draws, ...) {
  sigma <- get_dpar(draws, "sigma", i = i)
  alpha <- get_dpar(draws, "alpha", i = i)
  mu <- get_dpar(draws, "mu", i = i)
  args <- nlist(mu, sigma, alpha)
  rng_continuous(
    nrng = draws$nsamples, dist = "skew_normal", args = args, 
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_gaussian_mv <- function(i, draws, ...) {
  Mu <- get_Mu(draws, i = i)
  Sigma <- get_Sigma(draws, i = i)
  rmn <- function(s) {
    rmulti_normal(1, mu = Mu[s, ], Sigma = Sigma[s, , ])
  }
  do.call(rbind, lapply(1:draws$nsamples, rmn))
}

predict_student_mv <- function(i, draws, ...) {
  nu <- get_dpar(draws, "nu", i = i)
  Mu <- get_Mu(draws, i = i)
  Sigma <- get_Sigma(draws, i = i)
  rmst <- function(s) {
    rmulti_student_t(1, df = nu[s], mu = Mu[s, ], Sigma = Sigma[s, , ])
  }
  do.call(rbind, lapply(1:draws$nsamples, rmst))
}

predict_gaussian_cov <- function(i, draws, ...) {
  # currently, only ARMA1 processes are implemented
  obs <- with(draws$ac, begin_tg[i]:(begin_tg[i] + nobs_tg[i] - 1))
  mu <- as.matrix(get_dpar(draws, "mu", i = obs))
  args <- list(
    sigma = get_dpar(draws, "sigma", i = obs),
    se = draws$data$se[obs], nrows = length(obs)
  )
  ar <- draws$ac$ar
  ma <- draws$ac$ma
  if (!is.null(ar) && is.null(ma)) {
    args$ar <- ar
    Sigma <- do.call(get_cov_matrix_ar1, args)
  } else if (is.null(ar) && !is.null(ma)) {
    args$ma <- ma
    Sigma <- do.call(get_cov_matrix_ma1, args)
  } else if (!is.null(ar) && !is.null(ma)) {
    args[c("ar", "ma")] <- draws$ac[c("ar", "ma")]
    Sigma <- do.call(get_cov_matrix_arma1, args)
  } else {
    Sigma <- do.call(get_cov_matrix_ident, args)
  }
  .fun <- function(s) {
    rmulti_normal(1, mu = mu[s, ], Sigma = Sigma[s, , ])
  }
  do.call(rbind, lapply(1:draws$nsamples, .fun))
}

predict_student_cov <- function(i, draws, ...) {
  # currently, only ARMA1 processes are implemented
  obs <- with(draws$ac, begin_tg[i]:(begin_tg[i] + nobs_tg[i] - 1))
  args <- list(
    sigma = get_dpar(draws, "sigma", i = obs), 
    se = draws$data$se[obs], nrows = length(obs)
  )
  ar <- draws$ac$ar
  ma <- draws$ac$ma
  if (!is.null(ar) && is.null(ma)) {
    args$ar <- ar
    Sigma <- do.call(get_cov_matrix_ar1, args)
  } else if (is.null(ar) && !is.null(ma)) {
    args$ma <- ma
    Sigma <- do.call(get_cov_matrix_ma1, args)
  } else if (!is.null(ar) && !is.null(ma)) {
    args[c("ar", "ma")] <- draws$ac[c("ar", "ma")]
    Sigma <- do.call(get_cov_matrix_arma1, args)
  } else {
    Sigma <- do.call(get_cov_matrix_ident, args)
  }
  mu <- as.matrix(get_dpar(draws, "mu", i = obs))
  nu <- as.matrix(get_dpar(draws, "nu", i = obs))
  .fun <- function(s) {
    rmulti_student_t(1, df = nu[s, ], mu = mu[s, ], Sigma = Sigma[s, , ])
  }
  do.call(rbind, lapply(1:draws$nsamples, .fun))
}

predict_gaussian_lagsar <- function(i, draws, ...) {
  stopifnot(i == 1)
  .predict_gaussian_lagsar <- function(s) {
    W_new <- with(draws, diag(nobs) - ac$lagsar[s] * ac$W)
    mu <- as.numeric(solve(W_new) %*% mu[s, ])
    Sigma <- solve(crossprod(W_new)) * sigma[s]^2
    rmulti_normal(1, mu = mu, Sigma = Sigma)
  }
  mu <- get_dpar(draws, "mu")
  sigma <- get_dpar(draws, "sigma")
  do.call(rbind, lapply(1:draws$nsamples, .predict_gaussian_lagsar))
}

predict_student_lagsar <- function(i, draws, ...) {
  stopifnot(i == 1)
  .predict_student_lagsar <- function(s) {
    W_new <- with(draws, diag(nobs) - ac$lagsar[s] * ac$W)
    mu <- as.numeric(solve(W_new) %*% mu[s, ])
    Sigma <- solve(crossprod(W_new)) * sigma[s]^2
    rmulti_student_t(1, df = nu[s], mu = mu, Sigma = Sigma)
  }
  mu <- get_dpar(draws, "mu")
  sigma <- get_dpar(draws, "sigma")
  nu <- get_dpar(draws, "nu")
  do.call(rbind, lapply(1:draws$nsamples, .predict_student_lagsar))
}

predict_gaussian_errorsar <- function(i, draws, ...) {
  stopifnot(i == 1)
  .predict_gaussian_errorsar <- function(s) {
    W_new <- with(draws, diag(nobs) - ac$errorsar[s] * ac$W)
    Sigma <- solve(crossprod(W_new)) * sigma[s]^2
    rmulti_normal(1, mu = mu[s, ], Sigma = Sigma)
  }
  mu <- get_dpar(draws, "mu")
  sigma <- get_dpar(draws, "sigma")
  do.call(rbind, lapply(1:draws$nsamples, .predict_gaussian_errorsar))
}

predict_student_errorsar <- function(i, draws, ...) {
  stopifnot(i == 1)
  .predict_student_errorsar <- function(s) {
    W_new <- with(draws, diag(nobs) - ac$errorsar[s] * ac$W)
    Sigma <- solve(crossprod(W_new)) * sigma[s]^2
    rmulti_student_t(1, df = nu[s], mu = mu[s, ], Sigma = Sigma)
  }
  mu <- get_dpar(draws, "mu")
  sigma <- get_dpar(draws, "sigma")
  nu <- get_dpar(draws, "nu")
  do.call(rbind, lapply(1:draws$nsamples, .predict_student_errorsar))
}

predict_gaussian_fixed <- function(i, draws, ...) {
  stopifnot(i == 1)
  mu <- as.matrix(get_dpar(draws, "mu"))
  .fun <- function(s) {
    rmulti_normal(1, mu = mu[s, ], Sigma = draws$ac$V)
  }
  do.call(rbind, lapply(1:draws$nsamples, .fun))
}

predict_student_fixed <- function(i, draws, ...) {
  stopifnot(i == 1)
  mu <- as.matrix(get_dpar(draws, "mu"))
  nu <- as.matrix(get_dpar(draws, "nu"))
  .fun <- function(s) {
    rmulti_student_t(1, df = nu[s, ], mu = mu[s, ], Sigma = draws$ac$V)
  }
  do.call(rbind, lapply(1:draws$nsamples, .fun))
}

predict_binomial <- function(i, draws, ntrys = 5, ...) {
  args <- list(
    size = draws$data$trials[i], 
    prob = get_dpar(draws, "mu", i = i)
  )
  rng_discrete(
    nrng = draws$nsamples, dist = "binom", args = args, 
    lb = draws$data$lb[i], ub = draws$data$ub[i], 
    ntrys = ntrys
  )
}

predict_bernoulli <- function(i, draws, ...) {
  # truncation not useful
  mu <- get_dpar(draws, "mu", i = i)
  rbinom(length(mu), size = 1, prob = mu)
}

predict_poisson <- function(i, draws, ntrys = 5, ...) {
  args <- list(lambda = get_dpar(draws, "mu", i = i))
  rng_discrete(
    nrng = draws$nsamples, dist = "pois", args = args, 
    lb = draws$data$lb[i], ub = draws$data$ub[i],
    ntrys = ntrys
  )
}

predict_negbinomial <- function(i, draws, ntrys = 5, ...) {
  args <- list(
    mu = get_dpar(draws, "mu", i = i), 
    size = get_dpar(draws, "shape", i = i)
  )
  rng_discrete(
    nrng = draws$nsamples, dist = "nbinom", args = args, 
    lb = draws$data$lb[i], ub = draws$data$ub[i],
    ntrys = ntrys
  )
}

predict_geometric <- function(i, draws, ntrys = 5, ...) {
  args <- list(
    mu = get_dpar(draws, "mu", i = i), 
    size = 1
  )
  rng_discrete(
    nrng = draws$nsamples, dist = "nbinom", args = args, 
    lb = draws$data$lb[i], ub = draws$data$ub[i], 
    ntrys = ntrys
  )
}

predict_exponential <- function(i, draws, ...) {
  args <- list(rate = 1 / get_dpar(draws, "mu", i = i))
  rng_continuous(
    nrng = draws$nsamples, dist = "exp", args = args, 
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_gamma <- function(i, draws, ...) {
  shape <- get_dpar(draws, "shape", i = i)
  args <- nlist(shape, scale = get_dpar(draws, "mu", i = i) / shape)
  rng_continuous(
    nrng = draws$nsamples, dist = "gamma", args = args, 
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_weibull <- function(i, draws, ...) {
  shape <- get_dpar(draws, "shape", i = i)
  mu <- get_dpar(draws, "mu", i = i) 
  scale <- ilink(mu / shape, draws$f$link)
  args <- list(shape = shape, scale = scale)
  rng_continuous(
    nrng = draws$nsamples, dist = "weibull", args = args, 
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_frechet <- function(i, draws, ...) {
  nu <- get_dpar(draws, "nu", i = i)
  scale <- get_dpar(draws, "mu", i = i) / gamma(1 - 1 / nu)
  args <- list(scale = scale, shape = nu)
  rng_continuous(
    nrng = draws$nsamples, dist = "frechet", args = args, 
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_gen_extreme_value <- function(i, draws, ...) {
  sigma <- get_dpar(draws, "sigma", i = i)
  xi <- get_dpar(draws, "xi", i = i)
  mu <- get_dpar(draws, "mu", i = i)
  args <- nlist(mu, sigma, xi)
  rng_continuous(
    nrng = draws$nsamples, dist = "gen_extreme_value", 
    args = args, lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_inverse.gaussian <- function(i, draws, ...) {
  args <- list(
    mu = get_dpar(draws, "mu", i = i), 
    shape = get_dpar(draws, "shape", i = i)
  )
  rng_continuous(
    nrng = draws$nsamples, dist = "inv_gaussian", args = args, 
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_exgaussian <- function(i, draws, ...) {
  args <- list(
    mu = get_dpar(draws, "mu", i = i), 
    sigma = get_dpar(draws, "sigma", i = i),
    beta = get_dpar(draws, "beta", i = i)
  )
  rng_continuous(
    nrng = draws$nsamples, dist = "exgaussian", args = args, 
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_wiener <- function(i, draws, negative_rt = FALSE, ...) {
  args <- list(
    delta = get_dpar(draws, "mu", i = i), 
    alpha = get_dpar(draws, "bs", i = i),
    tau = get_dpar(draws, "ndt", i = i),
    beta = get_dpar(draws, "bias", i = i),
    types = if (negative_rt) c("q", "resp") else "q"
  )
  out <- rng_continuous(
    nrng = 1, dist = "wiener", args = args, 
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
  if (negative_rt) {
    # code lower bound responses as negative RTs
    out <- out[["q"]] * ifelse(out[["resp"]], 1, -1)
  }
  out
}

predict_beta <- function(i, draws, ...) {
  mu <- get_dpar(draws, "mu", i = i)
  phi <- get_dpar(draws, "phi", i = i)
  args <- list(shape1 = mu * phi, shape2 = (1 - mu) * phi)
  rng_continuous(
    nrng = draws$nsamples, dist = "beta", args = args, 
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_von_mises <- function(i, draws, ...) {
  args <- list(
    mu = get_dpar(draws, "mu", i = i), 
    kappa = get_dpar(draws, "kappa", i = i)
  )
  rng_continuous(
    nrng = draws$nsamples, dist = "von_mises", args = args,
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_asym_laplace <- function(i, draws, ...) {
  args <- list(
    mu = get_dpar(draws, "mu", i = i), 
    sigma = get_dpar(draws, "sigma", i = i),
    quantile = get_dpar(draws, "quantile", i = i)
  )
  rng_continuous(
    nrng = draws$nsamples, dist = "asym_laplace", args = args, 
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_hurdle_poisson <- function(i, draws, ...) {
  # theta is the bernoulli hurdle parameter
  theta <- get_dpar(draws, "hu", i = i) 
  lambda <- get_dpar(draws, "mu", i = i)
  ndraws <- draws$nsamples
  # compare with theta to incorporate the hurdle process
  hu <- runif(ndraws, 0, 1)
  # sample from a truncated poisson distribution
  # by adjusting lambda and adding 1
  t = -log(1 - runif(ndraws) * (1 - exp(-lambda)))
  ifelse(hu < theta, 0, rpois(ndraws, lambda = lambda - t) + 1)
}

predict_hurdle_negbinomial <- function(i, draws, ...) {
  # theta is the bernoulli hurdle parameter
  theta <- get_dpar(draws, "hu", i = i)
  mu <- get_dpar(draws, "mu", i = i)
  ndraws <- draws$nsamples
  # compare with theta to incorporate the hurdle process
  hu <- runif(ndraws, 0, 1)
  # sample from an approximate(!) truncated negbinomial distribution
  # by adjusting mu and adding 1
  t = -log(1 - runif(ndraws) * (1 - exp(-mu)))
  shape <- get_dpar(draws, "shape", i = i)
  ifelse(hu < theta, 0, rnbinom(ndraws, mu = mu - t, size = shape) + 1)
}

predict_hurdle_gamma <- function(i, draws, ...) {
  # theta is the bernoulli hurdle parameter
  theta <- get_dpar(draws, "hu", i = i)
  shape <- get_dpar(draws, "shape", i = i)
  scale <- get_dpar(draws, "mu", i = i) / shape
  ndraws <- draws$nsamples
  # compare with theta to incorporate the hurdle process
  hu <- runif(ndraws, 0, 1)
  ifelse(hu < theta, 0, rgamma(ndraws, shape = shape, scale = scale))
}

predict_hurdle_lognormal <- function(i, draws, ...) {
  # theta is the bernoulli hurdle parameter
  theta <- get_dpar(draws, "hu", i = i)
  mu <- get_dpar(draws, "mu", i = i)
  sigma <- get_dpar(draws, "sigma", i = i)
  ndraws <- draws$nsamples
  # compare with theta to incorporate the hurdle process
  hu <- runif(ndraws, 0, 1)
  ifelse(hu < theta, 0, rlnorm(ndraws, meanlog = mu, sdlog = sigma))
}

predict_zero_inflated_beta <- function(i, draws, ...) {
  # theta is the bernoulli hurdle parameter
  theta <- get_dpar(draws, "zi", i = i)
  mu <- get_dpar(draws, "mu", i = i)
  phi <- get_dpar(draws, "phi", i = i)
  # compare with theta to incorporate the hurdle process
  hu <- runif(draws$nsamples, 0, 1)
  ifelse(
    hu < theta, 0, 
    rbeta(draws$nsamples, shape1 = mu * phi, shape2 = (1 - mu) * phi)
  )
}

predict_zero_one_inflated_beta <- function(i, draws, ...) {
  zoi <- get_dpar(draws, "zoi", i)
  coi <- get_dpar(draws, "coi", i)
  mu <- get_dpar(draws, "mu", i = i)
  phi <- get_dpar(draws, "phi", i = i)
  hu <- runif(draws$nsamples, 0, 1)
  one_or_zero <- runif(draws$nsamples, 0, 1)
  ifelse(hu < zoi, 
    ifelse(one_or_zero < coi, 1, 0),
    rbeta(draws$nsamples, shape1 = mu * phi, shape2 = (1 - mu) * phi)
  )
}

predict_zero_inflated_poisson <- function(i, draws, ...) {
  # theta is the bernoulli zero-inflation parameter
  theta <- get_dpar(draws, "zi", i = i)
  lambda <- get_dpar(draws, "mu", i = i)
  ndraws <- draws$nsamples
  # compare with theta to incorporate the zero-inflation process
  zi <- runif(ndraws, 0, 1)
  ifelse(zi < theta, 0, rpois(ndraws, lambda = lambda))
}

predict_zero_inflated_negbinomial <- function(i, draws, ...) {
  # theta is the bernoulli zero-inflation parameter
  theta <- get_dpar(draws, "zi", i = i)
  mu <- get_dpar(draws, "mu", i = i)
  shape <- get_dpar(draws, "shape", i = i)
  ndraws <- draws$nsamples
  # compare with theta to incorporate the zero-inflation process
  zi <- runif(ndraws, 0, 1)
  ifelse(zi < theta, 0, rnbinom(ndraws, mu = mu, size = shape))
}

predict_zero_inflated_binomial <- function(i, draws, ...) {
  # theta is the bernoulii zero-inflation parameter
  theta <- get_dpar(draws, "zi", i = i)
  trials <- draws$data$trials[i]
  prob <- get_dpar(draws, "mu", i = i)
  ndraws <- draws$nsamples
  # compare with theta to incorporate the zero-inflation process
  zi <- runif(ndraws, 0, 1)
  ifelse(zi < theta, 0, rbinom(ndraws, size = trials, prob = prob))
}

predict_categorical <- function(i, draws, ...) {
  eta <- sapply(names(draws$dpars), get_dpar, draws = draws, i = i)
  p <- pcategorical(seq_len(draws$data$ncat), eta = eta)
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
  ncat <- draws$data$ncat
  disc <- get_dpar(draws, "disc", i = i)
  eta <- (disc * get_dpar(draws, "mu", i = i))
  p <- pordinal(
    seq_len(ncat), eta = eta, ncat = ncat, 
    family = family, link = draws$f$link
  )
  first_greater(p, target = runif(draws$nsamples, min = 0, max = 1))
}

predict_mixture <- function(i, draws, ...) {
  families <- family_names(draws$f)
  theta <- get_theta(draws, i = i)
  smix <- rng_mix(theta)
  out <- rep(NA, draws$nsamples)
  for (j in seq_along(families)) {
    sample_ids <- which(smix == j)
    if (length(sample_ids)) {
      predict_fun <- paste0("predict_", families[j])
      predict_fun <- get(predict_fun, asNamespace("brms"))
      tmp_draws <- pseudo_draws_for_mixture(draws, j, sample_ids)
      out[sample_ids] <- predict_fun(i, tmp_draws, ...)
    }
  }
  out
}

#---------------predict helper-functions----------------------------

rng_continuous <- function(nrng, dist, args, lb = NULL, ub = NULL) {
  # random numbers from (possibly truncated) continuous distributions
  # Args:
  #   nrng: number of random values to generate
  #   dist: name of a distribution for which the functions
  #         p<dist>, q<dist>, and r<dist> are available
  #   args: dditional arguments passed to the distribution functions
  # Returns:
  #   a vector of random values draws from the distribution
  if (is.null(lb) && is.null(ub)) {
    # sample as usual
    rdist <- paste0("r", dist)
    out <- do.call(rdist, c(nrng, args))
  } else {
    # sample from truncated distribution
    if (is.null(lb)) lb <- -Inf
    if (is.null(ub)) ub <- Inf
    pdist <- paste0("p", dist)
    qdist <- paste0("q", dist)
    plb <- do.call(pdist, c(list(lb), args))
    pub <- do.call(pdist, c(list(ub), args))
    rng <- list(runif(nrng, min = plb, max = pub))
    out <- do.call(qdist, c(rng, args))
    # remove infinte values caused by numerical imprecision
    out[out %in% c(-Inf, Inf)] <- NA
  }
  out
}

rng_discrete <- function(nrng, dist, args, lb = NULL, ub = NULL, ntrys = 5) {
  # random numbers from (possibly truncated) discrete distributions
  # currently rejection sampling is used for truncated distributions
  # Args:
  #   nrng: number of random values to generate
  #   dist: name of a distribution for which the functions
  #         p<dist>, q<dist>, and r<dist> are available
  #   args: dditional arguments passed to the distribution functions
  #   lb: optional lower truncation bound
  #   ub: optional upper truncation bound
  #   ntrys: number of trys in rejection sampling for truncated models
  # Returns:
  #   a vector of random values draws from the distribution
  rdist <- get(paste0("r", dist), mode = "function")
  if (is.null(lb) && is.null(ub)) {
    # sample as usual
    do.call(rdist, c(nrng, args))
  } else {
    # sample from truncated distribution via rejection sampling
    if (is.null(lb)) lb <- -Inf
    if (is.null(ub)) ub <- Inf
    rng <- matrix(do.call(rdist, c(nrng * ntrys, args)), ncol = ntrys)
    apply(rng, 1, extract_valid_sample, lb = lb, ub = ub)
  }
}

rng_mix <- function(theta) {
  # sample the ID of the mixture component
  stopifnot(is.matrix(theta))
  mix_comp <- seq_len(ncol(theta))
  ulapply(seq_len(nrow(theta)), function(s)
    sample(mix_comp, 1, prob = theta[s, ])
  )
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

get_pct_invalid <- function(x, lb = NULL, ub = NULL) {
  # percentage of invalid predictions of truncated discrete models
  # Args:
  #   x: matrix of predicted values
  #   lb: optional lower truncation bound
  #   ub: optional upper truncation bound
  if (!(is.null(lb) && is.null(ub))) {
    if (is.null(lb)) lb <- -Inf
    if (is.null(ub)) ub <- Inf
    # ensures correct comparison with vector bounds
    x <- c(t(x))
    sum(x <= lb | x > ub, na.rm = TRUE) / length(x[!is.na(x)]) 
  } else {
    0
  }
}
