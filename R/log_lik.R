log_lik_internal <- function(draws, ...) {
  # for use in non-pointwise evaluation only
  UseMethod("log_lik_internal")
}

#' @export
log_lik_internal.mvbrmsdraws <- function(draws, combine = TRUE, ...) {
  if (length(draws$mvpars$rescor)) {
    draws$mvpars$Mu <- get_Mu(draws)
    draws$mvpars$Sigma <- get_Sigma(draws)
    out <- log_lik_internal.brmsdraws(draws, ...)
  } else {
    out <- lapply(draws$resps, log_lik_internal, ...)
    if (combine) {
      out <- Reduce("+", out)
    } else {
      along <- ifelse(length(out) > 1L, 3, 2)
      out <- do.call(abind, c(out, along = along))
    }
  }
  out
}

#' @export
log_lik_internal.brmsdraws <- function(draws, ...) {
  log_lik_fun <- paste0("log_lik_", draws$f$fun)
  log_lik_fun <- get(log_lik_fun, asNamespace("brms"))
  for (dp in names(draws$dpars)) {
    draws$dpars[[dp]] <- get_dpar(draws, dpar = dp)
  }
  N <- choose_N(draws)
  log_lik <- do.call(cbind, lapply(seq_len(N), log_lik_fun, draws = draws))
  old_order <- draws$old_order
  sort <- isTRUE(ncol(log_lik) != length(old_order))
  reorder_obs(log_lik, old_order, sort = sort)
}

log_lik_pointwise <- function(data_i, draws) {
  # for use in pointwise evaluation only
  # cannot be made an S3 methods since i must be the first argument
  i <- data_i$i
  if (is.mvbrmsdraws(draws) && !length(draws$mvpars$rescor)) {
    out <- lapply(draws$resps, log_lik_pointwise, i = i)
    out <- Reduce("+", out)
  } else {
    log_lik_fun <- paste0("log_lik_", draws$f$fun)
    log_lik_fun <- get(log_lik_fun, asNamespace("brms"))
    out <- log_lik_fun(i, draws)
  }
  out
}

# All log_lik_<family> functions have the same arguments structure
# Args:
#  i: the column of draws to use i.e. the ith obervation 
#     in the initial data.frame 
#  draws: A named list returned by extract_draws containing 
#         all required data and samples
#  data: ignored; included for compatibility with loo::loo.function      
# Returns:
#   A vector of length draws$nsamples containing the pointwise 
#   log-likelihood fo the ith observation 
log_lik_gaussian <- function(i, draws, data = data.frame()) {
  args <- list(
    mean = get_dpar(draws, "mu", i = i), 
    sd = get_dpar(draws, "sigma", i = i)
  )
  # log_lik_censor computes the conventional log_lik in case of no censoring 
  out <- log_lik_censor(dist = "norm", args = args, i = i, data = draws$data)
  out <- log_lik_truncate(
    out, cdf = pnorm, args = args, i = i, data = draws$data
  )
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_student <- function(i, draws, data = data.frame()) {
  args <- list(
    df = get_dpar(draws, "nu", i = i), 
    mu = get_dpar(draws, "mu", i = i), 
    sigma = get_dpar(draws, "sigma", i = i)
  )
  out <- log_lik_censor(
    dist = "student_t", args = args, i = i, data = draws$data
  )
  out <- log_lik_truncate(
    out, cdf = pstudent_t, args = args, i = i, data = draws$data
  )
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_lognormal <- function(i, draws, data = data.frame()) {
  sigma <- get_dpar(draws, "sigma", i = i)
  args <- list(meanlog = get_dpar(draws, "mu", i), sdlog = sigma)
  out <- log_lik_censor(dist = "lnorm", args = args, i = i, data = draws$data)
  out <- log_lik_truncate(
    out, cdf = plnorm, args = args, i = i, data = draws$data
  )
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_shifted_lognormal <- function(i, draws, data = data.frame()) {
  sigma <- get_dpar(draws, "sigma", i = i)
  ndt <- get_dpar(draws, "ndt", i = i)
  args <- list(meanlog = get_dpar(draws, "mu", i), sdlog = sigma, shift = ndt)
  out <- log_lik_censor("shifted_lnorm", args, i = i, data = draws$data)
  out <- log_lik_truncate(out, pshifted_lnorm, args, i = i, data = draws$data)
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_skew_normal <- function(i, draws, data = data.frame()) {
  sigma <- get_dpar(draws, "sigma", i = i)
  alpha <- get_dpar(draws, "alpha", i = i)
  mu <- get_dpar(draws, "mu", i)
  args <- nlist(mu, sigma, alpha)
  out <- log_lik_censor(
    dist = "skew_normal", args = args, i = i, data = draws$data
  )
  out <- log_lik_truncate(
    out, cdf = pskew_normal, args = args, i = i, data = draws$data
  )
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_gaussian_mv <- function(i, draws, data = data.frame()) {
  Mu <- get_Mu(draws, i = i)
  Sigma <- get_Sigma(draws, i = i)
  dmn <- function(s) {
    dmulti_normal(
      draws$data$Y[i, ], mu = Mu[s, ],
      Sigma = Sigma[s, , ], log = TRUE
    )
  }
  out <- sapply(1:draws$nsamples, dmn)
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_student_mv <- function(i, draws, data = data.frame()) {
  nu <- get_dpar(draws, "nu", i = i)
  Mu <- get_Mu(draws, i = i)
  Sigma <- get_Sigma(draws, i = i)
  dmst <- function(s) {
    dmulti_student_t(
      draws$data$Y[i, ], df = nu[s], mu = Mu[s, ],
      Sigma = Sigma[s, , ], log = TRUE
    )
  }
  out <- sapply(1:draws$nsamples, dmst)
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_gaussian_cov <- function(i, draws, data = data.frame()) {
  # currently, only ARMA1 processes are implemented
  obs <- with(draws$ac, begin_tg[i]:(begin_tg[i] + nobs_tg[i] - 1))
  args <- list(
    sigma = get_dpar(draws, "sigma", obs),
    se = draws$data$se[obs], nrows = length(obs)
  )
  if (!is.null(draws$ac$ar) && is.null(draws$ac$ma)) {
    args$ar <- draws$ac$ar
    Sigma <- do.call(get_cov_matrix_ar1, args)
  } else if (is.null(draws$ac$ar) && !is.null(draws$ac$ma)) {
    args$ma <- draws$ac$ma
    Sigma <- do.call(get_cov_matrix_ma1, args)
  } else {
    args[c("ar", "ma")] <- draws$ac[c("ar", "ma")]
    Sigma <- do.call(get_cov_matrix_arma1, args)
  }
  # make sure par[s, ] is valid even if obs is of length 1
  mu <- as.matrix(get_dpar(draws, "mu", obs))
  # weights, truncation and censoring not yet allowed
  sapply(1:draws$nsamples, function(s)
    dmulti_normal(
      draws$data$Y[obs], mu = mu[s, ], 
      Sigma = Sigma[s, , ], log = TRUE
    )
  )
}

log_lik_student_cov <- function(i, draws, data = data.frame()) {
  # currently, only ARMA1 processes are implemented
  obs <- with(draws$ac, begin_tg[i]:(begin_tg[i] + nobs_tg[i] - 1))
  args <- list(
    sigma = get_dpar(draws, "sigma", obs),
    se = draws$data$se[obs], nrows = length(obs)
  )
  if (!is.null(draws$ac$ar) && is.null(draws$ac$ma)) {
    args$ar <- draws$ac$ar
    Sigma <- do.call(get_cov_matrix_ar1, args)
  } else if (is.null(draws$ac$ar) && !is.null(draws$ac$ma)) {
    args$ma <- draws$ac$ma
    Sigma <- do.call(get_cov_matrix_ma1, args)
  } else {
    args[c("ar", "ma")] <- draws$ac[c("ar", "ma")]
    Sigma <- do.call(get_cov_matrix_arma1, args)
  }
  # make sure par[s, ] is valid even if obs is of length 1
  mu <- as.matrix(get_dpar(draws, "mu", obs))
  nu <- as.matrix(get_dpar(draws, "nu", obs))
  # weights, truncation and censoring not yet allowed
  sapply(1:draws$nsamples, function(s)
    dmulti_student_t(
      draws$data$Y[obs], df = nu[s, ], mu = mu[s, ], 
      Sigma = Sigma[s, , ], log = TRUE
    )
  )
}

log_lik_gaussian_lagsar <- function(i, draws, data = data.frame()) {
  stopifnot(i == 1)
  .log_lik_gaussian_lagsar <- function(s) {
    W_new <- with(draws, diag(nobs) - ac$lagsar[s, ] * ac$W)
    mu <- as.numeric(solve(W_new) %*% mu[s, ])
    Sigma <- solve(crossprod(W_new)) * sigma[s]^2
    dmulti_normal(draws$data$Y, mu = mu, Sigma = Sigma, log = TRUE)
  }
  mu <- get_dpar(draws, "mu")
  sigma <- get_dpar(draws, "sigma")
  # weights, truncation and censoring not yet allowed
  sapply(1:draws$nsamples, .log_lik_gaussian_lagsar)
}

log_lik_student_lagsar <- function(i, draws, data = data.frame()) {
  stopifnot(i == 1)
  .log_lik_student_lagsar <- function(s) {
    W_new <- with(draws, diag(nobs) - ac$lagsar[s, ] * ac$W)
    mu <- as.numeric(solve(W_new) %*% mu[s, ])
    Sigma <- solve(crossprod(W_new)) * sigma[s]^2
    dmulti_student_t(
      draws$data$Y, df = nu[s], mu = mu, 
      Sigma = Sigma, log = TRUE
    )
  }
  mu <- get_dpar(draws, "mu")
  sigma <- get_dpar(draws, "sigma")
  nu <- get_dpar(draws, "nu")
  # weights, truncation and censoring not yet allowed
  sapply(1:draws$nsamples, .log_lik_student_lagsar)
}

log_lik_gaussian_errorsar <- function(i, draws, data = data.frame()) {
  stopifnot(i == 1)
  .log_lik_gaussian_errorsar <- function(s) {
    W_new <- with(draws, diag(nobs) - ac$errorsar[s, ] * ac$W)
    Sigma <- solve(crossprod(W_new)) * sigma[s]^2
    dmulti_normal(draws$data$Y, mu = mu[s, ], Sigma = Sigma, log = TRUE)
  }
  mu <- get_dpar(draws, "mu")
  sigma <- get_dpar(draws, "sigma")
  # weights, truncation and censoring not yet allowed
  sapply(1:draws$nsamples, .log_lik_gaussian_errorsar)
}

log_lik_student_errorsar <- function(i, draws, data = data.frame()) {
  stopifnot(i == 1)
  .log_lik_student_errorsar <- function(s) {
    W_new <- with(draws, diag(nobs) - ac$errorsar[s, ] * ac$W)
    Sigma <- solve(crossprod(W_new)) * sigma[s]^2
    dmulti_student_t(
      draws$data$Y, df = nu[s], mu = mu[s, ], 
      Sigma = Sigma, log = TRUE
    )
  }
  mu <- get_dpar(draws, "mu")
  sigma <- get_dpar(draws, "sigma")
  nu <- get_dpar(draws, "nu")
  # weights, truncation and censoring not yet allowed
  sapply(1:draws$nsamples, .log_lik_student_errorsar)
}

log_lik_gaussian_fixed <- function(i, draws, data = data.frame()) {
  stopifnot(i == 1)
  mu <- as.matrix(get_dpar(draws, "mu"))
  ulapply(1:draws$nsamples, function(s) 
    dmulti_normal(
      draws$data$Y, mu = mu[s, ], 
      Sigma = draws$ac$V, log = TRUE
    )
  )
}

log_lik_student_fixed <- function(i, draws, data = data.frame()) {
  stopifnot(i == 1)
  mu <- as.matrix(get_dpar(draws, "mu"))
  nu <- as.matrix(get_dpar(draws, "nu"))
  sapply(1:draws$nsamples, function(s) 
    dmulti_student_t(
      draws$data$Y, df = nu[s, ], mu = mu[s, ],
      Sigma = draws$ac$V, log = TRUE
    )
  )
}

log_lik_binomial <- function(i, draws, data = data.frame()) {
  trials <- draws$data$trials[i]
  args <- list(size = trials, prob = get_dpar(draws, "mu", i))
  out <- log_lik_censor(
    dist = "binom", args = args, i = i, data = draws$data
  )
  out <- log_lik_truncate(
    out, cdf = pbinom, args = args, i = i, data = draws$data
  )
  log_lik_weight(out, i = i, data = draws$data)
}  

log_lik_bernoulli <- function(i, draws, data = data.frame()) {
  args <- list(size = 1, prob = get_dpar(draws, "mu", i))
  out <- log_lik_censor(
    dist = "binom", args = args, i = i, data = draws$data
  )
  # no truncation allowed
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_poisson <- function(i, draws, data = data.frame()) {
  args <- list(lambda = get_dpar(draws, "mu", i))
  out <- log_lik_censor(
    dist = "pois", args = args, i = i, data = draws$data
  )
  out <- log_lik_truncate(
    out, cdf = ppois, args = args, i = i, data = draws$data
  )
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_negbinomial <- function(i, draws, data = data.frame()) {
  shape <- get_dpar(draws, "shape", i = i)
  args <- list(mu = get_dpar(draws, "mu", i), size = shape)
  out <- log_lik_censor(
    dist = "nbinom", args = args, i = i, data = draws$data
  )
  out <- log_lik_truncate(
    out, cdf = pnbinom, args = args, i = i, data = draws$data
  )
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_geometric <- function(i, draws, data = data.frame()) {
  args <- list(mu = get_dpar(draws, "mu", i), size = 1)
  out <- log_lik_censor(
    dist = "nbinom", args = args, i = i, data = draws$data
  )
  out <- log_lik_truncate(
    out, cdf = pnbinom, args = args, i = i, data = draws$data
  )
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_exponential <- function(i, draws, data = data.frame()) {
  args <- list(rate = 1 / get_dpar(draws, "mu", i))
  out <- log_lik_censor(dist = "exp", args = args, i = i, data = draws$data)
  out <- log_lik_truncate(
    out, cdf = pexp, args = args, i = i, data = draws$data
  )
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_gamma <- function(i, draws, data = data.frame()) {
  shape <- get_dpar(draws, "shape", i = i)
  args <- list(
    shape = shape, 
    scale = get_dpar(draws, "mu", i) / shape
  )
  out <- log_lik_censor(dist = "gamma", args = args, i = i, data = draws$data)
  out <- log_lik_truncate(
    out, cdf = pgamma, args = args, i = i, data = draws$data
  )
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_weibull <- function(i, draws, data = data.frame()) {
  shape <- get_dpar(draws, "shape", i = i)
  scale <- get_dpar(draws, "mu", i = i) / gamma(1 + 1 / shape)
  args <- list(shape = shape, scale = scale)
  out <- log_lik_censor(
    dist = "weibull", args = args, i = i, data = draws$data
  )
  out <- log_lik_truncate(
    out, cdf = pweibull, args = args, i = i, data = draws$data
  )
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_frechet <- function(i, draws, data = data.frame()) {
  nu <- get_dpar(draws, "nu", i = i)
  scale <- get_dpar(draws, "mu", i = i) / gamma(1 - 1 / nu)
  args <- list(scale = scale, shape = nu)
  out <- log_lik_censor(
    dist = "frechet", args = args, i = i, data = draws$data
  )
  out <- log_lik_truncate(
    out, cdf = pfrechet, args = args, i = i, data = draws$data
  )
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_gen_extreme_value <- function(i, draws, data = data.frame()) {
  sigma <- get_dpar(draws, "sigma", i = i)
  xi <- get_dpar(draws, "xi", i = i)
  mu <- get_dpar(draws, "mu", i)
  args <- nlist(mu, sigma, xi)
  out <- log_lik_censor(dist = "gen_extreme_value", args = args, 
                       i = i, data = draws$data)
  out <- log_lik_truncate(out, cdf = pgen_extreme_value, 
                         args = args, i = i, data = draws$data)
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_inverse.gaussian <- function(i, draws, data = data.frame()) {
  args <- list(mu = get_dpar(draws, "mu", i), 
               shape = get_dpar(draws, "shape", i = i))
  out <- log_lik_censor(dist = "inv_gaussian", args = args, 
                       i = i, data = draws$data)
  out <- log_lik_truncate(out, cdf = pinv_gaussian, args = args,
                         i = i, data = draws$data)
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_exgaussian <- function(i, draws, data = data.frame()) {
  args <- list(mu = get_dpar(draws, "mu", i), 
               sigma = get_dpar(draws, "sigma", i = i),
               beta = get_dpar(draws, "beta", i = i))
  out <- log_lik_censor(dist = "exgaussian", args = args, 
                       i = i, data = draws$data)
  out <- log_lik_truncate(out, cdf = pexgaussian, args = args,
                         i = i, data = draws$data)
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_wiener <- function(i, draws, data = data.frame()) {
  args <- list(
    delta = get_dpar(draws, "mu", i), 
    alpha = get_dpar(draws, "bs", i = i),
    tau = get_dpar(draws, "ndt", i = i),
    beta = get_dpar(draws, "bias", i = i),
    resp = draws$data[["dec"]][i]
  )
  out <- do.call("dwiener", c(draws$data$Y[i], args, log = TRUE))
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_beta <- function(i, draws, data = data.frame()) {
  mu <- get_dpar(draws, "mu", i)
  phi <- get_dpar(draws, "phi", i)
  args <- list(shape1 = mu * phi, shape2 = (1 - mu) * phi)
  out <- log_lik_censor(dist = "beta", args = args, i = i, data = draws$data)
  out <- log_lik_truncate(
    out, cdf = pbeta, args = args, i = i, data = draws$data
  )
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_von_mises <- function(i, draws, data = data.frame()) {
  args <- list(
    mu = get_dpar(draws, "mu", i), 
    kappa = get_dpar(draws, "kappa", i = i)
  )
  out <- log_lik_censor(
    dist = "von_mises", args = args, i = i, data = draws$data
  )
  out <- log_lik_truncate(
    out, cdf = pvon_mises, args = args, i = i, data = draws$data
  )
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_asym_laplace <- function(i, draws, ...) {
  args <- list(mu = get_dpar(draws, "mu", i), 
               sigma = get_dpar(draws, "sigma", i = i),
               quantile = get_dpar(draws, "quantile", i = i))
  out <- log_lik_censor(dist = "asym_laplace", args = args, 
                       i = i, data = draws$data)
  out <- log_lik_truncate(out, cdf = pvon_mises, args = args,
                         i = i, data = draws$data)
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_hurdle_poisson <- function(i, draws, data = data.frame()) {
  theta <- get_dpar(draws, "hu", i)
  args <- list(lambda = get_dpar(draws, "mu", i))
  out <- log_lik_hurdle_discrete(pdf = dpois, theta = theta, 
                                args = args, i = i, data = draws$data)
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_hurdle_negbinomial <- function(i, draws, data = data.frame()) {
  theta <- get_dpar(draws, "hu", i)
  args <- list(mu = get_dpar(draws, "mu", i), 
               size = get_dpar(draws, "shape", i = i))
  out <- log_lik_hurdle_discrete(pdf = dnbinom, theta = theta, 
                                args = args, i = i, data = draws$data)
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_hurdle_gamma <- function(i, draws, data = data.frame()) {
  theta <- get_dpar(draws, "hu", i)
  shape <- get_dpar(draws, "shape", i = i)
  args <- list(
    shape = shape, 
    scale = get_dpar(draws, "mu", i) / shape
  )
  out <- log_lik_hurdle_continuous(
    pdf = dgamma, theta = theta, args = args, i = i, data = draws$data
  )
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_hurdle_lognormal <- function(i, draws, data = data.frame()) {
  theta <- get_dpar(draws, "hu", i)
  sigma <- get_dpar(draws, "sigma", i = i)
  args <- list(meanlog = get_dpar(draws, "mu", i), sdlog = sigma)
  out <- log_lik_hurdle_continuous(pdf = dlnorm, theta = theta, 
                                  args = args, i = i, data = draws$data)
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_zero_inflated_poisson <- function(i, draws, data = data.frame()) {
  theta <- get_dpar(draws, "zi", i)
  args <- list(lambda = get_dpar(draws, "mu", i))
  out <- log_lik_zero_inflated(pdf = dpois, theta = theta, 
                              args = args, i = i, data = draws$data)
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_zero_inflated_negbinomial <- function(i, draws, data = data.frame()) {
  theta <- get_dpar(draws, "zi", i)
  args <- list(
    mu = get_dpar(draws, "mu", i), 
    size = get_dpar(draws, "shape", i = i)
  )
  out <- log_lik_zero_inflated(
    pdf = dnbinom, theta = theta, args = args, i = i, data = draws$data
  )
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_zero_inflated_binomial <- function(i, draws, data = data.frame()) {
  trials <- draws$data$trials[i] 
  theta <- get_dpar(draws, "zi", i)
  args <- list(size = trials, prob = get_dpar(draws, "mu", i))
  out <- log_lik_zero_inflated(pdf = dbinom, theta = theta, 
                              args = args, i = i, data = draws$data)
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_zero_inflated_beta <- function(i, draws, data = data.frame()) {
  theta <- get_dpar(draws, "zi", i)
  mu <- get_dpar(draws, "mu", i)
  phi <- get_dpar(draws, "phi", i)
  args <- list(shape1 = mu * phi, shape2 = (1 - mu) * phi)
  # zi_beta is technically a hurdle model
  out <- log_lik_hurdle_continuous(
    pdf = dbeta, theta = theta, args = args, i = i, data = draws$data)
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_zero_one_inflated_beta <- function(i, draws, data = data.frame()) {
  zoi <- get_dpar(draws, "zoi", i)
  coi <- get_dpar(draws, "coi", i)
  if (draws$data$Y[i] %in% c(0, 1)) {
    out <- dbinom(1, size = 1, prob = zoi, log = TRUE) + 
      dbinom(draws$data$Y[i], size = 1, prob = coi, log = TRUE)
  } else {
    phi <- get_dpar(draws, "phi", i)
    mu <- get_dpar(draws, "mu", i)
    args <- list(shape1 = mu * phi, shape2 = (1 - mu) * phi)
    out <- dbinom(0, size = 1, prob = zoi, log = TRUE) + 
      do.call(dbeta, c(draws$data$Y[i], args, log = TRUE))
  }
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_categorical <- function(i, draws, data = data.frame()) {
  stopifnot(draws$f$link == "logit")
  eta <- sapply(names(draws$dpars), get_dpar, draws = draws, i = i)
  out <- dcategorical(draws$data$Y[i], eta = eta, log = TRUE)
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_cumulative <- function(i, draws, data = data.frame()) {
  ncat <- draws$data$ncat
  eta <- get_dpar(draws, "disc", i = i) * get_dpar(draws, "mu", i = i)
  y <- draws$data$Y[i]
  if (y == 1) { 
    out <- log(ilink(eta[, 1], draws$f$link))
  } else if (y == ncat) {
    out <- log(1 - ilink(eta[, y - 1], draws$f$link)) 
  } else {
    out <- log(
      ilink(eta[, y], draws$f$link) - ilink(eta[, y - 1], draws$f$link)
    )
  }
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_sratio <- function(i, draws, data = data.frame()) {
  ncat <- draws$data$ncat
  eta <- get_dpar(draws, "disc", i = i) * get_dpar(draws, "mu", i = i)
  y <- draws$data$Y[i]
  q <- sapply(1:min(y, ncat - 1), 
    function(k) 1 - ilink(eta[, k], draws$f$link)
  )
  if (y == 1) {
    out <- log(1 - q[, 1]) 
  } else if (y == 2) {
    out <- log(1 - q[, 2]) + log(q[, 1])
  } else if (y == ncat) {
    out <- rowSums(log(q))
  } else {
    out <- log(1 - q[, y]) + rowSums(log(q[, 1:(y - 1)]))
  }
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_cratio <- function(i, draws, data = data.frame()) {
  ncat <- draws$data$ncat
  eta <- get_dpar(draws, "disc", i = i) * get_dpar(draws, "mu", i = i)
  y <- draws$data$Y[i]
  q <- sapply(1:min(y, ncat-1), function(k) ilink(eta[, k], draws$f$link))
  if (y == 1) {
    out <- log(1 - q[, 1])
  }  else if (y == 2) {
    out <- log(1 - q[, 2]) + log(q[, 1])
  } else if (y == ncat) {
    out <- rowSums(log(q))
  } else {
    out <- log(1 - q[, y]) + rowSums(log(q[, 1:(y - 1)]))
  }
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_acat <- function(i, draws, data = data.frame()) {
  ncat <- draws$data$ncat
  eta <- get_dpar(draws, "disc", i = i) * get_dpar(draws, "mu", i = i)
  y <- draws$data$Y[i]
  if (draws$f$link == "logit") { # more efficient calculation 
    q <- sapply(1:(ncat - 1), function(k) eta[, k])
    p <- cbind(rep(0, nrow(eta)), q[, 1], 
               matrix(0, nrow = nrow(eta), ncol = ncat - 2))
    if (ncat > 2) {
      p[, 3:ncat] <- sapply(3:ncat, function(k) rowSums(q[, 1:(k - 1)]))
    }
    out <- p[, y] - log(rowSums(exp(p)))
  } else {
    q <- sapply(1:(ncat - 1), function(k) 
      ilink(eta[, k], draws$f$link))
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
  log_lik_weight(out, i = i, data = draws$data)
}

log_lik_custom <- function(i, draws, data = data.frame()) {
  log_lik_fun <- paste0("log_lik_", draws$f$name)
  log_lik_fun <- get(log_lik_fun, draws$f$env)
  log_lik_fun(i = i, draws = draws)
}

log_lik_mixture <- function(i, draws, data = data.frame()) {
  families <- family_names(draws$f)
  theta <- get_theta(draws, i = i)
  out <- array(NA, dim = dim(theta))
  for (j in seq_along(families)) {
    log_lik_fun <- paste0("log_lik_", families[j])
    log_lik_fun <- get(log_lik_fun, asNamespace("brms"))
    tmp_draws <- pseudo_draws_for_mixture(draws, j)
    out[, j] <- exp(log(theta[, j]) + log_lik_fun(i, tmp_draws))
  }
  if (isTRUE(draws[["pp_mixture"]])) {
    out <- log(out) - log(rowSums(out))
  } else {
    out <- log(rowSums(out))
  }
  out
}

# ----------- log_lik helper-functions -----------

log_lik_censor <- function(dist, args, i, data) {
  # compute (possibly censored) log_lik values
  # Args:
  #   dist: name of a distribution for which the functions
  #         d<dist> (pdf) and p<dist> (cdf) are available
  #   args: additional arguments passed to pdf and cdf
  #   data: data initially passed to Stan
  # Returns:
  #   vector of log_lik values
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

log_lik_truncate <- function(x, cdf, args, i, data) {
  # adjust log_lik in truncated models
  # Args:
  #   x: vector of log_lik values
  #   cdf: a cumulative distribution function 
  #   args: arguments passed to cdf
  #   i: observation number
  #   data: data initally passed to Stan
  # Returns:
  #   vector of log_lik values
  if (!(is.null(data$lb) && is.null(data$ub))) {
    lb <- if (is.null(data$lb)) -Inf else data$lb[i]
    ub <- if (is.null(data$ub)) -Inf else data$ub[i]
    x - log(do.call(cdf, c(ub, args)) - do.call(cdf, c(lb, args)))
  } else {
    x
  }
}

log_lik_weight <- function(x, i, data) {
  # weight log_lik values according to defined weights
  # Args:
  #   x: vector of log_lik values
  #   i: observation number
  #   data: data initially passed to Stan
  # Returns:
  #   vector of log_lik values
  if ("weights" %in% names(data)) {
    x * data$weights[i]
  } else {
    x
  }
}

log_lik_hurdle_discrete <- function(pdf, theta, args, i, data) {
  # log_lik values for discrete hurdle models
  # Args:
  #  pdf: a probability density function 
  #  theta: bernoulli hurdle parameter
  #  args: arguments passed to pdf
  #  data: data initially passed to Stan
  # Returns:
  #   vector of log_lik values
  if (data$Y[i] == 0) {
    dbinom(1, size = 1, prob = theta, log = TRUE)
  } else {
    dbinom(0, size = 1, prob = theta, log = TRUE) + 
      do.call(pdf, c(data$Y[i], args, log = TRUE)) -
      log(1 - do.call(pdf, c(0, args)))
  }
}

log_lik_hurdle_continuous <- function(pdf, theta, args, i, data) {
  # log_lik values for continuous hurdle models
  # does not call log(1 - do.call(pdf, c(0, args)))
  # Args:
  #   same as log_lik_hurdle_discrete
  if (data$Y[i] == 0) {
    dbinom(1, size = 1, prob = theta, log = TRUE)
  } else {
    dbinom(0, size = 1, prob = theta, log = TRUE) + 
      do.call(pdf, c(data$Y[i], args, log = TRUE))
  }
}

log_lik_zero_inflated <- function(pdf, theta, args, i, data) {
  # log_lik values for zero-inflated models
  # Args:
  #  pdf: a probability density function 
  #  theta: bernoulli zero-inflation parameter
  #  args: arguments passed to pdf
  #  data: data initially passed to Stan
  # Returns:
  #   vector of log_lik values
  if (data$Y[i] == 0) {
    log(dbinom(1, size = 1, prob = theta) + 
        dbinom(0, size = 1, prob = theta) *
          do.call(pdf, c(0, args)))
  } else {
    dbinom(0, size = 1, prob = theta, log = TRUE) +
      do.call(pdf, c(data$Y[i], args, log = TRUE))
  }
}
