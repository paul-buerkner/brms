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
    out <- do_call(abind, c(out, along = along))
  }
  out
}

#' @export
predict_internal.brmsdraws <- function(draws, summary = TRUE, transform = NULL,
                                       sort = FALSE, robust = FALSE, 
                                       probs = c(0.025, 0.975), ...) {
  for (nlp in names(draws$nlpars)) {
    draws$nlpars[[nlp]] <- get_nlpar(draws, nlpar = nlp)
  }
  for (dp in names(draws$dpars)) {
    draws$dpars[[dp]] <- get_dpar(draws, dpar = dp)
  }
  predict_fun <- paste0("predict_", draws$family$fun)
  predict_fun <- get(predict_fun, asNamespace("brms"))
  N <- choose_N(draws)
  out <- lapply(seq_len(N), predict_fun, draws = draws, ...)
  if (grepl("_mv$", draws$family$fun)) {
    out <- do_call(abind, c(out, along = 3))
    out <- aperm(out, perm = c(1, 3, 2))
    dimnames(out)[[3]] <- names(draws$resps)
  } else if (has_multicol(draws$family)) {
    out <- do_call(abind, c(out, along = 3))
    out <- aperm(out, perm = c(1, 3, 2))
    dimnames(out)[[3]] <- draws$data$cats
  } else {
    out <- do_call(cbind, out) 
  }
  colnames(out) <- NULL
  if (use_int(draws$family)) {
    out <- check_discrete_trunc_bounds(
      out, lb = draws$data$lb, ub = draws$data$ub
    )
  }
  out <- reorder_obs(out, draws$old_order, sort = sort)
  # transform predicted response samples before summarizing them 
  if (!is.null(transform)) {
    out <- do_call(transform, list(out))
  }
  attr(out, "levels") <- draws$data$cats
  if (summary) {
    if (is_ordinal(draws$family)) {
      levels <- seq_len(max(draws$data$nthres) + 1)
      out <- posterior_table(out, levels = levels)
    } else if (is_categorical(draws$family)) {
      levels <- seq_len(draws$data$ncat)
      out <- posterior_table(out, levels = levels)
    } else {
      out <- posterior_summary(out, probs = probs, robust = robust)
    }
  }
  out
}

# All predict_<family> functions have the same arguments structure
# @param i the column of draws to use that is the ith obervation 
#   in the initial data.frame 
# @param draws A named list returned by extract_draws containing 
#   all required data and samples
# @param ... ignored arguments
# @param A vector of length draws$nsamples containing samples
#   from the posterior predictive distribution
predict_gaussian <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "norm",
    mean = get_dpar(draws, "mu", i = i), 
    sd = get_dpar(draws, "sigma", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_student <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "student_t", 
    df = get_dpar(draws, "nu", i = i), 
    mu = get_dpar(draws, "mu", i = i), 
    sigma = get_dpar(draws, "sigma", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_lognormal <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "lnorm",
    meanlog = get_dpar(draws, "mu", i = i), 
    sdlog = get_dpar(draws, "sigma", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_shifted_lognormal <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "shifted_lnorm",
    meanlog = get_dpar(draws, "mu", i = i), 
    sdlog = get_dpar(draws, "sigma", i = i),
    shift = get_dpar(draws, "ndt", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_skew_normal <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "skew_normal",
    mu = get_dpar(draws, "mu", i = i),
    sigma = get_dpar(draws, "sigma", i = i),
    alpha = get_dpar(draws, "alpha", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_gaussian_mv <- function(i, draws, ...) {
  Mu <- get_Mu(draws, i = i)
  Sigma <- get_Sigma(draws, i = i)
  .predict <- function(s) {
    rmulti_normal(1, mu = Mu[s, ], Sigma = Sigma[s, , ])
  }
  rblapply(seq_len(draws$nsamples), .predict)
}

predict_student_mv <- function(i, draws, ...) {
  nu <- get_dpar(draws, "nu", i = i)
  Mu <- get_Mu(draws, i = i)
  Sigma <- get_Sigma(draws, i = i)
  .predict <- function(s) {
    rmulti_student_t(1, df = nu[s], mu = Mu[s, ], Sigma = Sigma[s, , ])
  }
  rblapply(seq_len(draws$nsamples), .predict)
}

predict_gaussian_cov <- function(i, draws, ...) {
  obs <- with(draws$ac, begin_tg[i]:end_tg[i])
  mu <- as.matrix(get_dpar(draws, "mu", i = obs))
  Sigma <- get_cov_matrix_autocor(draws, obs)
  .predict <- function(s) {
    rmulti_normal(1, mu = mu[s, ], Sigma = Sigma[s, , ])
  }
  rblapply(seq_len(draws$nsamples), .predict)
}

predict_student_cov <- function(i, draws, ...) {
  obs <- with(draws$ac, begin_tg[i]:end_tg[i])
  nu <- as.matrix(get_dpar(draws, "nu", i = obs))
  mu <- as.matrix(get_dpar(draws, "mu", i = obs))
  Sigma <- get_cov_matrix_autocor(draws, obs)
  .predict <- function(s) {
    rmulti_student_t(1, df = nu[s, ], mu = mu[s, ], Sigma = Sigma[s, , ])
  }
  rblapply(seq_len(draws$nsamples), .predict)
}

predict_gaussian_lagsar <- function(i, draws, ...) {
  stopifnot(i == 1)
  .predict <- function(s) {
    W_new <- with(draws, diag(nobs) - ac$lagsar[s] * ac$W)
    mu <- as.numeric(solve(W_new) %*% mu[s, ])
    Sigma <- solve(crossprod(W_new)) * sigma[s]^2
    rmulti_normal(1, mu = mu, Sigma = Sigma)
  }
  mu <- get_dpar(draws, "mu")
  sigma <- get_dpar(draws, "sigma")
  rblapply(seq_len(draws$nsamples), .predict)
}

predict_student_lagsar <- function(i, draws, ...) {
  stopifnot(i == 1)
  .predict <- function(s) {
    W_new <- with(draws, diag(nobs) - ac$lagsar[s] * ac$W)
    mu <- as.numeric(solve(W_new) %*% mu[s, ])
    Sigma <- solve(crossprod(W_new)) * sigma[s]^2
    rmulti_student_t(1, df = nu[s], mu = mu, Sigma = Sigma)
  }
  mu <- get_dpar(draws, "mu")
  sigma <- get_dpar(draws, "sigma")
  nu <- get_dpar(draws, "nu")
  rblapply(seq_len(draws$nsamples), .predict)
}

predict_gaussian_errorsar <- function(i, draws, ...) {
  stopifnot(i == 1)
  .predict <- function(s) {
    W_new <- with(draws, diag(nobs) - ac$errorsar[s] * ac$W)
    Sigma <- solve(crossprod(W_new)) * sigma[s]^2
    rmulti_normal(1, mu = mu[s, ], Sigma = Sigma)
  }
  mu <- get_dpar(draws, "mu")
  sigma <- get_dpar(draws, "sigma")
  rblapply(seq_len(draws$nsamples), .predict)
}

predict_student_errorsar <- function(i, draws, ...) {
  stopifnot(i == 1)
  .predict <- function(s) {
    W_new <- with(draws, diag(nobs) - ac$errorsar[s] * ac$W)
    Sigma <- solve(crossprod(W_new)) * sigma[s]^2
    rmulti_student_t(1, df = nu[s], mu = mu[s, ], Sigma = Sigma)
  }
  mu <- get_dpar(draws, "mu")
  sigma <- get_dpar(draws, "sigma")
  nu <- get_dpar(draws, "nu")
  rblapply(seq_len(draws$nsamples), .predict)
}

predict_gaussian_fixed <- function(i, draws, ...) {
  stopifnot(i == 1)
  mu <- as.matrix(get_dpar(draws, "mu"))
  .predict <- function(s) {
    rmulti_normal(1, mu = mu[s, ], Sigma = draws$ac$V)
  }
  rblapply(seq_len(draws$nsamples), .predict)
}

predict_student_fixed <- function(i, draws, ...) {
  stopifnot(i == 1)
  mu <- as.matrix(get_dpar(draws, "mu"))
  nu <- as.matrix(get_dpar(draws, "nu"))
  .predict <- function(s) {
    rmulti_student_t(1, df = nu[s, ], mu = mu[s, ], Sigma = draws$ac$V)
  }
  rblapply(seq_len(draws$nsamples), .predict)
}

predict_binomial <- function(i, draws, ntrys = 5, ...) {
  rdiscrete(
    n = draws$nsamples, dist = "binom", 
    size = draws$data$trials[i], 
    prob = get_dpar(draws, "mu", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i], 
    ntrys = ntrys
  )
}

predict_bernoulli <- function(i, draws, ...) {
  mu <- get_dpar(draws, "mu", i = i)
  rbinom(length(mu), size = 1, prob = mu)
}

predict_poisson <- function(i, draws, ntrys = 5, ...) {
  rdiscrete(
    n = draws$nsamples, dist = "pois",
    lambda = get_dpar(draws, "mu", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i],
    ntrys = ntrys
  )
}

predict_negbinomial <- function(i, draws, ntrys = 5, ...) {
  rdiscrete(
    n = draws$nsamples, dist = "nbinom",
    mu = get_dpar(draws, "mu", i = i), 
    size = get_dpar(draws, "shape", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i],
    ntrys = ntrys
  )
}

predict_geometric <- function(i, draws, ntrys = 5, ...) {
  rdiscrete(
    n = draws$nsamples, dist = "nbinom",
    mu = get_dpar(draws, "mu", i = i), size = 1,
    lb = draws$data$lb[i], ub = draws$data$ub[i], 
    ntrys = ntrys
  )
}

predict_discrete_weibull <- function(i, draws, ntrys = 5, ...) {
  rdiscrete(
    n = draws$nsamples, dist = "discrete_weibull",
    mu = get_dpar(draws, "mu", i = i), 
    shape = get_dpar(draws, "shape", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i],
    ntrys = ntrys
  )
}

predict_com_poisson <- function(i, draws, ntrys = 5, ...) {
  rdiscrete(
    n = draws$nsamples, dist = "com_poisson",
    mu = get_dpar(draws, "mu", i = i), 
    shape = get_dpar(draws, "shape", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i],
    ntrys = ntrys
  )
}

predict_exponential <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "exp",
    rate = 1 / get_dpar(draws, "mu", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_gamma <- function(i, draws, ...) {
  shape <- get_dpar(draws, "shape", i = i)
  scale <- get_dpar(draws, "mu", i = i) / shape
  rcontinuous(
    n = draws$nsamples, dist = "gamma",
    shape = shape, scale = scale,
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_weibull <- function(i, draws, ...) {
  shape <- get_dpar(draws, "shape", i = i)
  scale <- get_dpar(draws, "mu", i = i) / gamma(1 + 1 / shape) 
  rcontinuous(
    n = draws$nsamples, dist = "weibull",
    shape = shape, scale = scale,
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_frechet <- function(i, draws, ...) {
  nu <- get_dpar(draws, "nu", i = i)
  scale <- get_dpar(draws, "mu", i = i) / gamma(1 - 1 / nu)
  rcontinuous(
    n = draws$nsamples, dist = "frechet",
    scale = scale, shape = nu,
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_gen_extreme_value <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "gen_extreme_value", 
    sigma = get_dpar(draws, "sigma", i = i),
    xi = get_dpar(draws, "xi", i = i),
    mu = get_dpar(draws, "mu", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_inverse.gaussian <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "inv_gaussian",
    mu = get_dpar(draws, "mu", i = i), 
    shape = get_dpar(draws, "shape", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_exgaussian <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "exgaussian",
    mu = get_dpar(draws, "mu", i = i), 
    sigma = get_dpar(draws, "sigma", i = i),
    beta = get_dpar(draws, "beta", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_wiener <- function(i, draws, negative_rt = FALSE, ...) {
  out <- rcontinuous(
    n = 1, dist = "wiener", 
    delta = get_dpar(draws, "mu", i = i), 
    alpha = get_dpar(draws, "bs", i = i),
    tau = get_dpar(draws, "ndt", i = i),
    beta = get_dpar(draws, "bias", i = i),
    types = if (negative_rt) c("q", "resp") else "q",
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
  rcontinuous(
    n = draws$nsamples, dist = "beta", 
    shape1 = mu * phi, shape2 = (1 - mu) * phi,
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_von_mises <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "von_mises",
    mu = get_dpar(draws, "mu", i = i), 
    kappa = get_dpar(draws, "kappa", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_asym_laplace <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "asym_laplace",
    mu = get_dpar(draws, "mu", i = i), 
    sigma = get_dpar(draws, "sigma", i = i),
    quantile = get_dpar(draws, "quantile", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

predict_zero_inflated_asym_laplace <- function(i, draws, ...) {
  zi <- get_dpar(draws, "zi", i = i)
  tmp <- runif(draws$nsamples, 0, 1)
  ifelse(
    tmp < zi, 0, 
    rcontinuous(
      n = draws$nsamples, dist = "asym_laplace",
      mu = get_dpar(draws, "mu", i = i), 
      sigma = get_dpar(draws, "sigma", i = i),
      quantile = get_dpar(draws, "quantile", i = i),
      lb = draws$data$lb[i], ub = draws$data$ub[i]
    )
  )
}

predict_cox <- function(i, draws, ...) {
  stop2("Cannot sample from the posterior predictive ",
        "distribution for family 'cox'.")
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
  eta <- insert_refcat(eta, family = draws$family)
  p <- pcategorical(seq_len(draws$data$ncat), eta = eta)
  first_greater(p, target = runif(draws$nsamples, min = 0, max = 1))
}

predict_multinomial <- function(i, draws, ...) {
  eta <- sapply(names(draws$dpars), get_dpar, draws = draws, i = i)
  eta <- insert_refcat(eta, family = draws$family)
  p <- pcategorical(seq_len(draws$data$ncat), eta = eta)
  size <- draws$data$trials[i]
  out <- lapply(seq_rows(p), function(s) t(rmultinom(1, size, p[s, ])))
  do_call(rbind, out)
}

predict_dirichlet <- function(i, draws, ...) {
  mu_dpars <- str_subset(names(draws$dpars), "^mu")
  eta <- sapply(mu_dpars, get_dpar, draws = draws, i = i)
  eta <- insert_refcat(eta, family = draws$family)
  phi <- get_dpar(draws, "phi", i = i)
  cats <- seq_len(draws$data$ncat)
  alpha <- dcategorical(cats, eta = eta) * phi
  rdirichlet(draws$nsamples, alpha = alpha)
}

predict_cumulative <- function(i, draws, ...) {
  predict_ordinal(i = i, draws = draws)
}

predict_sratio <- function(i, draws, ...) {
  predict_ordinal(i = i, draws = draws)
}

predict_cratio <- function(i, draws, ...) {
  predict_ordinal(i = i, draws = draws)
}

predict_acat <- function(i, draws, ...) {
  predict_ordinal(i = i, draws = draws)
}  

predict_ordinal <- function(i, draws, ...) {
  thres <- subset_thres(draws, i)
  nthres <- NCOL(thres)
  p <- pordinal(
    seq_len(nthres + 1), 
    eta = get_dpar(draws, "mu", i = i), 
    disc = get_dpar(draws, "disc", i = i),
    thres = thres,
    family = draws$family$family, 
    link = draws$family$link
  )
  first_greater(p, target = runif(draws$nsamples, min = 0, max = 1))
}

predict_custom <- function(i, draws, ...) {
  predict_fun <- draws$family$predict
  if (!is.function(predict_fun)) {
    predict_fun <- paste0("predict_", draws$family$name)
    predict_fun <- get(predict_fun, draws$family$env)
  }
  predict_fun(i = i, draws = draws, ...)
}

predict_mixture <- function(i, draws, ...) {
  families <- family_names(draws$family)
  theta <- get_theta(draws, i = i)
  smix <- sample_mixture_ids(theta)
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

# ------------ predict helper-functions ----------------------
# random numbers from (possibly truncated) continuous distributions
# @param n number of random values to generate
# @param dist name of a distribution for which the functions
#   p<dist>, q<dist>, and r<dist> are available
# @param ... additional arguments passed to the distribution functions
# @return vector of random values draws from the distribution
rcontinuous <- function(n, dist, ..., lb = NULL, ub = NULL) {
  args <- list(...)
  if (is.null(lb) && is.null(ub)) {
    # sample as usual
    rdist <- paste0("r", dist)
    out <- do_call(rdist, c(list(n), args))
  } else {
    # sample from truncated distribution
    if (is.null(lb)) lb <- -Inf
    if (is.null(ub)) ub <- Inf
    pdist <- paste0("p", dist)
    qdist <- paste0("q", dist)
    plb <- do_call(pdist, c(list(lb), args))
    pub <- do_call(pdist, c(list(ub), args))
    out <- runif(n, min = plb, max = pub)
    out <- do_call(qdist, c(list(out), args))
    # remove infinte values caused by numerical imprecision
    out[out %in% c(-Inf, Inf)] <- NA
  }
  out
}

# random numbers from (possibly truncated) discrete distributions
# currently rejection sampling is used for truncated distributions
# @param n number of random values to generate
# @param dist name of a distribution for which the functions
#   p<dist>, q<dist>, and r<dist> are available
# @param ... additional arguments passed to the distribution functions
# @param lb optional lower truncation bound
# @param ub optional upper truncation bound
# @param ntrys number of trys in rejection sampling for truncated models
# @return a vector of random values draws from the distribution
rdiscrete <- function(n, dist, ..., lb = NULL, ub = NULL, ntrys = 5) {
  args <- list(...)
  rdist <- paste0("r", dist)
  if (is.null(lb) && is.null(ub)) {
    # sample as usual
    out <- do_call(rdist, c(list(n), args))
  } else {
    # sample from truncated distribution via rejection sampling
    if (is.null(lb)) lb <- -Inf
    if (is.null(ub)) ub <- Inf
    out <- matrix(do_call(rdist, c(list(n * ntrys), args)), ncol = ntrys)
    out <- apply(out, 1, extract_valid_sample, lb = lb, ub = ub)
  }
  out
}

# sample from the IDs of the mixture components
sample_mixture_ids <- function(theta) {
  stopifnot(is.matrix(theta))
  mix_comp <- seq_cols(theta)
  ulapply(seq_rows(theta), function(s)
    sample(mix_comp, 1, prob = theta[s, ])
  )
}

# extract the first valid predicted value per Stan sample per observation 
# @param x draws to be check against truncation boundaries
# @param lb vector of lower bounds
# @param ub vector of upper bound
# @return a valid truncated sample or else the closest boundary
extract_valid_sample <- function(x, lb, ub) {
  valid <- match(TRUE, x >= lb & x <= ub)
  if (is.na(valid)) {
    # no valid truncated value found
    # set sample to lb or ub
    # 1e-10 is only to identify the invalid draws later on
    out <- ifelse(max(x) < lb, lb - 1e-10, ub + 1e-10)
  } else {
    out <- x[valid]
  }
  out
}

# check for invalid predictions of truncated discrete models
# @param x matrix of predicted values
# @param lb optional lower truncation bound
# @param ub optional upper truncation bound
# @param thres threshold (in %) of invalid values at which to warn the user
# @return rounded values of 'x'
check_discrete_trunc_bounds <- function(x, lb = NULL, ub = NULL, thres = 0.01) {
  if (is.null(lb) && is.null(ub)) {
    return(x)
  }
  if (is.null(lb)) lb <- -Inf
  if (is.null(ub)) ub <- Inf
  thres <- as_one_numeric(thres)
  # ensure correct comparison with vector bounds
  y <- as.vector(t(x))
  pct_invalid <- mean(y < lb | y > ub, na.rm = TRUE)
  if (pct_invalid >= thres) {
    warning2(
      round(pct_invalid * 100), "% of all predicted values ", 
      "were invalid. Increasing argument 'ntrys' may help."
    )
  }
  round(x)
}
