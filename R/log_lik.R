#' Compute the Pointwise Log-Likelihood
#' 
#' @aliases log_lik logLik.brmsfit
#' 
#' @param object A fitted model object of class \code{brmsfit}. 
#' @inheritParams posterior_predict.brmsfit
#' @param combine Only relevant in multivariate models.
#'   Indicates if the log-likelihoods of the submodels should
#'   be combined per observation (i.e. added together; the default) 
#'   or if the log-likelihoods should be returned separately.
#' @param pointwise A flag indicating whether to compute the full
#'   log-likelihood matrix at once (the default), or just return
#'   the likelihood function along with all data and samples
#'   required to compute the log-likelihood separately for each
#'   observation. The latter option is rarely useful when
#'   calling \code{log_lik} directly, but rather when computing
#'   \code{\link{waic}} or \code{\link{loo}}.
#' @param add_point_draws For internal use only. Ensures compatibility with the
#'   \code{\link{loo_subsample}} method.
#' 
#' @return Usually, an S x N matrix containing the pointwise log-likelihood
#'  samples, where S is the number of samples and N is the number 
#'  of observations in the data. For multivariate models and if 
#'  \code{combine} is \code{FALSE}, an S x N x R array is returned, 
#'  where R is the number of response variables.
#'  If \code{pointwise = TRUE}, the output is a function
#'  with a \code{draws} attribute containing all relevant
#'  data and posterior samples.
#' 
#' @aliases log_lik
#' @method log_lik brmsfit
#' @export
#' @export log_lik
#' @importFrom rstantools log_lik
log_lik.brmsfit <- function(object, newdata = NULL, re_formula = NULL,
                            resp = NULL, nsamples = NULL, subset = NULL, 
                            pointwise = FALSE, combine = TRUE,
                            add_point_draws = FALSE, ...) {
  pointwise <- as_one_logical(pointwise)
  combine <- as_one_logical(combine)
  add_point_draws <- as_one_logical(add_point_draws)
  contains_samples(object)
  object <- restructure(object)
  draws <- extract_draws(
    object, newdata = newdata, re_formula = re_formula, resp = resp, 
    subset = subset, nsamples = nsamples, check_response = TRUE, ...
  )
  if (add_point_draws) {
    # required for the loo_subsample method
    # Computing a point estimate based on the full draws object is too
    # difficult due to its highly nested structure. As an alternative, a second
    # draws object is created from the point estimates of the samples directly.
    attr(draws, "point_draws") <- extract_draws(
      object, newdata = newdata, re_formula = re_formula, resp = resp, 
      subset = subset, nsamples = nsamples, check_response = TRUE, 
      point = "median", ...
    )
  }
  if (pointwise) {
    stopifnot(combine)
    log_lik <- log_lik_pointwise
    attr(log_lik, "draws") <- draws
    attr(log_lik, "data") <- data.frame(i = seq_len(choose_N(draws)))
  } else {
    log_lik <- log_lik(draws, combine = combine)
    if (anyNA(log_lik)) {
      warning2(
        "NAs were found in the log-likelihood. Possibly this is because ",
        "some of your responses contain NAs. If you use 'mi' terms, try ", 
        "setting 'resp' to those response variables without missing values. ",
        "Alternatively, use 'newdata' to predict only complete cases."
      )
    }
  }
  log_lik
}

#' @export
logLik.brmsfit <- function(object, newdata = NULL, re_formula = NULL,
                           resp = NULL, nsamples = NULL, subset = NULL, 
                           pointwise = FALSE, combine = TRUE, ...) {
  cl <- match.call()
  cl[[1]] <- quote(log_lik)
  eval(cl, parent.frame())
}

#' @export
log_lik.mvbrmsdraws <- function(object, combine = TRUE, ...) {
  if (length(object$mvpars$rescor)) {
    object$mvpars$Mu <- get_Mu(object)
    object$mvpars$Sigma <- get_Sigma(object)
    out <- log_lik.brmsdraws(object, ...)
  } else {
    out <- lapply(object$resps, log_lik, ...)
    if (combine) {
      out <- Reduce("+", out)
    } else {
      along <- ifelse(length(out) > 1L, 3, 2)
      out <- do_call(abind, c(out, along = along))
    }
  }
  out
}

#' @export
log_lik.brmsdraws <- function(object, ...) {
  log_lik_fun <- paste0("log_lik_", object$family$fun)
  log_lik_fun <- get(log_lik_fun, asNamespace("brms"))
  for (nlp in names(object$nlpars)) {
    object$nlpars[[nlp]] <- get_nlpar(object, nlpar = nlp)
  }
  for (dp in names(object$dpars)) {
    object$dpars[[dp]] <- get_dpar(object, dpar = dp)
  }
  N <- choose_N(object)
  out <- cblapply(seq_len(N), log_lik_fun, draws = object)
  colnames(out) <- NULL
  old_order <- object$old_order
  sort <- isTRUE(ncol(out) != length(old_order))
  reorder_obs(out, old_order, sort = sort)
}

# evaluate log_lik in a pointwise manner
# cannot be an S3 method since 'i' must be the first argument
log_lik_pointwise <- function(data_i, draws, ...) {
  i <- data_i$i
  if (is.mvbrmsdraws(draws) && !length(draws$mvpars$rescor)) {
    out <- lapply(draws$resps, log_lik_pointwise, i = i)
    out <- Reduce("+", out)
  } else {
    log_lik_fun <- paste0("log_lik_", draws$family$fun)
    log_lik_fun <- get(log_lik_fun, asNamespace("brms"))
    out <- log_lik_fun(i, draws)
  }
  out
}

# All log_lik_<family> functions have the same arguments structure
# @param i the column of draws to use the the ith obervation 
#   in the initial data.frame 
# @param draws A named list returned by extract_draws containing 
#   all required data and samples
# @param data (ignored) included for compatibility with loo::loo.function      
# @return a vector of length draws$nsamples containing the pointwise 
#   log-likelihood fo the ith observation 
log_lik_gaussian <- function(i, draws, data = data.frame()) {
  args <- list(
    mean = get_dpar(draws, "mu", i = i), 
    sd = get_dpar(draws, "sigma", i = i)
  )
  # log_lik_censor computes the conventional log_lik in case of no censoring 
  out <- log_lik_censor(dist = "norm", args = args, i = i, draws = draws)
  out <- log_lik_truncate(
    out, cdf = pnorm, args = args, i = i, draws = draws
  )
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_student <- function(i, draws, data = data.frame()) {
  args <- list(
    df = get_dpar(draws, "nu", i = i), 
    mu = get_dpar(draws, "mu", i = i), 
    sigma = get_dpar(draws, "sigma", i = i)
  )
  out <- log_lik_censor(
    dist = "student_t", args = args, i = i, draws = draws
  )
  out <- log_lik_truncate(
    out, cdf = pstudent_t, args = args, i = i, draws = draws
  )
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_lognormal <- function(i, draws, data = data.frame()) {
  sigma <- get_dpar(draws, "sigma", i = i)
  args <- list(meanlog = get_dpar(draws, "mu", i), sdlog = sigma)
  out <- log_lik_censor(dist = "lnorm", args = args, i = i, draws = draws)
  out <- log_lik_truncate(
    out, cdf = plnorm, args = args, i = i, draws = draws
  )
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_shifted_lognormal <- function(i, draws, data = data.frame()) {
  sigma <- get_dpar(draws, "sigma", i = i)
  ndt <- get_dpar(draws, "ndt", i = i)
  args <- list(meanlog = get_dpar(draws, "mu", i), sdlog = sigma, shift = ndt)
  out <- log_lik_censor("shifted_lnorm", args, i = i, draws = draws)
  out <- log_lik_truncate(out, pshifted_lnorm, args, i = i, draws = draws)
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_skew_normal <- function(i, draws, data = data.frame()) {
  sigma <- get_dpar(draws, "sigma", i = i)
  alpha <- get_dpar(draws, "alpha", i = i)
  mu <- get_dpar(draws, "mu", i)
  args <- nlist(mu, sigma, alpha)
  out <- log_lik_censor(
    dist = "skew_normal", args = args, i = i, draws = draws
  )
  out <- log_lik_truncate(
    out, cdf = pskew_normal, args = args, i = i, draws = draws
  )
  log_lik_weight(out, i = i, draws = draws)
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
  log_lik_weight(out, i = i, draws = draws)
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
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_gaussian_time <- function(i, draws, data = data.frame()) {
  obs <- with(draws$ac, begin_tg[i]:end_tg[i])
  Y <- as.numeric(draws$data$Y[obs])
  mu <- as.matrix(get_dpar(draws, "mu", i = obs))
  Sigma <- get_cov_matrix_ac(draws, obs)
  .log_lik <- function(s) {
    C <- as.matrix(Sigma[s, , ])
    Cinv <- solve(C)
    e <- Y - mu[s, ]
    g <- solve(C, e)
    cbar <- diag(Cinv)
    yloo <- Y - g / cbar
    sdloo <- sqrt(1 / cbar)
    ll <- dnorm(Y, yloo, sdloo, log = TRUE)
    return(as.numeric(ll))
  }
  rblapply(seq_len(draws$nsamples), .log_lik)
}

log_lik_student_time <- function(i, draws, data = data.frame()) {
  obs <- with(draws$ac, begin_tg[i]:end_tg[i])
  Y <- as.numeric(draws$data$Y[obs])
  nu <- as.matrix(get_dpar(draws, "nu", i = obs))
  mu <- as.matrix(get_dpar(draws, "mu", i = obs))
  Sigma <- get_cov_matrix_ac(draws, obs)
  .log_lik <- function(s) {
    df <- nu[s, ]
    C <- as.matrix(Sigma[s, , ])
    Cinv <- solve(C)
    e <- Y - mu[s, ]
    g <- solve(C, e)
    cbar <- diag(Cinv)
    yloo <- Y - g / cbar
    sdloo <- sqrt(1 / cbar * student_t_cov_factor(df, Cinv, e))
    dfloo <- df + nrow(Cinv) - 1
    ll <- dstudent_t(Y, dfloo, yloo, sdloo, log = TRUE)
    return(as.numeric(ll))
  }
  rblapply(seq_len(draws$nsamples), .log_lik)
}

log_lik_gaussian_lagsar <- function(i, draws, data = data.frame()) {
  mu <- get_dpar(draws, "mu")
  sigma <- get_dpar(draws, "sigma")
  Y <- as.numeric(draws$data$Y)
  I <- diag(draws$nobs)
  stopifnot(i == 1)
  # see http://mc-stan.org/loo/articles/loo2-non-factorizable.html
  .log_lik <- function(s) {
    IB <- I - with(draws$ac, lagsar[s, ] * Msar)
    Cinv <- t(IB) %*% IB / sigma[s]^2
    e <- Y - solve(IB, mu[s, ])
    g <- Cinv %*% e
    cbar <- diag(Cinv)
    yloo <- Y - g / cbar
    sdloo <- sqrt(1 / cbar)
    ll <- dnorm(Y, yloo, sdloo, log = TRUE)
    return(as.numeric(ll))
  }
  rblapply(seq_len(draws$nsamples), .log_lik)
}

log_lik_student_lagsar <- function(i, draws, data = data.frame()) {
  nu <- get_dpar(draws, "nu")
  mu <- get_dpar(draws, "mu")
  sigma <- get_dpar(draws, "sigma")
  Y <- as.numeric(draws$data$Y)
  I <- diag(draws$nobs)
  stopifnot(i == 1)
  # see http://mc-stan.org/loo/articles/loo2-non-factorizable.html
  .log_lik <- function(s) {
    df <- nu[s]
    IB <- I - with(draws$ac, lagsar[s, ] * Msar)
    Cinv <- t(IB) %*% IB / sigma[s]^2
    e <- Y - solve(IB, mu[s, ])
    g <- Cinv %*% e
    cbar <- diag(Cinv)
    yloo <- Y - g / cbar
    sdloo <- sqrt(1 / cbar * student_t_cov_factor(df, Cinv, e))
    dfloo <- df + nrow(Cinv) - 1
    ll <- dstudent_t(Y, dfloo, yloo, sdloo, log = TRUE)
    return(as.numeric(ll))
  }
  rblapply(seq_len(draws$nsamples), .log_lik)
}

log_lik_gaussian_errorsar <- function(i, draws, data = data.frame()) {
  stopifnot(i == 1)
  mu <- get_dpar(draws, "mu")
  sigma <- get_dpar(draws, "sigma")
  Y <- as.numeric(draws$data$Y)
  I <- diag(draws$nobs)
  .log_lik <- function(s) {
    IB <- I - with(draws$ac, errorsar[s, ] * Msar)
    Cinv <- t(IB) %*% IB / sigma[s]^2
    e <- Y - mu[s, ]
    g <- Cinv %*% e
    cbar <- diag(Cinv)
    yloo <- Y - g / cbar
    sdloo <- sqrt(1 / cbar)
    ll <- dnorm(Y, yloo, sdloo, log = TRUE)
    return(as.numeric(ll))
  }
  rblapply(seq_len(draws$nsamples), .log_lik)
}

log_lik_student_errorsar <- function(i, draws, data = data.frame()) {
  stopifnot(i == 1)
  nu <- get_dpar(draws, "nu")
  mu <- get_dpar(draws, "mu")
  sigma <- get_dpar(draws, "sigma")
  Y <- as.numeric(draws$data$Y)
  I <- diag(draws$nobs)
  .log_lik <- function(s) {
    df <- nu[s]
    IB <- I - with(draws$ac, errorsar[s, ] * Msar)
    Cinv <- t(IB) %*% IB / sigma[s]^2
    e <- Y - mu[s, ]
    g <- Cinv %*% e
    cbar <- diag(Cinv)
    yloo <- Y - g / cbar
    sdloo <- sqrt(1 / cbar * student_t_cov_factor(df, Cinv, e))
    dfloo <- df + nrow(Cinv) - 1
    ll <- dstudent_t(Y, dfloo, yloo, sdloo, log = TRUE)
    return(as.numeric(ll))
  }
  rblapply(seq_len(draws$nsamples), .log_lik)
}

log_lik_gaussian_fcor <- function(i, draws, data = data.frame()) {
  stopifnot(i == 1)
  Y <- as.numeric(draws$data$Y)
  mu <- get_dpar(draws, "mu")
  Sigma <- get_cov_matrix_ac(draws)
  .log_lik <- function(s) {
    C <- as.matrix(Sigma[s, , ])
    Cinv <- solve(C)
    e <- Y - mu[s, ]
    g <- solve(C, e)
    cbar <- diag(Cinv)
    yloo <- Y - g / cbar
    sdloo <- sqrt(1 / cbar)
    ll <- dnorm(Y, yloo, sdloo, log = TRUE)
    return(as.numeric(ll))
  }
  rblapply(seq_len(draws$nsamples), .log_lik)
}

log_lik_student_fcor <- function(i, draws, data = data.frame()) {
  stopifnot(i == 1)
  Y <- as.numeric(draws$data$Y)
  nu <- get_dpar(draws, "nu")
  mu <- get_dpar(draws, "mu")
  Sigma <- get_cov_matrix_ac(draws)
  .log_lik <- function(s) {
    df <- nu[s]
    C <- as.matrix(Sigma[s, , ])
    Cinv <- solve(C)
    e <- Y - mu[s, ]
    g <- solve(C, e)
    cbar <- diag(Cinv)
    yloo <- Y - g / cbar
    sdloo <- sqrt(1 / cbar * student_t_cov_factor(df, Cinv, e))
    dfloo <- df + nrow(Cinv) - 1
    ll <- dstudent_t(Y, dfloo, yloo, sdloo, log = TRUE)
    return(as.numeric(ll))
  }
  rblapply(seq_len(draws$nsamples), .log_lik)
}

log_lik_binomial <- function(i, draws, data = data.frame()) {
  trials <- draws$data$trials[i]
  args <- list(size = trials, prob = get_dpar(draws, "mu", i))
  out <- log_lik_censor(
    dist = "binom", args = args, i = i, draws = draws
  )
  out <- log_lik_truncate(
    out, cdf = pbinom, args = args, i = i, draws = draws
  )
  log_lik_weight(out, i = i, draws = draws)
}  

log_lik_bernoulli <- function(i, draws, data = data.frame()) {
  args <- list(size = 1, prob = get_dpar(draws, "mu", i))
  out <- log_lik_censor(
    dist = "binom", args = args, i = i, draws = draws
  )
  # no truncation allowed
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_poisson <- function(i, draws, data = data.frame()) {
  args <- list(lambda = get_dpar(draws, "mu", i))
  out <- log_lik_censor(
    dist = "pois", args = args, i = i, draws = draws
  )
  out <- log_lik_truncate(
    out, cdf = ppois, args = args, i = i, draws = draws
  )
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_negbinomial <- function(i, draws, data = data.frame()) {
  shape <- get_dpar(draws, "shape", i = i)
  args <- list(mu = get_dpar(draws, "mu", i), size = shape)
  out <- log_lik_censor(
    dist = "nbinom", args = args, i = i, draws = draws
  )
  out <- log_lik_truncate(
    out, cdf = pnbinom, args = args, i = i, draws = draws
  )
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_geometric <- function(i, draws, data = data.frame()) {
  args <- list(mu = get_dpar(draws, "mu", i), size = 1)
  out <- log_lik_censor(
    dist = "nbinom", args = args, i = i, draws = draws
  )
  out <- log_lik_truncate(
    out, cdf = pnbinom, args = args, i = i, draws = draws
  )
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_discrete_weibull <- function(i, draws, data = data.frame()) {
  args <- list(
    mu = get_dpar(draws, "mu", i), 
    shape = get_dpar(draws, "shape", i = i)
  )
  out <- log_lik_censor(
    dist = "discrete_weibull", args = args, i = i, draws = draws
  )
  out <- log_lik_truncate(
    out, cdf = pdiscrete_weibull, args = args, i = i, draws = draws
  )
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_com_poisson <- function(i, draws, data = data.frame()) {
  args <- list(
    mu = get_dpar(draws, "mu", i), 
    shape = get_dpar(draws, "shape", i = i)
  )
  # no censoring or truncation allowed yet
  out <- do_call(dcom_poisson, c(draws$data$Y[i], args, log = TRUE))
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_exponential <- function(i, draws, data = data.frame()) {
  args <- list(rate = 1 / get_dpar(draws, "mu", i))
  out <- log_lik_censor(dist = "exp", args = args, i = i, draws = draws)
  out <- log_lik_truncate(
    out, cdf = pexp, args = args, i = i, draws = draws
  )
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_gamma <- function(i, draws, data = data.frame()) {
  shape <- get_dpar(draws, "shape", i = i)
  scale <- get_dpar(draws, "mu", i) / shape
  args <- nlist(shape, scale)
  out <- log_lik_censor(dist = "gamma", args = args, i = i, draws = draws)
  out <- log_lik_truncate(
    out, cdf = pgamma, args = args, i = i, draws = draws
  )
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_weibull <- function(i, draws, data = data.frame()) {
  shape <- get_dpar(draws, "shape", i = i)
  scale <- get_dpar(draws, "mu", i = i) / gamma(1 + 1 / shape)
  args <- list(shape = shape, scale = scale)
  out <- log_lik_censor(
    dist = "weibull", args = args, i = i, draws = draws
  )
  out <- log_lik_truncate(
    out, cdf = pweibull, args = args, i = i, draws = draws
  )
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_frechet <- function(i, draws, data = data.frame()) {
  nu <- get_dpar(draws, "nu", i = i)
  scale <- get_dpar(draws, "mu", i = i) / gamma(1 - 1 / nu)
  args <- list(scale = scale, shape = nu)
  out <- log_lik_censor(
    dist = "frechet", args = args, i = i, draws = draws
  )
  out <- log_lik_truncate(
    out, cdf = pfrechet, args = args, i = i, draws = draws
  )
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_gen_extreme_value <- function(i, draws, data = data.frame()) {
  sigma <- get_dpar(draws, "sigma", i = i)
  xi <- get_dpar(draws, "xi", i = i)
  mu <- get_dpar(draws, "mu", i)
  args <- nlist(mu, sigma, xi)
  out <- log_lik_censor(dist = "gen_extreme_value", args = args, 
                       i = i, draws = draws)
  out <- log_lik_truncate(out, cdf = pgen_extreme_value, 
                         args = args, i = i, draws = draws)
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_inverse.gaussian <- function(i, draws, data = data.frame()) {
  args <- list(mu = get_dpar(draws, "mu", i), 
               shape = get_dpar(draws, "shape", i = i))
  out <- log_lik_censor(dist = "inv_gaussian", args = args, 
                       i = i, draws = draws)
  out <- log_lik_truncate(out, cdf = pinv_gaussian, args = args,
                         i = i, draws = draws)
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_exgaussian <- function(i, draws, data = data.frame()) {
  args <- list(mu = get_dpar(draws, "mu", i), 
               sigma = get_dpar(draws, "sigma", i = i),
               beta = get_dpar(draws, "beta", i = i))
  out <- log_lik_censor(dist = "exgaussian", args = args, 
                       i = i, draws = draws)
  out <- log_lik_truncate(out, cdf = pexgaussian, args = args,
                         i = i, draws = draws)
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_wiener <- function(i, draws, data = data.frame()) {
  args <- list(
    delta = get_dpar(draws, "mu", i), 
    alpha = get_dpar(draws, "bs", i = i),
    tau = get_dpar(draws, "ndt", i = i),
    beta = get_dpar(draws, "bias", i = i),
    resp = draws$data[["dec"]][i]
  )
  out <- do_call(dwiener, c(draws$data$Y[i], args, log = TRUE))
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_beta <- function(i, draws, data = data.frame()) {
  mu <- get_dpar(draws, "mu", i)
  phi <- get_dpar(draws, "phi", i)
  args <- list(shape1 = mu * phi, shape2 = (1 - mu) * phi)
  out <- log_lik_censor(dist = "beta", args = args, i = i, draws = draws)
  out <- log_lik_truncate(
    out, cdf = pbeta, args = args, i = i, draws = draws
  )
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_von_mises <- function(i, draws, data = data.frame()) {
  args <- list(
    mu = get_dpar(draws, "mu", i), 
    kappa = get_dpar(draws, "kappa", i = i)
  )
  out <- log_lik_censor(
    dist = "von_mises", args = args, i = i, draws = draws
  )
  out <- log_lik_truncate(
    out, cdf = pvon_mises, args = args, i = i, draws = draws
  )
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_asym_laplace <- function(i, draws, ...) {
  args <- list(
    mu = get_dpar(draws, "mu", i), 
    sigma = get_dpar(draws, "sigma", i),
    quantile = get_dpar(draws, "quantile", i)
  )
  out <- log_lik_censor(dist = "asym_laplace", args, i, draws)
  out <- log_lik_truncate(out, pasym_laplace, args, i, draws)
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_zero_inflated_asym_laplace <- function(i, draws, ...) {
  args <- list(
    mu = get_dpar(draws, "mu", i), 
    sigma = get_dpar(draws, "sigma", i),
    quantile = get_dpar(draws, "quantile", i),
    zi = get_dpar(draws, "zi", i)
  )
  out <- log_lik_censor(dist = "zero_inflated_asym_laplace", args, i, draws)
  out <- log_lik_truncate(out, pzero_inflated_asym_laplace, args, i, draws)
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_cox <- function(i, draws, ...) {
  args <- list(
    mu = get_dpar(draws, "mu", i),
    bhaz = draws$bhaz$bhaz[, i], 
    cbhaz = draws$bhaz$cbhaz[, i]
  )
  out <- log_lik_censor(dist = "cox", args = args, i = i, draws = draws)
  out <- log_lik_truncate(out, cdf = pcox, args = args, i = i, draws = draws)
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_hurdle_poisson <- function(i, draws, data = data.frame()) {
  hu <- get_dpar(draws, "hu", i)
  lambda <- get_dpar(draws, "mu", i)
  args <- nlist(lambda, hu)
  out <- log_lik_censor("hurdle_poisson", args, i, draws)
  out <- log_lik_truncate(out, phurdle_poisson, args, i, draws)
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_hurdle_negbinomial <- function(i, draws, data = data.frame()) {
  hu <- get_dpar(draws, "hu", i)
  mu <- get_dpar(draws, "mu", i)
  shape <- get_dpar(draws, "shape", i = i)
  args <- nlist(mu, shape, hu)
  out <- log_lik_censor("hurdle_negbinomial", args, i, draws)
  out <- log_lik_truncate(out, phurdle_negbinomial, args, i, draws)
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_hurdle_gamma <- function(i, draws, data = data.frame()) {
  hu <- get_dpar(draws, "hu", i)
  shape <- get_dpar(draws, "shape", i = i)
  scale <- get_dpar(draws, "mu", i) / shape
  args <- nlist(shape, scale, hu)
  out <- log_lik_censor("hurdle_gamma", args, i, draws)
  out <- log_lik_truncate(out, phurdle_gamma, args, i, draws)
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_hurdle_lognormal <- function(i, draws, data = data.frame()) {
  hu <- get_dpar(draws, "hu", i)
  mu <- get_dpar(draws, "mu", i)
  sigma <- get_dpar(draws, "sigma", i = i)
  args <- nlist(mu, sigma, hu)
  out <- log_lik_censor("hurdle_lognormal", args, i, draws)
  out <- log_lik_truncate(out, phurdle_lognormal, args, i, draws)
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_zero_inflated_poisson <- function(i, draws, data = data.frame()) {
  zi <- get_dpar(draws, "zi", i)
  lambda <- get_dpar(draws, "mu", i)
  args <- nlist(lambda, zi)
  out <- log_lik_censor("zero_inflated_poisson", args, i, draws)
  out <- log_lik_truncate(out, pzero_inflated_poisson, args, i, draws)
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_zero_inflated_negbinomial <- function(i, draws, data = data.frame()) {
  zi <- get_dpar(draws, "zi", i)
  mu <- get_dpar(draws, "mu", i)
  shape <- get_dpar(draws, "shape", i = i)
  args <- nlist(mu, shape, zi)
  out <- log_lik_censor("zero_inflated_negbinomial", args, i, draws)
  out <- log_lik_truncate(out, pzero_inflated_negbinomial, args, i, draws)
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_zero_inflated_binomial <- function(i, draws, data = data.frame()) {
  trials <- draws$data$trials[i] 
  mu <- get_dpar(draws, "mu", i) 
  zi <- get_dpar(draws, "zi", i)
  args <- list(size = trials, prob = mu, zi)
  out <- log_lik_censor("zero_inflated_binomial", args, i, draws)
  out <- log_lik_truncate(out, pzero_inflated_binomial, args, i, draws)
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_zero_inflated_beta <- function(i, draws, data = data.frame()) {
  zi <- get_dpar(draws, "zi", i)
  mu <- get_dpar(draws, "mu", i)
  phi <- get_dpar(draws, "phi", i)
  args <- nlist(shape1 = mu * phi, shape2 = (1 - mu) * phi, zi)
  out <- log_lik_censor("zero_inflated_beta", args, i, draws)
  out <- log_lik_truncate(out, pzero_inflated_beta, args, i, draws)
  log_lik_weight(out, i = i, draws = draws)
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
      do_call(dbeta, c(draws$data$Y[i], args, log = TRUE))
  }
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_categorical <- function(i, draws, data = data.frame()) {
  stopifnot(draws$family$link == "logit")
  eta <- sapply(names(draws$dpars), get_dpar, draws = draws, i = i)
  eta <- insert_refcat(eta, family = draws$family)
  out <- dcategorical(draws$data$Y[i], eta = eta, log = TRUE)
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_multinomial <- function(i, draws, data = data.frame()) {
  stopifnot(draws$family$link == "logit")
  eta <- sapply(names(draws$dpars), get_dpar, draws = draws, i = i)
  eta <- insert_refcat(eta, family = draws$family)
  out <- dmultinomial(draws$data$Y[i, ], eta = eta, log = TRUE)
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_dirichlet <- function(i, draws, data = data.frame()) {
  stopifnot(draws$family$link == "logit")
  mu_dpars <- str_subset(names(draws$dpars), "^mu")
  eta <- sapply(mu_dpars, get_dpar, draws = draws, i = i)
  eta <- insert_refcat(eta, family = draws$family)
  phi <- get_dpar(draws, "phi", i = i)
  cats <- seq_len(draws$data$ncat)
  alpha <- dcategorical(cats, eta = eta) * phi
  out <- ddirichlet(draws$data$Y[i, ], alpha = alpha, log = TRUE)
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_cumulative <- function(i, draws, data = data.frame()) {
  disc <- get_dpar(draws, "disc", i = i)
  mu <- get_dpar(draws, "mu", i = i)
  thres <- subset_thres(draws, i)
  nthres <- NCOL(thres)
  eta <- disc * (thres - mu)
  y <- draws$data$Y[i]
  if (y == 1) { 
    out <- log(ilink(eta[, 1], draws$family$link))
  } else if (y == nthres + 1) {
    out <- log(1 - ilink(eta[, y - 1], draws$family$link)) 
  } else {
    out <- log(
      ilink(eta[, y], draws$family$link) - 
        ilink(eta[, y - 1], draws$family$link)
    )
  }
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_sratio <- function(i, draws, data = data.frame()) {
  disc <- get_dpar(draws, "disc", i = i)
  mu <- get_dpar(draws, "mu", i = i)
  thres <- subset_thres(draws, i)
  nthres <- NCOL(thres)
  eta <- disc * (thres - mu)
  y <- draws$data$Y[i]
  q <- sapply(seq_len(min(y, nthres)), 
    function(k) 1 - ilink(eta[, k], draws$family$link)
  )
  if (y == 1) {
    out <- log(1 - q[, 1]) 
  } else if (y == 2) {
    out <- log(1 - q[, 2]) + log(q[, 1])
  } else if (y == nthres + 1) {
    out <- rowSums(log(q))
  } else {
    out <- log(1 - q[, y]) + rowSums(log(q[, 1:(y - 1)]))
  }
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_cratio <- function(i, draws, data = data.frame()) {
  disc <- get_dpar(draws, "disc", i = i)
  mu <- get_dpar(draws, "mu", i = i)
  thres <- subset_thres(draws, i)
  nthres <- NCOL(thres)
  eta <- disc * (mu - thres)
  y <- draws$data$Y[i]
  q <- sapply(seq_len(min(y, nthres)), 
    function(k) ilink(eta[, k], draws$family$link)
  )
  if (y == 1) {
    out <- log(1 - q[, 1])
  }  else if (y == 2) {
    out <- log(1 - q[, 2]) + log(q[, 1])
  } else if (y == nthres + 1) {
    out <- rowSums(log(q))
  } else {
    out <- log(1 - q[, y]) + rowSums(log(q[, 1:(y - 1)]))
  }
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_acat <- function(i, draws, data = data.frame()) {
  disc <- get_dpar(draws, "disc", i = i)
  mu <- get_dpar(draws, "mu", i = i)
  thres <- subset_thres(draws, i)
  nthres <- NCOL(thres)
  eta <- disc * (mu - thres)
  y <- draws$data$Y[i]
  if (draws$family$link == "logit") { # more efficient calculation 
    q <- sapply(1:nthres, function(k) eta[, k])
    p <- cbind(rep(0, nrow(eta)), q[, 1], 
               matrix(0, nrow = nrow(eta), ncol = nthres - 1))
    if (nthres > 1L) {
      p[, 3:(nthres + 1)] <- 
        sapply(3:(nthres + 1), function(k) rowSums(q[, 1:(k - 1)]))
    }
    out <- p[, y] - log(rowSums(exp(p)))
  } else {
    q <- sapply(1:nthres, function(k) 
      ilink(eta[, k], draws$family$link))
    p <- cbind(apply(1 - q[, 1:nthres], 1, prod), 
               matrix(0, nrow = nrow(eta), ncol = nthres))
    if (nthres > 1L) {
      p[, 2:nthres] <- sapply(2:nthres, function(k) 
        apply(as.matrix(q[, 1:(k - 1)]), 1, prod) * 
          apply(as.matrix(1 - q[, k:nthres]), 1, prod))
    }
    p[, nthres + 1] <- apply(q[, 1:nthres], 1, prod)
    out <- log(p[, y]) - log(apply(p, 1, sum))
  }
  log_lik_weight(out, i = i, draws = draws)
}

log_lik_custom <- function(i, draws, data = data.frame()) {
  log_lik_fun <- draws$family$log_lik
  if (!is.function(log_lik_fun)) {
    log_lik_fun <- paste0("log_lik_", draws$family$name)
    log_lik_fun <- get(log_lik_fun, draws$family$env)
  }
  log_lik_fun(i = i, draws = draws)
}

log_lik_mixture <- function(i, draws, data = data.frame()) {
  families <- family_names(draws$family)
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
  log_lik_weight(out, i = i, draws = draws)
}

# ----------- log_lik helper-functions -----------
# compute (possibly censored) log_lik values
# @param dist name of a distribution for which the functions
#   d<dist> (pdf) and p<dist> (cdf) are available
# @param args additional arguments passed to pdf and cdf
# @param draws a brmsdraws object
# @return vector of log_lik values
log_lik_censor <- function(dist, args, i, draws) {
  pdf <- get(paste0("d", dist), mode = "function")
  cdf <- get(paste0("p", dist), mode = "function")
  y <- draws$data$Y[i]
  cens <- draws$data$cens[i]
  if (is.null(cens) || cens == 0) {
    x <- do_call(pdf, c(y, args, log = TRUE))
  } else if (cens == 1) {
    x <- do_call(cdf, c(y, args, lower.tail = FALSE, log.p = TRUE))
  } else if (cens == -1) {
    x <- do_call(cdf, c(y, args, log.p = TRUE))
  } else if (cens == 2) {
    rcens <- draws$data$rcens[i]
    x <- log(do_call(cdf, c(rcens, args)) - do_call(cdf, c(y, args)))
  }
  x
}

# adjust log_lik in truncated models
# @param x vector of log_lik values
# @param cdf a cumulative distribution function 
# @param args arguments passed to cdf
# @param i observation number
# @param draws a brmsdraws object
# @return vector of log_lik values
log_lik_truncate <- function(x, cdf, args, i, draws) {
  lb <- draws$data$lb[i]
  ub <- draws$data$ub[i]
  if (!(is.null(lb) && is.null(ub))) {
    if (is.null(lb)) lb <- -Inf
    if (is.null(ub)) ub <- Inf
    x - log(do_call(cdf, c(ub, args)) - do_call(cdf, c(lb, args)))
  }
  x
}

# weight log_lik values according to defined weights
# @param x vector of log_lik values
# @param i observation number
# @param draws a brmsdraws object
# @return vector of log_lik values
log_lik_weight <- function(x, i, draws) {
  weight <- draws$data$weights[i]
  if (!is.null(weight)) {
    x <- x * weight
  }
  x
}

# after some discussion with Aki Vehtari and Daniel Simpson,
# I disallowed computation of log-likelihood values for some models
# until pointwise solutions are implemented
stop_no_pw <- function() {
  stop2("Cannot yet compute pointwise log-likelihood for this model ",
        "because the observations are not conditionally independent.")
}

# multiplicate factor for conditional student-t models
# see http://proceedings.mlr.press/v33/shah14.pdf
# note that brms parameterizes C instead of Cov(y) = df / (df - 2) * C
# @param df degrees of freedom parameter
# @param Cinv inverse of the full matrix
# @param e vector of error terms, that is, y - mu
student_t_cov_factor <- function(df, Cinv, e) {
  beta1 <- ulapply(seq_rows(Cinv), student_t_beta1_i, Cinv, e)
  (df + beta1) / (df + nrow(Cinv) - 1)
}

# beta1 in equation (6) of http://proceedings.mlr.press/v33/shah14.pdf
# @param i observation index to exclude in the submatrix
# @param Cinv inverse of the full matrix
# @param e vector of error terms, that is, y - mu
# @param vector of length one
student_t_beta1_i <- function(i, Cinv, e) {
  sub_Cinv_i <- sub_inverse_symmetric(Cinv, i)
  t(e[-i]) %*% sub_Cinv_i %*% e[-i]
}

# efficient submatrix inverse for a symmetric matrix
# see http://www.scielo.org.mx/pdf/cys/v20n2/1405-5546-cys-20-02-00251.pdf
# @param Cinv inverse of the full matrix
# @param i observation index to exclude in the submatrix
# @return inverse of the submatrix after removing observation i
sub_inverse_symmetric <- function(Cinv, i) {
  csub <- Cinv[i, -i]
  D <- outer(csub, csub)
  Cinv[-i, -i] - D / Cinv[i, i]
}
