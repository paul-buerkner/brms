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
#'   the likelihood function along with all data and draws
#'   required to compute the log-likelihood separately for each
#'   observation. The latter option is rarely useful when
#'   calling \code{log_lik} directly, but rather when computing
#'   \code{\link{waic}} or \code{\link{loo}}.
#' @param add_point_estimate For internal use only. Ensures compatibility
#'   with the \code{\link{loo_subsample}} method.
#'
#' @return Usually, an S x N matrix containing the pointwise log-likelihood
#'  draws, where S is the number of draws and N is the number
#'  of observations in the data. For multivariate models and if
#'  \code{combine} is \code{FALSE}, an S x N x R array is returned,
#'  where R is the number of response variables.
#'  If \code{pointwise = TRUE}, the output is a function
#'  with a \code{draws} attribute containing all relevant
#'  data and posterior draws.
#'
#' @template details-newdata-na
#' @template details-allow_new_levels
#'
#' @aliases log_lik
#' @method log_lik brmsfit
#' @export
#' @export log_lik
#' @importFrom rstantools log_lik
log_lik.brmsfit <- function(object, newdata = NULL, re_formula = NULL,
                            resp = NULL, ndraws = NULL, draw_ids = NULL,
                            pointwise = FALSE, combine = TRUE,
                            add_point_estimate = FALSE,
                            cores = NULL, ...) {
  pointwise <- as_one_logical(pointwise)
  combine <- as_one_logical(combine)
  add_point_estimate <- as_one_logical(add_point_estimate)
  contains_draws(object)
  object <- restructure(object)
  prep <- prepare_predictions(
    object, newdata = newdata, re_formula = re_formula, resp = resp,
    ndraws = ndraws, draw_ids = draw_ids, check_response = TRUE, ...
  )
  if (add_point_estimate) {
    # required for the loo_subsample method
    # Computing a point estimate based on the full prep object is too
    # difficult due to its highly nested structure. As an alternative, a second
    # prep object is created from the point estimates of the draws directly.
    attr(prep, "point_estimate") <- prepare_predictions(
      object, newdata = newdata, re_formula = re_formula, resp = resp,
      ndraws = ndraws, draw_ids = draw_ids, check_response = TRUE,
      point_estimate = "median", ...
    )
  }
  if (pointwise) {
    stopifnot(combine)
    log_lik <- log_lik_pointwise
    # names need to be 'data' and 'draws' as per ?loo::loo.function
    attr(log_lik, "data") <- data.frame(i = seq_len(choose_N(prep)))
    attr(log_lik, "draws") <- prep
  } else {
    log_lik <- log_lik(prep, combine = combine, cores = cores)
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
                           resp = NULL, ndraws = NULL, draw_ids = NULL,
                           pointwise = FALSE, combine = TRUE,
                           cores = NULL, ...) {
  cl <- match.call()
  cl[[1]] <- quote(log_lik)
  eval(cl, parent.frame())
}

#' @export
log_lik.mvbrmsprep <- function(object, combine = TRUE, ...) {
  if (length(object$mvpars$rescor)) {
    object$mvpars$Mu <- get_Mu(object)
    object$mvpars$Sigma <- get_Sigma(object)
    out <- log_lik.brmsprep(object, ...)
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
log_lik.brmsprep <- function(object, cores = NULL, ...) {
  cores <- validate_cores_post_processing(cores)
  log_lik_fun <- paste0("log_lik_", object$family$fun)
  log_lik_fun <- get(log_lik_fun, asNamespace("brms"))
  if (is.customfamily(object$family)) {
    # ensure that the method can be found during parallel execution
    object$family$log_lik <- custom_family_method(object$family, "log_lik")
  }
  for (nlp in names(object$nlpars)) {
    object$nlpars[[nlp]] <- get_nlpar(object, nlpar = nlp)
  }
  for (dp in names(object$dpars)) {
    object$dpars[[dp]] <- get_dpar(object, dpar = dp)
  }
  N <- choose_N(object)
  out <- plapply(seq_len(N), log_lik_fun, cores = cores, prep = object)
  out <- do_call(cbind, out)
  colnames(out) <- NULL
  old_order <- object$old_order
  sort <- isTRUE(ncol(out) != length(old_order))
  reorder_obs(out, old_order, sort = sort)
}

# evaluate log_lik in a pointwise manner
# cannot be an S3 method since 'data_i' must be the first argument
# names must be 'data_i' and 'draws' as per ?loo::loo.function
log_lik_pointwise <- function(data_i, draws, ...) {
  i <- data_i$i
  if (is.mvbrmsprep(draws) && !length(draws$mvpars$rescor)) {
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
# @param i index of the observatio for which to compute log-lik values
# @param prep A named list returned by prepare_predictions containing
#   all required data and posterior draws
# @return a vector of length prep$ndraws containing the pointwise
#   log-likelihood for the ith observation
log_lik_gaussian <- function(i, prep) {
  mu <- get_dpar(prep, "mu", i = i)
  sigma <- get_dpar(prep, "sigma", i = i)
  sigma <- add_sigma_se(sigma, prep, i = i)
  args <- list(mean = mu, sd = sigma)
  # log_lik_censor computes the conventional log_lik in case of no censoring
  out <- log_lik_censor(dist = "norm", args = args, i = i, prep = prep)
  out <- log_lik_truncate(
    out, cdf = pnorm, args = args, i = i, prep = prep
  )
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_student <- function(i, prep) {
  nu <- get_dpar(prep, "nu", i = i)
  mu <- get_dpar(prep, "mu", i = i)
  sigma <- get_dpar(prep, "sigma", i = i)
  sigma <- add_sigma_se(sigma, prep, i = i)
  args <- list(df = nu, mu = mu, sigma = sigma)
  out <- log_lik_censor(
    dist = "student_t", args = args, i = i, prep = prep
  )
  out <- log_lik_truncate(
    out, cdf = pstudent_t, args = args, i = i, prep = prep
  )
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_lognormal <- function(i, prep) {
  sigma <- get_dpar(prep, "sigma", i = i)
  args <- list(meanlog = get_dpar(prep, "mu", i), sdlog = sigma)
  out <- log_lik_censor(dist = "lnorm", args = args, i = i, prep = prep)
  out <- log_lik_truncate(
    out, cdf = plnorm, args = args, i = i, prep = prep
  )
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_shifted_lognormal <- function(i, prep) {
  sigma <- get_dpar(prep, "sigma", i = i)
  ndt <- get_dpar(prep, "ndt", i = i)
  args <- list(meanlog = get_dpar(prep, "mu", i), sdlog = sigma, shift = ndt)
  out <- log_lik_censor("shifted_lnorm", args, i = i, prep = prep)
  out <- log_lik_truncate(out, pshifted_lnorm, args, i = i, prep = prep)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_skew_normal <- function(i, prep) {
  mu <- get_dpar(prep, "mu", i)
  sigma <- get_dpar(prep, "sigma", i = i)
  sigma <- add_sigma_se(sigma, prep, i = i)
  alpha <- get_dpar(prep, "alpha", i = i)
  args <- nlist(mu, sigma, alpha)
  out <- log_lik_censor(
    dist = "skew_normal", args = args, i = i, prep = prep
  )
  out <- log_lik_truncate(
    out, cdf = pskew_normal, args = args, i = i, prep = prep
  )
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_gaussian_mv <- function(i, prep) {
  Mu <- get_Mu(prep, i = i)
  Sigma <- get_Sigma(prep, i = i)
  dmn <- function(s) {
    dmulti_normal(
      prep$data$Y[i, ], mu = Mu[s, ],
      Sigma = Sigma[s, , ], log = TRUE
    )
  }
  out <- sapply(1:prep$ndraws, dmn)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_student_mv <- function(i, prep) {
  nu <- get_dpar(prep, "nu", i = i)
  Mu <- get_Mu(prep, i = i)
  Sigma <- get_Sigma(prep, i = i)
  dmst <- function(s) {
    dmulti_student_t(
      prep$data$Y[i, ], df = nu[s], mu = Mu[s, ],
      Sigma = Sigma[s, , ], log = TRUE
    )
  }
  out <- sapply(1:prep$ndraws, dmst)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_gaussian_time <- function(i, prep) {
  obs <- with(prep$ac, begin_tg[i]:end_tg[i])
  Jtime <- prep$ac$Jtime_tg[i, ]
  Y <- as.numeric(prep$data$Y[obs])
  mu <- as.matrix(get_dpar(prep, "mu", i = obs))
  Sigma <- get_cov_matrix_ac(prep, obs, Jtime = Jtime)
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
  rblapply(seq_len(prep$ndraws), .log_lik)
}

log_lik_student_time <- function(i, prep) {
  obs <- with(prep$ac, begin_tg[i]:end_tg[i])
  Jtime <- prep$ac$Jtime_tg[i, ]
  Y <- as.numeric(prep$data$Y[obs])
  nu <- as.matrix(get_dpar(prep, "nu", i = obs))
  mu <- as.matrix(get_dpar(prep, "mu", i = obs))
  Sigma <- get_cov_matrix_ac(prep, obs, Jtime = Jtime)
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
  rblapply(seq_len(prep$ndraws), .log_lik)
}

log_lik_gaussian_lagsar <- function(i, prep) {
  mu <- get_dpar(prep, "mu")
  sigma <- get_dpar(prep, "sigma")
  Y <- as.numeric(prep$data$Y)
  I <- diag(prep$nobs)
  stopifnot(i == 1)
  # see http://mc-stan.org/loo/articles/loo2-non-factorizable.html
  .log_lik <- function(s) {
    IB <- I - with(prep$ac, lagsar[s, ] * Msar)
    Cinv <- t(IB) %*% IB / sigma[s]^2
    e <- Y - solve(IB, mu[s, ])
    g <- Cinv %*% e
    cbar <- diag(Cinv)
    yloo <- Y - g / cbar
    sdloo <- sqrt(1 / cbar)
    ll <- dnorm(Y, yloo, sdloo, log = TRUE)
    return(as.numeric(ll))
  }
  rblapply(seq_len(prep$ndraws), .log_lik)
}

log_lik_student_lagsar <- function(i, prep) {
  nu <- get_dpar(prep, "nu")
  mu <- get_dpar(prep, "mu")
  sigma <- get_dpar(prep, "sigma")
  Y <- as.numeric(prep$data$Y)
  I <- diag(prep$nobs)
  stopifnot(i == 1)
  # see http://mc-stan.org/loo/articles/loo2-non-factorizable.html
  .log_lik <- function(s) {
    df <- nu[s]
    IB <- I - with(prep$ac, lagsar[s, ] * Msar)
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
  rblapply(seq_len(prep$ndraws), .log_lik)
}

log_lik_gaussian_errorsar <- function(i, prep) {
  stopifnot(i == 1)
  mu <- get_dpar(prep, "mu")
  sigma <- get_dpar(prep, "sigma")
  Y <- as.numeric(prep$data$Y)
  I <- diag(prep$nobs)
  .log_lik <- function(s) {
    IB <- I - with(prep$ac, errorsar[s, ] * Msar)
    Cinv <- t(IB) %*% IB / sigma[s]^2
    e <- Y - mu[s, ]
    g <- Cinv %*% e
    cbar <- diag(Cinv)
    yloo <- Y - g / cbar
    sdloo <- sqrt(1 / cbar)
    ll <- dnorm(Y, yloo, sdloo, log = TRUE)
    return(as.numeric(ll))
  }
  rblapply(seq_len(prep$ndraws), .log_lik)
}

log_lik_student_errorsar <- function(i, prep) {
  stopifnot(i == 1)
  nu <- get_dpar(prep, "nu")
  mu <- get_dpar(prep, "mu")
  sigma <- get_dpar(prep, "sigma")
  Y <- as.numeric(prep$data$Y)
  I <- diag(prep$nobs)
  .log_lik <- function(s) {
    df <- nu[s]
    IB <- I - with(prep$ac, errorsar[s, ] * Msar)
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
  rblapply(seq_len(prep$ndraws), .log_lik)
}

log_lik_gaussian_fcor <- function(i, prep) {
  stopifnot(i == 1)
  Y <- as.numeric(prep$data$Y)
  mu <- get_dpar(prep, "mu")
  Sigma <- get_cov_matrix_ac(prep)
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
  rblapply(seq_len(prep$ndraws), .log_lik)
}

log_lik_student_fcor <- function(i, prep) {
  stopifnot(i == 1)
  Y <- as.numeric(prep$data$Y)
  nu <- get_dpar(prep, "nu")
  mu <- get_dpar(prep, "mu")
  Sigma <- get_cov_matrix_ac(prep)
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
  rblapply(seq_len(prep$ndraws), .log_lik)
}

log_lik_binomial <- function(i, prep) {
  trials <- prep$data$trials[i]
  args <- list(size = trials, prob = get_dpar(prep, "mu", i))
  out <- log_lik_censor(
    dist = "binom", args = args, i = i, prep = prep
  )
  out <- log_lik_truncate(
    out, cdf = pbinom, args = args, i = i, prep = prep
  )
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_bernoulli <- function(i, prep) {
  args <- list(size = 1, prob = get_dpar(prep, "mu", i))
  out <- log_lik_censor(
    dist = "binom", args = args, i = i, prep = prep
  )
  # no truncation allowed
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_beta_binomial <- function(i, prep) {
  trials <- prep$data$trials[i]
  mu <- get_dpar(prep, "mu", i)
  phi <- get_dpar(prep, "phi", i)
  args <- nlist(size = trials, mu, phi)
  out <- log_lik_censor("beta_binomial", args, i, prep)
  out <- log_lik_truncate(out, pbeta_binomial, args, i, prep)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_poisson <- function(i, prep) {
  mu <- get_dpar(prep, "mu", i)
  mu <- multiply_dpar_rate_denom(mu, prep, i = i)
  args <- list(lambda = mu)
  out <- log_lik_censor(
    dist = "pois", args = args, i = i, prep = prep
  )
  out <- log_lik_truncate(
    out, cdf = ppois, args = args, i = i, prep = prep
  )
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_negbinomial <- function(i, prep) {
  mu <- get_dpar(prep, "mu", i)
  mu <- multiply_dpar_rate_denom(mu, prep, i = i)
  shape <- get_dpar(prep, "shape", i)
  shape <- multiply_dpar_rate_denom(shape, prep, i = i)
  args <- list(mu = mu, size = shape)
  out <- log_lik_censor(
    dist = "nbinom", args = args, i = i, prep = prep
  )
  out <- log_lik_truncate(
    out, cdf = pnbinom, args = args, i = i, prep = prep
  )
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_negbinomial2 <- function(i, prep) {
  mu <- get_dpar(prep, "mu", i)
  mu <- multiply_dpar_rate_denom(mu, prep, i = i)
  sigma <- get_dpar(prep, "sigma", i)
  shape <- multiply_dpar_rate_denom(1 / sigma, prep, i = i)
  args <- list(mu = mu, size = shape)
  out <- log_lik_censor(
    dist = "nbinom", args = args, i = i, prep = prep
  )
  out <- log_lik_truncate(
    out, cdf = pnbinom, args = args, i = i, prep = prep
  )
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_geometric <- function(i, prep) {
  mu <- get_dpar(prep, "mu", i)
  mu <- multiply_dpar_rate_denom(mu, prep, i = i)
  shape <- 1
  shape <- multiply_dpar_rate_denom(shape, prep, i = i)
  args <- list(mu = mu, size = shape)
  out <- log_lik_censor(
    dist = "nbinom", args = args, i = i, prep = prep
  )
  out <- log_lik_truncate(
    out, cdf = pnbinom, args = args, i = i, prep = prep
  )
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_discrete_weibull <- function(i, prep) {
  args <- list(
    mu = get_dpar(prep, "mu", i),
    shape = get_dpar(prep, "shape", i = i)
  )
  out <- log_lik_censor(
    dist = "discrete_weibull", args = args, i = i, prep = prep
  )
  out <- log_lik_truncate(
    out, cdf = pdiscrete_weibull, args = args, i = i, prep = prep
  )
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_com_poisson <- function(i, prep) {
  args <- list(
    mu = get_dpar(prep, "mu", i),
    shape = get_dpar(prep, "shape", i = i)
  )
  # no censoring or truncation allowed yet
  out <- do_call(dcom_poisson, c(prep$data$Y[i], args, log = TRUE))
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_exponential <- function(i, prep) {
  args <- list(rate = 1 / get_dpar(prep, "mu", i))
  out <- log_lik_censor(dist = "exp", args = args, i = i, prep = prep)
  out <- log_lik_truncate(
    out, cdf = pexp, args = args, i = i, prep = prep
  )
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_gamma <- function(i, prep) {
  shape <- get_dpar(prep, "shape", i = i)
  scale <- get_dpar(prep, "mu", i) / shape
  args <- nlist(shape, scale)
  out <- log_lik_censor(dist = "gamma", args = args, i = i, prep = prep)
  out <- log_lik_truncate(
    out, cdf = pgamma, args = args, i = i, prep = prep
  )
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_weibull <- function(i, prep) {
  shape <- get_dpar(prep, "shape", i = i)
  scale <- get_dpar(prep, "mu", i = i) / gamma(1 + 1 / shape)
  args <- list(shape = shape, scale = scale)
  out <- log_lik_censor(
    dist = "weibull", args = args, i = i, prep = prep
  )
  out <- log_lik_truncate(
    out, cdf = pweibull, args = args, i = i, prep = prep
  )
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_frechet <- function(i, prep) {
  nu <- get_dpar(prep, "nu", i = i)
  scale <- get_dpar(prep, "mu", i = i) / gamma(1 - 1 / nu)
  args <- list(scale = scale, shape = nu)
  out <- log_lik_censor(
    dist = "frechet", args = args, i = i, prep = prep
  )
  out <- log_lik_truncate(
    out, cdf = pfrechet, args = args, i = i, prep = prep
  )
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_gen_extreme_value <- function(i, prep) {
  sigma <- get_dpar(prep, "sigma", i = i)
  xi <- get_dpar(prep, "xi", i = i)
  mu <- get_dpar(prep, "mu", i)
  args <- nlist(mu, sigma, xi)
  out <- log_lik_censor(dist = "gen_extreme_value", args = args,
                       i = i, prep = prep)
  out <- log_lik_truncate(out, cdf = pgen_extreme_value,
                         args = args, i = i, prep = prep)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_inverse.gaussian <- function(i, prep) {
  args <- list(mu = get_dpar(prep, "mu", i),
               shape = get_dpar(prep, "shape", i = i))
  out <- log_lik_censor(dist = "inv_gaussian", args = args,
                       i = i, prep = prep)
  out <- log_lik_truncate(out, cdf = pinv_gaussian, args = args,
                         i = i, prep = prep)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_exgaussian <- function(i, prep) {
  args <- list(mu = get_dpar(prep, "mu", i),
               sigma = get_dpar(prep, "sigma", i = i),
               beta = get_dpar(prep, "beta", i = i))
  out <- log_lik_censor(dist = "exgaussian", args = args,
                       i = i, prep = prep)
  out <- log_lik_truncate(out, cdf = pexgaussian, args = args,
                         i = i, prep = prep)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_wiener <- function(i, prep) {
  args <- list(
    delta = get_dpar(prep, "mu", i),
    alpha = get_dpar(prep, "bs", i = i),
    tau = get_dpar(prep, "ndt", i = i),
    beta = get_dpar(prep, "bias", i = i),
    resp = prep$data[["dec"]][i]
  )
  out <- do_call(dwiener, c(prep$data$Y[i], args, log = TRUE))
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_beta <- function(i, prep) {
  mu <- get_dpar(prep, "mu", i)
  phi <- get_dpar(prep, "phi", i)
  args <- list(shape1 = mu * phi, shape2 = (1 - mu) * phi)
  out <- log_lik_censor(dist = "beta", args = args, i = i, prep = prep)
  out <- log_lik_truncate(
    out, cdf = pbeta, args = args, i = i, prep = prep
  )
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_von_mises <- function(i, prep) {
  args <- list(
    mu = get_dpar(prep, "mu", i),
    kappa = get_dpar(prep, "kappa", i = i)
  )
  out <- log_lik_censor(
    dist = "von_mises", args = args, i = i, prep = prep
  )
  out <- log_lik_truncate(
    out, cdf = pvon_mises, args = args, i = i, prep = prep
  )
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_asym_laplace <- function(i, prep, ...) {
  args <- list(
    mu = get_dpar(prep, "mu", i),
    sigma = get_dpar(prep, "sigma", i),
    quantile = get_dpar(prep, "quantile", i)
  )
  out <- log_lik_censor(dist = "asym_laplace", args, i, prep)
  out <- log_lik_truncate(out, pasym_laplace, args, i, prep)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_zero_inflated_asym_laplace <- function(i, prep, ...) {
  args <- list(
    mu = get_dpar(prep, "mu", i),
    sigma = get_dpar(prep, "sigma", i),
    quantile = get_dpar(prep, "quantile", i),
    zi = get_dpar(prep, "zi", i)
  )
  out <- log_lik_censor(dist = "zero_inflated_asym_laplace", args, i, prep)
  out <- log_lik_truncate(out, pzero_inflated_asym_laplace, args, i, prep)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_cox <- function(i, prep, ...) {
  args <- list(
    mu = get_dpar(prep, "mu", i),
    bhaz = prep$bhaz$bhaz[, i],
    cbhaz = prep$bhaz$cbhaz[, i]
  )
  out <- log_lik_censor(dist = "cox", args = args, i = i, prep = prep)
  out <- log_lik_truncate(out, cdf = pcox, args = args, i = i, prep = prep)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_hurdle_poisson <- function(i, prep) {
  hu <- get_dpar(prep, "hu", i)
  lambda <- get_dpar(prep, "mu", i)
  args <- nlist(lambda, hu)
  out <- log_lik_censor("hurdle_poisson", args, i, prep)
  out <- log_lik_truncate(out, phurdle_poisson, args, i, prep)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_hurdle_negbinomial <- function(i, prep) {
  hu <- get_dpar(prep, "hu", i)
  mu <- get_dpar(prep, "mu", i)
  shape <- get_dpar(prep, "shape", i = i)
  args <- nlist(mu, shape, hu)
  out <- log_lik_censor("hurdle_negbinomial", args, i, prep)
  out <- log_lik_truncate(out, phurdle_negbinomial, args, i, prep)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_hurdle_gamma <- function(i, prep) {
  hu <- get_dpar(prep, "hu", i)
  shape <- get_dpar(prep, "shape", i = i)
  scale <- get_dpar(prep, "mu", i) / shape
  args <- nlist(shape, scale, hu)
  out <- log_lik_censor("hurdle_gamma", args, i, prep)
  out <- log_lik_truncate(out, phurdle_gamma, args, i, prep)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_hurdle_lognormal <- function(i, prep) {
  hu <- get_dpar(prep, "hu", i)
  mu <- get_dpar(prep, "mu", i)
  sigma <- get_dpar(prep, "sigma", i = i)
  args <- nlist(mu, sigma, hu)
  out <- log_lik_censor("hurdle_lognormal", args, i, prep)
  out <- log_lik_truncate(out, phurdle_lognormal, args, i, prep)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_hurdle_cumulative <- function(i, prep) {
  mu <- get_dpar(prep, "mu", i = i)
  hu <- get_dpar(prep, "hu", i = i)
  disc <- get_dpar(prep, "disc", i = i)
  thres <- subset_thres(prep, i)
  nthres <- NCOL(thres)
  eta <- disc * (thres - mu)
  y <- prep$data$Y[i]
  if (y == 0L) {
    out <- dbinom(1, size = 1, prob = hu, log = TRUE)
  } else if (y == 1L) {
    out <- log_cdf(eta[, 1L], prep$family$link) +
      dbinom(0, size = 1, prob = hu, log = TRUE)
  } else if (y == nthres + 1L) {
    out <- log_ccdf(eta[, y - 1L], prep$family$link) +
      dbinom(0, size = 1, prob = hu, log = TRUE)
  } else {
    out <- log_diff_exp(
      log_cdf(eta[, y], prep$family$link),
      log_cdf(eta[, y - 1L], prep$family$link)
    ) + dbinom(0, size = 1, prob = hu, log = TRUE)
  }
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_mixcure_lognormal <- function(i, prep) {
  mu <- get_dpar(prep, "mu", i)
  sigma <- get_dpar(prep, "sigma", i = i)
  inc <- get_dpar(prep, "inc", i)
  args <- nlist(mu = mu, sigma = sigma, inc = inc)
  out <- log_lik_censor("mixcure_lognormal", args, i, prep)
  out <- log_lik_truncate(out, pmixcure_lognormal, args, i, prep)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_mixcure_weibull <- function(i, prep) {
  shape <- get_dpar(prep, "shape", i = i)
  scale <- get_dpar(prep, "mu", i = i) / gamma(1 + 1 / shape)
  inc <- get_dpar(prep, "inc", i)
  args <- list(shape = shape, scale = scale, inc = inc)
  out <- log_lik_censor(
    dist = "mixcure_weibull", args = args, i = i, prep = prep
  )
  out <- log_lik_truncate(
    out, cdf = pmixcure_weibull, args = args, i = i, prep = prep
  )
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_zero_inflated_poisson <- function(i, prep) {
  zi <- get_dpar(prep, "zi", i)
  lambda <- get_dpar(prep, "mu", i)
  args <- nlist(lambda, zi)
  out <- log_lik_censor("zero_inflated_poisson", args, i, prep)
  out <- log_lik_truncate(out, pzero_inflated_poisson, args, i, prep)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_zero_inflated_negbinomial <- function(i, prep) {
  zi <- get_dpar(prep, "zi", i)
  mu <- get_dpar(prep, "mu", i)
  shape <- get_dpar(prep, "shape", i = i)
  args <- nlist(mu, shape, zi)
  out <- log_lik_censor("zero_inflated_negbinomial", args, i, prep)
  out <- log_lik_truncate(out, pzero_inflated_negbinomial, args, i, prep)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_zero_inflated_binomial <- function(i, prep) {
  trials <- prep$data$trials[i]
  mu <- get_dpar(prep, "mu", i)
  zi <- get_dpar(prep, "zi", i)
  args <- list(size = trials, prob = mu, zi = zi)
  out <- log_lik_censor("zero_inflated_binomial", args, i, prep)
  out <- log_lik_truncate(out, pzero_inflated_binomial, args, i, prep)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_zero_inflated_beta_binomial <- function(i, prep) {
  trials <- prep$data$trials[i]
  mu <- get_dpar(prep, "mu", i)
  phi <- get_dpar(prep, "phi", i)
  zi <- get_dpar(prep, "zi", i)
  args <- nlist(size = trials, mu, phi, zi)
  out <- log_lik_censor("zero_inflated_beta_binomial", args, i, prep)
  out <- log_lik_truncate(out, pzero_inflated_beta_binomial, args, i, prep)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_zero_inflated_beta <- function(i, prep) {
  zi <- get_dpar(prep, "zi", i)
  mu <- get_dpar(prep, "mu", i)
  phi <- get_dpar(prep, "phi", i)
  args <- nlist(shape1 = mu * phi, shape2 = (1 - mu) * phi, zi)
  out <- log_lik_censor("zero_inflated_beta", args, i, prep)
  out <- log_lik_truncate(out, pzero_inflated_beta, args, i, prep)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_zero_one_inflated_beta <- function(i, prep) {
  zoi <- get_dpar(prep, "zoi", i)
  coi <- get_dpar(prep, "coi", i)
  if (prep$data$Y[i] %in% c(0, 1)) {
    out <- dbinom(1, size = 1, prob = zoi, log = TRUE) +
      dbinom(prep$data$Y[i], size = 1, prob = coi, log = TRUE)
  } else {
    phi <- get_dpar(prep, "phi", i)
    mu <- get_dpar(prep, "mu", i)
    args <- list(shape1 = mu * phi, shape2 = (1 - mu) * phi)
    out <- dbinom(0, size = 1, prob = zoi, log = TRUE) +
      do_call(dbeta, c(prep$data$Y[i], args, log = TRUE))
  }
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_categorical <- function(i, prep) {
  stopifnot(prep$family$link == "logit")
  eta <- get_Mu(prep, i = i)
  eta <- insert_refcat(eta, refcat = prep$refcat)
  out <- dcategorical(prep$data$Y[i], eta = eta, log = TRUE)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_multinomial <- function(i, prep) {
  stopifnot(prep$family$link == "logit")
  eta <- get_Mu(prep, i = i)
  eta <- insert_refcat(eta, refcat = prep$refcat)
  out <- dmultinomial(prep$data$Y[i, ], eta = eta, log = TRUE)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_dirichlet <- function(i, prep) {
  stopifnot(prep$family$link == "logit")
  eta <- get_Mu(prep, i = i)
  eta <- insert_refcat(eta, refcat = prep$refcat)
  phi <- get_dpar(prep, "phi", i = i)
  cats <- seq_len(prep$data$ncat)
  alpha <- dcategorical(cats, eta = eta) * phi
  out <- ddirichlet(prep$data$Y[i, ], alpha = alpha, log = TRUE)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_dirichlet2 <- function(i, prep) {
  mu <- get_Mu(prep, i = i)
  out <- ddirichlet(prep$data$Y[i, ], alpha = mu, log = TRUE)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_logistic_normal <- function(i, prep, ...) {
  mu <- get_Mu(prep, i = i)
  Sigma <- get_Sigma(prep, i = i, cor_name = "lncor")
  dlmn <- function(s) {
    dlogistic_normal(
      prep$data$Y[i, ], mu = mu[s, ], Sigma = Sigma[s, , ],
      refcat = prep$refcat, log = TRUE
    )
  }
  out <- sapply(1:prep$ndraws, dlmn)
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_cumulative <- function(i, prep) {
  disc <- get_dpar(prep, "disc", i = i)
  mu <- get_dpar(prep, "mu", i = i)
  thres <- subset_thres(prep, i)
  nthres <- NCOL(thres)
  eta <- disc * (thres - mu)
  y <- prep$data$Y[i]
  if (y == 1L) {
    out <- log_cdf(eta[, 1L], prep$family$link)
  } else if (y == nthres + 1L) {
    out <- log_ccdf(eta[, y - 1L], prep$family$link)
  } else {
    out <- log_diff_exp(
      log_cdf(eta[, y], prep$family$link),
      log_cdf(eta[, y - 1L], prep$family$link)
    )
  }
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_sratio <- function(i, prep) {
  disc <- get_dpar(prep, "disc", i = i)
  mu <- get_dpar(prep, "mu", i = i)
  thres <- subset_thres(prep, i)
  nthres <- NCOL(thres)
  eta <- disc * (thres - mu)
  y <- prep$data$Y[i]
  q <- sapply(seq_len(min(y, nthres)),
    function(k) log_ccdf(eta[, k], prep$family$link)
  )
  if (y == 1L) {
    out <- log1m_exp(q[, 1L])
  } else if (y == 2L) {
    out <- log1m_exp(q[, 2L]) + q[, 1L]
  } else if (y == nthres + 1L) {
    out <- rowSums(q)
  } else {
    out <- log1m_exp(q[, y]) + rowSums(q[, 1L:(y - 1L)])
  }
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_cratio <- function(i, prep) {
  disc <- get_dpar(prep, "disc", i = i)
  mu <- get_dpar(prep, "mu", i = i)
  thres <- subset_thres(prep, i)
  nthres <- NCOL(thres)
  eta <- disc * (mu - thres)
  y <- prep$data$Y[i]
  q <- sapply(seq_len(min(y, nthres)),
    function(k) log_cdf(eta[, k], prep$family$link)
  )
  if (y == 1L) {
    out <- log1m_exp(q[, 1L])
  } else if (y == 2L) {
    out <- log1m_exp(q[, 2L]) + q[, 1L]
  } else if (y == nthres + 1L) {
    out <- rowSums(q)
  } else {
    out <- log1m_exp(q[, y]) + rowSums(q[, 1L:(y - 1L)])
  }
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_acat <- function(i, prep) {
  disc <- get_dpar(prep, "disc", i = i)
  mu <- get_dpar(prep, "mu", i = i)
  thres <- subset_thres(prep, i)
  nthres <- NCOL(thres)
  eta <- disc * (mu - thres)
  y <- prep$data$Y[i]
  # TODO: check if computation can be made more numerically stable
  if (prep$family$link == "logit") {
    # more efficient computation for logit link
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
      inv_link(eta[, k], prep$family$link))
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
  log_lik_weight(out, i = i, prep = prep)
}

log_lik_custom <- function(i, prep) {
  custom_family_method(prep$family, "log_lik")(i, prep)
}

log_lik_mixture <- function(i, prep) {
  families <- family_names(prep$family)
  theta <- get_theta(prep, i = i)
  out <- array(NA, dim = dim(theta))
  for (j in seq_along(families)) {
    log_lik_fun <- paste0("log_lik_", families[j])
    log_lik_fun <- get(log_lik_fun, asNamespace("brms"))
    tmp_draws <- pseudo_prep_for_mixture(prep, j)
    out[, j] <- exp(log(theta[, j]) + log_lik_fun(i, tmp_draws))
  }
  if (isTRUE(prep[["pp_mixture"]])) {
    out <- log(out) - log(rowSums(out))
  } else {
    out <- log(rowSums(out))
  }
  log_lik_weight(out, i = i, prep = prep)
}

# ----------- log_lik helper-functions -----------
# compute (possibly censored) log_lik values
# @param dist name of a distribution for which the functions
#   d<dist> (pdf) and p<dist> (cdf) are available
# @param args additional arguments passed to pdf and cdf
# @param prep a brmsprep object
# @return vector of log_lik values
log_lik_censor <- function(dist, args, i, prep) {
  pdf <- get(paste0("d", dist), mode = "function")
  cdf <- get(paste0("p", dist), mode = "function")
  y <- prep$data$Y[i]
  cens <- prep$data$cens[i]
  if (is.null(cens) || cens == 0) {
    x <- do_call(pdf, c(y, args, log = TRUE))
  } else if (cens == 1) {
    x <- do_call(cdf, c(y, args, lower.tail = FALSE, log.p = TRUE))
  } else if (cens == -1) {
    x <- do_call(cdf, c(y, args, log.p = TRUE))
  } else if (cens == 2) {
    rcens <- prep$data$rcens[i]
    x <- log(do_call(cdf, c(rcens, args)) - do_call(cdf, c(y, args)))
  }
  x
}

# adjust log_lik in truncated models
# @param x vector of log_lik values
# @param cdf a cumulative distribution function
# @param args arguments passed to cdf
# @param i observation number
# @param prep a brmsprep object
# @return vector of log_lik values
log_lik_truncate <- function(x, cdf, args, i, prep) {
  lb <- prep$data[["lb"]][i]
  ub <- prep$data[["ub"]][i]
  if (is.null(lb) && is.null(ub)) {
    return(x)
  }
  if (!is.null(lb)) {
    log_cdf_lb <- do_call(cdf, c(lb, args, log.p = TRUE))
  } else {
    log_cdf_lb <- rep(-Inf, length(x))
  }
  if (!is.null(ub)) {
    log_cdf_ub <- do_call(cdf, c(ub, args, log.p = TRUE))
  } else {
    log_cdf_ub <- rep(0, length(x))
  }
  x - log_diff_exp(log_cdf_ub, log_cdf_lb)
}

# weight log_lik values according to defined weights
# @param x vector of log_lik values
# @param i observation number
# @param prep a brmsprep object
# @return vector of log_lik values
log_lik_weight <- function(x, i, prep) {
  weight <- prep$data$weights[i]
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
