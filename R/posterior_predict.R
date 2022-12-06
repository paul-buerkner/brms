#' Draws from the Posterior Predictive Distribution
#'
#' Compute posterior draws of the posterior predictive distribution. Can be
#' performed for the data used to fit the model (posterior predictive checks) or
#' for new data. By definition, these draws have higher variance than draws
#' of the expected value of the posterior predictive distribution computed by
#' \code{\link{posterior_epred.brmsfit}}. This is because the residual error
#' is incorporated in \code{posterior_predict}. However, the estimated means of
#' both methods averaged across draws should be very similar.
#'
#' @inheritParams prepare_predictions
#' @param object An object of class \code{brmsfit}.
#' @param re.form Alias of \code{re_formula}.
#' @param transform (Deprecated) A function or a character string naming
#'   a function to be applied on the predicted responses
#'   before summary statistics are computed.
#' @param negative_rt Only relevant for Wiener diffusion models.
#'   A flag indicating whether response times of responses
#'   on the lower boundary should be returned as negative values.
#'   This allows to distinguish responses on the upper and
#'   lower boundary. Defaults to \code{FALSE}.
#' @param sort Logical. Only relevant for time series models.
#'  Indicating whether to return predicted values in the original
#'  order (\code{FALSE}; default) or in the order of the
#'  time series (\code{TRUE}).
#' @param ntrys Parameter used in rejection sampling
#'   for truncated discrete models only
#'   (defaults to \code{5}). See Details for more information.
#' @param cores Number of cores (defaults to \code{1}). On non-Windows systems,
#'   this argument can be set globally via the \code{mc.cores} option.
#' @param ... Further arguments passed to \code{\link{prepare_predictions}}
#'   that control several aspects of data validation and prediction.
#'
#' @return An \code{array} of draws. In univariate models,
#'   the output is as an S x N matrix, where S is the number of posterior
#'   draws and N is the number of observations. In multivariate models, an
#'   additional dimension is added to the output which indexes along the
#'   different response variables.
#'
#' @template details-newdata-na
#' @template details-allow_new_levels
#' @details For truncated discrete models only: In the absence of any general
#'   algorithm to sample from truncated discrete distributions, rejection
#'   sampling is applied in this special case. This means that values are
#'   sampled until a value lies within the defined truncation boundaries. In
#'   practice, this procedure may be rather slow (especially in \R). Thus, we
#'   try to do approximate rejection sampling by sampling each value
#'   \code{ntrys} times and then select a valid value. If all values are
#'   invalid, the closest boundary is used, instead. If there are more than a
#'   few of these pathological cases, a warning will occur suggesting to
#'   increase argument \code{ntrys}.
#'
#' @examples
#' \dontrun{
#' ## fit a model
#' fit <- brm(time | cens(censored) ~ age + sex + (1 + age || patient),
#'            data = kidney, family = "exponential", init = "0")
#'
#' ## predicted responses
#' pp <- posterior_predict(fit)
#' str(pp)
#'
#' ## predicted responses excluding the group-level effect of age
#' pp <- posterior_predict(fit, re_formula = ~ (1 | patient))
#' str(pp)
#'
#' ## predicted responses of patient 1 for new data
#' newdata <- data.frame(
#'   sex = factor(c("male", "female")),
#'   age = c(20, 50),
#'   patient = c(1, 1)
#' )
#' pp <- posterior_predict(fit, newdata = newdata)
#' str(pp)
#' }
#'
#' @aliases posterior_predict
#' @method posterior_predict brmsfit
#' @importFrom rstantools posterior_predict
#' @export
#' @export posterior_predict
posterior_predict.brmsfit <- function(
  object, newdata = NULL, re_formula = NULL, re.form = NULL,
  transform = NULL, resp = NULL, negative_rt = FALSE,
  ndraws = NULL, draw_ids = NULL, sort = FALSE, ntrys = 5,
  cores = NULL, ...
) {
  cl <- match.call()
  if ("re.form" %in% names(cl) && !missing(re.form)) {
    re_formula <- re.form
  }
  contains_draws(object)
  object <- restructure(object)
  prep <- prepare_predictions(
    object, newdata = newdata, re_formula = re_formula, resp = resp,
    ndraws = ndraws, draw_ids = draw_ids, check_response = FALSE, ...
  )
  posterior_predict(
    prep, transform = transform, sort = sort, ntrys = ntrys,
    negative_rt = negative_rt, cores = cores, summary = FALSE
  )
}

#' @export
posterior_predict.mvbrmsprep <- function(object, ...) {
  if (length(object$mvpars$rescor)) {
    object$mvpars$Mu <- get_Mu(object)
    object$mvpars$Sigma <- get_Sigma(object)
    out <- posterior_predict.brmsprep(object, ...)
  } else {
    out <- lapply(object$resps, posterior_predict, ...)
    along <- ifelse(length(out) > 1L, 3, 2)
    out <- do_call(abind, c(out, along = along))
  }
  out
}

#' @export
posterior_predict.brmsprep <- function(object, transform = NULL, sort = FALSE,
                                       summary = FALSE, robust = FALSE,
                                       probs = c(0.025, 0.975),
                                       cores = NULL, ...) {
  summary <- as_one_logical(summary)
  cores <- validate_cores_post_processing(cores)
  if (is.customfamily(object$family)) {
    # ensure that the method can be found during parallel execution
    object$family$posterior_predict <-
      custom_family_method(object$family, "posterior_predict")
  }
  for (nlp in names(object$nlpars)) {
    object$nlpars[[nlp]] <- get_nlpar(object, nlpar = nlp)
  }
  for (dp in names(object$dpars)) {
    object$dpars[[dp]] <- get_dpar(object, dpar = dp)
  }
  pp_fun <- paste0("posterior_predict_", object$family$fun)
  pp_fun <- get(pp_fun, asNamespace("brms"))
  N <- choose_N(object)
  out <- plapply(seq_len(N), pp_fun, cores = cores, prep = object, ...)
  if (grepl("_mv$", object$family$fun)) {
    out <- do_call(abind, c(out, along = 3))
    out <- aperm(out, perm = c(1, 3, 2))
    dimnames(out)[[3]] <- names(object$resps)
  } else if (has_multicol(object$family)) {
    out <- do_call(abind, c(out, along = 3))
    out <- aperm(out, perm = c(1, 3, 2))
    dimnames(out)[[3]] <- object$cats
  } else {
    out <- do_call(cbind, out)
  }
  colnames(out) <- rownames(out) <- NULL
  if (use_int(object$family)) {
    out <- check_discrete_trunc_bounds(
      out, lb = object$data$lb, ub = object$data$ub
    )
  }
  out <- reorder_obs(out, object$old_order, sort = sort)
  # transform predicted response draws before summarizing them
  if (!is.null(transform)) {
    # deprecated as of brms 2.12.3
    warning2("Argument 'transform' is deprecated ",
             "and will be removed in the future.")
    out <- do_call(transform, list(out))
  }
  attr(out, "levels") <- object$cats
  if (summary) {
    # only for compatibility with the 'predict' method
    if (is_ordinal(object$family)) {
      levels <- seq_len(max(object$data$nthres) + 1)
      out <- posterior_table(out, levels = levels)
    } else if (is_categorical(object$family)) {
      levels <- seq_len(object$data$ncat)
      out <- posterior_table(out, levels = levels)
    } else {
      out <- posterior_summary(out, probs = probs, robust = robust)
    }
  }
  out
}

#' Draws from the Posterior Predictive Distribution
#'
#' This method is an alias of \code{\link{posterior_predict.brmsfit}}
#' with additional arguments for obtaining summaries of the computed draws.
#'
#' @inheritParams posterior_predict.brmsfit
#' @param summary Should summary statistics be returned
#'  instead of the raw values? Default is \code{TRUE}.
#' @param robust If \code{FALSE} (the default) the mean is used as
#'  the measure of central tendency and the standard deviation as
#'  the measure of variability. If \code{TRUE}, the median and the
#'  median absolute deviation (MAD) are applied instead.
#'  Only used if \code{summary} is \code{TRUE}.
#' @param probs  The percentiles to be computed by the \code{quantile}
#'  function. Only used if \code{summary} is \code{TRUE}.
#'
#' @return An \code{array} of predicted response values.
#'   If \code{summary = FALSE} the output resembles those of
#'   \code{\link{posterior_predict.brmsfit}}.
#'
#'   If \code{summary = TRUE} the output depends on the family: For categorical
#'   and ordinal families, the output is an N x C matrix, where N is the number
#'   of observations, C is the number of categories, and the values are
#'   predicted category probabilities. For all other families, the output is a N
#'   x E matrix where E = \code{2 + length(probs)} is the number of summary
#'   statistics: The \code{Estimate} column contains point estimates (either
#'   mean or median depending on argument \code{robust}), while the
#'   \code{Est.Error} column contains uncertainty estimates (either standard
#'   deviation or median absolute deviation depending on argument
#'   \code{robust}). The remaining columns starting with \code{Q} contain
#'   quantile estimates as specified via argument \code{probs}.
#'
#' @seealso \code{\link{posterior_predict.brmsfit}}
#'
#' @examples
#' \dontrun{
#' ## fit a model
#' fit <- brm(time | cens(censored) ~ age + sex + (1 + age || patient),
#'            data = kidney, family = "exponential", init = "0")
#'
#' ## predicted responses
#' pp <- predict(fit)
#' head(pp)
#'
#' ## predicted responses excluding the group-level effect of age
#' pp <- predict(fit, re_formula = ~ (1 | patient))
#' head(pp)
#'
#' ## predicted responses of patient 1 for new data
#' newdata <- data.frame(
#'   sex = factor(c("male", "female")),
#'   age = c(20, 50),
#'   patient = c(1, 1)
#' )
#' predict(fit, newdata = newdata)
#' }
#'
#' @export
predict.brmsfit <- function(object, newdata = NULL, re_formula = NULL,
                            transform = NULL, resp = NULL,
                            negative_rt = FALSE, ndraws = NULL, draw_ids = NULL,
                            sort = FALSE, ntrys = 5, cores = NULL, summary = TRUE,
                            robust = FALSE, probs = c(0.025, 0.975), ...) {
  contains_draws(object)
  object <- restructure(object)
  prep <- prepare_predictions(
    object, newdata = newdata, re_formula = re_formula, resp = resp,
    ndraws = ndraws, draw_ids = draw_ids, check_response = FALSE, ...
  )
  posterior_predict(
    prep, transform = transform, ntrys = ntrys, negative_rt = negative_rt,
    sort = sort, cores = cores, summary = summary, robust = robust,
    probs = probs
  )
}

#' Predictive Intervals
#'
#' Compute intervals from the posterior predictive distribution.
#'
#' @aliases predictive_interval
#'
#' @param object An \R object of class \code{brmsfit}.
#' @param prob A number p (0 < p < 1) indicating the desired probability mass to
#'   include in the intervals. Defaults to \code{0.9}.
#' @param ... Further arguments passed to \code{\link{posterior_predict}}.
#'
#' @return A matrix with 2 columns for the lower and upper bounds of the
#'   intervals, respectively, and as many rows as observations being predicted.
#'
#' @examples
#' \dontrun{
#' fit <- brm(count ~ zBase, data = epilepsy, family = poisson())
#' predictive_interval(fit)
#' }
#'
#' @importFrom rstantools predictive_interval
#' @export predictive_interval
#' @export
predictive_interval.brmsfit <- function(object, prob = 0.9, ...) {
  out <- posterior_predict(object, ...)
  predictive_interval(out, prob = prob)
}

# validate method name to obtain posterior predictions
# @param method name of the method
# @return validated name of the method
validate_pp_method <- function(method) {
  method <- as_one_character(method)
  if (method %in% c("posterior_predict", "predict", "pp")) {
    method <- "posterior_predict"
  } else if (method %in% c("posterior_epred", "fitted", "pp_expect")) {
    method <- "posterior_epred"
  } else if (method %in% c("posterior_linpred")) {
    method <- "posterior_linpred"
  } else if (method %in% c("predictive_error", "residuals")) {
    method <- "predictive_error"
  } else {
    stop2("Posterior predictive method '", method, "' it not supported.")
  }
  method
}

# ------------------- family specific posterior_predict methods ---------------------
# All posterior_predict_<family> functions have the same arguments structure
# @param i index of the observatio for which to compute pp values
# @param prep A named list returned by prepare_predictions containing
#   all required data and posterior draws
# @param ... ignored arguments
# @param A vector of length prep$ndraws containing draws
#   from the posterior predictive distribution
posterior_predict_gaussian <- function(i, prep, ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i = i)
  sigma <- get_dpar(prep, "sigma", i = i)
  sigma <- add_sigma_se(sigma, prep, i = i)
  rcontinuous(
    n = prep$ndraws, dist = "norm",
    mean = mu, sd = sigma,
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_student <- function(i, prep, ntrys = 5, ...) {
  nu <- get_dpar(prep, "nu", i = i)
  mu <- get_dpar(prep, "mu", i = i)
  sigma <- get_dpar(prep, "sigma", i = i)
  sigma <- add_sigma_se(sigma, prep, i = i)
  rcontinuous(
    n = prep$ndraws, dist = "student_t",
    df = nu, mu = mu, sigma = sigma,
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_lognormal <- function(i, prep, ntrys = 5, ...) {
  rcontinuous(
    n = prep$ndraws, dist = "lnorm",
    meanlog = get_dpar(prep, "mu", i = i),
    sdlog = get_dpar(prep, "sigma", i = i),
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_shifted_lognormal <- function(i, prep, ntrys = 5, ...) {
  rcontinuous(
    n = prep$ndraws, dist = "shifted_lnorm",
    meanlog = get_dpar(prep, "mu", i = i),
    sdlog = get_dpar(prep, "sigma", i = i),
    shift = get_dpar(prep, "ndt", i = i),
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_skew_normal <- function(i, prep, ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i = i)
  sigma <- get_dpar(prep, "sigma", i = i)
  sigma <- add_sigma_se(sigma, prep, i = i)
  alpha <- get_dpar(prep, "alpha", i = i)
  rcontinuous(
    n = prep$ndraws, dist = "skew_normal",
    mu = mu, sigma = sigma, alpha = alpha,
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_gaussian_mv <- function(i, prep, ...) {
  Mu <- get_Mu(prep, i = i)
  Sigma <- get_Sigma(prep, i = i)
  .predict <- function(s) {
    rmulti_normal(1, mu = Mu[s, ], Sigma = Sigma[s, , ])
  }
  rblapply(seq_len(prep$ndraws), .predict)
}

posterior_predict_student_mv <- function(i, prep, ...) {
  nu <- get_dpar(prep, "nu", i = i)
  Mu <- get_Mu(prep, i = i)
  Sigma <- get_Sigma(prep, i = i)
  .predict <- function(s) {
    rmulti_student_t(1, df = nu[s], mu = Mu[s, ], Sigma = Sigma[s, , ])
  }
  rblapply(seq_len(prep$ndraws), .predict)
}

posterior_predict_gaussian_time <- function(i, prep, ...) {
  obs <- with(prep$ac, begin_tg[i]:end_tg[i])
  Jtime <- prep$ac$Jtime_tg[i, ]
  mu <- as.matrix(get_dpar(prep, "mu", i = obs))
  Sigma <- get_cov_matrix_ac(prep, obs, Jtime = Jtime)
  .predict <- function(s) {
    rmulti_normal(1, mu = mu[s, ], Sigma = Sigma[s, , ])
  }
  rblapply(seq_len(prep$ndraws), .predict)
}

posterior_predict_student_time <- function(i, prep, ...) {
  obs <- with(prep$ac, begin_tg[i]:end_tg[i])
  Jtime <- prep$ac$Jtime_tg[i, ]
  nu <- as.matrix(get_dpar(prep, "nu", i = obs))
  mu <- as.matrix(get_dpar(prep, "mu", i = obs))
  Sigma <- get_cov_matrix_ac(prep, obs, Jtime = Jtime)
  .predict <- function(s) {
    rmulti_student_t(1, df = nu[s, ], mu = mu[s, ], Sigma = Sigma[s, , ])
  }
  rblapply(seq_len(prep$ndraws), .predict)
}

posterior_predict_gaussian_lagsar <- function(i, prep, ...) {
  stopifnot(i == 1)
  .predict <- function(s) {
    M_new <- with(prep, diag(nobs) - ac$lagsar[s] * ac$Msar)
    mu <- as.numeric(solve(M_new) %*% mu[s, ])
    Sigma <- solve(crossprod(M_new)) * sigma[s]^2
    rmulti_normal(1, mu = mu, Sigma = Sigma)
  }
  mu <- get_dpar(prep, "mu")
  sigma <- get_dpar(prep, "sigma")
  rblapply(seq_len(prep$ndraws), .predict)
}

posterior_predict_student_lagsar <- function(i, prep, ...) {
  stopifnot(i == 1)
  .predict <- function(s) {
    M_new <- with(prep, diag(nobs) - ac$lagsar[s] * ac$Msar)
    mu <- as.numeric(solve(M_new) %*% mu[s, ])
    Sigma <- solve(crossprod(M_new)) * sigma[s]^2
    rmulti_student_t(1, df = nu[s], mu = mu, Sigma = Sigma)
  }
  mu <- get_dpar(prep, "mu")
  sigma <- get_dpar(prep, "sigma")
  nu <- get_dpar(prep, "nu")
  rblapply(seq_len(prep$ndraws), .predict)
}

posterior_predict_gaussian_errorsar <- function(i, prep, ...) {
  stopifnot(i == 1)
  .predict <- function(s) {
    M_new <- with(prep, diag(nobs) - ac$errorsar[s] * ac$Msar)
    Sigma <- solve(crossprod(M_new)) * sigma[s]^2
    rmulti_normal(1, mu = mu[s, ], Sigma = Sigma)
  }
  mu <- get_dpar(prep, "mu")
  sigma <- get_dpar(prep, "sigma")
  rblapply(seq_len(prep$ndraws), .predict)
}

posterior_predict_student_errorsar <- function(i, prep, ...) {
  stopifnot(i == 1)
  .predict <- function(s) {
    M_new <- with(prep, diag(nobs) - ac$errorsar[s] * ac$Msar)
    Sigma <- solve(crossprod(M_new)) * sigma[s]^2
    rmulti_student_t(1, df = nu[s], mu = mu[s, ], Sigma = Sigma)
  }
  mu <- get_dpar(prep, "mu")
  sigma <- get_dpar(prep, "sigma")
  nu <- get_dpar(prep, "nu")
  rblapply(seq_len(prep$ndraws), .predict)
}

posterior_predict_gaussian_fcor <- function(i, prep, ...) {
  stopifnot(i == 1)
  mu <- as.matrix(get_dpar(prep, "mu"))
  Sigma <- get_cov_matrix_ac(prep)
  .predict <- function(s) {
    rmulti_normal(1, mu = mu[s, ], Sigma = Sigma[s, , ])
  }
  rblapply(seq_len(prep$ndraws), .predict)
}

posterior_predict_student_fcor <- function(i, prep, ...) {
  stopifnot(i == 1)
  nu <- as.matrix(get_dpar(prep, "nu"))
  mu <- as.matrix(get_dpar(prep, "mu"))
  Sigma <- get_cov_matrix_ac(prep)
  .predict <- function(s) {
    rmulti_student_t(1, df = nu[s, ], mu = mu[s, ], Sigma = Sigma[s, , ])
  }
  rblapply(seq_len(prep$ndraws), .predict)
}

posterior_predict_binomial <- function(i, prep, ntrys = 5, ...) {
  rdiscrete(
    n = prep$ndraws, dist = "binom",
    size = prep$data$trials[i],
    prob = get_dpar(prep, "mu", i = i),
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_beta_binomial <- function(i, prep, ntrys = 5, ...) {
  rdiscrete(
    n = prep$ndraws, dist = "beta_binomial",
    size = prep$data$trials[i],
    mu = get_dpar(prep, "mu", i = i),
    phi = get_dpar(prep, "phi", i = i),
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_bernoulli <- function(i, prep, ...) {
  mu <- get_dpar(prep, "mu", i = i)
  rbinom(length(mu), size = 1, prob = mu)
}

posterior_predict_poisson <- function(i, prep, ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i)
  mu <- multiply_dpar_rate_denom(mu, prep, i = i)
  rdiscrete(
    n = prep$ndraws, dist = "pois", lambda = mu,
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_negbinomial <- function(i, prep, ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i)
  mu <- multiply_dpar_rate_denom(mu, prep, i = i)
  shape <- get_dpar(prep, "shape", i)
  shape <- multiply_dpar_rate_denom(shape, prep, i = i)
  rdiscrete(
    n = prep$ndraws, dist = "nbinom",
    mu = mu, size = shape,
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_negbinomial2 <- function(i, prep, ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i)
  mu <- multiply_dpar_rate_denom(mu, prep, i = i)
  sigma <- get_dpar(prep, "sigma", i)
  shape <- multiply_dpar_rate_denom(1 / sigma, prep, i = i)
  rdiscrete(
    n = prep$ndraws, dist = "nbinom",
    mu = mu, size = shape,
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_geometric <- function(i, prep, ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i)
  mu <- multiply_dpar_rate_denom(mu, prep, i = i)
  shape <- 1
  shape <- multiply_dpar_rate_denom(shape, prep, i = i)
  rdiscrete(
    n = prep$ndraws, dist = "nbinom",
    mu = mu, size = shape,
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_discrete_weibull <- function(i, prep, ntrys = 5, ...) {
  rdiscrete(
    n = prep$ndraws, dist = "discrete_weibull",
    mu = get_dpar(prep, "mu", i = i),
    shape = get_dpar(prep, "shape", i = i),
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_com_poisson <- function(i, prep, ntrys = 5, ...) {
  rdiscrete(
    n = prep$ndraws, dist = "com_poisson",
    mu = get_dpar(prep, "mu", i = i),
    shape = get_dpar(prep, "shape", i = i),
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_exponential <- function(i, prep, ntrys = 5, ...) {
  rcontinuous(
    n = prep$ndraws, dist = "exp",
    rate = 1 / get_dpar(prep, "mu", i = i),
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_gamma <- function(i, prep, ntrys = 5, ...) {
  shape <- get_dpar(prep, "shape", i = i)
  scale <- get_dpar(prep, "mu", i = i) / shape
  rcontinuous(
    n = prep$ndraws, dist = "gamma",
    shape = shape, scale = scale,
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_weibull <- function(i, prep, ntrys = 5, ...) {
  shape <- get_dpar(prep, "shape", i = i)
  scale <- get_dpar(prep, "mu", i = i) / gamma(1 + 1 / shape)
  rcontinuous(
    n = prep$ndraws, dist = "weibull",
    shape = shape, scale = scale,
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_frechet <- function(i, prep, ntrys = 5, ...) {
  nu <- get_dpar(prep, "nu", i = i)
  scale <- get_dpar(prep, "mu", i = i) / gamma(1 - 1 / nu)
  rcontinuous(
    n = prep$ndraws, dist = "frechet",
    scale = scale, shape = nu,
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_gen_extreme_value <- function(i, prep, ntrys = 5, ...) {
  rcontinuous(
    n = prep$ndraws, dist = "gen_extreme_value",
    sigma = get_dpar(prep, "sigma", i = i),
    xi = get_dpar(prep, "xi", i = i),
    mu = get_dpar(prep, "mu", i = i),
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_inverse.gaussian <- function(i, prep, ntrys = 5, ...) {
  rcontinuous(
    n = prep$ndraws, dist = "inv_gaussian",
    mu = get_dpar(prep, "mu", i = i),
    shape = get_dpar(prep, "shape", i = i),
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_exgaussian <- function(i, prep, ntrys = 5, ...) {
  rcontinuous(
    n = prep$ndraws, dist = "exgaussian",
    mu = get_dpar(prep, "mu", i = i),
    sigma = get_dpar(prep, "sigma", i = i),
    beta = get_dpar(prep, "beta", i = i),
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_wiener <- function(i, prep, negative_rt = FALSE, ntrys = 5,
                                     ...) {
  out <- rcontinuous(
    n = 1, dist = "wiener",
    delta = get_dpar(prep, "mu", i = i),
    alpha = get_dpar(prep, "bs", i = i),
    tau = get_dpar(prep, "ndt", i = i),
    beta = get_dpar(prep, "bias", i = i),
    types = if (negative_rt) c("q", "resp") else "q",
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
  if (negative_rt) {
    # code lower bound responses as negative RTs
    out <- out[["q"]] * ifelse(out[["resp"]], 1, -1)
  }
  out
}

posterior_predict_beta <- function(i, prep, ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i = i)
  phi <- get_dpar(prep, "phi", i = i)
  rcontinuous(
    n = prep$ndraws, dist = "beta",
    shape1 = mu * phi, shape2 = (1 - mu) * phi,
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_von_mises <- function(i, prep, ntrys = 5, ...) {
  rcontinuous(
    n = prep$ndraws, dist = "von_mises",
    mu = get_dpar(prep, "mu", i = i),
    kappa = get_dpar(prep, "kappa", i = i),
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_asym_laplace <- function(i, prep, ntrys = 5, ...) {
  rcontinuous(
    n = prep$ndraws, dist = "asym_laplace",
    mu = get_dpar(prep, "mu", i = i),
    sigma = get_dpar(prep, "sigma", i = i),
    quantile = get_dpar(prep, "quantile", i = i),
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_zero_inflated_asym_laplace <- function(i, prep, ntrys = 5,
                                                         ...) {
  zi <- get_dpar(prep, "zi", i = i)
  tmp <- runif(prep$ndraws, 0, 1)
  ifelse(
    tmp < zi, 0,
    rcontinuous(
      n = prep$ndraws, dist = "asym_laplace",
      mu = get_dpar(prep, "mu", i = i),
      sigma = get_dpar(prep, "sigma", i = i),
      quantile = get_dpar(prep, "quantile", i = i),
      lb = prep$data$lb[i], ub = prep$data$ub[i],
      ntrys = ntrys
    )
  )
}

posterior_predict_cox <- function(i, prep, ...) {
  stop2("Cannot sample from the posterior predictive ",
        "distribution for family 'cox'.")
}

posterior_predict_hurdle_poisson <- function(i, prep, ...) {
  # theta is the bernoulli hurdle parameter
  theta <- get_dpar(prep, "hu", i = i)
  lambda <- get_dpar(prep, "mu", i = i)
  ndraws <- prep$ndraws
  # compare with theta to incorporate the hurdle process
  hu <- runif(ndraws, 0, 1)
  # sample from a truncated poisson distribution
  # by adjusting lambda and adding 1
  t = -log(1 - runif(ndraws) * (1 - exp(-lambda)))
  ifelse(hu < theta, 0, rpois(ndraws, lambda = lambda - t) + 1)
}

posterior_predict_hurdle_negbinomial <- function(i, prep, ...) {
  # theta is the bernoulli hurdle parameter
  theta <- get_dpar(prep, "hu", i = i)
  mu <- get_dpar(prep, "mu", i = i)
  ndraws <- prep$ndraws
  # compare with theta to incorporate the hurdle process
  hu <- runif(ndraws, 0, 1)
  # sample from an approximate(!) truncated negbinomial distribution
  # by adjusting mu and adding 1
  t = -log(1 - runif(ndraws) * (1 - exp(-mu)))
  shape <- get_dpar(prep, "shape", i = i)
  ifelse(hu < theta, 0, rnbinom(ndraws, mu = mu - t, size = shape) + 1)
}

posterior_predict_hurdle_gamma <- function(i, prep, ...) {
  # theta is the bernoulli hurdle parameter
  theta <- get_dpar(prep, "hu", i = i)
  shape <- get_dpar(prep, "shape", i = i)
  scale <- get_dpar(prep, "mu", i = i) / shape
  ndraws <- prep$ndraws
  # compare with theta to incorporate the hurdle process
  hu <- runif(ndraws, 0, 1)
  ifelse(hu < theta, 0, rgamma(ndraws, shape = shape, scale = scale))
}

posterior_predict_hurdle_lognormal <- function(i, prep, ...) {
  # theta is the bernoulli hurdle parameter
  theta <- get_dpar(prep, "hu", i = i)
  mu <- get_dpar(prep, "mu", i = i)
  sigma <- get_dpar(prep, "sigma", i = i)
  ndraws <- prep$ndraws
  # compare with theta to incorporate the hurdle process
  hu <- runif(ndraws, 0, 1)
  ifelse(hu < theta, 0, rlnorm(ndraws, meanlog = mu, sdlog = sigma))
}

posterior_predict_zero_inflated_beta <- function(i, prep, ...) {
  # theta is the bernoulli hurdle parameter
  theta <- get_dpar(prep, "zi", i = i)
  mu <- get_dpar(prep, "mu", i = i)
  phi <- get_dpar(prep, "phi", i = i)
  # compare with theta to incorporate the hurdle process
  hu <- runif(prep$ndraws, 0, 1)
  ifelse(
    hu < theta, 0,
    rbeta(prep$ndraws, shape1 = mu * phi, shape2 = (1 - mu) * phi)
  )
}

posterior_predict_zero_one_inflated_beta <- function(i, prep, ...) {
  zoi <- get_dpar(prep, "zoi", i)
  coi <- get_dpar(prep, "coi", i)
  mu <- get_dpar(prep, "mu", i = i)
  phi <- get_dpar(prep, "phi", i = i)
  hu <- runif(prep$ndraws, 0, 1)
  one_or_zero <- runif(prep$ndraws, 0, 1)
  ifelse(hu < zoi,
    ifelse(one_or_zero < coi, 1, 0),
    rbeta(prep$ndraws, shape1 = mu * phi, shape2 = (1 - mu) * phi)
  )
}

posterior_predict_zero_inflated_poisson <- function(i, prep, ...) {
  # theta is the bernoulli zero-inflation parameter
  theta <- get_dpar(prep, "zi", i = i)
  lambda <- get_dpar(prep, "mu", i = i)
  ndraws <- prep$ndraws
  # compare with theta to incorporate the zero-inflation process
  zi <- runif(ndraws, 0, 1)
  ifelse(zi < theta, 0, rpois(ndraws, lambda = lambda))
}

posterior_predict_zero_inflated_negbinomial <- function(i, prep, ...) {
  # theta is the bernoulli zero-inflation parameter
  theta <- get_dpar(prep, "zi", i = i)
  mu <- get_dpar(prep, "mu", i = i)
  shape <- get_dpar(prep, "shape", i = i)
  ndraws <- prep$ndraws
  # compare with theta to incorporate the zero-inflation process
  zi <- runif(ndraws, 0, 1)
  ifelse(zi < theta, 0, rnbinom(ndraws, mu = mu, size = shape))
}

posterior_predict_zero_inflated_binomial <- function(i, prep, ...) {
  # theta is the bernoulli zero-inflation parameter
  theta <- get_dpar(prep, "zi", i = i)
  trials <- prep$data$trials[i]
  prob <- get_dpar(prep, "mu", i = i)
  ndraws <- prep$ndraws
  # compare with theta to incorporate the zero-inflation process
  zi <- runif(ndraws, 0, 1)
  ifelse(zi < theta, 0, rbinom(ndraws, size = trials, prob = prob))
}

posterior_predict_zero_inflated_beta_binomial <- function(i, prep, ...) {
  # theta is the bernoulli zero-inflation parameter
  theta <- get_dpar(prep, "zi", i = i)
  trials <- prep$data$trials[i]
  mu <- get_dpar(prep, "mu", i = i)
  phi <- get_dpar(prep, "phi", i = i)
  ndraws <- prep$ndraws
  draws <- rbeta_binomial(ndraws, size = trials, mu = mu, phi = phi)
  # compare with theta to incorporate the zero-inflation process
  zi <- runif(ndraws, 0, 1)
  draws[zi < theta] <- 0
  draws
}

posterior_predict_categorical <- function(i, prep, ...) {
  eta <- get_Mu(prep, i = i)
  eta <- insert_refcat(eta, refcat = prep$refcat)
  p <- pcategorical(seq_len(prep$data$ncat), eta = eta)
  first_greater(p, target = runif(prep$ndraws, min = 0, max = 1))
}

posterior_predict_multinomial <- function(i, prep, ...) {
  eta <- get_Mu(prep, i = i)
  eta <- insert_refcat(eta, refcat = prep$refcat)
  p <- dcategorical(seq_len(prep$data$ncat), eta = eta)
  size <- prep$data$trials[i]
  rblapply(seq_rows(p), function(s) t(rmultinom(1, size, p[s, ])))
}

posterior_predict_dirichlet <- function(i, prep, ...) {
  eta <- get_Mu(prep, i = i)
  eta <- insert_refcat(eta, refcat = prep$refcat)
  phi <- get_dpar(prep, "phi", i = i)
  cats <- seq_len(prep$data$ncat)
  alpha <- dcategorical(cats, eta = eta) * phi
  rdirichlet(prep$ndraws, alpha = alpha)
}

posterior_predict_dirichlet2 <- function(i, prep, ...) {
  mu <- get_Mu(prep, i = i)
  rdirichlet(prep$ndraws, alpha = mu)
}

posterior_predict_logistic_normal <- function(i, prep, ...) {
  mu <- get_Mu(prep, i = i)
  Sigma <- get_Sigma(prep, i = i, cor_name = "lncor")
  .predict <- function(s) {
    rlogistic_normal(1, mu = mu[s, ], Sigma = Sigma[s, , ],
                     refcat = prep$refcat)
  }
  rblapply(seq_len(prep$ndraws), .predict)
}

posterior_predict_cumulative <- function(i, prep, ...) {
  posterior_predict_ordinal(i = i, prep = prep)
}

posterior_predict_sratio <- function(i, prep, ...) {
  posterior_predict_ordinal(i = i, prep = prep)
}

posterior_predict_cratio <- function(i, prep, ...) {
  posterior_predict_ordinal(i = i, prep = prep)
}

posterior_predict_acat <- function(i, prep, ...) {
  posterior_predict_ordinal(i = i, prep = prep)
}

posterior_predict_ordinal <- function(i, prep, ...) {
  thres <- subset_thres(prep, i)
  nthres <- NCOL(thres)
  p <- pordinal(
    seq_len(nthres + 1),
    eta = get_dpar(prep, "mu", i = i),
    disc = get_dpar(prep, "disc", i = i),
    thres = thres,
    family = prep$family$family,
    link = prep$family$link
  )
  first_greater(p, target = runif(prep$ndraws, min = 0, max = 1))
}

posterior_predict_custom <- function(i, prep, ...) {
  custom_family_method(prep$family, "posterior_predict")(i, prep, ...)
}

posterior_predict_mixture <- function(i, prep, ...) {
  families <- family_names(prep$family)
  theta <- get_theta(prep, i = i)
  smix <- sample_mixture_ids(theta)
  out <- rep(NA, prep$ndraws)
  for (j in seq_along(families)) {
    draw_ids <- which(smix == j)
    if (length(draw_ids)) {
      pp_fun <- paste0("posterior_predict_", families[j])
      pp_fun <- get(pp_fun, asNamespace("brms"))
      tmp_prep <- pseudo_prep_for_mixture(prep, j, draw_ids)
      out[draw_ids] <- pp_fun(i, tmp_prep, ...)
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
# @param ntrys number of trys in rejection sampling for truncated models
# @return vector of random values prep from the distribution
rcontinuous <- function(n, dist, ..., lb = NULL, ub = NULL, ntrys = 5) {
  args <- list(...)
  if (is.null(lb) && is.null(ub)) {
    # sample as usual
    rdist <- paste0("r", dist)
    out <- do_call(rdist, c(list(n), args))
  } else {
    # sample from truncated distribution
    pdist <- paste0("p", dist)
    qdist <- paste0("q", dist)
    if (!exists(pdist, mode = "function") || !exists(qdist, mode = "function")) {
      # use rejection sampling as CDF or quantile function are not available
      out <- rdiscrete(n, dist, ..., lb = lb, ub = ub, ntrys = ntrys)
    } else {
      if (is.null(lb)) lb <- -Inf
      if (is.null(ub)) ub <- Inf
      plb <- do_call(pdist, c(list(lb), args))
      pub <- do_call(pdist, c(list(ub), args))
      out <- runif(n, min = plb, max = pub)
      out <- do_call(qdist, c(list(out), args))
      # infinite values may be caused by numerical imprecision
      out[out %in% c(-Inf, Inf)] <- NA
    }
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
    out <- vector("list", ntrys)
    for (i in seq_along(out)) {
      # loop of the trys to prevent a mismatch between 'n'
      # and length of the parameter vectors passed as arguments
      out[[i]] <- as.vector(do_call(rdist, c(list(n), args)))
    }
    out <- do_call(cbind, out)
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
