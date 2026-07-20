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
#' @param output The type of output to return. Can be \code{"random"},
#'   \code{"probability"}, \code{"pit"}, \code{"density"}, or
#'   \code{"quantile"}. Defaults to \code{"random"}. In case of continuous
#'   distributions, \code{"probability"} is equivalent to \code{"pit"}.
#'   Not all families support outputs other than \code{"random"} yet;
#'   requesting an unsupported combination raises an error.
#' @param q Custom quantile for computing probability, PIT, or density values.
#'   It defaults to NULL in which case \code{prep$data$Y[i]} is used.
#' @param p Custom probability for computing quantile values.
#' @param lower.tail logical for computing probability or quantile values. It
#'   defaults to TRUE in which case probabilities are P(X < x) otherwise P(X > x).
#' @param log.p logical for computing probability or quantile values. It
#'   defaults to FALSE, if TRUE probabilities p are given as log(p).
#' @param log logical for computing density values. It defaults to FALSE.
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
  output = "random", q = NULL, p = NULL, cores = NULL, lower.tail = TRUE,
  log.p = FALSE, log = FALSE, ...
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
    negative_rt = negative_rt, cores = cores, summary = FALSE, output = output,
    q = q, p = p, lower.tail = lower.tail, log.p = log.p, log = log
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
                                       cores = NULL, output = "random", ...) {
  output <- as_one_character(output)
  output <- rlang::arg_match(
    output, values = c("random", "probability", "pit", "density", "quantile")
  )
  validate_pp_output_support(object$family$fun, output)

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
  out <- plapply(seq_len(N), pp_fun, .cores = cores, prep = object,
    output = output, ...)
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
                            robust = FALSE, probs = c(0.025, 0.975),
                            output = "random", ...) {
  contains_draws(object)
  object <- restructure(object)
  prep <- prepare_predictions(
    object, newdata = newdata, re_formula = re_formula, resp = resp,
    ndraws = ndraws, draw_ids = draw_ids, check_response = FALSE, ...
  )
  posterior_predict(
    prep, transform = transform, ntrys = ntrys, negative_rt = negative_rt,
    sort = sort, cores = cores, summary = summary, robust = robust,
    probs = probs, output = output
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
# @param i index of the observation for which to compute pp values
# @param prep A named list returned by prepare_predictions containing
#   all required data and posterior draws
# @param ... ignored arguments
# @param A vector of length prep$ndraws containing draws
#   from the posterior predictive distribution
posterior_predict_gaussian <- function(i, prep, output = "random", 
                                       ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i = i)
  sigma <- get_dpar(prep, "sigma", i = i)
  sigma <- add_sigma_se(sigma, prep, i = i)

  predict_continuous_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "norm", mean = mu, sd = sigma, ...
  )
}

posterior_predict_student <- function(i, prep, output = "random", 
                                      ntrys = 5, ...) {
  nu <- get_dpar(prep, "nu", i = i)
  mu <- get_dpar(prep, "mu", i = i)
  sigma <- get_dpar(prep, "sigma", i = i)
  sigma <- add_sigma_se(sigma, prep, i = i)

  predict_continuous_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "student_t", df = nu, mu = mu, sigma = sigma, ...
  )
}

posterior_predict_lognormal <- function(i, prep, output = "random", 
                                        ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i = i)
  sigma <- get_dpar(prep, "sigma", i = i)

  predict_continuous_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "lnorm", meanlog = mu, sdlog = sigma, ...
  )
}

posterior_predict_shifted_lognormal <- function(i, prep, output = "random",
                                                ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i = i)
  sigma <- get_dpar(prep, "sigma", i = i)
  ndt <- get_dpar(prep, "ndt", i = i)
  
  predict_continuous_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "shifted_lnorm", meanlog = mu, sdlog = sigma, shift = ndt, ...
  )
}

posterior_predict_skew_normal <- function(i, prep, output = "random",
                                          ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i = i)
  sigma <- get_dpar(prep, "sigma", i = i)
  sigma <- add_sigma_se(sigma, prep, i = i)
  alpha <- get_dpar(prep, "alpha", i = i)

  predict_continuous_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "skew_normal", mu = mu, sigma = sigma, alpha = alpha, ...
  )
}

posterior_predict_gaussian_mv <- function(i, prep, output = "random", ...) {
  validate_pp_output_support("gaussian_mv", output)
  Mu <- get_Mu(prep, i = i)
  Sigma <- get_Sigma(prep, i = i)
  .predict <- function(s) {
    rmulti_normal(1, mu = Mu[s, ], Sigma = Sigma[s, , ])
  }
  rblapply(seq_len(prep$ndraws), .predict)
}

posterior_predict_student_mv <- function(i, prep, output = "random", ...) {
  validate_pp_output_support("student_mv", output)
  nu <- get_dpar(prep, "nu", i = i)
  Mu <- get_Mu(prep, i = i)
  Sigma <- get_Sigma(prep, i = i)
  .predict <- function(s) {
    rmulti_student_t(1, df = nu[s], mu = Mu[s, ], Sigma = Sigma[s, , ])
  }
  rblapply(seq_len(prep$ndraws), .predict)
}

posterior_predict_gaussian_time <- function(i, prep, output = "random", ...) {
  validate_pp_output_support("gaussian_time", output)
  obs <- with(prep$ac, begin_tg[i]:end_tg[i])
  Jtime <- prep$ac$Jtime_tg[i, ]
  mu <- as.matrix(get_dpar(prep, "mu", i = obs))
  Sigma <- get_cov_matrix_ac(prep, obs, Jtime = Jtime)
  .predict <- function(s) {
    rmulti_normal(1, mu = mu[s, ], Sigma = Sigma[s, , ])
  }
  rblapply(seq_len(prep$ndraws), .predict)
}

posterior_predict_student_time <- function(i, prep, output = "random", ...) {
  validate_pp_output_support("student_time", output)
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

posterior_predict_gaussian_lagsar <- function(i, prep, output = "random", ...) {
  validate_pp_output_support("gaussian_lagsar", output)
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

posterior_predict_student_lagsar <- function(i, prep, output = "random", ...) {
  validate_pp_output_support("student_lagsar", output)
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

posterior_predict_gaussian_errorsar <- function(i, prep, output = "random",
                                                ...) {
  validate_pp_output_support("gaussian_errorsar", output)
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

posterior_predict_student_errorsar <- function(i, prep, output = "random",
                                               ...) {
  validate_pp_output_support("student_errorsar", output)
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

posterior_predict_gaussian_fcor <- function(i, prep, output = "random", ...) {
  validate_pp_output_support("gaussian_fcor", output)
  stopifnot(i == 1)
  mu <- as.matrix(get_dpar(prep, "mu"))
  Sigma <- get_cov_matrix_ac(prep)
  .predict <- function(s) {
    rmulti_normal(1, mu = mu[s, ], Sigma = Sigma[s, , ])
  }
  rblapply(seq_len(prep$ndraws), .predict)
}

posterior_predict_student_fcor <- function(i, prep, output = "random", ...) {
  validate_pp_output_support("student_fcor", output)
  stopifnot(i == 1)
  nu <- as.matrix(get_dpar(prep, "nu"))
  mu <- as.matrix(get_dpar(prep, "mu"))
  Sigma <- get_cov_matrix_ac(prep)
  .predict <- function(s) {
    rmulti_student_t(1, df = nu[s, ], mu = mu[s, ], Sigma = Sigma[s, , ])
  }
  rblapply(seq_len(prep$ndraws), .predict)
}

posterior_predict_binomial <- function(i, prep, output = "random", 
                                       ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i = i)
  size <- prep$data$trials[i]

  predict_discrete_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "binom", prob = mu, size = size, ...
  )
}

posterior_predict_beta_binomial <- function(i, prep, output = "random",
                                            ntrys = 5, ...) {
  size <- prep$data$trials[i]
  mu <- get_dpar(prep, "mu", i = i)
  phi <- get_dpar(prep, "phi", i = i)

  predict_discrete_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "beta_binomial", size = size, mu = mu, phi = phi, ...
  )
}

posterior_predict_bernoulli <- function(i, prep, output = "random", 
                                        ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i = i)

  predict_discrete_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "binom", prob = mu, size = 1, ...
  )
}

posterior_predict_poisson <- function(i, prep, output = "random", 
                                      ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i)
  mu <- multiply_dpar_rate_denom(mu, prep, i = i)

  predict_discrete_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "pois", lambda = mu, ...
  )
}

posterior_predict_negbinomial <- function(i, prep, output = "random",
                                          ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i)
  mu <- multiply_dpar_rate_denom(mu, prep, i = i)
  shape <- get_dpar(prep, "shape", i)
  shape <- multiply_dpar_rate_denom(shape, prep, i = i)

  predict_discrete_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "nbinom", mu = mu, size = shape, ...
  )
}

posterior_predict_negbinomial2 <- function(i, prep, output = "random",
                                           ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i)
  mu <- multiply_dpar_rate_denom(mu, prep, i = i)
  sigma <- get_dpar(prep, "sigma", i)
  shape <- multiply_dpar_rate_denom(1 / sigma, prep, i = i)

  predict_discrete_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "nbinom", mu = mu, size = shape, ...
  )
}

posterior_predict_geometric <- function(i, prep, output = "random",
                                        ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i)
  mu <- multiply_dpar_rate_denom(mu, prep, i = i)
  shape <- 1
  shape <- multiply_dpar_rate_denom(shape, prep, i = i)

  predict_discrete_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "nbinom", mu = mu, size = shape, ...
  )
}

posterior_predict_discrete_weibull <- function(i, prep, output = "random",
                                               ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i = i)
  shape <- get_dpar(prep, "shape", i = i)
  
  predict_discrete_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "discrete_weibull", mu = mu, shape = shape, ...
  )
}

posterior_predict_com_poisson <- function(i, prep, output = "random",
                                          ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i = i)
  shape <- get_dpar(prep, "shape", i = i)

  predict_discrete_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "com_poisson", mu = mu, shape = shape, ...
  )
}

posterior_predict_exponential <- function(i, prep, output = "random",
                                          ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i = i)
  rate <- 1 / mu
  
  predict_continuous_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "exp", rate = rate, ...
  )
}

posterior_predict_gamma <- function(i, prep, output = "random", 
                                    ntrys = 5, ...) {
  shape <- get_dpar(prep, "shape", i = i)
  scale <- get_dpar(prep, "mu", i = i) / shape

  predict_continuous_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "gamma", shape = shape, scale = scale, ...
  )
}

posterior_predict_weibull <- function(i, prep, output = "random", 
                                      ntrys = 5, ...) {
  shape <- get_dpar(prep, "shape", i = i)
  scale <- get_dpar(prep, "mu", i = i) / gamma(1 + 1 / shape)

  predict_continuous_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "weibull", shape = shape, scale = scale, ...
  )
}

posterior_predict_frechet <- function(i, prep, output = "random",
                                      ntrys = 5, ...) {
  nu <- get_dpar(prep, "nu", i = i)
  mu <- get_dpar(prep, "mu", i = i)
  scale <- mu / gamma(1 - 1 / nu)

  predict_continuous_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "frechet", scale = scale, shape = nu, ...
  )
}

posterior_predict_gen_extreme_value <- function(i, prep, output = "random",
                                                ntrys = 5, ...) {
  sigma <- get_dpar(prep, "sigma", i = i)
  xi <- get_dpar(prep, "xi", i = i)
  mu <- get_dpar(prep, "mu", i = i)

  predict_continuous_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "gen_extreme_value", sigma = sigma, xi = xi, mu = mu, ...
  )
}

posterior_predict_inverse.gaussian <- function(i, prep, output = "random",
                                               ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i = i)
  shape <- get_dpar(prep, "shape", i = i)

  predict_continuous_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "inv_gaussian", mu = mu, shape = shape, ...
  )
}

posterior_predict_exgaussian <- function(i, prep, output = "random", ntrys = 5,
                                         ...) {
  mu <- get_dpar(prep, "mu", i = i)
  sigma <- get_dpar(prep, "sigma", i = i)
  beta <- get_dpar(prep, "beta", i = i)

  predict_continuous_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "exgaussian", mu = mu, sigma = sigma, beta = beta, ...
  )
}

posterior_predict_wiener <- function(i, prep, output = "random",
                                     negative_rt = FALSE, ntrys = 5, ...) {
  validate_pp_output_support("wiener", output)
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

posterior_predict_beta <- function(i, prep, output = "random", ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i = i)
  phi <- get_dpar(prep, "phi", i = i)

  predict_continuous_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "beta", shape1 = mu * phi, shape2 = (1 - mu) * phi, ...
  )
}

posterior_predict_xbeta <- function(i, prep, output = "random", ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i = i)
  phi <- get_dpar(prep, "phi", i = i)
  kappa <- get_dpar(prep, "kappa", i = i)

  predict_continuous_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "xbeta", mu = mu, phi = phi, nu = kappa, ...
  )
}

posterior_predict_von_mises <- function(i, prep, output = "random", ntrys = 5,
                                        ...) {
  mu <- get_dpar(prep, "mu", i = i)
  kappa <- get_dpar(prep, "kappa", i = i)

  predict_continuous_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "von_mises", mu = mu, kappa = kappa, ...
  )
}

posterior_predict_asym_laplace <- function(i, prep, output = "random",
                                           ntrys = 5, ...) {
  mu <- get_dpar(prep, "mu", i = i)
  sigma <- get_dpar(prep, "sigma", i = i)
  quantile <- get_dpar(prep, "quantile", i = i)
  
  predict_continuous_helper(
    i = i, prep = prep, output = output, ntrys = ntrys,
    dist = "asym_laplace", mu = mu, sigma = sigma, quantile = quantile, ...
  )
}

posterior_predict_zero_inflated_asym_laplace <- function(i, prep,
                                                         output = "random",
                                                         ntrys = 5, ...) {
  zi <- get_dpar(prep, "zi", i = i)
  mu <- get_dpar(prep, "mu", i = i)
  sigma <- get_dpar(prep, "sigma", i = i)
  quantile <- get_dpar(prep, "quantile", i = i)

  if (output == "random") {
    out <- predict_continuous_helper(
      i = i, prep = prep, output = output, ntrys = ntrys,
      dist = "asym_laplace", mu = mu, sigma = sigma, quantile = quantile, ...
    )
    tmp <- runif(prep$ndraws, 0, 1)
    out <- ifelse(tmp < zi, 0, out)
  } else {
    out <- predict_continuous_helper(
      i = i, prep = prep, output = output, ntrys = ntrys,
      dist = "zero_inflated_asym_laplace",
      mu = mu, sigma = sigma, quantile = quantile, zi = zi, ...
    )
  }
  out
}

posterior_predict_cox <- function(i, prep, output = "random", ...) {
  stop2("Cannot sample from the posterior predictive ",
        "distribution for family 'cox'.")
}

posterior_predict_hurdle_poisson <- function(i, prep, output = "random",
                                             ntrys = 5, ...) {
  hu <- get_dpar(prep, "hu", i = i)
  lambda <- get_dpar(prep, "mu", i = i)
  lambda <- multiply_dpar_rate_denom(lambda, prep, i = i)

  if (output == "random") {
    if (!is.null(prep$data$lb[i]) || !is.null(prep$data$ub[i])) {
      warning2(
        "Truncated random sampling is not yet implemented for hurdle_poisson."
      )
    }
    qhurdle_poisson(runif(prep$ndraws), lambda = lambda, hu = hu)
  } else {
    predict_discrete_helper(
      i = i, prep = prep, output = output, ntrys = ntrys,
      dist = "hurdle_poisson", lambda = lambda, hu = hu, ...
    )
  }
}

posterior_predict_hurdle_negbinomial <- function(i, prep, output = "random",
                                                 ntrys = 5, ...) {
  hu <- get_dpar(prep, "hu", i = i)
  mu <- get_dpar(prep, "mu", i = i)
  mu <- multiply_dpar_rate_denom(mu, prep, i = i)
  shape <- get_dpar(prep, "shape", i = i)

  if (output == "random") {
    if (!is.null(prep$data$lb[i]) || !is.null(prep$data$ub[i])) {
      warning2(
        "Truncated random sampling is not yet implemented for hurdle_negbinomial."
      )
    }
    qhurdle_negbinomial(runif(prep$ndraws), mu = mu, shape = shape, hu = hu)
  } else {
    predict_discrete_helper(
      i = i, prep = prep, output = output, ntrys = ntrys,
      dist = "hurdle_negbinomial", mu = mu, shape = shape, hu = hu, ...
    )
  }
}

posterior_predict_hurdle_gamma <- function(i, prep, output = "random",
                                           ntrys = 5, ...) {
  hu <- get_dpar(prep, "hu", i = i)
  shape <- get_dpar(prep, "shape", i = i)
  scale <- get_dpar(prep, "mu", i = i) / shape

  if (output == "random") {
    if (!is.null(prep$data$lb[i]) || !is.null(prep$data$ub[i])) {
      warning2(
        "Truncated random sampling is not yet implemented for hurdle_gamma."
      )
    }
    qhurdle_gamma(runif(prep$ndraws), shape = shape, scale = scale, hu = hu)
  } else {
    predict_continuous_helper(
      i = i, prep = prep, output = output, ntrys = ntrys,
      dist = "hurdle_gamma", shape = shape, scale = scale, hu = hu, ...
    )
  }
}

posterior_predict_hurdle_lognormal <- function(i, prep, output = "random",
                                               ntrys = 5, ...) {
  hu <- get_dpar(prep, "hu", i = i)
  mu <- get_dpar(prep, "mu", i = i)
  sigma <- get_dpar(prep, "sigma", i = i)

  if (output == "random") {
    if (!is.null(prep$data$lb[i]) || !is.null(prep$data$ub[i])) {
      warning2(
        "Truncated random sampling is not yet implemented for hurdle_lognormal."
      )
    }
    qhurdle_lognormal(runif(prep$ndraws), mu = mu, sigma = sigma, hu = hu)
  } else {
    predict_continuous_helper(
      i = i, prep = prep, output = output, ntrys = ntrys,
      dist = "hurdle_lognormal", mu = mu, sigma = sigma, hu = hu, ...
    )
  }
}

posterior_predict_hurdle_cumulative <- function(i, prep, output = "random",
                                                ...) {
  predict_discrete_helper(
    i = i, prep = prep, output = output,
    dist = "hurdle_cumulative",
    eta = get_dpar(prep, "mu", i = i),
    disc = get_dpar(prep, "disc", i = i),
    hu = get_dpar(prep, "hu", i = i),
    thres = subset_thres(prep, i),
    link = prep$family$link,
    ...
  )
}

posterior_predict_zero_inflated_beta <- function(i, prep, output = "random",
                                                 ntrys = 5, ...) {
  zi <- get_dpar(prep, "zi", i = i)
  mu <- get_dpar(prep, "mu", i = i)
  phi <- get_dpar(prep, "phi", i = i)
  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi

  if (output == "random") {
    out <- predict_continuous_helper(
      i = i, prep = prep, output = output, ntrys = ntrys,
      dist = "beta", shape1 = shape1, shape2 = shape2, ...
    )
    tmp <- runif(prep$ndraws, 0, 1)
    out <- ifelse(tmp < zi, 0, out)
  } else {
    out <- predict_continuous_helper(
      i = i, prep = prep, output = output, ntrys = ntrys,
      dist = "zero_inflated_beta", shape1 = shape1, shape2 = shape2, zi = zi, ...
    )
  }
  out
}

posterior_predict_zero_one_inflated_beta <- function(i, prep, output = "random",
                                                     ntrys = 5, ...) {
  zoi <- get_dpar(prep, "zoi", i = i)
  coi <- get_dpar(prep, "coi", i = i)
  mu <- get_dpar(prep, "mu", i = i)
  phi <- get_dpar(prep, "phi", i = i)
  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi

  if (output == "random") {
    tmp <- runif(prep$ndraws, 0, 1)
    one_or_zero <- runif(prep$ndraws, 0, 1)
    out <- ifelse(
      tmp < zoi,
      ifelse(one_or_zero < coi, 1, 0),
      rbeta(prep$ndraws, shape1 = shape1, shape2 = shape2)
    )
  } else {
    out <- predict_continuous_helper(
      i = i, prep = prep, output = output, ntrys = ntrys,
      dist = "zero_one_inflated_beta",
      shape1 = shape1, shape2 = shape2, zoi = zoi, coi = coi, ...
    )
  }
  out
}

posterior_predict_zero_inflated_poisson <- function(i, prep, output = "random",
                                                    ntrys = 5, ...) {
  # zi is the bernoulli zero-inflation parameter
  zi <- get_dpar(prep, "zi", i = i)
  lambda <- get_dpar(prep, "mu", i = i)
  lambda <- multiply_dpar_rate_denom(lambda, prep, i = i)

  if (output == "random") {
    out <- predict_discrete_helper(
      i = i, prep = prep, output = output, ntrys = ntrys,
      dist = "pois", lambda = lambda, ...
    )
    tmp <- runif(prep$ndraws, 0, 1)
    out <- ifelse(tmp < zi, 0L, out)
  } else {
    out <- predict_discrete_helper(
      i = i, prep = prep, output = output, ntrys = ntrys,
      dist = "zero_inflated_poisson", lambda = lambda, zi = zi, ...
    )
  }
  out
}

posterior_predict_zero_inflated_negbinomial <- function(i, prep,
                                                        output = "random", ...) {
  zi <- get_dpar(prep, "zi", i = i)
  mu <- get_dpar(prep, "mu", i = i)
  shape <- get_dpar(prep, "shape", i = i)

  if (output == "random") {
    out <- predict_discrete_helper(i = i, prep = prep, output = output,
      dist = "nbinom", mu = mu, size = shape, ...)
    tmp <- runif(prep$ndraws, 0, 1)
    out <- ifelse(tmp < zi, 0L, out)
  } else {
    out <- predict_discrete_helper(i = i, prep = prep, output = output,
      dist = "zero_inflated_negbinomial", mu = mu, shape = shape, zi = zi, ...)
  }
  out
}

posterior_predict_zero_inflated_binomial <- function(i, prep, output = "random",
                                                     ntrys = 5, ...) {
  zi <- get_dpar(prep, "zi", i = i)
  trials <- prep$data$trials[i]
  prob <- get_dpar(prep, "mu", i = i)

  if (output == "random") {
    out <- predict_discrete_helper(
      i = i, prep = prep, output = output, ntrys = ntrys,
      dist = "binom", size = trials, prob = prob, ...
    )
    tmp <- runif(prep$ndraws, 0, 1)
    out <- ifelse(tmp < zi, 0L, out)
  } else {
    out <- predict_discrete_helper(
      i = i, prep = prep, output = output, ntrys = ntrys,
      dist = "zero_inflated_binomial", size = trials, prob = prob, zi = zi, ...
    )
  }
  out
}

posterior_predict_zero_inflated_beta_binomial <- function(i, prep,
                                                          output = "random", 
                                                          ntrys = 5, ...) {
  zi <- get_dpar(prep, "zi", i = i)
  trials <- prep$data$trials[i]
  mu <- get_dpar(prep, "mu", i = i)
  phi <- get_dpar(prep, "phi", i = i)

  if (output == "random") {
    out <- predict_discrete_helper(
      i = i, prep = prep, output = output, ntrys = ntrys,
      dist = "beta_binomial", size = trials, mu = mu, phi = phi, ...
    )
    tmp <- runif(prep$ndraws, 0, 1)
    out <- ifelse(tmp < zi, 0L, out)
  } else {
    out <- predict_discrete_helper(
      i = i, prep = prep, output = output, ntrys = ntrys,
      dist = "zero_inflated_beta_binomial",
      size = trials, mu = mu, phi = phi, zi = zi, ...
    )
  }
  out
}

posterior_predict_categorical <- function(i, prep, output = "random", ...) {
  eta <- get_Mu(prep, i = i)
  eta <- insert_refcat(eta, refcat = prep$refcat)
  predict_discrete_helper(
    i = i, prep = prep, output = output,
    dist = "categorical",
    eta = eta,
    ...
  )
}

posterior_predict_multinomial <- function(i, prep, output = "random", ...) {
  validate_pp_output_support("multinomial", output)
  eta <- get_Mu(prep, i = i)
  eta <- insert_refcat(eta, refcat = prep$refcat)
  p <- dcategorical(seq_len(prep$data$ncat), eta = eta)
  size <- prep$data$trials[i]
  rblapply(seq_rows(p), function(s) t(rmultinom(1, size, p[s, ])))
}

posterior_predict_dirichlet_multinomial <- function(i, prep, output = "random",
                                                    ...) {
  validate_pp_output_support("dirichlet_multinomial", output)
  eta <- get_Mu(prep, i = i)
  eta <- insert_refcat(eta, refcat = prep$refcat)
  phi <- get_dpar(prep, "phi", i = i)
  alpha <- dcategorical(seq_len(prep$data$ncat), eta = eta) * phi
  p <- rdirichlet(prep$ndraws, alpha = alpha)
  size <- prep$data$trials[i]
  rblapply(seq_rows(p), function(s) t(rmultinom(1, size, p[s, ])))
}

posterior_predict_dirichlet <- function(i, prep, output = "random", ...) {
  validate_pp_output_support("dirichlet", output)
  eta <- get_Mu(prep, i = i)
  eta <- insert_refcat(eta, refcat = prep$refcat)
  phi <- get_dpar(prep, "phi", i = i)
  cats <- seq_len(prep$data$ncat)
  alpha <- dcategorical(cats, eta = eta) * phi
  rdirichlet(prep$ndraws, alpha = alpha)
}

posterior_predict_dirichlet2 <- function(i, prep, output = "random", ...) {
  validate_pp_output_support("dirichlet2", output)
  mu <- get_Mu(prep, i = i)
  rdirichlet(prep$ndraws, alpha = mu)
}

posterior_predict_logistic_normal <- function(i, prep, output = "random", ...) {
  validate_pp_output_support("logistic_normal", output)
  mu <- get_Mu(prep, i = i)
  Sigma <- get_Sigma(prep, i = i, cor_name = "lncor")
  .predict <- function(s) {
    rlogistic_normal(1, mu = mu[s, ], Sigma = Sigma[s, , ],
                     refcat = prep$refcat)
  }
  rblapply(seq_len(prep$ndraws), .predict)
}

posterior_predict_cumulative <- function(i, prep, output = "random", ...) {
  posterior_predict_ordinal(i = i, prep = prep, output = output, ...)
}

posterior_predict_sratio <- function(i, prep, output = "random", ...) {
  posterior_predict_ordinal(i = i, prep = prep, output = output, ...)
}

posterior_predict_cratio <- function(i, prep, output = "random", ...) {
  posterior_predict_ordinal(i = i, prep = prep, output = output, ...)
}

posterior_predict_acat <- function(i, prep, output = "random", ...) {
  posterior_predict_ordinal(i = i, prep = prep, output = output, ...)
}

posterior_predict_ordinal <- function(i, prep, output = "random", ...) {
  thres <- subset_thres(prep, i)
  eta <- get_dpar(prep, "mu", i = i)
  disc <- get_dpar(prep, "disc", i = i)

  predict_discrete_helper(
    i = i, prep = prep, output = output,
    dist = "ordinal",
    eta = eta,
    disc = disc,
    thres = thres,
    family = prep$family$family,
    link = prep$family$link,
    ...
  )
}

posterior_predict_custom <- function(i, prep, output = "random", ...) {
  fun <- custom_family_method(prep$family, "posterior_predict")
  args <- names(formals(fun))
  if ("output" %in% args || "..." %in% args) {
    fun(i, prep, output = output, ...)
  } else {
    # older custom methods may not accept output yet
    if (!identical(output, "random")) {
      stop2(
        "Output '", output, "' is not yet implemented for this custom family. ",
        "Only output = 'random' is currently supported."
      )
    }
    fun(i, prep, ...)
  }
}

posterior_predict_mixture <- function(i, prep, output = "random", ...) {
  validate_pp_output_support("mixture", output)
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
      out[draw_ids] <- pp_fun(i, tmp_prep, output = output, ...)
    }
  }
  out
}

# ------------ predict helper-functions ----------------------

# families that currently only support output = "random"
# (multivariate, time-series, compositional, and diffusion models)
pp_output_random_only_families <- c(
  "gaussian_mv", "student_mv",
  "gaussian_time", "student_time",
  "gaussian_lagsar", "student_lagsar",
  "gaussian_errorsar", "student_errorsar",
  "gaussian_fcor", "student_fcor",
  "wiener",
  "multinomial", "dirichlet_multinomial",
  "dirichlet", "dirichlet2", "logistic_normal",
  "mixture"
)

# error if a non-random output is requested for an unsupported family
# @param family_fun name of the family method suffix (object$family$fun)
# @param output validated output type
# @noRd
validate_pp_output_support <- function(family_fun, output) {
  if (output == "random") {
    return(invisible(NULL))
  }
  family_fun <- as_one_character(family_fun)
  if (family_fun %in% pp_output_random_only_families) {
    stop2(
      "Output '", output, "' is not yet implemented for family '",
      family_fun, "'. Only output = 'random' is currently supported."
    )
  }
  invisible(NULL)
}

# random numbers from (possibly truncated) continuous distributions
# @param n number of random values to generate
# @param distribution name of a distribution for which the functions
#   p<distribution>, q<distribution>, and r<distribution> are available
# @param ... additional arguments passed to the distribution functions
# @param ntrys number of trys in rejection sampling for truncated models
# @return vector of random values prep from the distribution
rcontinuous <- function(n, distribution, ..., lb = NULL, ub = NULL, ntrys = 5) {
  args <- validate_distribution_args(distribution, ...)

  if (is.null(lb) && is.null(ub)) {
    # sample as usual
    rdist <- paste0("r", distribution)
    out <- do_call(rdist, c(list(n), args))
  } else {
    # sample from truncated distribution
    pdist <- paste0("p", distribution)
    qdist <- paste0("q", distribution)
    if (!exists(pdist, mode = "function") || !exists(qdist, mode = "function")) {
      # use rejection sampling as CDF or quantile function are not available
      out <- rdiscrete(n, distribution, ..., lb = lb, ub = ub, ntrys = ntrys)
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
# @param distribution name of a distribution for which the functions
#   p<distribution>, q<distribution>, and r<distribution> are available
# @param ... additional arguments passed to the distribution functions
# @param lb optional lower truncation bound
# @param ub optional upper truncation bound
# @param ntrys number of trys in rejection sampling for truncated models
# @return a vector of random values draws from the distribution
rdiscrete <- function(n, distribution, ..., lb = NULL, ub = NULL, ntrys = 5) {
  args <- validate_distribution_args(distribution, ...)
  rdist <- paste0("r", distribution)
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

# predict random numbers, CDF/PIT, density, or quantiles from
# continuous distributions
# @param i index of the observation for which to compute pp values
# @param prep A named list returned by prepare_predictions containing
#   all required data and posterior draws
# @param output "random", "probability", "pit", "density", or "quantile"
# @param distribution name of the distribution
# @param ntrys number of trys in rejection sampling for truncated models
# @param q optional custom quantile value; if NULL, the default is
#   prep$data$Y[i]
# @param p optional custom probability value for quantile output
# @param ... additional arguments passed to the distribution functions
# @return a vector of draws
predict_continuous_helper <- function(
  i, prep, output, distribution, ntrys, q = NULL, p = NULL, ...
) {
  lb <- prep$data$lb[i]
  ub <- prep$data$ub[i]
  if (output %in% c("probability", "pit", "density") && is.null(q)) q <- prep$data$Y[i]
  if (output == "quantile" && is.null(p)) {
    stop2("Argument 'p' must be specified when output = 'quantile'.")
  }

  switch(output,
    "random" = {
      dots <- list(...)
      dots[c("log.p", "lower.tail", "log", "p")] <- NULL
      do.call(rcontinuous, c(list(n = prep$ndraws, distribution = distribution,
                             lb = lb, ub = ub, ntrys = ntrys), dots))
    },
    # empty "probability" = , is a "fall-through", it means if the value is
    # "probability", do nothing and execute the next case's code block instead.
    "probability" = ,
    "pit" = {
      compute_cdf(q = q, distribution = distribution, lb = lb, ub = ub,
                  randomized = FALSE, ...)
    },
    "density" = {
      compute_density(q = q, distribution = distribution, lb = lb, ub = ub, ...)
    },
    "quantile" = {
      compute_quantile(p = p, distribution = distribution, lb = lb, ub = ub, ...)
    }
  )
}

# predict random numbers, CDF/PIT, density, or quantiles from discrete distributions
# @param i index of the observation for which to compute pp values
# @param prep A named list returned by prepare_predictions containing
#   all required data and posterior draws
# @param output "random", "probability", "pit", "density", or "quantile"
# @param distribution name of the distribution
# @param ntrys number of trys in rejection sampling for truncated models
# @param q optional custom quantile value; if NULL, the default is
#   prep$data$Y[i]
# @param p optional custom probability value for quantile output
# @param ... additional arguments passed to the distribution functions
# @return a vector of draws
predict_discrete_helper <- function(i, prep, output, distribution, ntrys = NULL,
                                    q = NULL, p = NULL, ...) {
  lb <- prep$data$lb[i]
  ub <- prep$data$ub[i]
  if (output %in% c("probability", "pit", "density") && is.null(q)) q <- prep$data$Y[i]
  if (output == "quantile" && is.null(p)) {
    stop2("Argument 'p' must be specified when output = 'quantile'.")
  }

  switch(output,
    "random" = {
      dots <- list(...)
      dots[c("log.p", "lower.tail", "log", "p")] <- NULL
      do.call(rdiscrete, c(list(n = prep$ndraws, distribution = distribution,
                                lb = lb, ub = ub, ntrys = ntrys), dots))
    },
    "probability" = {
      compute_cdf(q = q, distribution = distribution, lb = lb, ub = ub,
                  randomized = FALSE, ...)
    },
    "pit" = {
      compute_cdf(q = q, distribution = distribution, lb = lb, ub = ub,
                  randomized = TRUE, ...)
    },
    "density" = {
      compute_density(q = q, distribution = distribution, lb = lb, ub = ub, ...)
    },
    "quantile" = {
      compute_quantile(p = p, distribution = distribution, lb = lb, ub = ub, ...)
    }
  )
}

# compute cdf dependent on whether the distribution is truncated or not
# and whether to use the randomized PIT
# @param q quantile value(s) for which to compute the CDF
# @param distribution name of a distribution for which the functions
# @param lb optional lower truncation bound
# @param ub optional upper truncation bound
# @param randomized logical indicating whether to use the randomized PIT
# @param lower.tail logical; if TRUE (default) probabilities are P(X < x)
# otherwise, P(X > x)
# @param log.p logical; if TRUE probabilities p are given as log(p)
# @param ... additional arguments passed to the distribution functions
# @return a vector of probability values
# @noRd
compute_cdf <- function(q, distribution, lb, ub, randomized, lower.tail = TRUE,
                        log.p = FALSE, ...) {
  args <- validate_distribution_args(distribution, fun_prefix = "p", ...)
  pdist <- paste0("p", distribution)
  # prepare computation of (non-)truncated cdf
  F_internal <- function(q) {
    if (is.null(lb) && is.null(ub)) {
      do_call(pdist, c(list(q), args))
    } else {
      denom <- do_call(pdist, c(list(ub), args)) - do_call(pdist, c(list(lb), args))
      if (any(denom == 0)) stop("Division by zero")
      (do_call(pdist, c(list(q), args)) - do_call(pdist, c(list(lb), args))) / denom
    }
  }
  # randomized PIT specifically for discrete data (see, e.g.,
  # Czado, C., Gneiting, T., Held, L.: Predictive model
  # assessment for count data. Biometrics 65(4), 1254–1261 (2009).)
  # F(y-1) + V * [F(y) - F(y-1)] with V ~ Unif(0,1)
  if (isTRUE(randomized)) {
    v <- runif(length(q))
    probs <- F_internal(q - 1) + v * (F_internal(q) - F_internal(q - 1))
  } else if (isFALSE(randomized)) {
    probs <- F_internal(q)
  }
  if (isFALSE(lower.tail)) probs <- 1 - probs
  if (isTRUE(log.p)) probs <- log(probs)
  return(probs)
}

# compute density dependent on whether the distribution is truncated or not
# @param q quantile value(s) for which to compute the density
# @param distribution name of a distribution
# @param lb optional lower truncation bound
# @param ub optional upper truncation bound
# @param log logical; if TRUE densities are given as log(d)
# @param ... additional arguments passed to the distribution functions
# @return a vector of density values
# @noRd
compute_density <- function(q, distribution, lb, ub, log = FALSE, ...) {
  dargs <- validate_distribution_args(distribution, fun_prefix = "d", ...)
  pargs <- validate_distribution_args(distribution, fun_prefix = "p", ...)
  ddist <- paste0("d", distribution)
  pdist <- paste0("p", distribution)
  if (is.null(lb) && is.null(ub)) {
    return(do_call(ddist, c(list(q), dargs, log = log)))
  }
  if (is.null(lb)) {
    cdf_lb <- rep(0, length(q))
  } else {
    cdf_lb <- do_call(pdist, c(list(lb), pargs))
  }
  if (is.null(ub)) {
    cdf_ub <- rep(1, length(q))
  } else {
    cdf_ub <- do_call(pdist, c(list(ub), pargs))
  }
  denom <- cdf_ub - cdf_lb
  if (any(denom == 0)) stop("Division by zero")
  dens <- do_call(ddist, c(list(q), dargs, log = FALSE)) / denom
  if (!is.null(lb)) dens[q < lb] <- 0
  if (!is.null(ub)) dens[q > ub] <- 0
  if (isTRUE(log)) dens <- log(dens)
  dens
}

# compute quantile dependent on whether the distribution is truncated or not
# @param p probability value(s) for which to compute the quantile
# @param distribution name of a distribution
# @param lb optional lower truncation bound
# @param ub optional upper truncation bound
# @param lower.tail logical; if TRUE (default) probabilities are P(X < x)
# otherwise, P(X > x)
# @param log.p logical; if TRUE probabilities p are given as log(p)
# @param ... additional arguments passed to the distribution functions
# @return a vector of quantile values
# @noRd
compute_quantile <- function(p, distribution, lb, ub, lower.tail = TRUE,
                             log.p = FALSE, ...) {
  qargs <- validate_distribution_args(distribution, fun_prefix = "q", ...)
  pargs <- validate_distribution_args(distribution, fun_prefix = "p", ...)
  qdist <- paste0("q", distribution)
  pdist <- paste0("p", distribution)
  if (is.null(lb) && is.null(ub)) {
    return(do_call(
      qdist, c(list(p), qargs, lower.tail = lower.tail, log.p = log.p)
    ))
  }
  p <- validate_p_dist(p, lower.tail = lower.tail, log.p = log.p)
  if (is.null(lb)) {
    cdf_lb <- rep(0, length(p))
  } else {
    cdf_lb <- do_call(pdist, c(list(lb), pargs))
  }
  if (is.null(ub)) {
    cdf_ub <- rep(1, length(p))
  } else {
    cdf_ub <- do_call(pdist, c(list(ub), pargs))
  }
  denom <- cdf_ub - cdf_lb
  if (any(denom == 0)) stop("Division by zero")
  p_internal <- p * denom + cdf_lb
  do_call(qdist, c(list(p_internal), qargs))
}

# ensure that only arguments that are accepted by distribution functions are passed
validate_distribution_args <- function(distribution, fun_prefix = "p", ...) {
  args <- list(...)
  fun_prefix <- as_one_character(fun_prefix)
  rdist <- paste0(fun_prefix, distribution)
  rdist_fun <- match.fun(rdist)
  rdist_formals <- names(formals(rdist_fun))

  if (!is.null(rdist_formals) && !("..." %in% rdist_formals)) {
    args <- args[names(args) %in% rdist_formals]
  }

  args
}
