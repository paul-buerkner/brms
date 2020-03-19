#' Samples from the Posterior Predictive Distribution
#' 
#' Compute posterior samples of the posterior predictive distribution. Can be
#' performed for the data used to fit the model (posterior predictive checks) or
#' for new data. By definition, these samples have higher variance than samples
#' of the means of the posterior predictive distribution computed by
#' \code{\link{posterior_predict.brmsfit}}. This is because the residual error
#' is incorporated in \code{posterior_predict}. However, the estimated means of
#' both methods averaged across samples should be very similar.
#' 
#' @inheritParams extract_draws
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
#' @param ... Further arguments passed to \code{\link{extract_draws}}
#'   that control several aspects of data validation and prediction.
#' 
#' @return An \code{array} of predicted response values. If \code{summary =
#'   FALSE}, the output is as an S x N matrix, where S is the number of
#'   posterior samples and N is the number of observations. In multivariate
#'   models, an additional dimension is added to the output which indexes along
#'   the different response variables.
#' 
#' @details \code{NA} values within factors in \code{newdata}, 
#'   are interpreted as if all dummy variables of this factor are 
#'   zero. This allows, for instance, to make predictions of the grand mean 
#'   when using sum coding.  
#' 
#'   For truncated discrete models only: In the absence of any general algorithm
#'   to sample from truncated discrete distributions, rejection sampling is
#'   applied in this special case. This means that values are sampled until a
#'   value lies within the defined truncation boundaries. In practice, this
#'   procedure may be rather slow (especially in \R). Thus, we try to do
#'   approximate rejection sampling by sampling each value \code{ntrys} times
#'   and then select a valid value. If all values are invalid, the closest
#'   boundary is used, instead. If there are more than a few of these
#'   pathological cases, a warning will occur suggesting to increase argument
#'   \code{ntrys}.
#' 
#' @examples 
#' \dontrun{
#' ## fit a model
#' fit <- brm(time | cens(censored) ~ age + sex + (1 + age || patient), 
#'            data = kidney, family = "exponential", inits = "0")
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
  nsamples = NULL, subset = NULL, sort = FALSE, ntrys = 5, ...
) {
  cl <- match.call()
  if ("re.form" %in% names(cl)) {
    re_formula <- re.form
  }
  contains_samples(object)
  object <- restructure(object)
  draws <- extract_draws(
    object, newdata = newdata, re_formula = re_formula, resp = resp, 
    nsamples = nsamples, subset = subset, check_response = FALSE, ...
  )
  posterior_predict(
    draws, transform = transform, sort = sort, ntrys = ntrys, 
    negative_rt = negative_rt, summary = FALSE
  )
}

#' @export
posterior_predict.mvbrmsdraws <- function(object, ...) {
  if (length(object$mvpars$rescor)) {
    object$mvpars$Mu <- get_Mu(object)
    object$mvpars$Sigma <- get_Sigma(object)
    out <- posterior_predict.brmsdraws(object, ...)
  } else {
    out <- lapply(object$resps, posterior_predict, ...)
    along <- ifelse(length(out) > 1L, 3, 2)
    out <- do_call(abind, c(out, along = along))
  }
  out
}

#' @export
posterior_predict.brmsdraws <- function(object, transform = NULL, sort = FALSE,
                                        summary = FALSE, robust = FALSE, 
                                        probs = c(0.025, 0.975), ...) {
  for (nlp in names(object$nlpars)) {
    object$nlpars[[nlp]] <- get_nlpar(object, nlpar = nlp)
  }
  for (dp in names(object$dpars)) {
    object$dpars[[dp]] <- get_dpar(object, dpar = dp)
  }
  pp_fun <- paste0("posterior_predict_", object$family$fun)
  pp_fun <- get(pp_fun, asNamespace("brms"))
  N <- choose_N(object)
  out <- lapply(seq_len(N), pp_fun, draws = object, ...)
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
  colnames(out) <- NULL
  if (use_int(object$family)) {
    out <- check_discrete_trunc_bounds(
      out, lb = object$data$lb, ub = object$data$ub
    )
  }
  out <- reorder_obs(out, object$old_order, sort = sort)
  # transform predicted response samples before summarizing them 
  if (!is.null(transform)) {
    # deprecated as of brms 2.12.3
    warning2("Argument 'transform' is deprecated ", 
             "and will be removed in the future.")
    out <- do_call(transform, list(out))
  }
  attr(out, "levels") <- object$cats
  summary <- as_one_logical(summary)
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

#' Samples from the Posterior Predictive Distribution
#' 
#' This method is an alias of \code{\link{posterior_predict.brmsfit}}
#' with additional arguments for obtaining summaries of the computed samples.
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
#'   predicted category probabilites. For all other families, the output is a N
#'   x E matrix where E = \code{2 + length(probs)} is the number of summary
#'   statistics: The \code{Estimate} column contains point estimates (either
#'   mean or median depending on argument \code{robust}), while the
#'   \code{Est.Error} column contains uncertainty estimates (either standard
#'   deviation or median absolute deviation depending on argument
#'   \code{robust}). The remaining columns starting with \code{Q} contain
#'   quantile estimates as specifed via argument \code{probs}.  
#'   
#' @seealso \code{\link{posterior_predict.brmsfit}} 
#'  
#' @examples 
#' \dontrun{
#' ## fit a model
#' fit <- brm(time | cens(censored) ~ age + sex + (1 + age || patient), 
#'            data = kidney, family = "exponential", inits = "0")
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
                            negative_rt = FALSE, nsamples = NULL, 
                            subset = NULL, sort = FALSE, ntrys = 5, 
                            summary = TRUE, robust = FALSE,
                            probs = c(0.025, 0.975), ...) {
  contains_samples(object)
  object <- restructure(object)
  draws <- extract_draws(
    object, newdata = newdata, re_formula = re_formula, resp = resp, 
    nsamples = nsamples, subset = subset, check_response = FALSE, ...
  )
  posterior_predict(
    draws, transform = transform, ntrys = ntrys, negative_rt = negative_rt, 
    sort = sort, summary = summary, robust = robust, probs = probs
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
  } else if (method %in% c("pp_expect", "fitted")) {
    method <- "pp_expect"
  } else if (method %in% c("predictive_error", "residuals")) {
    method <- "predictive_error"
  } else if (method %in% c("posterior_linpred")) {
    method <- "posterior_linpred"
  } else {
    stop2("Posterior predictive method '", method, "' it not supported.")
  }
  method
}

# ------------------- family specific posterior_predict methods ---------------------
# All posterior_predict_<family> functions have the same arguments structure
# @param i the column of draws to use that is the ith obervation 
#   in the initial data.frame 
# @param draws A named list returned by extract_draws containing 
#   all required data and samples
# @param ... ignored arguments
# @param A vector of length draws$nsamples containing samples
#   from the posterior predictive distribution
posterior_predict_gaussian <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "norm",
    mean = get_dpar(draws, "mu", i = i), 
    sd = get_dpar(draws, "sigma", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

posterior_predict_student <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "student_t", 
    df = get_dpar(draws, "nu", i = i), 
    mu = get_dpar(draws, "mu", i = i), 
    sigma = get_dpar(draws, "sigma", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

posterior_predict_lognormal <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "lnorm",
    meanlog = get_dpar(draws, "mu", i = i), 
    sdlog = get_dpar(draws, "sigma", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

posterior_predict_shifted_lognormal <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "shifted_lnorm",
    meanlog = get_dpar(draws, "mu", i = i), 
    sdlog = get_dpar(draws, "sigma", i = i),
    shift = get_dpar(draws, "ndt", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

posterior_predict_skew_normal <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "skew_normal",
    mu = get_dpar(draws, "mu", i = i),
    sigma = get_dpar(draws, "sigma", i = i),
    alpha = get_dpar(draws, "alpha", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

posterior_predict_gaussian_mv <- function(i, draws, ...) {
  Mu <- get_Mu(draws, i = i)
  Sigma <- get_Sigma(draws, i = i)
  .predict <- function(s) {
    rmulti_normal(1, mu = Mu[s, ], Sigma = Sigma[s, , ])
  }
  rblapply(seq_len(draws$nsamples), .predict)
}

posterior_predict_student_mv <- function(i, draws, ...) {
  nu <- get_dpar(draws, "nu", i = i)
  Mu <- get_Mu(draws, i = i)
  Sigma <- get_Sigma(draws, i = i)
  .predict <- function(s) {
    rmulti_student_t(1, df = nu[s], mu = Mu[s, ], Sigma = Sigma[s, , ])
  }
  rblapply(seq_len(draws$nsamples), .predict)
}

posterior_predict_gaussian_time <- function(i, draws, ...) {
  obs <- with(draws$ac, begin_tg[i]:end_tg[i])
  mu <- as.matrix(get_dpar(draws, "mu", i = obs))
  Sigma <- get_cov_matrix_ac(draws, obs)
  .predict <- function(s) {
    rmulti_normal(1, mu = mu[s, ], Sigma = Sigma[s, , ])
  }
  rblapply(seq_len(draws$nsamples), .predict)
}

posterior_predict_student_time <- function(i, draws, ...) {
  obs <- with(draws$ac, begin_tg[i]:end_tg[i])
  nu <- as.matrix(get_dpar(draws, "nu", i = obs))
  mu <- as.matrix(get_dpar(draws, "mu", i = obs))
  Sigma <- get_cov_matrix_ac(draws, obs)
  .predict <- function(s) {
    rmulti_student_t(1, df = nu[s, ], mu = mu[s, ], Sigma = Sigma[s, , ])
  }
  rblapply(seq_len(draws$nsamples), .predict)
}

posterior_predict_gaussian_lagsar <- function(i, draws, ...) {
  stopifnot(i == 1)
  .predict <- function(s) {
    M_new <- with(draws, diag(nobs) - ac$lagsar[s] * ac$Msar)
    mu <- as.numeric(solve(M_new) %*% mu[s, ])
    Sigma <- solve(crossprod(M_new)) * sigma[s]^2
    rmulti_normal(1, mu = mu, Sigma = Sigma)
  }
  mu <- get_dpar(draws, "mu")
  sigma <- get_dpar(draws, "sigma")
  rblapply(seq_len(draws$nsamples), .predict)
}

posterior_predict_student_lagsar <- function(i, draws, ...) {
  stopifnot(i == 1)
  .predict <- function(s) {
    M_new <- with(draws, diag(nobs) - ac$lagsar[s] * ac$Msar)
    mu <- as.numeric(solve(M_new) %*% mu[s, ])
    Sigma <- solve(crossprod(M_new)) * sigma[s]^2
    rmulti_student_t(1, df = nu[s], mu = mu, Sigma = Sigma)
  }
  mu <- get_dpar(draws, "mu")
  sigma <- get_dpar(draws, "sigma")
  nu <- get_dpar(draws, "nu")
  rblapply(seq_len(draws$nsamples), .predict)
}

posterior_predict_gaussian_errorsar <- function(i, draws, ...) {
  stopifnot(i == 1)
  .predict <- function(s) {
    M_new <- with(draws, diag(nobs) - ac$errorsar[s] * ac$Msar)
    Sigma <- solve(crossprod(M_new)) * sigma[s]^2
    rmulti_normal(1, mu = mu[s, ], Sigma = Sigma)
  }
  mu <- get_dpar(draws, "mu")
  sigma <- get_dpar(draws, "sigma")
  rblapply(seq_len(draws$nsamples), .predict)
}

posterior_predict_student_errorsar <- function(i, draws, ...) {
  stopifnot(i == 1)
  .predict <- function(s) {
    M_new <- with(draws, diag(nobs) - ac$errorsar[s] * ac$Msar)
    Sigma <- solve(crossprod(M_new)) * sigma[s]^2
    rmulti_student_t(1, df = nu[s], mu = mu[s, ], Sigma = Sigma)
  }
  mu <- get_dpar(draws, "mu")
  sigma <- get_dpar(draws, "sigma")
  nu <- get_dpar(draws, "nu")
  rblapply(seq_len(draws$nsamples), .predict)
}

posterior_predict_gaussian_fcor <- function(i, draws, ...) {
  stopifnot(i == 1)
  mu <- as.matrix(get_dpar(draws, "mu"))
  Sigma <- get_cov_matrix_ac(draws)
  .predict <- function(s) {
    rmulti_normal(1, mu = mu[s, ], Sigma = Sigma[s, , ])
  }
  rblapply(seq_len(draws$nsamples), .predict)
}

posterior_predict_student_fcor <- function(i, draws, ...) {
  stopifnot(i == 1)
  nu <- as.matrix(get_dpar(draws, "nu"))
  mu <- as.matrix(get_dpar(draws, "mu"))
  Sigma <- get_cov_matrix_ac(draws)
  .predict <- function(s) {
    rmulti_student_t(1, df = nu[s, ], mu = mu[s, ], Sigma = Sigma[s, , ])
  }
  rblapply(seq_len(draws$nsamples), .predict)
}

posterior_predict_binomial <- function(i, draws, ntrys = 5, ...) {
  rdiscrete(
    n = draws$nsamples, dist = "binom", 
    size = draws$data$trials[i], 
    prob = get_dpar(draws, "mu", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i], 
    ntrys = ntrys
  )
}

posterior_predict_bernoulli <- function(i, draws, ...) {
  mu <- get_dpar(draws, "mu", i = i)
  rbinom(length(mu), size = 1, prob = mu)
}

posterior_predict_poisson <- function(i, draws, ntrys = 5, ...) {
  rdiscrete(
    n = draws$nsamples, dist = "pois",
    lambda = get_dpar(draws, "mu", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_negbinomial <- function(i, draws, ntrys = 5, ...) {
  rdiscrete(
    n = draws$nsamples, dist = "nbinom",
    mu = get_dpar(draws, "mu", i = i), 
    size = get_dpar(draws, "shape", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_geometric <- function(i, draws, ntrys = 5, ...) {
  rdiscrete(
    n = draws$nsamples, dist = "nbinom",
    mu = get_dpar(draws, "mu", i = i), size = 1,
    lb = draws$data$lb[i], ub = draws$data$ub[i], 
    ntrys = ntrys
  )
}

posterior_predict_discrete_weibull <- function(i, draws, ntrys = 5, ...) {
  rdiscrete(
    n = draws$nsamples, dist = "discrete_weibull",
    mu = get_dpar(draws, "mu", i = i), 
    shape = get_dpar(draws, "shape", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_com_poisson <- function(i, draws, ntrys = 5, ...) {
  rdiscrete(
    n = draws$nsamples, dist = "com_poisson",
    mu = get_dpar(draws, "mu", i = i), 
    shape = get_dpar(draws, "shape", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i],
    ntrys = ntrys
  )
}

posterior_predict_exponential <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "exp",
    rate = 1 / get_dpar(draws, "mu", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

posterior_predict_gamma <- function(i, draws, ...) {
  shape <- get_dpar(draws, "shape", i = i)
  scale <- get_dpar(draws, "mu", i = i) / shape
  rcontinuous(
    n = draws$nsamples, dist = "gamma",
    shape = shape, scale = scale,
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

posterior_predict_weibull <- function(i, draws, ...) {
  shape <- get_dpar(draws, "shape", i = i)
  scale <- get_dpar(draws, "mu", i = i) / gamma(1 + 1 / shape) 
  rcontinuous(
    n = draws$nsamples, dist = "weibull",
    shape = shape, scale = scale,
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

posterior_predict_frechet <- function(i, draws, ...) {
  nu <- get_dpar(draws, "nu", i = i)
  scale <- get_dpar(draws, "mu", i = i) / gamma(1 - 1 / nu)
  rcontinuous(
    n = draws$nsamples, dist = "frechet",
    scale = scale, shape = nu,
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

posterior_predict_gen_extreme_value <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "gen_extreme_value", 
    sigma = get_dpar(draws, "sigma", i = i),
    xi = get_dpar(draws, "xi", i = i),
    mu = get_dpar(draws, "mu", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

posterior_predict_inverse.gaussian <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "inv_gaussian",
    mu = get_dpar(draws, "mu", i = i), 
    shape = get_dpar(draws, "shape", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

posterior_predict_exgaussian <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "exgaussian",
    mu = get_dpar(draws, "mu", i = i), 
    sigma = get_dpar(draws, "sigma", i = i),
    beta = get_dpar(draws, "beta", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

posterior_predict_wiener <- function(i, draws, negative_rt = FALSE, ...) {
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

posterior_predict_beta <- function(i, draws, ...) {
  mu <- get_dpar(draws, "mu", i = i)
  phi <- get_dpar(draws, "phi", i = i)
  rcontinuous(
    n = draws$nsamples, dist = "beta", 
    shape1 = mu * phi, shape2 = (1 - mu) * phi,
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

posterior_predict_von_mises <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "von_mises",
    mu = get_dpar(draws, "mu", i = i), 
    kappa = get_dpar(draws, "kappa", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

posterior_predict_asym_laplace <- function(i, draws, ...) {
  rcontinuous(
    n = draws$nsamples, dist = "asym_laplace",
    mu = get_dpar(draws, "mu", i = i), 
    sigma = get_dpar(draws, "sigma", i = i),
    quantile = get_dpar(draws, "quantile", i = i),
    lb = draws$data$lb[i], ub = draws$data$ub[i]
  )
}

posterior_predict_zero_inflated_asym_laplace <- function(i, draws, ...) {
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

posterior_predict_cox <- function(i, draws, ...) {
  stop2("Cannot sample from the posterior predictive ",
        "distribution for family 'cox'.")
}

posterior_predict_hurdle_poisson <- function(i, draws, ...) {
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

posterior_predict_hurdle_negbinomial <- function(i, draws, ...) {
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

posterior_predict_hurdle_gamma <- function(i, draws, ...) {
  # theta is the bernoulli hurdle parameter
  theta <- get_dpar(draws, "hu", i = i)
  shape <- get_dpar(draws, "shape", i = i)
  scale <- get_dpar(draws, "mu", i = i) / shape
  ndraws <- draws$nsamples
  # compare with theta to incorporate the hurdle process
  hu <- runif(ndraws, 0, 1)
  ifelse(hu < theta, 0, rgamma(ndraws, shape = shape, scale = scale))
}

posterior_predict_hurdle_lognormal <- function(i, draws, ...) {
  # theta is the bernoulli hurdle parameter
  theta <- get_dpar(draws, "hu", i = i)
  mu <- get_dpar(draws, "mu", i = i)
  sigma <- get_dpar(draws, "sigma", i = i)
  ndraws <- draws$nsamples
  # compare with theta to incorporate the hurdle process
  hu <- runif(ndraws, 0, 1)
  ifelse(hu < theta, 0, rlnorm(ndraws, meanlog = mu, sdlog = sigma))
}

posterior_predict_zero_inflated_beta <- function(i, draws, ...) {
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

posterior_predict_zero_one_inflated_beta <- function(i, draws, ...) {
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

posterior_predict_zero_inflated_poisson <- function(i, draws, ...) {
  # theta is the bernoulli zero-inflation parameter
  theta <- get_dpar(draws, "zi", i = i)
  lambda <- get_dpar(draws, "mu", i = i)
  ndraws <- draws$nsamples
  # compare with theta to incorporate the zero-inflation process
  zi <- runif(ndraws, 0, 1)
  ifelse(zi < theta, 0, rpois(ndraws, lambda = lambda))
}

posterior_predict_zero_inflated_negbinomial <- function(i, draws, ...) {
  # theta is the bernoulli zero-inflation parameter
  theta <- get_dpar(draws, "zi", i = i)
  mu <- get_dpar(draws, "mu", i = i)
  shape <- get_dpar(draws, "shape", i = i)
  ndraws <- draws$nsamples
  # compare with theta to incorporate the zero-inflation process
  zi <- runif(ndraws, 0, 1)
  ifelse(zi < theta, 0, rnbinom(ndraws, mu = mu, size = shape))
}

posterior_predict_zero_inflated_binomial <- function(i, draws, ...) {
  # theta is the bernoulii zero-inflation parameter
  theta <- get_dpar(draws, "zi", i = i)
  trials <- draws$data$trials[i]
  prob <- get_dpar(draws, "mu", i = i)
  ndraws <- draws$nsamples
  # compare with theta to incorporate the zero-inflation process
  zi <- runif(ndraws, 0, 1)
  ifelse(zi < theta, 0, rbinom(ndraws, size = trials, prob = prob))
}

posterior_predict_categorical <- function(i, draws, ...) {
  eta <- sapply(names(draws$dpars), get_dpar, draws = draws, i = i)
  eta <- insert_refcat(eta, family = draws$family)
  p <- pcategorical(seq_len(draws$data$ncat), eta = eta)
  first_greater(p, target = runif(draws$nsamples, min = 0, max = 1))
}

posterior_predict_multinomial <- function(i, draws, ...) {
  eta <- sapply(names(draws$dpars), get_dpar, draws = draws, i = i)
  eta <- insert_refcat(eta, family = draws$family)
  p <- dcategorical(seq_len(draws$data$ncat), eta = eta)
  size <- draws$data$trials[i]
  out <- lapply(seq_rows(p), function(s) t(rmultinom(1, size, p[s, ])))
  do_call(rbind, out)
}

posterior_predict_dirichlet <- function(i, draws, ...) {
  mu_dpars <- str_subset(names(draws$dpars), "^mu")
  eta <- sapply(mu_dpars, get_dpar, draws = draws, i = i)
  eta <- insert_refcat(eta, family = draws$family)
  phi <- get_dpar(draws, "phi", i = i)
  cats <- seq_len(draws$data$ncat)
  alpha <- dcategorical(cats, eta = eta) * phi
  rdirichlet(draws$nsamples, alpha = alpha)
}

posterior_predict_cumulative <- function(i, draws, ...) {
  posterior_predict_ordinal(i = i, draws = draws)
}

posterior_predict_sratio <- function(i, draws, ...) {
  posterior_predict_ordinal(i = i, draws = draws)
}

posterior_predict_cratio <- function(i, draws, ...) {
  posterior_predict_ordinal(i = i, draws = draws)
}

posterior_predict_acat <- function(i, draws, ...) {
  posterior_predict_ordinal(i = i, draws = draws)
}  

posterior_predict_ordinal <- function(i, draws, ...) {
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

posterior_predict_custom <- function(i, draws, ...) {
  pp_fun <- draws$family$predict
  if (!is.function(pp_fun)) {
    pp_fun <- paste0("posterior_predict_", draws$family$name)
    pp_fun <- get(pp_fun, draws$family$env)
  }
  pp_fun(i = i, draws = draws, ...)
}

posterior_predict_mixture <- function(i, draws, ...) {
  families <- family_names(draws$family)
  theta <- get_theta(draws, i = i)
  smix <- sample_mixture_ids(theta)
  out <- rep(NA, draws$nsamples)
  for (j in seq_along(families)) {
    sample_ids <- which(smix == j)
    if (length(sample_ids)) {
      pp_fun <- paste0("posterior_predict_", families[j])
      pp_fun <- get(pp_fun, asNamespace("brms"))
      tmp_draws <- pseudo_draws_for_mixture(draws, j, sample_ids)
      out[sample_ids] <- pp_fun(i, tmp_draws, ...)
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
