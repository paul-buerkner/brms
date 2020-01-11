#' Compute Weighted Expectations Using LOO
#' 
#' These functions are wrappers around the \code{\link[loo]{E_loo}} 
#' function of the \pkg{loo} package.
#'
#' @aliases loo_predict loo_linpred loo_predictive_interval
#' 
#' @param object An object of class \code{brmsfit}.
#' @param type The statistic to be computed on the results. 
#'   Can by either \code{"mean"} (default), \code{"var"}, or
#'   \code{"quantile"}.
#' @param probs A vector of quantiles to compute. 
#'   Only used if \code{type = quantile}.
#' @param prob For \code{loo_predictive_interval}, a scalar in \eqn{(0,1)}
#'   indicating the desired probability mass to include in the intervals. The
#'   default is \code{prob = 0.9} (\eqn{90}\% intervals).
#' @param psis_object An optional object returend by \code{\link[loo]{psis}}. 
#'   If \code{psis_object} is missing then \code{\link[loo]{psis}} is executed 
#'   internally, which may be time consuming for models fit to very large datasets. 
#' @param ... Optional arguments passed to the underlying methods that is 
#'   \code{\link[brms:log_lik.brmsfit]{log_lik}}, as well as
#'   \code{\link[brms:posterior_predict.brmsfit]{posterior_predict}} or
#'   \code{\link[brms:posterior_linpred.brmsfit]{posterior_linpred}}. 
#' @inheritParams posterior_predict.brmsfit
#'   
#' @return \code{loo_predict} and \code{loo_linpred} return a vector with one 
#'   element per observation. The only exception is if \code{type = "quantile"} 
#'   and \code{length(probs) >= 2}, in which case a separate vector for each 
#'   element of \code{probs} is computed and they are returned in a matrix with 
#'   \code{length(probs)} rows and one column per observation.
#'   
#'   \code{loo_predictive_interval} returns a matrix with one row per 
#'   observation and two columns. 
#'   \code{loo_predictive_interval(..., prob = p)} is equivalent to 
#'   \code{loo_predict(..., type = "quantile", probs = c(a, 1-a))} with 
#'   \code{a = (1 - p)/2}, except it transposes the result and adds informative 
#'   column names.
#'   
#' @examples
#' \dontrun{
#' ## data from help("lm")
#' ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
#' trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
#' d <- data.frame(
#'   weight = c(ctl, trt), 
#'   group = gl(2, 10, 20, labels = c("Ctl", "Trt"))
#' ) 
#' fit <- brm(weight ~ group, data = d)
#' loo_predictive_interval(fit, prob = 0.8)
#' 
#' ## optionally log-weights can be pre-computed and reused
#' psis <- loo::psis(-log_lik(fit), cores = 2)
#' loo_predictive_interval(fit, prob = 0.8, psis_object = psis)
#' loo_predict(fit, type = "var", psis_object = psis)
#' }
#' 
#' @method loo_predict brmsfit
#' @importFrom rstantools loo_predict
#' @export loo_predict
#' @export 
loo_predict.brmsfit <- function(object, type = c("mean", "var", "quantile"), 
                                probs = 0.5, psis_object = NULL, resp = NULL,
                                ...) {
  type <- match.arg(type)
  stopifnot_resp(object, resp)
  if (is.null(psis_object)) {
    message("Running PSIS to compute weights")
    psis_object <- compute_loo(object, criterion = "psis", resp = resp, ...)
  }
  preds <- posterior_predict(object, resp = resp, ...)
  loo::E_loo(preds, psis_object, type = type, probs = probs)$value
}

#' @rdname loo_predict.brmsfit
#' @method loo_linpred brmsfit
#' @importFrom rstantools loo_linpred
#' @export loo_linpred
#' @export 
loo_linpred.brmsfit <- function(object, type = c("mean", "var", "quantile"), 
                                probs = 0.5, psis_object = NULL, resp = NULL,
                                ...) {
  type <- match.arg(type)
  stopifnot_resp(object, resp)
  family <- family(object, resp = resp)
  if (is_ordinal(family) || is_categorical(family)) {
    stop2("Method 'loo_linpred' is not implemented ", 
          "for categorical or ordinal models")
  }
  if (is.null(psis_object)) {
    message("Running PSIS to compute weights")
    psis_object <- compute_loo(object, criterion = "psis", resp = resp, ...)
  }
  preds <- posterior_linpred(object, resp = resp, ...)
  loo::E_loo(preds, psis_object, type = type, probs = probs)$value
}

#' @rdname loo_predict.brmsfit
#' @method loo_predictive_interval brmsfit
#' @importFrom rstantools loo_predictive_interval
#' @export loo_predictive_interval
#' @export
loo_predictive_interval.brmsfit <- function(object, prob = 0.9,
                                            psis_object = NULL, ...) {
  if (length(prob) != 1L) {
    stop2("Argument 'prob' should be of length 1.")
  }
  alpha <- (1 - prob) / 2
  probs <- c(alpha, 1 - alpha)
  labs <- paste0(100 * probs, "%")
  intervals <- loo_predict(
    object, type = "quantile", probs = probs, 
    psis_object = psis_object, ...
  )
  rownames(intervals) <- labs
  t(intervals)
}

#' Compute a LOO-adjusted R-squared for regression models
#' 
#' @aliases loo_R2
#' 
#' @inheritParams posterior_predict.brmsfit
#' @param ... Further arguments passed to 
#'   \code{\link[brms:pp_expect.brmsfit]{pp_expect}} and
#'   \code{\link[brms:log_lik.brmsfit]{log_lik}},
#'   which are used in the computation of the R-squared values.
#' 
#' @return A real value per response variable indicating 
#' the LOO-adjusted R-squared.
#'  
#' @examples
#' \dontrun{
#' fit <- brm(mpg ~ wt + cyl, data = mtcars)
#' summary(fit)
#' loo_R2(fit)
#' 
#' # compute R2 with new data
#' nd <- data.frame(mpg = c(10, 20, 30), wt = c(4, 3, 2), cyl = c(8, 6, 4))
#' loo_R2(fit, newdata = nd)
#' }
#' 
#' @method loo_R2 brmsfit
#' @export
loo_R2.brmsfit <- function(object, resp = NULL, ...) {
  contains_samples(object)
  object <- restructure(object)
  resp <- validate_resp(resp, object)
  R2 <- get_criterion(object, "loo_R2")
  if (is.vector(R2)) {
    take <- names(R2) %in% paste0("R2", resp)
    R2 <- R2[take]
    return(R2)
  } 
  family <- family(object, resp = resp)
  if (conv_cats_dpars(family)) {
    stop2("'loo_R2' is not defined for unordered categorical models.")
  }
  if (is_ordinal(family)) {
    warning2(
      "Predictions are treated as continuous variables in ",
      "'loo_R2' which is likely invalid for ordinal families."
    )
  }
  # see http://discourse.mc-stan.org/t/stan-summary-r2-or-adjusted-r2/4308/4
  .loo_R2 <- function(y, ypred, ll) {
    r_eff <- r_eff_log_lik(ll, object)
    psis_object <- loo::psis(log_ratios = -ll, r_eff = r_eff)
    ypredloo <- loo::E_loo(ypred, psis_object, log_ratios = -ll)$value
    eloo <- ypredloo - y
    return(1 - var(eloo) / var(y))
  }
  args_y <- list(object, warn = TRUE, ...)
  args_ypred <- list(object, sort = TRUE, ...)
  R2 <- named_list(paste0("R2", resp))
  for (i in seq_along(R2)) {
    # assumes expectations of different responses to be independent
    args_ypred$resp <- args_y$resp <- resp[i]
    y <- do_call(get_y, args_y)
    ypred <- do.call(pp_expect, args_ypred)
    ll <- do_call(log_lik, args_ypred)
    if (is_ordinal(family(object, resp = resp[i]))) {
      ypred <- ordinal_probs_continuous(ypred)
    }
    R2[[i]] <- .loo_R2(y, ypred, ll)
  }
  R2 <- unlist(R2)
  names(R2) <- paste0("R2", resp)
  R2
}

# temporary generic until available in loo
#' @export
loo_R2 <- function(object, ...) {
  UseMethod("loo_R2")
}
