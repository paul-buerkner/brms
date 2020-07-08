#' Posterior Samples of Predictive Errors
#' 
#' Compute posterior samples of predictive errors, that is, observed minus
#' predicted responses. Can be performed for the data used to fit the model
#' (posterior predictive checks) or for new data. 
#' 
#' @inheritParams posterior_predict.brmsfit
#' 
#' @return An S x N \code{array} of predictive error samples, where S is the
#'   number of posterior samples and N is the number of observations.
#' 
#' @examples 
#' \dontrun{
#' ## fit a model
#' fit <- brm(rating ~ treat + period + carry + (1|subject), 
#'            data = inhaler, cores = 2)
#' 
#' ## extract predictive errors
#' pe <- predictive_error(fit)
#' str(pe)
#' }
#' 
#' @aliases predictive_error
#' @method predictive_error brmsfit
#' @importFrom rstantools predictive_error
#' @export
#' @export predictive_error
predictive_error.brmsfit <- function(
  object, newdata = NULL, re_formula = NULL, re.form = NULL, 
  resp = NULL, nsamples = NULL, subset = NULL, sort = FALSE, ...
) {
  cl <- match.call()
  if ("re.form" %in% names(cl)) {
    re_formula <- re.form
  }
  .predictive_error(
    object, newdata = newdata, re_formula = re_formula,
    method = "posterior_predict", type = "ordinary", resp = resp, 
    nsamples = nsamples, subset = subset, sort = sort, ...
  )
}

#' Posterior Samples of Residuals/Predictive Errors
#' 
#' This method is an alias of \code{\link{predictive_error.brmsfit}}
#' with additional arguments for obtaining summaries of the computed samples.
#' 
#' @inheritParams predictive_error.brmsfit
#' @param method Method use to obtain predictions. Either
#'  \code{"pp_expect"} (the default) or \code{"posterior_predict"}.
#'  Using \code{"posterior_predict"} is recommended
#'  but \code{"pp_expect"} is the current default for 
#'  reasons of backwards compatibility.
#' @param type The type of the residuals, 
#'  either \code{"ordinary"} or \code{"pearson"}. 
#'  More information is provided under 'Details'.
#' @param summary Should summary statistics be returned
#'  instead of the raw values? Default is \code{TRUE}..
#' @param robust If \code{FALSE} (the default) the mean is used as 
#'  the measure of central tendency and the standard deviation as 
#'  the measure of variability. If \code{TRUE}, the median and the 
#'  median absolute deviation (MAD) are applied instead.
#'  Only used if \code{summary} is \code{TRUE}.
#' @param probs The percentiles to be computed by the \code{quantile} 
#'  function. Only used if \code{summary} is \code{TRUE}. 
#'  
#' @return An \code{array} of predictive error/residual samples. If
#'   \code{summary = FALSE} the output resembles those of
#'   \code{\link{predictive_error.brmsfit}}. If \code{summary = TRUE} the output
#'   is an N x E matrix, where N is the number of observations and E denotes
#'   the summary statistics computed from the samples.
#'  
#' @details Residuals of type \code{'ordinary'} are of the form \eqn{R = Y -
#'   Yrep}, where \eqn{Y} is the observed and \eqn{Yrep} is the predicted response.
#'   Residuals of type \code{pearson} are of the form \eqn{R = (Y - Yrep) /
#'   SD(Y)}, where \eqn{SD(Y)} is an estimation of the standard deviation of
#'   \eqn{Y}.
#'   
#' @examples 
#' \dontrun{
#' ## fit a model
#' fit <- brm(rating ~ treat + period + carry + (1|subject), 
#'            data = inhaler, cores = 2)
#' 
#' ## extract residuals/predictive errors
#' res <- residuals(fit)
#' head(res)
#' }
#'
#' @export
residuals.brmsfit <- function(object, newdata = NULL, re_formula = NULL, 
                              method = "pp_expect",
                              type = c("ordinary", "pearson"),
                              resp = NULL, nsamples = NULL,
                              subset = NULL, sort = FALSE, 
                              summary = TRUE, robust = FALSE, 
                              probs = c(0.025, 0.975), ...) {
  summary <- as_one_logical(summary)
  out <- .predictive_error(
    object, newdata = newdata, re_formula = re_formula,
    method = method, type = type, resp = resp, 
    nsamples = nsamples, subset = subset, sort = sort, ...
  )
  if (summary) {
    out <- posterior_summary(out, probs = probs, robust = robust)
  }
  out
}

# internal function doing the work for predictive_error.brmsfit
.predictive_error <- function(object, newdata, re_formula, method, type,  
                              resp, nsamples, subset, sort, ...) {
  contains_samples(object)
  object <- restructure(object)
  method <- validate_pp_method(method)
  type <- match.arg(type, c("ordinary", "pearson"))
  resp <- validate_resp(resp, object)
  family <- family(object, resp = resp)
  if (is_polytomous(family)) {
    stop2("Predictive errors are not defined for ordinal or categorical models.")
  }
  subset <- subset_samples(object, subset, nsamples)
  pred_args <- nlist(
    object, newdata, re_formula, resp, subset, 
    summary = FALSE, sort = TRUE, ...
  )
  yrep <- do_call(method, pred_args)
  y <- get_y(object, resp, newdata = newdata, warn = TRUE, ...)
  old_order <- attr(y, "old_order")
  if (length(dim(yrep)) == 3L) {
    # multivariate model
    y <- lapply(seq_cols(y), function(i) y[, i])
    y <- lapply(y, as_draws_matrix, dim = dim(yrep)[1:2])
    y <- abind(y, along = 3)
    dimnames(y)[[3]] <- dimnames(yrep)[[3]]
  } else {
    y <- as_draws_matrix(y, dim = dim(yrep))
  }
  out <- y - yrep
  remove(y, yrep)
  if (type == "pearson") {
    # deprecated as of brms 2.10.6
    warning2("Type 'pearson' is deprecated and will be removed in the future.")
    # get predicted standard deviation for each observation
    pred_args$summary <- TRUE
    pred <- do_call("predict", pred_args)
    if (length(dim(pred)) == 3L) {
      sd_pred <- array2list(pred[, 2, ])
      sd_pred <- lapply(sd_pred, as_draws_matrix, dim = dim(out)[1:2])
      sd_pred <- abind(sd_pred, along = 3)
    } else {
      sd_pred <- as_draws_matrix(pred[, 2], dim = dim(out))
    }
    out <- out / sd_pred
  }
  reorder_obs(out, old_order, sort = sort)
}
