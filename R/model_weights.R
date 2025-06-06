#' Model Weighting Methods
#'
#' Compute model weights in various ways, for instance, via
#' stacking of posterior predictive distributions, Akaike weights,
#' or marginal likelihoods.
#'
#' @inheritParams loo.brmsfit
#' @param weights Name of the criterion to compute weights from. Should be one
#'   of \code{"loo"}, \code{"waic"}, \code{"kfold"}, \code{"stacking"} (current
#'   default), \code{"bma"}, or \code{"pseudobma"}. For the former three
#'   options, Akaike weights will be computed based on the information criterion
#'   values returned by the respective methods. For \code{"stacking"} and
#'   \code{"pseudobma"}, method \code{\link{loo_model_weights}} will be used to
#'   obtain weights. For \code{"bma"}, method \code{\link{post_prob}} will be
#'   used to compute Bayesian model averaging weights based on log marginal
#'   likelihood values (make sure to specify reasonable priors in this case).
#'   For some methods, \code{weights} may also be a numeric vector of
#'   pre-specified weights.
#'
#' @return A numeric vector of weights for the models.
#'
#' @examples
#' \dontrun{
#' # model with 'treat' as predictor
#' fit1 <- brm(rating ~ treat + period + carry, data = inhaler)
#' summary(fit1)
#'
#' # model without 'treat' as predictor
#' fit2 <- brm(rating ~ period + carry, data = inhaler)
#' summary(fit2)
#'
#' # obtain Akaike weights based on the WAIC
#' model_weights(fit1, fit2, weights = "waic")
#' }
#'
#' @export
model_weights.brmsfit <- function(x, ..., weights = "stacking",
                                  model_names = NULL) {
  weights <- validate_weights_method(weights)
  args <- split_dots(x, ..., model_names = model_names)
  models <- args$models
  args$models <- NULL
  model_names <- names(models)
  if (weights %in% c("loo", "waic", "kfold")) {
    # Akaike weights based on information criteria
    ics <- rep(NA, length(models))
    for (i in seq_along(ics)) {
      args$x <- models[[i]]
      args$model_names <- names(models)[i]
      ics[i] <- SW(do_call(weights, args))$estimates[3, 1]
    }
    ic_diffs <- ics - min(ics)
    out <- exp(-ic_diffs / 2)
  } else if (weights %in% c("stacking", "pseudobma")) {
    args <- c(unname(models), args)
    args$method <- weights
    out <- do_call("loo_model_weights", args)
  } else if (weights %in% "bma") {
    args <- c(unname(models), args)
    out <- do_call("post_prob", args)
  }
  out <- as.numeric(out)
  out <- out / sum(out)
  names(out) <- model_names
  out
}

#' @rdname model_weights.brmsfit
#' @export
model_weights <- function(x, ...) {
  UseMethod("model_weights")
}

# validate name of the applied weighting method
validate_weights_method <- function(method) {
  method <- as_one_character(method)
  method <- tolower(method)
  if (method == "loo2") {
    warning2("Weight method 'loo2' is deprecated. Use 'stacking' instead.")
    method <- "stacking"
  }
  if (method == "marglik") {
    warning2("Weight method 'marglik' is deprecated. Use 'bma' instead.")
    method <- "bma"
  }
  options <- c("loo", "waic", "kfold", "stacking", "pseudobma", "bma")
  match.arg(method, options)
}

#' Posterior predictive draws averaged across models
#'
#' Compute posterior predictive draws averaged across models.
#' Weighting can be done in various ways, for instance using
#' Akaike weights based on information criteria or
#' marginal likelihoods.
#'
#' @inheritParams model_weights.brmsfit
#' @param method Method used to obtain predictions to average over. Should be
#'   one of \code{"posterior_predict"} (default), \code{"posterior_epred"},
#'   \code{"posterior_linpred"} or \code{"predictive_error"}.
#' @param control Optional \code{list} of further arguments
#'   passed to the function specified in \code{weights}.
#' @param ndraws Total number of posterior draws to use.
#' @param nsamples Deprecated alias of \code{ndraws}.
#' @param seed A single numeric value passed to \code{\link{set.seed}}
#'   to make results reproducible.
#' @param summary Should summary statistics
#'   (i.e. means, sds, and 95\% intervals) be returned
#'  instead of the raw values? Default is \code{TRUE}.
#' @param robust If \code{FALSE} (the default) the mean is used as
#'  the measure of central tendency and the standard deviation as
#'  the measure of variability. If \code{TRUE}, the median and the
#'  median absolute deviation (MAD) are applied instead.
#'  Only used if \code{summary} is \code{TRUE}.
#' @param probs  The percentiles to be computed by the \code{quantile}
#'  function. Only used if \code{summary} is \code{TRUE}.
#'
#' @return Same as the output of the method specified
#'   in argument \code{method}.
#'
#' @details Weights are computed with the \code{\link{model_weights}} method.
#'
#' @seealso \code{\link{model_weights}}, \code{\link{posterior_average}}
#'
#' @examples
#' \dontrun{
#' # model with 'treat' as predictor
#' fit1 <- brm(rating ~ treat + period + carry, data = inhaler)
#' summary(fit1)
#'
#' # model without 'treat' as predictor
#' fit2 <- brm(rating ~ period + carry, data = inhaler)
#' summary(fit2)
#'
#' # compute model-averaged predicted values
#' (df <- unique(inhaler[, c("treat", "period", "carry")]))
#' pp_average(fit1, fit2, newdata = df)
#'
#' # compute model-averaged fitted values
#' pp_average(fit1, fit2, method = "fitted", newdata = df)
#' }
#'
#' @export
pp_average.brmsfit <- function(
    x, ..., weights = "stacking", method = "posterior_predict",
    ndraws = NULL, nsamples = NULL, summary = TRUE, probs = c(0.025, 0.975),
    robust = FALSE, model_names = NULL, control = list(), seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  method <- validate_pp_method(method)
  ndraws <- use_alias(ndraws, nsamples)
  if (any(c("draw_ids", "subset") %in% names(list(...)))) {
    stop2("Cannot use argument 'draw_ids' in pp_average.")
  }
  args <- split_dots(x, ..., model_names = model_names)
  args$summary <- FALSE
  models <- args$models
  args$models <- NULL
  if (!match_response(models)) {
    stop2("Can only average models predicting the same response.")
  }
  if (is.null(ndraws)) {
    ndraws <- ndraws(models[[1]])
  }
  ndraws <- as_one_integer(ndraws)
  weights <- validate_weights(weights, models, control)
  ndraws <- round_largest_remainder(weights * ndraws)
  names(weights) <- names(ndraws) <- names(models)
  out <- named_list(names(models))
  for (i in seq_along(out)) {
    if (ndraws[i] > 0) {
      args$object <- models[[i]]
      args$ndraws <- ndraws[i]
      out[[i]] <- do_call(method, args)
    }
  }
  out <- do_call(rbind, out)
  if (summary) {
    out <- posterior_summary(out, probs = probs, robust = robust)
  }
  attr(out, "weights") <- weights
  attr(out, "ndraws") <- ndraws
  out
}

#' @rdname pp_average.brmsfit
#' @export
pp_average <- function(x, ...) {
  UseMethod("pp_average")
}

# validate weights passed to model averaging functions
# see pp_average.brmsfit for more documentation
validate_weights <- function(weights, models, control = list()) {
  if (!is.numeric(weights)) {
    weight_args <- c(unname(models), control)
    weight_args$weights <- weights
    weights <- do_call(model_weights, weight_args)
  } else {
    if (length(weights) != length(models)) {
      stop2(
        "If numeric, 'weights' must have the same length ",
        "as the number of models."
      )
    }
    if (any(weights < 0)) {
      stop2("If numeric, 'weights' must be positive.")
    }
  }
  weights / sum(weights)
}

#' Posterior draws of parameters averaged across models
#'
#' Extract posterior draws of parameters averaged across models.
#' Weighting can be done in various ways, for instance using
#' Akaike weights based on information criteria or
#' marginal likelihoods.
#'
#' @inheritParams pp_average.brmsfit
#' @param variable Names of variables (parameters) for which to average across
#'   models. Only those variables can be averaged that appear in every model.
#'   Defaults to all overlapping variables.
#' @param pars Deprecated alias of \code{variable}.
#' @param missing An optional numeric value or a named list of numeric values
#'   to use if a model does not contain a variable for which posterior draws
#'   should be averaged. Defaults to \code{NULL}, in which case only those
#'   variables can be averaged that are present in all of the models.
#'
#' @return A \code{data.frame} of posterior draws.
#'
#' @details Weights are computed with the \code{\link{model_weights}} method.
#'
#' @seealso \code{\link{model_weights}}, \code{\link{pp_average}}
#'
#' @examples
#' \dontrun{
#' # model with 'treat' as predictor
#' fit1 <- brm(rating ~ treat + period + carry, data = inhaler)
#' summary(fit1)
#'
#' # model without 'treat' as predictor
#' fit2 <- brm(rating ~ period + carry, data = inhaler)
#' summary(fit2)
#'
#' # compute model-averaged posteriors of overlapping parameters
#' posterior_average(fit1, fit2, weights = "waic")
#' }
#'
#' @export
posterior_average.brmsfit <- function(
    x, ..., variable = NULL, pars = NULL, weights = "stacking", ndraws = NULL,
    nsamples = NULL, missing = NULL, model_names = NULL, control = list(),
    seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  variable <- use_alias(variable, pars)
  ndraws <- use_alias(ndraws, nsamples)
  models <- split_dots(x, ..., model_names = model_names, other = FALSE)
  vars_list <- lapply(models, variables)
  all_vars <- unique(unlist(vars_list))
  if (is.null(missing)) {
    common_vars <- lapply(vars_list, function(x) all_vars %in% x)
    common_vars <- all_vars[Reduce("&", common_vars)]
    if (is.null(variable)) {
      variable <- setdiff(common_vars, "lp__")
    }
    variable <- as.character(variable)
    inv_vars <- setdiff(variable, common_vars)
    if (length(inv_vars)) {
      inv_vars <- collapse_comma(inv_vars)
      stop2(
        "Parameters ", inv_vars, " cannot be found in all ",
        "of the models. Consider using argument 'missing'."
      )
    }
  } else {
    if (is.null(variable)) {
      variable <- setdiff(all_vars, "lp__")
    }
    variable <- as.character(variable)
    inv_vars <- setdiff(variable, all_vars)
    if (length(inv_vars)) {
      inv_vars <- collapse_comma(inv_vars)
      stop2("Parameters ", inv_vars, " cannot be found in any of the models.")
    }
    if (is.list(missing)) {
      all_miss_vars <- unique(ulapply(
        models, function(m) setdiff(variable, variables(m))
      ))
      inv_vars <- setdiff(all_miss_vars, names(missing))
      if (length(inv_vars)) {
        stop2(
          "Argument 'missing' has no value for parameters ",
          collapse_comma(inv_vars), "."
        )
      }
      missing <- lapply(missing, as_one_numeric, allow_na = TRUE)
    } else {
      missing <- as_one_numeric(missing, allow_na = TRUE)
      missing <- named_list(variable, missing)
    }
  }
  if (is.null(ndraws)) {
    ndraws <- ndraws(models[[1]])
  }
  ndraws <- as_one_integer(ndraws)
  weights <- validate_weights(weights, models, control)
  ndraws <- round_largest_remainder(weights * ndraws)
  names(weights) <- names(ndraws) <- names(models)
  out <- named_list(names(models))
  for (i in seq_along(out)) {
    if (ndraws[i] > 0) {
      draw <- sample(seq_len(ndraws(models[[i]])), ndraws[i])
      draw <- sort(draw)
      found_vars <- intersect(variable, variables(models[[i]]))
      if (length(found_vars)) {
        out[[i]] <- as.data.frame(
          models[[i]],
          variable = found_vars, draw = draw
        )
      } else {
        out[[i]] <- as.data.frame(matrix(
          numeric(0),
          nrow = ndraws[i], ncol = 0
        ))
      }
      if (!is.null(missing)) {
        miss_vars <- setdiff(variable, names(out[[i]]))
        if (length(miss_vars)) {
          out[[i]][miss_vars] <- missing[miss_vars]
        }
      }
    }
  }
  out <- do_call(rbind, out)
  rownames(out) <- NULL
  attr(out, "weights") <- weights
  attr(out, "ndraws") <- ndraws
  out
}

#' @rdname posterior_average.brmsfit
#' @export
posterior_average <- function(x, ...) {
  UseMethod("posterior_average")
}
