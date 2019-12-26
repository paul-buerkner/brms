#' Model Weighting Methods
#' 
#' Compute model weights in various ways, for instance via
#' stacking of predictive distributions, Akaike weights, or
#' marginal likelihoods.
#' 
#' @inheritParams loo.brmsfit
#' @param weights Name of the criterion to compute weights from. Should be one
#'   of \code{"loo"}, \code{"waic"}, \code{"kfold"}, \code{"stacking"} (current
#'   default), or \code{"bma"}, \code{"pseudobma"}, For the former three
#'   options, Akaike weights will be computed based on the information criterion
#'   values returned by the respective methods. For \code{"stacking"} and
#'   \code{"pseudobma"} method \code{\link{loo_model_weights}} will be used to
#'   obtain weights. For \code{"bma"}, method \code{\link{post_prob}} will be
#'   used to compute Bayesian model averaging weights based on log marginal
#'   likelihood values (make sure to specify reasonable priors in this case).
#'   Some some method, \code{weights} may also be be a numeric vector of
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

#' Posterior samples of parameters averaged across models
#' 
#' Extract posterior samples of parameters averaged across models.
#' Weighting can be done in various ways, for instance using
#' Akaike weights based on information criteria or 
#' marginal likelihoods.
#' 
#' @inheritParams pp_average.brmsfit
#' @param pars Names of parameters for which to average across models.
#'   Only those parameters can be averaged that appear in every model.
#'   Defaults to all overlapping parameters.
#' @param missing An optional numeric value or a named list of numeric values 
#'   to use if a model does not contain a parameter for which posterior samples 
#'   should be averaged. Defaults to \code{NULL}, in which case only those
#'   parameters can be averaged that are present in all of the models.
#' 
#' @return A \code{data.frame} of posterior samples. Samples are rows
#'   and parameters are columns.
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
  x, ..., pars = NULL, weights = "stacking", nsamples = NULL,
  missing = NULL, model_names = NULL, control = list(),
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed) 
  }
  models <- split_dots(x, ..., model_names = model_names, other = FALSE)
  pars_list <- lapply(models, parnames)
  all_pars <- unique(unlist(pars_list))
  if (is.null(missing)) {
    common_pars <- lapply(pars_list, function(x) all_pars %in% x)
    common_pars <- all_pars[Reduce("&", common_pars)]
    if (is.null(pars)) {
      pars <- setdiff(common_pars, "lp__")
    }
    pars <- as.character(pars)
    inv_pars <- setdiff(pars, common_pars)
    if (length(inv_pars)) {
      inv_pars <- collapse_comma(inv_pars)
      stop2(
        "Parameters ", inv_pars, " cannot be found in all ", 
        "of the models. Consider using argument 'missing'."
      )
    }
  } else {
    if (is.null(pars)) {
      pars <- setdiff(all_pars, "lp__")
    }
    pars <- as.character(pars)
    inv_pars <- setdiff(pars, all_pars)
    if (length(inv_pars)) {
      inv_pars <- collapse_comma(inv_pars)
      stop2("Parameters ", inv_pars, " cannot be found in any of the models.")
    }
    if (is.list(missing)) {
      all_miss_pars <- unique(ulapply(
        models, function(m) setdiff(pars, parnames(m))
      ))
      inv_pars <- setdiff(all_miss_pars, names(missing))
      if (length(inv_pars)) {
        stop2("Argument 'missing' has no value for parameters ",
              collapse_comma(inv_pars), ".")
      }
      missing <- lapply(missing, as_one_numeric, allow_na = TRUE)
    } else {
      missing <- as_one_numeric(missing, allow_na = TRUE)
      missing <- named_list(pars, missing)
    }
  }
  if (is.null(nsamples)) {
    nsamples <- nsamples(models[[1]])
  }
  weights <- validate_weights(weights, models, control)
  nsamples <- round_largest_remainder(weights * nsamples)
  names(weights) <- names(nsamples) <- names(models)
  out <- named_list(names(models))
  for (i in seq_along(out)) {
    if (nsamples[i] > 0) {
      subset <- sample(seq_len(nsamples(models[[i]])), nsamples[i])
      subset <- sort(subset)
      ps <- posterior_samples(
        models[[i]], pars = pars, 
        subset = subset, exact_match = TRUE
      )
      if (!is.null(ps)) {
        out[[i]] <- ps
      } else {
        out[[i]] <- as.data.frame(matrix(
          numeric(0), nrow = nsamples[i], ncol = 0
        ))
      }
      if (!is.null(missing)) {
        miss_pars <- setdiff(pars, names(out[[i]]))
        if (length(miss_pars)) {
          out[[i]][miss_pars] <- missing[miss_pars]
        }
      }
    }
  }
  out <- do_call(rbind, out)
  rownames(out) <- NULL
  attr(out, "weights") <- weights
  attr(out, "nsamples") <- nsamples
  out
}

#' @rdname posterior_average.brmsfit
#' @export
posterior_average <- function(x, ...) {
  UseMethod("posterior_average")
}

#' Posterior predictive samples averaged across models
#' 
#' Compute posterior predictive samples averaged across models.
#' Weighting can be done in various ways, for instance using
#' Akaike weights based on information criteria or 
#' marginal likelihoods.
#' 
#' @inheritParams model_weights.brmsfit
#' @param method Method used to obtain predictions to average over. Should be
#'   one of \code{"posterior_predict"} (default), \code{"pp_expect"}, or
#'   \code{"predictive_error"}.
#' @param control Optional \code{list} of further arguments 
#'   passed to the function specified in \code{weights}.
#' @param nsamples Total number of posterior samples to use.
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
  nsamples = NULL, summary = TRUE, probs = c(0.025, 0.975), robust = FALSE,
  model_names = NULL, control = list(), seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed) 
  }
  method <- validate_pp_method(method)
  if ("subset" %in% names(list(...))) {
    stop2("Cannot use argument 'subset' in pp_average.")
  }
  args <- split_dots(x, ..., model_names = model_names)
  args$summary <- FALSE
  models <- args$models
  args$models <- NULL
  if (!match_response(models)) {
    stop2("Can only average models predicting the same response.")
  }
  if (is.null(nsamples)) {
    nsamples <- nsamples(models[[1]])
  }
  weights <- validate_weights(weights, models, control)
  nsamples <- round_largest_remainder(weights * nsamples)
  names(weights) <- names(nsamples) <- names(models)
  out <- named_list(names(models))
  for (i in seq_along(out)) {
    if (nsamples[i] > 0) {
      args$object <- models[[i]]
      args$nsamples <- nsamples[i]
      out[[i]] <- do_call(method, args)
    }
  }
  out <- do_call(rbind, out)
  if (summary) {
    out <- posterior_summary(out, probs = probs, robust = robust) 
  }
  attr(out, "weights") <- weights
  attr(out, "nsamples") <- nsamples
  out
}

#' @rdname pp_average.brmsfit
#' @export
pp_average <- function(x, ...) {
  UseMethod("pp_average")
}
