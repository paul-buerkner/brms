#' Extract Diagnostic Quantities of \pkg{brms} Models
#'
#' Extract quantities that can be used to diagnose sampling behavior
#' of the algorithms applied by \pkg{Stan} at the back-end of \pkg{brms}.
#'
#' @name diagnostic-quantities
#' @aliases log_posterior nuts_params rhat neff_ratio
#'
#' @param object A \code{brmsfit} object.
#' @param pars An optional character vector of parameter names.
#'   For \code{nuts_params} these will be NUTS sampler parameter
#'   names rather than model parameters. If pars is omitted
#'   all parameters are included.
#' @param ... Arguments passed to individual methods.
#'
#' @return The exact form of the output depends on the method.
#'
#' @details For more details see
#'   \code{\link[bayesplot:bayesplot-extractors]{bayesplot-extractors}}.
#'
#' @examples
#' \dontrun{
#' fit <- brm(time ~ age * sex, data = kidney)
#'
#' lp <- log_posterior(fit)
#' head(lp)
#'
#' np <- nuts_params(fit)
#' str(np)
#' # extract the number of divergence transitions
#' sum(subset(np, Parameter == "divergent__")$Value)
#'
#' head(rhat(fit))
#' head(neff_ratio(fit))
#' }
NULL

#' @rdname diagnostic-quantities
#' @importFrom bayesplot log_posterior
#' @export log_posterior
#' @export
log_posterior.brmsfit <- function(object, ...) {
  contains_draws(object)
  bayesplot::log_posterior(object$fit, ...)
}

#' @rdname diagnostic-quantities
#' @importFrom bayesplot nuts_params
#' @export nuts_params
#' @export
nuts_params.brmsfit <- function(object, pars = NULL, ...) {
  contains_draws(object)
  bayesplot::nuts_params(object$fit, pars = pars, ...)
}

#' @rdname diagnostic-quantities
#' @importFrom bayesplot rhat
#' @export rhat
#' @export
rhat.brmsfit <- function(object, pars = NULL, ...) {
  contains_draws(object)
  # bayesplot uses outdated rhat code from rstan
  # bayesplot::rhat(object$fit, pars = pars, ...)
  draws <- as_draws_array(object, variable = pars, ...)
  tmp <- posterior::summarize_draws(draws, rhat = posterior::rhat)
  rhat <- tmp$rhat
  names(rhat) <- tmp$variable
  rhat
}

#' @rdname diagnostic-quantities
#' @importFrom bayesplot neff_ratio
#' @export neff_ratio
#' @export
neff_ratio.brmsfit <- function(object, pars = NULL, ...) {
  contains_draws(object)
  # bayesplot uses outdated ess code from rstan
  # bayesplot::neff_ratio(object$fit, pars = pars, ...)
  draws <- as_draws_array(object, variable = pars, ...)
  # currently uses ess_bulk as ess estimate for the central tendency
  tmp <- posterior::summarize_draws(draws, ess = posterior::ess_bulk)
  ess <- tmp$ess
  names(ess) <- tmp$variable
  ess / ndraws(draws)
}

#' Extract Control Parameters of the NUTS Sampler
#'
#' Extract control parameters of the NUTS sampler such as
#' \code{adapt_delta} or \code{max_treedepth}.
#'
#' @param x An \R object
#' @param pars Optional names of the control parameters to be returned.
#'  If \code{NULL} (the default) all control parameters are returned.
#'  See \code{\link[rstan:stan]{stan}} for more details.
#' @param ... Currently ignored.
#'
#' @return A named \code{list} with control parameter values.
#'
#' @export
control_params <- function(x, ...) {
  UseMethod("control_params")
}

#' @rdname control_params
#' @export
control_params.brmsfit <- function(x, pars = NULL, ...) {
  contains_draws(x)
  if (is_equal(x$backend, "cmdstanr")) {
    out <- attr(x$fit, "metadata")$metadata
  } else {
    out <- attr(x$fit@sim$samples[[1]], "args")$control
  }
  if (!is.null(pars)) {
    out <- out[pars]
  }
  out
}
