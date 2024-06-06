#' Prior sensitivity: Create priorsense data
#'
#' The \code{create_priorsense_data.brmsfit} method can be used to
#' create the data structure needed by the \pkg{priorsense} package
#' for performing power-scaling sensitivity analysis. This method is
#' called automatically when performing powerscaling via
#' \code{\link[priorsense:powerscale]{powerscale}} or other related
#' functions, so you will rarely need to call it manually yourself.
#'
#' @param x A \code{brmsfit} object.
#' @param ... Currently unused.
#'
#' @return A \code{priorsense_data} object to be used in conjunction
#' with the \pkg{priorsense} package.
#'
#' @examples
#' \dontrun{
#' # fit a simple model with non-uniform priors
#' fit <- brm(count ~ zAge + zBase * Trt,
#'            data = epilepsy, family = poisson(),
#'            prior = prior(normal(0, 1), class = "b"))
#' summary(fit)
#'
#' # The following code requires the 'priorsense' package to be installed:
#' library(priorsense)
#'
#' # perform powerscaling of the prior
#' powerscale(fit, alpha = 1.5, component = "prior")
#'
#' # perform powerscaling sensitivity checks
#' powerscale_sensitivity(fit)
#'
#' # create powerscaling sensitivity plots
#' powerscale_plot_dens(fit)
#' }
#'
#' @exportS3Method priorsense::create_priorsense_data brmsfit
create_priorsense_data.brmsfit <- function(x, ...) {
  priorsense::create_priorsense_data.default(
    x = get_draws_ps(x),
    fit = x,
    log_prior = log_prior_ps(x),
    log_lik = log_lik_ps(x),
    log_prior_fn = log_prior_ps,
    log_lik_fn = log_lik_ps,
    log_ratio_fn = powerscale_log_ratio,
    ...
  )
}

log_lik_ps <- function(x) {
  log_lik <- log_lik(x)
  log_lik <- posterior::as_draws_array(log_lik)
  posterior::variables(log_lik) <- paste0("log_lik[", 1:nvariables(log_lik), "]")
  log_lik
}

log_prior_ps <- function(x, log_prior_name = "lprior") {
  posterior::subset_draws(
    posterior::as_draws_array(x),
    variable = log_prior_name
  )
}

get_draws_ps <- function(x, variable = NULL, regex = FALSE,
                         log_prior_name = "lprior") {
  excluded_variables <- c(log_prior_name, "lp__")
  draws <- posterior::as_draws_df(x, regex = regex)
  if (is.null(variable)) {
    # remove unnecessary variables
    variable <- posterior::variables(x)
    variable <- variable[!(variable %in% excluded_variables)]
    draws <- posterior::subset_draws(draws, variable = variable)
  }
  draws
}

powerscale_log_ratio <- function(draws, fit, alpha, component_fn) {
  component_draws <- component_fn(fit)
  component_draws <- rowsums_draws(component_draws)
  component_draws * (alpha - 1)
}

rowsums_draws <- function(x) {
  posterior::draws_array(
    sum = rowSums(
      posterior::as_draws_array(x),
      dims = 2
    ),
    .nchains = posterior::nchains(x)
  )
}
