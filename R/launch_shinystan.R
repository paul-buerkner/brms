#' Interface to \pkg{shinystan}
#'
#' Provide an interface to \pkg{shinystan} for models fitted with \pkg{brms}
#'
#' @aliases launch_shinystan
#'
#' @param object A fitted model object typically of class \code{brmsfit}.
#' @param rstudio Only relevant for RStudio users.
#' The default (\code{rstudio=FALSE}) is to launch the app
#' in the default web browser rather than RStudio's pop-up Viewer.
#' Users can change the default to \code{TRUE}
#' by setting the global option \cr \code{options(shinystan.rstudio = TRUE)}.
#' @param ... Optional arguments to pass to \code{\link[shiny:runApp]{runApp}}
#'
#' @return An S4 shinystan object
#'
#' @examples
#' \dontrun{
#' fit <- brm(rating ~ treat + period + carry + (1|subject),
#'            data = inhaler, family = "gaussian")
#' launch_shinystan(fit)
#' }
#'
#' @seealso \code{\link[shinystan:launch_shinystan]{launch_shinystan}}
#'
#' @exportS3Method shinystan::launch_shinystan brmsfit
launch_shinystan.brmsfit <- function(
  object, rstudio = getOption("shinystan.rstudio"), ...
) {
  contains_draws(object)
  if (object$algorithm != "sampling") {
    return(shinystan::launch_shinystan(object$fit, rstudio = rstudio, ...))
  }
  inc_warmup <- isTRUE(object$fit@sim$n_save[1] > niterations(object))
  draws <- as.array(object, inc_warmup = inc_warmup)
  warmup <- if (inc_warmup) nwarmup(object) else 0
  sampler_params <- rstan::get_sampler_params(object$fit, inc_warmup = inc_warmup)
  control <- object$fit@stan_args[[1]]$control
  if (is.null(control)) {
    max_td <- 10
  } else {
    max_td <- control$max_treedepth
    if (is.null(max_td)) {
      max_td <- 10
    }
  }
  sso <- shinystan::as.shinystan(
    X = draws,
    model_name = object$fit@model_name,
    warmup = warmup,
    sampler_params = sampler_params,
    max_treedepth = max_td,
    algorithm = "NUTS"
  )
  shinystan::launch_shinystan(sso, rstudio = rstudio, ...)
}
