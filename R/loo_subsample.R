#' Efficient approximate leave-one-out cross-validation (LOO) using subsampling
#'
#' @aliases loo_subsample
#'
#' @inheritParams loo.brmsfit
#'
#' @details More details can be found on
#' \code{\link[loo:loo_subsample]{loo_subsample}}.
#'
#' @examples
#' \dontrun{
#' # model with population-level effects only
#' fit1 <- brm(rating ~ treat + period + carry,
#'   data = inhaler
#' )
#' (loo1 <- loo_subsample(fit1))
#'
#' # model with an additional varying intercept for subjects
#' fit2 <- brm(rating ~ treat + period + carry + (1 | subject),
#'   data = inhaler
#' )
#' (loo2 <- loo_subsample(fit2))
#'
#' # compare both models
#' loo_compare(loo1, loo2)
#' }
#'
#' @importFrom loo loo_subsample
#' @export loo_subsample
#' @export
loo_subsample.brmsfit <- function(x, ..., compare = TRUE, resp = NULL,
                                  model_names = NULL) {
  args <- split_dots(x, ..., model_names = model_names)
  c(args) <- nlist(
    criterion = "loo_subsample", compare, resp,
    add_point_estimate = TRUE
  )
  do_call(compute_loolist, args)
}

# compute 'loo_subsample' criterion using the 'loo' package
# @param model_name ignored but included to avoid being passed to '...'
.loo_subsample <- function(x, newdata, resp, model_name, ...) {
  loo_args <- prepare_loo_args(
    x,
    newdata = newdata, resp = resp,
    pointwise = TRUE, ...
  )
  do_call("loo_subsample", loo_args, pkg = "loo")
}

# methods required in loo_subsample
#' @importFrom loo .ndraws
#' @export
.ndraws.brmsprep <- function(x) {
  x$ndraws
}

#' @export
.ndraws.mvbrmsprep <- function(x) {
  x$ndraws
}

#' @importFrom loo .thin_draws
#' @export
.thin_draws.brmsprep <- function(draws, loo_approximation_draws) {
  # brmsprep objects are too complex to implement a post-hoc subsetting method
  if (length(loo_approximation_draws)) {
    stop2("'loo_approximation_draws' is not supported for brmsfit objects.")
  }
  draws
}

#' @export
.thin_draws.mvbrmsprep <- function(draws, loo_approximation_draws) {
  if (length(loo_approximation_draws)) {
    stop2("'loo_approximation_draws' is not supported for brmsfit objects.")
  }
  draws
}

#' @importFrom loo .compute_point_estimate
#' @export
.compute_point_estimate.brmsprep <- function(draws) {
  # point estimates are stored in the form of an attribute rather
  # than computed on the fly due to the complexity of brmsprep objects
  attr(draws, "point_estimate")
}

#' @export
.compute_point_estimate.mvbrmsprep <- function(draws) {
  attr(draws, "point_estimate")
}
