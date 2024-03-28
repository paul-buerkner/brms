#' Compute exact cross-validation for problematic observations
#'
#' Compute exact cross-validation for problematic observations for which
#' approximate leave-one-out cross-validation may return incorrect results.
#' Models for problematic observations can be run in parallel using the
#' \pkg{future} package.
#'
#' @inheritParams predict.brmsfit
#' @param x An \R object of class \code{brmsfit} or \code{loo} depending
#'   on the method.
#' @param loo An \R object of class \code{loo}. If \code{NULL},
#'   brms will try to extract a precomputed \code{loo} object
#'   from the fitted model, added there via \code{\link{add_criterion}}.
#' @param fit An \R object of class \code{brmsfit}.
#' @param k_threshold The threshold at which Pareto \eqn{k}
#'   estimates are treated as problematic. Defaults to \code{0.7}.
#'   See \code{\link[loo:pareto-k-diagnostic]{pareto_k_ids}}
#'   for more details.
#' @param check Logical; If \code{TRUE} (the default), some checks
#'   check are performed if the \code{loo} object was generated
#'   from the \code{brmsfit} object passed to argument \code{fit}.
#' @param recompile Logical, indicating whether the Stan model should be
#'   recompiled. This may be necessary if you are running \code{reloo} on
#'   another machine than the one used to fit the model.
#' @param future_args A list of further arguments passed to
#'   \code{\link[future:future]{future}} for additional control over parallel
#'   execution if activated.
#' @param ... Further arguments passed to
#'   \code{\link{update.brmsfit}} and \code{\link{log_lik.brmsfit}}.
#'
#' @return An object of the class \code{loo}.
#'
#' @details
#' Warnings about Pareto \eqn{k} estimates indicate observations
#' for which the approximation to LOO is problematic (this is described in
#' detail in Vehtari, Gelman, and Gabry (2017) and the
#' \pkg{\link[loo:loo-package]{loo}} package documentation).
#' If there are \eqn{J} observations with \eqn{k} estimates above
#' \code{k_threshold}, then \code{reloo} will refit the original model
#' \eqn{J} times, each time leaving out one of the \eqn{J}
#' problematic observations. The pointwise contributions of these observations
#' to the total ELPD are then computed directly and substituted for the
#' previous estimates from these \eqn{J} observations that are stored in the
#' original \code{loo} object.
#'
#' @seealso \code{\link{loo}}, \code{\link{kfold}}
#'
#' @examples
#' \dontrun{
#' fit1 <- brm(count ~ zAge + zBase * Trt + (1|patient),
#'             data = epilepsy, family = poisson())
#' # throws warning about some pareto k estimates being too high
#' (loo1 <- loo(fit1))
#' (reloo1 <- reloo(fit1, loo = loo1, chains = 1))
#' }
#'
#' @export
reloo.brmsfit <- function(x, loo = NULL, k_threshold = 0.7, newdata = NULL,
                          resp = NULL, check = TRUE, recompile = NULL,
                          future_args = list(), ...) {
  stopifnot(is.brmsfit(x), is.list(future_args))
  if (is.brmsfit_multiple(x)) {
    warn_brmsfit_multiple(x)
    class(x) <- "brmsfit"
  }
  loo <- loo %||% x$criteria[["loo"]]
  if (is.null(loo)) {
    stop2("No 'loo' object was provided and none is stored within the model.")
  } else if (!is.loo(loo)) {
    stop2("Inputs to the 'loo' argument must be of class 'loo'.")
  }
  if (is.null(newdata)) {
    mf <- model.frame(x)
  } else {
    mf <- as.data.frame(newdata)
  }
  mf <- rm_attr(mf, c("terms", "brmsframe"))
  if (NROW(mf) != NROW(loo$pointwise)) {
    stop2("Number of observations in 'loo' and 'x' do not match.")
  }
  check <- as_one_logical(check)
  if (check) {
    yhash_loo <- attr(loo, "yhash")
    yhash_fit <- hash_response(x, newdata = newdata)
    if (!is_equal(yhash_loo, yhash_fit)) {
      stop2(
        "Response values used in 'loo' and 'x' do not match. ",
        "If this is a false positive, please set 'check' to FALSE."
      )
    }
  }
  if (is.null(loo$diagnostics$pareto_k)) {
    stop2("No Pareto k estimates found in the 'loo' object.")
  }
  obs <- loo::pareto_k_ids(loo, k_threshold)
  J <- length(obs)
  if (J == 0L) {
    message(
      "No problematic observations found. ",
      "Returning the original 'loo' object."
    )
    return(loo)
  }

  # ensure that the model can be run in the current R session
  x <- recompile_model(x, recompile = recompile)

  # split dots for use in log_lik and update
  dots <- list(...)
  ll_arg_names <- arg_names("log_lik")
  ll_arg_names <- intersect(names(dots), ll_arg_names)
  ll_args <- dots[ll_arg_names]
  ll_args$allow_new_levels <- TRUE
  ll_args$sample_new_levels <-
    first_not_null(ll_args$sample_new_levels, "gaussian")
  ll_args$resp <- resp
  ll_args$combine <- TRUE
  # cores is used in both log_lik and update
  up_arg_names <- setdiff(names(dots), setdiff(ll_arg_names, "cores"))
  up_args <- dots[up_arg_names]
  up_args$object <- x
  up_args$refresh <- 0

  .reloo <- function(j) {
    message(
      "\nFitting model ", j, " out of ", J,
      " (leaving out observation ", obs[j], ")"
    )
    omitted <- obs[j]
    mf_omitted <- mf[-omitted, , drop = FALSE]
    up_args$newdata <- mf_omitted
    up_args$data2 <- subset_data2(x$data2, -omitted)
    fit_j <- SW(do_call(update, up_args))
    ll_args$object <- fit_j
    ll_args$newdata <- mf[omitted, , drop = FALSE]
    ll_args$newdata2 <- subset_data2(x$data2, omitted)
    return(do_call(log_lik, ll_args))
  }

  message(
    J, " problematic observation(s) found.",
    "\nThe model will be refit ", J, " times."
  )
  # TODO: separate parallel and non-parallel code to enable better printing?
  future_args$X <- seq_len(J)
  future_args$FUN <- .reloo
  future_args$future.seed <- TRUE
  lls <- do_call("future_lapply", future_args, pkg = "future.apply")

  # most of the following code is taken from rstanarm:::reloo
  # compute elpd_{loo,j} for each of the held out observations
  elpd_loo <- ulapply(lls, log_mean_exp)
  # compute \hat{lpd}_j for each of the held out observations (using log-lik
  # matrix from full posterior, not the leave-one-out posteriors)
  mf_obs <- mf[obs, , drop = FALSE]
  data2_obs <- subset_data2(x$data2, obs)
  ll_x <- log_lik(x, newdata = mf_obs, newdata2 = data2_obs)
  hat_lpd <- apply(ll_x, 2, log_mean_exp)
  # compute effective number of parameters
  p_loo <- hat_lpd - elpd_loo
  # replace parts of the loo object with these computed quantities
  sel <- c("elpd_loo", "p_loo", "looic")
  loo$pointwise[obs, sel] <- cbind(elpd_loo, p_loo, -2 * elpd_loo)
  new_pw <- loo$pointwise[, sel, drop = FALSE]
  loo$estimates[, 1] <- colSums(new_pw)
  loo$estimates[, 2] <- sqrt(nrow(loo$pointwise) * apply(new_pw, 2, var))
  # what should we do about pareto-k? for now setting them to 0
  loo$diagnostics$pareto_k[obs] <- 0
  loo
}

#' @rdname reloo.brmsfit
#' @export
reloo.loo <- function(x, fit, ...) {
  reloo(fit, loo = x, ...)
}

# the generic will eventually be moved to 'loo'
#' @rdname reloo.brmsfit
#' @export
reloo <- function(x, ...) {
  UseMethod("reloo")
}
