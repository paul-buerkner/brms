#' Compute a Bayesian version of R-squared for regression models
#'
#' @aliases bayes_R2
#'
#' @inheritParams predict.brmsfit
#' @param ... Further arguments passed to
#'   \code{\link[brms:posterior_epred.brmsfit]{posterior_epred}},
#'   which is used in the computation of the R-squared values.
#'
#' @return If \code{summary = TRUE}, an M x C matrix is returned
#'  (M = number of response variables and c = \code{length(probs) + 2})
#'  containing summary statistics of the Bayesian R-squared values.
#'  If \code{summary = FALSE}, the posterior draws of the Bayesian
#'  R-squared values are returned in an S x M matrix (S is the number of draws).
#'
#' @details For an introduction to the approach, see Gelman et al. (2018)
#'  and \url{https://github.com/jgabry/bayes_R2/}.
#'
#' @references Andrew Gelman, Ben Goodrich, Jonah Gabry & Aki Vehtari. (2018).
#'   R-squared for Bayesian regression models, \emph{The American Statistician}.
#'   \code{10.1080/00031305.2018.1549100} (Preprint available at
#'   \url{https://stat.columbia.edu/~gelman/research/published/bayes_R2_v3.pdf})
#'
#' @examples
#' \dontrun{
#' fit <- brm(mpg ~ wt + cyl, data = mtcars)
#' summary(fit)
#' bayes_R2(fit)
#'
#' # compute R2 with new data
#' nd <- data.frame(mpg = c(10, 20, 30), wt = c(4, 3, 2), cyl = c(8, 6, 4))
#' bayes_R2(fit, newdata = nd)
#' }
#'
#' @method bayes_R2 brmsfit
#' @importFrom rstantools bayes_R2
#' @export bayes_R2
#' @export
bayes_R2.brmsfit <- function(object, resp = NULL, summary = TRUE,
                             robust = FALSE, probs = c(0.025, 0.975), ...) {
  contains_draws(object)
  object <- restructure(object)
  resp <- validate_resp(resp, object)
  summary <- as_one_logical(summary)
  # check for precomputed values
  R2 <- get_criterion(object, "bayes_R2")
  has_stored <- is.matrix(R2)
  use_stored <- as.logical(length(list(...)))
  if (has_stored && !use_stored) {
    message("Recomputing 'bayes_R2'")
  }
  if (has_stored && use_stored) {
    # assumes unsummarized 'R2' as ensured by 'add_criterion'
    take <- colnames(R2) %in% paste0("R2", resp)
    R2 <- R2[, take, drop = FALSE]
    if (summary) {
      R2 <- posterior_summary(R2, probs = probs, robust = robust)
    }
    return(R2)
  }
  family <- family(object, resp = resp)
  if (conv_cats_dpars(family)) {
    stop2("'bayes_R2' is not defined for unordered categorical models.")
  }
  if (is_ordinal(family)) {
    warning2(
      "Predictions are treated as continuous variables in ",
      "'bayes_R2' which is likely invalid for ordinal families."
    )
  }
  args_y <- list(object, warn = TRUE, ...)
  args_ypred <- list(object, sort = TRUE, ...)
  R2 <- named_list(paste0("R2", resp))
  for (i in seq_along(R2)) {
    # assumes expectations of different responses to be independent
    args_ypred$resp <- args_y$resp <- resp[i]
    y <- do_call(get_y, args_y)
    ypred <- do_call(posterior_epred, args_ypred)
    if (is_ordinal(family(object, resp = resp[i]))) {
      ypred <- ordinal_probs_continuous(ypred)
    }
    R2[[i]] <- .bayes_R2(y, ypred)
  }
  R2 <- do_call(cbind, R2)
  colnames(R2) <- paste0("R2", resp)
  if (summary) {
    R2 <- posterior_summary(R2, probs = probs, robust = robust)
  }
  R2
}

# internal function of bayes_R2.brmsfit
# see https://github.com/jgabry/bayes_R2/blob/master/bayes_R2.pdf
.bayes_R2 <- function(y, ypred, ...) {
  e <- -1 * sweep(ypred, 2, y)
  var_ypred <- matrixStats::rowVars(ypred)
  var_e <- matrixStats::rowVars(e)
  as.matrix(var_ypred / (var_ypred + var_e))
}
