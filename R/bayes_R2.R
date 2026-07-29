#' Compute a Bayesian version of R-squared for regression models
#'
#' @aliases bayes_R2
#'
#' @inheritParams predict.brmsfit
#' @param method Character string specifying how R2 is computed. Either
#'   `"model"` to use the model-based variances, or `"residual"` to
#'   use the residual-based variances. If `NULL` (default), the
#'   model-based R2 is computed where possible, falling back to residual-based
#'   R2 for families without current support of model-based variances.
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
#' @details For an introduction to the approach, see Gelman et al. (2019)
#'  and \url{https://github.com/jgabry/bayes_R2/}.
#'  For Gaussian and Bernoulli models, \code{bayes_R2} uses model-based residual
#'  variances as proposed in Gelman et al. (2019), with Bernoulli models using 
#'  Tjur's pseudo-variance. For Gaussian models with heteroscedastic
#'  sigma, the mean residual variance is used as an approximation (see Tjur 
#'  (2009) for discussion on this approximation). For other families, 
#'  \code{bayes_R2} warns and falls back to residual-based variances.
#'
#' @references Andrew Gelman, Ben Goodrich, Jonah Gabry & Aki Vehtari. (2019).
#'   R-squared for Bayesian regression models, \emph{The American Statistician},
#'   73(3):307-309. \code{10.1080/00031305.2018.1549100} (Preprint available at
#'   \url{https://acris.aalto.fi/ws/portalfiles/portal/34206843/bayes_R2_v3.pdf})
#' 
#'   Tue Tjur. (2009). Coefficient of determination in logistic regression 
#'   models - A new proposal: The coefficient of discrimination, \emph{The 
#'   American Statistician}, 63:366-372. \code{10.1198/tast.2009.08210}
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
bayes_R2.brmsfit <- function(object, resp = NULL, method = NULL, summary = TRUE,
                             robust = FALSE, probs = c(0.025, 0.975), ...) {
  contains_draws(object)
  object <- restructure(object)
  resp <- validate_resp(resp, object)
  summary <- as_one_logical(summary)
  if (!is.null(method)) {
    rlang::arg_match(method, values = c("model", "residual"))
  }
  # check for precomputed values
  R2 <- get_criterion(object, "bayes_R2")
  has_stored <- is.matrix(R2)
  further_arg_names <- c("resp")
  use_stored <- !length(list(...)) &&
    !any(further_arg_names %in% names(match.call()))
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
  args_ypred <- list(object, sort = TRUE, ...)
  R2 <- named_list(paste0("R2", resp))
  warned_families <- character(0)
  # TODO: find supported families automatically once more are supported
  model_variance_families <- c("gaussian", "bernoulli")
  
  for (i in seq_along(R2)) {
    # assumes expectations of different responses to be independent
    args_ypred$resp <- resp[i]
    ypred <- do_call(posterior_epred, args_ypred)
    if (is_ordinal(family(object, resp = resp[i]))) {
      ypred <- ordinal_probs_continuous(ypred)
    }
    family_name <- family(object, resp = resp[i])$family

    method_i <- method
    if (is.null(method_i)) {
      method_i <- str_if(
        family_name %in% model_variance_families, 
        "model", "residual"
      )
      if (method_i == "residual") {
        warning2(
          "No model-based residual variance is currently implemented for ",
          "family '", family_name, "'\nin 'bayes_R2'. ",
          "Falling back to residual-based R2 computation."
        )
        warned_families <- c(warned_families, family_name)
      }
    }
    
    if (method_i == "model") {
      R2[[i]] <- bayes_R2_model(ypred, family_name, args_ypred = args_ypred, ...)
    } else if (method_i == "residual") {
      R2[[i]] <- bayes_R2_residual(ypred, object, resp = resp[i], ...)
    }
  }

  R2 <- do_call(cbind, R2)
  colnames(R2) <- paste0("R2", resp)
  if (summary) {
    R2 <- posterior_summary(R2, probs = probs, robust = robust)
  }
  R2
}

bayes_R2_model <- function(ypred, family_name, ...) {
  var_res_fun <- get(paste0(".var_res_", family_name), mode = "function")
  var_res <- var_res_fun(ypred, ...)
  .bayes_R2(ypred, var_res)
}

bayes_R2_residual <- function(ypred, object, resp, ...) {
  y <- do_call(get_y, c(list(object, warn = TRUE, resp = resp), list(...)))
  res <- -1 * sweep(ypred, 2, y)
  var_res <- matrixStats::rowVars(res)
  .bayes_R2(ypred, var_res)
}

.bayes_R2 <- function(ypred, var_res) {
  var_ypred <- matrixStats::rowVars(ypred)
  as.matrix(var_ypred / (var_ypred + var_res))
}

# ------------------- family-specific residual variances -----------------
.var_res_gaussian <- function(ypred, args_ypred, ...) {
  args_sigma <- args_ypred
  args_sigma$dpar <- "sigma"
  sigma <- do_call(posterior_epred, args_sigma)
  # use mean of heteroscedastic sigma as approximate 
  # (see Tjur (2009) for discussion)
  matrixStats::rowMeans2(sigma)^2
}

.var_res_bernoulli <- function(ypred, ...) {
  matrixStats::rowMeans2(ypred * (1 - ypred))
}

