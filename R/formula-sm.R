# This file contains functions dealing with the extended
# formula syntax to specify smooth terms via mgcv

#' Defining smooths in \pkg{brms} formulas
#'
#' Functions used in definition of smooth terms within a model formulas.
#' The function does not evaluate a (spline) smooth - it exists purely
#' to help set up a model using spline based smooths.
#'
#' @param ... Arguments passed to \code{\link[mgcv:s]{mgcv::s}} or
#'  \code{\link[mgcv:t2]{mgcv::t2}}.
#'
#' @details The function defined here are just simple wrappers of the respective
#'   functions of the \pkg{mgcv} package. When using them, please cite the
#'   appropriate references obtained via \code{citation("mgcv")}.
#'
#'  \pkg{brms} uses the "random effects" parameterization of smoothing splines
#'  as explained in \code{\link[mgcv:gamm]{mgcv::gamm}}. A nice tutorial on this
#'  topic can be found in Pedersen et al. (2019). The answers provided in this
#'  \href{https://discourse.mc-stan.org/t/better-priors-non-flat-for-gams-brms/23012/4}{Stan discourse post}
#'  may also be helpful.
#'
#' @seealso \code{\link{brmsformula}},
#'   \code{\link[mgcv:s]{mgcv::s}}, \code{\link[mgcv:t2]{mgcv::t2}}
#'
#' @references
#' Pedersen, E. J., Miller, D. L., Simpson, G. L., & Ross, N. (2019).
#' Hierarchical generalized additive models in ecology: an introduction with
#' mgcv. PeerJ.
#'
#' @examples
#' \dontrun{
#' # simulate some data
#' dat <- mgcv::gamSim(1, n = 200, scale = 2)
#'
#' # fit univariate smooths for all predictors
#' fit1 <- brm(y ~ s(x0) + s(x1) + s(x2) + s(x3),
#'             data = dat, chains = 2)
#' summary(fit1)
#' plot(conditional_smooths(fit1), ask = FALSE)
#'
#' # fit a more complicated smooth model
#' fit2 <- brm(y ~ t2(x0, x1) + s(x2, by = x3),
#'             data = dat, chains = 2)
#' summary(fit2)
#' plot(conditional_smooths(fit2), ask = FALSE)
#' }
#'
#' @export
s <- function(...) {
  mgcv::s(...)
}

#' @rdname s
#' @export
t2 <- function(...) {
  mgcv::t2(...)
}

# extract information about smooth terms
# @param x either a formula or a list containing an element "sm"
# @param data data.frame containing the covariates
tidy_smef <- function(x, data = NULL) {
  if (is.formula(x)) {
    x <- brmsterms(x, check_response = FALSE)$dpars$mu
  }
  form <- x[["sm"]]
  if (!is.formula(form)) {
    return(empty_data_frame())
  }
  out <- data.frame(term = all_terms(form), stringsAsFactors = FALSE)
  nterms <- nrow(out)
  out$sfun <- get_matches("^[^\\(]+", out$term)
  out$vars <- out$byvars <- out$covars <- vector("list", nterms)
  for (i in seq_len(nterms)) {
    sm <- eval2(out$term[i])
    out$covars[[i]] <- sm$term
    if (sm$by != "NA") {
      out$byvars[[i]] <- sm$by
    }
    out$vars[[i]] <- c(out$covars[[i]], out$byvars[[i]])
  }
  out$label <- paste0(out$sfun, rename(ulapply(out$vars, collapse)))
  # prepare information inferred from the data
  sdata <- x$sdata$sm
  if (is.null(sdata)) {
    # TODO: remove once refactor is complete
    sdata <- data_sm(x, data)
  }
  bylevels <- attr(sdata$Xs, "bylevels")
  nby <- lengths(bylevels)
  tmp <- vector("list", nterms)
  for (i in seq_len(nterms)) {
    tmp[[i]] <- out[i, , drop = FALSE]
    tmp[[i]]$termnum <- i
    if (nby[i] > 0L) {
      tmp[[i]] <- do_call(rbind, repl(tmp[[i]], nby[i]))
      tmp[[i]]$bylevel <- rm_wsp(bylevels[[i]])
      tmp[[i]]$byterm <- paste0(tmp[[i]]$term, tmp[[i]]$bylevel)
      str_add(tmp[[i]]$label) <- rename(tmp[[i]]$bylevel)
    } else {
      tmp[[i]]$bylevel <- NA
      tmp[[i]]$byterm <- tmp[[i]]$term
    }
  }
  out <- do_call(rbind, tmp)
  out$knots <- sdata[grepl("^knots_", names(sdata))]
  out$nbases <- lengths(out$knots)
  attr(out, "Xs_names") <- colnames(sdata$Xs)
  rownames(out) <- NULL
  out
}

# check if smooths are present in the model
has_smooths <- function(bterms) {
  length(get_effect(bterms, target = "sm")) > 0L
}
