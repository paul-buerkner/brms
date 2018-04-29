#' Display Smooth Terms
#' 
#' Display smooth \code{s} and \code{t2} terms of models
#' fitted with \pkg{brms}.
#' 
#' @inheritParams marginal_effects
#' @param smooths Optional character vector of smooth terms
#'   to display. If \code{NULL} (the default) all smooth terms
#'   are shown.
#' @param subset A numeric vector specifying
#'  the posterior samples to be used. 
#'  If \code{NULL} (the default), all samples are used.
#' @param nsamples Positive integer indicating how many 
#'  posterior samples should be used. 
#'  If \code{NULL} (the default) all samples are used.
#'  Ignored if \code{subset} is not \code{NULL}.
#' @param ... Currently ignored.
#'   
#' @return For the \code{brmsfit} method, 
#' an object of class \code{brmsMarginalEffects}. See
#' \code{\link{marginal_effects}} for 
#' more details and documentation of the related plotting function.
#' 
#' @details Two-dimensional smooth terms will be visualized using
#'   either contour or raster plots.
#'   
#' @examples 
#' \dontrun{
#' set.seed(0) 
#' dat <- mgcv::gamSim(1, n = 200, scale = 2)
#' fit <- brm(y ~ s(x0) + s(x1) + s(x2) + s(x3), data = dat)
#' # show all smooth terms
#' plot(marginal_smooths(fit), rug = TRUE, ask = FALSE)
#' # show only the smooth term s(x2)
#' plot(marginal_smooths(fit, smooths = "s(x2)"), ask = FALSE)
#' 
#' # fit and plot a two-dimensional smooth term
#' fit2 <- brm(y ~ t2(x0, x2), data = dat)
#' ms <- marginal_smooths(fit2)
#' plot(ms, stype = "contour")
#' plot(ms, stype = "raster")
#' }
#' 
#' @export
marginal_smooths <- function(x, ...) {
  UseMethod("marginal_smooths")
}

marginal_smooths_internal <- function(x, ...) {
  # compute predictions based on the smooths terms only
  UseMethod("marginal_smooths_internal")
}

#' @export
marginal_smooths_internal.default <- function(x, ...) {
  NULL
}

#' @export
marginal_smooths_internal.mvbrmsterms <- function(x, ...) {
  out <- list()
  for (r in names(x$terms)) {
    out <- c(out, marginal_smooths_internal(x$terms[[r]], ...))
  }
  out
}

#' @export
marginal_smooths_internal.brmsterms <- function(x, ...) {
  out <- list()
  for (dp in names(x$dpars)) {
    if (is.btl(x$dpars[[dp]])) {
      btl <- x$dpars[[dp]]
      out <- c(out, marginal_smooths_internal(btl, ...))
    } else if (is.btnl(x$dpars[[dp]])) {
      for (nlp in names(x$dpars[[dp]]$nlpars)) {
        btl <- x$dpars[[dp]]$nlpars[[nlp]]
        out <- c(out, marginal_smooths_internal(btl, ...))
      }
    }
  }
  out
}

#' @export
marginal_smooths_internal.btl <- function(x, fit, samples, smooths,
                                          conditions, int_conditions, 
                                          probs, resolution, too_far,
                                          spaghetti, ...) {
  # Args:
  #   fit: brmsfit object
  #   samples: extract posterior samples
  #   smooths: optional names of smooth terms to plot
  #   conditions: output of prepare_conditions
  #   int_conditions: values of by-vars at which to evalute smooths
  #   ...: currently ignored
  stopifnot(is.brmsfit(fit))
  out <- list()
  mf <- model.frame(fit)
  too_many_covars <- FALSE
  smef <- tidy_smef(x, mf)
  smterms <- unique(smef$term) 
  for (i in seq_along(smterms)) {
    # loop over smooth terms and compute their predictions
    term <- smterms[i]
    sub_smef <- subset2(smef, term = term)
    byvars <- sub_smef$byvars[[1]]
    covars <- sub_smef$covars[[1]]
    vars <- sub_smef$vars[[1]]
    ncovars <- length(covars)
    if (ncovars > 2L) {
      too_many_covars <- TRUE
    }
    include_smooth <- !length(smooths) || term %in% smooths
    if (include_smooth && !too_many_covars) {
      values <- named_list(vars)
      for (cv in covars) {
        if (is.numeric(mf[[cv]])) {
          values[[cv]] <- seq(
            min(mf[[cv]]), max(mf[[cv]]), 
            length.out = resolution
          )
        } else {
          values[[cv]] <- levels(factor(mf[[cv]]))
        }
      }
      for (cv in byvars) {
        if (cv %in% names(int_conditions)) {
          values[[cv]] <- int_conditions[[cv]]
        } else if (is.numeric(mf[[cv]])) {
          mean2 <- mean(mf[[cv]], na.rm = TRUE)
          sd2 <- sd(mf[[cv]], na.rm = TRUE)
          values[[cv]] <- (-1:1) * sd2 + mean2
        } else {
          values[[cv]] <- levels(factor(mf[[cv]]))
        }
      }
      newdata <- expand.grid(values)
      if (ncovars == 2L && too_far > 0) {
        # exclude prediction grid points too far from data
        ex_too_far <- mgcv::exclude.too.far(
          g1 = newdata[[covars[1]]], 
          g2 = newdata[[covars[2]]], 
          d1 = mf[, covars[1]],
          d2 = mf[, covars[2]],
          dist = too_far
        )
        newdata <- newdata[!ex_too_far, ]  
      }
      other_vars <- setdiff(names(conditions), vars)
      newdata[, other_vars] <- conditions[1, other_vars]
      sdata <- standata(
        fit, newdata, re_formula = NA, 
        internal = TRUE, check_response = FALSE
      )
      draws_args <- nlist(x, samples, sdata, data = mf, smooths_only = TRUE)
      draws <- do.call(extract_draws, draws_args)
      J <- which(smef$termnum == i)
      scs <- unlist(attr(draws$fe$X, "smcols")[J])
      draws$fe$X <- draws$fe$X[, scs, drop = FALSE]
      draws$fe$b <- draws$fe$b[, scs, drop = FALSE]
      draws$sm <- draws$sm[J]
      eta <- predictor(draws, i = NULL)
      spa_data <- NULL
      if (spaghetti && ncovars == 1L) {
        sample <- rep(seq_len(nrow(eta)), each = ncol(eta))
        spa_data <- data.frame(as.numeric(t(eta)), factor(sample))
        colnames(spa_data) <- c("estimate__", "sample__")
        spa_data <- cbind(newdata[, vars, drop = FALSE], spa_data)
      }
      eta <- posterior_summary(eta, robust = TRUE, probs = probs)
      colnames(eta) <- c("estimate__", "se__", "lower__", "upper__")
      eta <- cbind(newdata[, vars, drop = FALSE], eta)
      if (length(byvars)) {
        # byvars will be plotted as facets
        bynumeric <- byvars[ulapply(eta[byvars], is.numeric)]
        for (cv in bynumeric) {
          eta[[cv]] <- round(eta[[cv]], 2)
        }
        eta$cond__ <- Reduce(paste_colon, eta[, byvars, drop = FALSE]) 
      }
      response <- combine_prefix(x, keep_mu = TRUE)
      response <- paste0(response, ": ", term)
      attr(eta, "response") <- response
      attr(eta, "effects") <- covars
      attr(eta, "surface") <- ncovars == 2L
      attr(eta, "spaghetti") <- spa_data
      attr(eta, "points") <- mf[, vars, drop = FALSE]
      attr(eta, "too_many_covars") <- too_many_covars
      out[[response]] <- eta
    }
  }
  out
}
