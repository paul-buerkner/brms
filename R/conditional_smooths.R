#' Display Smooth Terms
#' 
#' Display smooth \code{s} and \code{t2} terms of models
#' fitted with \pkg{brms}.
#' 
#' @aliases marginal_smooths marginal_smooths.brmsfit
#' 
#' @inheritParams conditional_effects.brmsfit
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
#' an object of class \code{brms_conditional_effects}. See
#' \code{\link{conditional_effects}} for 
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
#' plot(conditional_smooths(fit), rug = TRUE, ask = FALSE)
#' # show only the smooth term s(x2)
#' plot(conditional_smooths(fit, smooths = "s(x2)"), ask = FALSE)
#' 
#' # fit and plot a two-dimensional smooth term
#' fit2 <- brm(y ~ t2(x0, x2), data = dat)
#' ms <- conditional_smooths(fit2)
#' plot(ms, stype = "contour")
#' plot(ms, stype = "raster")
#' }
#' 
#' @export
conditional_smooths.brmsfit <- function(x, smooths = NULL,
                                        int_conditions = NULL,
                                        prob = 0.95, spaghetti = FALSE,
                                        resolution = 100, too_far = 0,
                                        subset = NULL, nsamples = NULL,
                                        probs = NULL, ...) {
  probs <- validate_ci_bounds(prob, probs = probs)
  spaghetti <- as_one_logical(spaghetti)
  contains_samples(x)
  x <- restructure(x)
  x <- exclude_terms(x, incl_autocor = FALSE)
  smooths <- rm_wsp(as.character(smooths))
  conditions <- prepare_conditions(x)
  subset <- subset_samples(x, subset, nsamples)
  # call as.matrix only once to save time and memory
  samples <- as.matrix(x, subset = subset)
  bterms <- brmsterms(exclude_terms(x$formula, smooths_only = TRUE))
  out <- conditional_smooths(
    bterms, fit = x, samples = samples, smooths = smooths, 
    conditions = conditions, int_conditions = int_conditions, 
    too_far = too_far, resolution = resolution, probs = probs, 
    spaghetti = spaghetti
  )
  if (!length(out)) {
    stop2("No valid smooth terms found in the model.")
  }
  structure(out, class = "brms_conditional_effects", smooths_only = TRUE)
}

#' @rdname conditional_smooths.brmsfit
#' @export
conditional_smooths <- function(x, ...) {
  UseMethod("conditional_smooths")
}

#' @export
conditional_smooths.default <- function(x, ...) {
  NULL
}

#' @export
conditional_smooths.mvbrmsterms <- function(x, ...) {
  out <- list()
  for (r in names(x$terms)) {
    c(out) <- conditional_smooths(x$terms[[r]], ...)
  }
  out
}

#' @export
conditional_smooths.brmsterms <- function(x, ...) {
  out <- list()
  for (dp in names(x$dpars)) {
    c(out) <- conditional_smooths(x$dpars[[dp]], ...)
  }
  for (nlp in names(x$nlpars)) {
    c(out) <- conditional_smooths(x$nlpars[[nlp]], ...)
  }
  out
}

# conditional smooths for a single predicted parameter
# @param fit brmsfit object
# @param samples extract posterior samples
# @param smooths optional names of smooth terms to plot
# @param conditions output of prepare_conditions
# @param int_conditions values of by-vars at which to evalute smooths
# @param ...: currently ignored
# @return a named list with one element per smooth term
#' @export
conditional_smooths.btl <- function(x, fit, samples, smooths,
                                    conditions, int_conditions, 
                                    probs, resolution, too_far,
                                    spaghetti, ...) {
  stopifnot(is.brmsfit(fit))
  out <- list()
  mf <- model.frame(fit)
  smef <- tidy_smef(x, mf)
  smterms <- unique(smef$term)
  if (!length(smooths)) {
    I <- seq_along(smterms)
  } else {
    I <- which(smterms %in% smooths)
  }
  for (i in I) {
    # loop over smooth terms and compute their predictions
    term <- smterms[i]
    sub_smef <- subset2(smef, term = term)
    # extract raw variable names before transformations
    covars <- all_vars(sub_smef$covars[[1]])
    byvars <- all_vars(sub_smef$byvars[[1]])
    ncovars <- length(covars)
    if (ncovars > 2L) {
      byvars <- c(covars[3:ncovars], byvars)
      covars <- covars[1:2]
      ncovars <- 2L
    }
    vars <- c(covars, byvars)
    values <- named_list(vars)
    is_numeric <- setNames(rep(FALSE, ncovars), covars)
    for (cv in covars) {
      if (is.numeric(mf[[cv]])) {
        is_numeric[cv] <- TRUE
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
        int_cond <- int_conditions[[cv]]
        if (is.function(int_cond)) {
          int_cond <- int_cond(mf[[cv]])
        }
        values[[cv]] <- int_cond
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
    for (v in other_vars) {
      cval <- conditions[1, v]
      if (length(dim(cval)) == 2L) {
        # matrix columns don't have automatic broadcasting apparently
        cval <- matrix(cval, nrow(newdata), ncol(cval), byrow = TRUE)
      }
      newdata[[v]] <- cval
    }
    sdata <- standata(
      fit, newdata, re_formula = NA, 
      internal = TRUE, check_response = FALSE
    )
    prep_args <- nlist(x, samples, sdata, data = mf)
    prep <- do_call(prepare_predictions, prep_args)
    J <- which(smef$termnum == i)
    scs <- unlist(attr(prep$sm$fe$Xs, "smcols")[J])
    prep$sm$fe$Xs <- prep$sm$fe$Xs[, scs, drop = FALSE]
    prep$sm$fe$bs <- prep$sm$fe$bs[, scs, drop = FALSE]
    prep$sm$re <- prep$sm$re[J]
    prep$family <- brmsfamily("gaussian")
    eta <- predictor(prep, i = NULL)
    effects <- na.omit(sub_smef$covars[[1]][1:2])
    marg_data <- add_effects__(newdata[, vars, drop = FALSE], effects)
    if (length(byvars)) {
      # byvars will be plotted as facets
      marg_data$cond__ <- rows2labels(marg_data[, byvars, drop = FALSE]) 
    } else {
      marg_data$cond__ <- factor(1)
    }
    spa_data <- NULL
    if (spaghetti && ncovars == 1L && is_numeric[1]) {
      sample <- rep(seq_rows(eta), each = ncol(eta))
      spa_data <- data.frame(as.numeric(t(eta)), factor(sample))
      colnames(spa_data) <- c("estimate__", "sample__")
      spa_data <- cbind(marg_data, spa_data)
    }
    eta <- posterior_summary(eta, robust = TRUE, probs = probs)
    colnames(eta) <- c("estimate__", "se__", "lower__", "upper__")
    eta <- cbind(marg_data, eta)
    response <- combine_prefix(x, keep_mu = TRUE)
    response <- paste0(response, ": ", term)
    points <- mf[, vars, drop = FALSE]
    points <- add_effects__(points, covars)
    attr(eta, "response") <- response
    attr(eta, "effects") <- effects
    attr(eta, "surface") <- all(is_numeric) && ncovars == 2L
    attr(eta, "spaghetti") <- spa_data
    attr(eta, "points") <- points
    out[[response]] <- eta
  }
  out
}

# the name 'marginal_smooths' is deprecated as of brms 2.10.3
# do not remove it eventually as it has been used in the brms papers
#' @export
marginal_smooths <- function(x, ...) {
  UseMethod("marginal_smooths")
}

#' @export
marginal_smooths.brmsfit <- function(x, ...) {
  warning2("Method 'marginal_smooths' is deprecated. ",
           "Please use 'conditional_smooths' instead.")
  conditional_smooths.brmsfit(x, ...)
}

