# Evaluate the fitted baseline without extrapolating its spline basis.
# Paired evaluation associates each time with one row of simplex draws;
# common-time evaluation shares a time grid across all draws.
bhaz_at <- function(t, sbhaz, basis, tail = "constant", integrate = FALSE,
                    paired = TRUE) {
  tail <- match.arg(tail, c("constant", "loglinear", "none"))
  stopifnot(inherits(basis, "MSpline"), is.matrix(sbhaz))
  bounds <- attr(basis, "Boundary.knots")
  t <- as.numeric(t)
  if (paired && nrow(sbhaz) == 1L && length(t) > 1L) {
    sbhaz <- sbhaz[rep(1L, length(t)), , drop = FALSE]
  }
  S <- nrow(sbhaz)
  if (paired) {
    if (!length(t) %in% c(1L, S)) {
      stop2("Times must have length 1 or the number of draws.")
    }
    t <- rep_len(t, S)
    out <- numeric(S)
  } else {
    out <- matrix(0, S, length(t))
  }
  inner <- which(t >= bounds[1] & t <= bounds[2])
  if (length(inner)) {
    Z <- bhaz_basis_matrix(t[inner], basis = basis, integrate = integrate)
    if (paired) {
      out[inner] <- matrixStats::rowSums2(
        Z * sbhaz[inner, , drop = FALSE]
      )
    } else {
      out[, inner] <- tcrossprod(sbhaz, Z)
    }
  }
  if (integrate) {
    out <- pmin(pmax(out, 0), 1)
  }
  outer <- which(t > bounds[2])
  if (length(outer)) {
    h <- as.vector(sbhaz %*%
      as.vector(bhaz_basis_matrix(bounds[2], basis = basis)))
    if (tail != "none" && any(!is.finite(h) | h <= 0)) {
      stop2("The Cox tail requires a finite positive boundary hazard.")
    }
    if (paired) {
      h <- h[outer]
      tt <- t[outer]
    } else {
      h <- matrix(h, S, length(outer))
      tt <- data2draws(t[outer], dim(h))
    }
    value <- switch(tail,
      constant = if (integrate) 1 + (tt - bounds[2]) * h else h,
      loglinear = if (integrate) {
        exp(bounds[2] * h * log(tt / bounds[2]))
      } else {
        h * (tt / bounds[2])^(bounds[2] * h - 1)
      },
      none = h * NA_real_
    )
    if (paired) out[outer] <- value else out[, outer] <- value
  }
  if (paired) out[is.na(t)] <- NA_real_ else out[, is.na(t)] <- NA_real_
  out
}

# Conditional cumulative hazards avoid subtraction of saturated CDFs.
cox_truncation <- function(mu, sbhaz, basis, tail, lb = NULL, ub = NULL) {
  Hlo <- if (is.null(lb) || identical(lb, -Inf)) {
    rep(0, nrow(sbhaz))
  } else {
    mu * bhaz_at(lb, sbhaz, basis, tail, integrate = TRUE)
  }
  Hup <- if (is.null(ub) || identical(ub, Inf)) {
    rep(Inf, nrow(sbhaz))
  } else {
    mu * bhaz_at(ub, sbhaz, basis, tail, integrate = TRUE)
  }
  if (anyNA(Hlo) || anyNA(Hup)) {
    stop2("Cox truncation bounds must lie within the spline support ",
          "when tail = 'none'.")
  }
  width <- Hup - Hlo
  if (any(!is.finite(Hlo) | width <= 0 | is.na(width))) {
    stop2("Invalid truncation bounds: no probability mass between 'lb' and 'ub'.")
  }
  nlist(Hlo, Hup, width)
}

# Invert H0 using bracketed Newton steps, shrinking the active set.
qcox_baseline <- function(target, sbhaz, basis, tail, tol = 1e-10,
                          maxit = 100L, lb = NULL, ub = NULL) {
  bounds <- attr(basis, "Boundary.knots")
  tol <- rep_len(tol, length(target))
  out <- rep(NA_real_, length(target))
  out[which(target == 0)] <- bounds[1]
  out[which(target == 1)] <- bounds[2]
  beyond <- which(target > 1)
  if (length(beyond)) {
    h <- bhaz_at(bounds[2], sbhaz, basis)[beyond]
    if (tail != "none" && any(!is.finite(h) | h <= 0)) {
      stop2("The Cox tail requires a finite positive boundary hazard.")
    }
    out[beyond] <- switch(tail,
      constant = bounds[2] + (target[beyond] - 1) / h,
      loglinear = bounds[2] *
        exp(log(target[beyond]) / (bounds[2] * h)),
      none = -bounds[2]
    )
  }
  active <- which(target > 0 & target < 1)
  if (length(active)) {
    grid <- seq(bounds[1], bounds[2], length.out = 25L)
    Hg <- bhaz_at(grid, sbhaz, basis, integrate = TRUE, paired = FALSE)
    pos <- max.col(Hg[active, , drop = FALSE] < target[active],
                   ties.method = "last")
    pos <- pmin(pmax(pos, 1L), length(grid) - 1L)
    lo <- rep(bounds[1], length(target))
    hi <- rep(bounds[2], length(target))
    lo[active] <- grid[pos]
    hi[active] <- grid[pos + 1L]
    Hlo <- Hg[cbind(active, pos)]
    Hhi <- Hg[cbind(active, pos + 1L)]
    if (!is.null(lb)) {
      lo[active] <- pmax(lo[active], rep_len(lb, length(target))[active])
    }
    if (!is.null(ub)) {
      hi[active] <- pmin(hi[active], rep_len(ub, length(target))[active])
    }
    if (!is.null(lb) || !is.null(ub)) {
      sb <- sbhaz[active, , drop = FALSE]
      Hlo <- bhaz_at(lo[active], sb, basis, integrate = TRUE)
      Hhi <- bhaz_at(hi[active], sb, basis, integrate = TRUE)
    }
    w <- (target[active] - Hlo) / pmax(Hhi - Hlo, .Machine$double.xmin)
    out[active] <- lo[active] +
      pmin(pmax(w, 0), 1) * (hi[active] - lo[active])
    for (iter in seq_len(maxit + 1L)) {
      sb <- sbhaz[active, , drop = FALSE]
      f <- bhaz_at(out[active], sb, basis, integrate = TRUE) - target[active]
      # Adjacent representable times limit the attainable probability accuracy.
      resolved <- hi[active] - lo[active] <= .Machine$double.eps *
        pmax(abs(lo[active]), abs(hi[active]))
      keep <- abs(f) > tol[active] & !resolved
      active <- active[keep]
      if (!length(active)) break
      if (iter > maxit) {
        stop2("Cox quantile inversion did not converge after ", maxit,
              " iterations.")
      }
      f <- f[keep]
      x <- out[active]
      lo[active] <- ifelse(f < 0, x, lo[active])
      hi[active] <- ifelse(f > 0, x, hi[active])
      step <- x - f / bhaz_at(x, sbhaz[active, , drop = FALSE], basis)
      bad <- !is.finite(step) | step <= lo[active] | step >= hi[active]
      step[bad] <- (lo[active][bad] + hi[active][bad]) / 2
      out[active] <- step
    }
  }
  out
}

# Keep old likelihood preparation usable while requiring the canonical
# response-level basis and fitting version for new time-aware features.
check_cox_prep <- function(prep) {
  if (!is_cox(prep$family) || is.mixfamily(prep$family) ||
      isTRUE(prep$bhaz$mixture)) {
    stop2("This method requires a non-mixture Cox response.")
  }
  if (!isTRUE(prep$bhaz$version >= "2.23.2") ||
      !inherits(prep$bhaz$basis, "MSpline")) {
    stop2("Cox predictions require a model fitted with brms 2.23.2 or later ",
          "and its stored baseline basis. Please refit your model.")
  }
  invisible(NULL)
}

posterior_predict_cox_trunc <- function(i, prep, output, args, q = NULL,
                                       p = NULL, lower.tail = TRUE,
                                       log.p = FALSE, log = FALSE, ...) {
  args$lb <- prep$data$lb[i]
  args$ub <- prep$data$ub[i]
  if (is.null(q)) q <- prep$data$Y[i]
  if (output == "quantile" && is.null(p)) {
    stop2("Argument 'p' must be specified when output = 'quantile'.")
  }
  switch(output,
    random = do_call(rcox, c(list(n = prep$ndraws), args)),
    probability = ,
    pit = do_call(pcox, c(list(q = q, lower.tail = lower.tail, log.p = log.p),
                          args)),
    density = do_call(dcox, c(list(x = q, log = log), args)),
    quantile = do_call(qcox, c(list(p = p, lower.tail = lower.tail,
                                  log.p = log.p), args))
  )
}

prepare_cox_predictions <- function(object, ..., check_response = FALSE) {
  if (is.brmsprep(object)) {
    if (length(list(...))) {
      stop2("Prediction arguments must be supplied when preparing the object.")
    }
    prep <- object
  } else {
    prep <- prepare_predictions(object, check_response = check_response, ...)
  }
  if (!is.brmsprep(prep)) {
    stop2("Please select a single Cox response via 'resp'.")
  }
  check_cox_prep(prep)
  prep$dpars$mu <- get_dpar(prep, "mu")
  if (is.null(dim(prep$dpars$mu))) {
    prep$dpars$mu <- matrix(prep$dpars$mu, prep$ndraws, prep$nobs)
  }
  if (any(!is.finite(prep$dpars$mu) | prep$dpars$mu <= 0)) {
    stop2("Cox predictions require finite positive values of 'mu'.")
  }
  prep
}

#' Cox Baseline and Conditional Survival Quantities
#'
#' Evaluate the fitted baseline or conditional hazard, cumulative hazard,
#' and survival function without changing the fitted spline basis.
#'
#' @param t Times, either a scalar or one per prediction observation.
#'   Baseline functions accept a common vector of times.
#' @param object A Cox \code{brmsfit}, or a single-response
#'   \code{brmsprep} returned by \code{\link{prepare_predictions}}.
#' @param stratum Index of the baseline hazard stratum.
#' @param sort See \code{\link{posterior_predict.brmsfit}}.
#' @param ... For fitted objects, arguments passed to
#'   \code{\link{prepare_predictions}}, including \code{newdata},
#'   \code{resp}, \code{re_formula}, \code{draw_ids}, and \code{tail}.
#' @return Baseline functions return draws by times matrices.
#'   Conditional functions return draws by observations matrices.
#' @details The tail policy is described in
#'   \code{\link{posterior_predict.brmsfit}}. Conditional quantities are for
#'   the untruncated event-time distribution; RMST and predictive CDFs
#'   additionally condition on any model truncation bounds.
#' @name cox_survival
NULL

#' @rdname cox_survival
#' @export
cox_bhaz <- function(t, object, stratum = 1L, ...) {
  prep <- prepare_cox_predictions(object, ...)
  stratum <- as_one_integer(stratum)
  if (stratum < 1L || stratum > length(prep$bhaz$sbhaz)) {
    stop2("Invalid baseline hazard stratum.")
  }
  bhaz_at(t, prep$bhaz$sbhaz[[stratum]], prep$bhaz$basis,
           prep$bhaz$tail, paired = FALSE)
}

#' @rdname cox_survival
#' @export
cox_cumbhaz <- function(t, object, stratum = 1L, ...) {
  prep <- prepare_cox_predictions(object, ...)
  stratum <- as_one_integer(stratum)
  if (stratum < 1L || stratum > length(prep$bhaz$sbhaz)) {
    stop2("Invalid baseline hazard stratum.")
  }
  bhaz_at(t, prep$bhaz$sbhaz[[stratum]], prep$bhaz$basis,
           prep$bhaz$tail, integrate = TRUE, paired = FALSE)
}

cox_quantity <- function(t, prep, integrate, sort) {
  if (!length(t) %in% c(1L, prep$nobs)) {
    stop2("'t' must have length 1 or the number of observations.")
  }
  t <- rep_len(t, prep$nobs)
  if (!sort && length(prep$old_order)) t <- t[order(prep$old_order)]
  out <- matrix(NA_real_, prep$ndraws, prep$nobs)
  for (s in seq_along(prep$bhaz$sbhaz)) {
    take <- which(prep$bhaz$Jgrbhaz == s)
    if (!length(take)) next
    out[, take] <- bhaz_at(
      t[take], prep$bhaz$sbhaz[[s]], prep$bhaz$basis, prep$bhaz$tail,
      integrate = integrate, paired = FALSE
    ) * prep$dpars$mu[, take, drop = FALSE]
  }
  reorder_obs(out, prep$old_order, sort = sort)
}

#' @rdname cox_survival
#' @export
cox_hazard <- function(t, object, sort = FALSE, ...) {
  cox_quantity(t, prepare_cox_predictions(object, ...), FALSE, sort)
}

#' @rdname cox_survival
#' @export
cox_cumhaz <- function(t, object, sort = FALSE, ...) {
  cox_quantity(t, prepare_cox_predictions(object, ...), TRUE, sort)
}

#' @rdname cox_survival
#' @export
cox_surv <- function(t, object, sort = FALSE, ...) {
  exp(-cox_cumhaz(t, object, sort = sort, ...))
}

#' @rdname cox_survival
#' @param times A common vector of finite nonnegative times.
#' @param probs Posterior quantiles for streaming survival summaries.
#' @return \code{cox_surv_curve} returns a draws by times by observations
#'   array; \code{cox_surv_summary} returns a times by quantiles by observations
#'   array without materializing the full curve array.
#' @export
cox_surv_curve <- function(times, object, sort = FALSE, ...) {
  prep <- prepare_cox_predictions(object, ...)
  cox_curve(times, prep, sort = sort)
}

#' @rdname cox_survival
#' @export
cox_surv_summary <- function(times, object, probs = c(0.025, 0.5, 0.975),
                              sort = FALSE, ...) {
  if (!length(probs) || anyNA(probs) || any(probs < 0 | probs > 1)) {
    stop2("'probs' must be a non-empty vector of probabilities.")
  }
  prep <- prepare_cox_predictions(object, ...)
  cox_curve(times, prep, probs = probs, sort = sort)
}

cox_curve <- function(times, prep, probs = NULL, sort = FALSE) {
  if (!length(times) || anyNA(times) || any(!is.finite(times) | times < 0)) {
    stop2("'times' must be a non-empty vector of finite nonnegative times.")
  }
  if (is.null(probs)) {
    out <- array(NA_real_, c(prep$ndraws, length(times), prep$nobs))
  } else {
    out <- array(NA_real_, c(length(probs), length(times), prep$nobs))
  }
  for (s in seq_along(prep$bhaz$sbhaz)) {
    H0 <- bhaz_at(times, prep$bhaz$sbhaz[[s]], prep$bhaz$basis,
                  prep$bhaz$tail, integrate = TRUE, paired = FALSE)
    for (i in which(prep$bhaz$Jgrbhaz == s)) {
      S <- exp(-H0 * prep$dpars$mu[, i])
      out[, , i] <- if (is.null(probs)) S else
        t(matrixStats::colQuantiles(S, probs = probs, drop = FALSE))
    }
  }
  # reorder_obs expects observations in the second dimension.
  out <- aperm(reorder_obs(aperm(out, c(1, 3, 2)), prep$old_order, sort),
               c(1, 3, 2))
  if (!is.null(probs)) out <- aperm(out, c(2, 1, 3))
  out
}

#' Apply Censoring to Cox Predictive Draws
#'
#' @param yrep A draws by observations matrix of sign-encoded predictions.
#' @param cens_times Positive censoring times: a common value, one per
#'   observation, a matrix matching \code{yrep}, or \code{"reverse_km"}
#'   to simulate them
#'   using \code{\link{cox_censor_times}}. \code{Inf} means no censoring.
#' @param y,status Observed times and binary event indicators, required for
#'   reverse Kaplan-Meier simulation.
#' @return A matrix with right-censored times encoded as negative values.
#' @details Censoring is applied explicitly, not inferred automatically by
#'   \code{posterior_predict}. Plotting these censored replications requires
#'   an external plotting function that understands their event indicators.
#' @export
cox_censor <- function(yrep, cens_times, y = NULL, status = NULL) {
  if (!is.matrix(yrep) || !is.numeric(yrep) || anyNA(yrep)) {
    stop2("'yrep' must be a numeric matrix without missing values.")
  }
  if (is.character(cens_times)) {
    match.arg(cens_times, "reverse_km")
    ct <- matrix(cox_censor_times(y, status, length(yrep)), nrow(yrep))
  } else if (is.matrix(cens_times)) {
    if (!identical(dim(cens_times), dim(yrep))) {
      stop2("'cens_times' must have the same dimensions as 'yrep'.")
    }
    ct <- cens_times
  } else {
    if (!length(cens_times) %in% c(1L, ncol(yrep))) {
      stop2("'cens_times' must have length 1 or the number of observations.")
    }
    ct <- matrix(cens_times, nrow(yrep), ncol(yrep), byrow = TRUE)
  }
  if (!is.numeric(ct) || anyNA(ct) || any(ct <= 0)) {
    stop2("Censoring times must be positive, with Inf for no censoring.")
  }
  time <- pmin(abs(yrep), ct)
  ifelse(yrep < 0 | abs(yrep) > ct, -time, time)
}

#' Simulate Censoring Times from a Reverse Kaplan-Meier Estimate
#'
#' @param y Positive finite observed times.
#' @param status Binary event indicators (1 for an event, 0 for censoring).
#' @param n Number of censoring times to simulate.
#' @return Positive censoring times, with \code{Inf} for draws beyond the
#'   observed censoring distribution.
#' @details Uses the optional \pkg{prodlim} package. No KM plotting
#'   implementation is provided by brms.
#' @export
cox_censor_times <- function(y, status, n = length(y)) {
  y <- as.numeric(y)
  status <- as.numeric(status)
  n <- as_one_integer(n)
  if (!length(y) || length(y) != length(status) || anyNA(status) ||
      any(!status %in% 0:1) || any(!is.finite(y) | y <= 0) || n < 0L) {
    stop2("Supply positive finite 'y', matching binary 'status', and n >= 0.")
  }
  if (all(status == 1)) return(rep(Inf, n))
  require_package("prodlim")
  dat <- data.frame(time = y, status = status)
  fit <- prodlim::prodlim(prodlim::Hist(time, status) ~ 1, data = dat,
                          reverse = TRUE)
  jumps <- sort(unique(y))
  G <- as.vector(predict(fit, times = jumps))
  if (anyNA(G)) {
    stop2("The reverse Kaplan-Meier estimate is undefined at observed times.")
  }
  index <- findInterval(-runif(n), -G) + 1L
  out <- jumps[pmin(index, length(jumps))]
  out[index > length(jumps)] <- Inf
  out
}

randomize_cox_pit <- function(p, cens, upper = NULL) {
  if (is.null(cens)) return(p)
  if (anyNA(cens) || any(!cens %in% c(-1, 0, 1, 2))) {
    stop2("Invalid Cox censoring indicators.")
  }
  cm <- data2draws(cens, dim(p))
  lo <- hi <- p
  lo[cm == -1] <- 0
  hi[cm == 1] <- 1
  if (any(cm == 2)) {
    if (is.null(upper)) {
      stop2("Interval-censored PIT requires upper endpoints 'rcens'.")
    }
    hi[cm == 2] <- upper[cm == 2]
  }
  if (anyNA(lo) || anyNA(hi) || any(lo > hi)) {
    stop2("The censoring interval has undefined or invalid CDF endpoints.")
  }
  lo + data2draws(runif(ncol(p)), dim(p)) * (hi - lo)
}

#' Censoring-Aware Cox Calibration Diagnostics
#'
#' @inheritParams cox_survival
#' @param randomize_cens Randomize CDFs within observed censoring intervals.
#' @param average Aggregate over posterior draws before randomization.
#'   If \code{FALSE}, return draws by observations conditional values.
#' @return \code{cox_pit} returns one PIT value per observation if
#'   \code{average = TRUE}, otherwise a draws by observations matrix.
#'   \code{cox_snell} returns a draws by observations residual matrix.
#' @details These are numerical diagnostics, not plotting functions.
#'   Unlike \code{posterior_predict(output = "pit")}, censoring-aware PIT
#'   randomizes right-, left-, and interval-censored observations within
#'   their possible CDF intervals. Randomization uses one uniform per
#'   observation, shared across draws; it must not be averaged away.
#'   Cox-Snell residuals subtract the cumulative hazard at a left-truncation
#'   bound. They are not provided for finite upper truncation.
#' @name cox_diagnostics
NULL

#' @rdname cox_diagnostics
#' @export
cox_pit <- function(object, randomize_cens = TRUE, average = TRUE,
                     sort = FALSE, ...) {
  prep <- prepare_cox_predictions(object, ..., check_response = TRUE)
  randomize_cens <- as_one_logical(randomize_cens)
  average <- as_one_logical(average)
  p <- posterior_predict(prep, output = "probability", sort = TRUE)
  cens <- if (randomize_cens) prep$data$cens else NULL
  upper <- NULL
  if (any(cens == 2)) {
    if (is.null(prep$data$rcens)) {
      stop2("Interval-censored PIT requires upper endpoints 'rcens'.")
    }
    upper <- vapply(seq_len(prep$nobs), function(i) {
      posterior_predict_cox(i, prep, output = "probability",
                            q = prep$data$rcens[i])
    }, numeric(prep$ndraws))
    dim(upper) <- c(prep$ndraws, prep$nobs)
  }
  if (average) {
    p <- matrix(colMeans(p), nrow = 1L)
    if (!is.null(upper)) upper <- matrix(colMeans(upper), nrow = 1L)
  }
  p <- randomize_cox_pit(p, cens, upper)
  p <- reorder_obs(p, prep$old_order, sort)
  if (average) as.vector(p) else p
}

#' @rdname cox_diagnostics
#' @export
cox_snell <- function(object, sort = FALSE, ...) {
  prep <- prepare_cox_predictions(object, ..., check_response = TRUE)
  if (any(is.finite(prep$data$ub))) {
    stop2("Cox-Snell residuals are not provided for finite upper truncation.")
  }
  out <- matrix(NA_real_, prep$ndraws, prep$nobs)
  for (i in seq_len(prep$nobs)) {
    mu <- prep$dpars$mu[, i]
    sb <- prep$bhaz$sbhaz[[prep$bhaz$Jgrbhaz[i]]]
    tr <- cox_truncation(mu, sb, prep$bhaz$basis, prep$bhaz$tail,
                         prep$data$lb[i], prep$data$ub[i])
    out[, i] <- pmax(mu * bhaz_at(prep$data$Y[i], sb, prep$bhaz$basis,
                                 prep$bhaz$tail, integrate = TRUE) - tr$Hlo, 0)
  }
  reorder_obs(out, prep$old_order, sort)
}

#' @rdname cox_diagnostics
#' @param loo_obj Optional \code{loo} result computed with
#'   \code{save_psis = TRUE} from the same fit, response, and all its draws
#'   in their original order. Subsets and reordered draws require recomputing
#'   PSIS rather than reusing incompatible weights. Precomputed weights
#'   require the original data: do not supply \code{newdata}.
#' @param ndraws,draw_ids See \code{\link{prepare_predictions}}.
#' @param newdata,resp See \code{\link{prepare_predictions}}.
#' @details \code{loo_pit_cox} returns PSIS-weighted CDFs randomized after
#'   weighting. Truncated and interval-censored fits are rejected because
#'   the existing Cox likelihood path does not evaluate arbitrary bounds.
#'   When computing weights, the likelihood and predictive CDFs use one
#'   shared preparation, including any simulated new-group effects.
#' @export
loo_pit_cox <- function(object, loo_obj = NULL, randomize_cens = TRUE,
                         ndraws = NULL, draw_ids = NULL,
                         newdata = NULL, resp = NULL, ...) {
  if (!is.brmsfit(object)) {
    stop2("'loo_pit_cox' requires a fitted brms model.")
  }
  if (!is.null(loo_obj) && !is.null(newdata)) {
    stop2("Precomputed PSIS weights require the original data; ",
          "do not supply 'newdata'.")
  }
  ids <- validate_draw_ids(object, draw_ids, ndraws) %||%
    seq_len(posterior::ndraws(object))
  prep <- prepare_cox_predictions(
    object, newdata = newdata, resp = resp, draw_ids = ids, ...,
    check_response = TRUE
  )
  if (is_trunc(prep) || any(prep$data$cens == 2)) {
    stop2("Cox LOO-PIT is not supported for truncated or interval-censored ",
          "models by the existing likelihood implementation.")
  }
  if (is.null(loo_obj)) {
    ll <- log_lik(prep, ...)
    r_eff <- r_eff_log_lik(ll[order(ids), , drop = FALSE], fit = object)
    loo_obj <- loo::loo(ll, r_eff = r_eff, save_psis = TRUE, ...)
  } else {
    if (!identical(ids, seq_len(posterior::ndraws(object)))) {
      stop2("Recompute PSIS for draw subsets or reordered draws; ",
            "do not supply 'loo_obj'.")
    }
    if (!identical(attr(loo_obj, "yhash"), hash_response(object, resp = resp))) {
      stop2("'loo_obj' does not match the requested response and observations.")
    }
  }
  if (is.null(loo_obj$psis_object)) {
    stop2("'loo_obj' must be computed with save_psis = TRUE.")
  }
  w <- weights(loo_obj$psis_object, normalize = TRUE, log = FALSE)
  p <- posterior_predict(prep, output = "probability", sort = FALSE)
  if (!identical(dim(w), dim(p))) {
    stop2("PSIS weights do not match the requested response and draws.")
  }
  p <- matrix(colSums(w * p), nrow = 1L)
  cens <- if (as_one_logical(randomize_cens) && !is.null(prep$data$cens)) {
    as.vector(reorder_obs(matrix(prep$data$cens, nrow = 1L), prep$old_order))
  } else {
    NULL
  }
  as.vector(randomize_cox_pit(p, cens))
}

#' Decode Censored Cox Predictions
#'
#' @param x A numeric vector or matrix of Cox predictive draws.
#' @return A list with \code{time = abs(x)} and event indicators
#'   \code{status = as.integer(x > 0)}, preserving the dimensions of \code{x}.
#' @details With \code{tail = "none"}, negative predictions indicate
#'   right censoring at the upper boundary of the fitted baseline spline.
#' @export
cox_decode <- function(x) {
  list(time = abs(x), status = (x > 0) * 1L)
}
