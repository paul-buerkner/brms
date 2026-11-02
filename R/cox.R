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
