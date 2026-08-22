# Helpers for the physical sum-to-zero parameterization of group effects

# Does a local linear predictor contain a sum-to-zero group-level block?
has_re_s2z <- function(x) {
  is.bframel(x) && has_rows(x$frame$re) &&
    "s2z" %in% names(x$frame$re) && any(x$frame$re$s2z)
}

# Return validated centering fractions for one S2Z block. Logical values from
# reframes created before partial centering was added remain valid endpoints.
re_s2z_center_values <- function(r) {
  stopifnot(is.reframe(r), has_rows(r))
  value <- r[["s2z_center"]]
  if (is.null(value)) {
    return(rep(1, nrow(r)))
  }
  if (!is.numeric(value) && !is.logical(value)) {
    stop2("Internal error: invalid sum-to-zero centering fractions.")
  }
  value <- as.numeric(value)
  if (length(value) != nrow(r) || anyNA(value) ||
      any(!is.finite(value)) || any(value < 0 | value > 1)) {
    stop2("Internal error: sum-to-zero centering fractions must be in [0, 1].")
  }
  value
}

# Return the per-coefficient auto flags, defaulting to FALSE for old reframes.
re_s2z_center_auto <- function(r) {
  stopifnot(is.reframe(r), has_rows(r))
  value <- r[["s2z_center_auto"]]
  if (is.null(value)) {
    return(rep(FALSE, nrow(r)))
  }
  if (!is.logical(value) || length(value) != nrow(r) || anyNA(value)) {
    stop2("Internal error: invalid automatic sum-to-zero centering flags.")
  }
  value
}

# Classify an S2Z block without changing the legacy endpoint code paths.
re_s2z_center_mode <- function(r) {
  rho <- re_s2z_center_values(r)
  if (any(re_s2z_center_auto(r)) || any(rho > 0 & rho < 1) ||
      length(unique(rho)) > 1L) {
    return("partial")
  }
  if (all(rho == 1)) "centered" else "noncentered"
}

# Return all local linear predictor frames contained in a brms frame.
all_bframel <- function(x) {
  if (is.bframel(x)) {
    return(list(x))
  }
  if (is.mvbrmsframe(x)) {
    return(unlist(lapply(x$terms, all_bframel), recursive = FALSE))
  }
  if (is.brmsframe(x)) {
    return(unlist(lapply(c(x$dpars, x$nlpars), all_bframel), recursive = FALSE))
  }
  list()
}

# Find the local linear predictor frame belonging to a group-level block.
re_s2z_bframel <- function(bframe, r) {
  stopifnot(is.anybrmsframe(bframe), is.reframe(r), has_rows(r))
  rp <- unique(combine_prefix(check_prefix(r)))
  out <- Filter(
    function(x) identical(combine_prefix(check_prefix(x)), rp),
    all_bframel(bframe)
  )
  if (length(out) != 1L) {
    stop2("Internal error while matching a sum-to-zero group-level block ",
          "to its linear predictor.")
  }
  out[[1]]
}

# Effective scalar population-level prior for one coefficient.
re_s2z_prior <- function(prior, bframe, class, coef = "") {
  stopifnot(is.brmsprior(prior), is.bframel(bframe))
  px <- check_prefix(bframe)
  p <- subset2(prior, class = class, ls = px)
  if (class == "b" && nzchar(coef)) {
    pcoef <- subset2(p, coef = coef)
    if (nrow(pcoef) == 1L && nzchar(pcoef$prior)) {
      ans <- pcoef
    } else {
      ans <- subset2(p, coef = "")
    }
  } else {
    ans <- subset2(p, coef = "")
  }
  if (!nrow(ans)) {
    return(list(dist = "flat", location = 0, scale = 1, df = NA_real_))
  }
  if (nrow(ans) > 1L) {
    # Follow stan_prior() and prefer the most specific populated base row.
    take <- nzchar(ans$prior) | nzchar(ans$lb) | nzchar(ans$ub) |
      nzchar(ans$tag)
    ans <- ans[take, , drop = FALSE]
    if (nrow(ans) > 1L) {
      specificity <- rowSums(nzchar(ans[, vars_prefix(), drop = FALSE]))
      ans <- ans[which.max(specificity), , drop = FALSE]
    }
  }
  if (nrow(ans) && (nzchar(ans$lb) || nzchar(ans$ub))) {
    stop2("The sum-to-zero parameterization does not yet ",
          "support bounded population-level priors (coefficient '", coef,
          "').")
  }
  if (nrow(ans) && nzchar(ans$tag)) {
    stop2("The sum-to-zero parameterization does not yet ",
          "support tagged population-level priors (coefficient '", coef,
          "').")
  }
  value <- if (nrow(ans)) ans$prior else ""
  parse_re_s2z_prior(value, coef = coef)
}

# Parse the conditionally Gaussian scalar priors supported by the fast path.
parse_re_s2z_prior <- function(prior, coef = "") {
  prior <- gsub("[[:space:]]+", "", prior)
  if (!nzchar(prior)) {
    return(list(dist = "flat", location = 0, scale = 1, df = NA_real_))
  }
  call <- try(str2lang(prior), silent = TRUE)
  if (inherits(call, "try-error") || !is.call(call)) {
    stop2("Prior '", prior, "' is not supported by the ",
          "sum-to-zero parameterization (coefficient '", coef, "').")
  }
  dist <- as.character(call[[1]])
  args <- as.list(call[-1])
  number <- function(x) {
    out <- suppressWarnings(as.numeric(deparse0(x)))
    if (length(out) != 1L || !is.finite(out)) {
      stop2("All arguments of population-level priors used with the ",
            "sum-to-zero parameterization must currently be ",
            "numeric constants (coefficient '", coef, "').")
    }
    out
  }
  if (dist == "std_normal" && !length(args)) {
    out <- list(dist = "normal", location = 0, scale = 1, df = NA_real_)
  } else if (dist == "normal" && length(args) == 2L) {
    out <- list(
      dist = "normal", location = number(args[[1]]),
      scale = number(args[[2]]), df = NA_real_
    )
  } else if (dist == "student_t" && length(args) == 3L) {
    out <- list(
      dist = "student", df = number(args[[1]]),
      location = number(args[[2]]), scale = number(args[[3]])
    )
  } else if (dist == "cauchy" && length(args) == 2L) {
    out <- list(
      dist = "student", df = 1, location = number(args[[1]]),
      scale = number(args[[2]])
    )
  } else {
    stop2(
      "Prior '", prior, "' is not supported by the ",
      "sum-to-zero parameterization. Supported population-level priors ",
      "are flat, normal, student_t, and cauchy (coefficient '", coef, "')."
    )
  }
  if (!(out$scale > 0) || out$dist == "student" && !(out$df > 0)) {
    stop2("Scale and degrees-of-freedom arguments must be positive in ",
          "population-level priors used with the sum-to-zero ",
          "parameterization (coefficient '", coef, "').")
  }
  out
}

# Check that an estimated scale parameter retains its required support.
validate_re_s2z_scale_bounds <- function(prior, class) {
  stopifnot(is.brmsprior(prior), length(class) == 1L)
  bounds <- stan_base_prior(prior, col = c("lb", "ub"))
  lb <- bounds[["lb"]][1L]
  ub <- bounds[["ub"]][1L]
  numeric_bound <- function(x) {
    suppressWarnings(as.numeric(x))
  }
  lb_num <- numeric_bound(lb)
  if (!nzchar(lb) || length(lb_num) != 1L || !is.finite(lb_num) ||
      lb_num < 0) {
    stop2("Class '", class, "' must have a finite non-negative lower ",
          "bound for gr(..., s2z = TRUE).")
  }
  if (nzchar(ub)) {
    ub_num <- numeric_bound(ub)
    if (length(ub_num) != 1L || !is.finite(ub_num) || ub_num <= 0 ||
        ub_num <= lb_num) {
      stop2("Class '", class, "' must have a finite positive upper bound ",
            "greater than its lower bound when one is specified for ",
            "gr(..., s2z = TRUE).")
    }
  }
  invisible(NULL)
}

# Check fixed group scales before the Gaussian precision is constructed.
validate_re_s2z_sd_prior <- function(prior, r) {
  stopifnot(is.brmsprior(prior), is.reframe(r), has_rows(r))
  px <- check_prefix(r)
  p <- subset2(
    prior, class = "sd", coef = c(r$coef, ""),
    group = c(r$group[1], ""), ls = px
  )
  validate_re_s2z_scale_bounds(p, class = "sd")
  base_prior <- stan_base_prior(p)
  for (coef in r$coef) {
    pcoef <- subset2(p, coef = coef)
    coef_prior <- pcoef$prior[nzchar(pcoef$prior)]
    stopifnot(length(coef_prior) <= 1L)
    value <- if (length(coef_prior)) {
      coef_prior[[1]]
    } else {
      base_prior
    }
    if (!stan_is_constant_prior(value)) {
      next
    }
    call <- try(str2lang(value), silent = TRUE)
    fixed <- if (inherits(call, "try-error") || length(call) < 2L) {
      NA_real_
    } else {
      suppressWarnings(as.numeric(deparse0(call[[2]])))
    }
    if (length(fixed) != 1L || !is.finite(fixed) || fixed <= 0) {
      stop2("Group-level standard deviations fixed with 'constant' must ",
            "be positive numeric scalars for gr(..., s2z = TRUE) ",
            "(coefficient '", coef, "').")
    }
  }
  invisible(NULL)
}

# Check fixed log-scale heterogeneity parameters. Zero is the exact
# shared-scale limit, so it is intentionally allowed here.
validate_re_s2z_sdlog_prior <- function(prior, r) {
  stopifnot(is.brmsprior(prior), is.reframe(r), has_rows(r))
  if (!identical(r$scale[1], "varying")) {
    return(invisible(NULL))
  }
  px <- check_prefix(r)
  p <- subset2(
    prior, class = "sdlog", coef = c(r$coef, ""),
    group = c(r$group[1], ""), ls = px
  )
  validate_re_s2z_scale_bounds(p, class = "sdlog")
  base_prior <- stan_base_prior(p)
  for (coef in r$coef) {
    pcoef <- subset2(p, coef = coef)
    coef_prior <- pcoef$prior[nzchar(pcoef$prior)]
    stopifnot(length(coef_prior) <= 1L)
    value <- if (length(coef_prior)) coef_prior[[1]] else base_prior
    if (!stan_is_constant_prior(value)) {
      next
    }
    call <- try(str2lang(value), silent = TRUE)
    fixed <- if (inherits(call, "try-error") || length(call) < 2L) {
      NA_real_
    } else {
      suppressWarnings(as.numeric(deparse0(call[[2]])))
    }
    if (length(fixed) != 1L || !is.finite(fixed) || fixed < 0) {
      stop2("Log-scale standard deviations fixed with 'constant' must ",
            "be non-negative numeric scalars for gr(..., scale = ",
            "\"varying\") (coefficient '", coef, "').")
    }
  }
  invisible(NULL)
}

# Describe one local S2Z block, including the fixed/group design mapping.
re_s2z_info <- function(bframe, prior = NULL, id = NULL) {
  stopifnot(is.bframel(bframe))
  if (!has_re_s2z(bframe)) {
    return(NULL)
  }
  r <- bframe$frame$re[bframe$frame$re$s2z, , drop = FALSE]
  ids <- unique(r$id)
  if (is.null(id) && length(ids) != 1L) {
    stop2("An explicit group-level ID is required when describing multiple ",
          "sum-to-zero blocks.")
  }
  if (!is.null(id)) {
    if (length(id) != 1L || !id %in% ids) {
      stop2("Invalid sum-to-zero group-level ID.")
    }
    ids <- id
  }
  r_id <- subset2(bframe$frame$re, id = ids)
  if (!all(r_id$s2z)) {
    stop2("All coefficients sharing a group-level ID must use the same ",
          "'s2z' setting.")
  }
  re_s2z_center_values(r_id)
  re_s2z_center_auto(r_id)
  r <- r_id
  qnames <- bframe$frame$fe$vars
  match_q <- match(r$coef, qnames)
  if (anyNA(match_q)) {
    missing <- collapse_comma(r$coef[is.na(match_q)])
    stop2("Every sum-to-zero group-level coefficient must have a matching ",
          "population-level design column. Missing: ", missing, ".")
  }
  out <- nlist(
    id = ids, r, qnames, match_q,
    center = bframe$frame$fe$center,
    fixef = bframe$frame$fe$vars_stan,
    p = usc(combine_prefix(check_prefix(bframe)))
  )
  if (!is.null(prior)) {
    specs <- vector("list", length(qnames))
    for (i in seq_along(qnames)) {
      if (out$center && qnames[i] == "Intercept") {
        specs[[i]] <- re_s2z_prior(prior, bframe, class = "Intercept")
      } else {
        specs[[i]] <- re_s2z_prior(
          prior, bframe, class = "b", coef = qnames[i]
        )
      }
    }
    out$prior <- specs
  }
  out
}

# Describe all local S2Z blocks in one linear predictor. Blocks retain their
# own covariance models but share one joint omitted-mean system.
re_s2z_infos <- function(bframe, prior = NULL) {
  stopifnot(is.bframel(bframe))
  if (!has_re_s2z(bframe)) {
    return(list())
  }
  ids <- unique(bframe$frame$re$id[bframe$frame$re$s2z])
  lapply(ids, function(id) re_s2z_info(bframe, prior = prior, id = id))
}

# Construct an S2Z design matrix in covariance-block coefficient order. A
# block may be assembled from multiple grouping terms that share an ID.
re_s2z_design_matrix <- function(bframe, data, id = NULL) {
  stopifnot(is.bframel(bframe))
  info <- re_s2z_info(bframe, id = id)
  if (is.null(info)) {
    return(matrix(numeric(), nrow = nrow(data), ncol = 0L))
  }
  r <- info$r
  out <- matrix(NA_real_, nrow = nrow(data), ncol = nrow(r))
  for (gn in unique(r$gn)) {
    take <- which(r$gn == gn)
    rg <- r[take, , drop = FALSE]
    Zg <- get_model_matrix(rg$form[[1]], data = data, rename = FALSE)
    if (ncol(Zg) != length(take)) {
      stop2("Internal mismatch in the sum-to-zero group-level design matrix.")
    }
    out[, take] <- Zg
  }
  colnames(out) <- r$coef
  out
}

# Data-driven centering fractions for one S2Z block. Each column is normalized
# by its global nonzero RMS. Within a group, pivoted QR removes the span of all
# other varying columns before squared residual exposure is accumulated. This
# makes the rule invariant to coefficient rescaling and conservative for
# locally unidentified slopes and interactions.
re_s2z_auto_rho <- function(bframe, data, threshold = 25, id = NULL) {
  stopifnot(is.bframel(bframe), is.numeric(threshold), length(threshold) == 1L,
            is.finite(threshold), threshold > 0)
  info <- re_s2z_info(bframe, id = id)
  if (is.null(info)) {
    return(matrix(numeric(), nrow = 0L, ncol = 0L))
  }
  r <- info$r
  Z <- re_s2z_design_matrix(bframe, data = data, id = info$id)
  M <- ncol(Z)
  levels <- get_levels(r)[[r$group[1]]]
  J <- length(levels)
  group_var <- r$gcall[[1]]$groups
  if (length(group_var) != 1L) {
    stop2("Automatic sum-to-zero centering requires one grouping variable.")
  }
  group_index <- match(get(group_var, data), levels)
  if (anyNA(group_index)) {
    stop2("Internal mismatch while computing automatic sum-to-zero centering.")
  }

  nonzero_rms <- vapply(seq_len(M), function(k) {
    zk <- Z[, k]
    keep <- is.finite(zk) & zk != 0
    if (!any(keep)) 1 else sqrt(mean(zk[keep]^2))
  }, numeric(1))
  Z <- sweep(Z, 2L, nonzero_rms, "/")
  exposure <- matrix(0, nrow = J, ncol = M)
  qr_tol <- sqrt(.Machine$double.eps)
  for (j in seq_len(J)) {
    rows <- which(group_index == j)
    if (!length(rows)) {
      next
    }
    Zj <- Z[rows, , drop = FALSE]
    for (k in seq_len(M)) {
      residual <- Zj[, k]
      if (M > 1L) {
        others <- Zj[, -k, drop = FALSE]
        keep <- colSums(abs(others)) > 0
        others <- others[, keep, drop = FALSE]
        if (ncol(others)) {
          residual <- lm.fit(
            x = others, y = residual, tol = qr_tol
          )$residuals
        }
      }
      value <- sum(residual^2)
      zero_tol <- qr_tol * max(1, sum(Zj[, k]^2))
      if (value < zero_tol) {
        value <- 0
      }
      exposure[j, k] <- value
    }
  }
  rho <- exposure / (exposure + threshold)
  colnames(rho) <- r$coef
  rownames(rho) <- levels
  rho
}

# Validate S2Z structure after formulas, data, and priors have been resolved.
validate_re_s2z <- function(bframe, prior) {
  stopifnot(is.anybrmsframe(bframe), is.brmsprior(prior))
  all_frames <- all_bframel(bframe)
  for (x in all_frames) {
    r <- x$frame$re
    varying <- rep(FALSE, nrow(r))
    if (has_rows(r) && "scale" %in% names(r)) {
      varying <- r$scale == "varying"
    }
    if (any(varying & !r$s2z)) {
      stop2("Group-varying scales currently require ",
            "gr(..., s2z = TRUE, scale = \"varying\").")
    }
  }
  # IDs may couple ordinary group effects across linear predictors. Such a
  # coupling is not compatible with the local S2Z omitted-mean systems, even
  # if only one occurrence of the shared ID requests S2Z. Check all terms
  # before filtering to S2Z-containing predictors so mixed conventional/S2Z
  # uses fail with the same clear error in either predictor order.
  s2z_ids <- unique(unlist(lapply(all_frames, function(x) {
    r <- x$frame$re
    if (!has_rows(r)) {
      return(integer())
    }
    r$id[r$s2z]
  }), use.names = FALSE))
  for (id in s2z_ids) {
    n_frames <- sum(vapply(all_frames, function(x) {
      r <- x$frame$re
      has_rows(r) && id %in% r$id
    }, logical(1)))
    if (n_frames > 1L) {
      stop2("A sum-to-zero group-level ID cannot span multiple ",
            "linear predictors.")
    }
  }
  frames <- Filter(has_re_s2z, all_frames)
  if (!length(frames)) {
    return(invisible(NULL))
  }
  if (get_sample_prior(prior) != "no") {
    stop2("Argument 'sample_prior' is not yet supported together with ",
          "gr(..., s2z = TRUE).")
  }
  ids <- integer()
  for (x in frames) {
    infos <- re_s2z_infos(x, prior = prior)
    for (info in infos) {
      r <- info$r
      if (info$id %in% ids) {
        stop2("A sum-to-zero group-level ID cannot span multiple ",
              "linear predictors.")
      }
      ids <- c(ids, info$id)
      if (length(unique(r$group)) != 1L) {
        stop2("All coefficients sharing a sum-to-zero group-level ID must ",
              "use the same grouping factor.")
      }
      if (length(unique(r$scale)) != 1L) {
        stop2("All coefficients sharing a group-level ID must use the same ",
              "'scale' setting.")
      }
      if (length(unique(r$dist)) != 1L) {
        stop2("All coefficients sharing a group-level ID must use the same ",
              "group-level distribution.")
      }
      if (length(unique(r$cor)) != 1L) {
        stop2("All coefficients sharing a group-level ID must use the same ",
              "'cor' setting.")
      }
      if (any(nzchar(r$gtype)) || any(nzchar(r$type))) {
        stop2("The sum-to-zero parameterization currently ",
              "supports only ordinary 'gr' terms.")
      }
      has_pw <- any(vapply(r$gcall, function(gcall) {
        isTRUE(nzchar(gcall$pw))
      }, logical(1)))
      if (any(nzchar(r$by), na.rm = TRUE) ||
          any(nzchar(r$cov), na.rm = TRUE) || has_pw) {
        stop2("Arguments 'by', 'cov', and 'pw' are not yet supported ",
              "together with gr(..., s2z = TRUE).")
      }
      if (length(get_levels(r)[[r$group[1]]]) < 2L) {
        stop2("The sum-to-zero parameterization requires at ",
              "least two observed grouping levels.")
      }
      if (is_ordinal(x$family)) {
        stop2("Ordinal thresholds are not yet supported together with ",
              "gr(..., s2z = TRUE).")
      }
      if (order_intercepts(x)) {
        stop2("Ordered mixture intercepts are not yet supported together ",
              "with gr(..., s2z = TRUE).")
      }
      if (x$frame$fe$sparse || x$frame$fe$decomp != "none") {
        stop2("Sparse and QR-decomposed population-level design matrices ",
              "are not yet supported together with gr(..., s2z = TRUE).")
      }
      if (has_special_prior(prior, check_prefix(r), class = "sd") ||
          has_special_prior(prior, check_prefix(r), class = "sdlog") ||
          has_special_prior(prior, x, class = "b")) {
        stop2("Special priors are not yet supported together with ",
              "gr(..., s2z = TRUE).")
      }
      validate_re_s2z_sd_prior(prior, r)
      validate_re_s2z_sdlog_prior(prior, r)
    }
  }
  invisible(NULL)
}

# Validate that matched fixed and group columns are numerically identical.
validate_re_s2z_design <- function(bframe, data) {
  stopifnot(is.bframel(bframe))
  if (!has_re_s2z(bframe)) {
    return(invisible(NULL))
  }
  X <- bframe$sdata$fe$X
  for (info in re_s2z_infos(bframe)) {
    Z <- re_s2z_design_matrix(bframe, data = data, id = info$id)
    if (ncol(Z) != nrow(info$r)) {
      stop2("Internal mismatch in the sum-to-zero group-level design matrix.")
    }
    for (j in seq_len(ncol(Z))) {
      same <- isTRUE(all.equal(
        unname(Z[, j]), unname(X[, info$match_q[j]]),
        tolerance = sqrt(.Machine$double.eps), check.attributes = FALSE
      ))
      if (!same) {
        stop2("Population- and group-level design columns must be identical ",
              "for gr(..., s2z = TRUE). Mismatch for coefficient '",
              info$r$coef[j], "'.")
      }
    }
  }
  invisible(NULL)
}
