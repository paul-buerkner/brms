# Helpers for the physical sum-to-zero parameterization of group effects

# Does a local linear predictor contain a sum-to-zero group-level block?
has_re_s2z <- function(x) {
  is.bframel(x) && has_rows(x$frame$re) &&
    "s2z" %in% names(x$frame$re) && any(x$frame$re$s2z)
}

# Does the ordinary location component of a brmsterms object contain an S2Z
# group term? Likelihood generation works from terms rather than frames, where
# the flag still lives in the parsed gr() calls.
has_re_s2z_terms <- function(x) {
  if (!is.brmsterms(x)) {
    return(FALSE)
  }
  mu <- x$dpars[["mu"]]
  re <- mu$re %||% NULL
  has_rows(re) && any(vapply(re$gcall, function(gcall) {
    isTRUE(gcall$s2z)
  }, logical(1)))
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

# A neutral placeholder used for fixed-only coordinates. Keeping the prior
# list full-length makes the descriptor backwards compatible with Stan code
# that still indexes it in population-coordinate order.
re_s2z_flat_prior <- function() {
  list(dist = "flat", location = 0, scale = 1, df = NA_real_)
}

# Collect the model coordinates needed to make an S2Z diagnostic actionable.
re_s2z_context <- function(bframe, r = NULL, coef = NULL, prior = NULL) {
  stopifnot(is.bframel(bframe))
  px <- check_prefix(bframe, keep_mu = TRUE)
  response <- px$resp %||% ""
  if ((!length(response) || !nzchar(response)) &&
      is.formula(bframe$respform)) {
    response_vars <- all.vars(bframe$respform)
    if (length(response_vars)) {
      response <- response_vars[1L]
    }
  }
  if ((!length(response) || !nzchar(response)) &&
      nzchar(bframe$frame$s2z_context$response %||% "")) {
    response <- bframe$frame$s2z_context$response
  }
  if (!length(response) || !nzchar(response)) {
    response <- "<unspecified>"
  }
  family <- bframe$family$family %||% ""
  if (!length(family) || !nzchar(family)) {
    family <- bframe$frame$s2z_context$family %||% ""
  }
  if (!length(family) || !nzchar(family)) {
    family <- "<unspecified>"
  }
  dpar <- px$dpar %||% ""
  nlpar <- px$nlpar %||% ""
  if (!nzchar(dpar)) {
    dpar <- "mu"
  }
  if (!nzchar(nlpar)) {
    nlpar <- "<none>"
  }
  group <- id <- "<unspecified>"
  if (!is.null(r) && has_rows(r)) {
    groups <- unique(r$group[nzchar(r$group)])
    if (length(groups)) {
      group <- paste(groups, collapse = ", ")
    }
    public_id <- r$gcall[[1L]]$id %||% NA
    if (length(public_id) == 1L && !is.na(public_id)) {
      public_id <- trimws(as.character(public_id))
    }
    if (length(public_id) == 1L && !is.na(public_id) &&
        nzchar(public_id)) {
      id <- public_id
    } else {
      ids <- unique(r$id)
      if (length(ids)) {
        id <- paste(ids, collapse = ", ")
      }
    }
    if (is.null(coef)) {
      coef <- unique(r$coef)
    }
  }
  coef <- unique(as.character(coef %||% character()))
  coef <- coef[nzchar(coef)]
  list(
    response = response, family = family, dpar = dpar, nlpar = nlpar,
    group = group, id = id, coef = coef, prior = prior
  )
}

format_re_s2z_context <- function(context) {
  stopifnot(is.list(context))
  out <- paste0(
    "response '", context$response, "', family '", context$family,
    "', dpar '", context$dpar, "', nlpar '", context$nlpar,
    "', group '", context$group, "', ID '", context$id, "'"
  )
  if (length(context$coef)) {
    out <- paste0(
      out, ", coefficient(s) '", paste(context$coef, collapse = ", "), "'"
    )
  }
  if (!is.null(context$prior) && nzchar(as.character(context$prior))) {
    out <- paste0(out, ", prior '", as.character(context$prior), "'")
  }
  out
}

# Suggest stable, predictor-local IDs when one public ID was reused across
# linear predictors. These labels are diagnostic examples only; they do not
# affect the generated model.
re_s2z_local_id_examples <- function(contexts) {
  stopifnot(is.list(contexts), length(contexts) > 1L)
  labels <- vapply(contexts, function(context) {
    predictor <- if (context$nlpar != "<none>") {
      context$nlpar
    } else {
      context$dpar
    }
    label <- paste(context$response, predictor, "s2z", sep = "_")
    label <- gsub("[^[:alnum:]_]+", "_", label)
    gsub("^_+|_+$", "", label)
  }, character(1))
  labels <- make.unique(labels, sep = "_")
  paste0('`id = "', labels, '"`', collapse = ", ")
}

stop_re_s2z <- function(context, capability, problem, remedy) {
  stop2(
    problem, "\nS2Z capability '", capability, "' is unavailable for ",
    format_re_s2z_context(context), ". Remedy: ", remedy
  )
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

# Effective scalar population-level prior for one active coefficient.
re_s2z_prior <- function(prior, bframe, class, coef = "", r = NULL) {
  stopifnot(is.brmsprior(prior), is.bframel(bframe))
  context <- re_s2z_context(bframe, r = r, coef = coef)
  px <- check_prefix(bframe)
  p <- subset2(prior, class = class, ls = px)
  if (class == "b" && nzchar(coef)) {
    pcoef <- subset2(p, coef = coef)
    populated <- nrow(pcoef) == 1L && any(
      nzchar(pcoef$prior) | nzchar(pcoef$lb) | nzchar(pcoef$ub) |
        nzchar(pcoef$tag)
    )
    if (populated) {
      ans <- pcoef
    } else {
      ans <- subset2(p, coef = "")
    }
  } else {
    ans <- subset2(p, coef = "")
  }
  if (!nrow(ans)) {
    return(re_s2z_flat_prior())
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
    stop_re_s2z(
      context, "active_prior_bounds",
      paste0(
        "The sum-to-zero parameterization does not yet support bounded ",
        "population-level priors (coefficient '", coef, "')."
      ),
      "remove the bound or make this coefficient fixed-only."
    )
  }
  if (nrow(ans) && nzchar(ans$tag)) {
    stop_re_s2z(
      context, "active_prior_tag",
      paste0(
        "The sum-to-zero parameterization does not yet support tagged ",
        "population-level priors (coefficient '", coef, "')."
      ),
      "remove the tag or make this coefficient fixed-only."
    )
  }
  value <- if (nrow(ans)) ans$prior else ""
  context$prior <- value
  parse_re_s2z_prior(value, coef = coef, context = context)
}

# Parse the exact scalar priors supported by the S2Z paths. Logistic priors
# use an explicit omitted-mean fallback; the remaining proper distributions
# are conditionally Gaussian.
parse_re_s2z_prior <- function(prior, coef = "", context = NULL) {
  prior <- gsub("[[:space:]]+", "", prior)
  if (!nzchar(prior)) {
    return(re_s2z_flat_prior())
  }
  fail <- function(capability, problem, remedy) {
    if (is.null(context)) {
      stop2(problem, " Remedy: ", remedy)
    }
    context$prior <- prior
    stop_re_s2z(context, capability, problem, remedy)
  }
  call <- try(str2lang(prior), silent = TRUE)
  if (inherits(call, "try-error") || !is.call(call)) {
    fail(
      "active_prior_distribution",
      paste0(
        "Prior '", prior, "' is not supported by the sum-to-zero ",
        "parameterization (coefficient '", coef, "')."
      ),
      "use flat, normal, student_t, cauchy, or logistic with numeric constants."
    )
  }
  dist <- as.character(call[[1]])
  args <- as.list(call[-1])
  number <- function(x) {
    out <- suppressWarnings(as.numeric(deparse0(x)))
    if (length(out) != 1L || !is.finite(out)) {
      fail(
        "active_prior_arguments",
        paste0(
          "All arguments of population-level priors used with the ",
          "sum-to-zero parameterization must currently be numeric ",
          "constants (coefficient '", coef, "')."
        ),
        "replace symbolic arguments with finite numeric constants."
      )
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
  } else if (dist == "logistic" && length(args) == 2L) {
    out <- list(
      dist = "logistic", location = number(args[[1]]),
      scale = number(args[[2]]), df = NA_real_
    )
  } else {
    fail(
      "active_prior_distribution",
      paste0(
        "Prior '", prior, "' is not supported by the sum-to-zero ",
        "parameterization. Supported population-level priors are flat, ",
        "normal, student_t, cauchy, and logistic (coefficient '", coef, "')."
      ),
      "choose one of the supported priors or make this coefficient fixed-only."
    )
  }
  if (!(out$scale > 0) || out$dist == "student" && !(out$df > 0)) {
    fail(
      "active_prior_arguments",
      paste0(
        "Scale and degrees-of-freedom arguments must be positive in ",
        "population-level priors used with the sum-to-zero ",
        "parameterization (coefficient '", coef, "')."
      ),
      "supply a positive finite scale and, for student_t, positive degrees of freedom."
    )
  }
  out
}

# Check fixed group scales before the Gaussian precision is constructed.
validate_re_s2z_sd_prior <- function(prior, r, bframe = NULL) {
  stopifnot(is.brmsprior(prior), is.reframe(r), has_rows(r))
  px <- check_prefix(r)
  p <- subset2(
    prior, class = "sd", coef = c(r$coef, ""),
    group = c(r$group[1], ""), ls = px
  )
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
      problem <- paste0(
        "Group-level standard deviations fixed with 'constant' must be ",
        "positive numeric scalars for gr(..., s2z = TRUE) ",
        "(coefficient '", coef, "')."
      )
      if (is.null(bframe)) {
        stop2(problem)
      }
      stop_re_s2z(
        re_s2z_context(bframe, r = r, coef = coef, prior = value),
        "fixed_sd_positive", problem,
        "use constant(value) with a finite value greater than zero, or estimate the scale."
      )
    }
  }
  invisible(NULL)
}

# Is this the location predictor of an ordinal response? Ordinal auxiliary
# predictors keep the ordinary scalar-coordinate map used by PR 1.
re_s2z_is_ordinal_location <- function(bframe) {
  stopifnot(is.bframel(bframe))
  if (!is_ordinal(bframe$family)) {
    return(FALSE)
  }
  px <- check_prefix(bframe)
  if (nzchar(px$nlpar %||% "")) {
    return(FALSE)
  }
  dpar <- dpar_class(px$dpar %||% "")
  !length(dpar) || !nzchar(dpar) || identical(dpar, "mu")
}

# Primitive ordinal threshold coordinates. Flexible thresholds retain one
# primitive per threshold. Equidistant thresholds retain only their translated
# first threshold; delta is a spacing and is not part of the omitted-location
# map. Group-specific vectors are kept separate even when their lengths agree.
re_s2z_threshold_frame <- function(bframe) {
  stopifnot(is.bframel(bframe), re_s2z_is_ordinal_location(bframe))
  groups <- as.character(get_thres_groups(bframe))
  if (!length(groups)) {
    groups <- ""
  }
  equidistant <- has_equidistant_thres(bframe)
  out <- lapply(seq_along(groups), function(i) {
    group <- groups[i]
    if (equidistant) {
      coef <- ""
      threshold_index <- 1L
      kind <- "first"
    } else {
      coef <- as.character(get_thres(bframe, group = group))
      threshold_index <- seq_along(coef)
      kind <- rep("flexible", length(coef))
    }
    data.frame(
      q = NA_integer_, group = rep(group, length(threshold_index)),
      group_index = rep.int(i, length(threshold_index)), coef = coef,
      threshold_index = threshold_index, kind = kind,
      stringsAsFactors = FALSE
    )
  })
  out <- do_call(rbind, out)
  out$q <- seq_rows(out)
  out$source <- ifelse(out$kind == "first", "first_Intercept", "Intercept")
  out$qname <- ifelse(
    out$kind == "first",
    ifelse(
      nzchar(out$group), paste0("first_Intercept[", out$group, "]"),
      "first_Intercept"
    ),
    ifelse(
      nzchar(out$group),
      paste0("Intercept[", out$group, ",", out$coef, "]"),
      paste0("Intercept[", out$coef, "]")
    )
  )
  out
}

# Center X exactly as the generated ordinal likelihood does. This R copy is
# used only to validate and record the affine map; Stan recomputes Xc from X.
re_s2z_likelihood_design <- function(bframe) {
  stopifnot(is.bframel(bframe))
  X <- bframe$sdata$fe$X
  if (is.null(X)) {
    n <- length(bframe$frame$resp$values %||% numeric())
    X <- matrix(numeric(), nrow = n, ncol = 0L)
  }
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (isTRUE(bframe$frame$fe$center) && ncol(X)) {
    X <- sweep(X, 2L, colMeans(X), FUN = "-")
  }
  X
}

# Build a name-preserving candidate slope map. The exact affine solver below
# keeps this selector when it already works, but may replace columns when a
# different fixed-effect basis spans the same varying design (for example,
# treatment-coded factors or algebraically rescaled terms). A varying
# intercept deliberately has no matching population coordinate in an ordinal
# model; its column is represented by a alone and shifts every threshold row.
re_s2z_ordinal_C <- function(r, slope_names) {
  stopifnot(is.reframe(r), is.character(slope_names))
  C <- matrix(
    0, nrow = length(slope_names), ncol = nrow(r),
    dimnames = list(slope_names, r$coef)
  )
  slope <- r$coef != "Intercept"
  matched <- match(r$coef[slope], slope_names)
  present <- !is.na(matched)
  if (any(present)) {
    C[cbind(matched[present], which(slope)[present])] <- 1
  }
  C
}

# Construct the raw group design without consulting a descriptor. This avoids
# a descriptor/design recursion while the exact affine map is being built.
.re_s2z_design_matrix_r <- function(r, data) {
  stopifnot(is.reframe(r), has_rows(r), is.data.frame(data))
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

# Build ordinal block descriptors around the exact affine identity
# Z = 1 a' + Xc C. q always consists of threshold primitives first and
# population slopes second. H follows q = theta - H m, hence its threshold
# rows are -a and its population-slope rows are +C.
.re_s2z_ordinal_descriptors <- function(bframe, data = NULL) {
  stopifnot(is.bframel(bframe), re_s2z_is_ordinal_location(bframe))
  re <- bframe$frame$re
  ids <- unique(re$id[re$s2z])
  threshold <- re_s2z_threshold_frame(bframe)
  threshold_q <- threshold$q
  slope_names <- as.character(
    bframe$frame$fe$vars_stan %||% character()
  )
  slope_q <- length(threshold_q) + seq_along(slope_names)
  qnames <- c(threshold$qname, slope_names)
  qtype <- c(rep("threshold", length(threshold_q)),
             rep("b", length(slope_q)))
  Xc <- re_s2z_likelihood_design(bframe)
  stopifnot(ncol(Xc) == length(slope_names))
  center <- bframe$frame$fe$center

  out <- lapply(ids, function(id) {
    r <- subset2(re, id = id)
    C <- re_s2z_ordinal_C(r, slope_names)
    match_slope <- match(r$coef, slope_names)
    match_slope[r$coef == "Intercept"] <- NA_integer_
    match_q <- ifelse(
      is.na(match_slope), NA_integer_, length(threshold_q) + match_slope
    )

    if (!is.null(data)) {
      Z <- .re_s2z_design_matrix_r(r, data = data)
      residual <- Z - Xc %*% C
      # Preserve exact name-matched selectors when possible. For every
      # remaining column, solve the general affine system against [1, Xc].
      # lm.fit() also supplies a deterministic estimable solution when the
      # fixed design is rank deficient; aliased coefficients are set to zero.
      candidate_a <- unname(residual[1L, ])
      candidate_deviation <- sweep(residual, 2L, candidate_a, FUN = "-")
      candidate_scale <- pmax(
        1,
        apply(abs(Z), 2L, max),
        apply(abs(Xc %*% C), 2L, max)
      )
      candidate_tolerance <- 64 * .Machine$double.eps * candidate_scale
      candidate_error <- apply(abs(candidate_deviation), 2L, max)
      solve_columns <- which(
        !is.finite(candidate_error) |
          candidate_error > candidate_tolerance
      )
      if (length(solve_columns)) {
        affine_design <- cbind(1, Xc)
        fit <- lm.fit(
          affine_design, Z[, solve_columns, drop = FALSE],
          tol = sqrt(.Machine$double.eps), singular.ok = TRUE
        )
        affine_coef <- fit$coefficients
        if (!is.matrix(affine_coef)) {
          affine_coef <- matrix(affine_coef, ncol = 1L)
        }
        affine_coef[is.na(affine_coef)] <- 0
        # Remove numerical dust so ordinary selectors retain literal 0/1
        # entries in generated Stan while nontrivial maps remain unchanged.
        for (value in c(-1, 0, 1)) {
          close <- abs(affine_coef - value) <=
            128 * .Machine$double.eps * pmax(1, abs(affine_coef))
          affine_coef[close] <- value
        }
        candidate_a[solve_columns] <- affine_coef[1L, ]
        if (nrow(C)) {
          C[, solve_columns] <- affine_coef[-1L, , drop = FALSE]
        }
        residual <- Z - Xc %*% C
      }
      # The first row is the affine constant used in the checked identity.
      # A scale-aware machine-precision envelope admits only arithmetic
      # roundoff, not a substantively approximate design match.
      a <- unname(residual[1L, ])
      deviation <- sweep(residual, 2L, a, FUN = "-")
      affine_error_by_coef <- if (ncol(deviation)) {
        apply(abs(deviation), 2L, max)
      } else {
        numeric()
      }
      names(affine_error_by_coef) <- r$coef
      affine_error <- if (length(affine_error_by_coef)) {
        max(affine_error_by_coef)
      } else {
        0
      }
      affine_scale <- max(1, abs(Z), abs(Xc %*% C), na.rm = TRUE)
      affine_tolerance <- 64 * .Machine$double.eps * affine_scale
      affine_ok_by_coef <- is.finite(affine_error_by_coef) &
        affine_error_by_coef <= affine_tolerance
      affine_ok <- all(affine_ok_by_coef)
    } else {
      # Frames created without source data should normally carry the cached
      # descriptor produced by validate_re_s2z_design(). This fallback is the
      # ordinary named-design map and remains useful for legacy frames.
      Xraw <- bframe$sdata$fe$X
      means_X <- if (isTRUE(center) && !is.null(Xraw) && ncol(Xraw)) {
        colMeans(Xraw)
      } else {
        rep(0, nrow(C))
      }
      a <- as.numeric(crossprod(means_X, C))
      a[r$coef == "Intercept"] <- 1
      affine_error <- NA_real_
      affine_error_by_coef <- setNames(rep(NA_real_, nrow(r)), r$coef)
      affine_tolerance <- NA_real_
      affine_ok <- NA
      affine_ok_by_coef <- setNames(rep(NA, nrow(r)), r$coef)
    }

    H <- matrix(
      0, nrow = length(qnames), ncol = nrow(r),
      dimnames = list(qnames, r$coef)
    )
    if (length(threshold_q)) {
      H[threshold_q, ] <- matrix(
        rep(-a, times = length(threshold_q)),
        nrow = length(threshold_q), byrow = TRUE
      )
    }
    if (length(slope_q)) {
      H[slope_q, ] <- C
    }
    block_active_q <- which(rowSums(abs(H)) > 0)
    nlist(
      id, r, ordinal = TRUE, qnames, qtype, threshold_q, slope_q,
      threshold, match_q, match_slope, center,
      fixef = slope_names,
      p = usc(combine_prefix(check_prefix(bframe))),
      context = re_s2z_context(bframe, r = r),
      a, C, H, block_active_q, affine_ok, affine_error,
      affine_error_by_coef, affine_ok_by_coef, affine_tolerance
    )
  })

  active_q <- sort(unique(unlist(lapply(out, `[[`, "block_active_q"),
                                    use.names = FALSE)))
  inactive_q <- setdiff(seq_along(qnames), active_q)
  active_names <- qnames[active_q]
  inactive_names <- qnames[inactive_q]
  active_index <- match(seq_along(qnames), active_q)
  lapply(out, function(info) {
    info$active_q <- active_q
    info$inactive_q <- inactive_q
    info$active_names <- active_names
    info$inactive_names <- inactive_names
    info$active_index <- active_index
    info$match_active <- match(info$match_q, active_q)
    info
  })
}

# Build validated-block descriptors without looking at priors. Active rows are
# the union of all rows touched by the local omitted-mean maps. In a centered
# fixed design, any matched non-intercept slope also structurally touches the
# temporary intercept row, even when its realized column mean happens to be 0.
.re_s2z_descriptors <- function(bframe, data = NULL) {
  stopifnot(is.bframel(bframe))
  if (!has_re_s2z(bframe)) {
    return(list())
  }
  if (re_s2z_is_ordinal_location(bframe)) {
    return(.re_s2z_ordinal_descriptors(bframe, data = data))
  }
  re <- bframe$frame$re
  ids <- unique(re$id[re$s2z])
  qnames <- bframe$frame$fe$vars
  center <- bframe$frame$fe$center
  out <- lapply(ids, function(id) {
    r <- subset2(re, id = id)
    match_q <- match(r$coef, qnames)
    nlist(
      id, r, ordinal = FALSE, qnames,
      qtype = rep("b", length(qnames)), match_q, center,
      fixef = bframe$frame$fe$vars_stan,
      p = usc(combine_prefix(check_prefix(bframe))),
      context = re_s2z_context(bframe, r = r)
    )
  })
  active_q <- sort(unique(unlist(lapply(out, function(info) {
    rows <- info$match_q[!is.na(info$match_q)]
    if (info$center && any(info$r$coef != "Intercept")) {
      rows <- c(rows, 1L)
    }
    rows
  }), use.names = FALSE)))
  inactive_q <- setdiff(seq_along(qnames), active_q)
  active_names <- qnames[active_q]
  inactive_names <- qnames[inactive_q]
  active_index <- match(seq_along(qnames), active_q)
  lapply(out, function(info) {
    info$block_active_q <- sort(unique(c(
      info$match_q[!is.na(info$match_q)],
      if (info$center && any(info$r$coef != "Intercept")) 1L
    )))
    info$active_q <- active_q
    info$inactive_q <- inactive_q
    info$active_names <- active_names
    info$inactive_names <- inactive_names
    info$active_index <- active_index
    info$match_active <- match(info$match_q, active_q)
    info
  })
}

# Extract and parse the ordinary prior for one primitive ordinal threshold.
# Threshold metadata stays on the parsed specification so code generation can
# retain distinct priors for grouped vectors and unequal threshold counts.
re_s2z_threshold_prior <- function(prior, bframe, threshold, r) {
  stopifnot(
    is.brmsprior(prior), is.bframel(bframe), is.data.frame(threshold),
    nrow(threshold) == 1L, is.reframe(r), has_rows(r)
  )
  px <- check_prefix(bframe)
  group <- threshold$group[[1L]]
  coef <- threshold$coef[[1L]]
  p <- subset2(
    prior, class = "Intercept", coef = c(coef, ""),
    group = c(group, ""), ls = px
  )
  populated <- function(x) {
    if (!nrow(x)) {
      return(logical())
    }
    nzchar(x$prior) | nzchar(x$lb) | nzchar(x$ub) | nzchar(x$tag)
  }
  choose_most_specific <- function(x) {
    if (nrow(x) <= 1L) {
      return(x)
    }
    exact_group <- x$group == group
    if (any(exact_group)) {
      x <- x[exact_group, , drop = FALSE]
    }
    if (nrow(x) > 1L) {
      prefix_vars <- intersect(vars_prefix(), names(x))
      specificity <- rowSums(nzchar(x[, prefix_vars, drop = FALSE]))
      x <- x[which.max(specificity), , drop = FALSE]
    }
    x
  }

  ans <- p[FALSE, , drop = FALSE]
  if (nzchar(coef)) {
    pcoef <- subset2(p, coef = coef)
    pcoef <- pcoef[populated(pcoef), , drop = FALSE]
    ans <- choose_most_specific(pcoef)
  }
  if (!nrow(ans)) {
    pbase <- subset2(p, coef = "")
    pbase <- pbase[populated(pbase), , drop = FALSE]
    ans <- choose_most_specific(pbase)
  }

  label <- if (threshold$kind[[1L]] == "first") {
    "first threshold"
  } else {
    paste0("threshold ", coef)
  }
  if (nzchar(group)) {
    label <- paste0(label, " in threshold group '", group, "'")
  }
  context <- re_s2z_context(bframe, r = r, coef = label)
  if (nrow(ans) && (nzchar(ans$lb) || nzchar(ans$ub))) {
    stop_re_s2z(
      context, "ordinal_threshold_prior_bounds",
      paste0("Bounded ordinal threshold priors are not supported by the ",
             "sum-to-zero parameterization (", label, ")."),
      "remove the bound or use s2z = FALSE for this ordinal predictor."
    )
  }
  if (nrow(ans) && nzchar(ans$tag)) {
    stop_re_s2z(
      context, "ordinal_threshold_prior_tag",
      paste0("Tagged ordinal threshold priors are not supported by the ",
             "sum-to-zero parameterization (", label, ")."),
      "remove the tag or use s2z = FALSE for this ordinal predictor."
    )
  }
  value <- if (nrow(ans)) ans$prior[[1L]] else ""
  context$prior <- value
  spec <- parse_re_s2z_prior(value, coef = label, context = context)
  if (identical(spec$dist, "logistic")) {
    stop_re_s2z(
      context, "ordinal_active_prior_distribution",
      paste0(
        "Logistic priors on S2Z-active ordinal coordinates are not yet ",
        "supported (", label, ")."
      ),
      "use a normal, student_t, or cauchy threshold prior."
    )
  }
  spec$class <- "Intercept"
  spec$group <- group
  spec$coef <- coef
  spec$q <- threshold$q[[1L]]
  spec$kind <- threshold$kind[[1L]]
  spec$threshold_index <- threshold$threshold_index[[1L]]
  spec$prior <- value
  spec
}

# Attach active-coordinate priors while leaving fixed-only slope slots as
# neutral placeholders. Every ordinal threshold primitive receives its own
# specification, including H=0 rows, so the joint dense system has one clear
# owner for all threshold densities.
.re_s2z_attach_priors <- function(infos, bframe, prior) {
  if (!length(infos)) {
    return(infos)
  }
  qnames <- infos[[1L]]$qnames
  active_q <- infos[[1L]]$active_q
  specs <- rep(list(re_s2z_flat_prior()), length(qnames))
  if (isTRUE(infos[[1L]]$ordinal)) {
    threshold <- infos[[1L]]$threshold
    for (i in seq_rows(threshold)) {
      qi <- threshold$q[i]
      specs[[qi]] <- re_s2z_threshold_prior(
        prior, bframe, threshold = threshold[i, , drop = FALSE],
        r = infos[[1L]]$r
      )
    }
    prior_q <- intersect(active_q, infos[[1L]]$slope_q)
  } else {
    prior_q <- active_q
  }
  for (i in prior_q) {
    # Attribute prior diagnostics to a block that actually touches this
    # population coordinate. In a multi-block predictor the first block need
    # not contain the failing coefficient.
    touching <- which(vapply(infos, function(info) {
      i %in% info$block_active_q
    }, logical(1)))
    r_context <- infos[[touching[1L]]]$r
    if (!isTRUE(infos[[1L]]$ordinal) &&
        infos[[1L]]$center && qnames[i] == "Intercept") {
      specs[[i]] <- re_s2z_prior(
        prior, bframe, class = "Intercept", r = r_context
      )
    } else {
      specs[[i]] <- re_s2z_prior(
        prior, bframe, class = "b", coef = qnames[i], r = r_context
      )
    }
    if (isTRUE(infos[[1L]]$ordinal) &&
        identical(specs[[i]]$dist, "logistic")) {
      context <- re_s2z_context(
        bframe, r = r_context, coef = qnames[i],
        prior = specs[[i]]$prior %||% "logistic"
      )
      stop_re_s2z(
        context, "ordinal_active_prior_distribution",
        paste0(
          "Logistic priors on S2Z-active ordinal coordinates are not yet ",
          "supported (coefficient '", qnames[i], "')."
        ),
        "use a normal, student_t, cauchy, or flat slope prior."
      )
    }
  }
  lapply(infos, function(info) {
    info$prior <- specs
    if (isTRUE(info$ordinal)) {
      info$threshold_prior <- unname(specs[info$threshold_q])
    }
    info
  })
}

# Describe all local S2Z blocks in one linear predictor. Blocks retain their
# own covariance models but share one active-coordinate omitted-mean system.
re_s2z_infos <- function(bframe, prior = NULL, data = NULL) {
  stopifnot(is.bframel(bframe))
  if (!has_re_s2z(bframe)) {
    return(list())
  }
  if (is.null(data) &&
      is.list(bframe$frame$re_s2z) && length(bframe$frame$re_s2z)) {
    out <- bframe$frame$re_s2z
  } else {
    out <- .re_s2z_descriptors(bframe, data = data)
  }
  if (!is.null(prior)) {
    out <- .re_s2z_attach_priors(out, bframe = bframe, prior = prior)
  }
  out
}

# Describe one local S2Z block, including the fixed/group design mapping.
re_s2z_info <- function(bframe, prior = NULL, id = NULL) {
  stopifnot(is.bframel(bframe))
  infos <- re_s2z_infos(bframe, prior = prior)
  if (!length(infos)) {
    return(NULL)
  }
  ids <- vapply(infos, `[[`, numeric(1), "id")
  if (is.null(id) && length(ids) != 1L) {
    stop2("An explicit group-level ID is required when describing multiple ",
          "sum-to-zero blocks.")
  }
  if (is.null(id)) {
    return(infos[[1L]])
  }
  if (length(id) != 1L || !id %in% ids) {
    stop2("Invalid sum-to-zero group-level ID.")
  }
  infos[[match(id, ids)]]
}

# Construct an S2Z design matrix in covariance-block coefficient order. A
# block may be assembled from multiple grouping terms that share an ID.
re_s2z_design_matrix <- function(bframe, data, id = NULL) {
  stopifnot(is.bframel(bframe))
  re <- bframe$frame$re
  ids <- unique(re$id[re$s2z])
  if (!length(ids)) {
    return(matrix(numeric(), nrow = nrow(data), ncol = 0L))
  }
  if (is.null(id) && length(ids) != 1L) {
    stop2("An explicit group-level ID is required when constructing a ",
          "design for multiple sum-to-zero blocks.")
  }
  if (is.null(id)) {
    id <- ids[[1L]]
  }
  if (length(id) != 1L || !id %in% ids) {
    stop2("Invalid sum-to-zero group-level ID.")
  }
  .re_s2z_design_matrix_r(subset2(re, id = id), data = data)
}

# Validate local structural eligibility before any fixed/group design match is
# attempted. The ordering here is intentional: family and term capabilities
# should diagnose themselves even when the corresponding ordinary fixed column
# is absent (for example, ordinal threshold models).
validate_re_s2z_structure <- function(bframe, data) {
  stopifnot(is.bframel(bframe), is.data.frame(data))
  if (!has_re_s2z(bframe)) {
    return(invisible(NULL))
  }
  re <- bframe$frame$re
  ids <- unique(re$id[re$s2z])
  first_r <- subset2(re, id = ids[1L])
  ordinal_location <- re_s2z_is_ordinal_location(bframe)
  if (ordinal_location) {
    context <- re_s2z_context(bframe, r = first_r)
    if (fix_intercepts(bframe)) {
      stop_re_s2z(
        context, "ordinal_fixed_thresholds",
        paste0(
          "Fixed or shared ordinal-mixture thresholds cannot absorb the ",
          "omitted sum-to-zero group-effect mean."
        ),
        "use component-specific flexible thresholds or set s2z = FALSE for this ordinal location predictor."
      )
    }
    if (has_sum_to_zero_thres(bframe)) {
      stop_re_s2z(
        context, "ordinal_sum_to_zero_thresholds",
        paste0(
          "Sum-to-zero ordinal thresholds do not retain the common ",
          "location coordinate required by gr(..., s2z = TRUE)."
        ),
        "use flexible or equidistant thresholds, or set s2z = FALSE."
      )
    }
    if (is.customfamily(bframe$family)) {
      stop_re_s2z(
        context, "ordinal_custom_family",
        paste0(
          "Custom ordinal likelihoods have not been validated for the ",
          "sum-to-zero threshold-location map."
        ),
        "use a built-in ordinal family or set s2z = FALSE."
      )
    }
  }
  if (order_intercepts(bframe)) {
    stop_re_s2z(
      re_s2z_context(bframe, r = first_r), "ordered_mixture",
      paste0(
        "Ordered mixture intercepts are not yet supported together with ",
        "gr(..., s2z = TRUE)."
      ),
      "disable ordered mixture intercepts or use a conventional group-level term."
    )
  }
  for (id in ids) {
    r <- subset2(re, id = id)
    context <- re_s2z_context(bframe, r = r)
    if (ordinal_location && any(r$type == "cs")) {
      stop_re_s2z(
        context, "ordinal_category_specific",
        paste0(
          "Category-specific sum-to-zero group effects require a ",
          "category-specific omitted-mean map."
        ),
        "use ordinary category-specific group effects with s2z = FALSE, or use a non-category-specific S2Z block."
      )
    }
    if (!all(r$s2z)) {
      stop_re_s2z(
        context, "same_id_s2z",
        paste0(
          "All coefficients sharing a group-level ID must use the same ",
          "'s2z' setting."
        ),
        "give every term with this ID the same s2z setting or use distinct IDs."
      )
    }
    if (length(unique(r$group)) != 1L) {
      stop_re_s2z(
        context, "same_id_group",
        paste0(
          "All coefficients sharing a sum-to-zero group-level ID must use ",
          "the same grouping factor."
        ),
        "use one grouping factor per S2Z ID or assign distinct IDs."
      )
    }
    if (length(unique(r$dist)) != 1L) {
      stop_re_s2z(
        context, "same_id_distribution",
        paste0(
          "All coefficients sharing a group-level ID must use the same ",
          "group-level distribution."
        ),
        "use one group-effect distribution per S2Z ID or assign distinct IDs."
      )
    }
    if (length(unique(r$cor)) != 1L) {
      stop_re_s2z(
        context, "same_id_correlation",
        paste0(
          "All coefficients sharing a group-level ID must use the same ",
          "'cor' setting."
        ),
        "use one correlation setting per S2Z ID or assign distinct IDs."
      )
    }
    if (any(nzchar(r$gtype), na.rm = TRUE) ||
        any(nzchar(r$type), na.rm = TRUE)) {
      stop_re_s2z(
        context, "ordinary_gr_only",
        paste0(
          "The sum-to-zero parameterization currently supports only ",
          "ordinary 'gr' terms."
        ),
        "remove the special group coefficient or use s2z = FALSE for this term."
      )
    }
    has_pw <- any(vapply(r$gcall, function(gcall) {
      isTRUE(nzchar(gcall$pw %||% ""))
    }, logical(1)))
    has_by <- any(nzchar(r$by), na.rm = TRUE)
    has_cov <- any(nzchar(r$cov), na.rm = TRUE)
    if (has_by || has_cov || has_pw) {
      capability <- if (has_by) "by" else if (has_cov) "cov" else "pw"
      remedy <- switch(
        capability,
        by = paste0(
          "remove 'by' or use s2z = FALSE to retain separate ",
          "variance-covariance matrices for its levels."
        ),
        cov = paste0(
          "remove 'cov' to use independent grouping levels, or use ",
          "s2z = FALSE to retain the supplied group covariance matrix."
        ),
        pw = paste0(
          "remove 'pw' or use s2z = FALSE to retain the group-level ",
          "prior weights."
        )
      )
      stop_re_s2z(
        context, capability,
        paste0(
          "Argument '", capability, "' is not yet supported together ",
          "with gr(..., s2z = TRUE)."
        ),
        remedy
      )
    }
    if (isTRUE(bframe$frame$fe$sparse)) {
      stop_re_s2z(
        context, "sparse",
        paste0(
          "Sparse and QR-decomposed population-level design matrices are ",
          "not yet supported together with gr(..., s2z = TRUE)."
        ),
        "disable sparse fixed-effect design generation for this predictor."
      )
    }
    if (!identical(bframe$frame$fe$decomp, "none")) {
      stop_re_s2z(
        context, "qr",
        paste0(
          "Sparse and QR-decomposed population-level design matrices are ",
          "not yet supported together with gr(..., s2z = TRUE)."
        ),
        "set decomp = 'none' for this predictor."
      )
    }
    # During model construction the physical constraint is defined over the
    # levels actually represented in the likelihood. Retained but unused
    # factor levels would otherwise introduce unobserved balancing effects.
    # Post-processing frames carry a fitted basis and may legitimately contain
    # only a subset of the fitted levels in newdata.
    levels <- get_levels(r)[[r$group[1L]]]
    if (is.null(bframe$basis)) {
      group_name <- r$gcall[[1L]]$groups
      group_data <- get(group_name, data)
      observed_levels <- unique(as.character(group_data[!is.na(group_data)]))
      if (length(observed_levels) < 2L) {
        stop_re_s2z(
          context, "minimum_levels",
          paste0(
            "The sum-to-zero parameterization requires at least two observed ",
            "grouping levels."
          ),
          "supply data with at least two observed levels or use a fixed effect."
        )
      }
      unused_levels <- setdiff(as.character(levels), observed_levels)
      if (length(unused_levels)) {
        stop_re_s2z(
          context, "unused_levels",
          paste0(
            "The sum-to-zero parameterization cannot retain unobserved ",
            "grouping levels: ", paste(unused_levels, collapse = ", "), "."
          ),
          "drop unused factor levels or set drop_unused_levels = TRUE."
        )
      }
    } else if (length(levels) < 2L) {
      stop_re_s2z(
        context, "minimum_levels",
        paste0(
          "The sum-to-zero parameterization requires at least two observed ",
          "grouping levels."
        ),
        "supply data with at least two observed levels or use a fixed effect."
      )
    }
  }
  invisible(NULL)
}

# Validate that matched fixed and group columns are numerically identical.
validate_re_s2z_design <- function(bframe, data) {
  stopifnot(is.bframel(bframe))
  if (!has_re_s2z(bframe)) {
    return(invisible(list()))
  }
  X <- bframe$sdata$fe$X
  infos <- re_s2z_infos(bframe, data = data)
  missing <- unique(unlist(lapply(infos, function(info) {
    if (isTRUE(info$ordinal)) {
      # A differently named varying column is valid when it lies in the
      # affine span of the ordinal likelihood design. Diagnose it as missing
      # only when both the name-based hint and the actual affine solve fail.
      info$r$coef[
        info$r$coef != "Intercept" & is.na(info$match_slope) &
          !info$affine_ok_by_coef
      ]
    } else {
      info$r$coef[is.na(info$match_q)]
    }
  }), use.names = FALSE))
  if (length(missing)) {
    first <- which(vapply(infos, function(info) {
      if (isTRUE(info$ordinal)) {
        any(
          info$r$coef != "Intercept" & is.na(info$match_slope) &
            !info$affine_ok_by_coef
        )
      } else {
        any(is.na(info$match_q))
      }
    }, logical(1)))[1L]
    level <- if (isTRUE(infos[[first]]$ordinal)) {
      "non-intercept sum-to-zero group-level slope"
    } else {
      "sum-to-zero group-level coefficient"
    }
    stop_re_s2z(
      re_s2z_context(bframe, r = infos[[first]]$r, coef = missing),
      "matching_name",
      paste0(
        "Every ", level, " must have a matching population-level design ",
        "column. Missing: ",
        paste(missing, collapse = ", "), "."
      ),
      "add population-level slope terms with these coefficient names or remove the unmatched group effects."
    )
  }
  if (isTRUE(infos[[1L]]$ordinal)) {
    bad <- which(!vapply(infos, `[[`, logical(1), "affine_ok"))
    if (length(bad)) {
      first <- infos[[bad[1L]]]
      bad_coef <- unique(unlist(lapply(infos[bad], function(info) {
        names(info$affine_ok_by_coef)[!info$affine_ok_by_coef]
      }), use.names = FALSE))
      stop_re_s2z(
        re_s2z_context(bframe, r = first$r, coef = bad_coef),
        "affine_identity",
        paste0(
          "The ordinal group design does not satisfy the required exact ",
          "affine identity Z = 1 a^T + Xc C for coefficient(s) '",
          paste(bad_coef, collapse = ", "), "'."
        ),
        paste0(
          "use the same named slope expressions and contrasts in the ",
          "population- and group-level designs; varying intercepts may map ",
          "directly to threshold location."
        )
      )
    }
    return(invisible(infos))
  }
  mismatches <- list()
  for (info in infos) {
    Z <- re_s2z_design_matrix(bframe, data = data, id = info$id)
    if (ncol(Z) != nrow(info$r)) {
      stop_re_s2z(
        info$context, "internal_design",
        "Internal mismatch in the sum-to-zero group-level design matrix.",
        "report this model and data as a brms issue."
      )
    }
    for (j in seq_len(ncol(Z))) {
      # The omitted-mean identity is exact only when the two likelihood design
      # columns are numerically identical; even a tiny tolerance would define
      # a different posterior.
      same <- identical(
        unname(Z[, j]), unname(X[, info$match_q[j]])
      )
      if (!same) {
        mismatches[[length(mismatches) + 1L]] <- list(
          info = info, coef = info$r$coef[j]
        )
      }
    }
  }
  if (length(mismatches)) {
    coef <- unique(vapply(mismatches, `[[`, character(1), "coef"))
    first <- mismatches[[1L]]$info
    stop_re_s2z(
      re_s2z_context(bframe, r = first$r, coef = coef),
      "matching_values",
      paste0(
        "Population- and group-level design columns must be identical for ",
        "gr(..., s2z = TRUE). Mismatch for coefficient(s) '",
        paste(coef, collapse = ", "), "'."
      ),
      "use the identical population- and group-level term expressions and contrasts."
    )
  }
  invisible(infos)
}

# Validate prior-dependent and cross-predictor constraints only after the
# complete effective prior table exists.
validate_re_s2z_prior_global <- function(bframe, prior) {
  stopifnot(is.anybrmsframe(bframe), is.brmsprior(prior))
  all_frames <- all_bframel(bframe)
  frames <- Filter(has_re_s2z, all_frames)
  if (!length(frames)) {
    return(invisible(NULL))
  }
  s2z_ids <- unique(unlist(lapply(frames, function(x) {
    x$frame$re$id[x$frame$re$s2z]
  }), use.names = FALSE))
  for (id in s2z_ids) {
    occurrences <- Filter(function(x) {
      r <- x$frame$re
      has_rows(r) && id %in% r$id
    }, all_frames)
    if (length(occurrences) > 1L) {
      contexts <- lapply(occurrences, function(x) {
        re_s2z_context(x, r = subset2(x$frame$re, id = id))
      })
      affected <- paste0(
        "Affected linear predictors:\n  - ",
        paste(
          vapply(contexts, format_re_s2z_context, character(1)),
          collapse = "\n  - "
        )
      )
      public_id <- contexts[[1L]]$id
      stop_re_s2z(
        contexts[[1L]], "cross_predictor_id",
        paste0(
          "A sum-to-zero group-level ID cannot span multiple linear ",
          "predictors.\n", affected
        ),
        paste0(
          "use a distinct ID in every listed linear predictor (for example, ",
          re_s2z_local_id_examples(contexts), "). For `mvbind(...)` or ",
          "category shorthand, omit the shared `| ", public_id,
          " |` tag or `id = \"", public_id, "\"` so brms allocates ",
          "predictor-local IDs; alternatively, expand the model into ",
          "separate `bf()` response formulas with those distinct IDs. ",
          "These rewrites do not retain cross-predictor group-effect ",
          "correlations; use s2z = FALSE if those correlations are required."
        )
      )
    }
  }
  if (get_sample_prior(prior) != "no") {
    first <- frames[[1L]]
    r <- first$frame$re[first$frame$re$s2z, , drop = FALSE]
    stop_re_s2z(
      re_s2z_context(first, r = r), "sample_prior",
      paste0(
        "Argument 'sample_prior' is not yet supported together with ",
        "gr(..., s2z = TRUE)."
      ),
      "set sample_prior = 'no' for this model."
    )
  }
  for (x in frames) {
    infos <- re_s2z_infos(x)
    active_q <- infos[[1L]]$active_q
    if (isTRUE(infos[[1L]]$ordinal)) {
      active_b_q <- intersect(active_q, infos[[1L]]$slope_q)
      active_b <- infos[[1L]]$qnames[active_b_q]
      active_intercept <- logical()
    } else {
      active_names <- infos[[1L]]$qnames[active_q]
      active_intercept <- infos[[1L]]$center & active_names == "Intercept"
      active_b_q <- active_q[!active_intercept]
      active_b <- active_names[!active_intercept]
    }
    if (length(active_b) && has_special_prior(prior, x, class = "b")) {
      touching <- which(vapply(infos, function(info) {
        active_b_q[1L] %in% info$block_active_q
      }, logical(1)))
      r_context <- infos[[touching[1L]]]$r
      stop_re_s2z(
        re_s2z_context(
          x, r = r_context, coef = active_b, prior = "special b prior"
        ),
        "active_special_prior",
        paste0(
          "Special population-level priors are not yet supported on ",
          "S2Z-active coordinates."
        ),
        "use a supported scalar prior on the active coordinates or make them fixed-only."
      )
    }
    active_threshold <- isTRUE(infos[[1L]]$ordinal) &&
      length(infos[[1L]]$threshold_q)
    if ((any(active_intercept) || active_threshold) &&
        has_special_prior(prior, x, class = "Intercept")) {
      coef <- if (active_threshold) {
        infos[[1L]]$qnames[infos[[1L]]$threshold_q]
      } else {
        "Intercept"
      }
      stop_re_s2z(
        re_s2z_context(
          x, r = infos[[1L]]$r, coef = coef,
          prior = "special Intercept prior"
        ),
        "active_special_prior",
        paste0(
          "Special intercept or threshold priors are not yet supported on ",
          "S2Z-active coordinates."
        ),
        "use supported scalar normal, student_t, or cauchy priors, or remove the active location map."
      )
    }
    for (info in infos) {
      if (has_special_prior(prior, check_prefix(info$r), class = "sd")) {
        stop_re_s2z(
          re_s2z_context(
            x, r = info$r, prior = "special group-scale prior"
          ),
          "active_special_prior",
          paste0(
            "Special priors are not yet supported on S2Z group-level ",
            "standard deviations."
          ),
          "use an ordinary sd prior or use s2z = FALSE for this block."
        )
      }
    }
    # Attaching priors performs bounds, tags, argument, and distribution
    # validation, restricted to the union of S2Z-active coordinates.
    re_s2z_infos(x, prior = prior)
    for (info in infos) {
      validate_re_s2z_sd_prior(prior, info$r, bframe = x)
    }
  }
  invisible(NULL)
}

# Compatibility entry point for downstream code while callers migrate to the
# explicitly named prior/global phase.
validate_re_s2z <- function(bframe, prior) {
  validate_re_s2z_prior_global(bframe, prior = prior)
}
