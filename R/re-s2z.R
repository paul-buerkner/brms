# Helpers for the physical sum-to-zero parameterization of group effects

# Does a local linear predictor contain a sum-to-zero group-level block?
has_re_s2z <- function(x) {
  is.bframel(x) && has_rows(x$frame$re) &&
    "s2z" %in% names(x$frame$re) && any(x$frame$re$s2z)
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

# Check fixed group scales before the Gaussian precision is constructed.
validate_re_s2z_sd_prior <- function(prior, r) {
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
      stop2("Group-level standard deviations fixed with 'constant' must ",
            "be positive numeric scalars for gr(..., s2z = TRUE) ",
            "(coefficient '", coef, "').")
    }
  }
  invisible(NULL)
}

# Describe one local S2Z block, including the fixed/group design mapping.
re_s2z_info <- function(bframe, prior = NULL) {
  stopifnot(is.bframel(bframe))
  if (!has_re_s2z(bframe)) {
    return(NULL)
  }
  r <- bframe$frame$re[bframe$frame$re$s2z, , drop = FALSE]
  ids <- unique(r$id)
  if (length(ids) != 1L) {
    stop2("Only one sum-to-zero group-level block is currently ",
          "supported per linear predictor.")
  }
  r_id <- subset2(bframe$frame$re, id = ids)
  if (!all(r_id$s2z)) {
    stop2("All coefficients sharing a group-level ID must use the same ",
          "'s2z' setting.")
  }
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

# Validate S2Z structure after formulas, data, and priors have been resolved.
validate_re_s2z <- function(bframe, prior) {
  stopifnot(is.anybrmsframe(bframe), is.brmsprior(prior))
  frames <- Filter(has_re_s2z, all_bframel(bframe))
  if (!length(frames)) {
    return(invisible(NULL))
  }
  if (get_sample_prior(prior) != "no") {
    stop2("Argument 'sample_prior' is not yet supported together with ",
          "gr(..., s2z = TRUE).")
  }
  ids <- integer()
  for (x in frames) {
    info <- re_s2z_info(x, prior = prior)
    r <- info$r
    if (info$id %in% ids) {
      stop2("A sum-to-zero group-level ID cannot span multiple ",
            "linear predictors.")
    }
    ids <- c(ids, info$id)
    if (r$gtype[1] != "" || any(nzchar(r$type))) {
      stop2("The sum-to-zero parameterization currently ",
            "supports only ordinary 'gr' terms.")
    }
    if (nzchar(r$by[1]) || nzchar(r$cov[1]) ||
        isTRUE(nzchar(r$gcall[[1]]$pw))) {
      stop2("Arguments 'by', 'cov', and 'pw' are not yet supported together ",
            "with gr(..., s2z = TRUE).")
    }
    if (length(unique(r$gn)) != 1L) {
      stop2("A sum-to-zero covariance block must be specified ",
            "in one group-level term.")
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
      stop2("Ordered mixture intercepts are not yet supported together with ",
            "gr(..., s2z = TRUE).")
    }
    if (x$frame$fe$sparse || x$frame$fe$decomp != "none") {
      stop2("Sparse and QR-decomposed population-level design matrices are ",
            "not yet supported together with gr(..., s2z = TRUE).")
    }
    if (has_special_prior(prior, check_prefix(r), class = "sd") ||
        has_special_prior(prior, x, class = "b")) {
      stop2("Special priors are not yet supported together with ",
            "gr(..., s2z = TRUE).")
    }
    validate_re_s2z_sd_prior(prior, r)
  }
  invisible(NULL)
}

# Validate that matched fixed and group columns are numerically identical.
validate_re_s2z_design <- function(bframe, data) {
  stopifnot(is.bframel(bframe))
  if (!has_re_s2z(bframe)) {
    return(invisible(NULL))
  }
  info <- re_s2z_info(bframe)
  Z <- get_model_matrix(info$r$form[[1]], data = data, rename = FALSE)
  X <- bframe$sdata$fe$X
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
  invisible(NULL)
}
