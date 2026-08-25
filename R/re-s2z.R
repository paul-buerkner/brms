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
    if (length(public_id) == 1L && !is.na(public_id) &&
        nzchar(as.character(public_id))) {
      id <- as.character(public_id)
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

# Build validated-block descriptors without looking at priors. Active rows are
# the union of all rows touched by the local omitted-mean maps. In a centered
# fixed design, any matched non-intercept slope also structurally touches the
# temporary intercept row, even when its realized column mean happens to be 0.
.re_s2z_descriptors <- function(bframe) {
  stopifnot(is.bframel(bframe))
  if (!has_re_s2z(bframe)) {
    return(list())
  }
  re <- bframe$frame$re
  ids <- unique(re$id[re$s2z])
  qnames <- bframe$frame$fe$vars
  center <- bframe$frame$fe$center
  out <- lapply(ids, function(id) {
    r <- subset2(re, id = id)
    match_q <- match(r$coef, qnames)
    nlist(
      id, r, qnames, match_q, center,
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

# Attach active-coordinate priors while leaving fixed-only slots as neutral
# placeholders for q-length consumers in the existing Stan generators.
.re_s2z_attach_priors <- function(infos, bframe, prior) {
  if (!length(infos)) {
    return(infos)
  }
  qnames <- infos[[1L]]$qnames
  active_q <- infos[[1L]]$active_q
  specs <- rep(list(re_s2z_flat_prior()), length(qnames))
  for (i in active_q) {
    # Attribute prior diagnostics to a block that actually touches this
    # population coordinate. In a multi-block predictor the first block need
    # not contain the failing coefficient.
    touching <- which(vapply(infos, function(info) {
      i %in% info$block_active_q
    }, logical(1)))
    r_context <- infos[[touching[1L]]]$r
    if (infos[[1L]]$center && qnames[i] == "Intercept") {
      specs[[i]] <- re_s2z_prior(
        prior, bframe, class = "Intercept", r = r_context
      )
    } else {
      specs[[i]] <- re_s2z_prior(
        prior, bframe, class = "b", coef = qnames[i], r = r_context
      )
    }
  }
  lapply(infos, function(info) {
    info$prior <- specs
    info
  })
}

# Describe all local S2Z blocks in one linear predictor. Blocks retain their
# own covariance models but share one active-coordinate omitted-mean system.
re_s2z_infos <- function(bframe, prior = NULL) {
  stopifnot(is.bframel(bframe))
  if (!has_re_s2z(bframe)) {
    return(list())
  }
  if (is.list(bframe$frame$re_s2z) && length(bframe$frame$re_s2z)) {
    out <- bframe$frame$re_s2z
  } else {
    out <- .re_s2z_descriptors(bframe)
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
  if (is_ordinal(bframe$family)) {
    stop_re_s2z(
      re_s2z_context(bframe, r = first_r), "ordinal_location",
      paste0(
        "Ordinal thresholds are not yet supported together with ",
        "gr(..., s2z = TRUE)."
      ),
      "use a conventional group-level term for this ordinal predictor, or wait for ordinal S2Z support."
    )
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
      stop_re_s2z(
        context, capability,
        paste0(
          "Arguments 'by', 'cov', and 'pw' are not yet supported together ",
          "with gr(..., s2z = TRUE)."
        ),
        "remove the unsupported argument or use s2z = FALSE for this term."
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
    return(invisible(NULL))
  }
  X <- bframe$sdata$fe$X
  infos <- re_s2z_infos(bframe)
  missing <- unique(unlist(lapply(infos, function(info) {
    info$r$coef[is.na(info$match_q)]
  }), use.names = FALSE))
  if (length(missing)) {
    first <- which(vapply(infos, function(info) {
      any(is.na(info$match_q))
    }, logical(1)))[1L]
    stop_re_s2z(
      re_s2z_context(bframe, r = infos[[first]]$r, coef = missing),
      "matching_name",
      paste0(
        "Every sum-to-zero group-level coefficient must have a matching ",
        "population-level design column. Missing: ",
        paste(missing, collapse = ", "), "."
      ),
      "add population-level terms with these coefficient names or remove the unmatched group effects."
    )
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
  invisible(NULL)
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
      first <- occurrences[[1L]]
      r <- subset2(first$frame$re, id = id)
      stop_re_s2z(
        re_s2z_context(first, r = r), "cross_predictor_id",
        paste0(
          "A sum-to-zero group-level ID cannot span multiple linear ",
          "predictors."
        ),
        "assign a distinct group-level ID within each linear predictor."
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
    active_names <- infos[[1L]]$qnames[active_q]
    active_intercept <- infos[[1L]]$center & active_names == "Intercept"
    active_b <- active_names[!active_intercept]
    if (length(active_b) && has_special_prior(prior, x, class = "b")) {
      active_b_q <- active_q[!active_intercept]
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
    if (any(active_intercept) &&
        has_special_prior(prior, x, class = "Intercept")) {
      stop_re_s2z(
        re_s2z_context(
          x, r = infos[[1L]]$r, coef = "Intercept",
          prior = "special Intercept prior"
        ),
        "active_special_prior",
        paste0(
          "Special population-level priors are not yet supported on ",
          "S2Z-active coordinates."
        ),
        "use a supported scalar intercept prior or remove the active intercept map."
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
