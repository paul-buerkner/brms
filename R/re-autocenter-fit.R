# Two-stage fitting helpers for automatic group-effect centering.

# Split automatic group-effect rows into the original gr() occurrences.  A
# resolved occurrence carries its fixed matrix on one of its expanded rows;
# frame_re() deliberately removes duplicate copies from the remaining rows.
re_autocenter_occurrences <- function(r) {
  stopifnot(is.reframe(r))
  if (!has_rows(r)) {
    return(list())
  }
  take_auto <- which(re_s2z_center_auto(r))
  if (!length(take_auto)) {
    return(list())
  }
  key <- paste(
    r$id, r$resp, r$dpar, r$nlpar, r$gn,
    sep = "\r"
  )
  key <- key[take_auto]
  groups <- lapply(unique(key), function(value) {
    take_auto[key == value]
  })
  lapply(groups, function(take) {
    supplied <- vapply(r$gcall[take], function(gcall) {
      !is.null(gcall[["s2z_center_data"]])
    }, logical(1))
    list(
      take = take,
      id = r$id[take[1L]],
      resp = r$resp[take[1L]],
      dpar = r$dpar[take[1L]],
      nlpar = r$nlpar[take[1L]],
      gn = r$gn[take[1L]],
      resolved = any(supplied)
    )
  })
}

# Return one descriptor per unresolved automatic covariance block.  Resolved
# matrices retain the internal auto marker so regenerated Stan code is
# identical, but their gr() calls already carry fixed data and need no pilot.
re_autocenter_specs <- function(bframe, unresolved = TRUE) {
  stopifnot(is.anybrmsframe(bframe))
  r_all <- bframe$frame$re
  if (!has_rows(r_all)) {
    return(list())
  }
  ids <- unique(r_all$id[re_s2z_center_auto(r_all)])
  out <- list()
  for (id in ids) {
    r <- subset2(r_all, id = id)
    occurrences <- re_autocenter_occurrences(r)
    resolved <- vapply(occurrences, `[[`, logical(1), "resolved")
    if (isTRUE(unresolved)) {
      if (length(resolved) && all(resolved)) {
        next
      }
      if (any(resolved)) {
        stop2(
          "Group-level ID ", id, " mixes resolved and unresolved automatic ",
          "centering occurrences. Request center = \"auto\" for every ",
          "occurrence sharing the ID, or supply fixed numeric matrices for ",
          "all of them."
        )
      }
    }
    levels <- get_levels(r)[[r$group[1L]]]
    latent <- all(r$s2z) && all(re_s2z_latent(r))
    cross <- all(r$s2z) && !latent && is_re_s2z_cross_id(bframe, id)
    if (latent) {
      dimension <- re_s2z_latent_dimension(r)
      M <- max(dimension)
      representatives <- match(seq_len(M), dimension)
      coefficients <- paste(
        r$nlpar[representatives], r$coef[representatives], sep = ":"
      )
    } else {
      dimension <- as.integer(r$cn)
      representatives <- seq_rows(r)
      M <- nrow(r)
      if (cross) {
        prefix <- combine_prefix(check_prefix(r))
        coefficients <- paste0(
          prefix, ":", r$coef, "[", dimension, "]"
        )
        coefficients <- coefficients[match(seq_len(M), dimension)]
      } else {
        coefficients <- r$coef[match(seq_len(M), dimension)]
      }
    }
    if (length(coefficients) != M || anyNA(coefficients)) {
      stop2("Internal mismatch in automatic group-centering dimensions.")
    }
    out[[as.character(id)]] <- nlist(
      id, r, levels, coefficients, dimension, latent, cross,
      G = length(levels), M
    )
  }
  out
}

has_re_autocenter <- function(bframe, unresolved = TRUE) {
  length(re_autocenter_specs(bframe, unresolved = unresolved)) > 0L
}

# Recover every fully resolved automatic matrix from a framed formula. This
# preserves centering_weights() metadata when an otherwise unrelated update
# recompiles a model whose automatic weights were already frozen.
resolved_re_autocenter_weights <- function(bframe) {
  stopifnot(is.anybrmsframe(bframe))
  specs <- re_autocenter_specs(bframe, unresolved = FALSE)
  if (!length(specs)) {
    return(list())
  }
  out <- list()
  for (id in names(specs)) {
    spec <- specs[[id]]
    occurrences <- re_autocenter_occurrences(spec$r)
    resolved <- vapply(occurrences, `[[`, logical(1), "resolved")
    if (!length(resolved) || !all(resolved)) {
      next
    }
    rho <- re_center_data_matrix(spec$r, levels = spec$levels)
    representatives <- match(seq_len(spec$M), spec$dimension)
    if (anyNA(representatives)) {
      stop2("Internal mismatch in resolved automatic centering dimensions.")
    }
    rho <- rho[, representatives, drop = FALSE]
    dimnames(rho) <- list(spec$levels, spec$coefficients)
    out[[id]] <- .as_brmsautocenter_resolved(rho)
  }
  out
}

# Name and validate matrices returned by the generic draw aggregator.
align_re_autocenter_summary <- function(summary, specs) {
  stopifnot(
    inherits(summary, "brms_re_autocenter_summary"), is.list(specs)
  )
  if (!setequal(names(summary$rho), names(specs))) {
    stop2("The precursor did not return every requested centering matrix.")
  }
  summary$rho <- lapply(names(specs), function(id) {
    spec <- specs[[id]]
    rho <- summary$rho[[id]]
    if (!identical(dim(rho), c(spec$G, spec$M))) {
      stop2("The precursor returned an invalid centering matrix for ID ", id,
            ".")
    }
    dimnames(rho) <- list(spec$levels, spec$coefficients)
    .as_brmsautocenter_resolved(rho)
  })
  names(summary$rho) <- names(specs)
  summary
}

# Convert full covariance-block matrices to one matrix per gr(center="auto")
# occurrence, in the same traversal order used by replace_re_autocenter().
re_autocenter_occurrence_weights <- function(bframe, weights, formula,
                                              aggregate = "median") {
  stopifnot(
    is.anybrmsframe(bframe), is.list(weights),
    is.formula(formula) || is.brmsformula(formula) ||
      is.mvbrmsformula(formula)
  )
  aggregate <- match.arg(aggregate, c("median", "mean"))
  r_all <- bframe$frame$re
  occurrences <- re_autocenter_occurrences(r_all)
  unresolved <- !vapply(occurrences, `[[`, logical(1), "resolved")
  occurrences <- occurrences[unresolved]
  if (!length(occurrences)) {
    return(list(center = list(), id = integer(), occurrence = integer()))
  }
  unresolved_ids <- unique(as.integer(vapply(
    occurrences, `[[`, numeric(1), "id"
  )))
  weight_ids <- names(weights)
  if (is.null(weight_ids) || !setequal(weight_ids, as.character(unresolved_ids))) {
    stop2("Automatic centering weights do not match the unresolved group IDs.")
  }

  center <- vector("list", length(occurrences))
  ids <- integer(length(occurrences))
  records <- data.frame(
    resp = character(length(occurrences)),
    dpar = character(length(occurrences)),
    nlpar = character(length(occurrences)),
    gn = integer(length(occurrences)),
    stringsAsFactors = FALSE
  )
  for (i in seq_along(occurrences)) {
    occurrence_info <- occurrences[[i]]
    r <- r_all[occurrence_info$take, , drop = FALSE]
    id <- as.character(occurrence_info$id)
    rho <- weights[[id]]
    if (is.null(rho)) {
      stop2("Missing automatic centering weights for group-level ID ", id,
            ".")
    }
    if (all(r$s2z) && all(re_s2z_latent(r))) {
      r_id <- subset2(r_all, id = r$id[1L])
      columns <- re_s2z_latent_dimension(r_id)[
        match(rownames(r), rownames(r_id))
      ]
      # Reframes do not guarantee informative row names. Fall back to the
      # stable latent key when necessary.
      if (anyNA(columns)) {
        columns <- match(
          re_s2z_latent_key(r),
          unique(re_s2z_latent_key(r_id))
        )
      }
    } else {
      columns <- as.integer(r$cn)
    }
    if (anyNA(columns) || any(columns < 1L | columns > ncol(rho))) {
      stop2("Internal mismatch while slicing automatic centering weights.")
    }
    value <- rho[, columns, drop = FALSE]
    colnames(value) <- r$coef
    center[[i]] <- .as_brmsautocenter_resolved(value)
    ids[i] <- as.integer(occurrence_info$id)
    records$resp[i] <- occurrence_info$resp
    records$dpar[i] <- occurrence_info$dpar
    records$nlpar[i] <- occurrence_info$nlpar
    records$gn[i] <- occurrence_info$gn
  }

  # terms_re() sorts grouping factors before frame IDs are assigned, whereas
  # replace_re_autocenter() necessarily walks the source formulas in their
  # written order. Recover that exact component order and then use `gn`, which
  # retains the original within-formula random-term position.
  components <- list()
  add_component <- function(resp, parameter = NULL, reserved = character()) {
    components[[length(components) + 1L]] <<-
      nlist(resp, parameter, reserved)
  }
  walk_components <- function(object, multivariate = FALSE) {
    if (is.mvbrmsformula(object)) {
      for (form in object$forms) {
        walk_components(form, multivariate = TRUE)
      }
      return(invisible(NULL))
    }
    if (is.brmsformula(object)) {
      resp <- if (multivariate) as.character(object$resp) else ""
      # An explicit mu formula replaces the main RHS in brmsterms(); without
      # one, the main formula is the (prefix-free) mu component.
      if (!"mu" %in% names(object$pforms)) {
        add_component(resp, reserved = names(object$pforms))
      }
      for (parameter in names(object$pforms)) {
        add_component(resp, parameter)
      }
      return(invisible(NULL))
    }
    if (is.formula(object)) {
      add_component("")
      return(invisible(NULL))
    }
    stop2("Automatic centering requires a formula or brms formula object.")
  }
  walk_components(formula)
  remaining <- seq_len(nrow(records))
  formula_groups <- list()
  for (component in components) {
    if (!length(remaining)) {
      break
    }
    same_response <- records$resp[remaining] == component$resp
    if (is.null(component$parameter)) {
      parameter <- ifelse(
        nzchar(records$dpar[remaining]),
        records$dpar[remaining], records$nlpar[remaining]
      )
      # Multicategory shorthand replicates one main-formula term into several
      # generated mean predictors. Exclude explicit pforms, then collapse
      # records with the same original term number to their one source call.
      same_parameter <- !parameter %in% component$reserved
    } else if (identical(component$parameter, "mu")) {
      same_parameter <- !nzchar(records$dpar[remaining]) &
        !nzchar(records$nlpar[remaining])
    } else {
      same_parameter <- records$dpar[remaining] == component$parameter |
        records$nlpar[remaining] == component$parameter
    }
    take <- remaining[same_response & same_parameter]
    if (length(take)) {
      take <- take[order(records$gn[take])]
      gn <- unique(records$gn[take])
      formula_groups <- c(
        formula_groups,
        lapply(gn, function(value) take[records$gn[take] == value])
      )
      remaining <- setdiff(remaining, take)
    }
  }
  if (length(remaining)) {
    stop2("Could not align automatic centering weights to formula occurrences.")
  }
  collapse_group <- function(index) {
    values <- center[index]
    template <- values[[1L]]
    compatible <- vapply(values, function(value) {
      identical(dim(value), dim(template)) &&
        identical(dimnames(value), dimnames(template))
    }, logical(1))
    if (!all(compatible)) {
      stop2(
        "Automatically expanded predictors require compatible centering ",
        "matrices for their shared source gr() term."
      )
    }
    if (length(values) == 1L) {
      return(template)
    }
    array <- simplify2array(lapply(values, unclass))
    if (identical(aggregate, "median")) {
      out <- apply(array, c(1L, 2L), stats::median)
    } else {
      out <- apply(array, c(1L, 2L), mean)
    }
    dim(out) <- dim(template)
    dimnames(out) <- dimnames(template)
    .as_brmsautocenter_resolved(out)
  }
  center <- lapply(formula_groups, collapse_group)
  ids <- vapply(
    formula_groups, function(index) as.integer(ids[index[1L]]), integer(1)
  )
  occurrence <- ave(
    seq_along(ids), ids, FUN = function(index) seq_along(index)
  )
  nlist(center, id = ids, occurrence = as.integer(occurrence))
}

# Fit a fully noncentered precursor and retain only its centering proposals.
# Pathfinder and HMC share the same contract; an HMC precursor is a separate
# short run and never mutates the coordinates of an active warmup.
.validate_re_autocenter_integer <- function(x, name, minimum = 1L) {
  x <- as_one_numeric(x)
  if (!is.finite(x) || x != floor(x) || x < minimum ||
      x > .Machine$integer.max) {
    stop2("Autocenter '", name, "' must be an integer not smaller than ",
          minimum, ".")
  }
  as.integer(x)
}

run_re_autocenter_pilot <- function(model, sdata, specs, control, backend,
                                    chains, cores, threads, opencl, init,
                                    seed, silent) {
  control <- validate_re_autocenter_control(control)
  stopifnot(is.list(sdata), is.list(specs), length(specs) > 0L)
  if (!identical(backend, "cmdstanr")) {
    stop2("Automatic group centering requires backend = \"cmdstanr\".")
  }
  pilot_data <- sdata
  candidates <- character(length(specs))
  for (i in seq_along(specs)) {
    spec <- specs[[i]]
    id <- spec$id
    pilot_data[[paste0("rho_s2z_", id)]] <- matrix(0, spec$G, spec$M)
    pilot_data[[paste0("compute_rho_center_candidate_", id)]] <- 1L
    candidates[i] <- paste0("rho_center_candidate_", id)
  }

  pilot_args <- control$pilot_args %||% list()
  if (identical(control$method, "pathfinder")) {
    pathfinder_aliases <- intersect(
      names(pilot_args), c("chains", "iter", "warmup")
    )
    if (length(pathfinder_aliases)) {
      stop2("Autocenter Pathfinder 'pilot_args' cannot contain ",
            collapse_comma(pathfinder_aliases), ".")
    }
    num_paths <- pilot_args$num_paths %||% min(max(1L, chains), 4L)
    pilot_args$num_paths <- NULL
    pilot_chains <- .validate_re_autocenter_integer(
      num_paths, "pilot_args$num_paths"
    )
    # Candidate generated quantities can be comparatively expensive.  A few
    # hundred approximate draws are ample for the default median aggregation.
    # CmdStan uses separate controls for one path and PSIS-combined paths;
    # interpret a user-supplied `draws` value for both unless they explicitly
    # distinguish `single_path_draws`.
    pilot_args$draws <- pilot_args$draws %||% 200L
    pilot_args$single_path_draws <-
      pilot_args$single_path_draws %||% pilot_args$draws
  } else {
    if ("num_paths" %in% names(pilot_args)) {
      stop2("Autocenter HMC 'pilot_args' cannot contain 'num_paths'.")
    }
    pilot_chains <- pilot_args$chains %||% min(max(1L, chains), 4L)
    pilot_iter <- pilot_args$iter %||% 500L
    pilot_warmup <- pilot_args$warmup %||% 400L
    pilot_args$chains <- pilot_args$iter <- pilot_args$warmup <- NULL
    pilot_chains <- .validate_re_autocenter_integer(
      pilot_chains, "pilot_args$chains"
    )
    pilot_iter <- .validate_re_autocenter_integer(
      pilot_iter, "pilot_args$iter"
    )
    pilot_warmup <- .validate_re_autocenter_integer(
      pilot_warmup, "pilot_args$warmup", minimum = 0L
    )
    if (pilot_warmup >= pilot_iter) {
      stop2("Autocenter HMC 'pilot_args$warmup' must be smaller than ",
            "'pilot_args$iter'.")
    }
  }
  pilot_control <- pilot_args$control %||% list()
  pilot_args$control <- NULL
  if (!is.list(pilot_control) ||
      length(pilot_control) &&
      (is.null(names(pilot_control)) || any(!nzchar(names(pilot_control))) ||
       anyDuplicated(names(pilot_control)))) {
    stop2("Autocenter HMC 'pilot_args$control' must be a named list.")
  }
  if (identical(control$method, "pathfinder") && length(pilot_control)) {
    stop2("Autocenter 'pilot_args$control' is available only for method = ",
          "\"hmc\".")
  }
  overlap <- intersect(names(pilot_control), names(pilot_args))
  if (length(overlap)) {
    stop2("Autocenter HMC controls cannot also appear directly in ",
          "'pilot_args': ", collapse_comma(overlap), ".")
  }
  invariant_control <- intersect(
    names(pilot_control),
    c("data", "seed", "init", "iter_sampling", "iter_warmup", "chains",
      "parallel_chains", "thin", "fixed_param", "threads_per_chain",
      "opencl_ids", "show_messages", "show_exceptions")
  )
  if (length(invariant_control)) {
    stop2("Autocenter HMC controls cannot override precursor orchestration ",
          "arguments: ", collapse_comma(invariant_control), ".")
  }
  if (length(intersect(
    names(pilot_args),
    c("model", "data", "algorithm", "backend", "exclude", "fixed_param")
  ))) {
    stop2("Reserved arguments cannot be supplied in autocenter 'pilot_args'.")
  }
  class(pilot_data) <- "list"
  if (isNA(seed)) {
    seed <- NULL
  }
  if (is_equal(init, "random")) {
    init <- NULL
  } else if (is_equal(init, "0")) {
    init <- 0
  }
  if (is.list(init)) {
    required_init <- if (identical(control$method, "pathfinder")) {
      1L
    } else {
      pilot_chains
    }
    if (length(init) < required_init) {
      stop2("The final-fit 'init' list has fewer entries than the ",
            control$method, " precursor requires.")
    }
    init <- init[seq_len(required_init)]
  }
  args <- nlist(data = pilot_data, seed, init)
  if (use_opencl(opencl)) {
    args$opencl_ids <- opencl$ids
  }
  threaded <- use_threading(threads, force = TRUE)
  if (identical(control$method, "pathfinder")) {
    c(args) <- list(
      num_paths = pilot_chains,
      show_messages = silent < 2,
      show_exceptions = silent == 0
    )
    if (threaded) {
      args$num_threads <- threads$threads
    }
    runner <- model$pathfinder
    algorithm <- "pathfinder"
  } else {
    c(args) <- list(
      iter_sampling = pilot_iter - pilot_warmup,
      iter_warmup = pilot_warmup,
      chains = pilot_chains,
      thin = 1L,
      parallel_chains = max(1L, min(cores, pilot_chains)),
      show_messages = silent < 2,
      show_exceptions = silent == 0
    )
    if (threaded) {
      args$threads_per_chain <- threads$threads
    }
    args[names(pilot_control)] <- pilot_control
    runner <- model$sample
    algorithm <- "HMC"
  }
  args[names(pilot_args)] <- pilot_args
  invalid_args <- setdiff(names(args), names(formals(runner)))
  if (length(invalid_args)) {
    stop2("Invalid autocenter ", algorithm, " 'pilot_args': ",
          collapse_comma(invalid_args), ".")
  }
  pilot_fit <- tryCatch(
    do_call(runner, args),
    error = function(error) {
      stop2("The automatic-centering ", algorithm,
            " precursor failed: ", conditionMessage(error))
    }
  )
  draws <- tryCatch(
    pilot_fit$draws(variables = candidates, format = "matrix"),
    error = function(error) {
      stop2("Could not read automatic-centering candidates from the ",
            algorithm, " precursor: ", conditionMessage(error))
    }
  )
  summary <- aggregate_re_autocenter_draws(draws, control = control)
  align_re_autocenter_summary(summary, specs)
}

#' @rdname autocenter_control
#' @export
centering_weights.brmsfit <- function(x, ...) {
  weights <- x$autocenter$weights %||% NULL
  if (is.null(weights)) {
    stop2("This brmsfit does not contain precursor centering weights.")
  }
  weights
}

#' @rdname autocenter_control
#' @export
centering_weights.brmsfit_multiple <- function(x, ...) {
  by_fit <- x$autocenter$by_fit %||% NULL
  if (is.null(by_fit)) {
    return(centering_weights.brmsfit(x, ...))
  }
  lapply(by_fit, function(autocenter) {
    weights <- autocenter$weights %||% NULL
    if (is.null(weights)) {
      stop2("A component fit does not contain precursor centering weights.")
    }
    weights
  })
}
