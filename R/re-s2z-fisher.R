# Analysis helpers for Stan-side Fisher centering of S2Z nonlinear parameters

.stop_re_s2z_fisher_nl <- function(...) {
  stop2("Stan-side Fisher centering for nonlinear S2Z effects ", ...)
}

# Does a nonlinear parameter contain coordinates that vary below the
# population-coefficient level? Such a parameter cannot enter the coordinate
# map of another latent block without making its restricted Jacobian depend on
# those latent coordinates.
.re_s2z_fisher_nlpar_is_latent <- function(x) {
  stopifnot(is.bframel(x))
  # Fixed offsets are population-only known quantities, not latent
  # coordinates. Mirror has_special_terms() for every genuinely varying term
  # while deliberately leaving offsets out of this classification.
  varying_terms <- c("sp", "sm", "gp", "ac", "cs")
  NROW(x[["re"]]) > 0L || any(lengths(x[varying_terms]))
}

# Validate the deliberately small expression language emitted by R's symbolic
# differentiator. Keeping this explicit makes unsupported outer nonlinear
# formulas fail during brms code generation rather than later in stanc.
.validate_re_s2z_fisher_derivative <- function(x) {
  if (is.numeric(x) && length(x) == 1L || is.name(x)) {
    return(invisible(NULL))
  }
  if (!is.call(x)) {
    .stop_re_s2z_fisher_nl("produced an unsupported derivative expression.")
  }
  fun <- as.character(x[[1L]])
  allowed <- c(
    "(", "+", "-", "*", "/", "^", "exp", "log", "sqrt",
    "sin", "cos", "tan", "sinh", "cosh", "tanh"
  )
  if (length(fun) != 1L || !fun %in% allowed) {
    .stop_re_s2z_fisher_nl(
      "produced unsupported derivative function '", fun, "'."
    )
  }
  lapply(as.list(x[-1L]), .validate_re_s2z_fisher_derivative)
  invisible(NULL)
}

# Recursively replace symbols in a derivative expression. Replacement values
# are language objects so the result can be rendered directly as Stan syntax.
.replace_re_s2z_fisher_symbols <- function(x, replacements) {
  if (is.name(x)) {
    key <- as.character(x)
    if (key %in% names(replacements)) {
      return(replacements[[key]])
    }
    return(x)
  }
  if (!is.call(x)) {
    return(x)
  }
  as.call(lapply(as.list(x), .replace_re_s2z_fisher_symbols,
                 replacements = replacements))
}

# Describe a population-only nonlinear predictor used by the outer derivative.
# The generated vector expression intentionally matches stan_fe's ordinary,
# dense, non-QR code path.
.re_s2z_fisher_dependency_info <- function(bframe, nlpar, id) {
  dep <- bframe$nlpars[[nlpar]]
  if (!is.bframel(dep)) {
    .stop_re_s2z_fisher_nl(
      "requires nonlinear dependency '", nlpar,
      "' to have a linear predictor."
    )
  }
  if (.re_s2z_fisher_nlpar_is_latent(dep)) {
    .stop_re_s2z_fisher_nl(
      "has a derivative depending on latent or group-varying nonlinear ",
      "parameter '", nlpar, "'."
    )
  }
  fe <- dep$frame$fe
  if (isTRUE(fe$center)) {
    .stop_re_s2z_fisher_nl(
      "does not yet support centered population-level predictors for ",
      "derivative dependency '", nlpar, "'."
    )
  }
  if (isTRUE(fe$sparse) || !identical(fe$decomp, "none")) {
    .stop_re_s2z_fisher_nl(
      "does not yet support sparse or QR predictors for derivative ",
      "dependency '", nlpar, "'."
    )
  }
  has_offset <- is.formula(dep$offset)
  if (!length(fe$vars_stan) && !has_offset) {
    .stop_re_s2z_fisher_nl(
      "requires derivative dependency '", nlpar,
      "' to contain population-level coefficients or a fixed offset."
    )
  }
  p <- usc(combine_prefix(check_prefix(dep)))
  # Qualify the private alias by the S2Z ID. The same population predictor may
  # supply loadings to more than one Fisher-centered factor in one model.
  vector_name <- paste0("fisher_s2z_", id, "_nlp", p)
  x_name <- paste0("X", p)
  coefficient_name <- paste0("b", p)
  vector_expression <- if (length(fe$vars_stan)) {
    paste(x_name, "*", coefficient_name)
  } else {
    paste0("offsets", p)
  }
  if (length(fe$vars_stan) && has_offset) {
    vector_expression <- paste0(
      "(", vector_expression, ") + offsets", p
    )
  }
  nlist(
    nlpar, prefix = p, frame = dep, x_name, coefficient_name,
    vector_name,
    vector_expression,
    observation_expression = paste0(vector_name, "[n]")
  )
}

# Describe the rows of one strict latent-score ID in a local nonlinear
# predictor. Conventional S2Z descriptors also encode an omitted population
# mean and therefore must not be used for strict latent scores.
.re_s2z_fisher_strict_rows <- function(bfl, id = NULL) {
  stopifnot(is.bframel(bfl))
  r_all <- bfl$frame$re
  if (!is.reframe(r_all) || !has_rows(r_all)) {
    .stop_re_s2z_fisher_nl(
      "requires a strict latent-score group-level block."
    )
  }
  is_strict <- r_all$s2z & re_s2z_latent(r_all)
  strict_ids <- unique(r_all$id[is_strict])
  if (is.null(id)) {
    if (!length(strict_ids)) {
      .stop_re_s2z_fisher_nl(
        "requires gr(..., s2z = TRUE, latent = TRUE, ",
        "center = \"fisher\")."
      )
    }
    if (length(strict_ids) != 1L) {
      .stop_re_s2z_fisher_nl(
        "requires an explicit group-level ID when a nonlinear parameter ",
        "contains multiple strict latent-score blocks."
      )
    }
    id <- strict_ids[[1L]]
  }
  r <- subset2(r_all, id = id)
  if (!has_rows(r) || !all(r$s2z) || !all(re_s2z_latent(r)) ||
      !identical(re_s2z_center_mode(r), "auto")) {
    .stop_re_s2z_fisher_nl(
      "requires gr(..., s2z = TRUE, latent = TRUE, ",
      "center = \"fisher\")."
    )
  }
  nlist(id = unique(r$id)[[1L]], r)
}

# Construct metadata for one Fisher-centered strict latent-score block living
# in a nonlinear response-location parameter. Both the outer derivative and
# its likelihood reference are made independent of the latent contrasts. Set
# strict = FALSE to obtain an unsupported-result object instead of throwing
# the diagnostic.
re_s2z_fisher_nl_info <- function(bframe, bfl, id = NULL, strict = TRUE) {
  strict <- as_one_logical(strict)
  analyze <- function() {
    if (!is.brmsframe(bframe) || !is.bframel(bfl)) {
      .stop_re_s2z_fisher_nl(
        "requires a univariate brms frame and a local linear predictor."
      )
    }
    info <- .re_s2z_fisher_strict_rows(bfl, id = id)
    target_nlpar <- unique(info$r$nlpar)
    if (length(target_nlpar) != 1L || !nzchar(target_nlpar)) {
      .stop_re_s2z_fisher_nl(
        "currently requires the S2Z block to belong to one nonlinear ",
        "parameter."
      )
    }
    if (!target_nlpar %in% names(bframe$nlpars) ||
        !identical(bframe$nlpars[[target_nlpar]], bfl)) {
      .stop_re_s2z_fisher_nl(
        "could not match nonlinear parameter '", target_nlpar,
        "' to its local predictor."
      )
    }
    family <- bframe$family
    response_family <- family$family
    is_gaussian_identity <- identical(response_family, "gaussian") &&
      identical(family$link, "identity")
    is_student_identity <- identical(response_family, "student") &&
      identical(family$link, "identity")
    is_location_identity <- is_gaussian_identity || is_student_identity
    outer <- bframe$dpars[["mu"]]
    if (!is.bframenl(outer) || !target_nlpar %in% outer$used_nlpars) {
      .stop_re_s2z_fisher_nl(
        "requires the S2Z nonlinear parameter to enter the outer mean formula."
      )
    }
    if (!isTRUE(outer$loop)) {
      .stop_re_s2z_fisher_nl(
        "currently requires an observation-local outer nonlinear formula ",
        "with loop = TRUE."
      )
    }

    outer_expression <- outer$formula[[2L]]
    .validate_re_s2z_fisher_derivative(outer_expression)
    derivative <- try(D(outer_expression, target_nlpar), silent = TRUE)
    if (inherits(derivative, "try-error")) {
      .stop_re_s2z_fisher_nl(
        "could not symbolically differentiate the outer mean with respect to '",
        target_nlpar, "'."
      )
    }
    .validate_re_s2z_fisher_derivative(derivative)
    dependencies <- all.vars(derivative)
    if (target_nlpar %in% dependencies) {
      .stop_re_s2z_fisher_nl(
        "has an outer derivative that depends on target nonlinear parameter '",
        target_nlpar, "'."
      )
    }
    nlpar_dependencies <- intersect(dependencies, names(bframe$nlpars))
    latent_nlpars <- names(Filter(
      .re_s2z_fisher_nlpar_is_latent, bframe$nlpars
    ))
    bad_latent <- intersect(nlpar_dependencies, latent_nlpars)
    if (length(bad_latent)) {
      .stop_re_s2z_fisher_nl(
        "has a derivative depending on latent or group-varying nonlinear ",
        "parameter '", bad_latent[1L], "'."
      )
    }
    covariates <- all.vars(outer$covars)
    unknown <- setdiff(dependencies, c(names(bframe$nlpars), covariates))
    if (length(unknown)) {
      .stop_re_s2z_fisher_nl(
        "has unsupported derivative dependency '", unknown[1L], "'."
      )
    }

    # The response-family information is evaluated at the population-only
    # nonlinear predictor. Retain fixed effects and offsets of every nonlinear
    # parameter, while removing all latent and group-varying contrasts. This
    # reference may depend on sampled population parameters, but never on the
    # strict score coordinates whose restricted Jacobian is used below.
    outer_dependencies <- all.vars(outer_expression)
    outer_nlpar_dependencies <- intersect(
      outer_dependencies, names(bframe$nlpars)
    )
    outer_covariate_dependencies <- intersect(outer_dependencies, covariates)
    outer_unknown <- setdiff(
      outer_dependencies, c(names(bframe$nlpars), covariates)
    )
    if (length(outer_unknown)) {
      .stop_re_s2z_fisher_nl(
        "has unsupported outer-reference dependency '",
        outer_unknown[1L], "'."
      )
    }
    population_nlpar_dependencies <- setdiff(
      outer_nlpar_dependencies, latent_nlpars
    )
    dependency_info <- lapply(population_nlpar_dependencies, function(nlpar) {
      .re_s2z_fisher_dependency_info(bframe, nlpar, id = info$id)
    })
    names(dependency_info) <- population_nlpar_dependencies
    population_replacements <- lapply(dependency_info, function(x) {
      str2lang(x$observation_expression)
    })
    outer_p <- usc(combine_prefix(check_prefix(outer)))
    covariate_replacements <- list()
    for (covar in outer_covariate_dependencies) {
      covar_index <- match(covar, covariates)
      covariate_replacements[[covar]] <- str2lang(
        paste0("C", outer_p, "_", covar_index, "[n]")
      )
    }
    derivative_replacements <- c(
      population_replacements[nlpar_dependencies],
      covariate_replacements[intersect(dependencies, covariates)]
    )
    derivative_stan <- .replace_re_s2z_fisher_symbols(
      derivative, derivative_replacements
    )
    derivative_stan <- deparse0(derivative_stan)
    latent_replacements <- lapply(
      intersect(outer_nlpar_dependencies, latent_nlpars),
      function(nlpar) {
        str2lang(stan_re_s2z_fisher_reference_eta(
          bframe$nlpars[[nlpar]], n = "n"
        ))
      }
    )
    names(latent_replacements) <- intersect(
      outer_nlpar_dependencies, latent_nlpars
    )
    outer_reference_replacements <- c(
      population_replacements, latent_replacements,
      covariate_replacements
    )
    outer_reference_stan <- .replace_re_s2z_fisher_symbols(
      outer_expression, outer_reference_replacements
    )
    outer_reference_stan <- deparse0(outer_reference_stan)
    covariate_dependencies <- intersect(dependencies, covariates)
    response_addition_terms <- names(Filter(is.formula, bframe$adforms))
    has_response_addition_terms <- length(response_addition_terms) > 0L
    sigma_dpar <- bframe$dpars[["sigma"]]
    sigma_fdpar <- bframe$fdpars[["sigma"]]
    sigma_kind <- if (!is.null(sigma_dpar)) {
      "predictor"
    } else if (!is.null(sigma_fdpar)) {
      "fixed"
    } else {
      "scalar_parameter"
    }
    sigma_is_scalar <- !identical(sigma_kind, "predictor")
    sigma_fixed_value <- if (identical(sigma_kind, "fixed")) {
      sigma_fdpar$value
    }
    sigma_parameter <- paste0("sigma", usc(bframe$resp))
    M <- nrow(info$r)
    cn <- unname(info$r$cn)
    frame_rows <- rownames(info$r)
    row_keys <- paste(info$r$id, info$r$nlpar, info$r$cn, sep = ":")
    row_metadata <- data.frame(
      frame_row = frame_rows,
      row_key = row_keys,
      id = unname(info$r$id),
      nlpar = unname(info$r$nlpar),
      coef = unname(info$r$coef),
      cn,
      stringsAsFactors = FALSE
    )
    structure(
      nlist(
        supported = TRUE, reason_unsupported = NULL, id = info$id,
        target_nlpar, outer_expression, derivative,
        outer_reference_stan, outer_dependencies,
        outer_nlpar_dependencies, outer_covariate_dependencies,
        obs_derivative_expr = derivative,
        obs_derivative_stan = derivative_stan,
        dzeta_dxi = rep(list(derivative), M),
        dzeta_dxi_stan = rep(derivative_stan, M),
        dependencies, nlpar_dependencies, covariate_dependencies,
        dependency_info, latent_nlpars, response_family,
        is_gaussian_identity, is_student_identity, is_location_identity,
        is_observation_local = TRUE, score_independent = TRUE,
        response_addition_terms, has_response_addition_terms,
        sigma_kind, sigma_is_scalar, sigma_fixed_value, sigma_parameter,
        coefficients = info$r$coef, cn, frame_rows, row_keys, row_metadata, M
      ),
      class = "re_s2z_fisher_nl_info"
    )
  }

  out <- tryCatch(analyze(), error = identity)
  if (!inherits(out, "error")) {
    return(out)
  }
  if (strict) {
    stop2(conditionMessage(out))
  }
  structure(
    list(supported = FALSE, reason_unsupported = conditionMessage(out)),
    class = "re_s2z_fisher_nl_info"
  )
}
