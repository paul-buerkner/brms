# Helpers for the physical sum-to-zero parameterization of group effects

# Does a local linear predictor contain a sum-to-zero group-level block?
has_re_s2z <- function(x) {
  is.bframel(x) && has_rows(x$frame$re) &&
    "s2z" %in% names(x$frame$re) && any(x$frame$re$s2z)
}

# Return the strict-latent flags, defaulting to FALSE for old reframes.
re_s2z_latent <- function(r) {
  stopifnot(is.reframe(r), has_rows(r))
  value <- r[["latent"]]
  if (is.null(value)) {
    return(rep(FALSE, nrow(r)))
  }
  if (!is.logical(value) || length(value) != nrow(r) || anyNA(value)) {
    stop2("Internal error: invalid strict-latent sum-to-zero flags.")
  }
  value
}

# Strict latent-score occurrences with the same nonlinear-parameter and
# coefficient names denote one shared score dimension, even when they appear
# in several response predictors. This is deliberately different from the
# usual same-ID semantics, which create distinct correlated group effects.
# Distinct nonlinear parameter names (eta1, eta2, ...) define the columns of a
# multivariate latent score.
re_s2z_latent_key <- function(r) {
  stopifnot(is.reframe(r), has_rows(r), all(re_s2z_latent(r)))
  paste(r$nlpar, r$coef, sep = "\r")
}

# Map every response-local occurrence to its shared covariance column.
re_s2z_latent_dimension <- function(r) {
  key <- re_s2z_latent_key(r)
  match(key, unique(key))
}

# One representative row per strict latent covariance dimension.
re_s2z_latent_dimensions <- function(r) {
  key <- re_s2z_latent_key(r)
  r[!duplicated(key), , drop = FALSE]
}

# Identify the covariance coordinates represented by a possibly mixed
# reframe. Ordinary group-level occurrences each retain their own coordinate.
# Strict latent occurrences sharing an ID and a nonlinear-parameter/
# coefficient key instead alias one coordinate across response predictors.
# The keys are stable under subsetting so prediction code can map a
# response-local occurrence back to the covariance parameter named from the
# full fitted reframe.
re_s2z_covariance_key <- function(r) {
  stopifnot(is.reframe(r))
  if (!has_rows(r)) {
    return(character())
  }
  key <- paste(
    "occurrence", r$id, r$resp, r$dpar, r$nlpar, r$cn, r$coef,
    sep = "\r"
  )
  for (id in unique(r$id)) {
    rows <- which(r$id == id)
    r_id <- r[rows, , drop = FALSE]
    if (all(r_id$s2z) && all(re_s2z_latent(r_id))) {
      key[rows] <- paste(
        "strict", id, re_s2z_latent_key(r_id), sep = "\r"
      )
    }
  }
  key
}

# Map every occurrence to the corresponding covariance coordinate.
re_s2z_covariance_dimension <- function(r, covariance_r = r) {
  stopifnot(is.reframe(r), is.reframe(covariance_r))
  key <- re_s2z_covariance_key(r)
  covariance_key <- re_s2z_covariance_key(
    re_s2z_covariance_dimensions(covariance_r)
  )
  match(key, covariance_key)
}

# Keep one representative row for each actual covariance coordinate while
# preserving the order in which those coordinates first occur.
re_s2z_covariance_dimensions <- function(r) {
  key <- re_s2z_covariance_key(r)
  r[!duplicated(key), , drop = FALSE]
}

# Does a local predictor contain conventional (omitted-mean) S2Z effects?
has_re_s2z_conventional <- function(x) {
  if (!has_re_s2z(x)) {
    return(FALSE)
  }
  r <- x$frame$re
  any(r$s2z & !re_s2z_latent(r))
}

# Return validated centering fractions for one group-effect block. Logical
# values from reframes created before partial centering was added remain valid
# endpoints. A missing column means the historical default: centered for S2Z
# and non-centered for ordinary group effects.
re_s2z_center_values <- function(r) {
  stopifnot(is.reframe(r), has_rows(r))
  value <- r[["s2z_center"]]
  if (is.null(value)) {
    return(as.numeric(r$s2z))
  }
  if (!is.numeric(value) && !is.logical(value)) {
    stop2("Internal error: invalid group-effect centering fractions.")
  }
  value <- as.numeric(value)
  if (length(value) != nrow(r) || anyNA(value) ||
      any(!is.finite(value)) || any(value < 0 | value > 1)) {
    stop2("Internal error: group-effect centering fractions must be in [0, 1].")
  }
  value
}

# Return the per-coefficient automatic-precursor flags, defaulting to FALSE
# for old reframes.
re_s2z_center_auto <- function(r) {
  stopifnot(is.reframe(r), has_rows(r))
  value <- r[["s2z_center_auto"]]
  if (is.null(value)) {
    return(rep(FALSE, nrow(r)))
  }
  if (!is.logical(value) || length(value) != nrow(r) || anyNA(value)) {
    stop2("Internal error: invalid automatic group-effect centering flags.")
  }
  value
}

# Does a group-effect block carry user-supplied level-specific centering
# fractions? Scalar fractions remain in the historical per-coefficient
# metadata; vectors and matrices are retained on the original gr() call until
# the grouping levels and coefficient columns are known.
re_center_has_data <- function(r) {
  stopifnot(is.reframe(r), has_rows(r))
  if (!"gcall" %in% names(r) || !is.list(r$gcall)) {
    return(FALSE)
  }
  any(vapply(r$gcall, function(x) {
    !is.null(x[["s2z_center_data"]])
  }, logical(1)))
}

# Align one slide-style level vector or its multivariate extension to a framed
# group-effect occurrence. The public `center` values are centering fractions
# rho = 1 - w, where w is the non-centering weight in the tutorial slides.
# A vector supplies one value per group level and is shared by the occurrence's
# coefficient columns. A matrix supplies level-by-coefficient values; a
# one-column matrix is shared across those columns.
.re_center_data_matrix <- function(value, levels, coef) {
  stopifnot(is.numeric(value), is.character(levels), is.character(coef))
  resolved_auto <- inherits(value, "brmsautocenter_resolved")
  if (!length(value) || anyNA(value) || any(!is.finite(value)) ||
      any(value < 0 | value > 1)) {
    stop2("Level-specific 'center' values must be finite numbers in [0, 1].")
  }
  if (is.null(dim(value))) {
    value_names <- names(value)
    if (is.null(value_names)) {
      if (length(value) != length(levels)) {
        stop2("An unnamed level-specific 'center' vector must have one ",
              "value per grouping level (expected ", length(levels),
              ", found ", length(value), ").")
      }
      aligned <- as.numeric(value)
    } else {
      if (any(!nzchar(value_names)) || anyDuplicated(value_names) ||
          !setequal(value_names, levels)) {
        stop2("Names of a level-specific 'center' vector must match the ",
              "grouping levels exactly.")
      }
      aligned <- as.numeric(value[match(levels, value_names)])
    }
    return(matrix(
      rep(aligned, length(coef)), nrow = length(levels), ncol = length(coef),
      dimnames = list(levels, coef)
    ))
  }

  if (!is.matrix(value) || !all(dim(value) > 0L)) {
    stop2("Level-specific 'center' values must be a numeric vector or matrix.")
  }
  value <- as.matrix(value)
  value_levels <- rownames(value)
  if (is.null(value_levels)) {
    if (nrow(value) != length(levels)) {
      stop2("An unnamed level-by-coefficient 'center' matrix must have one ",
            "row per grouping level (expected ", length(levels),
            ", found ", nrow(value), ").")
    }
  } else {
    invalid_levels <- any(!nzchar(value_levels)) ||
      anyDuplicated(value_levels)
    if (resolved_auto) {
      if (invalid_levels || !all(levels %in% value_levels)) {
        stop2(
          "Frozen automatic-centering weights do not contain every current ",
          "grouping level. Request center = \"auto\" again to run a new ",
          "precursor for changed levels."
        )
      }
    } else if (invalid_levels || !setequal(value_levels, levels)) {
      stop2("Row names of a level-by-coefficient 'center' matrix must match ",
            "the grouping levels exactly.")
    }
    value <- value[match(levels, value_levels), , drop = FALSE]
  }

  if (ncol(value) == 1L) {
    value <- matrix(
      rep(value[, 1L], length(coef)),
      nrow = length(levels), ncol = length(coef),
      dimnames = list(levels, coef)
    )
  } else {
    value_coef <- colnames(value)
    if (is.null(value_coef)) {
      if (ncol(value) != length(coef)) {
        stop2("An unnamed level-by-coefficient 'center' matrix must have one ",
              "column per group-level coefficient (expected ", length(coef),
              ", found ", ncol(value), ").")
      }
    } else {
      invalid_coef <- any(!nzchar(value_coef)) ||
        anyDuplicated(value_coef) || anyDuplicated(coef)
      if (resolved_auto) {
        if (invalid_coef || !all(coef %in% value_coef)) {
          stop2(
            "Frozen automatic-centering weights do not contain every ",
            "current group-level coefficient. Request center = \"auto\" ",
            "again to run a new precursor for changed coefficients."
          )
        }
      } else if (invalid_coef || !setequal(value_coef, coef)) {
        stop2("Column names of a level-by-coefficient 'center' matrix must ",
              "match the group-level coefficients exactly.")
      }
      value <- value[, match(coef, value_coef), drop = FALSE]
    }
    dimnames(value) <- list(levels, coef)
  }
  value
}

# Return the full level-by-coefficient fixed centering map for one covariance
# block. Separate gr() occurrences sharing an ID may each supply a level vector
# or a matrix for their own coefficient columns.
re_center_data_matrix <- function(r, levels) {
  stopifnot(is.reframe(r), has_rows(r), is.character(levels))
  rho <- matrix(
    rep(re_s2z_center_values(r), each = length(levels)),
    nrow = length(levels), ncol = nrow(r),
    dimnames = list(levels, r$coef)
  )
  if (!re_center_has_data(r)) {
    return(rho)
  }
  # `gn` identifies the shared covariance block, not a particular gr()
  # occurrence. A strict latent ID may intentionally collect terms from
  # several nonlinear predictors, each with its own vector or matrix.
  occurrence_label <- vapply(r$gcall, function(x) {
    as.character(x$label %||% "")
  }, character(1))
  occurrence_key <- paste(
    r$resp, r$dpar, r$nlpar, occurrence_label, sep = "\r"
  )
  occurrences <- split(seq_rows(r), occurrence_key)
  for (take in occurrences) {
    values <- lapply(r$gcall[take], `[[`, "s2z_center_data")
    supplied <- !vapply(values, is.null, logical(1))
    if (!any(supplied)) {
      next
    }
    supplied_values <- values[supplied]
    if (length(supplied_values) > 1L &&
        !all(vapply(supplied_values[-1L], identical, logical(1),
                    supplied_values[[1L]]))) {
      stop2("All coefficients from one group-level term must use the same ",
            "level-specific 'center' data.")
    }
    rho[, take] <- .re_center_data_matrix(
      supplied_values[[1L]], levels = levels, coef = r$coef[take]
    )
  }
  rho
}

# Build the population-location map used by the slide-style centering chart.
# If X is the raw population design and Z is the raw varying-coefficient
# design, each exactly representable column of X is written as Z %*% C.  This
# is more general than matching coefficient names: for example, treatment
# coding in X and cell-mean coding in Z imply locations (beta_0,
# beta_0 + beta_1). Population columns outside span(Z) have no corresponding
# group-level location and therefore receive a zero column in C.
.re_center_mean_matrix <- function(r, bframe, data) {
  stopifnot(
    is.reframe(r), has_rows(r), is.bframel(bframe), is.data.frame(data)
  )
  px <- check_prefix(bframe)
  p <- usc(combine_prefix(px))
  X <- bframe$sdata$fe[[paste0("X", p)]]
  if (is.null(X)) {
    stop2("Internal mismatch in the population design for group centering.")
  }
  X <- as.matrix(X)
  Z <- .re_s2z_design_matrix_r(r, data = data)
  if (nrow(X) != nrow(Z)) {
    stop2("Internal row mismatch in the population and group designs for ",
          "group centering.")
  }
  K <- ncol(X)
  M <- ncol(Z)
  C <- matrix(
    0, nrow = M, ncol = K,
    dimnames = list(r$coef, colnames(X))
  )
  if (!K || !M) {
    return(C)
  }

  # An ordinal location intercept is a threshold coordinate rather than a
  # population location. Slopes retain the ordinary design-span map.
  eligible <- rep(TRUE, K)
  if (re_s2z_is_ordinal_location(bframe)) {
    eligible[colnames(X) == "Intercept"] <- FALSE
  }
  if (!any(eligible)) {
    return(C)
  }
  fit <- lm.fit(
    x = Z, y = X[, eligible, drop = FALSE],
    tol = sqrt(.Machine$double.eps), singular.ok = TRUE
  )
  candidate <- fit$coefficients
  if (is.null(dim(candidate))) {
    candidate <- matrix(candidate, ncol = 1L)
  }
  candidate[is.na(candidate)] <- 0
  # Stabilize the serialized coordinate map and generated Stan expressions.
  for (value in c(-1, 0, 1)) {
    close <- abs(candidate - value) <=
      128 * .Machine$double.eps * pmax(1, abs(candidate))
    candidate[close] <- value
  }
  reconstructed <- Z %*% candidate
  target <- X[, eligible, drop = FALSE]
  error <- apply(abs(reconstructed - target), 2L, max)
  scale <- pmax(
    1, apply(abs(target), 2L, max), apply(abs(reconstructed), 2L, max)
  )
  exact <- is.finite(error) & error <= 64 * .Machine$double.eps * scale
  if (any(exact)) {
    C[, which(eligible)[exact]] <- candidate[, exact, drop = FALSE]
  }
  C
}

# Return one fitted design-basis map per ordinary centered covariance block.
# These matrices are cached as part of the fitted basis so newdata cannot
# silently redefine the meaning of saved latent coordinates.
frame_re_center_mean <- function(bframe, data, cached = NULL) {
  stopifnot(is.bframel(bframe), is.data.frame(data))
  r_all <- bframe$frame$re
  if (!has_rows(r_all)) {
    return(list())
  }
  requested <- !r_all$s2z &
    (re_s2z_center_auto(r_all) | re_s2z_center_values(r_all) > 0)
  ids <- unique(r_all$id[requested])
  if (!length(ids)) {
    return(list())
  }
  if (!is.null(cached)) {
    if (!is.list(cached)) {
      stop2("Invalid stored population-location maps for group centering.")
    }
    out <- list()
    for (id in ids) {
      r <- subset2(r_all, id = id)
      C <- cached[[as.character(id)]]
      if (!is.matrix(C) || !is.numeric(C) || anyNA(C) ||
          any(!is.finite(C)) || nrow(C) != nrow(r) ||
          ncol(C) != length(bframe$frame$fe$vars)) {
        stop2("Stored population-location map does not match group-level ID ",
              id, ".")
      }
      old_r <- rownames(C)
      old_x <- colnames(C)
      if (!is.null(old_r) && !identical(old_r, r$coef) ||
          !is.null(old_x) && !identical(old_x, bframe$frame$fe$vars)) {
        stop2("Stored population-location map does not match group-level ID ",
              id, ".")
      }
      dimnames(C) <- list(r$coef, bframe$frame$fe$vars)
      out[[as.character(id)]] <- C
    }
    return(out)
  }
  out <- list()
  for (id in ids) {
    r <- subset2(r_all, id = id)
    out[[as.character(id)]] <- .re_center_mean_matrix(
      r, bframe = bframe, data = data
    )
  }
  out
}

# Classify a group-effect block without changing the legacy endpoint paths.
re_s2z_center_mode <- function(r) {
  rho <- re_s2z_center_values(r)
  auto <- re_s2z_center_auto(r)
  if (any(auto)) {
    if (!all(auto)) {
      stop2("All coefficients sharing a group-level ID must use ",
            "automatic centering if any coefficient does.")
    }
    return("auto")
  }
  if (any(rho > 0 & rho < 1) || length(unique(rho)) > 1L) {
    return("partial")
  }
  if (all(rho == 1)) "centered" else "noncentered"
}

# Construct the orthonormal basis used by sum_to_zero_constrain_brms(). Its
# jth column has 1 / sqrt(j * (j + 1)) in rows 1 through j and the negative
# balancing weight in row j + 1. Keeping the R and Stan bases identical lets
# data-only covariance eigensystems replace the levelwise S2Z chart without a
# second change of coordinates.
re_s2z_basis <- function(n) {
  n <- as_one_integer(n)
  if (n < 2L) {
    stop2("A sum-to-zero basis requires at least two levels.")
  }
  K <- n - 1L
  out <- matrix(0, nrow = n, ncol = K)
  for (j in seq_len(K)) {
    weight <- 1 / sqrt(j * (j + 1))
    out[seq_len(j), j] <- weight
    out[j + 1L, j] <- -j * weight
  }
  out
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

# One scalar argument of an S2Z population prior. `code` is the original Stan
# expression and `value` is populated only when R can prove its current value.
# Keeping these separate lets the completed-square code retain arbitrary
# parameter dependence without pretending an unknown expression is numeric.
new_re_s2z_prior_arg <- function(value = NA_real_, code = "", role,
                                 dependencies = character(),
                                 dependency_blocks = character(),
                                 scalarized = FALSE,
                                 broadcasted = FALSE) {
  stopifnot(
    is.character(code), length(code) == 1L,
    is.character(role), length(role) == 1L,
    is.character(dependencies), is.character(dependency_blocks),
    is.logical(scalarized),
    length(scalarized) == 1L, is.logical(broadcasted),
    length(broadcasted) == 1L
  )
  known <- (is.double(value) || is.integer(value)) &&
    length(value) == 1L && is.finite(value)
  if (!known) {
    value <- NA_real_
  }
  structure(
    nlist(
      code, value = as.numeric(value), known, role,
      dependencies = unique(dependencies),
      dependency_blocks = unique(dependency_blocks),
      scalarized, broadcasted
    ),
    class = "re_s2z_prior_arg"
  )
}

is.re_s2z_prior_arg <- function(x) {
  inherits(x, "re_s2z_prior_arg")
}

# A neutral placeholder used for fixed-only coordinates. Keeping the prior
# list full-length makes the descriptor backwards compatible with Stan code
# that still indexes it in population-coordinate order.
re_s2z_flat_prior <- function() {
  list(
    dist = "flat",
    location = new_re_s2z_prior_arg(0, role = "location"),
    scale = new_re_s2z_prior_arg(1, role = "scale"),
    df = new_re_s2z_prior_arg(role = "degrees-of-freedom"),
    scalarized = FALSE, broadcasted = FALSE
  )
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
re_s2z_prior <- function(prior, bframe, class, coef = "", r = NULL,
                         stanvars = NULL) {
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
  broadcast_index <- broadcast_length <- NULL
  if (class == "b" && nzchar(coef) && nrow(ans) && !nzchar(ans$coef)) {
    fixef <- bframe$frame$fe$vars_stan
    broadcast_index <- match(coef, fixef)
    if (!is.na(broadcast_index)) {
      broadcast_length <- length(fixef)
    } else {
      broadcast_index <- NULL
    }
  }
  out <- parse_re_s2z_prior(
    value, coef = coef, context = context, stanvars = stanvars,
    broadcast_index = broadcast_index,
    broadcast_length = broadcast_length
  )
  out$prior <- value
  out
}

# Parse the exact scalar priors supported by the S2Z paths. Logistic priors
# use an explicit omitted-mean fallback; the remaining proper distributions
# are conditionally Gaussian.
parse_re_s2z_prior <- function(prior, coef = "", context = NULL,
                               stanvars = NULL, broadcast_index = NULL,
                               broadcast_length = NULL) {
  # Internal whitespace can separate operators from signed literals (for
  # example `x < -1`); removing it before str2lang() can retokenize the Stan
  # expression as an R assignment. Only trim the outer prior string.
  prior <- trimws(prior)
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
      paste0(
        "use flat, normal, student_t, cauchy, or logistic with scalar ",
        "Stan expressions."
      )
    )
  }
  dist <- as.character(call[[1]])
  args <- as.list(call[-1])
  stanvars <- if (is.stanvars(stanvars)) stanvars else empty_stanvars()
  data_stanvars <- if (is.stanvars(stanvars)) {
    subset_stanvars(stanvars, block = "data")
  } else {
    empty_stanvars()
  }
  stan_expression_code <- function(x) {
    is_ternary <- function(expr) {
      is.call(expr) && identical(as.character(expr[[1L]]), "?") &&
        length(expr) == 3L && is.call(expr[[3L]]) &&
        identical(as.character(expr[[3L]][[1L]]), ":") &&
        length(expr[[3L]]) == 3L
    }
    render <- function(expr) {
      if (is_ternary(expr)) {
        return(paste0(
          render(expr[[2L]]), " ? ", render(expr[[3L]][[2L]]),
          " : ", render(expr[[3L]][[3L]])
        ))
      }
      replacements <- character()
      mask_ternary <- function(node) {
        if (is_ternary(node)) {
          placeholder <- paste0(
            "s2zTernaryPlaceholder", length(replacements) + 1L
          )
          replacements[[placeholder]] <<- paste0("(", render(node), ")")
          return(as.name(placeholder))
        }
        if (is.call(node)) {
          return(as.call(lapply(as.list(node), mask_ternary)))
        }
        node
      }
      code <- deparse0(mask_ternary(expr))
      for (placeholder in names(replacements)) {
        code <- gsub(
          placeholder, replacements[[placeholder]], code, fixed = TRUE
        )
      }
      code
    }
    render(x)
  }
  expression_arg <- function(x, argument) {
    code <- stan_expression_code(x)
    numeric_code <- grepl(
      paste0(
        "^[+-]?((([0-9]+(\\.[0-9]*)?)|(\\.[0-9]+))",
        "([eE][+-]?[0-9]+)?)$"
      ),
      code, perl = TRUE
    )
    if (numeric_code) {
      value <- suppressWarnings(as.numeric(code))
      if (!is.finite(value)) {
        fail(
          "active_prior_arguments",
          paste0(
            "The ", argument, " argument in a population-level prior ",
            "used with the sum-to-zero parameterization must be finite ",
            "(coefficient '", coef, "')."
          ),
          "replace it with a finite scalar Stan expression."
        )
      }
      return(new_re_s2z_prior_arg(value, role = argument))
    }

    dependencies <- all.vars(x)
    strip_stan_comments <- function(code) {
      code <- paste(code, collapse = "\n")
      code <- gsub("(?s)/\\*.*?\\*/", " ", code, perl = TRUE)
      gsub("//[^\r\n]*", " ", code, perl = TRUE)
    }
    top_level_stan_code <- function(code) {
      chars <- strsplit(strip_stan_comments(code), "", fixed = TRUE)[[1L]]
      if (!length(chars)) {
        return("")
      }
      depth <- 0L
      quoted <- FALSE
      escaped <- FALSE
      for (i in seq_along(chars)) {
        char <- chars[[i]]
        if (quoted) {
          if (depth > 0L && char != "\n") {
            chars[[i]] <- " "
          }
          if (escaped) {
            escaped <- FALSE
          } else if (char == "\\") {
            escaped <- TRUE
          } else if (char == '"') {
            quoted <- FALSE
          }
          next
        }
        if (char == '"') {
          quoted <- TRUE
          if (depth > 0L) {
            chars[[i]] <- " "
          }
        } else if (char == "{") {
          depth <- depth + 1L
          chars[[i]] <- ";"
        } else if (char == "}") {
          depth <- max(0L, depth - 1L)
          chars[[i]] <- ";"
        } else if (depth > 0L && char != "\n") {
          chars[[i]] <- " "
        }
      }
      paste0(chars, collapse = "")
    }
    dependency_entry <- function(symbol) {
      if (symbol %in% names(stanvars) &&
          stanvars[[symbol]]$block == "data") {
        return(stanvars[[symbol]])
      }
      type_pattern <- paste0(
        "(real|int|complex|vector|row_vector|matrix|complex_vector|",
        "complex_row_vector|complex_matrix|simplex|unit_vector|ordered|",
        "positive_ordered|corr_matrix|cov_matrix|cholesky_factor_corr|",
        "cholesky_factor_cov)"
      )
      declaration <- paste0(
        "\\b", type_pattern, "\\b[[:space:]]*(<[^;>]+>)?",
        "[[:space:]]*(\\[[^]]+\\])?[[:space:]]+", symbol, "\\b"
      )
      matches <- which(vapply(stanvars, function(info) {
        info$block != "functions" && grepl(
          declaration, top_level_stan_code(info$scode), perl = TRUE
        )
      }, logical(1)))
      if (length(matches) == 1L) stanvars[[matches]] else NULL
    }
    dependency_info <- lapply(dependencies, dependency_entry)
    names(dependency_info) <- dependencies
    missing_dependencies <- dependencies[vapply(
      dependency_info, is.null, logical(1)
    )]
    if (length(missing_dependencies)) {
      fail(
        "active_prior_arguments",
        paste0(
          "Symbolic ", argument, " argument '", code, "' in a ",
          "population-level prior used with the sum-to-zero ",
          "parameterization refers to unknown variable(s): ",
          collapse_comma(missing_dependencies), " (coefficient '", coef,
          "')."
        ),
        paste0(
          "define each value with stanvar() in the data, transformed-data, ",
          "parameters, or transformed-parameters-start block."
        )
      )
    }

    if (length(dependencies)) {
      valid_dependency <- vapply(dependency_info, function(info) {
        info$block %in% c("data", "tdata", "parameters") ||
          info$block == "tparameters" && info$position == "start"
      }, logical(1))
      if (!all(valid_dependency)) {
        invalid <- dependencies[!valid_dependency]
        fail(
          "active_prior_arguments",
          paste0(
            "Symbolic ", argument, " argument '", code, "' in a ",
            "population-level prior used with the sum-to-zero ",
            "parameterization depends on value(s) declared too late: ",
            collapse_comma(invalid), " (coefficient '", coef, "')."
          ),
          paste0(
            "declare these values in data, transformed data, parameters, ",
            "or transformed parameters with position = \"start\"."
          )
        )
      }
    }
    dependency_blocks <- unique(vapply(
      dependency_info, `[[`, character(1), "block"
    ))

    call_heads <- function(expr) {
      if (!is.call(expr)) {
        return(character())
      }
      head <- expr[[1L]]
      out <- if (is.name(head)) as.character(head) else character()
      nested <- unlist(lapply(as.list(expr[-1L]), call_heads), use.names = FALSE)
      unique(c(out, nested))
    }
    function_names <- call_heads(x)
    invalid_function <- grepl("(_rng|_lp)$", function_names) |
      function_names == "target"
    if (any(invalid_function)) {
      fail(
        "active_prior_arguments",
        paste0(
          "Symbolic ", argument, " argument '", code, "' in a ",
          "population-level prior used with the sum-to-zero ",
          "parameterization contains an RNG or target-modifying function ",
          "(coefficient '", coef, "')."
        ),
        "use a deterministic scalar Stan expression."
      )
    }

    symbol <- if (is.name(x)) as.character(x) else ""
    if (nzchar(symbol) && symbol %in% names(data_stanvars)) {
      stanvar_info <- data_stanvars[[symbol]]
      value <- stanvar_info$sdata
      data_scode <- strip_stan_comments(stanvar_info$scode)
      numeric_value <- is.double(value) || is.integer(value) ||
        is.logical(value)
      scalar_declaration <- paste0(
        "^[[:space:]]*(real|int)[[:space:]]*",
        "(<[^>]+>)?[[:space:]]+", symbol,
        "[[:space:]]*;[[:space:]]*$"
      )
      vector_declaration <- paste0(
        "^[[:space:]]*(vector|row_vector)[[:space:]]*",
        "(<[^>]+>)?[[:space:]]*\\[[^]]+\\][[:space:]]+", symbol,
        "[[:space:]]*;[[:space:]]*$"
      )
      array_declaration <- paste0(
        "^[[:space:]]*array[[:space:]]*\\[[^]]+\\][[:space:]]+",
        "(real|int)[[:space:]]*(<[^>]+>)?[[:space:]]+", symbol,
        "[[:space:]]*;[[:space:]]*$"
      )
      scalar_scode <- length(stanvar_info$scode) == 1L &&
        grepl(scalar_declaration, data_scode, perl = TRUE)
      vector_scode <- length(stanvar_info$scode) == 1L &&
        (grepl(vector_declaration, data_scode, perl = TRUE) ||
         grepl(array_declaration, data_scode, perl = TRUE))
      selected_index <- NULL
      if (scalar_scode && numeric_value && length(value) == 1L) {
        selected_value <- as.numeric(value)
        selected_code <- symbol
      } else if (vector_scode && numeric_value) {
        if (!is.null(broadcast_index)) {
          if (length(value) != broadcast_length) {
            fail(
              "active_prior_arguments",
              paste0(
                "Vector-valued ", argument, " argument '", symbol,
                "' has length ", length(value), " but the population-level ",
                "coefficient vector has length ", broadcast_length,
                " (coefficient '", coef, "')."
              ),
              "supply a scalar or a vector with one value per coefficient."
            )
          }
          selected_index <- broadcast_index
        } else if (length(value) == 1L) {
          selected_index <- 1L
        }
        if (!is.null(selected_index)) {
          selected_value <- as.numeric(value[selected_index])
          selected_code <- paste0(symbol, "[", selected_index, "]")
        }
      }
      valid_value <- exists("selected_value", inherits = FALSE) &&
        length(selected_value) == 1L && is.finite(selected_value)
      if (!valid_value) {
        fail(
          "active_prior_arguments",
          paste0(
            "Symbolic ", argument, " argument '", symbol, "' in a ",
            "population-level prior used with the sum-to-zero ",
            "parameterization must resolve to one finite scalar value ",
            "(coefficient '", coef, "')."
          ),
          paste0(
            "use a scalar, index one element explicitly, or supply one ",
            "value per population-level coefficient."
          )
        )
      }
      return(new_re_s2z_prior_arg(
        selected_value, code = selected_code, role = argument,
        dependencies = dependencies, dependency_blocks = dependency_blocks,
        scalarized = !identical(selected_code, symbol),
        broadcasted = !is.null(broadcast_index)
      ))
    }
    if (nzchar(symbol) && length(dependencies) == 1L) {
      stanvar_info <- dependency_info[[symbol]]
      scode <- top_level_stan_code(stanvar_info$scode)
      declaration_boundary <- "(^|[;{}])[[:space:]]*"
      scalar_declaration <- paste0(
        declaration_boundary,
        "(real|int)[[:space:]]*(<[^;>]+>)?[[:space:]]+",
        symbol, "\\b"
      )
      vector_declaration <- paste0(
        declaration_boundary,
        "(vector|row_vector)[[:space:]]*(<[^;>]+>)?",
        "[[:space:]]*\\[([^]]+)\\][[:space:]]+", symbol,
        "\\b"
      )
      array_declaration <- paste0(
        declaration_boundary,
        "array[[:space:]]*\\[([^]]+)\\][[:space:]]+",
        "(real|int)[[:space:]]*(<[^;>]+>)?[[:space:]]+", symbol,
        "\\b"
      )
      scalar_scode <- grepl(scalar_declaration, scode, perl = TRUE)
      vector_scode <- grepl(vector_declaration, scode, perl = TRUE) ||
        grepl(array_declaration, scode, perl = TRUE)
      if (scalar_scode) {
        return(new_re_s2z_prior_arg(
          code = symbol, role = argument, dependencies = dependencies,
          dependency_blocks = dependency_blocks
        ))
      }
      if (vector_scode) {
        literal_1d_length <- function(pattern) {
          match <- regexec(pattern, scode, perl = TRUE)
          hit <- regmatches(scode, match)[[1L]]
          if (length(hit) < 3L) integer() else as.integer(trimws(hit[3L]))
        }
        vector_length_pattern <- paste0(
          declaration_boundary,
          "(?:vector|row_vector)[[:space:]]*(?:<[^;>]+>)?",
          "[[:space:]]*\\[([[:space:]]*[0-9]+[[:space:]]*)\\]",
          "[[:space:]]+", symbol, "\\b"
        )
        array_length_pattern <- paste0(
          declaration_boundary,
          "array[[:space:]]*\\[([[:space:]]*[0-9]+[[:space:]]*)\\]",
          "[[:space:]]+(?:real|int)\\b[[:space:]]*",
          "(?:<[^;>]+>)?[[:space:]]+", symbol, "\\b"
        )
        literal_length <- unique(c(
          literal_1d_length(vector_length_pattern),
          literal_1d_length(array_length_pattern)
        ))
        selected_index <- broadcast_index
        if (is.null(selected_index)) {
          if (length(literal_length) == 1L && literal_length == 1L) {
            selected_index <- 1L
          }
        }
        if (is.null(selected_index)) {
          fail(
            "active_prior_arguments",
            paste0(
              "Vector-valued ", argument, " argument '", symbol,
              "' must identify one scalar coordinate in a population-level ",
              "prior used with the sum-to-zero parameterization ",
              "(coefficient '", coef, "')."
            ),
            "index one element explicitly or use a scalar declaration."
          )
        }
        if (!is.null(broadcast_index) && length(literal_length) &&
            (length(literal_length) != 1L ||
             literal_length != broadcast_length)) {
          fail(
            "active_prior_arguments",
            paste0(
              "Vector-valued ", argument, " argument '", symbol,
              "' has declared length ", collapse_comma(literal_length),
              " but the population-level coefficient vector has length ",
              broadcast_length, " (coefficient '", coef, "')."
            ),
            "supply a scalar or a vector with one value per coefficient."
          )
        }
        selected_code <- paste0(symbol, "[", selected_index, "]")
        if (!is.null(broadcast_index)) {
          # Retain ordinary vectorized-prior semantics by checking the complete
          # runtime length before selecting a coordinate. Literal declarations
          # have already received the same earlier, contextual size check.
          selected_code <- paste0(
            "s2z_prior_coordinate_brms(", symbol, ", ", selected_index,
            ", ", broadcast_length, ")"
          )
        }
        return(new_re_s2z_prior_arg(
          code = selected_code,
          role = argument, dependencies = dependencies,
          dependency_blocks = dependency_blocks,
          scalarized = TRUE,
          broadcasted = !is.null(broadcast_index)
        ))
      }
      fail(
        "active_prior_arguments",
        paste0(
          "Symbolic ", argument, " argument '", symbol, "' must be a ",
          "scalar or one-dimensional Stan value (coefficient '", coef,
          "')."
        ),
        "use a real, int, vector, row_vector, or one-dimensional array."
      )
    }
    if (!is.null(broadcast_index)) {
      # The same global prior expression may evaluate to either a scalar
      # (ordinary Stan broadcasting) or a one-dimensional value. Let Stan's
      # overload resolution preserve that distinction, and require a vector
      # result to have exactly one entry per original population coefficient.
      code <- paste0(
        "s2z_prior_coordinate_brms(", code, ", ", broadcast_index,
        ", ", broadcast_length, ")"
      )
      return(new_re_s2z_prior_arg(
        code = code, role = argument, dependencies = dependencies,
        dependency_blocks = dependency_blocks,
        scalarized = TRUE, broadcasted = TRUE
      ))
    }
    new_re_s2z_prior_arg(
      code = code, role = argument, dependencies = dependencies,
      dependency_blocks = dependency_blocks
    )
  }
  if (dist == "std_normal" && !length(args)) {
    out <- list(
      dist = "normal",
      location = new_re_s2z_prior_arg(0, role = "location"),
      scale = new_re_s2z_prior_arg(1, role = "scale"),
      df = new_re_s2z_prior_arg(role = "degrees-of-freedom")
    )
  } else if (dist == "normal" && length(args) == 2L) {
    out <- list(
      dist = "normal",
      location = expression_arg(args[[1]], "location"),
      scale = expression_arg(args[[2]], "scale"),
      df = new_re_s2z_prior_arg(role = "degrees-of-freedom")
    )
  } else if (dist == "student_t" && length(args) == 3L) {
    out <- list(
      dist = "student",
      df = expression_arg(args[[1]], "degrees-of-freedom"),
      location = expression_arg(args[[2]], "location"),
      scale = expression_arg(args[[3]], "scale")
    )
  } else if (dist == "cauchy" && length(args) == 2L) {
    out <- list(
      dist = "student",
      df = new_re_s2z_prior_arg(1, role = "degrees-of-freedom"),
      location = expression_arg(args[[1]], "location"),
      scale = expression_arg(args[[2]], "scale")
    )
  } else if (dist == "logistic" && length(args) == 2L) {
    out <- list(
      dist = "logistic",
      location = expression_arg(args[[1]], "location"),
      scale = expression_arg(args[[2]], "scale"),
      df = new_re_s2z_prior_arg(role = "degrees-of-freedom")
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
  invalid_scale <- out$scale$known && !(out$scale$value > 0)
  invalid_df <- out$dist == "student" && out$df$known &&
    !(out$df$value > 0)
  if (invalid_scale || invalid_df) {
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
  out$scalarized <- any(vapply(
    out[c("location", "scale", "df")],
    function(x) is.re_s2z_prior_arg(x) && x$scalarized,
    logical(1)
  ))
  out$broadcasted <- any(vapply(
    out[c("location", "scale", "df")],
    function(x) is.re_s2z_prior_arg(x) && x$broadcasted,
    logical(1)
  ))
  out
}

# Check that an estimated scale parameter retains its required support.
validate_re_s2z_scale_bounds <- function(
    prior, class, setting = "gr(..., s2z = TRUE)"
) {
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
          "bound for ", setting, ".")
  }
  if (nzchar(ub)) {
    ub_num <- numeric_bound(ub)
    if (length(ub_num) != 1L || !is.finite(ub_num) || ub_num <= 0 ||
        ub_num <= lb_num) {
      stop2("Class '", class, "' must have a finite positive upper bound ",
            "greater than its lower bound when one is specified for ",
            setting, ".")
    }
  }
  invisible(NULL)
}

# Check fixed group scales before the Gaussian precision is constructed.
validate_re_s2z_sd_prior <- function(prior, r, bframe = NULL) {
  stopifnot(is.brmsprior(prior), is.reframe(r), has_rows(r))
  setting <- if (all(r$s2z)) {
    "gr(..., s2z = TRUE)"
  } else {
    "positive ordinary group-effect centering"
  }
  px <- check_prefix(r)
  p <- subset2(
    prior, class = "sd", coef = c(r$coef, ""),
    group = c(r$group[1], ""), ls = px
  )
  validate_re_s2z_scale_bounds(p, class = "sd", setting = setting)
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
        "positive numeric scalars for ", setting, " ",
        "(coefficient '", coef, "')."
      )
      if (is.null(bframe)) {
        stop2(problem)
      }
      if (!all(r$s2z)) {
        stop2(
          problem, " Use constant(value) with a finite value greater than ",
          "zero, or estimate the scale."
        )
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

# Resolve the scalar SD prior attached to one strict latent-score occurrence.
# Shared wide-response scores have only one covariance scale, so all aliases
# must imply the same prior even though brms's multivariate prior table is
# response-qualified.
re_s2z_effective_sd_prior <- function(prior, r) {
  stopifnot(
    is.brmsprior(prior), is.reframe(r), nrow(r) == 1L,
    isTRUE(r$s2z), isTRUE(re_s2z_latent(r))
  )
  px <- check_prefix(r)
  selected <- subset2(
    prior, class = "sd", coef = c(r$coef, ""),
    group = c(r$group, ""), ls = px
  )
  coefficient <- subset2(selected, coef = r$coef)
  coefficient <- coefficient[nzchar(coefficient$prior) |
    nzchar(coefficient$tag), , drop = FALSE]
  stopifnot(nrow(coefficient) <= 1L)
  value <- if (nrow(coefficient) && nzchar(coefficient$prior)) {
    coefficient$prior
  } else {
    stan_base_prior(selected)
  }
  tag <- if (nrow(coefficient) && nzchar(coefficient$tag)) {
    coefficient$tag
  } else {
    stan_base_prior(selected, col = "tag")
  }
  bounds <- stan_base_prior(selected, col = c("lb", "ub"))
  list(
    prior = as.character(value), tag = as.character(tag),
    lb = as.character(bounds$lb), ub = as.character(bounds$ub)
  )
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

# Validate the unrestricted centering chart for ordinary group effects. Plan 3
# deliberately covers predictor-local ordinary gr() covariance blocks only;
# structural combinations owned by later plans keep their established
# non-centered implementation.
validate_re_centered_group_effects <- function(bframe, prior) {
  stopifnot(is.anybrmsframe(bframe), is.brmsprior(prior))
  r_global <- bframe$frame$re
  stopifnot(is.reframe(r_global))
  if (!has_rows(r_global)) {
    return(invisible(NULL))
  }
  requested <- !r_global$s2z &
    (re_s2z_center_auto(r_global) | re_s2z_center_values(r_global) > 0)
  ids <- unique(r_global$id[requested])
  if (!length(ids)) {
    return(invisible(NULL))
  }
  frames <- all_bframel(bframe)
  for (id in ids) {
    r <- subset2(r_global, id = id)
    mode <- re_s2z_center_mode(r)
    if (any(r$s2z)) {
      stop2("All coefficients sharing a group-level ID must use the same ",
            "'s2z' setting.")
    }
    auto <- re_s2z_center_auto(r)
    if (any(auto) && !all(auto)) {
      stop2("All coefficients sharing a group-level ID must use automatic ",
            "centering if any coefficient does.")
    }
    if (length(unique(r$group)) != 1L) {
      stop2("Centered ordinary group effects sharing an ID must use the ",
            "same grouping factor.")
    }
    if (length(unique(r$dist)) != 1L ||
        !r$dist[1L] %in% c("gaussian", "student")) {
      distributions <- paste(unique(r$dist), collapse = ", ")
      stop2("Center mode '", mode, "' is not available for ordinary ",
            "group-effect distribution '", distributions, "'. Supported ",
            "distributions are 'gaussian' and 'student'; use center = 0 ",
            "for another distribution.")
    }
    if (length(unique(r$cor)) != 1L) {
      stop2("Centered ordinary group effects sharing an ID must use the ",
            "same 'cor' setting.")
    }
    if (length(unique(r$scale)) != 1L || r$scale[1L] != "shared") {
      stop2("Centered ordinary group effects currently require ",
            "scale = \"shared\".")
    }
    if (any(nzchar(r$gtype), na.rm = TRUE) ||
        any(nzchar(r$type), na.rm = TRUE)) {
      stop2("Center mode '", mode, "' for ordinary group effects currently ",
            "supports ordinary gr() coefficients only; multi-membership and ",
            "special group coefficients retain center = 0.")
    }
    has_pw <- any(vapply(r$gcall, function(gcall) {
      isTRUE(nzchar(gcall$pw))
    }, logical(1)))
    if (any(nzchar(r$by), na.rm = TRUE) ||
        any(nzchar(r$cov), na.rm = TRUE) || has_pw) {
      stop2("Center mode '", mode, "' for ordinary group effects is not yet ",
            "supported with 'by', 'cov', or 'pw'; use center = 0.")
    }
    id_frames <- Filter(function(x) {
      rx <- x$frame$re
      has_rows(rx) && id %in% rx$id
    }, frames)
    if (length(id_frames) != 1L) {
      stop2("Centered ordinary group-level IDs must be predictor-local; ",
            "use distinct IDs in separate response, distributional, or ",
            "nonlinear predictors, or use center = 0.")
    }
    if (has_re_s2z_conventional(id_frames[[1L]])) {
      stop2("Slide-style centering for ordinary group effects cannot yet be ",
            "combined with a conventional sum-to-zero group effect in the ",
            "same linear predictor; use center = 0 for the ordinary block ",
            "or place the blocks in separate predictors.")
    }
    if (identical(mode, "auto")) {
      if (re_s2z_is_ordinal_location(id_frames[[1L]])) {
        stop2("The automatic-centering precursor is not yet supported for ",
              "ordinary group effects in ordinal location predictors; use ",
              "a fixed center value in [0, 1].")
      }
      if (any(nzchar(r$nlpar))) {
        stop2("The automatic-centering precursor is not yet supported for ",
              "ordinary group effects in nonlinear predictors; use a fixed ",
              "center value in [0, 1].")
      }
    }
    validate_re_s2z_sd_prior(prior, r, bframe = id_frames[[1L]])
  }
  invisible(NULL)
}

# Preserve the group-effects branch's global validation of latent, centered,
# varying-scale, covariance, and cross-predictor S2Z blocks. PLAN-02's
# active-coordinate prior checks are layered on by validate_re_s2z_prior_global.
.validate_re_s2z_group_effects <- function(bframe, prior, stanvars = NULL) {
  stopifnot(is.anybrmsframe(bframe), is.brmsprior(prior))
  validate_re_centered_group_effects(bframe, prior = prior)
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
  # Conventional S2Z blocks usually have one omitted-mean system per linear
  # predictor. A supported multivariate-response block assembles those local
  # population coordinates into one global system. Strict
  # latent-score blocks omit that mean entirely, so one covariance block may
  # supply columns to several nonlinear predictors. Validate the distinction
  # globally before any local omitted-mean descriptors are constructed.
  r_global <- bframe$frame$re
  stopifnot(is.reframe(r_global))
  s2z_ids <- unique(r_global$id[r_global$s2z])
  for (id in s2z_ids) {
    r_id <- subset2(r_global, id = id)
    latent <- re_s2z_latent(r_id)
    if (any(latent) && !all(latent)) {
      stop2("All coefficients sharing a group-level ID must use the same ",
            "'latent' setting.")
    }
    n_frames <- sum(vapply(all_frames, function(x) {
      r <- x$frame$re
      has_rows(r) && id %in% r$id
    }, logical(1)))
    if (!all(latent) && n_frames > 1L) {
      id_frames <- Filter(function(x) {
        r <- x$frame$re
        has_rows(r) && id %in% r$id
      }, all_frames)
      if (!all(r_id$s2z)) {
        contexts <- lapply(id_frames, function(x) {
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
            "predictors unless every occurrence uses s2z = TRUE.\n",
            affected
          ),
          paste0(
            "use a distinct ID in every listed linear predictor (for ",
            "example, ", re_s2z_local_id_examples(contexts), "). For ",
            "`mvbind(...)` or category shorthand, omit the shared `| ",
            public_id, " |` tag or `id = \"", public_id, "\"` so brms ",
            "allocates predictor-local IDs; alternatively, expand the model ",
            "into separate `bf()` response formulas with those distinct ",
            "IDs. These rewrites do not retain cross-predictor group-effect ",
            "correlations; use s2z = FALSE if those correlations are required."
          )
        )
      }
      ordinal_frames <- Filter(re_s2z_is_ordinal_location, id_frames)
      if (length(ordinal_frames)) {
        ordinal_r <- subset2(ordinal_frames[[1L]]$frame$re, id = id)
        stop_re_s2z(
          re_s2z_context(ordinal_frames[[1L]], r = ordinal_r),
          "ordinal_cross_predictor_id",
          paste0(
            "An ordinal sum-to-zero group-level ID cannot yet span ",
            "multiple linear predictors."
          ),
          "use predictor-local S2Z IDs for ordinal location predictors."
        )
      }
      validate_re_s2z_cross_id(
        bframe, prior = prior, id = id, stanvars = stanvars
      )
    }
    if (!all(r_id$s2z)) {
      stop2("All coefficients sharing a strict or conventional sum-to-zero ",
            "group-level ID must use 's2z = TRUE'.")
    }
    if (!all(latent)) {
      next
    }
    if (n_frames > 1L && any(!nzchar(r_id$nlpar))) {
      stop2("A strict latent-score ID may span only nonlinear predictors.")
    }
    if (length(unique(r_id$group)) != 1L) {
      stop2("All coefficients sharing a strict latent-score ID must use ",
            "the same grouping factor.")
    }
    if (length(unique(r_id$cov)) != 1L) {
      stop2("All coefficients sharing a sum-to-zero group-level ID must use ",
            "the same 'cov' setting.")
    }
    if (length(unique(r_id$scale)) != 1L ||
        !identical(r_id$scale[1], "shared")) {
      stop2("Strict latent-score S2Z blocks currently require ",
            "scale = \"shared\".")
    }
    if (length(unique(r_id$dist)) != 1L ||
        !identical(r_id$dist[1], "gaussian")) {
      stop2("Strict latent-score S2Z blocks currently require Gaussian ",
            "group effects.")
    }
    if (length(unique(r_id$cor)) != 1L) {
      stop2("All coefficients sharing a strict latent-score ID must use ",
            "the same 'cor' setting.")
    }
    mode <- re_s2z_center_mode(r_id)
    if (!mode %in% c("centered", "noncentered", "partial", "auto")) {
      stop2("Strict latent-score S2Z blocks currently support only ",
            "centered, noncentered, fixed partial, and automatic centering ",
            "modes.")
    }
    if (any(nzchar(r_id$gtype)) || any(nzchar(r_id$type))) {
      stop2("Strict latent-score S2Z blocks currently support only ",
            "ordinary 'gr' terms.")
    }
    has_pw <- any(vapply(r_id$gcall, function(gcall) {
      isTRUE(nzchar(gcall$pw))
    }, logical(1)))
    if (any(nzchar(r_id$by), na.rm = TRUE) ||
        any(nzchar(r_id$cov), na.rm = TRUE) || has_pw) {
      stop2("Arguments 'by', 'cov', and 'pw' are not yet supported with ",
            "strict latent-score S2Z blocks.")
    }
    if (length(get_levels(r_id)[[r_id$group[1]]]) < 2L) {
      stop2("A strict latent-score S2Z block requires at least two ",
            "observed grouping levels.")
    }
    if (has_special_prior(prior, check_prefix(r_id), class = "sd")) {
      stop2("Special priors on group-level scales are not yet supported ",
            "with strict latent-score S2Z blocks.")
    }
    # Scale priors retain the nonlinear-parameter prefix of each occurrence.
    # The legacy validator expects one local prefix at a time, whereas a
    # strict score covariance block is intentionally allowed to span them.
    prefix_key <- combine_prefix(check_prefix(r_id))
    for (key in unique(prefix_key)) {
      validate_re_s2z_sd_prior(
        prior, r_id[prefix_key == key, , drop = FALSE]
      )
    }
    latent_key <- re_s2z_latent_key(r_id)
    for (key in unique(latent_key)) {
      rows <- which(latent_key == key)
      if (length(rows) < 2L) {
        next
      }
      specs <- lapply(rows, function(row) {
        re_s2z_effective_sd_prior(
          prior, r_id[row, , drop = FALSE]
        )
      })
      signatures <- vapply(
        specs, function(x) paste(unlist(x), collapse = "\r"), character(1)
      )
      if (length(unique(signatures)) != 1L) {
        stop2(
          "A strict latent score shared across responses has one covariance ",
          "scale and therefore requires identical 'sd' priors for nonlinear ",
          "parameter '", r_id$nlpar[rows[1L]], "' and coefficient '",
          r_id$coef[rows[1L]], "' in every response."
        )
      }
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
    infos <- re_s2z_infos(x)
    for (info in infos) {
      r <- info$r
      cross_id <- is_re_s2z_cross_id(bframe, info$id)
      if (info$id %in% ids && !cross_id) {
        stop2("A sum-to-zero group-level ID cannot span multiple ",
              "linear predictors.")
      }
      if (!info$id %in% ids) {
        ids <- c(ids, info$id)
      }
      if (length(unique(r$group)) != 1L) {
        stop2("All coefficients sharing a sum-to-zero group-level ID must ",
              "use the same grouping factor.")
      }
      if (length(unique(r$cov)) != 1L) {
        stop2("All coefficients sharing a sum-to-zero group-level ID must ",
              "use the same 'cov' setting.")
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
      auto <- re_s2z_center_auto(r)
      if (any(auto) && !all(auto)) {
        stop2("All coefficients sharing a sum-to-zero group-level ID must use ",
              "automatic centering if any coefficient does.")
      }
      if (any(nzchar(r$gtype)) || any(nzchar(r$type))) {
        stop2("The sum-to-zero parameterization currently ",
              "supports only ordinary 'gr' terms.")
      }
      has_pw <- any(vapply(r$gcall, function(gcall) {
        isTRUE(nzchar(gcall$pw))
      }, logical(1)))
      if (any(nzchar(r$by), na.rm = TRUE) || has_pw) {
        stop2("Arguments 'by' and 'pw' are not yet supported ",
              "together with gr(..., s2z = TRUE).")
      }
      if (length(get_levels(r)[[r$group[1]]]) < 2L) {
        stop2("The sum-to-zero parameterization requires at ",
              "least two observed grouping levels.")
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
          has_special_prior(prior, check_prefix(r), class = "sdlog")) {
        stop2("Special priors are not yet supported together with ",
              "gr(..., s2z = TRUE).")
      }
      validate_re_s2z_sd_prior(prior, r, bframe = x)
      validate_re_s2z_sdlog_prior(prior, r)
    }
  }
  invisible(NULL)
}

# Build ordinal block descriptors around the exact affine identity
# Z = 1 a' + Xc C. q always consists of threshold primitives first and
# population slopes second. H follows q = theta - H m, hence its threshold
# rows are -a and its population-slope rows are +C.
.re_s2z_ordinal_descriptors <- function(bframe, data = NULL) {
  stopifnot(is.bframel(bframe), re_s2z_is_ordinal_location(bframe))
  re <- bframe$frame$re
  ids <- unique(re$id[re$s2z & !re_s2z_latent(re)])
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
    re_s2z_center_values(r)
    re_s2z_center_auto(r)
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
  ids <- unique(re$id[re$s2z & !re_s2z_latent(re)])
  qnames <- bframe$frame$fe$vars
  center <- bframe$frame$fe$center
  out <- lapply(ids, function(id) {
    r <- subset2(re, id = id)
    re_s2z_center_values(r)
    re_s2z_center_auto(r)
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
re_s2z_threshold_prior <- function(prior, bframe, threshold, r,
                                   stanvars = NULL,
                                   broadcast_index = NULL,
                                   broadcast_length = NULL) {
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
  if (!nrow(ans) || nzchar(ans$coef)) {
    broadcast_index <- broadcast_length <- NULL
  }
  spec <- parse_re_s2z_prior(
    value, coef = label, context = context, stanvars = stanvars,
    broadcast_index = broadcast_index,
    broadcast_length = broadcast_length
  )
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
.re_s2z_attach_priors <- function(infos, bframe, prior, stanvars = NULL) {
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
      same_parameter <- threshold$group == threshold$group[i] &
        threshold$source == threshold$source[i]
      specs[[qi]] <- re_s2z_threshold_prior(
        prior, bframe, threshold = threshold[i, , drop = FALSE],
        r = infos[[1L]]$r, stanvars = stanvars,
        broadcast_index = threshold$threshold_index[i],
        broadcast_length = max(threshold$threshold_index[same_parameter])
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
        prior, bframe, class = "Intercept", r = r_context,
        stanvars = stanvars
      )
    } else {
      specs[[i]] <- re_s2z_prior(
        prior, bframe, class = "b", coef = qnames[i], r = r_context,
        stanvars = stanvars
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
re_s2z_infos <- function(bframe, prior = NULL, data = NULL, stanvars = NULL) {
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
    out <- .re_s2z_attach_priors(
      out, bframe = bframe, prior = prior, stanvars = stanvars
    )
  }
  out
}

# Describe one local S2Z block, including the fixed/group design mapping.
re_s2z_info <- function(bframe, prior = NULL, id = NULL, stanvars = NULL) {
  stopifnot(is.bframel(bframe))
  infos <- re_s2z_infos(bframe, prior = prior, stanvars = stanvars)
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

# The exact explicit-mean fallback for non-conditionally-Gaussian population
# priors composes with every supported centering chart for predictor-local
# blocks with shared scales and no known grouping covariance. Later structural
# extensions must not silently reuse its density or recovery algebra.
validate_re_s2z_logistic_chart <- function(infos, bframe, global_bframe) {
  stopifnot(
    is.list(infos), is.bframel(bframe), is.anybrmsframe(global_bframe)
  )
  if (!length(infos) || isTRUE(infos[[1L]]$ordinal)) {
    return(invisible(NULL))
  }
  specs <- infos[[1L]]$prior
  logistic_q <- which(vapply(specs, function(spec) {
    identical(spec$dist, "logistic")
  }, logical(1)))
  if (!length(logistic_q)) {
    return(invisible(NULL))
  }
  incompatible <- vapply(infos, function(info) {
    any(info$r$scale != "shared") || any(nzchar(info$r$cov)) ||
      is_re_s2z_cross_id(global_bframe, info$id)
  }, logical(1))
  if (!any(incompatible)) {
    return(invisible(NULL))
  }
  first <- infos[[which(incompatible)[1L]]]
  reasons <- unique(unlist(lapply(infos[incompatible], function(info) {
    c(
      if (any(info$r$scale != "shared")) "group-varying scales",
      if (any(nzchar(info$r$cov))) "a known grouping covariance",
      if (is_re_s2z_cross_id(global_bframe, info$id)) {
        "a cross-predictor ID"
      }
    )
  }), use.names = FALSE))
  prior_label <- unique(vapply(specs[logistic_q], function(spec) {
    spec$prior %||% "logistic"
  }, character(1)))
  stop_re_s2z(
    re_s2z_context(
      bframe, r = first$r, coef = infos[[1L]]$qnames[logistic_q],
      prior = paste(prior_label, collapse = ", ")
    ),
    "logistic_explicit_mean_chart",
    paste0(
      "Logistic priors on S2Z-active coordinates currently require a ",
      "predictor-local, shared-scale group block with no known grouping ",
      "covariance; found ",
      paste(reasons, collapse = ", "), "."
    ),
    paste0(
      "use scale = \"shared\", no cov argument, and a predictor-local ID; ",
      "otherwise use a normal, student_t, or cauchy prior."
    )
  )
}

# Construct an S2Z design matrix in covariance-block coefficient order. A
# block may be assembled from multiple grouping terms that share an ID.
re_s2z_design_matrix <- function(bframe, data, id = NULL) {
  stopifnot(is.bframel(bframe))
  re <- bframe$frame$re
  ids <- unique(re$id[re$s2z & !re_s2z_latent(re)])
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
  ids <- unique(re$id[re$s2z & !re_s2z_latent(re)])
  if (!length(ids)) {
    return(invisible(NULL))
  }
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
    if (ordinal_location && identical(re_s2z_center_mode(r), "auto")) {
      stop_re_s2z(
        context, "ordinal_fisher_centering",
        paste0(
          "The automatic-centering precursor is not yet supported for ordinal ",
          "sum-to-zero location predictors."
        ),
        "use a fixed center value in [0, 1], or set s2z = FALSE."
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
    if (has_by || has_pw) {
      capability <- if (has_by) "by" else "pw"
      remedy <- switch(
        capability,
        by = paste0(
          "remove 'by' or use s2z = FALSE to retain separate ",
          "variance-covariance matrices for its levels."
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
  if (!length(infos)) {
    return(invisible(list()))
  }
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
validate_re_s2z_prior_global <- function(bframe, prior, stanvars = NULL) {
  stopifnot(is.anybrmsframe(bframe), is.brmsprior(prior))
  .validate_re_s2z_group_effects(
    bframe, prior = prior, stanvars = stanvars
  )
  all_frames <- all_bframel(bframe)
  frames <- Filter(has_re_s2z, all_frames)
  if (!length(frames)) {
    return(invisible(NULL))
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
    if (!length(infos)) {
      # Strict latent-score blocks have no omitted population mean.
      next
    }
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
    prior_infos <- re_s2z_infos(
      x, prior = prior, stanvars = stanvars
    )
    validate_re_s2z_logistic_chart(
      prior_infos, bframe = x, global_bframe = bframe
    )
    for (info in infos) {
      validate_re_s2z_sd_prior(prior, info$r, bframe = x)
    }
  }
  invisible(NULL)
}

# Compatibility entry point for downstream code while callers migrate to the
# explicitly named prior/global phase.
validate_re_s2z <- function(bframe, prior, stanvars = NULL) {
  validate_re_s2z_prior_global(
    bframe, prior = prior, stanvars = stanvars
  )
}
