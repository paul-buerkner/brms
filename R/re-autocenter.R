# Helpers for estimating group-effect centering fractions in a precursor fit.
# The final Stan program only consumes the fixed numeric matrices produced by
# these helpers; none of the precursor settings are part of its target.

#' Configure and inspect a group-effect centering precursor
#'
#' Configure the separate precursor used by
#' \code{gr(..., center = "auto")} or extract the fixed centering matrices
#' stored in a fitted model.
#'
#' @param method Precursor algorithm. The default \code{"pathfinder"} uses
#'   CmdStanR Pathfinder. \code{"hmc"} uses a separate short HMC run; it does
#'   not modify the parameterization during warmup.
#' @param aggregate How valid candidate fractions are combined across precursor
#'   draws. The default is \code{"median"}; \code{"mean"} is also available.
#' @param fallback Value used if a candidate cell has no valid precursor draws:
#'   \code{"error"}, \code{"noncentered"}, \code{"centered"}, or a number in
#'   \code{[0, 1]}.
#' @param pilot_args A named list of additional arguments passed to the
#'   CmdStanR Pathfinder or HMC precursor. These arguments control only the
#'   precursor, not the final fit. Pathfinder uses multiple paths by default:
#'   \code{num_paths = min(max(1, chains), 4)}, where \code{chains} is the
#'   number of chains requested for the final fit. Thus, the usual four-chain
#'   fit uses four Pathfinder paths, while a one-chain fit uses one path unless
#'   \code{pilot_args$num_paths} is supplied explicitly. Pathfinder defaults
#'   both \code{draws} and \code{single_path_draws} to 200 to limit
#'   generated-quantities work. Supplying \code{draws} sets both defaults
#'   unless \code{single_path_draws} is supplied explicitly.
#'
#' @details Automatic centering is a two-stage workflow. The precursor is
#'   always fit in fully non-centered coordinates. At each precursor draw,
#'   level- and coefficient-specific candidate fractions are evaluated only in
#'   generated quantities, so they do not enter the precursor target density.
#'   Valid candidates are aggregated by the requested method. If a cell has no
#'   valid candidate, \code{fallback} either raises an error or supplies the
#'   requested fixed value.
#'   When multicategory shorthand expands one source \code{gr()} term into
#'   several generated mean predictors, their resolved proposals are combined
#'   elementwise by the same \code{aggregate} rule and the resulting fixed
#'   matrix is shared by those predictors.
#'
#'   The aggregated matrices are then passed as fixed data to a separate final
#'   HMC fit with a fresh warmup. No centering weight is sampled or updated
#'   within that fit. Automatic centering currently requires
#'   \code{backend = "cmdstanr"} and \code{algorithm = "sampling"} for the
#'   final fit. Use \code{centering_weights(fit)} to retrieve the matrices
#'   retained in the \code{brmsfit} object.
#'
#'   Once resolved, the fitted formula contains the fixed matrices. Ordinary
#'   refits therefore reuse them; refits on a subset of the fitted grouping
#'   levels align the stored rows by name. New levels or new coefficients need
#'   new weights. To run the stored request again, supply
#'   \code{center_control} to \code{update()}; alternatively, update with a
#'   formula that sets \code{center = "auto"} on every occurrence sharing the
#'   affected group ID.
#'   \code{brm_multiple()} compiles this common Stan program once and runs a
#'   separate precursor for each fitted data set.
#'
#' @return \code{autocenter_control()} returns a validated control list.
#'   \code{centering_weights()} returns named level-by-coefficient matrices.
#'   For a combined \code{brmsfit_multiple}, it returns one such named list per
#'   component data set.
#'
#' @examples
#' # Use eight Pathfinder paths for the automatic-centering precursor.
#' pathfinder_control <- autocenter_control(
#'   pilot_args = list(num_paths = 8)
#' )
#' @export
autocenter_control <- function(method = "pathfinder", aggregate = "median",
                               fallback = "error", pilot_args = list()) {
  validate_re_autocenter_control(nlist(
    method, aggregate, fallback, pilot_args
  ))
}

# Validate the small set of controls shared by Pathfinder and HMC precursors.
# A numeric fallback is a fixed centering fraction. The character endpoint
# aliases are retained in the validated object so diagnostics remain readable.
validate_re_autocenter_control <- function(control = NULL) {
  if (is.null(control)) {
    control <- list()
  }
  if (!is.list(control)) {
    stop2("Group-effect autocenter 'control' must be a list.")
  }
  control_names <- names(control)
  if (length(control) &&
      (is.null(control_names) || any(!nzchar(control_names)))) {
    stop2("Group-effect autocenter 'control' must be a named list.")
  }
  if (anyDuplicated(control_names)) {
    stop2("Group-effect autocenter 'control' names must be unique.")
  }
  valid_names <- c("method", "aggregate", "fallback", "pilot_args")
  invalid_names <- setdiff(control_names, valid_names)
  if (length(invalid_names)) {
    stop2(
      "Invalid group-effect autocenter control argument(s): ",
      collapse_comma(invalid_names), "."
    )
  }

  out <- list(
    method = control$method %||% "pathfinder",
    aggregate = control$aggregate %||% "median",
    fallback = control$fallback %||% "error",
    pilot_args = control$pilot_args %||% list()
  )
  out$method <- as_one_character(out$method)
  if (!out$method %in% c("pathfinder", "hmc")) {
    stop2("Autocenter 'method' must be \"pathfinder\" or \"hmc\".")
  }
  out$aggregate <- as_one_character(out$aggregate)
  if (!out$aggregate %in% c("median", "mean")) {
    stop2("Autocenter 'aggregate' must be \"median\" or \"mean\".")
  }
  if (is.character(out$fallback)) {
    out$fallback <- as_one_character(out$fallback)
    if (!out$fallback %in% c("error", "noncentered", "centered")) {
      stop2(
        "Character autocenter 'fallback' must be \"error\", ",
        "\"noncentered\", or \"centered\"."
      )
    }
  } else if (is.numeric(out$fallback)) {
    out$fallback <- as_one_numeric(out$fallback)
    if (!is.finite(out$fallback) ||
        out$fallback < 0 || out$fallback > 1) {
      stop2("Numeric autocenter 'fallback' must be in [0, 1].")
    }
  } else {
    stop2(
      "Autocenter 'fallback' must be \"error\", \"noncentered\", ",
      "\"centered\", or a number in [0, 1]."
    )
  }
  if (!is.list(out$pilot_args)) {
    stop2("Autocenter 'pilot_args' must be a list.")
  }
  pilot_names <- names(out$pilot_args)
  if (length(out$pilot_args) &&
      (is.null(pilot_names) || any(!nzchar(pilot_names)))) {
    stop2("Autocenter 'pilot_args' must be a named list.")
  }
  if (anyDuplicated(pilot_names)) {
    stop2("Autocenter 'pilot_args' names must be unique.")
  }
  class(out) <- c("brms_re_autocenter_control", "list")
  out
}

#' @rdname autocenter_control
#' @param x An object containing centering weights resolved by a precursor.
#' @param ... Reserved for methods.
#' @export
centering_weights <- function(x, ...) {
  UseMethod("centering_weights")
}

#' @export
centering_weights.brms_re_autocenter_summary <- function(x, ...) {
  x$rho
}

#' @export
centering_weights.default <- function(x, ...) {
  stop2("No 'centering_weights' method is available for this object.")
}

# Give resolved values a distinct class so gr() can distinguish a completed
# precursor from an unresolved center = "auto" request or a direct numeric
# partial-centering input.
.as_brmsautocenter_resolved <- function(x) {
  if (is.null(dim(x))) {
    x_names <- names(x)
    x <- matrix(
      as.numeric(x), ncol = 1L,
      dimnames = list(x_names, NULL)
    )
  }
  class(x) <- unique(c("brmsautocenter_resolved", class(x)))
  x
}

# Resolve a validated fallback to the value used for one candidate cell.
.re_autocenter_fallback <- function(fallback, variable) {
  if (identical(fallback, "error")) {
    stop2(
      "No valid precursor draws were available for autocenter candidate '",
      variable, "'."
    )
  }
  if (identical(fallback, "noncentered")) {
    return(0)
  }
  if (identical(fallback, "centered")) {
    return(1)
  }
  stopifnot(is.numeric(fallback), length(fallback) == 1L)
  fallback
}

# Parse the generated-quantity names without depending on a posterior draws
# class. Unrelated columns (for example .chain and .iteration) are ignored.
.re_autocenter_candidate_index <- function(variables) {
  stopifnot(is.character(variables))
  pattern <- paste0(
    "^rho_center_candidate_([A-Za-z0-9_]+)",
    "\\[([1-9][0-9]*),([1-9][0-9]*)\\]$"
  )
  matched <- regexec(pattern, variables, perl = TRUE)
  pieces <- regmatches(variables, matched)
  keep <- lengths(pieces) == 4L
  if (!any(keep)) {
    stop2(
      "No variables named 'rho_center_candidate_ID[j,k]' were found in ",
      "the precursor draws."
    )
  }
  pieces <- pieces[keep]
  out <- data.frame(
    column = which(keep),
    variable = vapply(pieces, `[[`, character(1), 1L),
    id = vapply(pieces, `[[`, character(1), 2L),
    level = as.integer(vapply(pieces, `[[`, character(1), 3L)),
    coefficient = as.integer(vapply(pieces, `[[`, character(1), 4L)),
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(out$variable)) {
    stop2("Autocenter candidate variable names must be unique.")
  }
  out
}

# Aggregate generated-quantity candidates into one fixed G by K matrix per
# framed group ID. Invalid draws are diagnosed and excluded. Values within
# floating-point tolerance of an endpoint are clipped to that endpoint; a
# fallback is used only when a cell has no valid draws at all.
aggregate_re_autocenter_draws <- function(draws, control = NULL) {
  control <- validate_re_autocenter_control(control)
  if (is.null(dim(draws)) || length(dim(draws)) != 2L || !nrow(draws)) {
    stop2("Autocenter precursor 'draws' must be a non-empty matrix-like object.")
  }
  variables <- colnames(draws)
  if (is.null(variables) || any(!nzchar(variables))) {
    stop2("Autocenter precursor 'draws' must have column names.")
  }
  index <- .re_autocenter_candidate_index(variables)
  if (is.data.frame(draws)) {
    candidate_draws <- as.matrix(draws[, index$column, drop = FALSE])
  } else {
    candidate_draws <- as.matrix(draws)[, index$column, drop = FALSE]
  }
  if (!is.numeric(candidate_draws)) {
    stop2("Autocenter candidate draws must be numeric.")
  }

  ids <- unique(index$id)
  rho <- setNames(vector("list", length(ids)), ids)
  diagnostics <- vector("list", nrow(index))
  tolerance <- sqrt(.Machine$double.eps)
  diagnostic_i <- 0L

  for (id in ids) {
    id_index <- index[index$id == id, , drop = FALSE]
    G <- max(id_index$level)
    K <- max(id_index$coefficient)
    expected <- expand.grid(
      level = seq_len(G), coefficient = seq_len(K),
      KEEP.OUT.ATTRS = FALSE
    )
    observed_key <- paste(id_index$level, id_index$coefficient, sep = ":")
    expected_key <- paste(expected$level, expected$coefficient, sep = ":")
    if (anyDuplicated(observed_key) || !setequal(observed_key, expected_key)) {
      stop2(
        "Autocenter candidates for group ID '", id,
        "' must form one complete G by K matrix."
      )
    }
    id_index <- id_index[match(expected_key, observed_key), , drop = FALSE]
    center <- matrix(
      NA_real_, nrow = G, ncol = K,
      dimnames = list(as.character(seq_len(G)), as.character(seq_len(K)))
    )

    for (i in seq_rows(id_index)) {
      values <- candidate_draws[, match(id_index$column[i], index$column)]
      finite <- is.finite(values)
      in_bounds <- finite & values >= -tolerance & values <= 1 + tolerance
      valid <- pmin(1, pmax(0, values[in_bounds]))
      fallback_used <- !length(valid)
      if (fallback_used) {
        estimate <- .re_autocenter_fallback(
          control$fallback, id_index$variable[i]
        )
      } else if (identical(control$aggregate, "median")) {
        estimate <- stats::median(valid)
      } else {
        estimate <- mean(valid)
      }
      estimate <- min(1, max(0, estimate))
      center[id_index$level[i], id_index$coefficient[i]] <- estimate

      finite_values <- values[finite]
      diagnostic_i <- diagnostic_i + 1L
      diagnostics[[diagnostic_i]] <- data.frame(
        id = id,
        level = id_index$level[i],
        coefficient = id_index$coefficient[i],
        variable = id_index$variable[i],
        n_draws = length(values),
        n_finite = sum(finite),
        n_valid = length(valid),
        n_out_of_bounds = sum(finite & !in_bounds),
        raw_min = if (length(finite_values)) min(finite_values) else NA_real_,
        raw_max = if (length(finite_values)) max(finite_values) else NA_real_,
        estimate = estimate,
        fallback_used = fallback_used,
        stringsAsFactors = FALSE
      )
    }
    rho[[id]] <- .as_brmsautocenter_resolved(center)
  }

  diagnostics <- do_call(rbind, diagnostics[seq_len(diagnostic_i)])
  rownames(diagnostics) <- NULL
  structure(
    list(rho = rho, diagnostics = diagnostics, control = control),
    class = c("brms_re_autocenter_summary", "list")
  )
}

# Is an expression a call to gr(), including explicitly namespaced calls?
.is_re_autocenter_gr_call <- function(expr) {
  if (!is.call(expr)) {
    return(FALSE)
  }
  head <- expr[[1L]]
  identical(head, as.name("gr")) ||
    (is.call(head) && length(head) == 3L &&
     as.character(head[[1L]]) %in% c("::", ":::") &&
     identical(as.character(head[[3L]]), "gr"))
}

# Does a center expression evaluate to the automatic-precursor marker?
.is_re_autocenter_marker <- function(expr, envir) {
  value <- tryCatch(eval(expr, envir = envir), error = function(e) NULL)
  is.character(value) && length(value) == 1L &&
    !is.na(value) && identical(value, "auto")
}

# Does a formula still contain an unresolved automatic request? Unlike Stan
# code comparison, this distinguishes a pending precursor from frozen weights.
has_unresolved_re_autocenter_formula <- function(x) {
  walk_formula <- function(formula) {
    formula_env <- environment(formula)
    if (is.null(formula_env)) {
      formula_env <- parent.frame()
    }
    walk_call <- function(expr) {
      if (!is.call(expr)) {
        return(FALSE)
      }
      if (.is_re_autocenter_gr_call(expr)) {
        center_position <- which(names(expr) == "center")
        return(
          length(center_position) == 1L &&
            .is_re_autocenter_marker(expr[[center_position]], formula_env)
        )
      }
      any(vapply(as.list(expr)[-1L], walk_call, logical(1)))
    }
    walk_call(formula)
  }
  walk_object <- function(object) {
    if (is.mvbrmsformula(object)) {
      return(any(vapply(object$forms, walk_object, logical(1))))
    }
    if (is.brmsformula(object)) {
      formulas <- object$pforms
      if (!"mu" %in% names(object$pforms)) {
        formulas <- c(list(object$formula), formulas)
      }
      return(any(vapply(formulas, walk_formula, logical(1))))
    }
    if (is.formula(object)) {
      return(walk_formula(object))
    }
    stop2("Automatic centering requires a formula or brms formula object.")
  }
  walk_object(x)
}

# Replace automatic gr() occurrences, in formula traversal order, by fixed
# numeric objects in private formula environments. IDs are the numeric IDs in
# the fitted reframe. Occurrence numbers disambiguate multiple gr() calls that
# contribute coefficient columns to the same framed ID.
replace_re_autocenter <- function(x, center, id = NULL,
                                  occurrence = NULL) {
  if (is.numeric(center)) {
    center <- list(center)
  }
  if (!is.list(center) || !length(center)) {
    stop2("Autocenter replacements must be a non-empty list of numeric values.")
  }
  center <- lapply(center, function(value) {
    value_dim <- dim(value)
    if (!is.numeric(value) || !length(value) ||
        (!is.null(value_dim) && length(value_dim) != 2L) ||
        anyNA(value) || any(!is.finite(value)) ||
        any(value < 0 | value > 1)) {
      stop2(
        "Each autocenter replacement must be a non-empty numeric vector or ",
        "matrix with finite values in [0, 1]."
      )
    }
    .as_brmsautocenter_resolved(value)
  })
  if (is.null(id)) {
    id <- names(center)
  }
  if (is.null(id) || length(id) != length(center)) {
    stop2("One framed group 'id' is required per autocenter replacement.")
  }
  id_integer <- SW(as.integer(id))
  if (anyNA(id_integer) || any(id_integer < 1L) ||
      any(as.character(id_integer) != as.character(id))) {
    stop2("Autocenter replacement 'id' values must be positive integers.")
  }
  if (is.null(occurrence)) {
    occurrence <- ave(
      seq_along(id_integer), id_integer,
      FUN = function(i) seq_along(i)
    )
  }
  occurrence_integer <- SW(as.integer(occurrence))
  if (length(occurrence_integer) != length(center) ||
      anyNA(occurrence_integer) || any(occurrence_integer < 1L) ||
      any(as.character(occurrence_integer) != as.character(occurrence))) {
    stop2("Autocenter replacement 'occurrence' values must be positive integers.")
  }
  replacement_key <- paste(id_integer, occurrence_integer, sep = ":")
  if (anyDuplicated(replacement_key)) {
    stop2("Framed autocenter ID/occurrence pairs must be unique.")
  }
  symbols <- paste0(
    ".brms_autocenter_id_", id_integer,
    "_occurrence_", occurrence_integer
  )
  used_symbols <- character(length(symbols))

  replacement_i <- 0L
  walk_formula <- function(formula) {
    source_env <- environment(formula)
    if (is.null(source_env)) {
      source_env <- parent.frame()
    }
    private_env <- new.env(parent = source_env)
    seen_terms <- new.env(parent = emptyenv())
    changed <- FALSE
    replace_duplicate <- function(expr, symbol) {
      if (!is.call(expr)) {
        return(expr)
      }
      if (.is_re_autocenter_gr_call(expr)) {
        center_position <- which(names(expr) == "center")
        if (length(center_position) == 1L &&
            .is_re_autocenter_marker(expr[[center_position]], source_env)) {
          expr[[center_position]] <- as.name(symbol)
          changed <<- TRUE
        }
        return(expr)
      }
      for (i in seq_along(expr)[-1L]) {
        expr[[i]] <- replace_duplicate(expr[[i]], symbol)
      }
      expr
    }
    walk_call <- function(expr) {
      if (!is.call(expr)) {
        return(expr)
      }
      head <- as.character(expr[[1L]])
      is_group_term <- length(head) == 1L && head %in% c("|", "||")
      if (is_group_term) {
        term_key <- deparse0(expr)
        if (exists(term_key, envir = seen_terms, inherits = FALSE)) {
          return(replace_duplicate(
            expr, get(term_key, envir = seen_terms, inherits = FALSE)
          ))
        }
        before <- replacement_i
        for (i in seq_along(expr)[-1L]) {
          expr[[i]] <- walk_call(expr[[i]])
        }
        if (replacement_i == before + 1L) {
          assign(term_key, used_symbols[replacement_i], envir = seen_terms)
        }
        return(expr)
      }
      if (.is_re_autocenter_gr_call(expr)) {
        center_position <- which(names(expr) == "center")
        if (length(center_position) > 1L) {
          stop2("A gr() call may contain only one 'center' argument.")
        }
        if (length(center_position) && .is_re_autocenter_marker(
          expr[[center_position]], source_env
        )) {
          replacement_i <<- replacement_i + 1L
          if (replacement_i > length(center)) {
            stop2(
              "More gr(center = \"auto\") occurrences were found than ",
              "autocenter replacements were supplied."
            )
          }
          symbol_base <- symbols[replacement_i]
          symbol <- symbol_base
          suffix <- 1L
          while (exists(symbol, envir = private_env, inherits = TRUE)) {
            suffix <- suffix + 1L
            symbol <- paste0(symbol_base, "_refit_", suffix)
          }
          used_symbols[replacement_i] <<- symbol
          expr[[center_position]] <- as.name(symbol)
          assign(symbol, center[[replacement_i]], envir = private_env)
          changed <<- TRUE
        }
        return(expr)
      }
      for (i in seq_along(expr)[-1L]) {
        expr[[i]] <- walk_call(expr[[i]])
      }
      expr
    }
    formula[] <- lapply(formula, walk_call)
    if (changed) {
      environment(formula) <- private_env
    }
    formula
  }
  walk_object <- function(object) {
    if (is.mvbrmsformula(object)) {
      object$forms <- lapply(object$forms, walk_object)
      return(object)
    }
    if (is.brmsformula(object)) {
      # brmsterms() uses an explicit mu formula in place of the main RHS.
      if (!"mu" %in% names(object$pforms)) {
        object$formula <- walk_formula(object$formula)
      }
      object$pforms <- lapply(object$pforms, walk_formula)
      return(object)
    }
    if (is.formula(object)) {
      return(walk_formula(object))
    }
    stop2("Autocenter replacement requires a formula or brms formula object.")
  }

  out <- walk_object(x)
  if (replacement_i < length(center)) {
    stop2(
      "Fewer gr(center = \"auto\") occurrences were found than autocenter ",
      "replacements were supplied."
    )
  }
  # Reuse the fitted-formula freezer so serialization and update() see the
  # same self-contained coordinate system as direct numeric center inputs.
  materialize_re_center(out)
}

# Reconstruct the automatic request represented by a fitted formula's private
# resolved matrices. This keeps the opt-in update(..., center_control = ...)
# rerun available after an unrelated recompilation.
restore_re_autocenter_request <- function(x) {
  found <- FALSE
  walk_formula <- function(formula) {
    formula_env <- environment(formula)
    if (is.null(formula_env)) {
      formula_env <- parent.frame()
    }
    walk_call <- function(expr) {
      if (!is.call(expr)) {
        return(expr)
      }
      if (.is_re_autocenter_gr_call(expr)) {
        center_position <- which(names(expr) == "center")
        if (length(center_position) == 1L) {
          value <- tryCatch(
            eval(expr[[center_position]], envir = formula_env),
            error = function(error) NULL
          )
          if (inherits(value, "brmsautocenter_resolved")) {
            expr[[center_position]] <- "auto"
            found <<- TRUE
          }
        }
        return(expr)
      }
      for (i in seq_along(expr)[-1L]) {
        expr[[i]] <- walk_call(expr[[i]])
      }
      expr
    }
    formula[] <- lapply(formula, walk_call)
    formula
  }
  walk_object <- function(object) {
    if (is.mvbrmsformula(object)) {
      object$forms <- lapply(object$forms, walk_object)
      return(object)
    }
    if (is.brmsformula(object)) {
      if (!"mu" %in% names(object$pforms)) {
        object$formula <- walk_formula(object$formula)
      }
      object$pforms <- lapply(object$pforms, walk_formula)
      return(object)
    }
    if (is.formula(object)) {
      return(walk_formula(object))
    }
    stop2("Automatic centering requires a formula or brms formula object.")
  }
  out <- walk_object(x)
  if (found) out else NULL
}
