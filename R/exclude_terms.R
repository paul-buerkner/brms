# exclude predictor terms from being evaluated
exclude_terms <- function(x, ...) {
  UseMethod("exclude_terms")
}

#' @export
exclude_terms.brmsfit <- function(x, ...) {
  x$formula <- exclude_terms(x$formula, ...)
  x
}

#' @export
exclude_terms.mvbrmsformula <- function(x, ...) {
  for (i in seq_along(x$forms)) {
    x$forms[[i]] <- exclude_terms(x$forms[[i]], ...)
  }
  x
}

#' @export
exclude_terms.brmsformula <- function(
  x, excl_term_types = NULL, incl_autocor = TRUE,
  smooths_only = FALSE, offset = TRUE, ...
) {
  excl_term_types <- as.character(excl_term_types)
  # TODO: deprecate the three arguments below?
  incl_autocor <- as_one_logical(incl_autocor)
  smooths_only <- as_one_logical(smooths_only)
  offset <- as_one_logical(offset)
  if (!incl_autocor) {
    c(excl_term_types) <- "ac"
  }
  if (!offset) {
    c(excl_term_types) <- "offset"
  }
  if (smooths_only) {
    excl_term_types <- setdiff(all_term_types(), "sm")
  }
  if (!length(excl_term_types)) {
    return(x)
  }
  invalid_types <- setdiff(excl_term_types, all_term_types())
  if (length(invalid_types)) {
    stop2("The following term types are invalid: ",
          collapse_comma(invalid_types))
  }
  attr(x$formula, "excl_term_types") <- excl_term_types
  for (i in seq_along(x$pforms)) {
    attr(x$pforms[[i]], "excl_term_types") <- excl_term_types
  }
  x
}

# extract names of excluded term types
excluded_term_types <- function(x) {
  as.character(attr(x, "excl_term_types", TRUE))
}
