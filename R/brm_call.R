#' Build – but **do not run** – a `brm()` call
#'
#' `create_brm_call()` is a wafer-thin wrapper around
#' \code{\link[brms]{brm}} that adds \code{call_only = TRUE}.
#' All other arguments are forwarded unchanged.
#'
#' @param ... Any arguments accepted by \code{\link[brms]{brm}}.
#'   They are stored unmodified in the returned object.  See the
#'   \pkg{brms} help page for the exhaustive list.
#'
#' @return A list of class \code{c("brm_call", "list")}.
#' @seealso \code{\link[brms]{brm}}
#' @export
create_brm_call <- function(...) {
  brm(..., call_only = TRUE)
}

#' @export
print.brm_call <- function(x, ...) {

  # helper: safe extraction by name (case-insensitive)
  get_arg <- function(nm) {
    hit <- match(tolower(nm), tolower(names(x)), nomatch = 0L)
    if (hit) x[[hit]] else NULL
  }

  ## ---- 1. headline fields ---------------------------------------
  formula <- get_arg("formula")
  data    <- get_arg("data")
  family  <- get_arg("family")

  cat("<brm_call>\n")

  if (!is.null(formula)) {
    cat("  Formula : ",
        paste(deparse(formula, nlines = 1L), collapse = " "),
        "\n", sep = "")
  }

  if (!is.null(data)) {
    dclass <- class(data)[1L]
    drows  <- tryCatch(NROW(data), error = function(e) NA_integer_)
    dcols  <- tryCatch(NCOL(data), error = function(e) NA_integer_)
    cat("  Data    : ", dclass,
        " [", drows, " x ", dcols, "]\n", sep = "")
  }

  if (!is.null(family)) {
    fam <- if (inherits(family, "family")) family$family else as.character(family)
    cat("  Family  : ", fam, "\n", sep = "")
  }

  ## ---- 2. everything else ---------------------------------------
  # which names have NOT been printed?
  printed <- c("formula", "data", "family")
  keep    <- !(tolower(names(x)) %in% printed)

  if (any(keep)) {
    cat("  Other arguments (", sum(keep), "):\n", sep = "")
    for (nm in names(x)[keep]) {
      val <- x[[nm]]

      # tiny preview of the value
      summary <- if (is.atomic(val) && length(val) == 1) {
        as.character(val)
      } else if (is.list(val)) {
        paste0("<list:", length(val), ">")
      } else {
        paste0("<", class(val)[1L], ">")
      }

      cat("    - ", nm, " = ", summary, "\n", sep = "")
    }
  }

  invisible(x)
}

#' @export
summary.brm_call <- function(object, ...) {
  cat("Summary of <brm_call>\n")
  str(object, max.level = 1, give.attr = FALSE, no.list = TRUE)
  invisible(object)
}
