#' Build – but **do not run** – a `brm()` call
#'
#' `create_brm_call()` is a wafer-thin wrapper around
#' \code{\link[brms]{brm}} that adds \code{call_only = TRUE}.
#' All other arguments are forwarded unchanged.
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

#' Collect `brm()` arguments into a tidy **brm_call** object
#'
#' Internal helper used at the very top of `brm()`.
#' It separates *formal* `brm()` arguments from any extra
#' Stan-backend tuning options that a user might pass through `...`,
#' then stores everything in a lightweight list with class
#' **`brm_call`**.
#' 
#' *Implementation notes*
#' * We grab the names of the **current** `brm()` formals at run-time
#'   (`names(formals(brms::brm))`) so the helper automatically stays in
#'   sync with upstream changes in **brms**.
#' * Any argument not in that set is treated as an
#'   *extra* Stan argument and saved under `dot_args`.
#'
#' @param ... Arguments passed from the public `brm()` wrapper.
#' @return A list of class `c("brm_call", "list")`.
#' @noRd
.create_brm_call <- function(...) {
  call_env  <- parent.frame()
  arg_names <- names(formals(brm))
  ## 1. drop the literal "..." from arg_names
  arg_names <- arg_names[arg_names != "..."]
  ## 2. capture every formal (already evaluated inside brm())
  brm_call <- setNames(
    lapply(arg_names, function(a) get(a, envir = call_env)),
    arg_names
  )
  ## 3. stash the dot-args for later splicing
  brm_call$dot_args <- list(...)
  class(brm_call) <- c("brm_call" , "list")
  brm_call
}

#' Checks if argument is a \code{brm_call} object
#'
#' @param x An \R object
#'
#' @export
is.brm_call <- function(x) {
  inherits(x, "brm_call")
}

#' Compare two brm_call if they are identical using hash value whcih was created
#' @noRd
identical_brm_calls <- function(c1, c2){
  if(isFALSE(is.brm_call(c1) && is.brm_call(c2)))  {
    stop2("Cannot compare other types than `brm_call`", .subclass = "brms_invalid_brm_call")
  }
  if(isTRUE(is.null(c1$model_hash) || is.null(c2$model_hash))){
    # warning2("Using `all.equal.list` fallback because `model_hash` was not computed.")
    return(all.equal.list(c1, c2))
  }
  # for this PR since we have not implemented hash methods yet
  # fall back to list check
  return(all.equal.list(c1, c2))
}

#' Checks if two brm_call objects are equal
#'
#' @param target A \code{brm_call} object
#' @param current Another \code{brm_call} object
#' @param ... Currently ignored
#'
#' @export
all.equal.brm_call <- function(target, current, ...) {
  identical_brm_calls(target, current)
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
  utils::str(object, max.level = 1, give.attr = FALSE, no.list = TRUE)
  invisible(object)
}
