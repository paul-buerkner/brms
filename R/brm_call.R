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

#' @noRd
is.brm_call <- function(x) inherits(x, "brm_call")

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
#' Get a package’s version as a string
#'
#' Silently returns `"NA"` when the package is not installed.
#' @noRd
v_package <- function(pkg) {
  if (identical(pkg, "mock")) {
    return("mock-1.0.0")
  }
  tryCatch(
    as.character(utils::packageVersion(pkg)),
    error = function(...) NA_character_
  )
}


# experimental lots of
# TODO here
hash_model_signature <- function(call) {
  # Build *just enough* to get Stan code & bframe ------------------------
  bframe  <- brmsterms(call$formula)
  prior   <- .validate_prior(call$prior, bframe)
  stanvars <- validate_stanvars(call$stanvars, call$stan_funs)

  scode <- .stancode(bframe,
                     prior     = prior,
                     stanvars  = stanvars,
                     backend   = call$backend,
                     threads   = call$threads,
                     opencl    = call$opencl,
                     normalize = call$normalize)

  # Canonicalise: drop comments, collapse whitespace --------------------
  scode <- gsub("^//.*?$|/\\*.*?\\*/", "", scode, perl = TRUE)    # rm comments
  scode <- gsub("[[:space:]]+", " ", scode)                       # one space
  scode <- trimws(scode)

  data = call$data
  if(is.null(data)){
    data_sig <- 'NULL'
  }else{
    data_sig = list(
      n    = NROW(call$data),
      vars = sort(names(call$data))
    )
  }

  ## map backend → package  -------------------------------------------------
  packages_backend <- list(rstan   = "rstan",
                           cmdstanr = "cmdstanr",
                           mock     = "mock")

  backend_pkg      <- packages_backend[[call$backend]]
  backend_version  <- v_package(backend_pkg)
  brms_version     <- v_package("brms")

  versions <- c(brms  = brms_version,
                backend = backend_version)


  compact <- function(x) x[!vapply(x, function(z) length(z) == 0 || is.na(z),
                                   logical(1))]

  payload <- list(
    versions  = versions,
    code_hash = digest::digest(scode, algo = "xxhash64", serialize = FALSE),
    data_sig  = data_sig,
    backend   = call$backend
  )

  payload <- c(payload,
               compact(list(threads  = call$threads,
                            opencl   = call$opencl,
                            normalize = call$normalize)))

  ## final model signature ---------------------------------------------------
  model_hash <- digest::digest(payload, algo = "xxhash64", serialize = TRUE)
  attr(brm_call, "model_hash") <- model_hash


  digest::digest(payload, algo = "xxhash64", serialize = TRUE)
}


