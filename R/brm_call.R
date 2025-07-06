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
  c1$model_hash == c2$model_hash
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

#' Get a package’s version as a string
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

#' local_digest lots of
#' @noRd
local_digest <- function(object, algo = "xxhash64", serialize = TRUE) {
  require_package('digest')
  digest::digest(object, algo = algo, serialize = serialize)
}

#' hash_model_signature
#' currently takes into account stancode produced, data hash with a threshold
#' and versions of backend and brms packages.
#' @noRd
hash_model_signature <- function(call) {
  if(!call$file_auto)
    return(call)
  data_row_threshold <- 1e5
  if(call$backend == 'mock'){
    call$model_hash <- 'mock_hash'
    return(call)
  }
  fit_args <- .create_fit_args(call)
  # check for model code in mock case or empty case maybe not able to produce
  if( !is.list(fit_args) || is.null(fit_args$model)){
    return(call)
  }
  if(call$backend == 'rstan') {
    scode <- fit_args$model@model_code
  } else { # cmdstanr
    scode <- fit_args$model$code()
  }
  scode <- paste(scode, collapse="\n")
  data = call$data
  data_sig <- "NULL"
  if(is.data.frame(data)){
    nrows <- NROW(data)
    if(nrows < data_row_threshold) {
      data_sig <- local_digest(data)
    }else{
      data_sig = list(
        n    = NROW(call$data),
        vars = sort(names(call$data))
      )
    }
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
  payload <- list(
    versions  = versions,
    code_hash = local_digest(scode, algo = "xxhash64", serialize = FALSE),
    data_sig  = data_sig,
    backend   = call$backend
  )
  payload <- c(payload,
               list(threads  = call$threads,
                            opencl   = call$opencl,
                            normalize = call$normalize))
  ## final model signature ---------------------------------------------------
  model_hash <- local_digest(payload, algo = "xxhash64", serialize = TRUE)
  call$model_hash <- model_hash
  call$file <- paste0("cache_brmsfit_", model_hash, ".rds")
  call$file_refit <- 'on_change'
  call
}
