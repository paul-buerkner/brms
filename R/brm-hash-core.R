#' Normalize data when NULL
#'
#' If `data` is NULL, return a fallback dataset used for testing.
#' Otherwise, return the original `data` as-is.
get_same_data_if_null <- function(data) {
  if (is.null(data)) {
    # Use a subset of epilepsy as a fallback test dataset
    correct_data <- epilepsy[-c(1, 8, 11, 25, 29, 12), ]
  } else {
    correct_data <- data
  }
  correct_data
}

#' Hash or fingerprint a data.frame for caching
#'
#' * **mode = "auto"**   (default)
#'   ‒ Hashes the full data if `nrow * ncol <= threshold`;
#'   ‒ otherwise hashes only its dimensions to stay fast.
#'
#' * **mode = "signature"**
#'   ‒ Creates a compact, human-readable signature (row/col count, names, types, small
#'      value summaries) and hashes that signature.
#'
#' * **mode = "full"**
#'   ‒ Always hashes the complete data.frame (may be slow for very large data).
#'
#' @param data        A data.frame to fingerprint.
#' @param threshold   Cell count above which the *auto* mode falls back to dim-only hashing.
#' @param mode        One of `"auto"`, `"signature"`, or `"full"`.
#' @param algo        Hash algorithm passed to \link[digest]{digest}.
#' @param include_summary  Logical; if `TRUE` the signature includes basic
#'   column summaries (used only when `mode = "signature"`).
#'
#' @return A character hash.
#' @export
hash_data <- function(data,
                      threshold = 1e7,
                      mode = c("auto", "signature", "full"),
                      algo = "xxhash64",
                      include_summary = FALSE) {

  if (!is.data.frame(data)) {
    stop("hash_data() expects a data.frame.")
  }

  require_package('digest')
  mode <- match.arg(mode)

  ## -- Helper: compact signature ------------------------------------------
  data_signature <- function(df) {
    list(
      nrow     = nrow(df),
      ncol     = ncol(df),
      colnames = names(df),
      types    = vapply(df, function(col) class(col)[1], character(1)),
      summary  = if (include_summary) {
        lapply(df, function(col) {
          if (is.numeric(col)) {
            round(mean(col, na.rm = TRUE), 4)
          } else if (is.factor(col) || is.character(col)) {
            uniq <- unique(as.character(col))
            uniq[1:min(3, length(uniq))]   # first few values/levels
          } else {
            NULL
          }
        })
      } else NULL
    )
  }

  ## -- Decide what to hash -------------------------------------------------
  if (mode == "full") {
    return(digest::digest(data, algo = algo))
  }

  if (mode == "signature") {
    sig <- data_signature(data)
    return(digest::digest(sig, algo = algo))
  }

  ## mode == "auto"
  size <- nrow(data) * ncol(data)
  if (size > threshold) {
    return(digest::digest(dim(data), algo = algo))
  }

  digest::digest(data, algo = algo)
}



#' #' Remove attached environments from object (recursive)
#' #'
#' #' Cleans up environment-related attributes from nested structures
#' #' to ensure consistent hashing and avoid environment pollution.
#' remove_env_attrs <- function(obj) {
#'   # Remove .Environment attribute if present
#'   if (!is.null(attr(obj, ".Environment"))) {
#'     attr(obj, ".Environment") <- NULL
#'   }
#'   # Recursively handle lists and pairlists
#'   if (is.list(obj) || is.pairlist(obj)) {
#'     obj <- lapply(obj, remove_env_attrs)
#'   }
#'   # Remove environment from formulas
#'   if (inherits(obj, "formula")) {
#'     environment(obj) <- NULL
#'   }
#'   obj
#' }