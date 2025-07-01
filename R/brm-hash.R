# hash_brms_call
#
# Computes a reproducible hash for a given brms call argument list.
# Designed to support caching by ensuring the same statistical inputs
# always return the same hash, regardless of argument order or environments.
#
# Args:
#   args_list: A named list of arguments to brm(), such as formula, data, family, prior, etc.
#   algo:      Digest algorithm to use (default: "xxhash64" for speed and reliability).
#   data_policy: Controls how the 'data' element is included in the hash.
#                "names"     → only column names and row count.
#                "signature" → object size of the data.
#                "hash"      → full digest of the data object.
#                "none"      → ignore data in the hash completely.
#
# Returns:
#   A character string representing the hash.
hash_brms_call <- function(args_list,
                           algo        = "xxhash64",
                           data_policy = c("names", "signature", "hash", "none")) {

  require_package("digest")

  ## ─── 2.  Sanitiser (recursive) ──────────────────────────────────────────
  clean <- function(x) {
    if (inherits(x, "formula")) {
      environment(x) <- emptyenv()
      return(as.character(x))
    }

    if (inherits(x, "brmsformula")) {
      # Strip environments from its internal formulas
      x$formula <- clean(x$formula)
      if (!is.null(x$pforms)) {
        x$pforms <- lapply(x$pforms, clean)
      }
      if (!is.null(x$nlpars)) {
        x$nlpars <- sort(x$nlpars)  # to be safe
      }
      return(x)
    }

    if (inherits(x, "mvbrmsformula")) {
      # Multi-response models: list of brmsformula
      x$forms <- lapply(x$forms, clean)
      return(x)
    }

    if (inherits(x, "family")) {
      return(list(family = x$family, link = x$link))
    }

    if (is.function(x)) {
      return(deparse(body(x), width.cutoff = 500L))
    }

    if (is.call(x) || is.language(x) || is.expression(x)) {
      return(deparse(x, width.cutoff = 500L))
    }

    if (is.list(x) && !inherits(x, "data.frame")) {
      x <- x[sort(names(x))]
      return(lapply(x, clean))
    }

    x
  }

  ## ─── 3.  Apply sanitiser to everything ──────────────────────────────────
  args_list <- lapply(args_list, clean)

  ## ─── 4.  Optional data fingerprint ──────────────────────────────────────
  policy <- match.arg(data_policy)
  if ("data" %in% names(args_list) && !is.null(args_list$data) && policy != "none") {
    d <- args_list$data
    args_list$data <- switch(
      policy,
      names      = list(colnames = names(d), nrow = nrow(d)),
      signature  = utils::object.size(d),
      hash       = digest::digest(d, algo = algo, serialize = TRUE)
    )
  }

  ## ─── 5.  Final digest ───────────────────────────────────────────────────
  digest::digest(args_list, algo = algo, serialize = TRUE)
}


# hash_dots
#
# Convenience wrapper for hash_brms_call that accepts ... instead of a named list.
# Useful when manually specifying formula, data, family, etc.
#
# Example:
#   hash_dots(formula = y ~ x, data = df, family = gaussian())
hash_dots  <- function(...){
  dots <- nlist(...)
  hash_brms_call( dots  )

}
