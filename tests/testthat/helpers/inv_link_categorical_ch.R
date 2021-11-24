# Very similar to inv_link_categorical(), but iterates over the observations and
# always assumes the first category to be the reference category:
inv_link_categorical_ch <- function(x, log = FALSE, refcat_ins = TRUE) {
  if (refcat_ins) {
    zeros_arr <- array(0, dim = c(head(dim(x), -1), 1))
    x <- abind::abind(zeros_arr, x)
  }
  ndim <- length(dim(x))
  # For testing purposes, only allow 3-dimensional arrays here:
  if (ndim <= 1) {
    x <- array(x, dim = c(1, 1, length(x)))
    ndim <- length(dim(x))
    need_drop <- TRUE
  } else if (ndim == 2) {
    x <- array(x, dim = c(dim(x)[1], 1, dim(x)[2]))
    ndim <- length(dim(x))
    need_drop <- TRUE
  } else if (ndim > 3) {
    stop("At most 3 dimensions are allowed here.")
  } else {
    need_drop <- FALSE
  }
  ndraws <- dim(x)[1]
  nobsv <- dim(x)[2]
  ncat <- dim(x)[3]
  .softmax <- if (log) {
    log_softmax
  } else {
    softmax
  }
  out <- aperm(
    array(
      sapply(seq_len(nobsv), function(i) {
        .softmax(slice(x, 2, i))
      }, simplify = "array"),
      dim = c(ndraws, ncat, nobsv)
    ),
    perm = c(1, 3, 2)
  )
  # Quick-and-dirty solution to drop the margin for a single observation (but
  # only if the input object was not a 3-dimensional array):
  if (need_drop) {
    return(slice(out, 2, 1))
  }
  out
}
environment(inv_link_categorical_ch) <- as.environment(asNamespace("brms"))
