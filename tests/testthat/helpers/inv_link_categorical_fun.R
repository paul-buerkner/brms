inv_link_categorical_ch <- function(x, log = FALSE) {
  ndim <- length(dim(x))
  # For testing purposes, only allow 3-dimensional arrays here:
  if (ndim <= 1) {
    x <- array(x, dim = c(1, 1, length(x)))
  } else if (ndim == 2) {
    x <- array(x, dim = c(dim(x)[1], 1, dim(x)[2]))
  } else if (ndim > 3) {
    stop("At most 3 dimensions are allowed here.")
  }
  # For testing purposes, iterate over the observations:
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
  # Quick-and-dirty solution to drop the margin for a single observation:
  if (nobsv == 1) {
    return(slice(out, 2, 1))
  }
  out
}
