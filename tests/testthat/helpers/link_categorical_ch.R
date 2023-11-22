# Very similar to link_categorical(), but iterates over the observations:
link_categorical_ch <- function(x, refcat = 1, return_refcat = FALSE) {
  # For testing purposes, only allow 3-dimensional arrays here:
  stopifnot(length(dim(x)) == 3)
  x_tosweep <- if (return_refcat) {
    x
  } else {
    slice(x, 3, -refcat, drop = FALSE)
  }
  ndraws <- dim(x)[1]
  nobsv <- dim(x)[2]
  ncat <- dim(x)[3]
  log(aperm(
    array(
      sapply(seq_len(nobsv), function(i) {
        slice(x_tosweep, 2, i) / slice(slice(x, 2, i), 2, refcat)
      }, simplify = "array"),
      dim = c(ndraws, ncat - !return_refcat, nobsv)
    ),
    perm = c(1, 3, 2)
  ))
}
environment(link_categorical_ch) <- as.environment(asNamespace("brms"))
