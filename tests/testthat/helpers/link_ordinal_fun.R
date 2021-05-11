link_ch <- function(x, link) {
  # switch() would be more straightforward, but for testing purposes, use if ()
  # here:
  if (link == "logit") {
    return(qlogis(x))
  } else if (link == "probit") {
    return(qnorm(x))
  } else if (link == "cauchit") {
    return(qcauchy(x))
  } else if (link == "cloglog") {
    return(log(-log(1 - x)))
  } else {
    stop("Unknown link.")
  }
}

link_cumulative_ch <- function(x, link) {
  # For testing purposes, only allow 3-dimensional arrays here:
  stopifnot(length(dim(x)) == 3)
  # For testing purposes, iterate over the observations:
  ndraws <- dim(x)[1]
  nobsv <- dim(x)[2]
  ncat <- dim(x)[3]
  x_cumsum <- aperm(
    array(
      sapply(seq_len(nobsv), function(i) {
        apply(x[, i, -ncat, drop = FALSE], 1, cumsum)
      }, simplify = "array"),
      dim = c(ncat - 1, ndraws, nobsv)
    ),
    perm = c(2, 3, 1)
  )
  link_ch(x_cumsum, link = link)
}

link_sratio_ch <- function(x, link) {
  # TODO
}

link_cratio_ch <- function(x, link) {
  # TODO
}

link_acat_ch <- function(x, link) {
  # TODO
}
