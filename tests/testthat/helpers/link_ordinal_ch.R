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

# Very similar to link_cumulative(), but iterates over the observations:
link_cumulative_ch <- function(x, link) {
  # For testing purposes, only allow 3-dimensional arrays here:
  stopifnot(length(dim(x)) == 3)
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

# The same as link_sratio(), but dropping margins:
link_sratio_ch <- function(x, link) {
  ndim <- length(dim(x))
  .F_k <- function(k) {
    if (k == 1) {
      prev_res <- list(F_k = NULL, S_km1_prod = 1)
    } else {
      prev_res <- .F_k(k - 1)
    }
    F_k <- slice(x, ndim, k) / prev_res$S_km1_prod
    return(list(
      F_k = abind::abind(prev_res$F_k, F_k, along = ndim),
      S_km1_prod = prev_res$S_km1_prod * (1 - F_k)
    ))
  }
  x <- .F_k(dim(x)[ndim] - 1)$F_k
  link_ch(x, link)
}
environment(link_sratio_ch) <- as.environment(asNamespace("brms"))

# The same as link_cratio(), but dropping margins:
link_cratio_ch <- function(x, link) {
  ndim <- length(dim(x))
  .F_k <- function(k) {
    if (k == 1) {
      prev_res <- list(F_k = NULL, F_km1_prod = 1)
    } else {
      prev_res <- .F_k(k - 1)
    }
    F_k <- 1 - slice(x, ndim, k) / prev_res$F_km1_prod
    return(list(
      F_k = abind::abind(prev_res$F_k, F_k, along = ndim),
      F_km1_prod = prev_res$F_km1_prod * F_k
    ))
  }
  x <- .F_k(dim(x)[ndim] - 1)$F_k
  link_ch(x, link)
}
environment(link_cratio_ch) <- as.environment(asNamespace("brms"))

# The same as link_acat(), but possibly dropping margins and not treating the
# logit link as a special case:
link_acat_ch <- function(x, link) {
  ndim <- length(dim(x))
  ncat <- dim(x)[ndim]
  dim_noncat <- dim(x)[-ndim]
  x <- slice(x, ndim, -1) / slice(x, ndim, -ncat)
  x <- inv_odds(x)
  array(link_ch(x, link), dim = c(dim_noncat, ncat - 1))
}
environment(link_acat_ch) <- as.environment(asNamespace("brms"))
