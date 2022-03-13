inv_link_cumulative_ch <- function(x, link) {
  x <- inv_link(x, link)
  ndim <- length(dim(x))
  ncat <- dim(x)[ndim] + 1
  out <- vector("list", ncat)
  out[[1]] <- slice(x, ndim, 1)
  if (ncat > 2) {
    .diff <- function(k) {
      slice(x, ndim, k) - slice(x, ndim, k - 1)
    }
    mid_cats <- 2:(ncat - 1)
    out[mid_cats] <- lapply(mid_cats, .diff)
  }
  out[[ncat]] <- 1 - slice(x, ndim, ncat - 1)
  abind::abind(out, along = ndim)
}
environment(inv_link_cumulative_ch) <- as.environment(asNamespace("brms"))

inv_link_sratio_ch <- function(x, link) {
  x <- inv_link(x, link)
  ndim <- length(dim(x))
  ncat <- dim(x)[ndim] + 1
  marg_noncat <- seq_along(dim(x))[-ndim]
  out <- vector("list", ncat)
  out[[1]] <- slice(x, ndim, 1)
  if (ncat > 2) {
    .condprod <- function(k) {
      slice(x, ndim, k) *
        apply(1 - slice(x, ndim, 1:(k - 1), drop = FALSE), marg_noncat, prod)
    }
    mid_cats <- 2:(ncat - 1)
    out[mid_cats] <- lapply(mid_cats, .condprod)
  }
  out[[ncat]] <- apply(1 - x, marg_noncat, prod)
  abind::abind(out, along = ndim)
}
environment(inv_link_sratio_ch) <- as.environment(asNamespace("brms"))

inv_link_cratio_ch <- function(x, link) {
  x <- inv_link(x, link)
  ndim <- length(dim(x))
  ncat <- dim(x)[ndim] + 1
  marg_noncat <- seq_along(dim(x))[-ndim]
  out <- vector("list", ncat)
  out[[1]] <- 1 - slice(x, ndim, 1)
  if (ncat > 2) {
    .condprod <- function(k) {
      (1 - slice(x, ndim, k)) *
        apply(slice(x, ndim, 1:(k - 1), drop = FALSE), marg_noncat, prod)
    }
    mid_cats <- 2:(ncat - 1)
    out[mid_cats] <- lapply(mid_cats, .condprod)
  }
  out[[ncat]] <- apply(x, marg_noncat, prod)
  abind::abind(out, along = ndim)
}
environment(inv_link_cratio_ch) <- as.environment(asNamespace("brms"))

inv_link_acat_ch <- function(x, link) {
  ndim <- length(dim(x))
  ncat <- dim(x)[ndim] + 1
  marg_noncat <- seq_along(dim(x))[-ndim]
  out <- vector("list", ncat)
  if (link == "logit") {
    # faster evaluation in this case
    out[[1]] <- array(1, dim = dim(x)[-ndim])
    out[[2]] <- exp(slice(x, ndim, 1))
    if (ncat > 2) {
      .catsum <- function(k) {
        exp(apply(slice(x, ndim, 1:(k - 1), drop = FALSE), marg_noncat, sum))
      }
      remaincats <- 3:ncat
      out[remaincats] <- lapply(remaincats, .catsum)
    }
  } else {
    x <- inv_link(x, link)
    out[[1]] <- apply(1 - x, marg_noncat, prod)
    if (ncat > 2) {
      .othercatprod <- function(k) {
        apply(slice(x, ndim, 1:(k - 1), drop = FALSE), marg_noncat, prod) *
          apply(slice(1 - x, ndim, k:(ncat - 1), drop = FALSE), marg_noncat, prod)
      }
      mid_cats <- 2:(ncat - 1)
      out[mid_cats] <- lapply(mid_cats, .othercatprod)
    }
    out[[ncat]] <- apply(x, marg_noncat, prod)
  }
  out <- abind::abind(out, along = ndim)
  catsum <- apply(out, marg_noncat, sum)
  sweep(out, marg_noncat, catsum, "/")
}
environment(inv_link_acat_ch) <- as.environment(asNamespace("brms"))
