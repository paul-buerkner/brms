### Only needed here in the unit tests (and only for testing in R CMD check):
seq_cols <- brms:::seq_cols
slice_col <- brms:::slice_col
### 

# Very similar to insert_refcat(), but iterates over the observations (if
# necessary):
insert_refcat_ch <- function(eta, family) {
  ndim <- length(dim(eta))
  if (ndim == 2) {
    return(insert_refcat_ch_i(eta, family = family))
  } else if (ndim == 3) {
    out <- abind::abind(lapply(seq_cols(eta), function(i) {
      insert_refcat_ch_i(slice_col(eta, i), family = family)
    }), along = 3)
    return(aperm(out, perm = c(1, 3, 2)))
  } else {
    stop2("eta has wrong dimensions.")
  }
}

# A matrix-only variant of insert_refcat() (used to be insert_refcat() before it
# was extended to arrays):
insert_refcat_ch_i <- function(eta, family) {
  stopifnot(is.matrix(eta), is.brmsfamily(family))
  if (!conv_cats_dpars(family) || isNA(family$refcat)) {
    return(eta)
  }
  # need to add zeros for the reference category
  zeros <- as.matrix(rep(0, nrow(eta)))
  if (is.null(family$refcat) || is.null(family$cats)) {
    # no information on the categories provided:
    # use the first category as the reference
    return(cbind(zeros, eta))
  }
  colnames(zeros) <- paste0("mu", family$refcat)
  iref <- match(family$refcat, family$cats)
  before <- seq_len(iref - 1)
  after <- setdiff(seq_cols(eta), before)
  cbind(eta[, before, drop = FALSE], zeros, eta[, after, drop = FALSE])
}
