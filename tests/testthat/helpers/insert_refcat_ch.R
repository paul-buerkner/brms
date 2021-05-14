insert_refcat_ch <- function(eta, family) {
  stopifnot(is.array(eta), is.brmsfamily(family))
  if (!conv_cats_dpars(family) || isNA(family$refcat)) {
    return(eta)
  }
  # need to add zeros for the reference category
  ndim <- length(dim(eta))
  dim_noncat <- dim(eta)[-ndim]
  zeros_arr <- array(0, dim = c(dim_noncat, 1))
  if (is.null(family$refcat) || is.null(family$cats)) {
    # no information on the categories provided:
    # use the first category as the reference
    return(abind::abind(zeros_arr, eta))
  }
  dimnames(zeros_arr)[[ndim]] <- paste0("mu", family$refcat)
  iref <- match(family$refcat, family$cats)
  before <- seq_len(iref - 1)
  after <- setdiff(seq_dim(eta, ndim), before)
  abind::abind(slice(eta, ndim, before, drop = FALSE),
               zeros_arr,
               slice(eta, ndim, after, drop = FALSE))
}
