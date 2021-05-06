### Only needed here in the unit tests:
ilink <- brms:::ilink
slice <- brms:::slice
### 

inv_link_cumulative_ch <- function(x, link) {
  x <- ilink(x, link)
  ndim <- length(dim(x))
  dim_noncat <- dim(x)[-ndim]
  ones_arr <- array(1, dim = c(dim_noncat, 1))
  zeros_arr <- array(0, dim = c(dim_noncat, 1))
  abind::abind(x, ones_arr) - abind::abind(zeros_arr, x)
}

inv_link_sratio_ch <- function(x, link) {
  x <- ilink(x, link)
  ndim <- length(dim(x))
  dim_noncat <- dim(x)[-ndim]
  marg_othdim <- seq_along(dim(x))[-ndim]
  ones_arr <- array(1, dim = c(dim_noncat, 1))
  Sx_cumprod <- aperm(apply(1 - x, marg_othdim, cumprod),
                      perm = c(marg_othdim + 1, 1))
  abind::abind(x, ones_arr) * abind::abind(ones_arr, Sx_cumprod)
}

inv_link_cratio_ch <- function(x, link) {
  x <- ilink(x, link)
  ndim <- length(dim(x))
  dim_noncat <- dim(x)[-ndim]
  marg_othdim <- seq_along(dim(x))[-ndim]
  ones_arr <- array(1, dim = c(dim_noncat, 1))
  x_cumprod <- aperm(apply(x, marg_othdim, cumprod),
                     perm = c(marg_othdim + 1, 1))
  abind::abind(1 - x, ones_arr) * abind::abind(ones_arr, x_cumprod)
}

inv_link_acat_ch <- function(x, link) {
  ndim <- length(dim(x))
  dim_noncat <- dim(x)[-ndim]
  marg_othdim <- seq_along(dim(x))[-ndim]
  ones_arr <- array(1, dim = c(dim_noncat, 1))
  if (link == "logit") { 
    # faster evaluation in this case
    exp_x_cumprod <- aperm(apply(exp(x), marg_othdim, cumprod),
                           perm = c(marg_othdim + 1, 1))
    out <- abind::abind(ones_arr, exp_x_cumprod)
  } else {
    x <- ilink(x, link)
    x_cumprod <- aperm(apply(x, marg_othdim, cumprod),
                       perm = c(marg_othdim + 1, 1))
    nthres <- dim(x)[ndim]
    Sx_cumprod_rev <- aperm(apply(
      1 - slice(x, ndim, rev(seq_len(nthres)), drop = FALSE),
      marg_othdim, cumprod
    ), perm = c(marg_othdim + 1, 1))
    Sx_cumprod_rev <- slice(
      Sx_cumprod_rev, ndim, rev(seq_len(nthres)), drop = FALSE
    )
    out <- abind::abind(ones_arr, x_cumprod) *
      abind::abind(Sx_cumprod_rev, ones_arr)
  }
  catsum <- apply(out, marg_othdim, sum)
  sweep(out, marg_othdim, catsum, "/")
}
