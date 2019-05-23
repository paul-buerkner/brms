# list parameters NOT to be saved by Stan
# @param bterms object of class brmsterms
# @param data data passed by the user
# @param ranef output of tidy_ranef
# @param save_ranef save varying effects per level?
# @param save_mevars save noise-free variables?
# @param save_all_pars save all variables from the 'parameters' block?
# @return a vector of parameters to be excluded
exclude_pars <- function(bterms, data = NULL, ranef = empty_ranef(),
                         save_ranef = TRUE, save_mevars = FALSE,
                         save_all_pars = FALSE) {
  save_ranef <- as_one_logical(save_ranef)
  save_mevars <- as_one_logical(save_mevars)
  save_all_pars <- as_one_logical(save_all_pars)
  out <- exclude_pars_internal(
    bterms, data = data, save_ranef = save_ranef, 
    save_all_pars = save_all_pars, save_mevars = save_mevars
  )
  meef <- tidy_meef(bterms, data)
  if (nrow(meef)) {
    I <- seq_along(unique(meef$grname))
    K <- seq_rows(meef)
    c(out) <- paste0(c("Xme", "Corme_"), I)
    if (!save_all_pars) {
      c(out) <- c(paste0("zme_", K), paste0("Lme_", I))
    }
    if (!save_mevars) {
      c(out) <- paste0("Xme_", K)
    }
  }
  if (has_rows(ranef)) {
    rm_re_pars <- c(if (!save_all_pars) c("z", "L"), "Cor", "r")
    for (id in unique(ranef$id)) {
      c(out) <- paste0(rm_re_pars, "_", id)
    }
    if (!save_ranef) {
      p <- usc(combine_prefix(ranef))
      c(out) <- paste0("r_", ranef$id, p, "_", ranef$cn)
    }
    tranef <- get_dist_groups(ranef, "student")
    if (!save_all_pars && has_rows(tranef)) {
      c(out) <- paste0(c("udf_", "dfm_"), tranef$ggn)
    }
  }
  att <- nlist(save_ranef, save_mevars, save_all_pars)
  do_call(structure, c(list(unique(out)), att))
}

exclude_pars_internal <- function(x, ...) {
  UseMethod("exclude_pars_internal")
}

#' @export
exclude_pars_internal.default <- function(x, ...) {
  NULL
}

#' @export
exclude_pars_internal.mvbrmsterms <- function(x, save_all_pars, ...) {
  out <- c("Rescor", "Sigma")
  if (!save_all_pars) {
    c(out) <- c("Lrescor", "LSigma")
  }
  for (i in seq_along(x$terms)) {
    c(out) <- exclude_pars_internal(x$terms[[i]], save_all_pars, ...)
  }
  out
}

#' @export
exclude_pars_internal.brmsterms <- function(x, save_ranef, save_all_pars, 
                                            save_mevars, ...) {
  p <- usc(combine_prefix(x))
  out <- paste0(c("chol_cov", names(x$dpars)), p)
  if (!save_all_pars) {
    c(out) <- c(
      paste0("temp", p, "_Intercept1"), 
      paste0("ordered", p, "_Intercept"),
      paste0("fixed", p, "_Intercept"),
      paste0(c("theta", "zcar", "zerr"), p)
    )
    for (dp in names(x$dpars)) {
      c(out) <- exclude_pars_internal(x$dpars[[dp]], ...)
    }
    for (nlp in names(x$nlpars)) {
      c(out) <- exclude_pars_internal(x$nlpars[[nlp]], ...)
    }
  }
  if (!save_mevars && is.formula(x$adforms$mi)) {
    c(out) <- paste0("Yl", p)
  }
  if (!save_ranef) {
    # latent residuals are like group-level effects
    c(out) <- paste0("err", p)
  }
  out
}

#' @export
exclude_pars_internal.btl <- function(x, data, ...) {
  p <- usc(combine_prefix(x))
  out <- c("bQ", "hs_global", "hs_local", "zb", "hs_localsp", "zbsp")
  out <- paste0(out, p)
  c(out) <- paste0("temp", p, "_Intercept")
  smef <- tidy_smef(x, data)
  for (i in seq_rows(smef)) {
    nb <- seq_len(smef$nbases[i])
    c(out) <- paste0("zs", p, "_", i, "_", nb)
  }
  out
}
