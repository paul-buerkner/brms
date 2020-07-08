# list parameters NOT to be saved by Stan
# @return a vector of parameter names to be excluded
exclude_pars <- function(x, ...) {
  UseMethod("exclude_pars")
}

#' @export
exclude_pars.default <- function(x, ...) {
  character(0)
}

# @param save_ranef save varying effects per level?
# @param save_mevars save noise-free variables?
# @param save_all_pars save all variables from the 'parameters' block?
#' @export
exclude_pars.brmsfit <- function(x, save_ranef = TRUE, save_mevars = FALSE,
                                 save_all_pars = FALSE, ...) {
  save_ranef <- as_one_logical(save_ranef)
  save_mevars <- as_one_logical(save_mevars)
  save_all_pars <- as_one_logical(save_all_pars)
  out <- character(0)
  bterms <- brmsterms(x$formula)
  c(out) <- exclude_pars(
    bterms, data = x$data, save_ranef = save_ranef, 
    save_all_pars = save_all_pars, save_mevars = save_mevars, ...
  )
  meef <- tidy_meef(bterms, x$data)
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
  ranef <- x$ranef
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
  out <- unique(out)
  att <- nlist(save_ranef, save_mevars, save_all_pars)
  attributes(out)[names(att)] <- att
  out
}

#' @export
exclude_pars.mvbrmsterms <- function(x, save_all_pars = FALSE, ...) {
  out <- c("Rescor", "Sigma")
  if (!save_all_pars) {
    c(out) <- c("Lrescor", "LSigma")
  }
  for (i in seq_along(x$terms)) {
    c(out) <- exclude_pars(x$terms[[i]], save_all_pars = save_all_pars, ...)
  }
  out
}

#' @export
exclude_pars.brmsterms <- function(x, save_ranef = TRUE, save_mevars = FALSE, 
                                   save_all_pars = FALSE, ...) {
  out <- character(0)
  p <- usc(combine_prefix(x))
  if (!save_all_pars) {
    par_classes <- c("ordered_Intercept", "fixed_Intercept", "theta")
    c(out) <- paste0(par_classes, p)
  }
  for (dp in names(x$dpars)) {
    c(out) <- exclude_pars(x$dpars[[dp]], save_all_pars = save_all_pars, ...)
  }
  for (nlp in names(x$nlpars)) {
    c(out) <- exclude_pars(x$nlpars[[nlp]], save_all_pars = save_all_pars, ...)
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
exclude_pars.btl <- function(x, data, save_all_pars = FALSE, ...) {
  out <- character(0)
  p <- usc(combine_prefix(x))
  c(out) <- paste0("chol_cor", p)
  if (!save_all_pars) {
    par_classes <- c(
      "bQ", "hs_global", "hs_local", "hs_slab", "zb", "hs_localsp", 
      "zbsp", "Intercept", "first_Intercept", "merged_Intercept",
      "zcar", "nszcar", "zerr"
    )
    c(out) <- paste0(par_classes, p)
    smef <- tidy_smef(x, data)
    for (i in seq_rows(smef)) {
      nb <- seq_len(smef$nbases[i])
      c(out) <- paste0("zs", p, "_", i, "_", nb)
    } 
  }
  out
}
