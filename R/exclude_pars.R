exclude_pars <- function(bterms, data = NULL, ranef = empty_ranef(),
                         save_ranef = TRUE, save_mevars = FALSE,
                         save_all_pars = FALSE) {
  # list parameters NOT to be saved by Stan
  # Args:
  #   bterms: object of class brmsterms
  #   data: data passed by the user
  #   ranef: output of tidy_ranef
  #   save_ranef: should group-level effects be saved?
  #   save_mevars: should samples of noise-free variables be saved?
  # Returns:
  #   a vector of parameters to be excluded
  save_ranef <- as_one_logical(save_ranef)
  save_mevars <- as_one_logical(save_mevars)
  save_all_pars <- as_one_logical(save_all_pars)
  out <- exclude_pars_internal(
    bterms, data = data, save_all_pars = save_all_pars,
    save_mevars = save_mevars
  )
  meef <- tidy_meef(bterms, data)
  if (nrow(meef)) {
    I <- seq_along(unique(meef$grname))
    K <- seq_len(nrow(meef))
    out <- c(out, paste0(c("Xme", "Corme_"), I))
    if (!save_all_pars) {
      out <- c(out, paste0("zme_", K), paste0("Lme_", I))
    }
    if (!save_mevars) {
      out <- c(out, paste0("Xme_", K))
    }
  }
  if (nrow(ranef)) {
    rm_re_pars <- c(if (!save_all_pars) c("z", "L"), "Cor", "r")
    for (id in unique(ranef$id)) {
      out <- c(out, paste0(rm_re_pars, "_", id))
    }
    if (!save_ranef) {
      p <- usc(combine_prefix(ranef))
      out <- c(out, paste0("r_", ranef$id, p, "_", ranef$cn))
    }
  }
  att <- nlist(save_ranef, save_mevars, save_all_pars)
  do.call(structure, c(list(unique(out)), att))
}

exclude_pars_internal <- function(x, ...) {
  UseMethod("exclude_pars_internal")
}

#' @export
exclude_pars_internal.mvbrmsterms <- function(x, save_all_pars, ...) {
  out <- c("Rescor", "Sigma")
  if (!save_all_pars) {
    out <- c(out, "Lrescor", "LSigma")
  }
  for (i in seq_along(x$terms)) {
    out <- c(out, exclude_pars_internal(x$terms[[i]], save_all_pars, ...))
  }
  out
}

#' @export
exclude_pars_internal.brmsterms <- function(x, save_all_pars, save_mevars, ...) {
  p <- usc(combine_prefix(x))
  out <- paste0(c("res_cov_matrix", names(x$dpars)), p)
  if (!save_all_pars) {
    out <- c(out,
      paste0("temp", p, "_Intercept1"), 
      paste0("ordered", p, "_Intercept"),
      paste0(c("theta", "zcar"), p)
    )
    for (dp in names(x$dpars)) {
      out <- c(out, exclude_pars_internal(x$dpars[[dp]], ...))
    }
  }
  if (!save_mevars && is.formula(x$adforms$mi)) {
    out <- c(out, paste0("Yl", p)) 
  }
  out
}

#' @export
exclude_pars_internal.btnl <- function(x, ...) {
  out <- NULL
  for (nlp in names(x$nlpars)) {
    out <- c(out, exclude_pars_internal(x$nlpars[[nlp]], ...))
  }
  out
}

#' @export
exclude_pars_internal.btl <- function(x, data, ...) {
  p <- usc(combine_prefix(x))
  out <- c(
    paste0("temp", p, "_Intercept"),
    paste0(c("hs_local", "hs_global", "zb"), p)
  )
  sms <- get_sm_labels(x, data)
  if (length(sms) && !is.null(data)) {
    for (i in seq_along(sms)) {
      nb <- seq_len(attr(sms, "nbases")[[i]])
      out <- c(out, paste0("zs", p, "_", i, "_", nb))
    } 
  }
  out
}
