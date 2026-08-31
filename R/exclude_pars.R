# list parameters NOT to be saved by Stan
# @return a vector of parameter names to be excluded
exclude_pars <- function(x, ...) {
  UseMethod("exclude_pars")
}

#' @export
exclude_pars.default <- function(x, ...) {
  character(0)
}

#' @export
exclude_pars.brmsfit <- function(x, bframe = NULL, ...) {
  out <- character(0)
  if (is.null(bframe)) {
    # needed for the moment until brmsframe is stored in brmsfit
    bframe <- brmsframe(x$formula, data = x$data)
  }
  stopifnot(is.anybrmsframe(bframe))
  c(out) <- exclude_pars(bframe, save_pars = x$save_pars, ...)
  c(out) <- exclude_pars_re(bframe, save_pars = x$save_pars, ...)
  c(out) <- exclude_pars_me(bframe, save_pars = x$save_pars, ...)
  out <- unique(out)
  out <- setdiff(out, x$save_pars$manual)
  out
}

#' @export
exclude_pars.mvbrmsframe <- function(x, save_pars, ...) {
  out <- c("Rescor", "Sigma")
  if (!save_pars$all) {
    c(out) <- c("Lrescor", "LSigma")
  }
  for (i in seq_along(x$terms)) {
    c(out) <- exclude_pars(x$terms[[i]], save_pars = save_pars, ...)
  }
  out
}

#' @export
exclude_pars.brmsframe <- function(x, save_pars, ...) {
  resp <- usc(combine_prefix(x))
  par_classes <- c("Lncor", "Cortime")
  out <- paste0(par_classes, resp)
  if (!save_pars$all) {
    par_classes <- c(
      "ordered_Intercept", "fixed_Intercept",
      "theta", "Llncor", "Lcortime"
    )
    c(out) <- paste0(par_classes, resp)
  }
  for (dp in names(x$dpars)) {
    c(out) <- exclude_pars(x$dpars[[dp]], save_pars = save_pars, ...)
  }
  for (nlp in names(x$nlpars)) {
    c(out) <- exclude_pars(x$nlpars[[nlp]], save_pars = save_pars, ...)
  }
  if (is.formula(x$adforms$mi)) {
    if (!(isTRUE(save_pars$latent) || x$resp %in% save_pars$latent)) {
      c(out) <- paste0("Yl", resp)
    }
  }
  if (!(isTRUE(save_pars$group) || ".err" %in% save_pars$group)) {
    # latent residuals are like group-level effects
    c(out) <- paste0("err", resp)
  }
  out
}

#' @export
exclude_pars.bframel <- function(x, save_pars, ...) {
  out <- character(0)
  p <- usc(combine_prefix(x))
  c(out) <- paste0("chol_cor", p)
  if (!save_pars$all) {
    # removed the "Intercept" and "first_Intercept" parameters from this list
    # to reduce the number of models that need refitting for moment matching
    par_classes <- c(
      "bQ", "zb", "zbsp", "zbs", "zar", "zma", "hs_local", "R2D2_phi",
      "scales", "merged_Intercept", "finite_Intercept", "zcar", "nszcar",
      "zerr"
    )
    c(out) <- paste0(par_classes, p)
    if (has_re_s2z_conventional(x)) {
      c(out) <- paste0(
        c(
          "theta_s2z", "z_theta_s2z", "theta_s2z_active", "fixed_s2z",
          "zfixed_s2z", "sdfixed_s2z"
        ),
        p
      )
      infos_s2z <- re_s2z_infos(x)
      if (length(infos_s2z)) {
        info_s2z <- infos_s2z[[1L]]
        inactive_q_s2z <- info_s2z$inactive_q
        if (isTRUE(info_s2z$ordinal)) {
          inactive_q_s2z <- intersect(inactive_q_s2z, info_s2z$slope_q)
          threshold <- info_s2z$threshold
          for (i in unique(threshold$group_index)) {
            row <- threshold[match(i, threshold$group_index), , drop = FALSE]
            threshold_class <- if (identical(row$kind, "first")) {
              "first_theta_Intercept"
            } else {
              "theta_Intercept"
            }
            threshold_group <- str_if(has_thres_groups(x), usc(i))
            c(out) <- paste0(threshold_class, p, threshold_group)
          }
        }
        n_inactive_s2z <- length(inactive_q_s2z)
        if (n_inactive_s2z) {
          c(out) <- paste0(
            "par_fixed_s2z", p, "_", seq_len(n_inactive_s2z)
          )
        }
        q <- length(info_s2z$qnames)
      } else {
        q <- length(x$frame$fe$vars)
      }
      c(out) <- paste0("udf_b_s2z", p, "_", seq_len(q))
    }
    smframe <- x$frame$sm
    for (i in seq_rows(smframe)) {
      nb <- seq_len(smframe$nbases[i])
      c(out) <- paste0("zs", p, "_", i, "_", nb)
    }
  }
  out
}

# exclude variables related to random effects
exclude_pars_re <- function(bframe, save_pars, ...) {
  reframe <- bframe$frame$re
  stopifnot(is.reframe(reframe))
  out <- list()
  if (!has_rows(reframe)) {
    return(out)
  }
  rm_re_pars <- c(if (!save_pars$all) c("z", "L"), "Cor", "r")
  for (id in unique(reframe$id)) {
    c(out) <- paste0(rm_re_pars, "_", id)
    r <- subset2(reframe, id = id)
    if (!save_pars$all && !isTRUE(r$s2z[1]) &&
        !identical(re_s2z_center_mode(r), "noncentered")) {
      c(out) <- paste0(c("L_center_re_", "log_jacobian_re_"), id)
      if (identical(re_s2z_center_mode(r), "auto")) {
        c(out) <- paste0(c("rho_s2z_", "mean_rho_s2z_"), id)
      }
    }
    if (!save_pars$all && isTRUE(r$s2z[1])) {
      s2z_classes <- c(
        "z_s2z", "z_mean_s2z", "r_s2z", "H_s2z",
        "q_explicit_s2z", "prior_mean_s2z",
        "prior_prec_s2z", "prior_scale_s2z",
        "W_matheron_s2z", "sqrt_W_matheron_s2z",
        "L_W_matheron_s2z", "theta_white_matheron_s2z",
        "group_scale_s2z", "group_prec_s2z",
        "L_Sigma_s2z", "L_mean_s2z", "Q_Sigma_s2z",
        "P_group_s2z", "group_info_s2z", "h_group_s2z",
        "P_s2z", "L_P_s2z", "H_joint_s2z", "h_joint_s2z",
        "D_s2z", "sqrt_D_s2z", "D_diag_s2z", "intercept_map_s2z",
        "rank1_info_s2z", "group_quad_s2z", "joint_quad_s2z",
        "log_det_partial_s2z",
        "mhat_s2z", "qhat_s2z", "white_s2z", "mean_r_s2z",
        "q_recovered_s2z", "z_sd_s2z", "z_sd_mean_s2z",
        "relative_sd_s2z", "reference_sd_s2z", "sd_level_s2z"
      )
      c(out) <- paste0(s2z_classes, "_", id)
      rp <- usc(combine_prefix(r))
      c(out) <- paste0("r_s2z_", r$id, rp, "_", r$cn)
      if (identical(re_s2z_center_mode(r), "auto")) {
        c(out) <- paste0(c("rho_s2z_", "mean_rho_s2z_"), id)
        if (all(re_s2z_latent(r))) {
          fisher_info <- stan_re_s2z_latent_fisher_info(
            id, r = r, bframe = bframe, threads = NULL
          )
          dependency_names <- vapply(
            fisher_info$dependency_info,
            function(x) x$vector_name,
            character(1)
          )
          c(out) <- dependency_names
        }
      }
    }
  }
  if (isFALSE(save_pars$group)) {
    p <- usc(combine_prefix(reframe))
    c(out) <- paste0("r_", reframe$id, p, "_", reframe$cn)
    varying_ids <- unique(reframe$id[reframe$scale == "varying"])
    c(out) <- paste0("sd_level_", varying_ids)
  } else if (is.character(save_pars$group)) {
    sub_reframe <- reframe[!reframe$group %in% save_pars$group, ]
    if (has_rows(sub_reframe)) {
      sub_p <- usc(combine_prefix(sub_reframe))
      c(out) <- paste0("r_", sub_reframe$id, sub_p, "_", sub_reframe$cn)
      varying_ids <- unique(sub_reframe$id[sub_reframe$scale == "varying"])
      c(out) <- paste0("sd_level_", varying_ids)
    }
  }
  reframe_t <- subset_reframe_dist(reframe, "student")
  if (!save_pars$all && has_rows(reframe_t)) {
    c(out) <- paste0(c("udf_", "dfm_"), reframe_t$ggn)
  }
  out
}

# exclude variables related to noise-free variables
exclude_pars_me <- function(bframe, save_pars, ...) {
  meframe <- bframe$frame$me
  stopifnot(is.meframe(meframe))
  out <- list()
  if (!has_rows(meframe)) {
    return(out)
  }
  I <- seq_along(unique(meframe$grname))
  K <- seq_rows(meframe)
  c(out) <- paste0(c("Corme_"), I)
  if (!save_pars$all) {
    c(out) <- c(paste0("zme_", K), paste0("Lme_", I))
  }
  if (isFALSE(save_pars$latent)) {
    c(out) <- paste0("Xme_", K)
  } else if (is.character(save_pars$latent)) {
    sub_K <- K[!meframe$xname %in% save_pars$latent]
    if (length(sub_K)) {
      c(out) <- paste0("Xme_", sub_K)
    }
  }
  out
}

#' Control Saving of Parameter Draws
#'
#' Control which (draws of) parameters should be saved in a \pkg{brms}
#' model. The output of this function is meant for usage in the
#' \code{save_pars} argument of \code{\link{brm}}.
#'
#' @param group A flag to indicate if group-level coefficients for
#'   each level of the grouping factors should be saved (default is
#'   \code{TRUE}). Set to \code{FALSE} to save memory. Alternatively,
#'   \code{group} may also be a character vector naming the grouping factors
#'   for which to save draws of coefficients.
#' @param latent A flag to indicate if draws of latent variables obtained by
#'   using \code{me} and \code{mi} terms should be saved (default is
#'   \code{FALSE}). Saving these draws allows to better use methods such as
#'   \code{posterior_predict} with the latent variables but leads to very large
#'   \R objects even for models of moderate size and complexity. Alternatively,
#'   \code{latent} may also be a character vector naming the latent variables
#'   for which to save draws.
#' @param all A flag to indicate if draws of all variables defined in Stan's
#'   \code{parameters} block should be saved (default is \code{FALSE}). Saving
#'   these draws is required in order to apply the certain methods such as
#'   \code{bridge_sampler} and \code{bayes_factor}.
#' @param manual A character vector naming Stan variable names which should be
#'   saved. These names should match the variable names inside the Stan code
#'   before renaming. This feature is meant for power users only and will rarely
#'   be useful outside of very special cases.
#'
#' @return A list of class \code{"save_pars"}.
#'
#' @examples
#' \dontrun{
#' # don't store group-level coefficients
#' fit <- brm(count ~ zAge + zBase * Trt + (1|patient),
#'            data = epilepsy, family = poisson(),
#'            save_pars = save_pars(group = FALSE))
#' variables(fit)
#' }
#'
#' @export
save_pars <- function(group = TRUE, latent = FALSE, all = FALSE,
                      manual = NULL) {
  out <- list()
  if (is.logical(group)) {
    out$group <- as_one_logical(group)
  } else {
    out$group <- as.character(group)
  }
  if (is.logical(latent)) {
    out$latent <- as_one_logical(latent)
  } else {
    out$latent <- as.character(latent)
  }
  out$all <- as_one_logical(all)
  out$manual <- as.character(manual)
  class(out) <- "save_pars"
  out
}

# validate 'save_pars' argument
# deprecated arguments:
# @param save_ranef save varying effects per level?
# @param save_mevars save noise-free variables?
# @param save_all_pars save all variables from the 'parameters' block?
# @return validated 'save_pars' argument
validate_save_pars <- function(save_pars, save_ranef = NULL, save_mevars = NULL,
                               save_all_pars = NULL) {
  if (is.null(save_pars)) {
    save_pars <- save_pars()
  }
  if (!is.save_pars(save_pars)) {
    stop2("Argument 'save_pars' needed to be created via 'save_pars()'.")
  }
  if (!is.null(save_ranef)) {
    warning2(
      "Argument 'save_ranef' is deprecated. Please use argument ",
      "'group' in function 'save_pars()' instead."
    )
    save_pars$group <- as_one_logical(save_ranef)
  }
  if (!is.null(save_mevars)) {
    warning2(
      "Argument 'save_mevars' is deprecated. Please use argument ",
      "'latent' in function 'save_pars()' instead."
    )
    save_pars$latent <- as_one_logical(save_mevars)
  }
  if (!is.null(save_all_pars)) {
    warning2(
      "Argument 'save_all_pars' is deprecated. Please use argument ",
      "'all' in function 'save_pars()' instead."
    )
    save_pars$all <- as_one_logical(save_all_pars)
  }
  save_pars
}

is.save_pars <- function(x) {
  inherits(x, "save_pars")
}
