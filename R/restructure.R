#' Retructure Old \code{brmsfit} Objects
#'
#' Restructure old \code{brmsfit} objects to work with 
#' the latest \pkg{brms} version. This function is called
#' internally when applying post-processing methods.
#' However, in order to avoid uncessary run time caused
#' by the restructuring, I recommend explicitely calling
#' \code{restructure} once per model after updating \pkg{brms}.
#' 
#' @param x An object of class \code{brmsfit}.
#' @param rstr_summary Logical; If \code{TRUE}, the cached summary
#'   stored by \pkg{rstan} is restructured as well.
#'   
#' @return A \code{brmsfit} object compatible with the latest version
#'   of \pkg{brms}.
#'   
#' @export
restructure <- function(x, rstr_summary = FALSE) {
  stopifnot(is.brmsfit(x))
  if (isTRUE(attr(x, "restructured"))) {
    return(x)  # already restructured
  }
  if (is.null(x$version)) {
    # this is the latest version without saving the version number
    x$version <- list(brms = package_version("0.9.1"))
  } else if (is.package_version(x$version)) {
    # also added the rstan version in brms 1.5.0
    x$version <- list(brms = x$version)
  }
  version <- x$version$brms
  if (version < utils::packageVersion("brms")) {
    # slot 'threshold' deprecated as of brms > 1.7.0
    if (is_ordinal(x$family) && is.null(x$family$threshold)) {
      x$family$threshold <- x$threshold
    }
    x$threshold <- NULL
    # slot 'nonlinear' deprecated as of brms > 0.9.1
    x$formula <- SW(amend_formula(
      formula(x), data = model.frame(x), 
      family = family(x), autocor = x$autocor,
      nonlinear = x$nonlinear
    ))
    x$nonlinear <- x$partial <- NULL
    x$formula[["old_mv"]] <- is_old_mv(x)
    bterms <- parse_bf(formula(x))
    if (!isTRUE(x$formula[["old_mv"]])) {
      x$data <- rm_attr(x$data, "brmsframe")
      x$data <- update_data(x$data, bterms) 
    }
    x$ranef <- tidy_ranef(bterms, model.frame(x))
    if ("prior_frame" %in% class(x$prior)) {
      class(x$prior) <- c("brmsprior", "data.frame") 
    }
    if (is(x$autocor, "cov_fixed")) {
      # deprecated as of brms 1.4.0
      class(x$autocor) <- "cor_fixed"
    }
    if (version <= "0.9.1") {
      # update gaussian("log") to lognormal() family
      nresp <- length(bterms$response)
      if (is_old_lognormal(x$family, nresp = nresp, version = version)) {
        x$family <- x$formula$family <- lognormal()
      }
    }
    if (version <= "0.10.0.9000") {
      if (length(bterms$dpars$mu$nlpars)) {
        # nlpar and group have changed positions
        change <- change_old_re(x$ranef, pars = parnames(x),
                                dims = x$fit@sim$dims_oi)
        x <- do_renaming_old(x, change)
      }
    }
    if (version < "1.0.0") {
      # double underscores were added to group-level parameters
      change <- change_old_re2(x$ranef, pars = parnames(x),
                               dims = x$fit@sim$dims_oi)
      x <- do_renaming_old(x, change)
    }
    if (version <= "1.0.1") {
      # names of spline parameters had to be changed after
      # allowing for multiple covariates in one spline term
      change <- change_old_sm(bterms, pars = parnames(x),
                              dims = x$fit@sim$dims_oi)
      x <- do_renaming_old(x, change)
    }
    if (version <= "1.2.0") {
      x$ranef$type[x$ranef$type == "mono"] <- "mo"
      x$ranef$type[x$ranef$type == "cse"] <- "cs"
    }
    if (version <= "1.8.0") {
      att <- attributes(x$exclude)
      if (is.null(att$save_ranef)) {
        attr(x$exclude, "save_ranef") <- 
          any(grepl("^r_", parnames(x))) || !nrow(x$ranef)
      }
      if (is.null(att$save_mevars)) {
        attr(x$exclude, "save_mevars") <- 
          any(grepl("^Xme_", parnames(x)))
      }
    }
    if (version <= "1.8.0.1") {
      x$prior[, c("resp", "dpar")] <- ""
    }
    if (version <= "1.9.0.3") {
      # names of monotonic parameters had to be changed after
      # allowing for interactions in monotonic terms
      change <- change_old_mo(bterms, x$data, pars = parnames(x))
      x <- do_renaming_old(x, change)
    }
    stan_env <- attributes(x$fit)$.MISC
    if (rstr_summary && exists("summary", stan_env)) {
      stan_summary <- get("summary", stan_env)
      old_parnames <- rownames(stan_summary$msd)
      if (!identical(old_parnames, parnames(x))) {
        # do not rename parameters in the cached summary
        # just let rstan compute the summary again
        remove("summary", pos = stan_env)
      }
    }
  }
  structure(x, restructured = TRUE)
}

change_old_re <- function(ranef, pars, dims) {
  # interchanges group and nlpar in names of group-level parameters
  # for brms <= 0.10.0.9000 only
  # Args:
  #   ranef: output of tidy_ranef
  #   pars: names of all parameters in the model
  #   dims: dimension of parameters
  # Returns:
  #   a list whose elements can be interpreted by do_renaming_old
  change <- list()
  for (id in unique(ranef$id)) {
    r <- subset2(ranef, id = id)
    g <- r$group[1]
    nlpar <- r$nlpar[1]
    stopifnot(nzchar(nlpar))
    # rename sd-parameters
    old_sd_names <- paste0("sd_", nlpar, "_", g, "_", r$coef)
    new_sd_names <- paste0("sd_", g, "_", nlpar, "_", r$coef)
    for (i in seq_along(old_sd_names)) {
      change <- lc(change, 
        change_simple(old_sd_names[i], new_sd_names[i], pars, dims)
      )
    }
    # rename cor-parameters
    new_cor_names <- get_cornames(
      paste0(nlpar, "_", r$coef), type = paste0("cor_", g),
      brackets = FALSE, sep = "_"
    )
    old_cor_names <- get_cornames(
      r$coef, brackets = FALSE, sep = "_",
      type = paste0("cor_", nlpar, "_", g)
    )
    for (i in seq_along(old_cor_names)) {
      change <- lc(change, 
        change_simple(old_cor_names[i], new_cor_names[i], pars, dims)
      )
    } 
    # rename r-parameters
    old_r_name <- paste0("r_", nlpar, "_", g)
    new_r_name <- paste0("r_", g, "_", nlpar)
    levels <- gsub("[ \t\r\n]", ".", attr(ranef, "levels")[[g]])
    index_names <- make_index_names(levels, r$coef, dim = 2)
    new_r_names <- paste0(new_r_name, index_names)
    change <- lc(change, 
      change_simple(
        old_r_name, new_r_names, pars, dims, 
        pnames = new_r_name
      )
    )
  }
  change
}

change_old_re2 <- function(ranef, pars, dims) {
  # add double underscore in group-level parameters
  # for brms < 1.0.0 only
  # assumes that group and nlpar are correctly ordered already
  # Args:
  #   ranef: output of tidy_ranef
  #   pars: names of all parameters in the model
  #   dims: dimension of parameters
  # Returns:
  #   a list whose elements can be interpreted by do_renaming_old
  change <- list()
  for (id in unique(ranef$id)) {
    r <- subset2(ranef, id = id)
    g <- r$group[1]
    nlpars_usc <- usc(r$nlpar, "suffix")
    # rename sd-parameters
    old_sd_names <- paste0("sd_", g, "_", nlpars_usc, r$coef)
    new_sd_names <- paste0("sd_", g, "__", nlpars_usc, r$coef)
    for (i in seq_along(old_sd_names)) {
      change <- lc(change, 
        change_simple(old_sd_names[i], new_sd_names[i], pars, dims)
      )
    }
    # rename cor-parameters
    new_cor_names <- get_cornames(
      paste0(nlpars_usc, r$coef), type = paste0("cor_", g),
      brackets = FALSE
    )
    old_cor_names <- get_cornames(
      paste0(nlpars_usc, r$coef), type = paste0("cor_", g),
      brackets = FALSE, sep = "_"
    )
    for (i in seq_along(old_cor_names)) {
      change <- lc(change, 
        change_simple(old_cor_names[i], new_cor_names[i], pars, dims))
    } 
    # rename r-parameters
    for (nlpar in unique(r$nlpar)) {
      sub_r <- r[r$nlpar == nlpar, ]
      old_r_name <- paste0("r_", g, usc(nlpar))
      new_r_name <- paste0("r_", g, usc(usc(nlpar)))
      levels <- gsub("[ \t\r\n]", ".", attr(ranef, "levels")[[g]])
      index_names <- make_index_names(levels, sub_r$coef, dim = 2)
      new_r_names <- paste0(new_r_name, index_names)
      change <- lc(change, 
        change_simple(
          old_r_name, new_r_names, pars, dims, 
          pnames = new_r_name
        )
      )
    }
  }
  change
}

change_old_sm <- function(bterms, pars, dims) {
  # change names of spline parameters fitted with brms <= 1.0.1
  # this became necessary after allowing smooths with multiple covariates
  stopifnot(is.brmsterms(bterms))
  .change_old_sm <- function(bt) {
    change <- list()
    sm_labels <- get_sm_labels(bt)
    if (length(sm_labels)) {
      p <- usc(combine_prefix(bt), "suffix")
      old_smooths <- rename(paste0(p, sm_labels))
      new_smooths <- rename(paste0(p, get_sm_labels(bt, covars = TRUE)))
      old_sds_pars <- paste0("sds_", old_smooths)
      new_sds_pars <- paste0("sds_", new_smooths, "_1")
      old_s_pars <- paste0("s_", old_smooths)
      new_s_pars <- paste0("s_", new_smooths, "_1")
      for (i in seq_along(old_smooths)) {
        change <- lc(change,
          change_simple(old_sds_pars[i], new_sds_pars[i], pars, dims)
        )
        indices <- seq_len(dims[[old_s_pars[i]]])
        new_s_par_indices <- paste0(new_s_pars[i], "[", indices, "]")
        change <- lc(change,
          change_simple(
            old_s_pars[i], new_s_par_indices, pars, dims,
            pnames = new_s_pars[i]
          )
        )
      }
    }
    return(change)
  }
  
  change <- list()
  if (length(bterms$response) > 1L) {
    for (r in bterms$response) {
      bterms$dpars$mu$resp <- r
      change <- c(change, .change_old_sm(bterms$dpars$mu))
    }
    bterms$dpars$mu <- NULL
  }
  for (dp in names(bterms$dpars)) {
    bt <- bterms$dpars[[dp]]
    if (length(bt$nlpars)) {
      for (nlp in names(bt$nlpars)) {
        change <- c(change, .change_old_sm(bt$nlpars[[nlp]]))
      }
    } else {
      change <- c(change, .change_old_sm(bt))
    }
  }
  change
}

change_old_mo <- function(bterms, data, pars) {
  # change names of monotonic effects fitted with brms <= 1.9.0
  # this became necessary after implementing monotonic interactions
  stopifnot(is.brmsterms(bterms))
  .change_old_mo <- function(bt) {
    change <- list()
    monef <- get_mo_labels(bt, data)
    if (!length(monef)) {
      return(change)
    }
    p <- usc(combine_prefix(bt))
    bmo_prefix <- paste0("bmo", p, "_")
    bmo_regex <- paste0("^", bmo_prefix, "[^_]+$")
    bmo_old <- pars[grepl(bmo_regex, pars)]
    bmo_new <- paste0(bmo_prefix, rename(monef))
    if (length(bmo_old) != length(bmo_new)) {
      stop2("Restructuring failed. Please refit your ", 
            "model with the latest version of brms.")
    }
    for (i in seq_along(bmo_old)) {
      pos <- grepl(paste0("^", bmo_old[i]), pars)
      change <- lc(change, nlist(pos, fnames = bmo_new[i]))
    }
    simo_regex <- paste0("^simplex", p, "_[^_]+$")
    simo_old_all <- pars[grepl(simo_regex, pars)]
    simo_index <- get_matches("\\[[[:digit:]]+\\]$", simo_old_all)
    simo_old <- unique(sub("\\[[[:digit:]]+\\]$", "", simo_old_all))
    simo_coef <- get_simo_labels(monef)
    for (i in seq_along(simo_old)) {
      regex_pos <- paste0("^", simo_old[i])
      pos <- grepl(regex_pos, pars)
      simo_new <- paste0("simo", p, "_", simo_coef[i])
      simo_index_part <- simo_index[grepl(regex_pos, simo_old_all)]
      simo_new <- paste0(simo_new, simo_index_part)
      change <- lc(change, nlist(pos, fnames = simo_new))
    }
    return(change)
  }
  
  change <- list()
  if (length(bterms$response) > 1L) {
    for (r in bterms$response) {
      bterms$dpars$mu$resp <- r
      change <- c(change, .change_old_mo(bterms$dpars$mu))
    }
    bterms$dpars$mu <- NULL
  }
  for (dp in names(bterms$dpars)) {
    bt <- bterms$dpars[[dp]]
    if (length(bt$nlpars)) {
      for (nlp in names(bt$nlpars)) {
        change <- c(change, .change_old_mo(bt$nlpars[[nlp]]))
      }
    } else {
      change <- c(change, .change_old_mo(bt))
    }
  }
  change
}

change_simple <- function(oldname, fnames, pars, dims,
                          pnames = fnames) {
  # helper function for very simple renaming
  # only used in renaming of old models
  pos <- grepl(paste0("^", oldname), pars)
  if (any(pos)) {
    out <- nlist(pos, oldname, pnames, fnames,
                 dims = dims[[oldname]])
  } else {
    out <- NULL
  }
  out
}

do_renaming_old <- function(x, change) {
  # perform actual renaming of Stan parameters
  # This is an old version of do_renaming allowing
  # for more complex changes in the stanfit object,
  # which are no longer needed for the renaming of new models
  # Args:
  #   change: a list of lists each element allowing
  #           to rename certain parameters
  #   x: An object of class brmsfit
  # Returns:
  #   A brmsfit object with updated parameter names
  .do_renaming_old <- function(x, change) {
    chains <- length(x$fit@sim$samples) 
    x$fit@sim$fnames_oi[change$pos] <- change$fnames
    for (i in seq_len(chains)) {
      names(x$fit@sim$samples[[i]])[change$pos] <- change$fnames
      if (!is.null(change$sort)) {
        x$fit@sim$samples[[i]][change$pos] <- 
          x$fit@sim$samples[[i]][change$pos][change$sort]
      }
    }
    if (!is.null(change$oldname)) {
      # Changing dims_oi interferes with rstan::unconstrain_pars(), 
      # which is used in the bridgesampling package. The code below is
      # only retained for backwards compatibility with older models.
      onp <- match(change$oldname, names(x$fit@sim$dims_oi))
      if (is.null(onp) || is.na(onp)) {
        warning2(
          "Parameter ", change$oldname, " could not be renamed. ",
          "This should not happen. \nPlease inform me so that ",
          "I can fix this problem."
        )
      } else {
        if (is.null(change$pnames)) {
          # only needed to collapse multiple r_<i> of the same grouping factor
          x$fit@sim$dims_oi[[onp]] <- NULL
        } else {
          # rename dims_oi to match names in fnames_oi
          dims <- x$fit@sim$dims_oi
          x$fit@sim$dims_oi <- c(
            if (onp > 1) dims[1:(onp - 1)],
            make_dims(change),
            dims[(onp + 1):length(dims)]
          )
        }
      }
      x$fit@sim$pars_oi <- names(x$fit@sim$dims_oi)
    }
    return(x)
  }
  for (i in seq_along(change)) {
    x <- .do_renaming_old(x, change[[i]])
  }
  x
}
