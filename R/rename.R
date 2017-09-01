rename_pars <- function(x) {
  # rename parameters (and possibly change their dimensions) 
  # within the stanfit object to ensure reasonable parameter names
  # Args:
  #   x: a brmsfit obejct
  # Returns:
  #   a brmfit object with adjusted parameter names and dimensions
  if (!length(x$fit@sim)) {
    return(x) 
  }
  x <- reorder_pars(x)
  bterms <- with(x, parse_bf(formula, family = family, autocor = autocor))
  family <- family(x)
  data <- model.frame(x)
  pars <- parnames(x)
  
  # find positions of parameters and define new names
  change <- list()
  resp <- bterms$response
  if (length(resp) > 1L) {
    # rename effects in multivariate models
    for (r in resp) {
      bterms$dpars[["mu"]]$resp <- r
      change_eff <- change_effects(
        bterms$dpars[["mu"]], data = data, 
        pars = pars, stancode = stancode(x)
      )
      change <- c(change, change_eff)
    }
    bterms$dpars[["mu"]] <- NULL
  }
  # rename effects of distributional parameters
  for (ap in names(bterms$dpars)) {
    change_eff <- change_effects(
      bterms$dpars[[ap]], data = data, 
      pars = pars, stancode = stancode(x)
    )
    change <- c(change, change_eff)
  }
  change <- c(change,
    change_re(x$ranef, pars = pars),
    change_Xme(bterms, pars = pars),
    change_autocor(bterms, data = data, pars = pars)
  )
  # rename residual parameters of multivariate linear models
  if (is_linear(family) && length(bterms$response) > 1L) {
    corfnames <- paste0("sigma_", bterms$response)
    change <- lc(change, 
      list(
        pos = grepl("^sigma\\[", pars), oldname = "sigma",
        pnames = corfnames, fnames = corfnames
      )
    )
    change <- c(change, 
      change_prior(
        class = "sigma", pars = pars, 
        names = bterms$response
      )
    )
    rescor_names <- get_cornames(
      bterms$response, type = "rescor", brackets = FALSE
    )
    change <- lc(change, 
      list(
        pos = grepl("^rescor\\[", pars), 
        oldname = "rescor", pnames = rescor_names,
        fnames = rescor_names
      )
    )
  }
  # perform the actual renaming in x$fit@sim
  x <- do_renaming(x, change)
  x <- compute_quantities(x)
  x
}

change_effects <- function(x, ...) {
  # helps in renaming parameters after model fitting
  UseMethod("change_effects")
}

#' @export
change_effects.btl <- function(x, data, pars, stancode = "", ...) {
  # helps in renaming various kinds of effects
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  c(
    change_fe(x, data, pars, stancode = stancode),
    change_sm(x, data, pars),
    change_cs(x, data, pars),
    change_mo(x, data, pars),
    change_me(x, data, pars),
    change_gp(x, data, pars)
  )
}

change_effects.btnl <- function(x, data, pars, ...) {
  # helps in renaming effects for non-linear parameters
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  change <- list()
  for (nlp in names(x$nlpars)) {
    change <- c(change, 
      change_effects(x$nlpars[[nlp]], data, pars, ...)
    )
  }
  change
}

change_fe <- function(bterms, data, pars, stancode = "") {
  # helps in renaming fixed effects parameters
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  change <- list()
  px <- check_prefix(bterms)
  fixef <- colnames(data_fe(bterms, data)$X)
  fixef <- rm_int_fe(fixef, stancode, px = px)
  if (length(fixef)) {
    b <- paste0("b", usc(combine_prefix(px), "prefix"))
    pos <- grepl(paste0("^", b, "\\["), pars)
    bnames <- paste0(b, "_", fixef)
    change <- lc(change, list(pos = pos, fnames = bnames))
    change <- c(change,
      change_prior(class = b, pars = pars, names = fixef)
    )
  }
  change
}

change_mo <- function(bterms, data, pars) {
  # helps in renaming monotonous effects parameters
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  change <- list()
  monef <- get_mo_labels(bterms, data)
  if (!length(monef)) {
    return(change) 
  }
  p <- usc(combine_prefix(bterms))
  simo_coef <- get_simo_labels(monef)
  monef <- rename(monef)
  bmo <- paste0("bmo", p)
  pos <- grepl(paste0("^", bmo, "\\["), pars)
  newnames <- paste0("bmo", p, "_", monef)
  change <- lc(change, nlist(pos, fnames = newnames))
  change <- c(change, 
    change_prior(class = bmo, pars = pars, names = monef)
  )
  for (i in seq_along(simo_coef)) {
    simo_old <- paste0("simo", p, "_", i)
    simo_new <- paste0("simo", p, "_", simo_coef[i])
    pos <- grepl(paste0("^", simo_old, "\\["), pars)
    simo_names <- paste0(simo_new, "[", seq_len(sum(pos)), "]")
    change <- lc(change, list(pos = pos, fnames = simo_names))
    change <- c(change,
      change_prior(
        class = simo_old, new_class = simo_new,
        pars = pars, is_vector = TRUE
      )
    )
  }
  change
}

change_cs <- function(bterms, data, pars) {
  # helps in renaming category specific effects parameters
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  change <- list()
  csef <- colnames(data_cs(bterms, data)$Xcs)
  if (length(csef)) {
    ncse <- length(csef)
    thres <- sum(grepl("^b_Intercept\\[", pars))
    csenames <- t(outer(csef, paste0("[", 1:thres, "]"), FUN = paste0))
    csenames <- paste0("bcs_", csenames)
    sort_cse <- ulapply(seq_len(ncse), seq, to = thres * ncse, by = ncse)
    change <- lc(change, 
      list(pos = grepl("^bcs\\[", pars), fnames = csenames, sort = sort_cse)
    )
    change <- c(change, 
      change_prior(class = "bcs", pars = pars, names = csef)
    )
  }
  change
}

change_me <- function(bterms, data, pars) {
  # helps in renaming parameters of noisy variables
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  change <- list()
  meef <- get_me_labels(bterms, data)
  if (length(meef)) {
    p <- usc(combine_prefix(bterms), "prefix")
    meef <- rename(meef)
    # rename coefficients of noise free terms
    bme <- paste0("bme", p)
    pos <- grepl(paste0("^", bme, "\\["), pars)
    bmenames <- paste0(bme, "_", meef)
    change <- lc(change, nlist(pos, fnames = bmenames))
    change <- c(change,
      change_prior(class = bme, pars = pars, names = meef)
    )
  }
  change
}

change_Xme <- function(bterms, pars) {
  # helps in renaming global noise free variables
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  stopifnot(is.brmsterms(bterms))
  change <- list()
  if (any(grepl("^Xme_", pars))) {
    uni_me <- get_uni_me(bterms)
    for (i in seq_along(uni_me)) {
      Xme <- paste0("Xme_", i)
      pos <- grepl(paste0("^", Xme, "\\["), pars)
      Xme_new <- paste0("Xme_", rename(uni_me[i]))
      fnames <- paste0(Xme_new, "[", seq_len(sum(pos)), "]")
      change <- lc(change, nlist(pos, fnames))
    }
  }
  change
}

change_gp <- function(bterms, data, pars) {
  # helps in renaming parameters of gaussian processes
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  change <- list()
  gpef <- get_gp_labels(bterms, data = data, covars = TRUE)
  p <- usc(combine_prefix(bterms), "prefix")
  for (i in seq_along(gpef)) {
    # rename GP hyperparameters
    by_levels = attr(gpef, "by_levels")[[i]]
    gp_names <- paste0(gpef[i], usc(by_levels))
    sdgp <- paste0("sdgp", p)
    sdgp_old <- paste0(sdgp, "_", i)
    sdgp_pos <- grepl(paste0("^", sdgp_old, "\\["), pars)
    sdgp_names <- paste0(sdgp, "_", gp_names)
    lscale <- paste0("lscale", p)
    lscale_old <- paste0(lscale, "_", i)
    lscale_pos <- grepl(paste0("^", lscale_old, "\\["), pars)
    lscale_names <- paste0(lscale, "_", gp_names)
    change <- lc(change, 
      nlist(pos = sdgp_pos, fnames = sdgp_names),
      nlist(pos = lscale_pos, fnames = lscale_names)
    )
    change <- c(change,
      change_prior(class = sdgp_old, pars = pars, 
                   names = gpef[i], new_class = sdgp),
      change_prior(class = lscale_old, pars = pars, 
                   names = gpef[i], new_class = lscale)
    )
    zgp <- paste0("zgp", p)
    zgp_old <- paste0(zgp, "_", i)
    zgp_pos <- grepl(paste0("^", zgp_old, "\\["), pars)
    if (any(zgp_pos)) {
      # users may choose not to save zgp
      zgp_new <- paste0(zgp, "_", gpef[i])
      fnames <- paste0(zgp_new, "[", seq_len(sum(zgp_pos)), "]")
      change_zgp <- nlist(pos = zgp_pos, fnames = fnames)
      change <- lc(change, change_zgp)
    }
  }
  change
}

change_sm <- function(bterms, data, pars) {
  # helps in renaming smoothing term parameters
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  change <- list()
  smooths <- get_sm_labels(bterms, data, covars = TRUE)
  if (length(smooths)) {
    stopifnot(!is.null(attr(smooths, "nbases")))
    p <- usc(combine_prefix(bterms), "prefix")
    sds <- paste0("sds", p)
    sds_names <- paste0(sds, "_", smooths)
    s <- paste0("s", p)
    s_names <- paste0(s, "_", smooths)
    for (i in seq_along(smooths)) {
      nb <- attr(smooths, "nbases")[[i]]
      for (j in seq_len(nb)) {
        ij <- paste0(i, "_", j)
        sds_pos <- grepl(paste0("^", sds, "_", ij), pars)
        change <- lc(change, 
          list(pos = sds_pos, fnames = paste0(sds_names[i], "_", j))
        )
        s_pos <- grepl(paste0("^", s, "_", ij), pars)
        s_fnames <- paste0(s_names[i], "_", j, "[", seq_len(sum(s_pos)), "]")
        change <- lc(change, list(pos = s_pos, fnames = s_fnames))
        new_prior_class <- paste0(sds, "_", smooths[i], "_", j)
        change <- c(change, 
          change_prior(
            class = paste0(sds, "_", ij), pars = pars, 
            new_class = new_prior_class
          )
        ) 
      }
    }
  }
  change
}
  

change_re <- function(ranef, pars) {
  # helps in renaming group-level parameters
  # Args:
  #   ranef: list returned by tidy_ranef
  #   pars: names of all model parameters
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  change <- list()
  if (nrow(ranef)) {
    for (id in unique(ranef$id)) {
      r <- subset2(ranef, id = id)
      g <- r$group[1]
      suffix <- paste0(usc(combine_prefix(r), "suffix"), r$coef)
      rfnames <- paste0("sd_", g, "__", suffix)
      rpos <- grepl(paste0("^sd_", id, "(\\[|$)"), pars)
      change <- lc(change, list(pos = rpos, fnames = rfnames))
      change <- c(change,
        change_prior(
          class = paste0("sd_", id), pars = pars,
          new_class = paste0("sd_", g), 
          names = paste0("_", suffix)
        )
      )
      # rename group-level correlations
      if (nrow(r) > 1L && isTRUE(r$cor[1])) {
        type <- paste0("cor_", g)
        cor_names <- get_cornames(
          suffix, type = paste0("cor_", g), brackets = FALSE
        )
        cor_pos <- grepl(paste0("^cor_", id, "(\\[|$)"), pars)
        change <- lc(change, list(pos = cor_pos, fnames = cor_names))
        change <- c(change,
          change_prior(
            class = paste0("cor_", id), pars = pars,
            new_class = paste0("cor_", g)
          )
        )
      }
    }
    if (any(grepl("^r_", pars))) {
      change <- c(change, change_re_levels(ranef, pars = pars))
    }
  }
  change
} 

change_re_levels <- function(ranef, pars)  {
  # helps in renaming random effects 'r_.' parameters
  # Args:
  #   ranef: output of tidy_ranef
  #   pars: names of all model parameters
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  change <- list()
  for (i in seq_len(nrow(ranef))) {
    r <- ranef[i, ]
    p <- usc(combine_prefix(r))
    r_parnames <- paste0("r_", r$id, p, "_", r$cn)
    r_regex <- paste0("^", r_parnames, "(\\[|$)")
    change_rl <- list(pos = grepl(r_regex, pars))
    r_new_parname <- paste0("r_", r$group, usc(p))
    # rstan doesn't like whitespaces in parameter names
    levels <- gsub("[ \t\r\n]", ".", attr(ranef, "levels")[[r$group]])
    index_names <- make_index_names(levels, r$coef, dim = 2)
    change_rl$fnames <- paste0(r_new_parname, index_names)
    change <- lc(change, change_rl)
  }
  change
}

change_autocor <- function(bterms, data, pars) {
  # helps in renaming autocor parameters
  change <- list()
  if (is.cor_bsts(bterms$autocor)) {
    data <- order_data(data, bterms = bterms)
    if (nzchar(bterms$time$group)) {
      group <- gsub("[ \t\r\n]", "", get(bterms$time$group, data))
    } else {
      group <- rep(1, nrow(data)) 
    }
    if (nzchar(bterms$time$time)) {
      time <- gsub("[ \t\r\n]", "", get(bterms$time$time, data))
    } else {
      time <- ulapply(unique(group), function(g) seq_len(sum(group == g)))
    }
    loclev_pars <- paste0("loclev[", group, ",", time, "]")
    change <- lc(change, 
      list(pos = grepl("^loclev\\[", pars), fnames = loclev_pars)             
    )
  }
  change
}

change_prior <- function(class, pars, names = NULL, new_class = class,
                         is_vector = FALSE) {
  # helps in renaming prior parameters
  # Args: 
  #   class: the class of the parameters for which prior names should be changed
  #   pars: names of all parameters in the model
  #   names: names to replace digits at the end of parameter names
  #   new_class: replacement of the orginal class name
  #   is_vector: indicate if the prior parameter is a vector
  # Return:
  #   a list whose elements can be interpreted by do_renaming
  change <- list()
  regex <- paste0("^prior_", class, "(_[[:digit:]]+|$|\\[)")
  pos_priors <- which(grepl(regex, pars))
  if (length(pos_priors)) {
    priors <- gsub(
      paste0("^prior_", class), 
      paste0("prior_", new_class), 
      pars[pos_priors]
    )
    if (is_vector) {
      change <- lc(change, list(pos = pos_priors, fnames = priors))
    } else {
      digits <- sapply(priors, function(prior) {
        d <- regmatches(prior, gregexpr("_[[:digit:]]+$", prior))[[1]]
        if (length(d)) as.numeric(substr(d, 2, nchar(d))) else 0
      })
      for (i in seq_along(priors)) {
        if (digits[i] && !is.null(names)) {
          priors[i] <- gsub("[[:digit:]]+$", names[digits[i]], priors[i])
        }
        if (pars[pos_priors[i]] != priors[i]) {
          change <- lc(change, list(pos = pos_priors[i], fnames = priors[i]))
        }
      }
    }
  }
  change
}

change_old_re <- function(ranef, pars, dims) {
  # interchanges group and nlpar in names of group-level parameters
  # for brms <= 0.10.0.9000 only
  # Args:
  #   ranef: output of tidy_ranef
  #   pars: names of all parameters in the model
  #   dims: dimension of parameters
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
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
  #   a list whose elements can be interpreted by do_renaming
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

change_simple <- function(oldname, fnames, pars, dims,
                          pnames = fnames) {
  # helper function for very simple renaming
  pos <- grepl(paste0("^", oldname), pars)
  if (any(pos)) {
    out <- nlist(pos, oldname, pnames, fnames,
                 dims = dims[[oldname]])
  } else {
    out <- NULL
  }
  out
}

rm_int_fe <- function(fixef, stancode, px = list()) {
  # identifies if the intercept has to be removed from fixef
  # and returns adjusted fixef names
  p <- usc(combine_prefix(px), "suffix")
  int <- paste0("b_", p, "Intercept = temp_", p, "Intercept")
  loclev <- "vector[N] loclev;"
  if (any(ulapply(c(int, loclev), grepl, stancode, fixed = TRUE))) {
    fixef <- setdiff(fixef, "Intercept")
  } 
  fixef
}

make_index_names <- function(rownames, colnames = NULL, dim = 1) {
  # compute index names in square brackets for indexing stan parameters
  # Args:
  #   rownames: a vector of row names
  #   colnames: a vector of columns 
  #   dim: The number of dimensions of the output either 1 or 2
  # Returns:
  #   all index pairs of rows and cols
  if (!dim %in% c(1, 2))
    stop("dim must be 1 or 2")
  if (dim == 1) {
    index_names <- paste0("[", rownames, "]")
  } else {
    temp <- outer(rownames, colnames, FUN = paste, sep = ",")
    index_names <- paste0("[", temp, "]")
  }
  index_names
}

combine_duplicates <- function(x, sep = NULL) {
  # combine elements of a list that have the same name
  # NOTE: not used in the renaming functions anymore
  # Args:
  #   x: a named list
  #   sep: names of attributes to also include
  #        in the duplication checks
  # Returns: 
  #   a list of possibly reduced length.
  # Examples:
  #   combine_duplicates(list(a = 1, a = c(2,3)))
  #   becomes list(a = c(1,2,3)) 
  if (!is.list(x) || is.null(names(x))) 
    stop("x must be a named list")
  dat <- data.frame(names = names(x))
  for (s in sep) {
    dat[[s]] <- lapply(x, function(y)
      if (is.null(attr(y, s))) NA else attr(y, s)
    )
  }
  unique_dat <- unique(dat)
  new_list <- vector("list", length = nrow(unique_dat))
  names(new_list) <- unique_dat$names
  for (i in seq_along(new_list)) {
    new_list[[i]] <- unname(unlist(x[matching_rows(i, dat)]))
    for (s in sep) {
      if (!all(is.na(unique_dat[[s]][[i]]))) {
        # NAs were introduced above to allow matching
        # and should not be kept as attributes
        attr(new_list[[i]], s) <- unique_dat[[s]][[i]]
      }
    }
  }
  new_list
}

matching_rows <- function(i, data, check.attributes = FALSE, ...) {
  if (!length(nrow(data))) {
    stop("data must have rows")
  }
  ulapply(1:nrow(data), function(k, ...) is_equal(data[i, ], data[k, ], ...),
          check.attributes = check.attributes, ...)
}

make_dims <- function(x) {
  # helper function to make correct dims for .@sims$dims_oi
  if (is.null(x$dim)) {
    x$dim <- numeric(0)
  }
  setNames(rep(list(x$dim), length(x$pnames)), x$pnames)
}

do_renaming <- function(x, change) {
  # perform actual renaming of Stan parameters
  # Args:
  #   change: a list of lists each element allowing
  #           to rename certain parameters
  #   x: An object of class brmsfit
  # Returns:
  #   A brmsfit object with updated parameter names
  .do_renaming <- function(x, change) {
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
    x <- .do_renaming(x, change[[i]])
  }
  x
}

reorder_pars <- function(x) {
  # order parameter samples after parameter class
  # Args:
  #   x: brmsfit object
  all_classes <- c(
    "b", "bmo", "bcs", "bme", "ar", "ma", "arr", "lagsar",
    "errorsar", "car", "sdcar", "sigmaLL", "sd", "cor", "sds", 
    "sdgp", "lscale", dpars(), "temp", "rescor", "delta", 
    "lasso", "simo", "r", "s", "zgp", "rcar", "loclev", 
    "Xme", "prior", "lp"
  )
  # reorder parameter classes
  class <- get_matches("^[^[:digit:]_]+", x$fit@sim$pars_oi)
  new_order <- order(
    factor(class, levels = all_classes),
    !grepl("_Intercept$", x$fit@sim$pars_oi)
  )
  x$fit@sim$dims_oi <- x$fit@sim$dims_oi[new_order]
  x$fit@sim$pars_oi <- names(x$fit@sim$dims_oi)
  # reorder single parameter names
  nsubpars <- ulapply(x$fit@sim$dims_oi, prod)
  num <- lapply(seq_along(new_order), function(x)
    as.numeric(paste0(x, ".", sprintf("%010d", seq_len(nsubpars[x]))))
  )
  new_order <- order(unlist(num[order(new_order)]))
  x$fit@sim$fnames_oi <- x$fit@sim$fnames_oi[new_order]
  chains <- length(x$fit@sim$samples)
  for (i in seq_len(chains)) {
    # subset_attr ensures that attributes are not removed
    x$fit@sim$samples[[i]] <- subset_attr(x$fit@sim$samples[[i]], new_order)
  }
  x
}

compute_quantities <- function(x) {
  # wrapper function to compute and store quantities in the stanfit 
  # object which were not computed / stored by Stan itself
  # Args:
  #   x: a brmsfit object
  stopifnot(is.brmsfit(x))
  x <- recreate_xi(x)
  x
}

recreate_xi <- function(x) {
  # helper function to recreate parameter xi, which is currently
  # defined in the Stan model block and thus not being stored
  # Args:
  #   x: a brmsfit object
  # Returns:
  #   a brmsfit object, with temp_xi being replaced by xi 
  stopifnot(is.brmsfit(x))
  if (!"temp_xi" %in% parnames(x)) {
    return(x)
  }
  draws <- try(extract_draws(x))
  if (is(draws, "try-error")) {
    warning2("Trying to compute 'xi' was unsuccessful. ",
             "Some S3 methods may not work as expected.")
    return(x)
  }
  mu <- ilink(get_eta(draws, i = NULL), draws$f$link)
  sigma <- get_sigma(draws$sigma, data = draws$data, i = NULL)
  y <- matrix(draws$data$Y, dim(mu)[1], dim(mu)[2], byrow = TRUE)
  bs <- - 1 / matrixStats::rowRanges((y - mu) / sigma)
  bs <- matrixStats::rowRanges(bs)
  temp_xi <- as.vector(as.matrix(x, pars = "temp_xi"))
  xi <- inv_logit(temp_xi) * (bs[, 2] - bs[, 1]) + bs[, 1]
  # write xi into stanfit object
  samp_chain <- length(xi) / x$fit@sim$chains
  for (i in seq_len(x$fit@sim$chains)) {
    xi_part <- xi[((i - 1) * samp_chain + 1):(i * samp_chain)]
    # add warmup samples not used anyway
    xi_part <- c(rep(0, x$fit@sim$warmup2[1]), xi_part)
    x$fit@sim$samples[[i]][["temp_xi"]] <- xi_part
  }
  change <- list(
    pos = grepl("^temp_xi$", parnames(x)), 
    oldname = "temp_xi", 
    fnames = "xi", 
    pnames = "xi"           
  )
  do_renaming(x, list(change))
}
