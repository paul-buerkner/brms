rename <- function(names, symbols = NULL, subs = NULL, 
                   fixed = TRUE, check_dup = FALSE) {
  # rename certain symbols in a character vector
  # Args:
  #   names: a character vector of names to be renamed
  #   symbols: the regular expressions in names to be replaced
  #   subs: the replacements
  #   fixed: same as for sub, grepl etc
  #   check_dup: logical; check for duplications in names after renaming
  # Returns: 
  #   renamed parameter vector of the same length as names
  symbols <- as.character(symbols)
  subs <- as.character(subs)
  if (!length(symbols)) {
    symbols <- c(" ", "(", ")", "[", "]", ",", "\"", "'", 
                 "+", "-", "*", "/", "^", "=", "!=")
  }
  if (!length(subs)) {
    subs <- c(rep("", 8), "P", "M", "MU", "D", "E", "EQ", "NEQ")
  }
  if (length(subs) == 1L) {
    subs <- rep(subs, length(symbols))
  }
  stopifnot(length(symbols) == length(subs))
  # avoid zero-length pattern error
  has_chars <- nzchar(symbols)
  symbols <- symbols[has_chars]
  subs <- subs[has_chars]
  new_names <- names
  for (i in seq_along(symbols)) {
    new_names <- gsub(symbols[i], subs[i], new_names, fixed = fixed)
  }
  dup <- duplicated(new_names)
  if (check_dup && any(dup)) {
    dup_names <- names[new_names %in% new_names[dup]]
    stop("Internal renaming of variables led to duplicated names. \n",
         "Occured for variables: ", paste(dup_names, collapse = ", "))
  }
  new_names
}

model_name <- function(family) {
  # create the name of the fitted stan model
  # Args:
  #   family: A family object
  if (!is(family, "family")) {
    mn <- "brms-model"
  } else {
    type <- ifelse(is.null(family$type), "", paste(",", family$type))
    mn <- paste0(family$family, "(",family$link, type, ") brms-model")
  }
  mn
}

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
  # some variables generally needed
  family <- family(x)
  ee <- extract_effects(x$formula, family = family)
  standata <- standata(x)
  
  # order parameter samples after parameter class
  chains <- length(x$fit@sim$samples) 
  all_classes <- c("b_Intercept", "b", "bm", "bcs", "ar", "ma", "arr", "sd", 
                   "cor", "sds", "sigma", "sigmaLL", "rescor", "nu", "shape",
                   "phi", "zi", "hu", "delta", "simplex", "r", "s", "loclev", 
                   "prior", "lp")
  class <- get_matches("^[^_\\[]+", x$fit@sim$fnames_oi)
  # make sure that the fixed effects intercept comes first
  if (length(ee$response) > 1L) {
    regex_resp <- paste0(usc(ee$response , "suffix"), collapse  = "|")
    regex_resp <- paste0("(", regex_resp, ")")
  } else {
    regex_resp <- ""
  }
  regex_intercept <- paste0("^b_", regex_resp, "Intercept($|\\[)")
  pos_intercept <- which(grepl(regex_intercept, x$fit@sim$fnames_oi))
  class[pos_intercept] <- "b_Intercept"
  ordered <- order(factor(class, levels = all_classes))
  x$fit@sim$fnames_oi <- x$fit@sim$fnames_oi[ordered]
  for (i in seq_len(chains)) {
    # subset_attr ensures that attributes are not removed
    x$fit@sim$samples[[i]] <- subset_attr(x$fit@sim$samples[[i]], ordered)
  }
  mclass <- regmatches(x$fit@sim$pars_oi, regexpr("^[^_]+", x$fit@sim$pars_oi))
  pos_intercept <- which(grepl(regex_intercept, x$fit@sim$pars_oi))
  mclass[pos_intercept] <- "b_Intercept"
  ordered <- order(factor(mclass, levels = all_classes))
  x$fit@sim$dims_oi <- x$fit@sim$dims_oi[ordered]
  x$fit@sim$pars_oi <- names(x$fit@sim$dims_oi)
  pars <- parnames(x)
  
  # find positions of parameters and define new names
  change <- list()
  if (length(ee$nonlinear)) {
    nlpars <- names(ee$nonlinear)
    for (p in nlpars) {
      splines <- get_spline_labels(ee$nonlinear[[p]], x$data, covars = TRUE)
      change_eff <- change_effects(
        pars = pars, dims = x$fit@sim$dims_oi,
        fixef = colnames(standata[[paste0("X_", p)]]),
        monef = colnames(standata[[paste0("Xm_", p)]]),
        splines = splines, nlpar = p)
      change <- c(change, change_eff)
    }
  } else {
    resp <- ee$response
    if (length(resp) > 1L) {
      for (r in resp) {
        fixef <- rm_int_fixef(colnames(standata[[paste0("X_", r)]]), 
                              stancode(x), nlpar = r)
        change_eff <- change_effects(
          pars = pars, dims = x$fit@sim$dims_oi, 
          fixef = fixef, monef = colnames(standata[[paste0("Xm_", r)]]),
          splines = get_spline_labels(ee, x$data, covars = TRUE),
          nlpar = r)
        change <- c(change, change_eff)
      }
    } else {
      fixef <- rm_int_fixef(colnames(standata[["X"]]), stancode(x))
      change_eff <- change_effects(
        pars = pars, dims = x$fit@sim$dims_oi, 
        fixef = fixef, monef = colnames(standata[["Xm"]]),
        splines = get_spline_labels(ee, x$data, covars = TRUE))
      change_csef <- change_csef(colnames(standata[["Xcs"]]), 
                                 pars = pars, ncat = standata$ncat)
      change <- c(change, change_eff, change_csef)
    }
  }
  # rename parameters related to auxilliary parameters
  for (ap in intersect(auxpars(), names(ee))) {
    change_eff <- change_effects(
      pars = pars, dims = x$fit@sim$dims_oi,
      fixef = colnames(standata[[paste0("X_", ap)]]),
      monef = colnames(standata[[paste0("Xm_", ap)]]),
      splines = get_spline_labels(ee[[ap]], x$data, covars = TRUE),
      nlpar = ap)
    change <- c(change, change_eff)
  }
  # rename group-level parameters separately
  change_ranef <- change_ranef(x$ranef, pars = pars, 
                               dims = x$fit@sim$dims_oi)
  change <- c(change, change_ranef)
  
  if (is.linear(family) && length(ee$response) > 1L) {
    # rename residual parameters of multivariate linear models
    corfnames <- paste0("sigma_", ee$response)
    change <- lc(change, 
      list(pos = grepl("^sigma\\[", pars), oldname = "sigma",
           pnames = corfnames, fnames = corfnames))
    change <- c(change, change_prior(class = "sigma", pars = pars, 
                                     names = ee$response))
    rescor_names <- get_cornames(ee$response, type = "rescor", 
                                  brackets = FALSE)
    change <- lc(change, list(pos = grepl("^rescor\\[", pars), 
                              oldname = "rescor", pnames = rescor_names,
                              fnames = rescor_names))
  }
  
  # perform the actual renaming in x$fit@sim
  for (i in seq_along(change)) {
    x <- do_renaming(change = change[[i]], x = x)
  }
  x$fit@sim$pars_oi <- names(x$fit@sim$dims_oi)
  x
}

change_effects <- function(pars, dims, fixef = NULL, monef = NULL, 
                           splines = NULL, nlpar = "") {
  # helps in renaming various kinds of effects
  change_fixef <- change_fixef(fixef, pars = pars, nlpar = nlpar)
  change_monef <- change_monef(monef, pars = pars, nlpar = nlpar)
  change_splines <- change_splines(splines, pars = pars, nlpar = nlpar)
  c(change_fixef, change_monef, change_splines)
}

change_fixef <- function(fixef, pars, nlpar = "") {
  # helps in renaming fixed effects parameters
  # Args:
  #   fixef: names of the fixed effects
  #   pars: names of all model parameters
  #   nlpar: optional string naming a non-linear parameter
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  change <- list()
  if (length(fixef)) {
    b <- paste0("b", usc(nlpar, "prefix"))
    pos <- grepl(paste0("^", b, "\\["), pars)
    bnames <- paste0(b, "_", fixef)
    change <- lc(change, 
      list(pos = pos, oldname = b, pnames = bnames, fnames = bnames))
    change <- c(change,
      change_prior(class = b, pars = pars, names = fixef))
  }
  change
}

change_monef <- function(monef, pars, nlpar = "") {
  # helps in renaming monotonous effects parameters
  # Args:
  #   monef: names of the monotonous effects
  #   pars: names of all model parameters
  #   nlpar: optional string naming a non-linear parameter
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  change <- list()
  if (length(monef)) {
    p <- usc(nlpar, "prefix")
    bm <- paste0("bm", p)
    newnames <- paste0("bm", p, "_", monef)
    change <- lc(change, list(pos = grepl(paste0("^", bm, "\\["), pars), 
                              oldname = bm, pnames = newnames, 
                              fnames = newnames))
    change <- c(change, change_prior(class = bm, pars = pars, names = monef))
    for (i in seq_along(monef)) {
      simplex <- paste0(paste0("simplex", p, "_"), c(i, monef[i]))
      pos <- grepl(paste0("^", simplex[1], "\\["), pars)
      change <- lc(change, 
        list(pos = pos, oldname = simplex[1], pnames = simplex[2], 
             fnames = paste0(simplex[2], "[", 1:sum(pos), "]"), 
             dim = sum(pos)))
      change <- c(change,
        change_prior(class = simplex[1], new_class = simplex[2],
                      pars = pars, is_vector = TRUE))
    }
  }
  change
}

change_csef <- function(csef, pars, ncat) {
  # helps in renaming category specific effects parameters
  # Args:
  #   cseef: names of the category specific effects
  #   pars: names of all model parameters
  #   ncat: number of response categories
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  change <- list()
  if (length(csef)) {
    ncse <- length(csef)
    thres <- ncat - 1
    csenames <- t(outer(csef, paste0("[", 1:thres, "]"), FUN = paste0))
    csenames <- paste0("bcs_", csenames)
    sort_cse <- ulapply(seq_len(ncse), seq, to = thres * ncse, by = ncse)
    change <- lc(change, list(pos = grepl("^bcs\\[", pars), oldname = "bcs", 
                              pnames = paste0("bcs_", csef), fnames = csenames,
                              sort = sort_cse, dim = thres))
    change <- c(change, change_prior(class = "bcs", pars = pars, names = csef))
  }
  change
}

change_splines <- function(splines, pars, nlpar = "") {
  # helps in renaming smoothing term parameters
  # Args:
  #   splines: smoothing terms
  #   pars: names of all model parameters
  #   nlpar: optional string naming a non-linear parameter
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  change <- list()
  if (length(splines)) {
    stopifnot(!is.null(attr(splines, "nbases")))
    p <- usc(nlpar, "prefix")
    splines <- rename(splines)
    sds <- paste0("sds", p)
    sds_names <- paste0(sds, "_", splines)
    s <- paste0("s", p)
    s_names <- paste0(s, "_", splines)
    for (i in seq_along(splines)) {
      nb <- attr(splines, "nbases")[[i]]
      for (j in seq_len(nb)) {
        ij <- paste0(i, "_", j)
        sds_pos <- grepl(paste0("^", sds, "_", ij), pars)
        change <- lc(change, 
          list(pos = sds_pos, oldname = paste0(sds, "_", ij), 
               pnames = paste0(sds_names[i], "_", j), 
               fnames = paste0(sds_names[i], "_", j)))
        s_pos <- grepl(paste0("^", s, "_", ij), pars)
        s_fnames <- paste0(s_names[i], "_", j, "[", seq_len(sum(s_pos)), "]")
        change <- lc(change, 
          list(pos = s_pos, oldname = paste0(s, "_", ij), 
               pnames = s_names[i], fnames = s_fnames, 
               dim = as.numeric(sum(s_pos))))
        new_prior_class <- paste0(sds, "_", splines[i], "_", j)
        change <- c(change, 
          change_prior(class = paste0(sds, "_", ij), pars = pars, 
                       new_class = new_prior_class)) 
      }
    }
  }
  change
}
  

change_ranef <- function(ranef, pars, dims) {
  # helps in renaming group-level parameters
  # Args:
  #   ranef: list returned by tidy_ranef
  #   pars: names of all model parameters
  #   dims: named list containing parameter dimensions
  #   nlpar: optional string naming a non-linear parameter
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  change <- list()
  if (nrow(ranef)) {
    for (id in unique(ranef$id)) {
      r <- ranef[ranef$id == id, ]
      g <- r$group[1]
      suffix <- paste0(usc(r$nlpar, "suffix"), r$coef)
      rfnames <- paste0("sd_", g, "__", suffix)
      change <- lc(change,
        list(pos = grepl(paste0("^sd_", id, "(\\[|$)"), pars),
             oldname = paste0("sd_", id), pnames = rfnames,
             fnames = rfnames))
      change <- c(change,
        change_prior(class = paste0("sd_", id), pars = pars,
                     new_class = paste0("sd_", g), names = suffix))
      # rename group-level correlations
      if (nrow(r) > 1L && isTRUE(r$cor[1])) {
        type <- paste0("cor_", g)
        cor_names <- get_cornames(suffix, type = paste0("cor_", g),
                                  brackets = FALSE)
        change <- lc(change,
          list(pos = grepl(paste0("^cor_", id, "(\\[|$)"), pars),
               oldname = paste0("cor_", id), pnames = cor_names,
               fnames = cor_names))
        change <- c(change,
          change_prior(class = paste0("cor_", id), pars = pars,
                       new_class = paste0("cor_", g)))
      }
    }
    if (any(grepl("^r_", pars))) {
      rc_args <- nlist(ranef, pars, dims)
      change <- c(change, do.call(change_ranef_levels, rc_args))
    }
  }
  change
} 

change_ranef_levels <- function(ranef, pars, dims)  {
  # helps in renaming random effects 'r_.' parameters
  # Args:
  #   ranef: output of tidy_ranef
  #   pars: names of all model parameters
  #   dims: named list containing parameter dimensions
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  change <- list()
  for (i in seq_len(nrow(ranef))) {
    r <- ranef[i, ]
    usc_nlpar <- usc(r$nlpar)
    r_parnames <- paste0("r_", r$id, usc_nlpar, "_", r$cn)
    r_regex <- paste0("^", r_parnames, "(\\[|$)")
    change_rl <- list(pos = grepl(r_regex, pars), 
                      oldname = r_parnames)
    r_new_parname <- paste0("r_", r$group, usc(usc_nlpar))
    # prepare for removal of redundant parameters r_<i>
    # and for combining group-level effects into one parameter matrix
    gf_matches <- which(ranef$group == r$group & ranef$nlpar == r$nlpar)
    nr <- length(gf_matches)
    old_dim <- dims[[change_rl$oldname]]
    if (gf_matches[1] == i) {
      # if this is the first RE term of this group
      # counted separately for each non-linear parameter
      change_rl$pnames <- r_new_parname
      change_rl$dim <- c(old_dim[1], nr)
    }
    # rstan doesn't like whitespaces in parameter names
    levels <- gsub("[ \t\r\n]", ".", attr(ranef, "levels")[[r$group]])
    index_names <- make_index_names(levels, r$coef, dim = 2)
    change_rl$fnames <- paste0(r_new_parname, index_names)
    change <- lc(change, change_rl)
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
  pos_priors <- which(grepl(paste0("^prior_", class, "(_|$|\\[)"), pars))
  if (length(pos_priors)) {
    priors <- gsub(paste0("^prior_", class), paste0("prior_", new_class), 
                   pars[pos_priors])
    if (is_vector) {
      change <- lc(change, 
                   list(pos = pos_priors, oldname = paste0("prior_", class),
                        pnames = priors, fnames = priors))
    } else {
      digits <- sapply(priors, function(prior) {
        d <- regmatches(prior, gregexpr("_[[:digit:]]+$", prior))[[1]]
        if (length(d)) 
          as.numeric(substr(d, 2, nchar(d))) 
        else 0
      })
      for (i in seq_along(priors)) {
        if (digits[i] && !is.null(names)) {
          priors[i] <- gsub("[[:digit:]]+$", names[digits[i]], priors[i])
        }
        if (pars[pos_priors[i]] != priors[i]) {
          change <- lc(change, 
            list(pos = pos_priors[i], oldname = pars[pos_priors[i]],
                 pnames = priors[i], fnames = priors[i]))
        }
      }
    }
  }
  change
}

change_old_ranef <- function(ranef, pars, dims) {
  # prepare for renaming of group-level parameters of models 
  # fitted with brms <= 0.10.0.9000
  # only relevant for non-linear models
  # Args:
  #   ranef: output of tidy_ranef
  #   pars: names of all parameters in the model
  #   dims: dimension of parameters
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  change <- list()
  for (id in unique(ranef$id)) {
    r <- ranef[ranef$id == id, ]
    g <- r$group[1]
    nlpar <- r$nlpar[1]
    # rename sd-parameters
    old_sd_names <- paste0("sd_", nlpar, "_", g, "_", r$coef)
    new_sd_names <- paste0("sd_", g, "__", nlpar, "_", r$coef)
    for (i in seq_along(old_sd_names)) {
      change <- lc(change, 
        change_simple(old_sd_names[i], new_sd_names[i], pars, dims))
    }
    # rename cor-parameters
    new_cor_names <- get_cornames(paste0(nlpar, "_", r$coef),
                                  type = paste0("cor_", g),
                                  brackets = FALSE)
    old_cor_names <- get_cornames(r$coef, brackets = FALSE, sep = "_",
                                  type = paste0("cor_", nlpar, "_", g))
    for (i in seq_along(old_cor_names)) {
      change <- lc(change, 
        change_simple(old_cor_names[i], new_cor_names[i], pars, dims))
    } 
    # rename r-parameters
    old_r_name <- paste0("r_", nlpar, "_", g)
    new_r_name <- paste0("r_", g, "__", nlpar)
    levels <- gsub("[ \t\r\n]", ".", attr(ranef, "levels")[[g]])
    index_names <- make_index_names(levels, r$coef, dim = 2)
    new_r_names <- paste0(new_r_name, index_names)
    change <- lc(change, 
      change_simple(old_r_name, new_r_names, pars, dims, 
                    pnames = new_r_name))
  }
  change
}

change_old_ranef2 <- function(ranef, pars, dims) {
  # prepare for renaming of group-level parameters of models 
  # fitted with brms > 0.10.0.9000 but < 1.0.0
  # Args:
  #   ranef: output of tidy_ranef
  #   pars: names of all parameters in the model
  #   dims: dimension of parameters
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  change <- list()
  for (id in unique(ranef$id)) {
    r <- ranef[ranef$id == id, ]
    g <- r$group[1]
    #nlpar <- r$nlpar[1]
    nlpars_usc <- usc(r$nlpar, "suffix")
    # rename sd-parameters
    old_sd_names <- paste0("sd_", g, "_", nlpars_usc, r$coef)
    new_sd_names <- paste0("sd_", g, "__", nlpars_usc, r$coef)
    for (i in seq_along(old_sd_names)) {
      change <- lc(change, 
        change_simple(old_sd_names[i], new_sd_names[i], pars, dims))
    }
    # rename cor-parameters
    new_cor_names <- get_cornames(paste0(nlpars_usc, r$coef),
                                  type = paste0("cor_", g),
                                  brackets = FALSE)
    old_cor_names <- get_cornames(paste0(nlpars_usc, r$coef),
                                  type = paste0("cor_", g),
                                  brackets = FALSE, sep = "_")
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
        change_simple(old_r_name, new_r_names, pars, dims, 
                      pnames = new_r_name))
    }
  }
  change
}

change_old_splines <- function(effects, pars, dims) {
  # change names of spline parameters fitted with brms <= 1.0.1
  # this became necessary after allowing splines with multiple covariates
  .change_old_splines <- function(e, nlpar = "") {
    change <- list()
    spline_labels <- get_spline_labels(e)
    if (length(spline_labels)) {
      p <- usc(nlpar, "suffix")
      old_splines <- rename(paste0(p, spline_labels))
      new_splines <- rename(paste0(p, get_spline_labels(e, covars = TRUE)))
      old_sds_pars <- paste0("sds_", old_splines)
      new_sds_pars <- paste0("sds_", new_splines, "_1")
      old_s_pars <- paste0("s_", old_splines)
      new_s_pars <- paste0("s_", new_splines, "_1")
      for (i in seq_along(old_splines)) {
        change <- lc(change,
          change_simple(old_sds_pars[i], new_sds_pars[i], pars, dims))
        indices <- seq_len(dims[[old_s_pars[i]]])
        new_s_par_indices <- paste0(new_s_pars[i], "[", indices, "]")
        change <- lc(change,
          change_simple(old_s_pars[i], new_s_par_indices, pars, dims,
                        pnames = new_s_pars[i]))
      }
    }
    return(change)
  }
  
  change <- list()
  spec_effects <- rmNULL(c(effects[auxpars()], effects$nonlinear))
  for (sp in names(spec_effects)) {
    se <- spec_effects[[sp]]
    change <- c(change, .change_old_splines(se, nlpar = sp))
  }
  resp <- effects$response
  if (length(resp) > 1L) {
    for (r in resp) {
      change <- c(change, .change_old_splines(effects, nlpar = r))
    }
  } else {
    change <- c(change, .change_old_splines(effects))
  }
  change
}

change_simple <- function(oldname, fnames, pars, dims,
                          pnames = fnames) {
  # helper function for very simple renaming
  # currently only used in the change_old_ranef functions
  pos <- grepl(paste0("^", oldname), pars)
  if (any(pos)) {
    out <- nlist(pos, oldname, pnames, fnames,
                 dims = dims[[oldname]])
  } else {
    out <- NULL
  }
  return(out)
}

rm_int_fixef <- function(fixef, stancode, nlpar = "") {
  # identifies if the intercept has to be removed from fixef
  # and returns adjusted fixef names
  p <- usc(nlpar, "suffix")
  regex <- paste0("b_", p, "Intercept = temp_", p, "Intercept")
  if (grepl(regex, stancode, fixed = TRUE)) {
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
      if (is.null(attr(y, s))) NA else attr(y, s))
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

do_renaming <- function(change, x) {
  # perform actual renaming of parameters
  # Args:
  #   change: A list containing all information to rename 
  #           a certain type of parameters (e.g., fixed effects)
  #   x: An object of class brmsfit
  # Returns:
  #   A brmsfit object with updated parameter names
  chains <- length(x$fit@sim$samples) 
  x$fit@sim$fnames_oi[change$pos] <- change$fnames
  for (i in 1:chains) {
    names(x$fit@sim$samples[[i]])[change$pos] <- change$fnames
    if (!is.null(change$sort)) {
      x$fit@sim$samples[[i]][change$pos] <- 
        x$fit@sim$samples[[i]][change$pos][change$sort]
    }
  }
  onp <- match(change$oldname, names(x$fit@sim$dims_oi))
  if (is.na(onp)) {
    warning("Parameter ", change$oldname, " could not be renamed. ",
            "This should not happen. \nPlease inform me so that ",
            "I can fix this problem.", call. = FALSE)
  } else {
    if (is.null(change$pnames)) {
      # only needed to collapse multiple r_<i> of the same grouping factor
      x$fit@sim$dims_oi[[onp]] <- NULL  
    } else { 
      # rename dims_oi to match names in fnames_oi
      dims <- x$fit@sim$dims_oi
      x$fit@sim$dims_oi <- c(if (onp > 1) dims[1:(onp - 1)], 
                             make_dims(change),
                             dims[(onp + 1):length(dims)])
    }
  }
  x
}
