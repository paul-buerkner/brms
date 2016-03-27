rename <- function(names, symbols = NULL, subs = NULL, 
                   fixed = TRUE, check_dup = FALSE) {
  # rename certain symbols in a character vector
  # 
  # Args:
  #   names: a character vector of names to be renamed
  #   symbols: the regular expressions in names to be replaced
  #   subs: the replacements
  #   fixed: same as for sub, grepl etc
  #   check_dup: logical; check for duplications in names after renaming
  # 
  # Returns: 
  #   renamed parameter vector of the same length as names
  if (is.null(symbols)) {
    symbols <- c(" ", "(", ")", "[", "]", ",", 
                 "+", "-", "*", "/", "^", "=", "!=")
  }
  if (is.null(subs)) {
    subs <- c(rep("", 6), "P", "M", "MU", "D", "E", "EQ", "NEQ")
  }
  if (length(subs) == 1) {
    subs <- rep(subs, length(symbols))
  }
  if (length(symbols) != length(subs)) {
    stop("length(symbols) != length(subs)")
  }
  new_names <- names
  for (i in seq_along(symbols)) {
    # avoid zero-length pattern error when nchar(symbols[i]) == 0
    if (nchar(symbols[i])) {
      new_names <- gsub(symbols[i], subs[i], new_names, fixed = fixed)
    }
  }
  dup <- duplicated(new_names)
  if (check_dup && any(dup)) 
    stop(paste0("Internal renaming of variables led to duplicated names. \n",
                "Occured for variables: ", 
                paste(names[which(new_names %in% new_names[dup])], 
                      collapse = ", ")))
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
  # rename parameters (and possibly change their dimensions) within the stanfit object 
  # to ensure reasonable parameter names for summary, plot, launch_shiny etc.
  # 
  # Args:
  #   x a brmsfit obejct
  #
  # Returns:
  #   a brmfit object with adjusted parameter names and dimensions
  if (!length(x$fit@sim)) return(x) 
  
  # order parameter samples after parameter class
  chains <- length(x$fit@sim$samples) 
  all_classes <- c("b_Intercept", "b", "bm", "bp", "ar", "ma", "arr", 
                   "sd", "cor", "sigma", "rescor", "nu", "shape", "phi",
                   "delta", "simplex", "r", "prior", "lp")
  class <- regmatches(x$fit@sim$fnames_oi, regexpr("^[^_\\[]+", x$fit@sim$fnames_oi))
  # make sure that the fixed effects intercept comes first
  pos_intercept <- which(grepl("^b_Intercept($|\\[)", x$fit@sim$fnames_oi))
  class[pos_intercept] <- "b_Intercept"
  ordered <- order(factor(class, levels = all_classes))
  x$fit@sim$fnames_oi <- x$fit@sim$fnames_oi[ordered]
  for (i in 1:chains) {
    # subset_attr ensures that attributes are not removed
    x$fit@sim$samples[[i]] <- subset_attr(x$fit@sim$samples[[i]], ordered)
  }
  mclass <- regmatches(x$fit@sim$pars_oi, regexpr("^[^_]+", x$fit@sim$pars_oi))
  # make sure that the fixed effects intercept comes first
  pos_intercept <- which(grepl("^b_Intercept($|\\[)", x$fit@sim$pars_oi))
  mclass[pos_intercept] <- "b_Intercept"
  ordered <- order(factor(mclass, levels = all_classes))
  x$fit@sim$dims_oi <- x$fit@sim$dims_oi[ordered]
  x$fit@sim$pars_oi <- names(x$fit@sim$dims_oi)
  
  # some variables generally needed
  pars <- parnames(x)
  family <- family(x)
  ee <- extract_effects(x$formula, family = family, nonlinear = x$nonlinear)
  change <- list()
  standata <- standata(x)
  
  # find positions of parameters and define new names
  # fixed effects
  if (length(x$nonlinear)) {
    nlpars <- paste0("_", names(ee$nonlinear))
    X_list <- standata[paste0("X", nlpars)]
  } else {
    nlpars <- ""
    X_list <- list(standata$X)
  }
  f <- lapply(X_list, colnames)
  for (i in seq_along(f)) {
    if (length(f[[i]])) {
      b <- paste0("b", nlpars[i])
      change <- lc(change, list(pos = grepl(paste0("^", b, "\\["), pars), 
                                oldname = b, pnames = paste0(b, "_", f[[i]]), 
                                fnames = paste0(b, "_", f[[i]])))
      change <- c(change, prior_changes(class = "b", pars = pars, names = f))
    }
  }
  # non-linear models currently don't use any special intercept handling
  intercepts <- names(get_intercepts(ee, data = x$data, family = family))
  if (length(intercepts) && !is_equal(intercepts, "Intercept")) {
    # for intercepts in models using multivariate formula syntax
    change <- lc(change, list(pos = grepl("^b_Intercept\\[", pars), 
                              oldname = "b_Intercept", 
                              pnames = paste0("b_", intercepts), 
                              fnames = paste0("b_", intercepts)))
  }
  # monotonous effects
  if (is.formula(ee$mono)) {
    monef <- colnames(standata$Xm)
    if (length(monef)) {
      change <- lc(change, list(pos = grepl("^bm\\[", pars), oldname = "bm", 
                                pnames = paste0("b_", monef), 
                                fnames = paste0("b_", monef)))
      change <- c(change, prior_changes(class = "bm", pars = pars, 
                                        names = monef))
      for (i in seq_along(monef)) {
        simplex <- paste0("simplex_", c(i, monef[i]))
        pos <- grepl(paste0("^", simplex[1], "\\["), pars)
        change <- lc(change, 
          list(pos = pos, oldname = simplex[1], pnames = simplex[2], 
               fnames = paste0(simplex[2], "[", 1:sum(pos), "]"), 
               dim = sum(pos)))
        change <- c(change,
          prior_changes(class = simplex[1], new_class = simplex[2],
                        pars = pars, is_vector = TRUE))
      }
    }
  } 
  # category specific effects
  if (is.formula(ee$cse)) {
    csef <- colnames(standata$Xp)
    ncse <- length(csef)
    if (ncse) {
      thres <- max(standata$max_obs) - 1
      csenames <- t(outer(csef, paste0("[",1:thres,"]"), FUN = paste0))
      csenames <- paste0("b_", csenames)
      sort_cse <- ulapply(1:ncse, seq, to = thres * ncse, by = ncse)
      change <- lc(change, list(pos = grepl("^bp\\[", pars), oldname = "bp", 
                                pnames = paste0("b_", csef), fnames = csenames,
                                sort = sort_cse, dim = thres))
      change <- c(change, prior_changes(class = "bp", pars = pars, 
                                        names = csef))
    }
  } 
  # random effects
  if (length(x$ranef)) {
    group <- names(x$ranef)
    gf <- make_group_frame(x$ranef)
    random <- get_random(ee)
    nlpar <- ulapply(x$ranef, get_nlpar, suffix = "_")
    for (i in seq_along(x$ranef)) {
      # k = i for ordinary (not non-linear) models
      k <- get_re_index(i, random)
      sd <- paste0("sd_", nlpar[i])
      rfnames <- paste0(sd, group[i], "_", x$ranef[[i]])
      change <- lc(change, 
        list(pos = grepl(paste0("^", sd, k,"(\\[|$)"), pars),
             oldname = paste0(sd, k), pnames = rfnames, 
             fnames = rfnames))
      change <- c(change, 
        prior_changes(class = paste0(sd, k), pars = pars, 
                      new_class = paste0(sd, group[i]),
                      names = x$ranef[[i]]))
      # rename random effects correlations
      if (length(x$ranef[[i]]) > 1 && random$cor[[i]]) {
        cor <- paste0("cor_", nlpar[i])
        cor_names <- get_cornames(x$ranef[[i]], brackets = FALSE,
                                  type = paste0(cor, group[i]))
        change <- lc(change, 
          list(pos = grepl(paste0("^", cor, k, "(\\[|$)"), pars),
               oldname = paste0(cor, k), pnames = cor_names,
               fnames = cor_names)) 
        change <- c(change, 
          prior_changes(class = paste0(cor, k), pars = pars, 
                        new_class = paste0(cor, group[i])))
      }
      if (any(grepl("^r_", pars))) {
        rc_args <- nlist(i, k, gf, pars, ranef = x$ranef, 
                         dims_oi = x$fit@sim$dims_oi)
        change <- lc(change, do.call(ranef_changes, rc_args))
      }  
    }
  }
  
  if (has_sigma(family, se = is.formula(ee$se), autocor = x$autocor)) {
    corfnames <- paste0("sigma_",ee$response)
    change <- lc(change, list(pos = grepl("^sigma", pars), oldname = "sigma",
                              pnames = corfnames, fnames = corfnames))
    change <- c(change, prior_changes(class = "sigma", pars = pars, 
                                      names = ee$response))
    # residual correlation paramaters
    if (is.linear(family) && length(ee$response) > 1) {
       rescor_names <- get_cornames(ee$response, type = "rescor", brackets = FALSE)
       change <- lc(change, list(pos = grepl("^rescor\\[", pars), oldname = "rescor",
                                 pnames = rescor_names, fnames = rescor_names))
    }
  } 
  
  # perform the actual renaming in x$fit@sim
  for (i in seq_along(change)) {
    x <- do_renaming(change = change[[i]], x = x)
  }
  x$fit@sim$pars_oi <- names(x$fit@sim$dims_oi)
  if (length(x$ranef)) { 
    # combines duplicated grouping factors to appear as if it was only one
    x$ranef <- combine_duplicates(x$ranef, sep = c("nlpar", "levels"))
  }
  x
}

make_group_frame <- function(ranef) {
  # make a little data.frame helping to rename and combine random effects
  # 
  # Args:
  #   ranef: a named list containing the random effects. The names are taken as grouping factors
  #
  # Returns: 
  #   A data.frame with length(ranef) rows and 3 columns: 
  #     g: the grouping factor of each terms
  #     nlpars: the non-linear parameter of each term (if present)
  #     first: a number corresponding to the first column for this term 
  #            in the final r_<gf> matrices
  #     last: a number corresponding to the last column for this term 
  #           in the final r_<gf> matrices
  group <- names(ranef)
  nlpar <- ulapply(ranef, function(y) attr(y, "nlpar"))
  if (is.null(nlpar)) nlpar <- rep("", length(ranef))
  out <- data.frame(g = group, nlp = nlpar, first = NA, last = NA)
  out[1, 3:4] <- c(1, length(ranef[[1]]))
  if (length(group) > 1) {
    for (i in 2:length(group)) {
      matches <- which(out$g[1:(i-1)] == group[i] & 
                       out$nlp[1:(i-1)] == nlpar[i])
      if (length(matches))
        out[i, 3:4] <- c(out$last[max(matches)] + 1, 
                         out$last[max(matches)] + length(ranef[[i]]))
      else out[i, 3:4] <- c(1, length(ranef[[i]]))
    }
  }
  out
}

make_index_names <- function(rownames, colnames = NULL, dim = 1) {
  # compute index names in square brackets for indexing stan parameters
  #
  # Args:
  #   rownames: a vector of row names
  #   colnames: a vector of columns 
  #   dim: The number of dimensions of the output either 1 or 2
  #
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
  #
  # Args:
  #   x: a named list
  #   sep: names of attributes to also include
  #        in the duplication checks
  #
  # Returns: 
  #   a list of possibly reduced length.
  # 
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

ranef_changes <- function(i, ranef, gf, dims_oi, pars, k = i)  {
  # helps in renaming random effects (r_) parameters
  # Args:
  #  i: the global index of the grouping factor under consideration
  #  ranef: output of gather_ranef
  #  gf: matrix as constructed by make_group_frame
  #  dims_oi: named list containing parameter dimensions
  #  k: the index in the Stan parameter names 
  #     may differ from i only for non-linear models
  # Returns:
  #  a list that can be interpreted by rename_pars
  group <- names(ranef)
  r <- "r_"
  nlpar <- attr(ranef[[i]], "nlpar")
  if (!is.null(nlpar)) {
    r <- paste0(r, nlpar, "_")
  } else nlpar <- ""
  r_parnames <- paste0("^", r, k, "(\\[|$)")
  change <- list(pos = grepl(r_parnames, pars), oldname = paste0(r, k))
  
  # prepare for removal of redundant parameters r_<i>
  # and for combining random effects into one paramater matrix
  # number of total REs for this grouping factor
  gf_matches <- which(gf$g == group[i] & gf$nlp == nlpar)
  n_ranefs <- max(gf$last[gf_matches]) 
  old_dim <- dims_oi[[change$oldname]]
  if (gf_matches[1] == i) {
    # if this is the first RE term of this group
    # counted separately for each non-linear parameter
    change$pnames <- paste0(r, group[i])
    change$dim <- if (n_ranefs == 1) old_dim else c(old_dim[1], n_ranefs) 
  } 
  # rstan doesn't like whitespaces in parameter names
  level_names <- gsub("[ \t\r\n]", ".", attr(ranef[[i]], "levels"))
  index_names <- make_index_names(rownames = level_names,
                                  colnames = ranef[[i]], dim = 2)
  change$fnames <- paste0(r, group[i], index_names)
  change
}

prior_changes <- function(class, pars, names = NULL, new_class = class,
                          is_vector = FALSE) {
  # helps in renaming prior parameters
  #
  # Args: 
  #   class: the class of the parameters for which prior names should be changed
  #   pars: all parameters in the model
  #   names: names to replace digits at the end of parameter names
  #   new_class: replacement of the orginal class name
  #   is_vector: indicate if the prior parameter is a vector
  #
  # Return:
  #   a list whose elements can be interpreted by rename_pars
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
      if (sum(abs(digits)) > 0 && is.null(names)) {
        stop("argument 'names' is missing")
      }
      for (i in 1:length(priors)) {
        if (digits[i]) {
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

make_dims <- function(x) {
  # helper function to make correct dims for .@sims$dims_oi
  if (is.null(x$dim)) 
    x$dim <- numeric(0)
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
    stop(paste("Parameter", change$oldname, "could not be renamed.",
               "Please report a bug."), call. = FALSE)
  }
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
  x
}