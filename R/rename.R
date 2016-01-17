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
                 "delta", "r", "prior", "lp")
  class <- regmatches(x$fit@sim$fnames_oi, regexpr("^[^_\\[]+", x$fit@sim$fnames_oi))
  # make sure that the fixed effects intercept comes first
  pos_intercept <- which(grepl("^b_Intercept($|\\[)", x$fit@sim$fnames_oi))
  class[pos_intercept] <- "b_Intercept"
  ordered <- order(factor(class, levels = all_classes))
  x$fit@sim$fnames_oi <- x$fit@sim$fnames_oi[ordered]
  for (i in 1:chains) {
    # keep_attr ensures that attributes are not removed
    x$fit@sim$samples[[i]] <- keep_attr(x$fit@sim$samples[[i]], ordered)
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
  ee <- extract_effects(x$formula, family = family)
  change <- list()
  standata <- standata(x)
  
  # find positions of parameters and define new names
  f <- colnames(standata$X)
  if (length(f) && !is.categorical(family)) {
    change <- lc(change, list(pos = grepl("^b\\[", pars), oldname = "b", 
                              pnames = paste0("b_",f), fnames = paste0("b_",f)))
    change <- c(change, prior_changes(class = "b", pars = pars, names = f))
    change <- c(change, prior_changes(class = "temp_Intercept", pars = pars, 
                                      new_class = "b_Intercept"))
  }
  
  if (is.formula(x$partial) || is.categorical(family)) {
    p <- colnames(standata$Xp)
    lp <- length(p)
    if (lp) {
      thres <- max(standata$max_obs) - 1
      pfnames <- paste0("b_",t(outer(p, paste0("[",1:thres,"]"), FUN = paste0)))
      change <- lc(change, list(pos = grepl("^bp\\[", pars), oldname = "bp", 
                                pnames = paste0("b_",p), fnames = pfnames,
                                sort = ulapply(1:lp, seq, to = thres*lp, by = lp),
                                dim = thres))
      change <- c(change, prior_changes(class = "bp", pars = pars, names = p))
    }
  } 
  
  if (length(x$ranef)) {
    group <- names(x$ranef)
    gf <- make_group_frame(x$ranef)
    for (i in seq_along(x$ranef)) {
      rfnames <- paste0("sd_", group[i],"_", x$ranef[[i]])
      change <- lc(change, list(pos = grepl(paste0("^sd_",i,"(\\[|$)"), pars),
                                oldname = paste0("sd_",i), pnames = rfnames, 
                                fnames = rfnames))
      change <- c(change, prior_changes(class = paste0("sd_",i), pars = pars, 
                                        names = x$ranef[[i]], 
                                        new_class = paste0("sd_",group[i])))
      
      if (length(x$ranef[[i]]) > 1 && ee$random$cor[[i]]) {
        cor_names <- get_cornames(x$ranef[[i]], type = paste0("cor_",group[i]), 
                                  brackets = FALSE)
        change <- lc(change, list(pos = grepl(paste0("^cor_",i,"(\\[|$)"), pars),
                                  oldname = paste0("cor_",i), pnames = cor_names,
                                  fnames = cor_names)) 
        change <- c(change, prior_changes(class = paste0("cor_",i), pars = pars, 
                                          new_class = paste0("cor_",group[i])))
      }
      if (any(grepl("^r_", pars))) {
        if (length(x$ranef[[i]]) == 1 || ee$random$cor[[i]]) {
          change <- lc(change, 
            ranef_changes(i = i, ranef = x$ranef, gf = gf, pars = pars,
                          dims_oi = x$fit@sim$dims_oi))
        } else {
          # multiple uncorrelated random effects
          for (j in seq_along(x$ranef[[i]])) {
            change <- lc(change, 
              ranef_changes(i = i, ranef = x$ranef, gf = gf, pars = pars,
                            dims_oi = x$fit@sim$dims_oi, j = j))
          }
        }
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
  # combines duplicated grouping factors to appear as if it was only one
  if (length(x$ranef)) 
    x$ranef <- combine_duplicates(x$ranef)
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
  #     first: a number corresponding to the first column for this term in the final r_<gf> matrices
  #     last: a number corresponding to the last column for this term in the final r_<gf> matrices
  group <- names(ranef)
  out <- data.frame(g = group, first = NA, last = NA)
  out[1, 2:3] <- c(1, length(ranef[[1]]))
  if (length(group) > 1) {
    for (i in 2:length(group)) {
      matches <- which(out$g[1:(i-1)] == group[i])
      if (length(matches))
        out[i, 2:3] <- c(out$last[max(matches)] + 1, 
                        out$last[max(matches)] + length(ranef[[i]]))
      else out[i, 2:3] <- c(1, length(ranef[[i]]))
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

combine_duplicates <- function(x) {
  # combine elements of a list that have the same name
  #
  # Args:
  #   x: a list
  #
  # Returns: 
  #   a list of possibly reduced length.
  # 
  # Examples:
  #   combine_duplicates(list(a = 1, a = c(2,3)))
  #   becomes list(a = c(1,2,3)) 
  if (!is.list(x)) 
    stop("x must be a list")
  if (is.null(names(x))) 
    stop("elements of x must be named")
  unique_names <- unique(names(x))
  new_list <- do.call(list, as.list(rep(NA, length(unique_names))))
  names(new_list) <- unique_names
  for (i in 1:length(unique_names)) {
    pos <- which(names(x) %in% unique_names[i])
    new_list[[unique_names[i]]] <- unname(unlist(x[pos]))
    # keep levels attribute
    attr(new_list[[unique_names[i]]], "levels") <- attr(x[[pos[1]]], "levels")
  }
  new_list
}

ranef_changes <- function(i, ranef, gf, dims_oi, pars, j = NULL)  {
  # helps in renaming random effects (r_) parameters
  # Args:
  #  i: the index of the grouping factor under consideration
  #  group: a vector of names of all grouping factors
  #  gf: matrix as constructed by make_group_frame
  #  dims_oi: named list containing parameter dimensions
  #  j: secondary indices used for uncorrelated random effects 
  # Returns:
  #  a list that can be interpreted by rename_pars
  group <- names(ranef)
  stopifnot(length(j) <= 1)
  r_index <- ifelse(is.null(j), i, paste0(i, "_", j))
  r_parnames <- paste0("^r_", r_index,"(\\[|$)")
  change <- list(pos = grepl(r_parnames, pars),
                 oldname = paste0("r_", r_index))
  
  # prepare for removal of redundant parameters r_<i>
  # and for combining random effects into one paramater matrix
  # number of total REs for this grouping factor
  n_ranefs <- max(gf$last[which(gf$g == group[i])]) 
  old_dim <- dims_oi[[change$oldname]]
  if (match(gf$g[i], group) == i && (is.null(j) || j == 1)) {
    change$pnames <- paste0("r_",group[i])
    change$dim <- if (n_ranefs == 1) old_dim 
                  else c(old_dim[1], n_ranefs) 
  } 
  # define index names of new parameter names
  colnames <- ranef[[i]]
  if (!is.null(j)) colnames <- colnames[j]
  index_names <- make_index_names(rownames = attr(ranef[[i]], "levels"),
                                  colnames = colnames, dim = 2)
  change$fnames <- paste0("r_", group[i], index_names)
  change
}

prior_changes <- function(class, pars, names = NULL, new_class = class) {
  # helps in renaming prior parameters
  #
  # Args: 
  #   class: the class of the parameters for which prior names should be changed
  #   pars: all parameters in the model
  #   names: names to replace digits at the end of parameter names
  #   new_class: replacment of the orginal class name
  #
  # Return:
  #   a list whose elements can be interpreted by rename_pars
  change <- list()
  pos_priors <- which(grepl(paste0("^prior_",class,"(_|$)"), pars))
  if (length(pos_priors)) {
    priors <- gsub(paste0("^prior_",class), paste0("prior_",new_class), 
                   pars[pos_priors])
    digits <- sapply(priors, function(prior) {
      d <- regmatches(prior, gregexpr("_[[:digit:]]+$", prior))[[1]]
      if (length(d)) 
        as.numeric(substr(d, 2, nchar(d))) 
      else 0
    })
    if (sum(abs(digits)) > 0 && is.null(names)) 
      stop("argument names is missing")
    for (i in 1:length(priors)) {
      if (digits[i]) 
        priors[i] <- gsub("[[:digit:]]+$", names[digits[i]], priors[i])
      if (pars[pos_priors[i]] != priors[i])
        change <- lc(change, list(pos = pos_priors[i], 
                                  oldname = pars[pos_priors[i]],
                                  pnames = priors[i],
                                  fnames = priors[i]))
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