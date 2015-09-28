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
  if (is.null(symbols))
    symbols <- c(" ", "(", ")", "[", "]", ",", "+", "-", "*", "/", "^", "=", "!=")
  if (is.null(subs))
    subs <- c(rep("", 6), "P", "M", "MU", "D", "E", "EQ", "NEQ")
  if (length(symbols) != length(subs)) 
    stop("length(symbols) != length(subs)")
  new_names <- names
  for (i in 1:length(symbols)) {
    new_names <- gsub(symbols[i], subs[i], new_names, fixed = fixed)
  }
  dup <- duplicated(new_names)
  if (check_dup && any(dup)) 
    stop(paste0("Internal renaming of variables led to duplicated names. \n",
                "Occured for variables: ", 
                paste(names[which(new_names %in% new_names[dup])], collapse = ", ")))
  new_names
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
  all_class <- c("b_Intercept", "b", "bp", "ar", "ma", "sd", "cor", 
                 "sigma", "rescor", "nu", "shape", "delta", "r", "prior", "lp")
  class <- regmatches(x$fit@sim$fnames_oi, regexpr("^[^_\\[]+", x$fit@sim$fnames_oi))
  # make sure that the fixed effects intercept comes first
  pos_intercept <- which(grepl("^b_Intercept($|\\[)", x$fit@sim$fnames_oi))
  class[pos_intercept] <- "b_Intercept"
  ordered <- order(factor(class, levels = all_class))
  x$fit@sim$fnames_oi <- x$fit@sim$fnames_oi[ordered]
  for (i in 1:chains) {
    # keep_attr ensures that attributes are not removed
    x$fit@sim$samples[[i]] <- keep_attr(x$fit@sim$samples[[i]], ordered)
  }
  mclass <- regmatches(x$fit@sim$pars_oi, regexpr("^[^_]+", x$fit@sim$pars_oi))
  # make sure that the fixed effects intercept comes first
  pos_intercept <- which(grepl("^b_Intercept($|\\[)", x$fit@sim$pars_oi))
  mclass[pos_intercept] <- "b_Intercept"
  ordered <- order(factor(mclass, levels = all_class))
  x$fit@sim$dims_oi <- x$fit@sim$dims_oi[ordered]
  x$fit@sim$pars_oi <- names(x$fit@sim$dims_oi)
  
  # some variables generally needed
  pars <- parnames(x)
  ee <- extract_effects(x$formula, family = x$family)
  change <- list()
  standata <- standata(x)
  
  # find positions of parameters and define new names
  f <- colnames(standata$X)
  if (length(f) && x$family != "categorical") {
    change <- lc(change, list(pos = grepl("^b\\[", pars), oldname = "b", 
                              pnames = paste0("b_",f), fnames = paste0("b_",f)))
    change <- c(change, prior_names(class = "b", pars = pars, names = f))
  }
  
  if (is.formula(x$partial) || x$family == "categorical") {
    if (x$family == "categorical") {
      p <- colnames(standata$X)
    } else {
      p <- colnames(standata$Xp)
    }
    lp <- length(p)
    thres <- max(standata$max_obs) - 1
    pfnames <- paste0("b_",t(outer(p, paste0("[",1:thres,"]"), FUN = paste0)))
    change <- lc(change, list(pos = grepl("^bp\\[", pars), oldname = "bp", 
                              pnames = paste0("b_",p), fnames = pfnames,
                              sort = ulapply(1:lp, seq, to = thres*lp, by = lp),
                              dim = thres))
    change <- c(change, prior_names(class = "bp", pars = pars, 
                                    names = p, new_class = "b"))
  }  
  
  if (length(x$ranef)) {
    group <- names(x$ranef)
    gf <- make_group_frame(x$ranef)
    for (i in 1:length(x$ranef)) {
      rfnames <- paste0("sd_",group[i],"_", x$ranef[[i]])
      change <- lc(change, list(pos = grepl(paste0("^sd_",i,"(\\[|$)"), pars),
                                oldname = paste0("sd_",i), pnames = rfnames, 
                                fnames = rfnames))
      change <- c(change, prior_names(class = paste0("sd_",i), pars = pars, 
                                      names = x$ranef[[i]], 
                                      new_class = paste0("sd_",group[i])))
      
      if (length(x$ranef[[i]]) > 1 && ee$cor[[i]]) {
        cor_names <- get_cornames(x$ranef[[i]], type = paste0("cor_",group[i]), 
                                  brackets = FALSE)
        change <- lc(change, list(pos = grepl(paste0("^cor_",i,"(\\[|$)"), pars),
                                  oldname = paste0("cor_",i), pnames = cor_names,
                                  fnames = cor_names)) 
        change <- c(change, prior_names(class = paste0("cor_",i), pars = pars, 
                                        new_class = paste0("cor_",group[i])))
      }
      if (any(grepl("^r_", pars))) {
        lc <- length(change) + 1
        change[[lc]] <- list(pos = grepl(paste0("^r_",i,"(\\[|$)"), pars),
                             oldname = paste0("r_",i))
        
        # prepare for removal of redundant parameters r_<i>
        # and for combining random effects into one paramater matrix
        # number of total REs for this grouping factor
        n_ranefs <- max(gf$last[which(gf$g == group[i])]) 
        old_dim <- x$fit@sim$dims_oi[[change[[lc]]$oldname]]
        indices <- make_indices(rows = 1:old_dim[1], cols = gf$first[i]:gf$last[i], 
                                dim = ifelse(n_ranefs == 1, 1, 2))
        if (match(gf$g[i], group) < i) 
          change[[lc]]$pnames <- NULL 
        else {
          change[[lc]]$pnames <- paste0("r_",group[i])
          change[[lc]]$dim <- if (n_ranefs == 1) old_dim 
                              else c(old_dim[1], n_ranefs) 
        } 
        change[[lc]]$fnames <- paste0("r_",group[i], indices)
      }  
    }
  }
  if (x$family %in% c("gaussian", "student", "cauchy") && !is.formula(ee$se)) {
    corfnames <- paste0("sigma_",ee$response)
    change <- lc(change, list(pos = grepl("^sigma", pars), oldname = "sigma",
                              pnames = corfnames, fnames = corfnames))
    change <- c(change, prior_names(class = "sigma", pars = pars, 
                                     names = ee$response))
    # residual correlation paramaters
    if (x$family == "gaussian" && length(ee$response) > 1) {
       rescor_names <- get_cornames(ee$response, type = "rescor", brackets = FALSE)
       change <- lc(change, list(pos = grepl("^rescor\\[", pars), oldname = "rescor",
                                 pnames = rescor_names, fnames = rescor_names))
    }
  } 
  
  # perform the actual renaming in x$fit@sim
  if (length(change)) {
    for (c in 1:length(change)) {
      x$fit@sim$fnames_oi[change[[c]]$pos] <- change[[c]]$fnames
      for (i in 1:chains) {
        names(x$fit@sim$samples[[i]])[change[[c]]$pos] <- change[[c]]$fnames
        if (!is.null(change[[c]]$sort)) {
          x$fit@sim$samples[[i]][change[[c]]$pos] <- 
            x$fit@sim$samples[[i]][change[[c]]$pos][change[[c]]$sort]
        }
      }
      onp <- match(change[[c]]$oldname, names(x$fit@sim$dims_oi))
      if (is.null(change[[c]]$pnames)) {
        # only needed to collapse multiple r_<i> of the same grouping factor
        x$fit@sim$dims_oi[[onp]] <- NULL  
      } else { 
        # rename dims_oi to match names in fnames_oi
        dims <- x$fit@sim$dims_oi
        x$fit@sim$dims_oi <- c(if (onp > 1) dims[1:(onp - 1)], 
                               make_dims(change[[c]]),
                               dims[(onp + 1):length(dims)])
      }
    }
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

make_indices <- function(rows, cols = NULL, dim = 1) {
  # make indices in square brackets for indexing stan parameters
  #
  # Args:
  #   rows: a vector of rows
  #   cols: a vector of columns
  #   dim: The number of dimensions of the output either 1 or 2
  #
  # Returns:
  #   all index pairs of rows and cols
  if (!dim %in% c(1,2))
    stop("dim must be 1 or 2")
  if (dim == 1) 
    indices <- paste0("[",rows,"]")
  else {
    indices <- paste0("[", outer(rows, cols, FUN = paste, sep = ","), "]")
  }
  indices
}

combine_duplicates <- function(x) {
  # combine elements of a list that have the same name
  #
  # Args:
  #   x: a list
  #
  # Returns: 
  #   a list of possibly reducte length.
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

prior_names <- function(class, pars, names = NULL, new_class = class) {
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