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
  bterms <- parse_bf(x$formula)
  data <- model.frame(x)
  meef <- tidy_meef(bterms, data)
  pars <- parnames(x)
  # find positions of parameters and define new names
  change <- c(
    change_effects(bterms, data = data, pars = pars, scode = stancode(x)),
    change_re(x$ranef, pars = pars),
    change_Xme(meef, pars = pars),
    change_autocor(bterms, data = data, pars = pars)
  )
  # perform the actual renaming in x$fit@sim
  x <- do_renaming(x, change)
  x <- compute_quantities(x)
  x <- reorder_pars(x)
  x
}

change_effects <- function(x, ...) {
  # helps in renaming parameters after model fitting
  UseMethod("change_effects")
}

#' @export
change_effects.mvbrmsterms <- function(x, pars, ...) {
  change <- list()
  for (i in seq_along(x$terms)) {
    change <- c(change, change_effects(x$terms[[i]], pars = pars, ...))
  }
  if (x$rescor) {
    rescor_names <- get_cornames(
      x$responses, type = "rescor", brackets = FALSE
    )
    change <- lc(change,
      list(pos = grepl("^rescor\\[", pars), fnames = rescor_names)
    )
  }
  change
}

#' @export
change_effects.brmsterms <- function(x, ...) {
  change <- list()
  for (dp in names(x$dpars)) {
    change <- c(change, change_effects(x$dpars[[dp]], ...))
  }
  if (is.formula(x$adforms$mi)) {
    change <- c(change, change_Ymi(x, ...))
  }
  change
}

#' @export
change_effects.btl <- function(x, data, pars, scode = "", ...) {
  # helps in renaming various kinds of effects
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  c(
    change_fe(x, data, pars, scode = scode),
    change_sm(x, data, pars),
    change_cs(x, data, pars),
    change_sp(x, data, pars),
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

change_fe <- function(bterms, data, pars, scode = "") {
  # helps in renaming fixed effects parameters
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  change <- list()
  px <- check_prefix(bterms)
  fixef <- colnames(data_fe(bterms, data)$X)
  fixef <- rm_int_fe(fixef, scode, px = px)
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

change_sp <- function(bterms, data, pars) {
  # helps in renaming special effects parameters
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  change <- list()
  spef <- tidy_spef(bterms, data)
  if (is.null(spef)) {
    return(change) 
  }
  p <- usc(combine_prefix(bterms))
  bsp <- paste0("bsp", p)
  pos <- grepl(paste0("^", bsp, "\\["), pars)
  newnames <- paste0("bsp", p, "_", spef$coef)
  change <- lc(change, nlist(pos, fnames = newnames))
  change <- c(change, 
    change_prior(class = bsp, pars = pars, names = spef$coef)
  )
  simo_coef <- get_simo_labels(spef)
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

change_Xme <- function(meef, pars) {
  # helps in renaming global noise free variables
  # Returns:
  #   a list whose elements can be interpreted by do_renaming
  stopifnot(is.meef_frame(meef))
  change <- list()
  levels <- attr(meef, "levels")
  groups <- unique(meef$grname)
  for (i in seq_along(groups)) {
    g <- groups[i]
    K <- which(meef$grname %in% g)
    # rename mean and sd parameters
    for (par in c("meanme", "sdme")) {
      hpar <- paste0(par, "_", i)
      pos <- grepl(paste0("^", hpar, "\\["), pars)
      hpar_new <- paste0(par, "_", meef$coef[K])
      change <- lc(change, nlist(pos, fnames = hpar_new))
      change <- c(change,
        change_prior(class = hpar, pars = pars, names = hpar_new)
      )
    }
    # rename latent variable parameters
    for (k in K) {
      if (any(grepl("^Xme_", pars))) {
        Xme <- paste0("Xme_", k)
        pos <- grepl(paste0("^", Xme, "\\["), pars)
        Xme_new <- paste0("Xme_", meef$coef[k])
        if (nzchar(g)) {
          indices <- gsub("[ \t\r\n]", ".", levels[[g]])
        } else {
          indices <- seq_len(sum(pos))
        }
        fnames <- paste0(Xme_new, "[", indices, "]") 
        change <- lc(change, nlist(pos, fnames))
      }
    }
    # rename correlation parameters
    if (meef$cor[K[1]] && length(K) > 1L) {
      cor_type <- paste0("corme", usc(g))
      cor_names <- get_cornames(meef$coef[K], cor_type, brackets = FALSE)
      cor_regex <- paste0("^corme_", i, "(\\[|$)")
      cor_pos <- grepl(cor_regex, pars)
      change <- lc(change, list(pos = cor_pos, fnames = cor_names))
      change <- c(change,
        change_prior(
          class = paste0("corme_", i), pars = pars,
          new_class = paste0("corme", usc(g))
        )
      )
    }
  }
  change
}

change_Ymi <- function(bterms, data, pars, ...) {
  # helps in renaming estimated missing values
  stopifnot(is.brmsterms(bterms))
  change <- list()
  if (is.formula(bterms$adforms$mi)) {
    resp <- usc(combine_prefix(bterms))
    resp_data <- data_response(bterms, data, check_response = FALSE)
    Ymi <- paste0("Ymi", resp)
    pos <- grepl(paste0("^", Ymi, "\\["), pars)
    if (any(pos)) {
      Jmi <- resp_data$Jmi
      fnames <- paste0(Ymi, "[", Jmi, "]")
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
      rnames <- get_rnames(r)
      sd_names <- paste0("sd_", g, "__", as.vector(rnames))
      sd_pos <- grepl(paste0("^sd_", id, "(\\[|$)"), pars)
      change <- lc(change, list(pos = sd_pos, fnames = sd_names))
      change <- c(change,
        change_prior(
          class = paste0("sd_", id), pars = pars,
          new_class = paste0("sd_", g), 
          names = paste0("_", as.vector(rnames))
        )
      )
      # rename group-level correlations
      if (nrow(r) > 1L && isTRUE(r$cor[1])) {
        type <- paste0("cor_", g)
        if (isTRUE(nzchar(r$by[1]))) {
          cor_names <- named_list(r$bylevels[[1]])
          for (j in seq_len(nrow(rnames))) {
            cor_names[[j]] <- get_cornames(
              rnames[, j], type, brackets = FALSE
            )
          }
          cor_names <- unlist(cor_names)
        } else {
          cor_names <- get_cornames(rnames, type, brackets = FALSE)
        }
        cor_regex <- paste0("^cor_", id, "(_[[:digit:]]+)?(\\[|$)")
        cor_pos <- grepl(cor_regex, pars)
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
    if (!is.null(bterms$time$group)) {
      group <- gsub("[ \t\r\n]", "", get(bterms$time$group, data))
    } else {
      group <- rep(1, nrow(data)) 
    }
    if (!is.null(bterms$time$time)) {
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

rm_int_fe <- function(fixef, scode, px = list()) {
  # identifies if the intercept has to be removed from fixef
  # and returns adjusted fixef names
  p <- usc(combine_prefix(px))
  regex <- paste0("(temp", p, "_Intercept)|(vector\\[N\\] loclev", p, ";)")
  if (grepl(regex, scode)) {
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
    "b", "bsp", "bcs", "ar", "ma", "arr", "lagsar",
    "errorsar", "car", "sdcar", "sigmaLL", "sd", "cor", "sds", 
    "sdgp", "lscale", dpars(), "temp", "rescor", "delta", 
    "lasso", "simo", "r", "s", "zgp", "rcar", "loclev", 
    "Ymi", "Yl", "meanme", "sdme", "corme", "Xme", "prior", "lp"
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
  has_subpars <- nsubpars > 0
  new_order <- new_order[has_subpars]
  nsubpars <- nsubpars[has_subpars]
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
  x <- compute_xi(x)
  x
}

compute_xi <- function(x, ...) {
  # helper function to compute parameter xi, which is currently
  # defined in the Stan model block and thus not being stored
  # Returns:
  #   All methods return a brmsfit object
  UseMethod("compute_xi")
}

#' @export
compute_xi.brmsfit <- function(x, ...) {
  if (!any(grepl("^temp_xi(_|$)", parnames(x)))) {
    return(x)
  }
  draws <- try(extract_draws(x))
  if (is(draws, "try-error")) {
    warning2("Trying to compute 'xi' was unsuccessful. ",
             "Some S3 methods may not work as expected.")
    return(x)
  }
  compute_xi(draws, fit = x, ...)
}

#' @export
compute_xi.mvbrmsdraws <- function(x, fit, ...) {
  stopifnot(is.brmsfit(fit))
  for (resp in names(x$resps)) {
    fit <- compute_xi(x$resps[[resp]], fit = fit, ...)
  }
  fit
}

#' @export
compute_xi.brmsdraws <- function(x, fit, ...) {
  stopifnot(is.brmsfit(fit))
  resp <- usc(x$resp)
  temp_xi_name <- paste0("temp_xi", resp)
  if (!temp_xi_name %in% parnames(fit)) {
    return(fit)
  }
  mu <- get_dpar(x, "mu")
  sigma <- get_dpar(x, "sigma")
  y <- matrix(x$data$Y, dim(mu)[1], dim(mu)[2], byrow = TRUE)
  bs <- - 1 / matrixStats::rowRanges((y - mu) / sigma)
  bs <- matrixStats::rowRanges(bs)
  temp_xi <- as.vector(as.matrix(fit, pars = temp_xi_name))
  xi <- inv_logit(temp_xi) * (bs[, 2] - bs[, 1]) + bs[, 1]
  # write xi into stanfit object
  xi_name <- paste0("xi", resp)
  samp_chain <- length(xi) / fit$fit@sim$chains
  for (i in seq_len(fit$fit@sim$chains)) {
    xi_part <- xi[((i - 1) * samp_chain + 1):(i * samp_chain)]
    # add warmup samples not used anyway
    xi_part <- c(rep(0, fit$fit@sim$warmup2[1]), xi_part)
    fit$fit@sim$samples[[i]][[xi_name]] <- xi_part
  }
  fit$fit@sim$pars_oi <- c(fit$fit@sim$pars_oi, xi_name)
  fit$fit@sim$dims_oi[[xi_name]] <- numeric(0)
  fit$fit@sim$fnames_oi <- c(fit$fit@sim$fnames_oi, xi_name)
  fit$fit@sim$n_flatnames <- fit$fit@sim$n_flatnames + 1
  fit
}
