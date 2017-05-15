message_new_method <- function(method, version) {
  message(
    "The '", method, "' method has been revised in brms ", version, ".\n",
    "To get the old output, set argument 'old' to TRUE."
  )
  invisible(NULL)
}

old_fixef_brmsfit <- function(object, estimate = "mean", ...) {
  # fixef.brmsfit as implemented prior to brms 1.7.0
  stopifnot(is.brmsfit(object))
  pars <- parnames(object)
  fpars <- pars[grepl(fixef_pars(), pars)]
  out <- as.matrix(object, pars = fpars, exact_match = TRUE)
  fpars <- gsub(fixef_pars(), "", fpars)
  out <- do.call(cbind, 
    lapply(estimate, get_estimate, samples = out, ...)
  )
  rownames(out) <- fpars
  out
}

old_ranef_brmsfit <- function(object, estimate = c("mean", "median"), 
                              var = FALSE, ...) {
  # ranef.brmsfit as implemented prior to brms 1.7.0
  .ranef <- function(group, nlpar = "") {
    # get group-level effects of a grouping factor
    # Args:
    #   group: name of a grouping factor
    #   nlpar: name of a non-linear parameter
    take <- object$ranef$group == group & object$ranef$nlpar == nlpar
    rnames <- object$ranef[take, "coef"]
    usc_nlpar <- usc(usc(nlpar))
    rpars <- pars[grepl(paste0("^r_", group, usc_nlpar, "\\["), pars)]
    if (!length(rpars)) {
      return(NULL)
    }
    rdims <- object$fit@sim$dims_oi[[paste0("r_", group, usc_nlpar)]]
    if (is.na(rdims[2])) rdims[2] <- 1
    levels <- attr(object$ranef, "levels")[[group]]
    if (is.null(levels)) {
      # avoid error in dimnames if levels are NULL 
      # for backwards compatibility with brms < 0.5.0 
      levels <- seq_len(rdims[1])
    }
    rs <- posterior_samples(object, pars = rpars, exact_match = TRUE)
    rs_array <- array(dim = c(rdims[1], rdims[2], nrow(rs)))
    k <- 0
    for (j in seq_len(rdims[2])) {
      for (l in seq_len(rdims[1])) {
        k <- k + 1
        rs_array[l, j, ] <- rs[, k]
      }
    }
    out <- get_estimate(estimate, samples = rs_array, margin = 1:2, ...)
    colnames(out) <- rnames
    if (var) {
      Var <- array(
        dim = c(rep(rdims[2], 2), rdims[1]), 
        dimnames = list(rnames, rnames, seq_len(rdims[1]))
      )
      if (rdims[2] == 1L) {
        for (j in seq_len(rdims[1])) {
          Var[, , j] <- var(rs_array[j, 1, ]) 
        }
      } else {
        for (j in seq_len(rdims[1])) {
          Var[, , j] <- cov(t(rs_array[j, , ]))
        }
      }
      dimnames(Var)[[3]] <- levels
      attr(out, "var") <- Var
    }
    rownames(out) <- levels
    if (nchar(nlpar)) {
      attr(out, "nlpar") <- nlpar
    }
    return(out)
  }
  
  stopifnot(is.brmsfit(object))
  estimate <- match.arg(estimate)
  pars <- parnames(object)
  group_nlpar <- unique(object$ranef[, c("group", "nlpar")])
  ranef <- named_list(group_nlpar$group)
  for (i in seq_along(ranef)) {
    ranef[[i]] <- .ranef(
      group = group_nlpar$group[i], 
      nlpar = group_nlpar$nlpar[i]
    )
  }
  rmNULL(ranef)
}

old_coef_brmsfit <- function(object, estimate = c("mean", "median"), ...) {
  # coef.brmsfit as implemented prior to brms 1.7.0
  .coef <- function(ranef, fixef) {
    # helper function to combine group and population-level effects
    all_ranef_names <- unique(ulapply(ranef, colnames))
    fixef_no_digits <- get_matches("^[^\\[]+", rownames(fixef))
    miss_fixef <- setdiff(all_ranef_names, rownames(fixef))
    miss_fixef_no_digits <- get_matches("^[^\\[]+", miss_fixef)
    new_fixef <- named_list(miss_fixef)
    for (k in seq_along(miss_fixef)) {
      # digits occur in ordinal models with category specific effects
      match_fixef <- match(miss_fixef_no_digits[k], rownames(fixef))
      if (!is.na(match_fixef)) {
        new_fixef[[k]] <- fixef[match_fixef, 1]
      } else if (!miss_fixef[k] %in% fixef_no_digits) {
        new_fixef[[k]] <- 0
      }
    }
    rm_fixef <- rownames(fixef) %in% miss_fixef_no_digits
    fixef <- fixef[!rm_fixef, , drop = FALSE]
    fixef <- do.call(rbind, c(list(fixef), rmNULL(new_fixef)))
    coef <- ranef
    for (i in seq_along(coef)) {
      ranef_names <- colnames(coef[[i]])
      ranef_no_digits <- get_matches("^[^\\[]+", ranef_names)
      miss_ranef <- setdiff(rownames(fixef), ranef_names)
      miss_ranef_no_digits <- get_matches("^[^\\[]+", miss_ranef)
      new_ranef <- named_list(miss_ranef)
      for (k in seq_along(miss_ranef)) {
        # digits occur in ordinal models with category specific effects
        match_ranef <- match(miss_ranef_no_digits[k], ranef_names)
        if (!is.na(match_ranef)) {
          new_ranef[[k]] <- coef[[i]][, match_ranef]
        } else if (!miss_ranef[k] %in% ranef_no_digits) {
          new_ranef[[k]] <- 0
        }
      }
      rm_ranef <- ranef_names %in% miss_ranef_no_digits
      coef[[i]] <- coef[[i]][, !rm_ranef, drop = FALSE]
      coef[[i]] <- do.call(cbind, c(list(coef[[i]]), rmNULL(new_ranef)))
      for (nm in colnames(coef[[i]])) {
        # correct the sign of thresholds in ordinal models
        sign <- ifelse(grepl("^Intercept\\[[[:digit:]]+\\]$", nm), -1, 1)
        coef[[i]][, nm] <- coef[[i]][, nm] + sign * fixef[nm, 1]
      }
    }
    return(coef)
  }
  
  stopifnot(is.brmsfit(object))
  estimate <- match.arg(estimate)
  fixef <- fixef(object, estimate = estimate, old = TRUE, ...)
  ranef <- ranef(object, estimate = estimate, old = TRUE, ...)
  nlpars <- ulapply(ranef, attr, "nlpar")
  if (length(nlpars)) {
    # do not combine effects of different nlpars
    unique_nlpars <- unique(nlpars)
    coef <- named_list(unique_nlpars)
    for (p in unique_nlpars) {
      ranef_temp <- ranef[nlpars %in% p]
      rx <- paste0("^", p, "_")
      take_rows <- grepl(rx, rownames(fixef))
      fixef_temp <- fixef[take_rows, , drop = FALSE]
      rownames(fixef_temp) <- sub(rx, "", rownames(fixef_temp))
      coef[[p]] <- .coef(ranef_temp, fixef_temp)
      for (i in seq_along(coef[[p]])) {
        attr(coef[[p]][[i]], "nlpar") <- p
      }
    }
    coef <- unlist(unname(coef), recursive = FALSE)
  } else {
    coef <- .coef(ranef, fixef)
  }
  coef
}

old_VarCorr_brmsfit <- function(x, estimate = "mean", ...) {
  # VarCorr.brmsfit as implemented prior to brms 1.7.0
  extract <- function(p) {
    # extracts samples for sd, cor and cov
    nr <- length(p$sd_pars)
    sd <- posterior_samples(x, pars = p$sd_pars, exact_match = TRUE)
    nsamples <- nrow(sd)
    out <- list(
      sd = do.call(cbind, lapply(estimate, get_estimate, samples = sd, ...))
    )
    rownames(out$sd) <- p$rnames 
    # calculate correlation and covariance matrices
    found_cor_pars <- intersect(p$cor_pars, parnames(x))
    if (length(found_cor_pars)) {
      cor <- posterior_samples(x, pars = paste0("^", found_cor_pars, "$"))
      if (length(found_cor_pars) < length(p$cor_pars)) { 
        # some correlations are missing and will be replaced by 0
        cor_all <- as.data.frame(
          matrix(0, nrow = nrow(cor), ncol = length(p$cor_pars))
        )
        names(cor_all) <- p$cor_pars
        for (i in seq_len(ncol(cor_all))) {
          found <- match(names(cor_all)[i], names(cor))
          if (!is.na(found)) {
            # correlation was estimated
            cor_all[, i] <- cor[, found]
          }
        }
        cor <- cor_all
      }
    } else {
      cor <- NULL
    } 
    # get_cov_matrix and array2list can be found in brmsfit-helpers.R
    matrices <- get_cov_matrix(sd = sd, cor = cor) 
    out$cor <- abind(lapply(
      estimate, get_estimate, samples = matrices$cor, 
      margin = c(2, 3), to.array = TRUE, ...
    ))
    out$cov <- abind(lapply(
      estimate, get_estimate, samples = matrices$cov, 
      margin = c(2,3), to.array = TRUE, ...
    )) 
    if (length(p$rnames) > 1) {
      dimnames(out$cor) <- list(p$rnames, p$rnames, dimnames(out$cor)[[3]])
      dimnames(out$cov) <- dimnames(out$cor)   
    }
    out$cor <- lapply(array2list(out$cor), function(x)
      if (is.null(dim(x))) 
        structure(matrix(x), dimnames = list(p$rnames, p$rnames)) 
      else x
    )
    out$cov <- lapply(array2list(out$cov), function(x)
      if (is.null(dim(x))) 
        structure(matrix(x), dimnames = list(p$rnames, p$rnames)) 
      else x
    )
    out
  }
  
  stopifnot(is.brmsfit(x))
  family <- family(x)
  bterms <- parse_bf(x$formula, family = family)
  if (nrow(x$ranef)) {
    get_names <- function(group) {
      # get names of group-level parameters
      r <- x$ranef[x$ranef$group == group, ]
      rnames <- paste0(usc(r$nlpar, "suffix"), r$coef)
      cor_type <- paste0("cor_", group)
      sd_pars <- paste0("sd_", group, "__", rnames)
      cor_pars <- get_cornames(rnames, type = cor_type, brackets = FALSE)
      nlist(rnames, type = cor_type, sd_pars, cor_pars)
    }
    group <- unique(x$ranef$group)
    p <- lapply(group, get_names)
  } else {
    p <- group <- NULL
  } 
  # special treatment of residuals variances in linear models
  has_sigma <- has_sigma(family, bterms, incmv = TRUE)
  if (has_sigma && !"sigma" %in% names(bterms$auxpars)) {
    cor_pars <- get_cornames(
      bterms$response, type = "rescor", brackets = FALSE
    )
    p <- lc(p, 
      list(rnames = bterms$response, cor_pars = cor_pars,
           sd_pars = c("sigma", paste0("sigma_", bterms$response)))
    )
    group <- c(group, "RESIDUAL")
  } 
  VarCorr <- lapply(p, extract)
  names(VarCorr) <- group
  class(VarCorr) <- "brmsVarCorr"
  VarCorr
}

#' @rdname VarCorr.brmsfit
#' @export
as.data.frame.brmsVarCorr <- function(x, ...) {
  estimates <- colnames(x[[1]]$sd)
  groups <- names(x)
  n_groups <- length(groups)
  names_coef <- lapply(x, function(y) rownames(y$sd))
  groups_col <- ulapply(1:n_groups, function(i) 
    c(groups[i], rep("", length(names_coef[[i]]) - 1)))
  max_cor <- max(ulapply(names_coef, length)) - 1
  # basic data.frame to be used in fill_base_frame
  base_frame <- as.data.frame(matrix(NA, nrow = length(groups_col),
                                     ncol = 4 + 2 * max_cor))
  names(base_frame) <- c("Group", "Name", "Std.Dev", rep("Cor", max_cor),
                         rep("Cov", max_cor + 1))
  base_frame[, 1:2] <- cbind(groups_col, unlist(names_coef))
  
  fill_base_frame <- function(estimate) {
    # fills the base_frame with SD and COR estimates
    # Args:
    #   estimate: The estimate being applied on the SD and COR parameters
    out <- base_frame
    pos <- 1
    for (i in 1:n_groups) {
      len <- length(names_coef[[i]])
      rows <- pos:(pos + len - 1)
      out[rows, "Std.Dev"] <- x[[i]]$sd[, estimate]
      if (len > 1) {
        # covariances and correlations present; add correlations
        cor_pos <- 4:(2 + len)
        cormat <- x[[i]]$cor[[estimate]][2:len, 1:(len-1), drop = FALSE]
        lt <- lower.tri(cormat, diag = TRUE)
        out[rows[2:length(rows)], cor_pos][lt] <- cormat[lt]
      }
      # add covariances
      cov_pos <- (4 + max_cor):(3 + max_cor + len)
      covmat <- x[[i]]$cov[[estimate]]
      lt <- lower.tri(covmat, diag = TRUE)
      out[rows, cov_pos][lt] <- covmat[lt]
      pos <- pos + len
    }
    out
  }
  
  out <- do.call(rbind, lapply(estimates, fill_base_frame))
  estimates_col <- ulapply(estimates, function(e)
    c(e, rep("", length(groups_col) - 1)))
  cbind(Estimate = estimates_col, out)
}

#' @export
print.brmsVarCorr <- function(x, digits = 2, ...) {
  dat <- as.data.frame(x)
  dat[, 4:ncol(dat)] <- round(as.matrix(dat[, 4:ncol(dat)]), digits = digits)
  dat[is.na(dat)] <- ""
  print(dat, row.names = FALSE, ...)
  invisible(x)
}
