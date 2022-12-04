#' Rename Parameters
#'
#' Rename parameters within the \code{stanfit} object after model fitting to
#' ensure reasonable parameter names. This function is usually called
#' automatically by \code{\link{brm}} and users will rarely be required to call
#' it themselves.
#'
#' @param x A brmsfit object.
#' @return A brmfit object with adjusted parameter names.
#'
#' @examples
#' \dontrun{
#' # fit a model manually via rstan
#' scode <- make_stancode(count ~ Trt, data = epilepsy)
#' sdata <- make_standata(count ~ Trt, data = epilepsy)
#' stanfit <- rstan::stan(model_code = scode, data = sdata)
#'
#' # feed the Stan model back into brms
#' fit <- brm(count ~ Trt, data = epilepsy, empty = TRUE)
#' fit$fit <- stanfit
#' fit <- rename_pars(fit)
#' summary(fit)
#' }
#'
#' @export
rename_pars <- function(x) {
  if (!length(x$fit@sim)) {
    return(x)
  }
  bterms <- brmsterms(x$formula)
  data <- model.frame(x)
  meef <- tidy_meef(bterms, data)
  pars <- variables(x)
  # find positions of parameters and define new names
  change <- c(
    change_effects(bterms, data = data, pars = pars),
    change_re(x$ranef, pars = pars),
    change_Xme(meef, pars = pars)
  )
  # perform the actual renaming in x$fit@sim
  x <- save_old_par_order(x)
  x <- do_renaming(x, change)
  x$fit <- repair_stanfit(x$fit)
  x <- compute_quantities(x)
  x <- reorder_pars(x)
  x
}

# helps in renaming parameters after model fitting
# @return a list whose elements can be interpreted by do_renaming
change_effects <- function(x, ...) {
  UseMethod("change_effects")
}

#' @export
change_effects.default <- function(x, ...) {
  NULL
}

#' @export
change_effects.mvbrmsterms <- function(x, data, pars, ...) {
  out <- list()
  for (i in seq_along(x$terms)) {
    c(out) <- change_effects(x$terms[[i]], data = data, pars = pars, ...)
  }
  if (x$rescor) {
    rescor_names <- get_cornames(
      x$responses, type = "rescor", brackets = FALSE
    )
    lc(out) <- clist(grepl("^rescor\\[", pars), rescor_names)
  }
  out
}

#' @export
change_effects.brmsterms <- function(x, data, ...) {
  data <- subset_data(data, x)
  out <- list()
  for (dp in names(x$dpars)) {
    c(out) <- change_effects(x$dpars[[dp]], data = data, ...)
  }
  for (nlp in names(x$nlpars)) {
    c(out) <- change_effects(x$nlpars[[nlp]], data = data, ...)
  }
  if (is.formula(x$adforms$mi)) {
    c(out) <- change_Ymi(x, data = data, ...)
  }
  c(out) <- change_thres(x, data = data, ...)
  c(out) <- change_family_cor_pars(x, data = data, ...)
  out
}

# helps in renaming parameters of additive predictor terms
# @param pars vector of all parameter names
#' @export
change_effects.btl <- function(x, ...) {
  c(change_fe(x, ...),
    change_sm(x, ...),
    change_cs(x, ...),
    change_sp(x, ...),
    change_gp(x, ...))
}

# helps in renaming fixed effects parameters
change_fe <- function(bterms, data, pars, ...) {
  out <- list()
  px <- check_prefix(bterms)
  fixef <- colnames(data_fe(bterms, data)$X)
  if (stan_center_X(bterms)) {
    fixef <- setdiff(fixef, "Intercept")
  }
  if (length(fixef)) {
    b <- paste0("b", usc(combine_prefix(px)))
    pos <- grepl(paste0("^", b, "\\["), pars)
    bnames <- paste0(b, "_", fixef)
    lc(out) <- clist(pos, bnames)
    c(out) <- change_prior(b, pars, names = fixef)
    c(out) <- change_special_prior_local(bterms, fixef, pars)
  }
  out
}

# helps in renaming special effects parameters
change_sp <- function(bterms, data, pars, ...) {
  out <- list()
  spef <- tidy_spef(bterms, data)
  if (!nrow(spef)) return(out)
  p <- usc(combine_prefix(bterms))
  bsp <- paste0("bsp", p)
  pos <- grepl(paste0("^", bsp, "\\["), pars)
  newnames <- paste0("bsp", p, "_", spef$coef)
  lc(out) <- clist(pos, newnames)
  c(out) <- change_prior(bsp, pars, names = spef$coef)
  simo_coef <- get_simo_labels(spef)
  for (i in seq_along(simo_coef)) {
    simo_old <- paste0("simo", p, "_", i)
    simo_new <- paste0("simo", p, "_", simo_coef[i])
    pos <- grepl(paste0("^", simo_old, "\\["), pars)
    simo_names <- paste0(simo_new, "[", seq_len(sum(pos)), "]")
    lc(out) <- clist(pos, simo_names)
    c(out) <- change_prior(
      simo_old, pars, new_class = simo_new, is_vector = TRUE
    )
  }
  out
}

# helps in renaming category specific effects parameters
change_cs <- function(bterms, data, pars, ...) {
  out <- list()
  csef <- colnames(data_cs(bterms, data)$Xcs)
  if (length(csef)) {
    p <- usc(combine_prefix(bterms))
    bcsp <- paste0("bcs", p)
    ncs <- length(csef)
    thres <- sum(grepl(paste0("^b", p, "_Intercept\\["), pars))
    csenames <- t(outer(csef, paste0("[", 1:thres, "]"), FUN = paste0))
    csenames <- paste0(bcsp, "_", csenames)
    sort_cse <- ulapply(seq_len(ncs), seq, to = thres * ncs, by = ncs)
    lc(out) <- clist(
      grepl(paste0("^", bcsp, "\\["), pars), csenames, sort = sort_cse
    )
    c(out) <- change_prior(bcsp, pars, names = csef)
  }
  out
}

# rename threshold parameters in ordinal models
change_thres <- function(bterms, pars, ...) {
  out <- list()
  # renaming is only required if multiple threshold were estimated
  if (!has_thres_groups(bterms)) {
    return(out)
  }
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  int <- paste0("b", p, "_Intercept")
  groups <- get_thres_groups(bterms)
  for (i in seq_along(groups)) {
    thres <- get_thres(bterms, groups[i])
    pos <- grepl(glue("^{int}_{i}\\["), pars)
    int_names <- glue("{int}[{groups[i]},{thres}]")
    lc(out) <- clist(pos, int_names)
  }
  out
}

# helps in renaming global noise free variables
# @param meef data.frame returned by 'tidy_meef'
change_Xme <- function(meef, pars, ...) {
  stopifnot(is.meef_frame(meef))
  out <- list()
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
      lc(out) <- clist(pos, hpar_new)
      c(out) <- change_prior(hpar, pars, names = hpar_new)
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
        lc(out) <- clist(pos, fnames)
      }
    }
    # rename correlation parameters
    if (meef$cor[K[1]] && length(K) > 1L) {
      cor_type <- paste0("corme", usc(g))
      cor_names <- get_cornames(meef$coef[K], cor_type, brackets = FALSE)
      cor_regex <- paste0("^corme_", i, "(\\[|$)")
      cor_pos <- grepl(cor_regex, pars)
      lc(out) <- clist(cor_pos, cor_names)
      c(out) <- change_prior(
        paste0("corme_", i), pars, new_class = paste0("corme", usc(g))
      )
    }
  }
  out
}

# helps in renaming estimated missing values
change_Ymi <- function(bterms, data, pars, ...) {
  stopifnot(is.brmsterms(bterms))
  out <- list()
  if (is.formula(bterms$adforms$mi)) {
    resp <- usc(combine_prefix(bterms))
    resp_data <- data_response(bterms, data, check_response = FALSE)
    Ymi <- paste0("Ymi", resp)
    pos <- grepl(paste0("^", Ymi, "\\["), pars)
    if (any(pos)) {
      Jmi <- resp_data$Jmi
      fnames <- paste0(Ymi, "[", Jmi, "]")
      lc(out) <- clist(pos, fnames)
    }
  }
  out
}

# helps in renaming parameters of gaussian processes
change_gp <- function(bterms, data, pars, ...) {
  out <- list()
  p <- usc(combine_prefix(bterms), "prefix")
  gpef <- tidy_gpef(bterms, data)
  for (i in seq_rows(gpef)) {
    # rename GP hyperparameters
    sfx1 <- gpef$sfx1[[i]]
    sfx2 <- as.vector(gpef$sfx2[[i]])
    sdgp <- paste0("sdgp", p)
    sdgp_old <- paste0(sdgp, "_", i)
    sdgp_pos <- grepl(paste0("^", sdgp_old, "\\["), pars)
    sdgp_names <- paste0(sdgp, "_", sfx1)
    lc(out) <- clist(sdgp_pos, sdgp_names)
    c(out) <- change_prior(sdgp_old, pars, names = sfx1, new_class = sdgp)

    lscale <- paste0("lscale", p)
    lscale_old <- paste0(lscale, "_", i)
    lscale_pos <- grepl(paste0("^", lscale_old, "\\["), pars)
    lscale_names <- paste0(lscale, "_", sfx2)
    lc(out) <- clist(lscale_pos, lscale_names)
    c(out) <- change_prior(lscale_old, pars, names = sfx2, new_class = lscale)

    zgp <- paste0("zgp", p)
    zgp_old <- paste0(zgp, "_", i)
    if (length(sfx1) > 1L) {
      # categorical 'by' variable
      for (j in seq_along(sfx1)) {
        zgp_old_sub <- paste0(zgp_old, "_", j)
        zgp_pos <- grepl(paste0("^", zgp_old_sub, "\\["), pars)
        if (any(zgp_pos)) {
          zgp_new <- paste0(zgp, "_", sfx1[j])
          fnames <- paste0(zgp_new, "[", seq_len(sum(zgp_pos)), "]")
          lc(out) <- clist(zgp_pos, fnames)
        }
      }
    } else {
      zgp_pos <- grepl(paste0("^", zgp_old, "\\["), pars)
      if (any(zgp_pos)) {
        zgp_new <- paste0(zgp, "_", sfx1)
        fnames <- paste0(zgp_new, "[", seq_len(sum(zgp_pos)), "]")
        lc(out) <- clist(zgp_pos, fnames)
      }
    }
  }
  out
}

# helps in renaming smoothing term parameters
change_sm <- function(bterms, data, pars, ...) {
  out <- list()
  smef <- tidy_smef(bterms, data)
  if (NROW(smef)) {
    p <- usc(combine_prefix(bterms), "prefix")
    Xs_names <- attr(smef, "Xs_names")
    if (length(Xs_names)) {
      bs <- paste0("bs", p)
      pos <- grepl(paste0("^", bs, "\\["), pars)
      bsnames <- paste0(bs, "_", Xs_names)
      lc(out) <- clist(pos, bsnames)
      c(out) <- change_prior(bs, pars, names = Xs_names)
    }
    sds <- paste0("sds", p)
    sds_names <- paste0(sds, "_", smef$label)
    s <- paste0("s", p)
    snames <- paste0(s, "_", smef$label)
    for (i in seq_rows(smef)) {
      for (j in seq_len(smef$nbases[i])) {
        ij <- paste0(i, "_", j)
        sds_pos <- grepl(paste0("^", sds, "_", ij), pars)
        lc(out) <- clist(sds_pos, paste0(sds_names[i], "_", j))
        spos <- grepl(paste0("^", s, "_", ij), pars)
        sfnames <- paste0(snames[i], "_", j, "[", seq_len(sum(spos)), "]")
        lc(out) <- clist(spos, sfnames)
        new_prior_class <- paste0(sds, "_", smef$label[i], "_", j)
        c(out) <- change_prior(
          paste0(sds, "_", ij), pars, new_class = new_prior_class
        )
      }
    }
  }
  out
}

# helps in renaming group-level parameters
# @param ranef: data.frame returned by 'tidy_ranef'
change_re <- function(ranef, pars, ...) {
  out <- list()
  if (has_rows(ranef)) {
    for (id in unique(ranef$id)) {
      r <- subset2(ranef, id = id)
      g <- r$group[1]
      rnames <- get_rnames(r)
      sd_names <- paste0("sd_", g, "__", as.vector(rnames))
      sd_pos <- grepl(paste0("^sd_", id, "(\\[|$)"), pars)
      lc(out) <- clist(sd_pos, sd_names)
      c(out) <- change_prior(
        paste0("sd_", id), pars, new_class = paste0("sd_", g),
        names = paste0("_", as.vector(rnames))
      )
      # rename group-level correlations
      if (nrow(r) > 1L && isTRUE(r$cor[1])) {
        type <- paste0("cor_", g)
        if (isTRUE(nzchar(r$by[1]))) {
          cor_names <- named_list(r$bylevels[[1]])
          for (j in seq_along(cor_names)) {
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
        lc(out) <- clist(cor_pos, cor_names)
        c(out) <- change_prior(
          paste0("cor_", id), pars, new_class = paste0("cor_", g)
        )
      }
    }
    if (any(grepl("^r_", pars))) {
      c(out) <- change_re_levels(ranef, pars = pars)
    }
    tranef <- get_dist_groups(ranef, "student")
    for (i in seq_rows(tranef)) {
      df_pos <- grepl(paste0("^df_", tranef$ggn[i], "$"), pars)
      df_name <- paste0("df_", tranef$group[i])
      lc(out) <- clist(df_pos, df_name)
    }
  }
  out
}

# helps in renaming varying effects parameters per level
# @param ranef: data.frame returned by 'tidy_ranef'
change_re_levels <- function(ranef, pars, ...)  {
  out <- list()
  for (i in seq_rows(ranef)) {
    r <- ranef[i, ]
    p <- usc(combine_prefix(r))
    r_parnames <- paste0("r_", r$id, p, "_", r$cn)
    r_regex <- paste0("^", r_parnames, "(\\[|$)")
    r_new_parname <- paste0("r_", r$group, usc(p))
    # rstan doesn't like whitespaces in parameter names
    levels <- gsub("[ \t\r\n]", ".", attr(ranef, "levels")[[r$group]])
    index_names <- make_index_names(levels, r$coef, dim = 2)
    fnames <- paste0(r_new_parname, index_names)
    lc(out) <- clist(grepl(r_regex, pars), fnames)
  }
  out
}

# helps to rename correlation parameters of likelihoods
change_family_cor_pars <- function(x, pars, ...) {
  stopifnot(is.brmsterms(x))
  out <- list()
  if (is_logistic_normal(x$family)) {
    predcats <- get_predcats(x$family)
    lncor_names <- get_cornames(
      predcats, type = "lncor", brackets = FALSE
    )
    lc(out) <- clist(grepl("^lncor\\[", pars), lncor_names)
  }
  out
}

# rename parameters related to special priors
change_special_prior_local <- function(bterms, coef, pars, ...) {
  out <- list()
  p <- combine_prefix(bterms)
  # rename parameters related to the R2D2 prior
  pos_R2D2_phi <- grepl(paste0("^R2D2_phi", p), pars)
  if (any(pos_R2D2_phi)) {
    phi <- paste0("R2D2_phi", p)
    new_phi <- paste0(phi, "_", coef)
    lc(out) <- clist(pos_R2D2_phi, new_phi)
    c(out) <- change_prior(phi, pars, names = coef, is_vector = TRUE)
  }
  out
}

# helps in renaming prior parameters
# @param class the class of the parameters
# @param pars names of all parameters in the model
# @param names names to replace digits at the end of parameter names
# @param new_class optional replacement of the orginal class name
# @param is_vector indicate if the prior parameter is a vector
change_prior <- function(class, pars, names = NULL, new_class = class,
                         is_vector = FALSE) {
  out <- list()
  # 'stan_rngprior' adds '__' before the digits to disambiguate
  regex <- paste0("^prior_", class, "(__[[:digit:]]+|$|\\[)")
  pos_priors <- which(grepl(regex, pars))
  if (length(pos_priors)) {
    priors <- gsub(
      paste0("^prior_", class),
      paste0("prior_", new_class),
      pars[pos_priors]
    )
    if (is_vector) {
      if (!is.null(names)) {
        .names <- paste0("_", names)
        for (i in seq_along(priors)) {
          priors[i] <- gsub("\\[[[:digit:]]+\\]$", .names[i], priors[i])
        }
      }
      lc(out) <- clist(pos_priors, priors)
    } else {
      digits <- sapply(priors, function(prior) {
        d <- regmatches(prior, gregexpr("__[[:digit:]]+$", prior))[[1]]
        if (length(d)) as.numeric(substr(d, 3, nchar(d))) else 0
      })
      for (i in seq_along(priors)) {
        if (digits[i] && !is.null(names)) {
          priors[i] <- gsub("_[[:digit:]]+$", names[digits[i]], priors[i])
        }
        if (pars[pos_priors[i]] != priors[i]) {
          lc(out) <- clist(pos_priors[i], priors[i])
        }
      }
    }
  }
  out
}

# helper for change_* functions
clist <- function(pos, fnames, ...) {
  structure(nlist(pos, fnames, ...), class = c("clist", "list"))
}

is.clist <- function(x) {
  inherits(x, "clist")
}

# compute index names in square brackets for indexing stan parameters
# @param rownames a vector of row names
# @param colnames a vector of columns
# @param dim the number of output dimensions
# @return all index pairs of rows and cols
make_index_names <- function(rownames, colnames = NULL, dim = 1) {
  if (!dim %in% c(1, 2))
    stop("dim must be 1 or 2")
  if (dim == 1) {
    index_names <- paste0("[", rownames, "]")
  } else {
    tmp <- outer(rownames, colnames, FUN = paste, sep = ",")
    index_names <- paste0("[", tmp, "]")
  }
  index_names
}

# save original order of the parameters in the stanfit object
save_old_par_order <- function(x) {
  x$fit@sim$pars_oi_old <- x$fit@sim$pars_oi
  x$fit@sim$dims_oi_old <- x$fit@sim$dims_oi
  x$fit@sim$fnames_oi_old <- x$fit@sim$fnames_oi
  x
}

# perform actual renaming of Stan parameters
# @param x a brmsfit object
# @param change a list of lists each element allowing
#   to rename certain parameters
# @return a brmsfit object with updated parameter names
do_renaming <- function(x, change) {
  .do_renaming <- function(x, change) {
    stopifnot(is.clist(change))
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
  chains <- length(x$fit@sim$samples)
  # temporary fix for issue #387 until fixed in rstan
  for (i in seq_len(chains)) {
    x$fit@sim$samples[[i]]$lp__.1 <- NULL
  }
  for (i in seq_along(change)) {
    x <- .do_renaming(x, change[[i]])
  }
  x
}

# order parameter draws after parameter class
# @param x brmsfit object
reorder_pars <- function(x) {
  all_classes <- unique(c(
    "b", "bs", "bsp", "bcs", "ar", "ma", "sderr", "lagsar", "errorsar", "car",
    "rhocar", "sdcar", "cosy", "correrr", "sd", "cor", "df", "sds", "sdgp",
    "lscale", valid_dpars(x), "lncor", "Intercept", "tmp", "rescor",
    "delta", "lasso", "simo", "r", "s", "zgp", "rcar", "sbhaz",
    "R2D2", "Ymi", "Yl", "meanme", "sdme", "corme", "Xme", "prior",
    "lprior", "lp"
  ))
  # reorder parameter classes
  class <- get_matches("^[^_]+", x$fit@sim$pars_oi)
  new_order <- order(
    factor(class, levels = all_classes),
    !grepl("_Intercept(_[[:digit:]]+)?$", x$fit@sim$pars_oi)
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
    # attributes of samples must be kept
    x$fit@sim$samples[[i]] <-
      subset_keep_attr(x$fit@sim$samples[[i]], new_order)
  }
  x
}

# wrapper function to compute and store quantities in the stanfit
# object which were not computed / stored by Stan itself
# @param x a brmsfit object
# @return a brmsfit object
compute_quantities <- function(x) {
  stopifnot(is.brmsfit(x))
  x <- compute_xi(x)
  x
}

# helper function to compute parameter xi, which is currently
# defined in the Stan model block and thus not being stored
# @param x a brmsfit object
# @return a brmsfit object
compute_xi <- function(x, ...) {
  UseMethod("compute_xi")
}

#' @export
compute_xi.brmsfit <- function(x, ...) {
  if (!any(grepl("^tmp_xi(_|$)", variables(x)))) {
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
compute_xi.mvbrmsprep <- function(x, fit, ...) {
  stopifnot(is.brmsfit(fit))
  for (resp in names(x$resps)) {
    fit <- compute_xi(x$resps[[resp]], fit = fit, ...)
  }
  fit
}

#' @export
compute_xi.brmsprep <- function(x, fit, ...) {
  stopifnot(is.brmsfit(fit))
  resp <- usc(x$resp)
  tmp_xi_name <- paste0("tmp_xi", resp)
  if (!tmp_xi_name %in% variables(fit)) {
    return(fit)
  }
  mu <- get_dpar(x, "mu")
  sigma <- get_dpar(x, "sigma")
  y <- matrix(x$data$Y, dim(mu)[1], dim(mu)[2], byrow = TRUE)
  bs <- -1 / matrixStats::rowRanges((y - mu) / sigma)
  bs <- matrixStats::rowRanges(bs)
  tmp_xi <- as.vector(as.matrix(fit, pars = tmp_xi_name))
  xi <- inv_logit(tmp_xi) * (bs[, 2] - bs[, 1]) + bs[, 1]
  # write xi into stanfit object
  xi_name <- paste0("xi", resp)
  samp_chain <- length(xi) / fit$fit@sim$chains
  for (i in seq_len(fit$fit@sim$chains)) {
    xi_part <- xi[((i - 1) * samp_chain + 1):(i * samp_chain)]
    # add warmup draws not used anyway
    xi_part <- c(rep(0, fit$fit@sim$warmup2[1]), xi_part)
    fit$fit@sim$samples[[i]][[xi_name]] <- xi_part
  }
  fit$fit@sim$pars_oi <- c(fit$fit@sim$pars_oi, xi_name)
  fit$fit@sim$dims_oi[[xi_name]] <- numeric(0)
  fit$fit@sim$fnames_oi <- c(fit$fit@sim$fnames_oi, xi_name)
  fit$fit@sim$n_flatnames <- fit$fit@sim$n_flatnames + 1
  fit
}
