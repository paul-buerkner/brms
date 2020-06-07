#' @export
#' @rdname prepare_predictions
prepare_predictions.brmsfit <- function(
  x, newdata = NULL, re_formula = NULL, 
  allow_new_levels = FALSE, sample_new_levels = "uncertainty",
  incl_autocor = TRUE, oos = NULL, resp = NULL, nsamples = NULL, 
  subset = NULL, nug = NULL, smooths_only = FALSE, offset = TRUE, 
  newdata2 = NULL, new_objects = NULL, point_estimate = NULL, ...
) {
  
  x <- restructure(x)
  # allows functions to fall back to old default behavior
  # which was used when originally fitting the model
  options(.brmsfit_version = x$version$brms)
  on.exit(options(.brmsfit_version = NULL))
  
  snl_options <- c("uncertainty", "gaussian", "old_levels")
  sample_new_levels <- match.arg(sample_new_levels, snl_options)
  warn_brmsfit_multiple(x, newdata = newdata)
  newdata2 <- use_alias(newdata2, new_objects)
  x <- exclude_terms(
    x, incl_autocor = incl_autocor, 
    offset = offset, smooths_only = smooths_only
  )
  resp <- validate_resp(resp, x)
  subset <- subset_samples(x, subset, nsamples)
  samples <- as.matrix(x, subset = subset)
  samples <- point_samples(samples, point_estimate)
  
  new_formula <- update_re_terms(x$formula, re_formula)
  bterms <- brmsterms(new_formula)
  ranef <- tidy_ranef(bterms, x$data)
  meef <- tidy_meef(bterms, x$data)
  new <- !is.null(newdata)
  sdata <- standata(
    x, newdata = newdata, re_formula = re_formula, 
    newdata2 = newdata2, resp = resp, 
    allow_new_levels = allow_new_levels, 
    internal = TRUE, ...
  )
  prep_ranef <- prepare_predictions_ranef(
    ranef = ranef, samples = samples, sdata = sdata, 
    resp = resp, old_ranef = x$ranef, 
    sample_new_levels = sample_new_levels,
  )
  prepare_predictions(
    bterms, samples = samples, sdata = sdata, data = x$data, 
    prep_ranef = prep_ranef, meef = meef, resp = resp, 
    sample_new_levels = sample_new_levels, nug = nug, 
    new = new, oos = oos, stanvars = x$stanvars
  )
}

prepare_predictions.mvbrmsterms <- function(x, samples, sdata, resp = NULL, ...) {
  resp <- validate_resp(resp, x$responses)
  if (length(resp) > 1) {
    if (has_subset(x)) {
      stop2("Argument 'resp' must be a single variable name ",
            "for models using addition argument 'subset'.")
    }
    out <- list(nsamples = nrow(samples), nobs = sdata$N)
    out$resps <- named_list(resp)
    out$old_order <- attr(sdata, "old_order")
    for (r in resp) {
      out$resps[[r]] <- prepare_predictions(
        x$terms[[r]], samples = samples, sdata = sdata, ...
      )
    }
    if (x$rescor) {
      out$family <- out$resps[[1]]$family
      out$family$fun <- paste0(out$family$family, "_mv")
      rescor <- get_cornames(resp, type = "rescor", brackets = FALSE)
      out$mvpars$rescor <- get_samples(samples, rescor, fixed = TRUE)
      if (out$family$family == "student") {
        # store in out$dpars so that get_dpar can be called on nu
        out$dpars$nu <- as.vector(get_samples(samples, "^nu$"))
      }
      out$data$N <- out$resps[[1]]$data$N
      out$data$weights <- out$resps[[1]]$data$weights
      Y <- lapply(out$resps, function(x) x$data$Y)
      out$data$Y <- do_call(cbind, Y)
    }
    out <- structure(out, class = "mvbrmsprep")
  } else {
    out <- prepare_predictions(
      x$terms[[resp]], samples = samples, sdata = sdata, ...
    )
  }
  out
}

#' @export
prepare_predictions.brmsterms <- function(x, samples, sdata, data, ...) {
  data <- subset_data(data, x)
  nsamples <- nrow(samples)
  nobs <- sdata[[paste0("N", usc(x$resp))]]
  resp <- usc(combine_prefix(x))
  out <- nlist(nsamples, nobs, resp = x$resp)
  out$family <- prepare_family(x)
  out$old_order <- attr(sdata, "old_order")
  valid_dpars <- valid_dpars(x)
  out$dpars <- named_list(valid_dpars)
  for (dp in valid_dpars) {
    dp_regex <- paste0("^", dp, resp, "$")
    if (is.btl(x$dpars[[dp]]) || is.btnl(x$dpars[[dp]])) {
      out$dpars[[dp]] <- prepare_predictions(
        x$dpars[[dp]], samples = samples, 
        sdata = sdata, data = data, ...
      )
    } else if (any(grepl(dp_regex, colnames(samples)))) {
      out$dpars[[dp]] <- as.vector(get_samples(samples, dp_regex))
    } else if (is.numeric(x$fdpars[[dp]]$value)) {
      # fixed dpars are stored as regular draws as of brms 2.12.9
      # so this manual extraction is only required for older models
      out$dpars[[dp]] <- x$fdpars[[dp]]$value
    }
  }
  out$nlpars <- named_list(names(x$nlpars))
  for (nlp in names(x$nlpars)) {
    out$nlpars[[nlp]] <- prepare_predictions(
      x$nlpars[[nlp]], samples = samples, 
      sdata = sdata, data = data, ...
    )
  }
  if (is.mixfamily(x$family)) {
    families <- family_names(x$family)
    thetas <- paste0("theta", seq_along(families))
    if (any(ulapply(out$dpars[thetas], is.list))) {
      # theta was predicted
      missing_id <- which(ulapply(out$dpars[thetas], is.null))
      out$dpars[[paste0("theta", missing_id)]] <- structure(
        as_draws_matrix(0, c(nsamples, nobs)), predicted = TRUE
      )
    } else {
      # theta was not predicted
      out$dpars$theta <- do_call(cbind, out$dpars[thetas])
      out$dpars[thetas] <- NULL
      if (nrow(out$dpars$theta) == 1L) {
        dim <- c(nrow(samples), ncol(out$dpars$theta))
        out$dpars$theta <- as_draws_matrix(out$dpars$theta, dim = dim)
      }
    }
  } 
  if (is_ordinal(x$family)) {
    # it is better to handle ordinal thresholds outside the
    # main predictor term in particular for use in custom families
    if (is.mixfamily(x$family)) {
      mu_pars <- str_subset(names(x$dpars), "^mu[[:digit:]]+")
      for (mu in mu_pars) {
        out$thres[[mu]] <- 
          prepare_predictions_thres(x$dpars[[mu]], samples, sdata, ...)
      }
    } else {
      out$thres <- prepare_predictions_thres(x$dpars$mu, samples, sdata, ...)
    }
  }
  if (is_cox(x$family)) {
    # prepare baseline hazard functions for the Cox model
    if (is.mixfamily(x$family)) {
      mu_pars <- str_subset(names(x$dpars), "^mu[[:digit:]]+")
      for (mu in mu_pars) {
        out$bhaz[[mu]] <- prepare_predictions_bhaz(
          x$dpars[[mu]], samples, sdata, ...
        )
      }
    } else {
      out$bhaz <- prepare_predictions_bhaz(x$dpars$mu, samples, sdata, ...)
    }
  }
  # response category names for categorical and ordinal models
  out$cats <- get_cats(x)
  # only include those autocor samples on the top-level 
  # of the output which imply covariance matrices on natural residuals
  out$ac <- prepare_predictions_ac(x$dpars$mu, samples, sdata, nat_cov = TRUE, ...)
  out$data <- prepare_predictions_data(x, sdata = sdata, data = data, ...)
  structure(out, class = "brmsprep")
}

#' @export
prepare_predictions.btnl <- function(x, samples, sdata, ...) {
  out <- list(
    family = x$family, nlform = x$formula[[2]],
    nsamples = nrow(samples), 
    nobs = sdata[[paste0("N", usc(x$resp))]],
    used_nlpars = x$used_nlpars
  )
  class(out) <- "bprepnl"
  p <- usc(combine_prefix(x))
  covars <- all.vars(x$covars)
  dim <- c(out$nsamples, out$nobs)
  for (i in seq_along(covars)) {
    cvalues <- sdata[[paste0("C", p, "_", i)]]
    out$C[[covars[i]]] <- as_draws_matrix(cvalues, dim = dim)
  }
  out
}

#' @export
prepare_predictions.btl <- function(x, samples, sdata, ...) {
  nsamples <- nrow(samples)
  nobs <- sdata[[paste0("N", usc(x$resp))]]
  out <- nlist(family = x$family, nsamples, nobs)
  class(out) <- "bprepl"
  out$fe <- prepare_predictions_fe(x, samples, sdata, ...)
  out$sp <- prepare_predictions_sp(x, samples, sdata, ...)
  out$cs <- prepare_predictions_cs(x, samples, sdata, ...)
  out$sm <- prepare_predictions_sm(x, samples, sdata, ...)
  out$gp <- prepare_predictions_gp(x, samples, sdata, ...)
  out$re <- prepare_predictions_re(x, sdata, ...)
  out$ac <- prepare_predictions_ac(x, samples, sdata, nat_cov = FALSE, ...)
  out$offset <- prepare_predictions_offset(x, sdata, ...)
  out
}

# prepare predictions of ordinary population-level effects
prepare_predictions_fe <- function(bterms, samples, sdata, ...) {
  out <- list()
  if (is.null(bterms[["fe"]])) {
    return(out)
  }
  p <- usc(combine_prefix(bterms))
  X <- sdata[[paste0("X", p)]]
  fixef <- colnames(X)
  if (length(fixef)) {
    out$X <- X
    b_pars <- paste0("b", p, "_", fixef)
    out$b <- get_samples(samples, b_pars, fixed = TRUE)
  }
  out
}

# prepare predictions of special effects terms
prepare_predictions_sp <- function(bterms, samples, sdata, data, 
                                   meef = empty_meef(), new = FALSE, ...) {
  out <- list()
  spef <- tidy_spef(bterms, data)
  if (!nrow(spef)) {
    return(out)
  }
  p <- usc(combine_prefix(bterms))
  resp <- usc(bterms$resp)
  # prepare calls evaluated in sp_predictor
  out$calls <- vector("list", nrow(spef))
  for (i in seq_along(out$calls)) {
    call <- spef$joint_call[[i]]
    if (!is.null(spef$calls_mo[[i]])) {
      new_mo <- paste0(".mo(simo_", spef$Imo[[i]], ", Xmo_", spef$Imo[[i]], ")")
      call <- rename(call, spef$calls_mo[[i]], new_mo)
    }
    if (!is.null(spef$calls_me[[i]])) {
      new_me <- paste0("Xme_", seq_along(meef$term))
      call <- rename(call, meef$term, new_me)
    }
    if (!is.null(spef$calls_mi[[i]])) {
      new_mi <- paste0("Yl_", spef$vars_mi[[i]])
      call <- rename(call, spef$calls_mi[[i]], new_mi)
    }
    if (spef$Ic[i] > 0) {
      str_add(call) <- paste0(" * Csp_", spef$Ic[i])
    }
    out$calls[[i]] <- parse(text = paste0(call))
  }
  # extract general data and parameters for special effects
  bsp_pars <- paste0("bsp", p, "_", spef$coef)
  out$bsp <- get_samples(samples, bsp_pars, fixed = TRUE)
  colnames(out$bsp) <- spef$coef
  # prepare predictions specific to monotonic effects
  simo_coef <- get_simo_labels(spef)
  Jmo <- sdata[[paste0("Jmo", p)]]
  out$simo <- out$Xmo <- named_list(simo_coef)
  for (i in seq_along(simo_coef)) {
    J <- seq_len(Jmo[i])
    simo_par <- paste0("simo", p, "_", simo_coef[i], "[", J, "]")
    out$simo[[i]] <- get_samples(samples, simo_par, fixed = TRUE)
    out$Xmo[[i]] <- sdata[[paste0("Xmo", p, "_", i)]]
  }
  # prepare predictions specific to noise-free effects
  warn_me <- FALSE
  if (nrow(meef)) {
    save_mevars <- any(grepl("^Xme_", colnames(samples)))
    warn_me <- warn_me || !new && !save_mevars
    out$Xme <- named_list(meef$coef)
    Xme_pars <- paste0("Xme_", escape_all(meef$coef), "\\[")
    Xn <- sdata[paste0("Xn_", seq_rows(meef))]
    noise <- sdata[paste0("noise_", seq_rows(meef))]
    groups <- unique(meef$grname)
    for (i in seq_along(groups)) {
      g <- groups[i]
      K <- which(meef$grname %in% g)
      if (nzchar(g)) {
        Jme <- sdata[[paste0("Jme_", i)]]
      }
      if (!new && save_mevars) {
        # extract original samples of latent variables
        for (k in K) {
          out$Xme[[k]] <- get_samples(samples, Xme_pars[k])
        }
      } else {
        # sample new values of latent variables
        if (nzchar(g)) {
          # TODO: reuse existing levels in predictions?
          # represent all indices between 1 and length(unique(Jme))
          Jme <- as.numeric(factor(Jme))
          me_dim <- c(nrow(out$bsp), max(Jme))
        } else {
          me_dim <- c(nrow(out$bsp), sdata$N)
        }
        for (k in K) {
          dXn <- as_draws_matrix(Xn[[k]], me_dim)
          dnoise <- as_draws_matrix(noise[[k]], me_dim)
          out$Xme[[k]] <- array(rnorm(prod(me_dim), dXn, dnoise), me_dim)
          remove(dXn, dnoise)
        }
      }
      if (nzchar(g)) {
        for (k in K) {
          out$Xme[[k]] <- out$Xme[[k]][, Jme, drop = FALSE]
        }
      }
    }
  }
  # prepare predictions specific to missing value variables
  dim <- c(nrow(out$bsp), sdata[[paste0("N", resp)]])
  vars_mi <- unique(unlist(spef$vars_mi))
  if (length(vars_mi)) {
    # we know at this point that the model is multivariate
    Yl_names <- paste0("Yl_", vars_mi)
    out$Yl <- named_list(Yl_names)
    for (i in seq_along(out$Yl)) {
      vmi <- vars_mi[i]
      Y <- as_draws_matrix(sdata[[paste0("Y_", vmi)]], dim)
      sdy <- sdata[[paste0("noise_", vmi)]]
      if (is.null(sdy)) {
        # missings only
        out$Yl[[i]] <- Y
        if (!new) {
          Ymi_pars <- paste0("Ymi_", vmi, "\\[")
          Ymi <- get_samples(samples, Ymi_pars)
          Jmi <- sdata[[paste0("Jmi_", vmi)]]
          out$Yl[[i]][, Jmi] <- Ymi
        }
      } else {
        # measurement-error in the response
        save_mevars <- any(grepl("^Yl_", colnames(samples)))
        if (save_mevars && !new) {
          Yl_pars <- paste0("Yl_", vmi, "\\[")
          out$Yl[[i]] <- get_samples(samples, Yl_pars)
        } else {
          warn_me <- warn_me || !new
          sdy <- as_draws_matrix(sdy, dim)
          out$Yl[[i]] <- rcontinuous(
            n = prod(dim), dist = "norm", 
            mean = Y, sd = sdy,
            lb = sdata[[paste0("lbmi_", vmi)]],
            ub = sdata[[paste0("ubmi_", vmi)]]
          )
          out$Yl[[i]] <- array(out$Yl[[i]], dim)
        }
      }
    }
  }
  if (warn_me) {
    warning2(
      "Noise-free variables were not saved. Please set ",
      "argument 'save_mevars' to TRUE when fitting your model. ",
      "Treating original data as if it was new data as a workaround."
    )
  }
  # prepare covariates
  ncovars <- max(spef$Ic)
  out$Csp <- vector("list", ncovars)
  for (i in seq_len(ncovars)) {
    out$Csp[[i]] <- sdata[[paste0("Csp", p, "_", i)]]
    out$Csp[[i]] <- as_draws_matrix(out$Csp[[i]], dim = dim)
  }
  out
}

# prepare predictions of category specific effects
prepare_predictions_cs <- function(bterms, samples, sdata, data, ...) {
  out <- list()
  if (!is_ordinal(bterms$family)) {
    return(out) 
  }
  resp <- usc(bterms$resp)
  out$nthres <- sdata[[paste0("nthres", resp)]]
  csef <- colnames(get_model_matrix(bterms$cs, data))
  if (length(csef)) {
    p <- usc(combine_prefix(bterms))
    cs_pars <- paste0("^bcs", p, "_", escape_all(csef), "\\[")
    out$bcs <- get_samples(samples, cs_pars)
    out$Xcs <- sdata[[paste0("Xcs", p)]]
  }
  out
}

# prepare predictions of smooth terms
prepare_predictions_sm <- function(bterms, samples, sdata, data, ...) {
  out <- list()
  smef <- tidy_smef(bterms, data)
  if (!NROW(smef)) {
    return(out)
  }
  p <- usc(combine_prefix(bterms))
  Xs_names <- attr(smef, "Xs_names")
  if (length(Xs_names)) {
    out$fe$Xs <- sdata[[paste0("Xs", p)]]
    # allow for "b_" prefix for compatibility with version <= 2.5.0
    bspars <- paste0("^bs?", p, "_", escape_all(Xs_names), "$")
    out$fe$bs <- get_samples(samples, bspars)
  }
  out$re <- named_list(smef$label)
  for (i in seq_rows(smef)) {
    sm <- list()
    for (j in seq_len(smef$nbases[i])) {
      sm$Zs[[j]] <- sdata[[paste0("Zs", p, "_", i, "_", j)]]
      spars <- paste0("^s", p, "_", smef$label[i], "_", j, "\\[")
      sm$s[[j]] <- get_samples(samples, spars)
    }
    out$re[[i]] <- sm
  }
  out
}

# prepare predictions for Gaussian processes
# @param new is new data used?
# @param nug small numeric value to avoid numerical problems in GPs
prepare_predictions_gp <- function(bterms, samples, sdata, data,
                                   new = FALSE, nug = NULL, ...) {
  gpef <- tidy_gpef(bterms, data)
  if (!nrow(gpef)) {
    return(list())
  }
  p <- usc(combine_prefix(bterms))
  if (is.null(nug)) {
    nug <- ifelse(new, 1e-8, 1e-11)
  }
  out <- named_list(gpef$label)
  for (i in seq_along(out)) {
    cons <- gpef$cons[[i]]
    if (length(cons)) {
      gp <- named_list(cons)
      for (j in seq_along(cons)) {
        gp[[j]] <- .prepare_predictions_gp(
          gpef, samples = samples, sdata = sdata,
          nug = nug, new = new, byj = j, p = p, i = i
        )
      }
      attr(gp, "byfac") <- TRUE
    } else {
      gp <- .prepare_predictions_gp(
        gpef, samples = samples, sdata = sdata,
        nug = nug, new = new, p = p, i = i
      )
    }
    out[[i]] <- gp
  }
  out
}

# prepare predictions for Gaussian processes
# @param gpef output of tidy_gpef
# @param p prefix created by combine_prefix()
# @param i indiex of the Gaussian process
# @param byj index for the contrast of a categorical 'by' variable
# @return a list to be evaluated by .predictor_gp()
.prepare_predictions_gp <- function(gpef, samples, sdata, nug, 
                                    new, p, i, byj = NULL) {
  sfx1 <- escape_all(gpef$sfx1[[i]])
  sfx2 <- escape_all(gpef$sfx2[[i]])
  if (is.null(byj)) {
    lvl <- ""
  } else {
    lvl <- gpef$bylevels[[i]][byj]
    sfx1 <- sfx1[byj]
    sfx2 <- sfx2[byj, ]
  }
  j <- usc(byj)
  pi <- paste0(p, "_", i)
  gp <- list()
  sdgp <- paste0("^sdgp", p, "_", sfx1, "$")
  gp$sdgp <- as.vector(get_samples(samples, sdgp))
  lscale <- paste0("^lscale", p, "_", sfx2, "$")
  gp$lscale <- get_samples(samples, lscale)
  zgp_regex <- paste0("^zgp", p, "_", sfx1, "\\[")
  gp$zgp <- get_samples(samples, zgp_regex)
  Xgp_name <- paste0("Xgp", pi, j)
  Igp_name <- paste0("Igp", pi, j)
  Jgp_name <- paste0("Jgp", pi, j)
  if (new && isNA(gpef$k[i])) {
    # in exact GPs old covariate values are required for predictions
    gp$x <- sdata[[paste0(Xgp_name, "_old")]]
    gp$nug <- 1e-11
    # computing GPs for new data requires the old GP terms
    gp$yL <- .predictor_gp(gp)
    gp$x_new <- sdata[[Xgp_name]]
    gp$Igp <- sdata[[Igp_name]]
  } else {
    gp$x <- sdata[[Xgp_name]]
    gp$Igp <- sdata[[Igp_name]]
    if (!isNA(gpef$k[i])) {
      gp$slambda <- sdata[[paste0("slambda", pi, j)]]
    }
  }
  gp$Jgp <- sdata[[Jgp_name]]
  # possible factor from 'by' variable
  gp$Cgp <- sdata[[paste0("Cgp", pi, j)]]
  gp$nug <- nug
  gp
}

# prepare predictions for all group level effects
# needs to be separate from 'prepare_predictions_re' to take correlations
# across responses and distributional parameters into account (#779)
# @param ranef output of 'tidy_ranef' based on the new formula and old data
# @param old_ranef same as 'ranef' but based on the original formula
# @return a named list with one element per group containing posterior draws 
#   of levels used in the data as well as additional meta-data
prepare_predictions_ranef <- function(ranef, samples, sdata, old_ranef, resp = NULL, 
                                      sample_new_levels = "uncertainty", ...) {
  if (!nrow(ranef)) {
    return(list())
  }
  # ensures subsetting 'ranef' by 'resp' works correctly
  resp <- resp %||% ""
  groups <- unique(ranef$group)
  out <- named_list(groups, list())
  for (g in groups) {
    # prepare general variables related to group g
    ranef_g <- subset2(ranef, group = g)
    old_ranef_g <- subset2(old_ranef, group = g)
    used_levels <- attr(sdata, "levels")[[g]]
    old_levels <- attr(old_ranef, "levels")[[g]]
    nlevels <- length(old_levels) 
    nranef <- nrow(ranef_g)
    # prepare samples of group-level effects
    rpars <- paste0("^r_", g, "(__.+)?\\[")
    rsamples <- get_samples(samples, rpars)
    if (!length(rsamples)) {
      stop2(
        "Group-level effects of group '", g, "' not found. ", 
        "Please set 'save_ranef' to TRUE when fitting your model."
      )
    }
    # only prepare predictions of effects specified in the new formula
    cols_match <- c("coef", "resp", "dpar", "nlpar")
    used_rpars <- which(find_rows(old_ranef_g, ls = ranef_g[cols_match]))
    used_rpars <- outer(seq_len(nlevels), (used_rpars - 1) * nlevels, "+")
    used_rpars <- as.vector(used_rpars)
    rsamples <- rsamples[, used_rpars, drop = FALSE]
    rsamples <- column_to_row_major_order(rsamples, nranef)
    # prepare data required for indexing parameters
    gtype <- ranef_g$gtype[1]
    resp_g <- intersect(ranef_g$resp, resp)[1]
    # any valid ID works here as J and W are independent of the ID
    id <- subset2(ranef_g, resp = resp)$id[1]
    idresp <- paste0(id, usc(resp_g))
    if (gtype == "mm") {
      ngf <- length(ranef_g$gcall[[1]]$groups)
      gf <- sdata[paste0("J_", idresp, "_", seq_len(ngf))]
      weights <- sdata[paste0("W_", idresp, "_", seq_len(ngf))]
    } else {
      gf <- sdata[paste0("J_", idresp)]
      weights <- list(rep(1, length(gf[[1]])))
    }
    # generate samples for new levels
    args_new_rsamples <- nlist(
      ranef = ranef_g, gf, used_levels, old_levels, 
      rsamples = rsamples, samples, sample_new_levels
    )
    new_rsamples <- do_call(get_new_rsamples, args_new_rsamples)
    max_level <- attr(new_rsamples, "max_level")
    gf <- attr(new_rsamples, "gf")
    rsamples <- cbind(rsamples, new_rsamples)
    # keep only those levels actually used in the current data
    levels <- unique(unlist(gf))
    rsamples <- subset_levels(rsamples, levels, nranef)
    # store all information required in 'prepare_predictions_re'
    out[[g]]$ranef <- ranef_g
    out[[g]]$rsamples <- rsamples
    out[[g]]$levels <- levels
    out[[g]]$nranef <- nranef
    out[[g]]$max_level <- max_level
    out[[g]]$gf <- gf
    out[[g]]$weights <- weights
  }
  out
}

# prepare predictions of group-level effects
# @param prep_ranef a named list with one element per group containing
#   posterior draws of levels as well as additional meta-data
prepare_predictions_re <- function(bterms, sdata, prep_ranef = list(),
                                   sample_new_levels = "uncertainty", ...) {
  out <- list()
  if (!length(prep_ranef)) {
    return(out) 
  }
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  ranef_px <- lapply(prep_ranef, "[[", "ranef")
  ranef_px <- do_call(rbind, ranef_px)
  ranef_px <- subset2(ranef_px, ls = px)
  if (!NROW(ranef_px)) {
    return(out)
  }
  groups <- unique(ranef_px$group)
  # assigning S4 objects requires initialisation of list elements
  out[c("Z", "Zsp", "Zcs")] <- list(named_list(groups))
  for (g in groups) {
    # extract variables specific to group 'g'
    ranef_g <- prep_ranef[[g]]$ranef
    ranef_g_px <- subset2(ranef_g, ls = px)
    rsamples <- prep_ranef[[g]]$rsamples
    nranef <- prep_ranef[[g]]$nranef
    levels <- prep_ranef[[g]]$levels
    max_level <- prep_ranef[[g]]$max_level
    gf <- prep_ranef[[g]]$gf
    weights <- prep_ranef[[g]]$weights
    # TODO: define 'select' according to parameter names not by position
    # store samples and corresponding data in the output
    # special group-level terms (mo, me, mi)
    ranef_g_px_sp <- subset2(ranef_g_px, type = "sp")
    if (nrow(ranef_g_px_sp)) {
      Z <- matrix(1, length(gf[[1]])) 
      out[["Zsp"]][[g]] <- prepare_Z(Z, gf, max_level, weights)
      for (co in ranef_g_px_sp$coef) {
        # select from all varying effects of that group
        select <- find_rows(ranef_g, ls = px) & 
          ranef_g$coef == co & ranef_g$type == "sp"
        select <- which(select)
        select <- select + nranef * (seq_along(levels) - 1)
        out[["rsp"]][[co]][[g]] <- rsamples[, select, drop = FALSE]
      }
    }
    # category specific group-level terms
    ranef_g_px_cs <- subset2(ranef_g_px, type = "cs")
    if (nrow(ranef_g_px_cs)) {
      # all categories share the same Z matrix
      ranef_g_px_cs_1 <- ranef_g_px_cs[grepl("\\[1\\]$", ranef_g_px_cs$coef), ]
      Znames <- paste0("Z_", ranef_g_px_cs_1$id, p, "_", ranef_g_px_cs_1$cn) 
      Z <- do_call(cbind, sdata[Znames])
      out[["Zcs"]][[g]] <- prepare_Z(Z, gf, max_level, weights)
      for (i in seq_len(sdata$nthres)) {
        index <- paste0("\\[", i, "\\]$")
        # select from all varying effects of that group
        select <- find_rows(ranef_g, ls = px) &
          grepl(index, ranef_g$coef) & ranef_g$type == "cs"
        select <- which(select)
        select <- as.vector(outer(select, nranef * (seq_along(levels) - 1), "+"))
        out[["rcs"]][[g]][[i]] <- rsamples[, select, drop = FALSE]
      }
    }
    # basic group-level terms
    ranef_g_px_basic <- subset2(ranef_g_px, type = c("", "mmc"))
    if (nrow(ranef_g_px_basic)) {
      Znames <- paste0("Z_", ranef_g_px_basic$id, p, "_", ranef_g_px_basic$cn)
      if (ranef_g_px_basic$gtype[1] == "mm") {
        ng <- length(ranef_g_px_basic$gcall[[1]]$groups)
        Z <- vector("list", ng)
        for (k in seq_len(ng)) {
          Z[[k]] <- do_call(cbind, sdata[paste0(Znames, "_", k)])
        }
      } else {
        Z <- do_call(cbind, sdata[Znames])
      }
      out[["Z"]][[g]] <- prepare_Z(Z, gf, max_level, weights)
      # select from all varying effects of that group
      select <- find_rows(ranef_g, ls = px) & ranef_g$type %in% c("", "mmc")
      select <- which(select)
      select <- as.vector(outer(select, nranef * (seq_along(levels) - 1), "+"))
      out[["r"]][[g]] <- rsamples[, select, drop = FALSE]
    }
  }
  out
}

# prepare predictions of autocorrelation parameters
# @param nat_cov extract terms for covariance matrices of natural residuals?
prepare_predictions_ac <- function(bterms, samples, sdata, oos = NULL, 
                                   nat_cov = FALSE, new = FALSE, ...) {
  out <- list()
  nat_cov <- as_one_logical(nat_cov)
  acef <- tidy_acef(bterms)
  acef <- subset2(acef, nat_cov = nat_cov)
  if (!NROW(acef)) {
    return(out)
  }
  out$acef <- acef
  p <- usc(combine_prefix(bterms))
  out$N_tg <- sdata[[paste0("N_tg", p)]]
  if (has_ac_class(acef, "arma")) {
    acef_arma <- subset2(acef, class = "arma")
    out$Y <- sdata[[paste0("Y", p)]]
    if (!is.null(oos)) {
      if (any(oos > length(out$Y))) {
        stop2("'oos' should not contain integers larger than N.")
      }
      # .predictor_arma has special behavior for NA responses 
      out$Y[oos] <- NA
    }
    out$J_lag <- sdata[[paste0("J_lag", p)]]
    if (acef_arma$p > 0) {
      out$ar <- get_samples(samples, paste0("^ar", p, "\\["))
    }
    if (acef_arma$q > 0) {
      out$ma <- get_samples(samples, paste0("^ma", p, "\\["))
    }
  }
  if (has_ac_class(acef, "cosy")) {
    out$cosy <-  get_samples(samples, paste0("^cosy", p, "$"))
  }
  if (use_ac_cov_time(acef)) {
    # prepare predictions for the covariance structures of time-series models
    out$begin_tg <- sdata[[paste0("begin_tg", p)]]
    out$end_tg <- sdata[[paste0("end_tg", p)]]
    if (has_cor_latent_residuals(bterms)) {
      regex_err <- paste0("^err", p, "\\[")
      has_err <- any(grepl(regex_err, colnames(samples)))
      if (has_err && !new) {
        out$err <- get_samples(samples, regex_err)
      } else {
        # need to sample correlated residuals
        out$err <- matrix(nrow = nrow(samples), ncol = length(out$Y))
        out$sderr <- get_samples(samples, paste0("^sderr", p, "$"))
        for (i in seq_len(out$N_tg)) {
          obs <- with(out, begin_tg[i]:end_tg[i])
          zeros <- rep(0, length(obs))
          cov <- get_cov_matrix_ac(list(ac = out), obs, latent = TRUE)
          .err <- function(s) rmulti_normal(1, zeros, Sigma = cov[s, , ])
          out$err[, obs] <- rblapply(seq_rows(samples), .err)
        }
      }
    }
  }
  if (has_ac_class(acef, "sar")) {
    out$lagsar <- get_samples(samples, paste0("^lagsar", p, "$"))
    out$errorsar <- get_samples(samples, paste0("^errorsar", p, "$"))
    out$Msar <- sdata[[paste0("Msar", p)]]
  }
  if (has_ac_class(acef, "car")) {
    acef_car <- subset2(acef, class = "car")
    if (new && acef_car$gr == "NA") {
      stop2("Without a grouping factor, CAR models cannot handle newdata.")
    }
    gcar <- sdata[[paste0("Jloc", p)]]
    Zcar <- matrix(rep(1, length(gcar)))
    out$Zcar <- prepare_Z(Zcar, list(gcar))
    rcar <- get_samples(samples, paste0("^rcar", p, "\\["))
    rcar <- rcar[, unique(gcar), drop = FALSE]
    out$rcar <- rcar
  }
  if (has_ac_class(acef, "fcor")) {
    out$Mfcor <- sdata[[paste0("Mfcor", p)]]
  }
  out
}

prepare_predictions_offset <- function(bterms, sdata, ...) {
  p <- usc(combine_prefix(bterms))
  sdata[[paste0("offsets", p)]]
}

# prepare predictions of ordinal thresholds
prepare_predictions_thres <- function(bterms, samples, sdata, ...) {
  out <- list()
  if (!is_ordinal(bterms$family)) {
    return(out)
  }
  resp <- usc(bterms$resp)
  out$nthres <- sdata[[paste0("nthres", resp)]]
  out$Jthres <- sdata[[paste0("Jthres", resp)]]
  p <- usc(combine_prefix(bterms))
  regex <- paste0("^b", p, "_Intercept\\[")
  out$thres <- get_samples(samples, regex)
  out
}

# prepare predictions of baseline functions for the cox model
prepare_predictions_bhaz <- function(bterms, samples, sdata, ...) {
  if (!is_cox(bterms$family)) {
    return(NULL)
  }
  out <- list()
  p <- usc(combine_prefix(bterms))
  sbhaz <- get_samples(samples, paste0("^sbhaz", p))
  Zbhaz <- sdata[[paste0("Zbhaz", p)]]
  out$bhaz <- tcrossprod(sbhaz, Zbhaz) 
  Zcbhaz <- sdata[[paste0("Zcbhaz", p)]]
  out$cbhaz <- tcrossprod(sbhaz, Zcbhaz)
  out
}

# extract data mainly related to the response variable
prepare_predictions_data <- function(bterms, sdata, data, stanvars = NULL, ...) {
  resp <- usc(combine_prefix(bterms))
  vars <- c(
    "Y", "trials", "ncat", "nthres", "se", "weights", 
    "denom", "dec", "cens", "rcens", "lb", "ub"
  )
  vars <- paste0(vars, resp)
  vars <- intersect(vars, names(sdata))
  # variables of variable length need to be handled via regular expression
  vl_vars <- c("vreal", "vint")
  vl_vars <- regex_or(vl_vars)
  vl_vars <- paste0("^", vl_vars, "[[:digit:]]+", escape_all(resp), "$")
  vl_vars <- str_subset(names(sdata), vl_vars)
  vars <- union(vars, vl_vars)
  out <- sdata[vars]
  if (length(stanvars)) {
    stopifnot(is.stanvars(stanvars))
    out[names(stanvars)] <- sdata[names(stanvars)]
  }
  out
}

# choose number of observations to be used in post-processing methods
choose_N <- function(prep) {
  stopifnot(is.brmsprep(prep) || is.mvbrmsprep(prep))
  if (!is.null(prep$ac$N_tg)) prep$ac$N_tg else prep$nobs
}

# create pseudo brmsprep objects for components of mixture models
# @param comp the mixture component number
# @param sample_ids see predict_mixture
pseudo_prep_for_mixture <- function(prep, comp, sample_ids = NULL) {
  stopifnot(is.brmsprep(prep), is.mixfamily(prep$family))
  if (!is.null(sample_ids)) {
    nsamples <- length(sample_ids)
  } else {
    nsamples <- prep$nsamples
  }
  out <- list(
    family = prep$family$mix[[comp]], nsamples = nsamples,
    nobs = prep$nobs, data = prep$data
  )
  out$family$fun <- out$family$family
  for (dp in valid_dpars(out$family)) {
    out$dpars[[dp]] <- prep$dpars[[paste0(dp, comp)]]
    if (length(sample_ids) && length(out$dpars[[dp]]) > 1L) {
      out$dpars[[dp]] <- p(out$dpars[[dp]], sample_ids, row = TRUE)
    }
  }
  if (is_ordinal(out$family)) {
    out$thres <- prep$thres[[paste0("mu", comp)]]
  }
  if (is_cox(out$family)) {
    out$bhaz <- prep$bhaz[[paste0("mu", comp)]]
  }
  # weighting should happen after computing the mixture
  out$data$weights <- NULL
  structure(out, class = "brmsprep")
}

# take relevant cols of a matrix of group-level terms
# if only a subset of levels is provided (for newdata)
# @param x a matrix typically samples of r or Z design matrices
#   samples need to be stored in row major order
# @param levels grouping factor levels to keep
# @param nranef number of group-level effects
subset_levels <- function(x, levels, nranef) {
  take_levels <- ulapply(levels, 
    function(l) ((l - 1) * nranef + 1):(l * nranef)
  )
  x[, take_levels, drop = FALSE]
}

# transform x from column to row major order
# rows represent levels and columns represent effects
# @param x a matrix of samples of group-level parameters
# @param nranef number of group-level effects
column_to_row_major_order <- function(x, nranef) {
  nlevels <- ncol(x) / nranef
  sort_levels <- ulapply(seq_len(nlevels),
    function(l) seq(l, ncol(x), by = nlevels)
  )
  x[, sort_levels, drop = FALSE]
}

# prepare group-level design matrices for use in 'predictor'
# @param Z (list of) matrices to be prepared
# @param gf (list of) vectors containing grouping factor values
# @param weights optional (list of) weights of the same length as gf
# @param max_level maximal level of 'gf'
# @return a sparse matrix representation of Z
prepare_Z <- function(Z, gf, max_level = NULL, weights = NULL) {
  if (!is.list(Z)) {
    Z <- list(Z)
  }
  if (!is.list(gf)) {
    gf <- list(gf)
  }
  if (is.null(weights)) {
    weights <- rep(1, length(gf[[1]]))
  }
  if (!is.list(weights)) {
    weights <- list(weights)
  }
  if (is.null(max_level)) {
    max_level <- max(unlist(gf))
  }
  levels <- unique(unlist(gf))
  nranef <- ncol(Z[[1]])
  Z <- mapply(
    expand_matrix, A = Z, x = gf, weights = weights,
    MoreArgs = nlist(max_level)
  )
  Z <- Reduce("+", Z)
  subset_levels(Z, levels, nranef)
}

# expand a matrix into a sparse matrix of higher dimension
# @param A matrix to be expanded
# @param x levels to expand the matrix
# @param max_level maximal number of levels that x can take on
# @param weights weights to apply to rows of A before expanding
# @param a sparse matrix of dimension nrow(A) x (ncol(A) * max_level)
expand_matrix <- function(A, x, max_level = max(x), weights = 1) {
  stopifnot(is.matrix(A))
  stopifnot(length(x) == nrow(A))
  stopifnot(all(is_wholenumber(x) & x > 0))
  stopifnot(length(weights) %in% c(1, nrow(A), prod(dim(A))))
  A <- A * as.vector(weights)
  K <- ncol(A)
  i <- rep(seq_along(x), each = K)
  make_j <- function(n, K, x) K * (x[n] - 1) + 1:K
  j <- ulapply(seq_along(x), make_j, K = K, x = x)
  Matrix::sparseMatrix(
    i = i, j = j, x = as.vector(t(A)),
    dims = c(nrow(A), ncol(A) * max_level)
  )
}

# generate samples for new group levels
# @param ranef 'ranef_frame' object of only a single grouping variable
# @param gf list of vectors of level indices in the current data
# @param rsamples matrix of group-level samples in row major order
# @param used_levels names of levels used in the current data
# @param old_levels names of levels used in the original data
# @param sample_new_levels specifies the way in which new samples are generated
# @param samples optional matrix of samples from all model parameters
# @return a matrix of samples for new group levels
get_new_rsamples <- function(ranef, gf, rsamples, used_levels, old_levels,
                             sample_new_levels, samples = NULL) {
  snl_options <- c("uncertainty", "gaussian", "old_levels")
  sample_new_levels <- match.arg(sample_new_levels, snl_options)
  g <- unique(ranef$group)
  stopifnot(length(g) == 1L)
  stopifnot(is.list(gf))
  used_by_per_level <- attr(used_levels, "by")
  old_by_per_level <- attr(old_levels, "by")
  new_levels <- setdiff(used_levels, old_levels)
  nranef <- nrow(ranef)
  nlevels <- length(old_levels)
  max_level <- nlevels
  
  out <- vector("list", length(gf))
  for (i in seq_along(gf)) {
    has_new_levels <- any(gf[[i]] > nlevels)
    if (has_new_levels) {
      if (sample_new_levels %in% c("old_levels", "gaussian")) {
        new_indices <- sort(setdiff(gf[[i]], seq_len(nlevels)))
        out[[i]] <- matrix(NA, nrow(rsamples), nranef * length(new_indices))
        if (sample_new_levels == "old_levels") {
          for (j in seq_along(new_indices)) {
            # choose an existing person to take the parameters from
            if (length(old_by_per_level)) {
              new_by <- used_by_per_level[used_levels == new_levels[j]]
              possible_levels <- old_levels[old_by_per_level == new_by]
              possible_levels <- which(old_levels %in% possible_levels)
              take_level <- sample(possible_levels, 1)
            } else {
              take_level <- sample(seq_len(nlevels), 1)
            }
            for (k in seq_len(nranef)) {
              take <- (take_level - 1) * nranef + k
              out[[i]][, (j - 1) * nranef + k] <- rsamples[, take]
            }
          }
        } else if (sample_new_levels == "gaussian") {
          if (any(!ranef$dist %in% "gaussian")) {
            stop2("Option sample_new_levels = 'gaussian' is not ",
                  "available for non-gaussian group-level effects.")
          }
          for (j in seq_along(new_indices)) {
            # extract hyperparameters used to compute the covariance matrix
            if (length(old_by_per_level)) {
              new_by <- used_by_per_level[used_levels == new_levels[j]]
              rnames <- as.vector(get_rnames(ranef, bylevels = new_by))
            } else {
              rnames <- get_rnames(ranef)
            }
            sd_pars <- paste0("sd_", g, "__", rnames)
            sd_samples <- get_samples(samples, sd_pars, fixed = TRUE)
            cor_type <- paste0("cor_", g)
            cor_pars <- get_cornames(rnames, cor_type, brackets = FALSE)
            cor_samples <- matrix(0, nrow(sd_samples), length(cor_pars))
            for (k in seq_along(cor_pars)) {
              if (cor_pars[k] %in% colnames(samples)) {
                cor_samples[, k] <- get_samples(
                  samples, cor_pars[k], fixed = TRUE
                )
              }
            }
            cov_matrix <- get_cov_matrix(sd_samples, cor_samples)
            # sample new levels from the normal distribution
            # implied by the covariance matrix
            indices <- ((j - 1) * nranef + 1):(j * nranef)
            out[[i]][, indices] <- t(apply(
              cov_matrix, 1, rmulti_normal, 
              n = 1, mu = rep(0, length(sd_pars))
            ))
          }
        }
        max_level <- max_level + length(new_indices)
      } else if (sample_new_levels == "uncertainty") {
        out[[i]] <- matrix(nrow = nrow(rsamples), ncol = nranef)
        # selected levels need to be the same for all varying effects
        # to correctly take their correlations into account
        sel_levels <- sample(seq_len(nlevels), NROW(rsamples), TRUE)
        for (k in seq_len(nranef)) {
          indices <- seq(k, nlevels * nranef, by = nranef)
          tmp <- rsamples[, indices, drop = FALSE]
          for (j in seq_rows(tmp)) {
            out[[i]][j, k] <- tmp[j, sel_levels[j]]
          }
        }
        max_level <- max_level + 1
        gf[[i]][gf[[i]] > nlevels] <- max_level
      }
    } else { 
      out[[i]] <- matrix(nrow = nrow(rsamples), ncol = 0)
    }
  }
  out <- do_call(cbind, out)
  structure(out, gf = gf, max_level = max_level)
}

# extract samples of selected parameters
get_samples <- function(x, pars, ...) {
  pars <- extract_pars(pars, all_pars = colnames(x), ...)
  x[, pars, drop = FALSE]
}

# compute point estimates of posterior samples
# currently used primarily for 'loo_subsample'
# @param samples matrix of posterior samples
# @param point optional name of the point estimate to be computed
# @return a matrix with one row and as many columns as parameters
point_samples <- function(samples, point_estimate = NULL) {
  if (is.null(point_estimate)) {
    return(samples)
  }
  point_estimate <- match.arg(point_estimate, c("mean", "median"))
  parnames <- colnames(samples)
  if (point_estimate == "mean") {
    samples <- matrixStats::colMeans2(samples)
  } else if (point_estimate == "median") {
    samples <- matrixStats::colMedians(samples)
  }
  samples <- t(samples)
  colnames(samples) <- parnames
  samples
}

is.brmsprep <- function(x) {
  inherits(x, "brmsprep")
}

is.mvbrmsprep <- function(x) {
  inherits(x, "mvbrmsprep")
}

is.bprepl <- function(x) {
  inherits(x, "bprepl")
}

is.bprepnl <- function(x) {
  inherits(x, "bprepnl")
}

#' Prepare Predictions
#' 
#' This method helps in preparing \pkg{brms} models for certin post-processing
#' tasks most notably various forms of predictions. Unless you are a package
#' developer, you will rarely need to call \code{prepare_predictions} directly.
#' 
#' @name prepare_predictions
#' @aliases prepare_predictions.brmsfit extract_draws
#'
#' @param x An \R object typically of class \code{'brmsfit'}.
#' @param newdata An optional data.frame for which to evaluate predictions. If
#'   \code{NULL} (default), the original data of the model is used.
#'   \code{NA} values within factors are interpreted as if all dummy
#'   variables of this factor are zero. This allows, for instance, to make
#'   predictions of the grand mean when using sum coding.
#' @param re_formula formula containing group-level effects to be considered in
#'   the prediction. If \code{NULL} (default), include all group-level effects;
#'   if \code{NA}, include no group-level effects.
#' @param allow_new_levels A flag indicating if new levels of group-level
#'   effects are allowed (defaults to \code{FALSE}). Only relevant if
#'   \code{newdata} is provided.
#' @param sample_new_levels Indicates how to sample new levels for grouping
#'   factors specified in \code{re_formula}. This argument is only relevant if
#'   \code{newdata} is provided and \code{allow_new_levels} is set to
#'   \code{TRUE}. If \code{"uncertainty"} (default), include group-level
#'   uncertainty in the predictions based on the variation of the existing
#'   levels. If \code{"gaussian"}, sample new levels from the (multivariate)
#'   normal distribution implied by the group-level standard deviations and
#'   correlations. This options may be useful for conducting Bayesian power
#'   analysis. If \code{"old_levels"}, directly sample new levels from the
#'   existing levels.
#' @param newdata2 A named \code{list} of objects containing new data, which
#'   cannot be passed via argument \code{newdata}. Required for some objects 
#'   used in autocorrelation structures, or \code{\link{stanvars}}.
#' @param new_objects Deprecated alias of \code{newdata2}.
#' @param incl_autocor A flag indicating if correlation structures originally
#'   specified via \code{autocor} should be included in the predictions.
#'   Defaults to \code{TRUE}.
#' @param offset Logical; Indicates if offsets should be included in the
#'   predictions. Defaults to \code{TRUE}.
#' @param oos Optional indices of observations for which to compute
#'   out-of-sample rather than in-sample predictions. Only required in models
#'   that make use of response values to make predictions, that is currently
#'   only ARMA models.
#' @param smooths_only Logical; If \code{TRUE} only predictions related to the
#'   computation of smooth terms will be prepared.
#' @param resp Optional names of response variables. If specified, predictions
#'   are performed only for the specified response variables.
#' @param subset A numeric vector specifying the posterior samples to be used.
#'   If \code{NULL} (the default), all samples are used.
#' @param nsamples Positive integer indicating how many posterior samples should
#'   be used. If \code{NULL} (the default) all samples are used. Ignored if
#'   \code{subset} is not \code{NULL}.
#' @param nug Small positive number for Gaussian process terms only. For
#'   numerical reasons, the covariance matrix of a Gaussian process might not be
#'   positive definite. Adding a very small number to the matrix's diagonal
#'   often solves this problem. If \code{NULL} (the default), \code{nug} is
#'   chosen internally.
#' @param point_estimate Shall the returned object contain only point estimates
#'   of the parameters instead of their posterior samples? Defaults to
#'   \code{NULL} in which case no point estimate is computed. Alternatively, may
#'   be set to \code{"mean"} or \code{"median"}. This argument is primarily
#'   implemented to ensure compatibility with the \code{\link{loo_subsample}}
#'   method.
#' @param ... Further arguments passed to \code{\link{validate_newdata}}.
#'
#' @return An object of class \code{'brmsprep'} or \code{'mvbrmsprep'},
#'   depending on whether a univariate or multivariate model is passed.
#'   
#' @export
prepare_predictions <- function(x, ...) {
  UseMethod("prepare_predictions")
}

#' @export
prepare_predictions.default <- function(x, ...) {
  NULL
}

# the name 'extract_draws' is deprecated as of brms 2.12.6
# remove it eventually in brms 3.0
#' @export
extract_draws <- function(x, ...) {
  warning2("Method 'extract_draws' is deprecated. ", 
           "Please use 'prepare_predictions' instead.")
  UseMethod("prepare_predictions")
}
