#' @export
#' @rdname extract_draws
extract_draws.brmsfit <- function(x, newdata = NULL, re_formula = NULL, 
                                  allow_new_levels = FALSE,
                                  sample_new_levels = "uncertainty",
                                  incl_autocor = TRUE, oos = NULL, resp = NULL,
                                  nsamples = NULL, subset = NULL, nug = NULL, 
                                  smooths_only = FALSE, offset = TRUE, 
                                  new_objects = list(), point = NULL, ...) {
  snl_options <- c("uncertainty", "gaussian", "old_levels")
  sample_new_levels <- match.arg(sample_new_levels, snl_options)
  warn_brmsfit_multiple(x, newdata = newdata)
  x <- restructure(x)
  if (!incl_autocor) {
    x <- remove_autocor(x) 
  }
  resp <- validate_resp(resp, x)
  subset <- subset_samples(x, subset, nsamples)
  samples <- as.matrix(x, subset = subset)
  samples <- point_samples(samples, point = point)
  
  # prepare (new) data and stan data 
  newdata <- validate_newdata(
    newdata, object = x, re_formula = re_formula, 
    resp = resp, allow_new_levels = allow_new_levels, 
    new_objects = new_objects, ...
  )
  new <- !isTRUE(attr(newdata, "old"))
  sdata <- standata(
    x, newdata = newdata, re_formula = re_formula, 
    resp = resp, allow_new_levels = allow_new_levels, 
    new_objects = new_objects, internal = TRUE, ...
  )
  new_formula <- update_re_terms(x$formula, re_formula)
  bterms <- parse_bf(new_formula)
  ranef <- tidy_ranef(bterms, x$data)
  meef <- tidy_meef(bterms, x$data)
  old_sdata <- trunc_bounds <- NULL
  if (new) {
    # extract_draws_re() also requires the levels from newdata
    # original level names are already passed via old_ranef
    used_levels <- attr(tidy_ranef(bterms, newdata), "levels")
    attr(ranef, "levels") <- used_levels
    if (length(get_effect(bterms, "gp"))) {
      # GPs for new data require the original data as well
      old_sdata <- standata(x, internal = TRUE, ...)
    }
    if (length(get_effect(bterms, "sp"))) {
      # truncation bounds for imputing missing values in new data
      trunc_bounds <- trunc_bounds(bterms, data = newdata, incl_family = TRUE)
    }
  }
  draws_ranef <- extract_draws_ranef(
    ranef = ranef, samples = samples, sdata = sdata, 
    resp = resp, old_ranef = x$ranef, 
    sample_new_levels = sample_new_levels,
  )
  extract_draws(
    bterms, samples = samples, sdata = sdata, data = x$data, 
    draws_ranef = draws_ranef, meef = meef, resp = resp, 
    sample_new_levels = sample_new_levels, nug = nug, 
    smooths_only = smooths_only, offset = offset, new = new, 
    oos = oos, stanvars = names(x$stanvars), old_sdata = old_sdata,
    trunc_bounds = trunc_bounds
  )
}

extract_draws.mvbrmsterms <- function(x, samples, sdata, resp = NULL, ...) {
  resp <- validate_resp(resp, x$responses)
  if (length(resp) > 1) {
    if (has_subset(x)) {
      stop2("Argument 'resp' must be a single variable name ",
            "for models using addition argument 'subset'.")
    }
    draws <- list(nsamples = nrow(samples), nobs = sdata$N)
    draws$resps <- named_list(resp)
    draws$old_order <- attr(sdata, "old_order")
    for (r in resp) {
      draws$resps[[r]] <- extract_draws(
        x$terms[[r]], samples = samples, sdata = sdata, ...
      )
    }
    if (x$rescor) {
      draws$family <- draws$resps[[1]]$family
      draws$family$fun <- paste0(draws$family$family, "_mv")
      rescor <- get_cornames(resp, type = "rescor", brackets = FALSE)
      draws$mvpars$rescor <- get_samples(samples, rescor, fixed = TRUE)
      if (draws$family$family == "student") {
        # store in draws$dpars so that get_dpar can be called on nu
        draws$dpars$nu <- as.vector(get_samples(samples, "^nu$"))
      }
      draws$data$N <- draws$resps[[1]]$data$N
      draws$data$weights <- draws$resps[[1]]$data$weights
      Y <- lapply(draws$resps, function(x) x$data$Y)
      draws$data$Y <- do_call(cbind, Y)
    }
    draws <- structure(draws, class = "mvbrmsdraws")
  } else {
    draws <- extract_draws(
      x$terms[[resp]], samples = samples, sdata = sdata, ...
    )
  }
  draws
}

#' @export
extract_draws.brmsterms <- function(x, samples, sdata, data, ...) {
  data <- subset_data(data, x)
  nsamples <- nrow(samples)
  nobs <- sdata[[paste0("N", usc(x$resp))]]
  resp <- usc(combine_prefix(x))
  draws <- nlist(family = prepare_family(x), nsamples, nobs, resp = x$resp)
  draws$old_order <- attr(sdata, "old_order")
  valid_dpars <- valid_dpars(x)
  draws$dpars <- named_list(valid_dpars)
  for (dp in valid_dpars) {
    dp_regex <- paste0("^", dp, resp, "$")
    if (is.btl(x$dpars[[dp]]) || is.btnl(x$dpars[[dp]])) {
      draws$dpars[[dp]] <- extract_draws(
        x$dpars[[dp]], samples = samples, 
        sdata = sdata, data = data, ...
      )
    } else if (is.numeric(x$fdpars[[dp]]$value)) {
      draws$dpars[[dp]] <- x$fdpars[[dp]]$value
    } else if (any(grepl(dp_regex, colnames(samples)))) {
      draws$dpars[[dp]] <- as.vector(get_samples(samples, dp_regex))
    }
  }
  draws$nlpars <- named_list(names(x$nlpars))
  for (nlp in names(x$nlpars)) {
    draws$nlpars[[nlp]] <- extract_draws(
      x$nlpars[[nlp]], samples = samples, 
      sdata = sdata, data = data, ...
    )
  }
  if (is.mixfamily(x$family)) {
    families <- family_names(x$family)
    thetas <- paste0("theta", seq_along(families))
    if (any(ulapply(draws$dpars[thetas], is.list))) {
      # theta was predicted
      missing_id <- which(ulapply(draws$dpars[thetas], is.null))
      draws$dpars[[paste0("theta", missing_id)]] <- structure(
        as_draws_matrix(0, c(nsamples, nobs)), predicted = TRUE
      )
    } else {
      # theta was not predicted
      draws$dpars$theta <- do_call(cbind, draws$dpars[thetas])
      draws$dpars[thetas] <- NULL
      if (nrow(draws$dpars$theta) == 1L) {
        dim <- c(nrow(samples), ncol(draws$dpars$theta))
        draws$dpars$theta <- as_draws_matrix(draws$dpars$theta, dim = dim)
      }
    }
  } 
  if (is_ordinal(x$family)) {
    # it is better to handle ordinal thresholds outside the
    # main predictor term in particular for use in custom families
    if (is.mixfamily(x$family)) {
      mu_pars <- str_subset(names(x$dpars), "^mu[[:digit:]]+")
      for (mu in mu_pars) {
        draws$thres[[mu]] <- 
          extract_draws_thres(x$dpars[[mu]], samples, sdata, ...)
      }
    } else {
      draws$thres <- extract_draws_thres(x$dpars$mu, samples, sdata, ...)
    }
  }
  if (is_cox(x$family)) {
    # prepare baseline hazard functions for the Cox model
    if (is.mixfamily(x$family)) {
      mu_pars <- str_subset(names(x$dpars), "^mu[[:digit:]]+")
      for (mu in mu_pars) {
        draws$bhaz[[mu]] <- extract_draws_bhaz(
          x$dpars[[mu]], samples, sdata, ...
        )
      }
    } else {
      draws$bhaz <- extract_draws_bhaz(x$dpars$mu, samples, sdata, ...)
    }
  }
  if (has_cor_natural_residuals(x)) {
    # only include autocor samples on the top-level of draws 
    # when using the covariance formulation of autocorrelations
    draws$ac <- extract_draws_autocor(x, samples, sdata, ...)
  }
  draws$data <- extract_draws_data(x, sdata = sdata, data = data, ...)
  structure(draws, class = "brmsdraws")
}

#' @export
extract_draws.btnl <- function(x, samples, sdata, ...) {
  draws <- list(
    family = x$family, nlform = x$formula[[2]],
    nsamples = nrow(samples), 
    nobs = sdata[[paste0("N", usc(x$resp))]],
    used_nlpars = x$used_nlpars
  )
  class(draws) <- "bdrawsnl"
  p <- usc(combine_prefix(x))
  C <- sdata[[paste0("C", p)]]
  for (cov in colnames(C)) {
    dim <- c(nrow(samples), nrow(C))
    draws$C[[cov]] <- as_draws_matrix(C[, cov], dim = dim)
  }
  draws
}

#' @export
extract_draws.btl <- function(x, samples, sdata, smooths_only = FALSE, 
                              offset = TRUE, ...) {
  smooths_only <- as_one_logical(smooths_only)
  offset <- as_one_logical(offset)
  nsamples <- nrow(samples)
  nobs <- sdata[[paste0("N", usc(x$resp))]]
  draws <- nlist(family = x$family, nsamples, nobs)
  class(draws) <- "bdrawsl"
  if (smooths_only) {
    # make sure only smooth terms will be included in draws
    draws$sm <- extract_draws_sm(x, samples, sdata, ...)
    return(draws)
  }
  draws$fe <- extract_draws_fe(x, samples, sdata, ...)
  draws$sp <- extract_draws_sp(x, samples, sdata, ...)
  draws$cs <- extract_draws_cs(x, samples, sdata, ...)
  draws$sm <- extract_draws_sm(x, samples, sdata, ...)
  draws$gp <- extract_draws_gp(x, samples, sdata, ...)
  draws$re <- extract_draws_re(x, sdata, ...)
  if (offset) {
    draws$offset <- extract_draws_offset(x, sdata, ...)
  }
  if (!has_cor_natural_residuals(x)) {
    draws$ac <- extract_draws_autocor(x, samples, sdata, ...)
  }
  draws
}

# extract draws of ordinary population-level effects
extract_draws_fe <- function(bterms, samples, sdata, ...) {
  draws <- list()
  p <- usc(combine_prefix(bterms))
  X <- sdata[[paste0("X", p)]]
  fixef <- colnames(X)
  if (length(fixef)) {
    draws$X <- X
    b_pars <- paste0("b", p, "_", fixef)
    draws$b <- get_samples(samples, b_pars, fixed = TRUE)
  }
  draws
}

# extract draws of special effects terms
extract_draws_sp <- function(bterms, samples, sdata, data, meef, 
                             new = FALSE, trunc_bounds = NULL, ...) {
  draws <- list()
  spef <- tidy_spef(bterms, data)
  if (!nrow(spef)) return(draws)
  p <- usc(combine_prefix(bterms))
  resp <- usc(bterms$resp)
  # prepare calls evaluated in sp_predictor
  draws$calls <- vector("list", nrow(spef))
  for (i in seq_along(draws$calls)) {
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
    draws$calls[[i]] <- parse(text = paste0(call))
  }
  # extract general data and parameters for special effects
  bsp_pars <- paste0("bsp", p, "_", spef$coef)
  draws$bsp <- get_samples(samples, bsp_pars, fixed = TRUE)
  # prepare draws specific to monotonic effects
  simo_coef <- get_simo_labels(spef)
  Jmo <- sdata[[paste0("Jmo", p)]]
  draws$simo <- draws$Xmo <- named_list(simo_coef)
  for (i in seq_along(simo_coef)) {
    J <- seq_len(Jmo[i])
    simo_par <- paste0("simo", p, "_", simo_coef[i], "[", J, "]")
    draws$simo[[i]] <- get_samples(samples, simo_par, fixed = TRUE)
    draws$Xmo[[i]] <- sdata[[paste0("Xmo", p, "_", i)]]
  }
  # prepare draws specific to noise-free effects
  warn_me <- FALSE
  if (nrow(meef)) {
    save_mevars <- any(grepl("^Xme_", colnames(samples)))
    warn_me <- warn_me || !new && !save_mevars
    draws$Xme <- named_list(meef$coef)
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
          draws$Xme[[k]] <- get_samples(samples, Xme_pars[k])
        }
      } else {
        # sample new values of latent variables
        if (nzchar(g)) {
          # represent all indices between 1 and length(unique(Jme))
          Jme <- as.numeric(factor(Jme))
          me_dim <- c(nrow(draws$bsp), max(Jme))
        } else {
          me_dim <- c(nrow(draws$bsp), sdata$N)
        }
        for (k in K) {
          dXn <- as_draws_matrix(Xn[[k]], me_dim)
          dnoise <- as_draws_matrix(noise[[k]], me_dim)
          draws$Xme[[k]] <- array(rnorm(prod(me_dim), dXn, dnoise), me_dim)
          remove(dXn, dnoise)
        }
      }
      if (nzchar(g)) {
        for (k in K) {
          draws$Xme[[k]] <- draws$Xme[[k]][, Jme, drop = FALSE]
        }
      }
    }
  }
  # prepare draws specific to missing value variables
  dim <- c(nrow(draws$bsp), sdata[[paste0("N", resp)]])
  vars_mi <- unique(unlist(spef$vars_mi))
  if (length(vars_mi)) {
    # we know at this point that the model is multivariate
    Yl_names <- paste0("Yl_", vars_mi)
    draws$Yl <- named_list(Yl_names)
    for (i in seq_along(draws$Yl)) {
      vmi <- vars_mi[i]
      Y <- as_draws_matrix(sdata[[paste0("Y_", vmi)]], dim)
      sdy <- sdata[[paste0("noise_", vmi)]]
      if (is.null(sdy)) {
        # missings only
        draws$Yl[[i]] <- Y
        if (!new) {
          Ymi_pars <- paste0("Ymi_", vmi, "\\[")
          Ymi <- get_samples(samples, Ymi_pars)
          Jmi <- sdata[[paste0("Jmi_", vmi)]]
          draws$Yl[[i]][, Jmi] <- Ymi
        }
      } else {
        # measurement-error in the response
        save_mevars <- any(grepl("^Yl_", colnames(samples)))
        if (save_mevars && !new) {
          Yl_pars <- paste0("Yl_", vmi, "\\[")
          draws$Yl[[i]] <- get_samples(samples, Yl_pars)
        } else {
          warn_me <- warn_me || !new
          sdy <- as_draws_matrix(sdy, dim)
          draws$Yl[[i]] <- rcontinuous(
            n = prod(dim), dist = "norm", 
            mean = Y, sd = sdy,
            lb = trunc_bounds[[vmi]]$lb,
            ub = trunc_bounds[[vmi]]$ub
          )
          draws$Yl[[i]] <- array(draws$Yl[[i]], dim)
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
  for (i in seq_len(ncovars)) {
    draws$Csp[[i]] <- sdata[[paste0("Csp", p, "_", i)]]
    draws$Csp[[i]] <- as_draws_matrix(draws$Csp[[i]], dim = dim)
  }
  draws
}

# extract draws of category specific effects
extract_draws_cs <- function(bterms, samples, sdata, data, ...) {
  draws <- list()
  if (is_ordinal(bterms$family)) {
    resp <- usc(bterms$resp)
    draws$nthres <- sdata[[paste0("nthres", resp)]]
    csef <- colnames(get_model_matrix(bterms$cs, data))
    if (length(csef)) {
      p <- usc(combine_prefix(bterms))
      cs_pars <- paste0("^bcs", p, "_", csef, "\\[")
      draws$bcs <- get_samples(samples, cs_pars)
      draws$Xcs <- sdata[[paste0("Xcs", p)]]
    }
  }
  draws
}

# extract draws of smooth terms
extract_draws_sm <- function(bterms, samples, sdata, data, ...) {
  draws <- list()
  smef <- tidy_smef(bterms, data)
  if (!NROW(smef)) {
    return(draws)
  }
  p <- usc(combine_prefix(bterms))
  Xs_names <- attr(smef, "Xs_names")
  if (length(Xs_names)) {
    draws$fe$Xs <- sdata[[paste0("Xs", p)]]
    # allow for "b_" prefix for compatibility with version <= 2.5.0
    bspars <- paste0("^bs?", p, "_", escape_all(Xs_names), "$")
    draws$fe$bs <- get_samples(samples, bspars)
  }
  draws$re <- named_list(smef$label)
  for (i in seq_rows(smef)) {
    sm <- list()
    for (j in seq_len(smef$nbases[i])) {
      sm$Zs[[j]] <- sdata[[paste0("Zs", p, "_", i, "_", j)]]
      spars <- paste0("^s", p, "_", smef$label[i], "_", j, "\\[")
      sm$s[[j]] <- get_samples(samples, spars)
    }
    draws$re[[i]] <- sm
  }
  draws
}

# extract draws for Gaussian processes
# @param new is new data used?
# @param nug small numeric value to avoid numerical problems in GPs
extract_draws_gp <- function(bterms, samples, sdata, data,
                             new = FALSE, nug = NULL,
                             old_sdata = NULL, ...) {
  gpef <- tidy_gpef(bterms, data)
  if (!nrow(gpef)) {
    return(list())
  }
  if (new) {
    stopifnot(!is.null(old_sdata))
  }
  p <- usc(combine_prefix(bterms))
  if (is.null(nug)) {
    nug <- ifelse(new, 1e-8, 1e-11)
  }
  draws <- named_list(gpef$label)
  for (i in seq_along(draws)) {
    cons <- gpef$cons[[i]]
    if (length(cons)) {
      gp <- named_list(cons)
      for (j in seq_along(cons)) {
        gp[[j]] <- .extract_draws_gp(
          gpef, samples = samples, sdata = sdata,
          old_sdata = old_sdata, nug = nug, new = new,
          byj = j, p = p, i = i
        )
      }
      attr(gp, "byfac") <- TRUE
    } else {
      gp <- .extract_draws_gp(
        gpef, samples = samples, sdata = sdata,
        old_sdata = old_sdata, nug = nug, new = new, 
        p = p, i = i
      )
    }
    draws[[i]] <- gp
  }
  draws
}

# extract draws for Gaussian processes
# @param gpef output of tidy_gpef
# @param p prefix created by combine_prefix()
# @param i indiex of the Gaussian process
# @param byj index for the contrast of a categorical 'by' variable
# @return a list to be evaluated by .predictor_gp()
.extract_draws_gp <- function(gpef, samples, sdata, old_sdata,
                              nug, new, p, i, byj = NULL) {
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
    # approximate GPs don't differentiate between new and old data
    gp$x <- old_sdata[[Xgp_name]]
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

# extract draws for all group level effects
# needs to be separate from 'extract_draws_re' to take correlations
# across responses and distributional parameters into account (#779)
# @param ranef output of 'tidy_ranef' based on the new formula 
#   and old data but storing levels obtained from new data
# @param old_ranef same as 'ranef' but based on the original formula
# @return a named list with one element per group containing posterior draws 
#   of levels used in the data as well as additional meta-data
extract_draws_ranef <- function(ranef, samples, sdata, resp, old_ranef, 
                                sample_new_levels = "uncertainty", ...) {
  if (!nrow(ranef)) {
    return(NULL)
  }
  groups <- unique(ranef$group)
  out <- named_list(groups, list())
  for (g in groups) {
    # prepare general variables related to group g
    ranef_g <- subset2(ranef, group = g)
    old_ranef_g <- subset2(old_ranef, group = g)
    used_levels <- attr(ranef, "levels")[[g]]
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
    # only extract draws of effects specified in the new formula
    used_rpars <- match(ranef_g$cn, old_ranef_g$cn)
    used_rpars <- outer(seq_len(nlevels), (used_rpars - 1) * nlevels, "+")
    used_rpars <- as.vector(used_rpars)
    rsamples <- rsamples[, used_rpars, drop = FALSE]
    rsamples <- column_to_row_major_order(rsamples, nranef)
    # prepare data required for indexing parameters
    gtype <- ranef_g$gtype[1]
    resp_g <- intersect(ranef_g$resp, resp)[1]
    # any valid ID works here as J and W are independent of the ID
    id <- ranef_g$id[1]
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
    # store all information required in 'extract_draws_re'
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

# extract draws and data of group-level effects
# @param draws_ranef a named list with one element per group containing
#   posterior draws of levels as well as additional meta-data
extract_draws_re <- function(bterms, sdata, draws_ranef,
                             sample_new_levels = "uncertainty", ...) {
  draws <- list()
  if (!length(draws_ranef)) {
    return(draws) 
  }
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  ranef_px <- lapply(draws_ranef, "[[", "ranef")
  ranef_px <- do_call(rbind, ranef_px)
  ranef_px <- subset2(ranef_px, ls = px)
  if (!NROW(ranef_px)) {
    return(draws)
  }
  groups <- unique(ranef_px$group)
  # assigning S4 objects requires initialisation of list elements
  draws[c("Z", "Zsp", "Zcs")] <- list(named_list(groups))
  for (g in groups) {
    # extract variables specific to group 'g'
    ranef_g <- draws_ranef[[g]]$ranef
    ranef_g_px <- subset2(ranef_g, ls = px)
    rsamples <- draws_ranef[[g]]$rsamples
    nranef <- draws_ranef[[g]]$nranef
    levels <- draws_ranef[[g]]$levels
    max_level <- draws_ranef[[g]]$max_level
    gf <- draws_ranef[[g]]$gf
    weights <- draws_ranef[[g]]$weights
    # store samples and corresponding data in 'draws'
    # special group-level terms (mo, me, mi)
    ranef_g_px_sp <- subset2(ranef_g_px, type = "sp")
    if (nrow(ranef_g_px_sp)) {
      Z <- matrix(1, length(gf[[1]])) 
      draws[["Zsp"]][[g]] <- prepare_Z(Z, gf, max_level, weights)
      for (co in ranef_g_px_sp$coef) {
        # select from all varying effects of that group
        select <- find_rows(ranef_g, ls = px) & 
          ranef_g$coef == co & ranef_g$type == "sp"
        select <- which(select)
        select <- select + nranef * (seq_along(levels) - 1)
        draws[["rsp"]][[co]][[g]] <- rsamples[, select, drop = FALSE]
      }
    }
    # category specific group-level terms
    ranef_g_px_cs <- subset2(ranef_g_px, type = "cs")
    if (nrow(ranef_g_px_cs)) {
      # all categories share the same Z matrix
      ranef_g_px_cs_1 <- ranef_g_px_cs[grepl("\\[1\\]$", ranef_g_px_cs$coef), ]
      Znames <- paste0("Z_", ranef_g_px_cs_1$id, p, "_", ranef_g_px_cs_1$cn) 
      Z <- do_call(cbind, sdata[Znames])
      draws[["Zcs"]][[g]] <- prepare_Z(Z, gf, max_level, weights)
      for (i in seq_len(sdata$nthres)) {
        index <- paste0("\\[", i, "\\]$")
        # select from all varying effects of that group
        select <- find_rows(ranef_g, ls = px) &
          grepl(index, ranef_g$coef) & ranef_g$type == "cs"
        select <- which(select)
        select <- as.vector(outer(select, nranef * (seq_along(levels) - 1), "+"))
        draws[["rcs"]][[g]][[i]] <- rsamples[, select, drop = FALSE]
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
      draws[["Z"]][[g]] <- prepare_Z(Z, gf, max_level, weights)
      # select from all varying effects of that group
      select <- find_rows(ranef_g, ls = px) & ranef_g$type %in% c("", "mmc")
      select <- which(select)
      select <- as.vector(outer(select, nranef * (seq_along(levels) - 1), "+"))
      draws[["r"]][[g]] <- rsamples[, select, drop = FALSE]
    }
  }
  draws
}

extract_draws_offset <- function(bterms, sdata, ...) {
  p <- usc(combine_prefix(bterms))
  sdata[[paste0("offset", p)]]
}

# extract draws of ordinal thresholds
extract_draws_thres <- function(bterms, samples, sdata, ...) {
  draws <- list()
  if (!is_ordinal(bterms$family)) {
    return(draws)
  }
  resp <- usc(bterms$resp)
  draws$nthres <- sdata[[paste0("nthres", resp)]]
  draws$Jthres <- sdata[[paste0("Jthres", resp)]]
  p <- usc(combine_prefix(bterms))
  regex <- paste0("^b", p, "_Intercept\\[")
  draws$thres <- get_samples(samples, regex)
  draws
}

# extract draws of baseline functions for the cox model
extract_draws_bhaz <- function(bterms, samples, sdata, ...) {
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

# extract draws of autocorrelation parameters
extract_draws_autocor <- function(bterms, samples, sdata, oos = NULL, 
                                  new = FALSE, ...) {
  draws <- list()
  autocor <- draws$autocor <- bterms$autocor
  p <- usc(combine_prefix(bterms))
  draws$N_tg <- sdata[[paste0("N_tg", p)]]
  if (is.cor_arma(autocor)) {
    draws$Y <- sdata[[paste0("Y", p)]]
    if (!is.null(oos)) {
      if (any(oos > length(draws$Y))) {
        stop2("'oos' should not contain integers larger than N.")
      }
      # .predictor_arma has special behavior for NA responses 
      draws$Y[oos] <- NA
    }
    draws$J_lag <- sdata[[paste0("J_lag", p)]]
    if (get_ar(autocor)) {
      draws$ar <- get_samples(samples, paste0("^ar", p, "\\["))
    }
    if (get_ma(autocor)) {
      draws$ma <- get_samples(samples, paste0("^ma", p, "\\["))
    }
  }
  if (is.cor_cosy(autocor)) {
    draws$cosy <-  get_samples(samples, paste0("^cosy", p, "$"))
  }
  if (use_cov(autocor)) {
    # draws for the 'covariance' versions of ARMA and COSY structures
    draws$begin_tg <- sdata[[paste0("begin_tg", p)]]
    draws$end_tg <- sdata[[paste0("end_tg", p)]]
    if (has_latent_residuals(bterms)) {
      regex_err <- paste0("^err", p, "\\[")
      has_err <- any(grepl(regex_err, colnames(samples)))
      if (has_err && !new) {
        draws$err <- get_samples(samples, regex_err)
      } else {
        # need to sample correlated residuals
        draws$err <- matrix(nrow = nrow(samples), ncol = length(draws$Y))
        draws$sderr <- get_samples(samples, paste0("^sderr", p, "$"))
        for (i in seq_len(draws$N_tg)) {
          obs <- with(draws, begin_tg[i]:end_tg[i])
          zeros <- rep(0, length(obs))
          cov <- get_cov_matrix_autocor(list(ac = draws), obs, latent = TRUE)
          .err <- function(s) rmulti_normal(1, zeros, Sigma = cov[s, , ])
          draws$err[, obs] <- rblapply(seq_rows(samples), .err)
        }
      }
    }
  }
  if (is.cor_sar(autocor)) {
    draws$lagsar <- get_samples(samples, paste0("^lagsar", p, "$"))
    draws$errorsar <- get_samples(samples, paste0("^errorsar", p, "$"))
    draws$W <- sdata[[paste0("W", p)]]
  }
  if (is.cor_car(autocor)) {
    if (new) {
      group <- parse_time(autocor)$group
      if (!isTRUE(nzchar(group))) {
        stop2("Without a grouping factor, CAR models cannot handle newdata.")
      }
    }
    gcar <- sdata[[paste0("Jloc", p)]]
    Zcar <- matrix(rep(1, length(gcar)))
    draws$Zcar <- prepare_Z(Zcar, list(gcar))
    rcar <- get_samples(samples, paste0("^rcar", p, "\\["))
    rcar <- rcar[, unique(gcar), drop = FALSE]
    draws$rcar <- rcar
  }
  draws
}

# extract data mainly related to the response variable
# @param stanvars: *names* of variables stored in slot 'stanvars'
extract_draws_data <- function(bterms, sdata, data, stanvars = NULL, ...) {
  resp <- usc(combine_prefix(bterms))
  vars <- c(
    "Y", "trials", "ncat", "nthres", "se", "weights", 
    "dec", "cens", "rcens", "lb", "ub"
  )
  vars <- paste0(vars, resp)
  vars <- intersect(vars, names(sdata))
  # variables of variable length need to be handled via regular expression
  vl_vars <- c("vreal", "vint")
  vl_vars <- regex_or(vl_vars)
  vl_vars <- paste0("^", vl_vars, "[[:digit:]]+", escape_all(resp), "$")
  vl_vars <- str_subset(names(sdata), vl_vars)
  vars <- union(vars, vl_vars)
  draws <- sdata[vars]
  if (length(stanvars)) {
    stopifnot(is.character(stanvars))
    draws[stanvars] <- sdata[stanvars]
  }
  if (has_cat(bterms)) {
    draws$cats <- get_cats(bterms)
  }
  draws
}

# choose number of observations to be used in post-processing methods
choose_N <- function(draws) {
  stopifnot(is.brmsdraws(draws) || is.mvbrmsdraws(draws))
  if (!is.null(draws$ac$N_tg)) draws$ac$N_tg else draws$nobs
}

# create pseudo brmsdraws objects for components of mixture models
# @param comp the mixture component number
# @param sample_ids see predict_mixture
pseudo_draws_for_mixture <- function(draws, comp, sample_ids = NULL) {
  stopifnot(is.brmsdraws(draws), is.mixfamily(draws$family))
  if (!is.null(sample_ids)) {
    nsamples <- length(sample_ids)
  } else {
    nsamples <- draws$nsamples
  }
  out <- list(
    family = draws$family$mix[[comp]], nsamples = nsamples,
    nobs = draws$nobs, data = draws$data
  )
  out$family$fun <- out$family$family
  for (dp in valid_dpars(out$family)) {
    out$dpars[[dp]] <- draws$dpars[[paste0(dp, comp)]]
    if (length(sample_ids) && length(out$dpars[[dp]]) > 1L) {
      out$dpars[[dp]] <- p(out$dpars[[dp]], sample_ids, row = TRUE)
    }
  }
  if (is_ordinal(out$family)) {
    out$thres <- draws$thres[[paste0("mu", comp)]]
  }
  if (is_cox(out$family)) {
    out$bhaz <- draws$bhaz[[paste0("mu", comp)]]
  }
  # weighting should happen after computing the mixture
  out$data$weights <- NULL
  structure(out, class = "brmsdraws")
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
point_samples <- function(samples, point = NULL) {
  if (!is.null(point)) {
    point <- match.arg(point, c("mean", "median"))
    parnames <- colnames(samples)
    if (point == "mean") {
      samples <- matrixStats::colMeans2(samples)
    } else if (point == "median") {
      samples <- matrixStats::colMedians(samples)
    }
    samples <- t(samples)
    colnames(samples) <- parnames
  }
  samples
}

is.brmsdraws <- function(x) {
  inherits(x, "brmsdraws")
}

is.mvbrmsdraws <- function(x) {
  inherits(x, "mvbrmsdraws")
}

is.bdrawsl <- function(x) {
  inherits(x, "bdrawsl")
}

is.bdrawsnl <- function(x) {
  inherits(x, "bdrawsnl")
}

#' Extract Data and Posterior Draws
#' 
#' This method helps in preparing \pkg{brms} models for certin post-processing
#' tasks most notably various forms of predictions. Unless you are a package
#' developer, you will rarely need to call \code{extract_draws} directly.
#' 
#' @name extract_draws
#' @aliases extract_draws.brmsfit
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
#' @param new_objects A named \code{list} of objects containing new data, which
#'   cannot be passed via argument \code{newdata}. Required for objects passed
#'   via \code{\link{stanvars}} and for \code{\link[brms:cor_sar]{cor_sar}} and
#'   \code{\link[brms:cor_fixed]{cor_fixed}} correlation structures.
#' @param incl_autocor A flag indicating if correlation structures originally
#'   specified via \code{autocor} should be included in the predictions.
#'   Defaults to \code{TRUE}.
#' @param offset Logical; Indicates if offsets should be included in the
#'   predictions. Defaults to \code{TRUE}.
#' @param oos Optional indices of observations for which to compute
#'   out-of-sample rather than in-sample predictions. Only required in models
#'   that make use of response values to make predictions, that is currently
#'   only ARMA models.
#' @param smooths_only Logical; If \code{TRUE} only draws related to the
#'   computation of smooth terms will be extracted.
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
#' @param point Shall the returned object contain only point estimates of the
#'   parameters instead of their posterior samples? Defaults to \code{NULL} in
#'   which case no point estimate is computed. Alternatively, may be set to
#'   \code{"mean"} or \code{"median"}. This argument is primarily implemented to
#'   ensure compatibility with the \code{\link{loo_subsample}} method.
#' @param ... Further arguments passed to \code{\link{validate_newdata}}.
#'
#' @return An object of class \code{'brmsdraws'} or \code{'mvbrmsdraws'},
#'   depending on whether a univariate or multivariate model is passed.
#'   
#' @export
extract_draws <- function(x, ...) {
  UseMethod("extract_draws")
}

#' @export
extract_draws.default <- function(x, ...) {
  NULL
}
