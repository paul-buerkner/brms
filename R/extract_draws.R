extract_draws <- function(x, ...) {
  # extract data and posterior draws
  UseMethod("extract_draws")
}

#' @export
extract_draws.brmsfit <- function(x, newdata = NULL, re_formula = NULL, 
                                  allow_new_levels = FALSE, 
                                  sample_new_levels = "uncertainty",
                                  new_objects = list(), incl_autocor = TRUE, 
                                  subset = NULL, nsamples = NULL, 
                                  nug = NULL, ...) {
  # extract all data and posterior draws required in (non)linear_predictor
  # Args:
  #   see doc of logLik.brmsfit
  #   ...: passed to amend_newdata
  # Returns:
  #   A named list to be interpreted by linear_predictor
  snl_options <- c("uncertainty", "gaussian", "old_levels")
  sample_new_levels <- match.arg(sample_new_levels, snl_options)
  x <- restructure(x)
  x <- remove_autocor(x, incl_autocor)
  bterms <- parse_bf(x$formula)
  subset <- subset_samples(x, subset, nsamples)
  nsamples <- nsamples(x, subset = subset)
  newd_args <- nlist(
    fit = x, newdata, re_formula, allow_new_levels, new_objects
  )
  draws <- nlist(
    f = prepare_family(x), nsamples = nsamples,
    data = do.call(amend_newdata, c(newd_args, list(...)))
  )
  args <- c(newd_args, nlist(
    sample_new_levels, subset, nsamples, nug
  ))
  keep <- !grepl("^(X|Z|J)", names(draws$data))
  draws$data <- subset_attr(draws$data, keep)
  
  # special treatment of multivariate models
  resp <- bterms$response
  if (length(resp) > 1L && !isTRUE(x$formula[["old_mv"]])) {
    draws$mu[["mv"]] <- named_list(resp)
    for (j in seq_along(resp)) {
      r <- resp[j]
      more_args <- list(x = bterms$dpars[["mu"]], mv = TRUE)
      more_args$x$resp <- r
      draws$mu[["mv"]][[r]] <- do.call(extract_draws, c(args, more_args))
      draws$mu[["mv"]][[r]][["f"]] <- bterms$dpars[["mu"]]$family
      if (isTRUE(ncol(draws$mu[["mv"]][[r]]$data$Y) > 1L)) {
        draws$mu[["mv"]][[r]]$data$Y <- draws$mu[["mv"]][[r]]$data$Y[, j]
      }
    }
    bterms$dpars[["mu"]] <- NULL
  }
  # extract draws of auxiliary parameters
  am_args <- nlist(x, subset)
  valid_dpars <- valid_dpars(x$family, bterms = bterms)
  for (ap in valid_dpars) {
    ap_regex <- paste0("^", ap, "($|_)")
    if (is.btl(bterms$dpars[[ap]]) || is.btnl(bterms$dpars[[ap]])) {
      p <- usc(combine_prefix(bterms$dpars[[ap]]))
      more_args <- nlist(
        x = bterms$dpars[[ap]], C = draws$data[[paste0("C", p)]]
      )
      draws[[ap]] <- do.call(extract_draws, c(args, more_args))
      draws[[ap]][["f"]] <- bterms$dpars[[ap]]$family
    } else if (is.numeric(bterms$fdpars[[ap]]$value)) {
      draws[[ap]] <- bterms$fdpars[[ap]]$value
    } else if (any(grepl(ap_regex, parnames(x)))) {
      draws[[ap]] <- do.call(as.matrix, c(am_args, pars = ap_regex))
    }
  }
  if (is.mixfamily(x$family)) {
    families <- family_names(x$family)
    thetas <- paste0("theta", seq_along(families))
    if (any(ulapply(draws[thetas], is.list))) {
      # theta was predicted
      missing_id <- which(ulapply(draws[thetas], is.null))
      draws[[paste0("theta", missing_id)]] <- structure(
        as_draws_matrix(0, c(draws$nsamples, draws$data$N)),
        predicted = TRUE
      )
    } else {
      # theta was not predicted
      draws$theta <- do.call(cbind, draws[thetas])
      draws[thetas] <- NULL
      if (nrow(draws$theta) == 1L) {
        draws$theta <- as_draws_matrix(
          draws$theta, dim = c(nsamples, ncol(draws$theta))
        )
      }
    }
  }
  if (is_linear(x$family) && length(bterms$response) > 1L) {
    # parameters for multivariate normal models
    draws$sigma <- do.call(as.matrix, c(am_args, pars = "^sigma($|_)"))
    draws$rescor <- do.call(as.matrix, c(am_args, pars = "^rescor_"))
    draws$Sigma <- get_cov_matrix(sd = draws$sigma, cor = draws$rescor)$cov
  }
  if (use_cov(x$autocor) || is.cor_sar(x$autocor)) {
    # only include autocor samples on the top-level of draws 
    # when using the covariance formulation of ARMA / SAR structures
    draws <- c(draws, 
      extract_draws_autocor(
        x, subset, newdata = newdata, sdata = draws$data
      )
    )
  }
  draws
}

#' @export
extract_draws.btnl <- function(x, C, ...) {
  # Args:
  #   C: matrix containing covariates
  #   ...: passed to extract_draws.btl
  nlpars <- names(x$nlpars)
  draws <- named_list(nlpars)
  for (nlp in nlpars) {
    rhs_formula <- x$nlpars[[nlp]]$formula
    draws$nlpars[[nlp]] <- extract_draws(x$nlpars[[nlp]], ...)
  }
  nsamples <- draws$nlpars[[nlpars[1]]]$nsamples
  stopifnot(is.matrix(C))
  for (cov in colnames(C)) {
    draws[["C"]][[cov]] <- as_draws_matrix(
      C[, cov], dim = c(nsamples, nrow(C))
    )
  }
  draws$nlform <- x$formula[[2]]
  draws
}

#' @export
extract_draws.btl <- function(x, fit, newdata = NULL, re_formula = NULL, 
                              allow_new_levels = FALSE, 
                              sample_new_levels = "uncertainty",
                              new_objects = list(), incl_autocor = TRUE, 
                              subset = NULL, smooths_only = FALSE, 
                              mv = FALSE, nug = NULL, ...) {
  # extract draws of all kinds of effects
  # Args:
  #   fit: a brmsfit object
  #   mv: is the model multivariate?
  #   ...: further elements to store in draws
  #   other arguments: see extract_draws.brmsfit
  # Returns:
  #   A named list to be interpreted by linear_predictor
  stopifnot(is.brmsfit(fit))
  dots <- list(...)
  fit <- remove_autocor(fit, incl_autocor)
  px <- check_prefix(x)
  p <- combine_prefix(px)
  nsamples <- nsamples(fit, subset = subset)
  all_group_vars <- get_all_group_vars(fit)
  fit$formula$formula <- update(fit$formula$formula, rhs(x$formula))
  # ensure that auxiliary parameters are not included (fixes #154)
  fit$formula$pforms <- fit$formula$pfix <- NULL
  fit$ranef <- tidy_ranef(parse_bf(fit$formula), fit$data)
  if (nzchar(p)) {
    # make sure not to evaluate family specific stuff
    fit$formula[["response"]] <- NA
    fit$family <- fit$formula$family <- .dpar_family()
  }
  new_formula <- update_re_terms(fit$formula, re_formula)
  bterms <- parse_bf(new_formula)
  new_ranef <- tidy_ranef(bterms, fit$data)
  newd_args <- nlist(
    fit, newdata, re_formula, allow_new_levels, 
    new_objects, all_group_vars, check_response = FALSE
  )
  draws <- list(
    # supresses messages of add_new_objects
    data = suppressMessages(do.call(amend_newdata, newd_args)), 
    old_cat = is_old_categorical(fit)
  )
  draws[names(dots)] <- dots
  args <- nlist(x = fit, subset)
  new <- !is.null(newdata)
  # deviate from the usual way of passing bterms and data
  # as bterms and px do not coincide in extract_draws.btl
  fixef <- colnames(draws$data[["X"]])
  smooths <- get_sm_labels(bterms$dpars$mu, fit$data, covars = TRUE)
  draws <- c(draws,
    extract_draws_fe(fixef, args, px = px, old_cat = draws$old_cat),
    extract_draws_sm(smooths, args, sdata = draws$data, px = px)
  )
  if (smooths_only) {
    # make sure only smooth terms will be included in draws
    return(structure(draws, predicted = TRUE))
  }
  monef <- get_mo_labels(bterms$dpars$mu, fit$data)
  csef <- colnames(draws$data[["Xcs"]])
  meef <- get_me_labels(bterms$dpars$mu, fit$data)
  gpef <- get_gp_labels(bterms$dpars$mu, fit$data, covars = TRUE)
  draws <- c(draws,
    extract_draws_mo(monef, args, sdata = draws$data, px = px),
    extract_draws_cs(csef, args, px = px, old_cat = draws$old_cat),
    extract_draws_me(meef, args, sdata = draws$data, px = px, new = new),
    extract_draws_gp(
      gpef, args, sdata = draws$data, px = px, 
      new = new, nug = nug, newd_args = newd_args
    ),
    extract_draws_re(
      new_ranef, args, sdata = draws$data, px = px,
      sample_new_levels = sample_new_levels
    )
  )
  if (!nzchar(p) || mv) {
    # only include autocorrelation parameters in draws of mu
    if (!(use_cov(fit$autocor) || is.cor_sar(fit$autocor))) {
      draws <- c(draws, 
        extract_draws_autocor(
          fit, subset, newdata = newdata, sdata = draws$data
        )
      )
    }
  }
  structure(draws, predicted = TRUE)
}

extract_draws_fe <- function(fixef, args, px = list(), old_cat = 0L) {
  # extract draws of ordinary population-level effects
  # Args:
  #   fixef: names of the population-level effects
  #   args: list of arguments passed to as.matrix.brmsfit
  #   nlpar: name of a non-linear parameter
  #   old_cat: see old_categorical
  # Returns: 
  #   A named list to be interpreted by linear_predictor
  stopifnot("x" %in% names(args))
  p <- usc(combine_prefix(px), "suffix")
  draws <- list()
  if (length(fixef) && old_cat != 1L) {
    b_pars <- paste0("b_", p, fixef)
    draws[["b"]] <- do.call(as.matrix, 
      c(args, list(pars = b_pars, exact = TRUE))
    )
  }
  draws
}

extract_draws_mo <- function(monef, args, sdata, px = list()) {
  # extract draws of monotonic effects
  # Args:
  #   monef: names of the monotonic effects
  #   args: list of arguments passed to as.matrix.brmsfit
  #   sdata: list returned by make_standata
  # Returns: 
  #   A named list to be interpreted by linear_predictor
  stopifnot("x" %in% names(args))
  p <- usc(combine_prefix(px), "suffix")
  draws <- list()
  if (!length(monef)) {
    return(draws) 
  }
  att <- attributes(monef)
  # prepare calls evaluated in 'mo_predictor()'
  calls <- rep(NA, length(monef))
  for (i in seq_along(calls)) {
    tmp <- paste0(".mo(simo_", att$Imo[[i]], ", Xmo_", att$Imo[[i]], ")")
    calls[i] <- paste0(tmp, collapse = " * ")
    if (att$not_one[i]) {
      str_add(calls[i]) <- paste0(" * Cmo_", att$Icmo[i])
    }
  }
  calls <- parse(text = paste0(calls))
  # extract data and parameters for monotonic effects
  monef <- rename(monef)
  draws$bmo <- named_list(monef)
  for (i in seq_along(monef)) {
    bmo_par <- paste0("bmo_", p, monef[i])
    draws[["bmo"]][[i]] <- do.call(as.matrix,
      c(args, list(pars = bmo_par, exact = TRUE))
    )
  }
  attr(draws[["bmo"]], "calls") <- calls
  simo_coef <- get_simo_labels(monef)
  draws$simo <- named_list(simo_coef)
  for (i in seq_along(simo_coef)) {
    simo_par <- paste0(
      "simo_", p, simo_coef[i], "[", seq_len(sdata$Jmo[i]), "]"
    )
    draws[["simo"]][[i]] <- do.call(as.matrix, 
      c(args, list(pars = simo_par, exact = TRUE))
    )
  }
  ncovars <- sum(att$not_one)
  dim <- c(nrow(draws[["bmo"]][[1]]), sdata[["N"]])
  for (i in seq_len(ncovars)) {
    cmo <- paste0("Cmo_", i)
    draws[["Cmo"]][[i]] <- as_draws_matrix(sdata[[cmo]], dim)
  }
  draws
}

extract_draws_cs <- function(csef, args, px = list(), old_cat = 0L) {
  # category specific effects 
  # extract draws of category specific effects
  # Args:
  #   csef: names of the category specific effects
  #   args: list of arguments passed to as.matrix.brmsfit
  #   nlpar: name of a non-linear parameter
  #   old_cat: see old_categorical
  # Returns: 
  #   A named list to be interpreted by linear_predictor
  stopifnot("x" %in% names(args))
  p <- usc(combine_prefix(px), "suffix")
  draws <- list()
  if (is_ordinal(args$x$family)) {
    draws[["Intercept"]] <- do.call(as.matrix, 
      c(args, list(pars = "^b_Intercept\\["))
    )
    if (length(csef)) {
      # as of brms > 1.0.1 the original prefix 'bcs' is used
      bcs <- ifelse(any(grepl("^bcs_", parnames(args$x))), "^bcs_", "^b_")
      cs_pars <- paste0(bcs, p, csef, "\\[")
      draws[["cs"]] <- do.call(as.matrix, c(args, list(pars = cs_pars)))
    }
  } else if (old_cat == 1L) {
    # old categorical models deprecated as of brms > 0.8.0
    draws[["cs"]] <- do.call(as.matrix, c(args, list(pars = "^b_")))
  }
  draws
}

extract_draws_me <- function(meef, args, sdata, px = list(),
                             new = FALSE) {
  # extract draws of noise-free effects
  # Args:
  #   meef: names of the noise free effects
  #   args: list of arguments passed to as.matrix.brmsfit
  #   sdata: list returned by make_standata
  #   nlpar: name of a non-linear parameter
  #   new: logical; are new data specified?
  # Returns: 
  #   A named list to be interpreted by linear_predictor
  stopifnot("x" %in% names(args))
  p <- usc(combine_prefix(px), "suffix")
  draws <- list()
  if (length(meef)) {
    if (new) {
      stop2("Predictions with noise-free variables are not yet ",
            "possible when passing new data.")
    }
    if (!any(grepl("^Xme_", parnames(args$x)))) {
      stop2("Noise-free variables were not saved. Please set ",
            "argument 'save_mevars' to TRUE when fitting your model.")
    }
    uni_me <- attr(meef, "uni_me")
    not_one <- attr(meef, "not_one")
    # prepare calls to evaluate noise-free data
    me_sp <- strsplit(rm_wsp(meef), ":")
    meef_terms <- rep(NA, length(me_sp))
    for (i in seq_along(me_sp)) {
      # remove non-me parts from the terms
      take <- grepl_expr("^me\\([^:]*\\)$", me_sp[[i]])
      me_sp[[i]] <- me_sp[[i]][take]
      meef_terms[i] <- paste0(me_sp[[i]], collapse = " * ")
    }
    new_me <- paste0("Xme_", seq_along(uni_me))
    meef_terms <- rename(meef_terms, uni_me, new_me)
    ci <- ulapply(seq_along(not_one), function(i) sum(not_one[1:i]))
    covars <- ifelse(not_one, paste0(" * Cme_", ci), "")
    meef_terms <- paste0(meef_terms, covars)
    # extract coefficient samples
    meef <- rename(meef)
    draws[["bme"]] <- named_list(meef)
    attr(draws[["bme"]], "calls") <- parse(text = meef_terms)
    for (j in seq_along(draws[["bme"]])) {
      bme_par <- paste0("bme_", p, meef[j])
      draws[["bme"]][[j]] <- do.call(as.matrix, 
        c(args, list(pars = bme_par, exact = TRUE))
      )
    }
    # extract noise-free variable samples
    uni_me <- rename(uni_me)
    draws[["Xme"]] <- named_list(uni_me)
    N <- seq_len(nobs(args$x))
    for (j in seq_along(draws[["Xme"]])) {
      Xme_pars <- paste0("Xme_", uni_me[j], "[", N, "]")
      draws[["Xme"]][[j]] <- do.call(as.matrix,
        c(args, list(pars = Xme_pars, exact = TRUE))
      )
    }
    # prepare covariates
    ncovars <- sum(not_one)
    for (j in seq_len(ncovars)) {
      cme <- paste0("Cme_", j)
      draws[["Cme"]][[j]] <- as_draws_matrix(
        sdata[[cme]], dim = dim(draws[["Xme"]][[1]])
      )
    }
  } 
  draws
}

extract_draws_sm <- function(smooths, args, sdata, px = list()) {
  # extract draws of smooth terms
  # Args:
  #   smooths: names of the smooth terms
  #   args: list of arguments passed to as.matrix.brmsfit
  #   sdata: list returned by make_standata
  #   nlpar: name of a non-linear parameter
  # Returns: 
  #   A named list to be interpreted by linear_predictor
  stopifnot("x" %in% names(args))
  p <- usc(combine_prefix(px), "suffix")
  draws <- list()
  if (length(smooths)) {
    draws[["Zs"]] <- draws[["s"]] <- named_list(smooths)
    for (i in seq_along(smooths)) {
      nb <- seq_len(attr(smooths, "nbases")[[i]])
      for (j in nb) {
        draws[["Zs"]][[smooths[i]]][[j]] <- sdata[[paste0("Zs_", i, "_", j)]]
        s_pars <- paste0("^s_", p, smooths[i], "_", j, "\\[")
        draws[["s"]][[smooths[i]]][[j]] <- 
          do.call(as.matrix, c(args, list(pars = s_pars)))
      }
    }
  }
  draws
}

extract_draws_gp <- function(gpef, args, sdata, px = list(), 
                             new = FALSE, nug = NULL,
                             newd_args = list()) {
  # extract draws for gaussian processes
  # Args:
  #   gpef: names of the gaussian process terms
  stopifnot("x" %in% names(args))
  sdata_old <- NULL
  if (length(gpef) && new) {
    newd_args["newdata"] <- list(NULL)
    sdata_old <- do.call(amend_newdata, newd_args)
  }
  p <- usc(combine_prefix(px), "suffix")
  draws <- list()
  if (length(gpef)) {
    draws[["gp"]] <- named_list(gpef)
    if (is.null(nug)) {
      nug <- ifelse(new, 1e-8, 1e-11)
    }
    for (i in seq_along(gpef)) {
      gp <- list()
      by_levels <- attr(gpef, "by_levels")[[i]]
      gp_names <- paste0(gpef[i], usc(by_levels))
      sdgp <- paste0("^sdgp_", p, gp_names)
      gp[["sdgp"]] <- do.call(as.matrix, c(args, list(pars = sdgp)))
      lscale <- paste0("^lscale_", p, gp_names)
      gp[["lscale"]] <- do.call(as.matrix, c(args, list(pars = lscale)))
      zgp <- paste0("^zgp_", p, gpef[i], "\\[")
      gp[["zgp"]] <- do.call(as.matrix, c(args, list(pars = zgp)))
      gp[["bynum"]] <- sdata[[paste0("Cgp_", i)]]
      Jgp_pos <- grepl(paste0("^Jgp_", i), names(sdata))
      if (new) {
        gp[["x"]] <- sdata_old[[paste0("Xgp_", i)]]
        gp[["x_new"]] <- sdata[[paste0("Xgp_", i)]]
        gp[["Jgp"]] <- sdata_old[Jgp_pos]
        gp[["Jgp_new"]] <- sdata[Jgp_pos]
        # computing GPs for new data requires the old GP terms
        gp[["yL"]] <- gp_predictor(
          x = gp[["x"]], sdgp = gp[["sdgp"]],
          lscale = gp[["lscale"]], zgp = gp[["zgp"]], 
          Jgp = gp[["Jgp"]]
        )
      } else {
        gp[["x"]] <- sdata[[paste0("Xgp_", i)]]
        gp[["Jgp"]] <- sdata[Jgp_pos]
      }
      gp[["nug"]] <- nug
      draws[["gp"]][[i]] <- gp
    }
  }
  draws
}

extract_draws_re <- function(ranef, args, sdata, px = list(),
                             sample_new_levels = "uncertainty") {
  # extract draws of group-level effects
  # Args:
  #   ranef: data.frame returned by tidy_ranef
  #   args: list of arguments passed to as.matrix.brmsfit
  #   sdata: list returned by make_standata
  #   nlpar: name of a non-linear parameter
  #   sample_new_levels: see help("predict.brmsfit")
  # Returns: 
  #   A named list to be interpreted by linear_predictor
  stopifnot("x" %in% names(args))
  p <- combine_prefix(px)
  draws <- list()
  groups <- unique(ranef$group)
  # assigning S4 objects requires initialisation of list elements
  draws[c("Z", "Zmo", "Zcs", "Zme")] <- list(named_list(groups))
  for (g in groups) {
    new_r <- subset2(ranef, group = g)
    r_pars <- paste0("^r_", g, usc(usc(p)), "\\[")
    r <- do.call(as.matrix, c(args, list(pars = r_pars)))
    if (is.null(r)) {
      stop2(
        "Group-level effects for each level of group ", 
        "'", g, "' not found. Please set 'save_ranef' to ", 
        "TRUE when fitting your model."
      )
    }
    nlevels <- ngrps(args$x)[[g]]
    old_r <- subset2(args$x$ranef, group = g)
    used_re <- match(new_r$coef, old_r$coef)
    used_re_pars <- outer(seq_len(nlevels), (used_re - 1) * nlevels, "+")
    used_re_pars <- as.vector(used_re_pars)
    r <- r[, used_re_pars, drop = FALSE]
    nranef <- nrow(new_r)
    gtype <- new_r$gtype[1]
    if (gtype == "mm") {
      ngf <- length(new_r$gcall[[1]]$groups)
      gf <- sdata[paste0("J_", new_r$id[1], "_", seq_len(ngf))]
      weights <- sdata[paste0("W_", new_r$id[1], "_", seq_len(ngf))]
    } else {
      gf <- sdata[paste0("J_", new_r$id[1])]
      weights <- list(rep(1, length(gf[[1]])))
    }
    # incorporate new gf levels (only if allow_new_levels is TRUE)
    new_r_levels <- vector("list", length(gf))
    max_level <- nlevels
    for (i in seq_along(gf)) {
      has_new_levels <- any(gf[[i]] > nlevels)
      if (has_new_levels) {
        if (sample_new_levels %in% c("old_levels", "gaussian")) {
          new_levels <- sort(setdiff(gf[[i]], seq_len(nlevels)))
          new_r_levels[[i]] <- matrix(
            nrow = nrow(r), ncol = nranef * length(new_levels)
          )
          if (sample_new_levels == "old_levels") {
            for (j in seq_along(new_levels)) {
              # choose a person to take the group-level effects from
              take_level <- sample(seq_len(nlevels), 1)
              for (k in seq_len(nranef)) {
                take <- (k - 1) * nlevels + take_level
                new_r_levels[[i]][, (j - 1) * nranef + k] <- r[, take]
              }
            }
          } else if (sample_new_levels == "gaussian") {
            # extract hyperparameters used to compute the covariance matrix
            sd_pars <- paste0("sd_", g, usc(usc(p)), "__", new_r$coef)
            sd_samples <- do.call(as.matrix, 
              c(args, list(pars = sd_pars, exact_match = TRUE))
            )
            cor_type <- paste0("cor_", g)
            cor_pars <- paste0(usc(p, "suffix"), new_r$coef)
            cor_pars <- get_cornames(cor_pars, type = cor_type, brackets = FALSE)
            cor_samples <- matrix(
              0, nrow = nrow(sd_samples), ncol = length(cor_pars)
            )
            for (j in seq_along(cor_pars)) {
              if (cor_pars[j] %in% parnames(args$x)) {
                cor_samples[, j] <- do.call(as.matrix, 
                  c(args, list(pars = cor_pars[j], exact_match = TRUE))
                )
              }
            }
            # compute the covariance matrix
            cov_matrix <- get_cov_matrix(sd_samples, cor_samples)$cov
            for (j in seq_along(new_levels)) {
              # sample new levels from the normal distribution
              # implied by the covariance matrix
              indices <- ((j - 1) * nranef + 1):(j * nranef)
              new_r_levels[[i]][, indices] <- t(apply(
                cov_matrix, 1, rmulti_normal, 
                n = 1, mu = rep(0, length(sd_pars))
              ))
            }
          }
          max_level <- max_level + length(new_levels)
        } else if (sample_new_levels == "uncertainty") {
          new_r_levels[[i]] <- matrix(nrow = nrow(r), ncol = nranef)
          for (k in seq_len(nranef)) {
            # sample values for the new level
            indices <- ((k - 1) * nlevels + 1):(k * nlevels)
            new_r_levels[[i]][, k] <- apply(r[, indices], 1, sample, size = 1)
          }
          max_level <- max_level + 1
          gf[[i]][gf[[i]] > nlevels] <- max_level
        }
      } else { 
        new_r_levels[[i]] <- matrix(nrow = nrow(r), ncol = 0)
      }
    }
    new_r_levels <- do.call(cbind, new_r_levels)
    # we need row major instead of column major order
    sort_levels <- ulapply(seq_len(nlevels),
      function(l) seq(l, ncol(r), nlevels)
    )
    r <- cbind(r[, sort_levels, drop = FALSE], new_r_levels)
    levels <- unique(unlist(gf))
    r <- subset_levels(r, levels, nranef)
    # monotonic group-level terms
    new_r_mo <- subset2(new_r, type = "mo")
    if (nrow(new_r_mo)) {
      Z <- matrix(1, length(gf[[1]])) 
      draws[["Zmo"]][[g]] <- prepare_Z(Z, gf, max_level, weights)
      draws[["rmo"]] <- named_list(new_r_mo$coef, list())
      for (co in names(draws[["rmo"]])) {
        take <- which(new_r$coef == co & new_r$type == "mo")
        take <- take + nranef * (seq_along(levels) - 1)
        draws[["rmo"]][[co]][[g]] <- r[, take, drop = FALSE]
      }
    }
    # noise-free group-level terms
    new_r_me <- subset2(new_r, type = "me")
    if (nrow(new_r_me)) {
      Z <- matrix(1, length(gf[[1]]))
      draws[["Zme"]][[g]] <- prepare_Z(Z, gf, max_level, weights)
      draws[["rme"]] <- named_list(new_r_me$coef, list())
      for (co in names(draws[["r_me"]])) {
        take <- which(new_r$coef == co & new_r$type == "me")
        take <- take + nranef * (seq_along(levels) - 1)
        draws[["rme"]][[co]][[g]] <- r[, take, drop = FALSE]
      }
    }
    # category specific group-level terms
    new_r_cs <- subset2(new_r, type = "cs")
    if (nrow(new_r_cs)) {
      Z <- do.call(cbind, lapply(unique(new_r_cs$gn), 
        function(k) sdata[[paste0("Z_", k)]]
      ))
      draws[["Zcs"]][[g]] <- prepare_Z(Z, gf, max_level, weights)
      draws[["rcs"]] <- named_list(seq_len(sdata$ncat - 1), list())
      for (i in names(draws[["rcs"]])) {
        index <- paste0("\\[", i, "\\]$")
        take <- which(grepl(index, new_r$coef) & new_r$type == "cs")
        take <- as.vector(outer(take, nranef * (seq_along(levels) - 1), "+"))
        draws[["rcs"]][[i]][[g]] <- r[, take, drop = FALSE]
      }
    }
    # basic group-level effects
    new_r_basic <- subset2(new_r, type = "")
    if (nrow(new_r_basic)) {
      Z <- do.call(cbind, lapply(unique(new_r_basic$gn), 
        function(k) sdata[[paste0("Z_", k)]]
      ))
      draws[["Z"]][[g]] <- prepare_Z(Z, gf, max_level, weights)
      take <- which(!nzchar(new_r$type))
      take <- as.vector(outer(take, nranef * (seq_along(levels) - 1), "+"))
      r <- r[, take, drop = FALSE]
      draws[["r"]][[g]] <- r
    }
  }
  draws
}

extract_draws_autocor <- function(fit, subset = NULL,
                                  newdata = NULL, sdata = NULL) {
  # extract draws of autocorrelation parameters
  # Args:
  #   sdata: output of make_standata
  stopifnot(is.brmsfit(fit))
  is_newdata <- !is.null(newdata)
  args <- nlist(x = fit, subset)
  draws <- list()
  if (get_ar(fit$autocor)) {
    draws[["ar"]] <- do.call(as.matrix, c(args, pars = "^ar\\["))
  }
  if (get_ma(fit$autocor)) {
    draws[["ma"]] <- do.call(as.matrix, c(args, pars = "^ma\\["))
  }
  if (get_arr(fit$autocor)) {
    draws[["arr"]] <- do.call(as.matrix, c(args, pars = "^arr\\["))
  }
  if (is.cor_sar(fit$autocor)) {
    draws[["lagsar"]] <- do.call(as.matrix, c(args, pars = "^lagsar$"))
    draws[["errorsar"]] <- do.call(as.matrix, c(args, pars = "^errorsar$"))
  }
  if (is.cor_car(fit$autocor)) {
    if (is_newdata) {
      group <- parse_time(fit$autocor)$group
      if (!isTRUE(nzchar(group))) {
        stop2("Without a grouping factor, CAR models cannot handle newdata.")
      }
    }
    gcar <- sdata[["Jloc"]]
    Zcar <- matrix(rep(1, length(gcar)))
    draws[["Zcar"]] <- prepare_Z(Zcar, list(gcar))
    rcar <- do.call(as.matrix, c(args, pars = "^rcar\\["))
    rcar <- rcar[, unique(gcar), drop = FALSE]
    draws[["rcar"]] <- rcar
  }
  if (is.cor_bsts(fit$autocor)) {
    if (is_newdata) {
      warning2("Local level terms are currently ignored ", 
               "when 'newdata' is specified.")
    } else {
      draws[["loclev"]] <- do.call(as.matrix, c(args, pars = "^loclev\\["))
    }
  }
  draws
}

subset_levels <- function(x, levels, nranef) {
  # take relevant cols of a matrix of group-level terms
  # if only a subset of levels is provided (for newdata)
  # Args:
  #   x: a matrix typically samples of r or Z design matrices
  #   levels: grouping factor levels to keep
  #   nranef: number of group-level effects
  take_levels <- ulapply(levels, 
    function(l) ((l - 1) * nranef + 1):(l * nranef)
  )
  x[, take_levels, drop = FALSE]
}

prepare_Z <- function(Z, gf, max_level = max(unlist(gf)),
                      weights = list(rep(1, length(gf[[1]])))) {
  # prepare group-level design matrices for use in linear_predictor
  # Args:
  #   Z: matrix to be prepared
  #   gf: list of vectors containing grouping factor values
  #   max_level: maximal level of gf
  #   weights: optional list of weights of the same length as gf
  nranef <- ncol(Z)
  levels <- unique(unlist(gf))
  Z <- mapply(
    expand_matrix, x = gf, weights = weights,
    MoreArgs = nlist(A = Z, max_level)
  )
  Z <- Reduce("+", Z)
  subset_levels(Z, levels, nranef)
}

expand_matrix <- function(A, x, max_level = max(x), weights = 1) {
  # expand a matrix into a sparse matrix of higher dimension
  # Args:
  #   A: matrix to be expanded
  #   x: levels to expand the matrix
  #   weights: weights to apply to rows of A before expanding
  #   nlevels: number of levels that x can take on
  # Returns:
  #   A sparse matrix of dimension nrow(A) x (ncol(A) * max_level)
  stopifnot(is.matrix(A))
  stopifnot(length(x) == nrow(A))
  stopifnot(all(is_wholenumber(x) & x > 0))
  stopifnot(length(weights) %in% c(1, nrow(A), prod(dim(A))))
  A <- A * weights
  K <- ncol(A)
  i <- rep(seq_along(x), each = K)
  make_j <- function(n, K, x) K * (x[n] - 1) + 1:K
  j <- ulapply(seq_along(x), make_j, K = K, x = x)
  Matrix::sparseMatrix(
    i = i, j = j, x = as.vector(t(A)),
    dims = c(nrow(A), ncol(A) * max_level)
  )
}
