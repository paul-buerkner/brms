extract_draws <- function(x, ...) {
  # extract data and posterior draws
  UseMethod("extract_draws")
}

#' @export
extract_draws.brmsfit <- function(x, newdata = NULL, re_formula = NULL, 
                                  sample_new_levels = "uncertainty",
                                  incl_autocor = TRUE, resp = NULL,
                                  subset = NULL, nsamples = NULL, nug = NULL, 
                                  smooths_only = FALSE, ...) {
  # extract all data and posterior draws required in (non)linear_predictor
  # Args:
  #   see doc of logLik.brmsfit
  #   ...: passed to validate_newdata
  # Returns:
  #   A named list to be interpreted by linear_predictor
  snl_options <- c("uncertainty", "gaussian", "old_levels")
  sample_new_levels <- match.arg(sample_new_levels, snl_options)
  x <- restructure(x)
  if (!incl_autocor) {
    x <- remove_autocor(x) 
  }
  sdata <- validate_newdata(
    newdata, x, re_formula = re_formula, resp = resp, ...
  )
  subset <- subset_samples(x, subset, nsamples)
  samples <- as.matrix(x, subset = subset)
  new_formula <- update_re_terms(x$formula, re_formula)
  bterms <- parse_bf(new_formula)
  ranef <- tidy_ranef(bterms, data = x$data)
  args <- nlist(
    x = bterms, samples, sdata, data = x$data,
    ranef, old_ranef = x$ranef, sample_new_levels,
    resp, nug, smooths_only, new = !is.null(newdata)
  )
  if (length(get_effect(bterms, "gp")) && !is.null(newdata)) {
    # GPs for new data require the original data as well
    args$old_sdata <- validate_newdata(
      newdata = NULL, fit = x, re_formula = re_formula, ...
    )
  }
  do.call(extract_draws, args)
}

extract_draws.mvbrmsterms <- function(x, samples, sdata, resp = NULL, ...) {
  resp <- validate_resp(resp, x$responses)
  if (length(resp) > 1) {
    draws <- list(nsamples = nrow(samples), nobs = sdata$N)
    draws$resps <- named_list(resp)
    draws$old_order <- attr(sdata, "old_order")
    for (r in resp) {
      draws$resps[[r]] <- extract_draws(
        x$terms[[r]], samples = samples, sdata = sdata, ...
      )
    }
    if (x$rescor) {
      draws$f <- draws$resps[[1]]$f
      draws$f$fun <- paste0(draws$f$family, "_mv")
      rescor <- get_cornames(resp, type = "rescor", brackets = FALSE)
      draws$mvpars$rescor <- get_samples(samples, rescor, exact = TRUE)
      if (draws$f$family == "student") {
        # store in draws$dpars so that get_dpar can be called on nu
        draws$dpars$nu <- as.vector(get_samples(samples, "^nu$"))
      }
      draws$data$N <- draws$resps[[1]]$data$N
      draws$data$weights <- draws$resps[[1]]$data$weights
      Y <- lapply(draws$resps, function(x) x$data$Y)
      draws$data$Y <- do.call(cbind, Y)
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
extract_draws.brmsterms <- function(x, samples, sdata, ...) {
  nsamples <- nrow(samples)
  nobs <- sdata$N
  resp <- usc(combine_prefix(x))
  draws <- nlist(f = prepare_family(x), nsamples, nobs, resp = x$resp)
  draws$old_order <- attr(sdata, "old_order")
  valid_dpars <- valid_dpars(x)
  draws$dpars <- named_list(valid_dpars)
  for (dp in valid_dpars) {
    dp_regex <- paste0("^", dp, resp, "$")
    if (is.btl(x$dpars[[dp]]) || is.btnl(x$dpars[[dp]])) {
      draws$dpars[[dp]] <- extract_draws(
        x$dpars[[dp]], samples = samples, sdata = sdata, ...
      )
    } else if (is.numeric(x$fdpars[[dp]]$value)) {
      draws$dpars[[dp]] <- x$fdpars[[dp]]$value
    } else if (any(grepl(dp_regex, colnames(samples)))) {
      draws$dpars[[dp]] <- as.vector(get_samples(samples, dp_regex))
    }
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
      draws$dpars$theta <- do.call(cbind, draws$dpars[thetas])
      draws$dpars[thetas] <- NULL
      if (nrow(draws$dpars$theta) == 1L) {
        dim <- c(nrow(samples), ncol(draws$dpars$theta))
        draws$dpars$theta <- as_draws_matrix(draws$dpars$theta, dim = dim)
      }
    }
  }
  if (use_cov(x$autocor) || is.cor_sar(x$autocor)) {
    # only include autocor samples on the top-level of draws 
    # when using the covariance formulation of ARMA / SAR structures
    draws$ac <- extract_draws_autocor(x, samples, sdata, ...)
  }
  draws$data <- extract_draws_data(x, sdata = sdata, ...)
  structure(draws, class = "brmsdraws")
}

#' @export
extract_draws.btnl <- function(x, samples, sdata, ...) {
  # Args:
  #   C: matrix containing covariates
  #   ...: passed to extract_draws.btl
  draws <- list(f = x$family, nsamples = nrow(samples), nobs = sdata$N)
  nlpars <- names(x$nlpars)
  for (nlp in nlpars) {
    draws$nlpars[[nlp]] <- extract_draws(
      x$nlpars[[nlp]], samples = samples, sdata = sdata, ...
    )
  }
  p <- usc(combine_prefix(x))
  C <- sdata[[paste0("C", p)]]
  stopifnot(is.matrix(C))
  for (cov in colnames(C)) {
    draws$C[[cov]] <- as_draws_matrix(
      C[, cov], dim = c(nrow(samples), nrow(C))
    )
  }
  draws$nlform <- x$formula[[2]]
  draws$f <- x$family
  draws
}

#' @export
extract_draws.btl <- function(x, samples, sdata, smooths_only = FALSE, ...) {
  # extract draws of all kinds of effects
  # Args:
  #   fit: a brmsfit object
  #   mv: is the model multivariate?
  #   ...: further elements to store in draws
  #   other arguments: see extract_draws.brmsfit
  # Returns:
  #   A named list to be interpreted by linear_predictor
  nsamples <- nrow(samples)
  draws <- nlist(f = x$family, nsamples, nobs = sdata$N)
  args <- nlist(bterms = x, samples, sdata, ...)
  draws$fe <- do.call(extract_draws_fe, args)
  draws$sm <- do.call(extract_draws_sm, args)
  if (smooths_only) {
    # make sure only smooth terms will be included in draws
    return(structure(draws, predicted = TRUE))
  }
  draws$sp <- do.call(extract_draws_sp, args)
  draws$cs <- do.call(extract_draws_cs, args)
  draws$gp <- do.call(extract_draws_gp, args)
  draws$re <- do.call(extract_draws_re, args)
  draws$offset <- do.call(extract_draws_offset, args)
  if (!(use_cov(x$autocor) || is.cor_sar(x$autocor))) {
    draws$ac <- do.call(extract_draws_autocor, args)
  }
  structure(draws, predicted = TRUE)
}

extract_draws_fe <- function(bterms, samples, sdata, ...) {
  # extract draws of ordinary population-level effects
  # Args:
  #   fixef: names of the population-level effects
  #   args: list of arguments passed to as.matrix.brmsfit
  #   nlpar: name of a non-linear parameter
  #   old_cat: see old_categorical
  # Returns: 
  #   A named list to be interpreted by linear_predictor
  # stopifnot("x" %in% names(args))
  draws <- list()
  p <- usc(combine_prefix(bterms))
  draws$X <- sdata[[paste0("X", p)]]
  fixef <- colnames(draws$X)
  if (length(fixef)) {
    b_pars <- paste0("b", p, "_", fixef)
    draws$b <- get_samples(samples, b_pars, exact = TRUE)
  }
  draws
}

extract_draws_sp <- function(bterms, samples, sdata, data, new = FALSE, ...) {
  # extract draws of special effects terms
  # Returns: 
  #   A named list to be interpreted by linear_predictor
  draws <- list()
  spef <- tidy_spef(bterms, data)
  if (is.null(spef)) {
    return(draws) 
  }
  p <- usc(combine_prefix(bterms))
  # prepare calls evaluated in sp_predictor
  draws$calls <- vector("list", nrow(spef))
  for (i in seq_along(draws$calls)) {
    call <- spef$call_prod[[i]]
    if (!is.null(spef$call_mo[[i]])) {
      new_mo <- paste0(
        ".mo(simo_", spef$Imo[[i]], ", Xmo_", spef$Imo[[i]], ")"
      )
      call <- rename(call, spef$call_mo[[i]], new_mo)
    }
    if (!is.null(spef$call_me[[i]])) {
      new_me <- paste0("Xme_", seq_along(spef$uni_me[[i]]))
      call <- rename(call, spef$uni_me[[i]], new_me)
    }
    if (!is.null(spef$call_mi[[i]])) {
      new_mi <- paste0("Yf_", spef$vars_mi[[i]])
      call <- rename(call, spef$call_mi[[i]], new_mi)
    }
    if (spef$Ic[i] > 0) {
      str_add(call) <- paste0(" * Csp_", spef$Ic[i])
    }
    draws$calls[[i]] <- parse(text = paste0(call))
  }
  # extract general data and parameters for special effects
  bsp_pars <- paste0("bsp", p, "_", spef$coef)
  draws$bsp <- get_samples(samples, bsp_pars, exact = TRUE)
  # prepare draws specific to monotonic effects
  simo_coef <- get_simo_labels(spef)
  draws$simo <- draws$Xmo <- named_list(simo_coef)
  for (i in seq_along(simo_coef)) {
    J <- seq_len(sdata$Jmo[i])
    simo_par <- paste0("simo", p, "_", simo_coef[i], "[", J, "]")
    draws$simo[[i]] <- get_samples(samples, simo_par, exact = TRUE)
    draws$Xmo[[i]] <- sdata[[paste0("Xmo", p, "_", i)]]
  }
  # prepare draws specific to noise-free effects
  uni_me <- rename(unique(unlist(spef$uni_me)))
  if (length(uni_me)) {
    if (new) {
      stop2("Predictions with noise-free variables are not yet ",
            "possible when passing new data.")
    }
    if (!any(grepl("^Xme_", colnames(samples)))) {
      stop2("Noise-free variables were not saved. Please set ",
            "argument 'save_mevars' to TRUE when fitting your model.")
    }
    draws$Xme <- named_list(uni_me)
    for (i in seq_along(draws$Xme)) {
      Xme_pars <- paste0("Xme", p, "_", uni_me[i], "\\[")
      draws$Xme[[i]] <- get_samples(samples, Xme_pars)
    }
  }
  # TODO: prepare draws specific to missing value effects
  ncovars <- max(spef$Ic)
  dim <- c(nrow(draws$bsp), sdata$N)
  for (i in seq_len(ncovars)) {
    draws$Csp[[i]] <- as_draws_matrix(
      sdata[[paste0("Csp", p, "_", i)]], dim = dim
    )
  }
  draws
}

extract_draws_cs <- function(bterms, samples, sdata, data, ...) {
  # category specific effects 
  # extract draws of category specific effects
  # Args:
  #   csef: names of the category specific effects
  #   args: list of arguments passed to as.matrix.brmsfit
  #   nlpar: name of a non-linear parameter
  #   old_cat: see old_categorical
  # Returns: 
  #   A named list to be interpreted by linear_predictor
  draws <- list()
  p <- usc(combine_prefix(bterms))
  int_regex <- paste0("^b", p, "_Intercept\\[")
  is_ordinal <- any(grepl(int_regex, colnames(samples))) 
  if (is_ordinal) {
    draws$ncat <- sdata[[paste0("ncat", p)]]
    draws$Intercept <- get_samples(samples, int_regex)
    csef <- colnames(get_model_matrix(bterms$cs, data))
    if (length(csef)) {
      # as of brms > 1.0.1 the original prefix 'bcs' is used
      bcs <- ifelse(any(grepl("^bcs_", colnames(samples))), "^bcs", "^b")
      cs_pars <- paste0(bcs, p, "_", csef, "\\[")
      draws$bcs <- get_samples(samples, cs_pars)
      draws$Xcs <- sdata[[paste0("Xcs", p)]]
    }
  }
  draws
}

extract_draws_sm <- function(bterms, samples, sdata, data, ...) {
  # extract draws of smooth terms
  # Args:
  #   smooths: names of the smooth terms
  #   args: list of arguments passed to as.matrix.brmsfit
  #   sdata: list returned by make_standata
  #   nlpar: name of a non-linear parameter
  # Returns: 
  #   A named list to be interpreted by linear_predictor
  smooths <- get_sm_labels(bterms, data, covars = TRUE)
  if (!length(smooths)) {
    return(list())
  }
  p <- usc(combine_prefix(bterms))
  draws <- named_list(smooths)
  for (i in seq_along(smooths)) {
    sm <- list()
    nb <- seq_len(attr(smooths, "nbases")[[i]])
    for (j in nb) {
      sm$Zs[[j]] <- sdata[[paste0("Zs", p, "_", i, "_", j)]]
      s_pars <- paste0("^s", p, "_", smooths[i], "_", j, "\\[")
      sm$s[[j]] <- get_samples(samples, s_pars)
    }
    draws[[i]] <- sm
  }
  draws
}

extract_draws_gp <- function(bterms, samples, sdata, data,
                             new = FALSE, nug = NULL,
                             old_sdata = NULL, ...) {
  # extract draws for gaussian processes
  # Args:
  #   gpef: names of the gaussian process terms
  gpef <- get_gp_labels(bterms, data, covars = TRUE)
  if (!length(gpef)) {
    return(list()) 
  }
  if (new) {
    stopifnot(!is.null(old_sdata))
  }
  p <- usc(combine_prefix(bterms))
  if (is.null(nug)) {
    nug <- ifelse(new, 1e-8, 1e-11)
  }
  draws <- named_list(gpef)
  for (i in seq_along(gpef)) {
    gp <- list()
    by_levels <- attr(gpef, "by_levels")[[i]]
    gp_names <- paste0(gpef[i], usc(by_levels))
    sdgp <- paste0("^sdgp", p, "_", gp_names)
    gp$sdgp <- get_samples(samples, sdgp)
    lscale <- paste0("^lscale", p, "_", gp_names)
    gp$lscale <- get_samples(samples, lscale)
    zgp <- paste0("^zgp", p, "_", gpef[i], "\\[")
    gp$zgp <- get_samples(samples, zgp)
    gp$bynum <- sdata[[paste0("Cgp", p, "_", i)]]
    Jgp_pos <- grepl(paste0("^Jgp", p, "_", i), names(sdata))
    if (new) {
      gp$x <- old_sdata[[paste0("Xgp", p, "_", i)]]
      gp$x_new <- sdata[[paste0("Xgp", p, "_", i)]]
      gp$Jgp <- old_sdata[Jgp_pos]
      gp$Jgp_new <- sdata[Jgp_pos]
      # computing GPs for new data requires the old GP terms
      gp$yL <- .predictor_gp(
        x = gp[["x"]], sdgp = gp[["sdgp"]],
        lscale = gp[["lscale"]], zgp = gp[["zgp"]], 
        Jgp = gp[["Jgp"]]
      )
    } else {
      gp$x <- sdata[[paste0("Xgp", p, "_", i)]]
      gp$Jgp <- sdata[Jgp_pos]
    }
    gp$nug <- nug
    draws[[i]] <- gp
  }
  draws
}

extract_draws_re <- function(bterms, samples, sdata, data, ranef, old_ranef, 
                             sample_new_levels = "uncertainty", ...) {
  # extract draws of group-level effects
  # Args:,
  #   ranef: data.frame returned by tidy_ranef
  #   args: list of arguments passed to as.matrix.brmsfit
  #   sdata: list returned by make_standata
  #   sample_new_levels: see help("predict.brmsfit")
  # Returns: 
  #   A named list to be interpreted by linear_predictor
  draws <- list()
  px <- check_prefix(bterms)
  ranef <- subset2(ranef, ls = px)
  if (!nrow(ranef)) {
    return(draws)
  }
  p <- combine_prefix(px)
  groups <- unique(ranef$group)
  old_ranef <- subset2(old_ranef, ls = px)
  # assigning S4 objects requires initialisation of list elements
  draws[c("Z", "Zsp", "Zcs")] <- list(named_list(groups))
  for (g in groups) {
    new_r <- subset2(ranef, group = g)
    r_pars <- paste0("^r_", g, usc(usc(p)), "\\[")
    r <- get_samples(samples, r_pars)
    if (is.null(r)) {
      stop2(
        "Group-level effects for each level of group ", 
        "'", g, "' not found. Please set 'save_ranef' to ", 
        "TRUE when fitting your model."
      )
    }
    nlevels <- length(attr(old_ranef, "levels")[[g]]) 
    old_r <- subset2(old_ranef, group = g)
    used_re <- match(new_r$coef, old_r$coef)
    used_re_pars <- outer(seq_len(nlevels), (used_re - 1) * nlevels, "+")
    used_re_pars <- as.vector(used_re_pars)
    r <- r[, used_re_pars, drop = FALSE]
    nranef <- nrow(new_r)
    gtype <- new_r$gtype[1]
    id <- new_r$id[1]
    if (gtype == "mm") {
      ngf <- length(new_r$gcall[[1]]$groups)
      gf <- sdata[paste0("J_", id, "_", seq_len(ngf))]
      weights <- sdata[paste0("W_", id, "_", seq_len(ngf))]
    } else {
      gf <- sdata[paste0("J_", id)]
      weights <- list(rep(1, length(gf[[1]])))
    }
    # incorporate new gf levels
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
            sd_samples <- get_samples(samples, sd_pars, exact = TRUE)
            cor_type <- paste0("cor_", g)
            cor_pars <- paste0(usc(p, "suffix"), new_r$coef)
            cor_pars <- get_cornames(cor_pars, type = cor_type, brackets = FALSE)
            cor_samples <- matrix(
              0, nrow = nrow(sd_samples), ncol = length(cor_pars)
            )
            for (j in seq_along(cor_pars)) {
              if (cor_pars[j] %in% colnames(samples)) {
                cor_samples[, j] <- get_samples(
                  samples, cor_pars[j], exact = TRUE
                )
              }
            }
            # compute the covariance matrix
            cov_matrix <- get_cov_matrix(sd_samples, cor_samples)
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
            new_r_levels[[i]][, k] <- apply(
              r[, indices, drop = FALSE], 1, sample, size = 1
            )
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
    # special group-level terms (mo, me, mi)
    new_r_sp <- subset2(new_r, type = "sp")
    if (nrow(new_r_sp)) {
      Z <- matrix(1, length(gf[[1]])) 
      draws[["Zsp"]][[g]] <- prepare_Z(Z, gf, max_level, weights)
      for (co in new_r_sp$coef) {
        take <- which(new_r$coef == co & new_r$type == "sp")
        take <- take + nranef * (seq_along(levels) - 1)
        draws[["rsp"]][[co]][[g]] <- r[, take, drop = FALSE]
      }
    }
    # category specific group-level terms
    new_r_cs <- subset2(new_r, type = "cs")
    if (nrow(new_r_cs)) {
      # all categories share the same Z matrix
      cn1 <- new_r_cs$cn[grepl("\\[1\\]$", new_r_cs$coef)]
      Znames <- paste0("Z_", id, usc(p), "_", cn1)
      Z <- do.call(cbind, sdata[Znames])
      draws[["Zcs"]][[g]] <- prepare_Z(Z, gf, max_level, weights)
      for (i in seq_len(sdata$ncat - 1)) {
        index <- paste0("\\[", i, "\\]$")
        take <- which(grepl(index, new_r$coef) & new_r$type == "cs")
        take <- as.vector(outer(take, nranef * (seq_along(levels) - 1), "+"))
        draws[["rcs"]][[g]][[i]] <- r[, take, drop = FALSE]
      }
    }
    # basic group-level effects
    new_r_basic <- subset2(new_r, type = "")
    if (nrow(new_r_basic)) {
      Znames <- paste0("Z_", id, usc(p), "_", new_r_basic$cn)
      Z <- do.call(cbind, sdata[Znames])
      draws[["Z"]][[g]] <- prepare_Z(Z, gf, max_level, weights)
      take <- which(!nzchar(new_r$type))
      take <- as.vector(outer(take, nranef * (seq_along(levels) - 1), "+"))
      r <- r[, take, drop = FALSE]
      draws[["r"]][[g]] <- r
    }
  }
  draws
}

extract_draws_offset <- function(bterms, samples, sdata, ...) {
  p <- usc(combine_prefix(bterms))
  sdata[[paste0("offset", p)]]
}

extract_draws_autocor <- function(bterms, samples, sdata, new = FALSE, ...) {
  # extract draws of autocorrelation parameters
  draws <- list()
  autocor <- bterms$autocor
  p <- usc(combine_prefix(bterms))
  draws$N_tg <- sdata[[paste0("N_tg", p)]]
  if (get_ar(autocor) || get_ma(autocor)) {
    draws$Y <- sdata[[paste0("Y", p)]]
    draws$J_lag <- sdata[[paste0("J_lag", p)]]
    if (get_ar(autocor)) {
      draws$ar <- get_samples(samples, paste0("^ar", p, "\\["))
    }
    if (get_ma(autocor)) {
      draws$ma <- get_samples(samples, paste0("^ma", p, "\\["))
    }
    if (use_cov(autocor)) {
      draws$begin_tg <- sdata[[paste0("begin_tg", p)]]
      draws$nobs_tg <- sdata[[paste0("nobs_tg", p)]]
    }
  }
  if (get_arr(autocor)) {
    draws$arr <- get_samples(samples, paste0("^arr", p, "\\["))
    draws$Yarr <- sdata[[paste0("Yarr", p)]]
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
  if (is.cor_bsts(autocor)) {
    if (new) {
      warning2(
        "Local level terms are currently ignored ", 
        "when 'newdata' is specified."
      )
    } else {
      draws$loclev <- get_samples(samples, paste0("^loclev", p, "\\["))
    }
  }
  draws
}

extract_draws_data <- function(bterms, sdata, ...) {
  # extract data mainly related to the response variable
  vars <- c(
    "Y", "trials", "ncat", "se", "weights", 
    "dec", "cens", "rcens", "lb", "ub"
  )
  resp <- usc(combine_prefix(bterms))
  draws <- sdata[paste0(vars, resp)]
  rmNULL(draws, recursive = FALSE)
}

pseudo_draws_for_mixture <- function(draws, comp, sample_ids = NULL) {
  # create pseudo brmsdraws objects for components of mixture models
  # Args:
  #   comp: the mixture component number
  #   sample_ids: see predict_mixture
  stopifnot(is.brmsdraws(draws), is.mixfamily(draws$f))
  if (!is.null(sample_ids)) {
    nsamples <- length(sample_ids)
  } else {
    nsamples <- draws$nsamples
  }
  out <- list(
    f = draws$f$mix[[comp]], nsamples = nsamples,
    nobs = draws$nobs, data = draws$data
  )
  out$f$fun <- out$f$family
  for (dp in valid_dpars(out$f)) {
    out$dpars[[dp]] <- draws$dpars[[paste0(dp, comp)]]
    if (!is.null(sample_ids)) {
      out$dpars[[dp]] <- p(out$dpars[[dp]], sample_ids, row = TRUE)
    }
  }
  structure(out, class = "brmsdraws")
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

get_samples <- function(x, pars, ...) {
  pars <- extract_pars(pars, all_pars = colnames(x), ...)
  x[, pars, drop = FALSE]
}

is.brmsdraws <- function(x) {
  inherits(x, "brmsdraws")
}

is.mvbrmsdraws <- function(x) {
  inherits(x, "mvbrmsdraws")
}

