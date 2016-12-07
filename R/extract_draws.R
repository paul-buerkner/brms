extract_draws <- function(x, newdata = NULL, re_formula = NULL, 
                          allow_new_levels = FALSE, incl_autocor = TRUE, 
                          subset = NULL, nsamples = NULL, ...) {
  # extract all data and posterior samples required in (non)linear_predictor
  # Args:
  #   x: an object of class brmsfit
  #   incl_autocor: include autocorrelation parameters in the output?
  #   other arguments: see doc of logLik.brmsfit
  # Returns:
  #   A named list to be understood by linear_predictor.
  #   For non-linear models, every element of draws$nonlinear
  #   can itself be passed to linear_predictor
  ee <- extract_effects(formula(x), family = family(x))
  if (is.null(subset) && !is.null(nsamples)) {
    subset <- sample(Nsamples(x), nsamples)
  }
  nsamples <- Nsamples(x, subset = subset)
  newd_args <- nlist(newdata, re_formula, allow_new_levels, incl_autocor)
  standata <- do.call(amend_newdata, c(newd_args, list(fit = x, ...)))
  draws <- nlist(f = prepare_family(x), data = standata, 
                 nsamples = nsamples, autocor = x$autocor)
  args <- c(newd_args, nlist(x, subset, nsamples, f = draws$f))
  
  # extract draws of auxiliary parameters
  am_args <- list(x = x, subset = subset)
  valid_auxpars <- valid_auxpars(family(x), effects = ee, autocor = x$autocor)
  for (ap in valid_auxpars) {
    if (!is.null(ee[[ap]])) {
      more_args <- list(rhs_formula = attr(x$formula, ap),
                        nlpar = ap, ilink = ilink_auxpars(ap))
      draws[[ap]] <- do.call(.extract_draws, c(args, more_args))
    } else {
      regex <- paste0("^", ap, "($|_)")
      draws[[ap]] <- do.call(as.matrix, c(am_args, pars = regex))
    }
  }
  if (is.linear(family(x)) && length(ee$response) > 1L) {
    # parameters for multivariate models
    draws[["sigma"]] <- do.call(as.matrix, c(am_args, pars = "^sigma($|_)"))
    draws[["rescor"]] <- do.call(as.matrix, c(am_args, pars = "^rescor_"))
    draws[["Sigma"]] <- get_cov_matrix(sd = draws$sigma, cor = draws$rescor)$cov
  }
  if (incl_autocor) {
    # autocorrelation parameters
    if (get_ar(x$autocor)) {
      draws[["ar"]] <- do.call(as.matrix, c(am_args, pars = "^ar\\["))
    }
    if (get_ma(x$autocor)) {
      draws[["ma"]] <- do.call(as.matrix, c(am_args, pars = "^ma\\["))
    }
    if (get_arr(x$autocor)) {
      draws[["arr"]] <- do.call(as.matrix, c(am_args, pars = "^arr\\["))
    }
    if (is(x$autocor, "cor_bsts")) {
      if (is.null(newdata)) {
        draws[["loclev"]] <- do.call(as.matrix, c(am_args, pars = "^loclev\\["))
      } else {
        warning2("Local level terms are currently ignored when ",
                 "'newdata' is specified.")
      }
    }
  }
  # samples of the (non-linear) predictor
  nlpars <- names(ee$nonlinear)
  if (length(nlpars)) {
    for (i in seq_along(nlpars)) {
      rhs_formula <- attr(x$formula, "nonlinear")[[i]]
      more_args <- nlist(rhs_formula, nlpar = nlpars[i])
      draws$nonlinear[[nlpars[i]]] <- 
        do.call(.extract_draws, c(args, more_args))
    }
    covars <- all.vars(rhs(ee$covars))
    if (!all(covars %in% colnames(standata$C))) {
      stop("Covariate matrix is invalid. Please report a bug.")
    }
    for (i in seq_along(covars)) {
      draws[["C"]][[covars[i]]] <- 
        matrix(standata[["C"]][, covars[i]], nrow = nsamples, 
               ncol = nrow(standata[["C"]]), byrow = TRUE)
    }
    draws$nlform <- ee$fixed[[3]]
    # remove redudant information to save working memory
    keep <- !grepl("^(X|Z|J|C)", names(draws$data))
    draws$data <- subset_attr(draws$data, keep)
  } else {
    x$formula <- rm_attr(formula(x), auxpars())
    x$ranef <- tidy_ranef(extract_effects(formula(x)), 
                          data = model.frame(x))
    resp <- ee$response
    if (length(resp) > 1L && !isTRUE(attr(formula(x), "old_mv"))) {
      # new multivariate models
      draws[["mv"]] <- named_list(resp)
      for (r in resp) {
        more_args <- list(nlpar = r)
        draws[["mv"]][[r]] <- do.call(.extract_draws, c(args, more_args))
      }
    } else {
      # univariate models
      draws <- c(draws, do.call(.extract_draws, args))
    }
  }
  draws
}

.extract_draws <- function(x, newdata = NULL, re_formula = NULL,
                           allow_new_levels = FALSE, 
                           incl_autocor = TRUE, subset = NULL,
                           nlpar = "", smooths_only = FALSE,
                           rhs_formula = ~ ., ...) {
  # helper function for extract_draws to extract posterior samples
  # of all regression effects
  # Args:
  #   x: a brmsift object
  #   rhs_formula: used to update the rhs for x$formula
  #                to extract draws of non-linear parameters
  #   smooths_only: extract only smoothing terms?
  #   nlpar: optional name of a non-linear parameter
  #   ...: further elements to store in draws
  # Returns:
  #   a named list
  dots <- list(...)
  nsamples <- Nsamples(x, subset = subset)
  # always update formula to make sure that formulae of
  # non-linear and auxiliary parameters are not included
  # fixes issue #154
  x$formula <- update.formula(formula(x), rhs(rhs_formula))
  x$ranef <- tidy_ranef(extract_effects(formula(x)), data = model.frame(x))
  if (nzchar(nlpar)) {
    # relevant for multivariate models only
    attr(x$formula, "response") <- nlpar 
  }
  new_formula <- update_re_terms(formula(x), re_formula = re_formula)
  ee <- extract_effects(new_formula, family = family(x))
  new_ranef <- tidy_ranef(ee, model.frame(x))
  nlpar_usc <- usc(nlpar, "suffix")
  usc_nlpar <- usc(usc(nlpar))
  newd_args <- nlist(fit = x, newdata, re_formula, 
                     allow_new_levels, incl_autocor)
  draws <- list(data = do.call(amend_newdata, newd_args), 
                old_cat = is.old_categorical(x), ...)
  
  if (smooths_only) {
    # make sure only smooth terms will be included in draws
    keep_elements <- setdiff(names(draws$data), c("Xm", "Xcs", "offset"))
    draws$data <- draws$data[keep_elements]
    new_ranef <- empty_ranef()
  }
  
  args <- list(x = x, subset = subset)
  if (isTRUE(ncol(draws$data[["X"]]) > 0) && draws$old_cat != 1L) {
    b_pars <- paste0("b_", nlpar_usc, colnames(draws$data[["X"]]))
    draws[["b"]] <- 
      do.call(as.matrix, c(args, list(pars = b_pars, exact = TRUE)))
  }
  # monotonic effects
  if (isTRUE(ncol(draws$data[["Xmo"]]) > 0)) {
    monef <- colnames(draws$data[["Xmo"]])
    draws[["bmo"]] <- draws$simplex <- named_list(monef)
    # as of brms > 1.2.0 the prefix 'bmo' is used
    bmo <- ifelse(any(grepl("^bmo_", parnames(x))), "bmo_",
                  ifelse(any(grepl("^bm_", parnames(x))), "bm_", "b_"))
    for (i in seq_along(monef)) {
      bmo_par <- paste0(bmo, nlpar_usc, monef[i])
      draws[["bmo"]][[i]] <- 
        do.call(as.matrix, c(args, list(pars = bmo_par, exact = TRUE)))
      simplex_par <- paste0("simplex_", nlpar_usc, monef[i], 
                            "[", seq_len(draws$data$Jm[i]), "]")
      draws[["simplex"]][[i]] <- 
        do.call(as.matrix, c(args, list(pars = simplex_par, exact = TRUE)))
      
    }
  }
  # noise-free effects
  meef <- get_me_labels(ee, x$data)
  if (length(meef)) {
    if (!is.null(newdata)) {
      stop2("Predictions with noise-free variables are not yet ",
            "possible when passing new data.")
    }
    if (!any(grepl(paste0("Xme_", nlpar_usc), parnames(x)))) {
      stop2("Noise-free variables were not saved. Please set ",
            "argument 'save_mevars' to TRUE when calling 'brm'.")
    }
    uni_me <- attr(meef, "uni_me")
    not_one <- attr(meef, "not_one")
    # prepare calls to evaluate noise-free data
    me_sp <- strsplit(gsub("[[:space:]]", "", meef), ":")
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
      bme_par <- paste0("bme_", nlpar_usc, meef[j])
      draws[["bme"]][[j]] <- 
        do.call(as.matrix, c(args, list(pars = bme_par, exact = TRUE)))
    }
    # extract noise-free variable samples
    uni_me <- rename(uni_me)
    draws[["Xme"]] <- named_list(uni_me)
    for (j in seq_along(draws[["Xme"]])) {
      Xme_pars <- paste0("Xme_", nlpar_usc, uni_me[j],
                         "[", seq_len(nobs(x)), "]")
      draws[["Xme"]][[j]] <- 
        do.call(as.matrix, c(args, list(pars = Xme_pars, exact = TRUE)))
    }
    # prepare covariates
    ncovars <- sum(not_one)
    for (j in seq_len(ncovars)) {
      cme <- paste0("Cme_", j)
      draws[["Cme"]][[j]] <- 
        matrix(draws$data[[cme]], nrow = nsamples, 
               ncol = length(draws$data[[cme]]), byrow = TRUE)
    }
  }
  # category specific effects 
  if (is.ordinal(family(x))) {
    draws[["Intercept"]] <- 
      do.call(as.matrix, c(args, list(pars = "^b_Intercept\\[")))
    if (isTRUE(ncol(draws$data[["Xcs"]]) > 0)) {
      # as of brms > 1.0.1 the original prefix 'bcs' is used
      bcs <- ifelse(any(grepl("^bcs_", parnames(x))), "^bcs_", "^b_")
      cs_pars <- paste0(bcs, colnames(draws$data$Xcs), "\\[")
      draws[["cs"]] <- do.call(as.matrix, c(args, list(pars = cs_pars)))
    }
  } else if (draws$old_cat == 1L) {
    # old categorical models deprecated as of brms > 0.8.0
    if (!is.null(draws$data[["X"]])) {
      draws[["cs"]] <- do.call(as.matrix, c(args, list(pars = "^b_")))
    }
  }
  # splines
  splines <- rename(get_spline_labels(ee, x$data, covars = TRUE))
  if (length(splines)) {
    draws[["Zs"]] <- draws[["s"]] <- named_list(splines)
    for (i in seq_along(splines)) {
      nb <- seq_len(attr(splines, "nbases")[[i]])
      for (j in nb) {
        draws[["Zs"]][[splines[i]]][[j]] <- 
          draws$data[[paste0("Zs_", i, "_", j)]]
        s_pars <- paste0("^s_", nlpar_usc, splines[i], "_", j, "\\[")
        draws[["s"]][[splines[i]]][[j]] <- 
          do.call(as.matrix, c(args, list(pars = s_pars)))
      }
    }
  }
  
  # requires initialization to assign S4 objects of the Matrix package
  groups <- unique(new_ranef$group) 
  draws[["Z"]] <- draws[["Zmo"]] <- 
    draws[["Zcs"]] <- draws[["Zme"]] <- named_list(groups)
  for (g in groups) {
    new_r <- new_ranef[new_ranef$group == g, ]
    gf <- get(paste0("J_", new_r$id[1]), draws$data)
    r_pars <- paste0("^r_", g, usc_nlpar, "\\[")
    r <- do.call(as.matrix, c(args, list(pars = r_pars)))
    if (is.null(r)) {
      stop2("Group-level effects for each level of group ", 
            "'", g, "' not found. Please set ranef = TRUE ",
            "when calling brm.")
    }
    nlevels <- ngrps(x)[[g]]
    old_r <- x$ranef[x$ranef$group == g, ]
    used_re <- match(new_r$coef, old_r$coef)
    used_re_pars <- outer(seq_len(nlevels), (used_re - 1) * nlevels, "+")
    used_re_pars <- as.vector(used_re_pars)
    r <- r[, used_re_pars, drop = FALSE]
    nranef <- nrow(new_r)
    # incorporate new gf levels (only if allow_new_levels is TRUE)
    has_new_levels <- anyNA(gf)
    if (has_new_levels) {
      new_r_level <- matrix(nrow = nrow(r), ncol = nranef)
      for (k in seq_len(nranef)) {
        # sample values of the new level for each group-level effect
        indices <- ((k - 1) * nlevels + 1):(k * nlevels)
        new_r_level[, k] <- apply(r[, indices], 1, sample, size = 1)
      }
      gf[is.na(gf)] <- nlevels + 1
    } else { 
      new_r_level <- matrix(nrow = nrow(r), ncol = 0)
    }
    # we need row major instead of column major order
    sort_levels <- ulapply(seq_len(nlevels), 
      function(l) seq(l, ncol(r), nlevels))
    r <- cbind(r[, sort_levels, drop = FALSE], new_r_level)
    levels <- unique(gf)
    r <- subset_levels(r, levels, nranef)
    # monotonic group-level terms
    new_r_mo <- new_r[new_r$type == "mo", ]
    if (nrow(new_r_mo)) {
      Z_mo <- expand_matrix(matrix(1, length(gf)), gf)
      if (length(levels) < nlevels) {
        Z_mo <- Z_mo[, levels, drop = FALSE]
      }
      draws[["Zmo"]][[g]] <- Z_mo
      draws[["rmo"]] <- named_list(new_r_mo$coef, list())
      for (co in names(draws[["rmo"]])) {
        take <- which(new_r$coef == co & new_r$type == "mo")
        take <- take + nranef * (seq_along(levels) - 1)
        draws[["rmo"]][[co]][[g]] <- r[, take, drop = FALSE]
      }
    }
    # noise-free group-level terms
    new_r_me <- new_r[new_r$type == "me", ]
    if (nrow(new_r_me)) {
      Z_me <- expand_matrix(matrix(1, length(gf)), gf)
      if (length(levels) < nlevels) {
        Z_me <- Z_me[, levels, drop = FALSE]
      }
      draws[["Zme"]][[g]] <- Z_me
      draws[["rme"]] <- named_list(new_r_me$coef, list())
      for (co in names(draws[["r_me"]])) {
        take <- which(new_r$coef == co & new_r$type == "me")
        take <- take + nranef * (seq_along(levels) - 1)
        draws[["rme"]][[co]][[g]] <- r[, take, drop = FALSE]
      }
    }
    # category specific group-level terms
    new_r_cs <- new_r[new_r$type == "cs", ]
    if (nrow(new_r_cs)) {
      Z <- prepare_Z(new_r_cs, gf, draws$data)
      draws[["Zcs"]][[g]] <- Z
      draws[["rcs"]] <- named_list(seq_len(draws$data$ncat - 1), list())
      for (i in names(draws[["rcs"]])) {
        index <- paste0("\\[", i, "\\]$")
        take <- which(grepl(index, new_r$coef) & new_r$type == "cs")
        take <- as.vector(outer(take, nranef * (seq_along(levels) - 1), "+"))
        draws[["rcs"]][[i]][[g]] <- r[, take, drop = FALSE]
      }
    }
    # basic group-level effects
    new_r_basic <- new_r[!nzchar(new_r$type), ]
    if (nrow(new_r_basic)) {
      Z <- prepare_Z(new_r_basic, gf, draws$data)
      draws[["Z"]][[g]] <- Z
      take <- which(!nzchar(new_r$type))
      take <- as.vector(outer(take, nranef * (seq_along(levels) - 1), "+"))
      r <- r[, take, drop = FALSE]
      draws[["r"]][[g]] <- r
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
    function(l) ((l - 1) * nranef + 1):(l * nranef))
  x[, take_levels, drop = FALSE]
}

prepare_Z <- function(ranef, gf, standata) {
  # prepare group-level design matrices for use in linear_predictor
  # Args:
  #   ranef: output of tidy_ranef
  #   standata: output of make_standata
  #   gf: values of the grouping factor
  Z <- lapply(unique(ranef$gn), 
              function(k) standata[[paste0("Z_", k)]])
  Z <- do.call(cbind, Z)
  nranef <- ncol(Z)
  Z <- expand_matrix(Z, gf)
  subset_levels(Z, unique(gf), nranef)
}

expand_matrix <- function(A, x) {
  # expand a matrix into a sparse matrix of higher dimension
  # Args:
  #   A: matrix to be expanded
  #   x: levels to expand the matrix
  # Notes:
  #   used in extract_draws
  # Returns:
  #   A sparse matrix of dimension nrow(A) x (ncol(A) * length(x))
  stopifnot(is.matrix(A))
  stopifnot(length(x) == nrow(A))
  stopifnot(all(is.wholenumber(x) & x > 0))
  K <- ncol(A)
  i <- rep(seq_along(x), each = K)
  make_j <- function(n, K, x) K * (x[n] - 1) + 1:K
  j <- ulapply(seq_along(x), make_j, K = K, x = x)
  Matrix::sparseMatrix(i = i, j = j, x = as.vector(t(A)))
}
