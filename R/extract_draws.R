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
  standata <- amend_newdata(newdata, fit = x, re_formula = re_formula,
                            allow_new_levels = allow_new_levels, ...)
  draws <- nlist(f = prepare_family(x), data = standata, 
                 nsamples = nsamples, autocor = x$autocor)
  
  # extract draws of auxiliary parameters
  args <- list(x = x, subset = subset)
  ilinks <- c(sigma = "exp", shape = "exp", nu = "exp", phi = "exp", 
              kappa = "exp", zi = "inv_logit", hu = "inv_logit") 
  valid_auxpars <- valid_auxpars(family(x), effects = ee, autocor = x$autocor)
  for (ap in valid_auxpars) {
    if (!is.null(ee[[ap]])) {
      auxfit <- x
      auxfit$formula <- update.formula(x$formula, rhs(attr(x$formula, ap)))
      auxfit$ranef <- tidy_ranef(extract_effects(auxfit$formula), 
                                 data = model.frame(x))
      auxstandata <- amend_newdata(newdata, fit = auxfit, 
                                   re_formula = re_formula,
                                   allow_new_levels = allow_new_levels)
      draws[[ap]] <- .extract_draws(auxfit, standata = auxstandata, nlpar = ap,
                                    re_formula = re_formula, subset = subset)
      rm(auxfit)
      draws[[ap]][c("f", "data", "nsamples", "ilink")] <- 
        list(draws$f, auxstandata, draws$nsamples, ilinks[ap])
    } else {
      regex <- paste0("^", ap, "($|_)")
      draws[[ap]] <- do.call(as.matrix, c(args, pars = regex))
    }
  }
  # parameters for multivariate models
  if (is.linear(family(x)) && length(ee$response) > 1L) {
    draws[["sigma"]] <- do.call(as.matrix, c(args, pars = "^sigma($|_)"))
    draws[["rescor"]] <- do.call(as.matrix, c(args, pars = "^rescor_"))
    draws[["Sigma"]] <- get_cov_matrix(sd = draws$sigma, cor = draws$rescor)$cov
  }
  # autocorrelation parameters
  if (incl_autocor) {
    if (get_ar(x$autocor)) {
      draws[["ar"]] <- do.call(as.matrix, c(args, pars = "^ar\\["))
    }
    if (get_ma(x$autocor)) {
      draws[["ma"]] <- do.call(as.matrix, c(args, pars = "^ma\\["))
    }
    if (get_arr(x$autocor)) {
      draws[["arr"]] <- do.call(as.matrix, c(args, pars = "^arr\\["))
    }
    if (is(x$autocor, "cor_bsts")) {
      if (is.null(newdata)) {
        draws[["loclev"]] <- do.call(as.matrix, c(args, pars = "^loclev\\["))
      } else {
        warning(paste("Local level terms are currently ignored when", 
                      "'newdata' is specified."), call. = FALSE)
      }
    }
  }
  # samples of the (non-linear) predictor
  nlpars <- names(ee$nonlinear)
  if (length(nlpars)) {
    for (i in seq_along(nlpars)) {
      nlfit <- x
      nlfit$formula <- update.formula(x$formula, 
                         rhs(attr(x$formula, "nonlinear")[[i]]))
      nlfit$ranef <- tidy_ranef(extract_effects(nlfit$formula), 
                               data = model.frame(x))
      nlstandata <- amend_newdata(newdata, fit = nlfit, re_formula = re_formula, 
                                  allow_new_levels = allow_new_levels)
      draws$nonlinear[[nlpars[i]]] <- 
        .extract_draws(nlfit, standata = nlstandata, nlpar = nlpars[i], 
                       re_formula = re_formula, subset = subset)
      rm(nlfit)
      draws$nonlinear[[nlpars[i]]][c("f", "data", "nsamples")] <- 
        list(f = draws$f, data = nlstandata, draws$nsamples)
    }
    covars <- all.vars(rhs(ee$covars))
    if (!all(covars %in% colnames(standata$C))) {
      stop("Covariate matrix is invalid. Please report a bug.")
    }
    for (i in seq_along(covars)) {
      draws[["C"]][[covars[i]]] <- 
        matrix(standata$C[, covars[i]], nrow = nsamples, 
               ncol = nrow(standata$C), byrow = TRUE)
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
        mvfit <- x
        attr(mvfit$formula, "response") <- r
        mvstandata <- amend_newdata(newdata, fit = mvfit, 
                                    re_formula = re_formula, 
                                    allow_new_levels = allow_new_levels)
        draws[["mv"]][[r]] <- 
          .extract_draws(mvfit, standata = mvstandata, subset = subset,
                         re_formula = re_formula, nlpar = r)
        rm(mvfit)
        draws[["mv"]][[r]][c("f", "data", "nsamples")] <- 
          list(draws$f, mvstandata, draws$nsamples)
      }
    } else {
      # univariate models
      draws <- c(draws, 
        .extract_draws(x, standata = standata, subset = subset,
                       re_formula = re_formula))  
    }
  }
  draws
}

.extract_draws <- function(x, standata, nlpar = "", re_formula = NULL, 
                           subset = NULL) {
  # helper function for extract_draws to extract posterior samples
  # of all regression effects
  # Args:
  #   x: a brmsift object
  #   standata: a list returned by make_standata
  #   nlpar: optional name of a non-linear parameter
  # Returns:
  #   a named list
  new_formula <- update_re_terms(formula(x), re_formula = re_formula)
  ee <- extract_effects(new_formula, family = family(x))
  nlpar_usc <- ifelse(nchar(nlpar), paste0(nlpar, "_"), "")
  args <- list(x = x, subset = subset)
  
  draws <- list(old_cat = is.old_categorical(x))
  if (!is.null(standata$X) && ncol(standata$X) && draws$old_cat != 1L) {
    b_pars <- paste0("b_", nlpar_usc, colnames(standata$X))
    draws[["b"]] <- 
      do.call(as.matrix, c(args, list(pars = b_pars, exact = TRUE)))
  }
  if (!is.null(standata$Xm) && ncol(standata$Xm)) {
    monef <- colnames(standata$Xm)
    draws[["bm"]] <- draws$simplex <- vector("list", length(monef))
    for (i in 1:length(monef)) {
      bm_par <- paste0("b_", nlpar_usc, monef[i])
      draws[["bm"]][[i]] <- 
        do.call(as.matrix, c(args, list(pars = bm_par, exact = TRUE)))
      simplex_par <- paste0("simplex_", nlpar_usc, monef[i], 
                            "[", 1:standata$Jm[i], "]")
      draws[["simplex"]][[i]] <- 
        do.call(as.matrix, c(args, list(pars = simplex_par, exact = TRUE)))
      
    }
  }
  if (is.ordinal(family(x))) {
    draws[["Intercept"]] <- 
      do.call(as.matrix, c(args, list(pars = "^b_Intercept\\[")))
    if (!is.null(standata$Xp) && ncol(standata$Xp)) {
      cse_pars <- paste0("^b_", colnames(standata$Xp), "\\[")
      draws[["cse"]] <- do.call(as.matrix, c(args, list(pars = cse_pars)))
    }
  } else if (draws$old_cat == 1L) {
    # old categorical models deprecated as of brms > 0.8.0
    if (!is.null(standata$X)) {
      draws[["cse"]] <- do.call(as.matrix, c(args, list(pars = "^b_")))
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
  
  # group-level effects
  usc_nlpar <- usc(usc(nlpar))
  new_ranef <- tidy_ranef(ee, model.frame(x))
  groups <- unique(new_ranef$group)
  # requires initialization to assign S4 objects of the Matrix package
  draws[["Z"]] <- draws[["Z_mono"]] <- draws[["Z_cse"]] <- named_list(groups)
  for (g in groups) {
    new_r <- new_ranef[new_ranef$group == g, ]
    gf <- get(paste0("J_", new_r$id[1]), standata)
    r_pars <- paste0("^r_", g, usc_nlpar, "\\[")
    r <- do.call(as.matrix, c(args, list(pars = r_pars)))
    if (is.null(r)) {
      stop("Group-level effects for each level of group '",
           g, "' not found. Please set ranef = TRUE ",
           "when calling brm.", call. = FALSE)
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
    # monotonic group-level terms separately
    new_r_mono <- new_r[new_r$type == "mono", ]
    if (nrow(new_r_mono)) {
      Z_mono <- expand_matrix(matrix(1, length(gf)), gf)
      if (length(levels) < nlevels) {
        Z_mono <- Z_mono[, levels, drop = FALSE]
      }
      draws[["Z_mono"]][[g]] <- Z_mono
      draws[["r_mono"]] <- named_list(new_r_mono$coef, list())
      for (co in names(draws[["r_mono"]])) {
        take <- which(new_r$coef == co & new_r$type == "mono")
        take <- take + nranef * (seq_along(levels) - 1)
        draws[["r_mono"]][[co]][[g]] <- r[, take, drop = FALSE]
      }
    }
    # category specific group-level terms
    new_r_cse <- new_r[new_r$type == "cse", ]
    if (nrow(new_r_cse)) {
      Z <- prepare_Z(new_r_cse, gf, standata)
      draws[["Z_cse"]][[g]] <- Z
      draws[["r_cse"]] <- named_list(seq_len(standata$ncat - 1), list())
      for (i in names(draws[["r_cse"]])) {
        index <- paste0("\\[", i, "\\]$")
        take <- which(grepl(index, new_r$coef) & new_r$type == "cse")
        take <- as.vector(outer(take, nranef * (seq_along(levels) - 1), "+"))
        draws[["r_cse"]][[i]][[g]] <- r[, take, drop = FALSE]
      }
    }
    # basic group-level effects
    new_r_basic <- new_r[!nzchar(new_r$type), ]
    if (nrow(new_r_basic)) {
      Z <- prepare_Z(new_r_basic, gf, standata)
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
