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
  links <- c(sigma = "exp", shape = "exp", nu = "exp", 
             phi = "exp", zi = "inv_logit", hu = "inv_logit") 
  valid_auxpars <- valid_auxpars(family(x), effects = ee, autocor = x$autocor)
  for (ap in valid_auxpars) {
    if (!is.null(ee[[ap]])) {
      auxfit <- x
      auxfit$formula <- update.formula(x$formula, rhs(attr(x$formula, ap)))
      auxfit$ranef <- gather_ranef(extract_effects(auxfit$formula), 
                                   data = model.frame(x))
      auxstandata <- amend_newdata(newdata, fit = auxfit, 
                                   re_formula = re_formula,
                                   allow_new_levels = allow_new_levels)
      draws[[ap]] <- .extract_draws(auxfit, standata = auxstandata, nlpar = ap,
                                    re_formula = re_formula, subset = subset)
      rm(auxfit)
      draws[[ap]][c("f", "data", "nsamples", "link")] <- 
        list(draws$f, auxstandata, draws$nsamples, links[ap])
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
      nlfit$ranef <- gather_ranef(extract_effects(nlfit$formula), 
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
    x$ranef <- gather_ranef(extract_effects(formula(x)), 
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
  if (!is.null(standata$X) && ncol(standata$X) && !draws$old_cat) {
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
      draws[["p"]] <- do.call(as.matrix, c(args, list(pars = cse_pars)))
    }
  } else if (draws$old_cat) {
    # old categorical models deprecated as of brms > 0.8.0
    if (!is.null(standata$X)) {
      draws[["p"]] <- do.call(as.matrix, c(args, list(pars = "^b_")))
    }
  }
  # splines
  splines <- rename(get_spline_labels(ee))
  for (i in seq_along(splines)) {
    draws[["Zs"]][[splines[i]]] <- draws$data[[paste0("Zs_", i)]]
    s_pars <- paste0("^s_", nlpar_usc, splines[i], "\\[")
    draws[["s"]][[splines[i]]] <- 
      do.call(as.matrix, c(args, list(pars = s_pars)))
  }
  # group-level effects
  usc_nlpar <- usc(nlpar, "prefix")
  new_ranef <- gather_ranef(ee, model.frame(x))
  groups <- unique(new_ranef$group)
  # requires initialization to assign S4 objects of the Matrix package
  draws[["Z"]] <- named_list(groups)
  for (g in groups) {
    new_r <- new_ranef[new_ranef$group == g, ]
    # create a single RE design matrix for every grouping factor
    Z <- lapply(unique(new_r$gn), 
                function(k) get(paste0("Z_", k), standata))
    Z <- do.call(cbind, Z)
    gf <- get(paste0("J_", new_r$id[1]), standata)
    r_pars <- paste0("^r_", g, usc_nlpar, "\\[")
    r <- do.call(as.matrix, c(args, list(pars = r_pars)))
    if (is.null(r)) {
      stop(paste0("Group-level effects for each level of group '",
                 g, "' not found. Please set ranef = TRUE ",
                 "when calling brm."), call. = FALSE)
    }
    # match columns of Z with corresponding RE estimates
    n_levels <- ngrps(x)[[g]]
    old_r <- x$ranef[x$ranef$group == g, ]
    used_re <- ulapply(new_r$coef, match, old_r$coef)
    used_re_pars <- ulapply(used_re, function(k) 
      1:n_levels + (k - 1) * n_levels)
    r <- r[, used_re_pars, drop = FALSE]
    nranef <- ncol(Z)
    max_levels <- ncol(r) / nranef
    has_new_levels <- anyNA(gf)
    if (has_new_levels) {
      # if new levels are present (only if allow_new_levels is TRUE)
      new_r <- matrix(nrow = nrow(r), ncol = nranef)
      for (k in 1:nranef) {
        # sample values of the new level for each random effect
        indices <- ((k - 1) * max_levels + 1):(k * max_levels)
        new_r[, k] <- apply(r[, indices], MARGIN = 1, FUN = sample, size = 1)
      }
      gf[is.na(gf)] <- max_levels + 1
    } else { 
      new_r <- matrix(nrow = nrow(r), ncol = 0)
    }
    # sort levels because we need row major instead of column major order
    sort_levels <- ulapply(1:max_levels, function(l) 
      seq(l, ncol(r), max_levels))
    r <- cbind(r[, sort_levels, drop = FALSE], new_r)
    Z <- expand_matrix(Z, gf)
    levels <- unique(gf)
    if (has_new_levels) max_levels <- max_levels + 1
    if (length(levels) < max_levels) {
      # if only a subset of levels is provided (only for newdata)
      take_levels <- ulapply(levels, function(l) 
        ((l - 1) * nranef + 1):(l * nranef))
      r <- r[, take_levels, drop = FALSE]
      Z <- Z[, take_levels, drop = FALSE]
    }
    draws[["Z"]][[g]] <- Z
    draws[["r"]][[g]] <- r
  }
  draws
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
