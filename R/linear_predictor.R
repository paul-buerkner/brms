linear_predictor <- function(x, standata, re_formula = NULL,
                             subset = NULL, nlpar = "") {
  # compute the linear predictor (eta) for brms models
  #
  # Args:
  #   x: a brmsfit object
  #   standata: a list containing the output of make_standata
  #             which was possibly called with newdata
  #   re_formula: formula containing random effects 
  #               to be considered in the prediction
  #   subset: A numeric vector indicating the posterior samples to be used.
  #           If NULL, all samples are used.
  #   nlpar: optional non-linear parameter for which to compute
  #          the linear predictor
  #
  # Returns:
  #   usually, an S x N matrix where S is the number of samples
  #   and N is the number of observations in the data.
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  # the linear predictor will be based on an updated formula 
  # if re_formula is specified
  new_ranef <- check_re_formula(re_formula, old_ranef = x$ranef, 
                                data = x$data)
  new_formula <- update_re_terms(x$formula, re_formula = re_formula)
  family <- family(x)
  # do not use the nonlinear argument here
  ee <- extract_effects(new_formula, family = family)
  args <- list(x = x, as.matrix = TRUE, subset = subset)
  nsamples <- if (!is.null(subset)) length(subset) else Nsamples(x)
  nlpar <- if (nchar(nlpar)) paste0(nlpar, "_")
  
  old_cat <- is.old_categorical(x)
  eta <- matrix(0, nrow = nsamples, ncol = standata$N)
  if (!is.null(standata$X) && ncol(standata$X) && !old_cat) {
    b_pars <- paste0("b_", nlpar, colnames(standata$X))
    b <- do.call(posterior_samples, 
                 c(args, list(pars = b_pars, exact = TRUE)))
    eta <- eta + fixef_predictor(X = standata$X, b = b)  
  }
  if (!is.null(standata$offset)) {
    eta <- eta + matrix(rep(standata$offset, nsamples), 
                        ncol = standata$N, byrow = TRUE)
  }
  # incorporate monotonous effects
  if (!is.null(standata$Xm) && ncol(standata$Xm)) {
    monef <- colnames(standata$Xm)
    for (i in 1:ncol(standata$Xm)) {
      bm_par <- paste0("b_", monef[i])
      bm <- do.call(posterior_samples, 
                    c(args, list(pars = bm_par, exact = TRUE)))
      simplex_par <- paste0("simplex_", monef[i], "[", 1:standata$Jm[i], "]")
      simplex <- do.call(posterior_samples, 
                         c(args, list(pars = simplex_par, exact = TRUE)))
      eta <- eta + monef_predictor(Xm = standata$Xm[, i], bm = as.vector(bm), 
                                   simplex = simplex)
    }
  }
  
  # incorporate random effects
  group <- names(new_ranef)
  for (i in seq_along(group)) {
    if (any(grepl("^J_", names(standata)))) {  # implies brms > 0.4.1
      # create a single RE design matrix for every grouping factor
      Z <- lapply(which(ee$random$group == group[i]), 
                  function(k) get(paste0("Z_", k), standata))
      Z <- do.call(cbind, Z)
      id <- match(group[i], ee$random$group)
      gf <- get(paste0("J_", id), standata)
    } else {  # implies brms <= 0.4.1
      Z <- as.matrix(get(paste0("Z_", group[i]), standata))
      gf <- get(group[i], standata)
    }
    r_pars <- paste0("^r_", nlpar, group[i],"\\[")
    r <- do.call(posterior_samples, c(args, list(pars = r_pars)))
    if (is.null(r)) {
      stop(paste("Random effects for each level of grouping factor",
                 group[i], "not found. Please set ranef = TRUE",
                 "when calling brm."), call. = FALSE)
    }
    # match columns of Z with corresponding RE estimates
    n_levels <- ngrps(x)[[group[[i]]]]
    used_re <- ulapply(new_ranef[[group[i]]], match, x$ranef[[group[i]]])
    used_re_pars <- ulapply(used_re, function(j) 
      1:n_levels + (j - 1) * n_levels)
    r <- r[, used_re_pars, drop = FALSE]
    # add REs to linear predictor
    eta <- eta + ranef_predictor(Z = Z, gf = gf, r = r) 
    rm(r)
  }
  # indicates if the model was fitted with brms <= 0.5.0
  old_autocor <- is.null(x$autocor$r)
  if (get_arr(x$autocor)) {
    # incorporate ARR effects
    if (old_autocor) {
      Yarr <- as.matrix(standata$Yar)
      arr <- do.call(posterior_samples, c(args, list(pars = "^ar\\[")))
    } else {
      # brms > 0.5.0
      Yarr <- as.matrix(standata$Yarr)
      arr <- do.call(posterior_samples, c(args, list(pars = "^arr\\[")))
    }
    eta <- eta + fixef_predictor(X = Yarr, b = arr)
  }
  if ((get_ar(x$autocor) || get_ma(x$autocor)) && !use_cov(x$autocor)) {
    # only run when ARMA effects were modeled as part of eta
    if (old_autocor) {
      ar <- NULL
    } else {
      ar <- do.call(posterior_samples, c(args, list(pars = "^ar\\[")))
    }
    ma <- do.call(posterior_samples, c(args, list(pars = "^ma\\[")))
    eta <- arma_predictor(standata = standata, ar = ar, ma = ma, 
                          eta = eta, link = x$link)
  }
  
  # transform eta to to etap for ordinal and categorical models
  if (is.ordinal(family)) {
    Intercept <- do.call(posterior_samples, 
                         c(args, list(pars = "^b_Intercept\\[")))
    if (!is.null(standata$Xp) && ncol(standata$Xp)) {
      cse_pars <- paste0("^b_", colnames(standata$Xp), "\\[")
      p <- do.call(posterior_samples, c(args, list(pars = cse_pars)))
      eta <- cse_predictor(Xp = standata$Xp, p = p, eta = eta, 
                           ncat = standata$max_obs)
    } else {
      eta <- array(eta, dim = c(dim(eta), standata$max_obs - 1))
    } 
    for (k in 1:(standata$max_obs - 1)) {
      if (family$family %in% c("cumulative", "sratio")) {
        eta[, , k] <-  Intercept[, k] - eta[, , k]
      } else {
        eta[, , k] <- eta[, , k] - Intercept[, k]
      }
    }
  } else if (is.categorical(family)) {
    if (old_cat) {
      # deprecated as of brms > 0.8.0
      if (!is.null(standata$X)) {
        p <- do.call(posterior_samples, c(args, list(pars = "^b_")))
        eta <- cse_predictor(Xp = standata$X, p = p, eta = eta, 
                             ncat = standata$max_obs)
      } else {
        eta <- array(eta, dim = c(dim(eta), standata$max_obs - 1))
      }
    } else {
      ncat1 <- standata$ncat - 1 
      eta <- array(eta, dim = c(nrow(eta), ncol(eta) / ncat1, ncat1))
    }
  }
  eta
}

nonlinear_predictor <- function(x, newdata = NULL, C = NULL, 
                                re_formula = NULL, subset = NULL,
                                allow_new_levels = FALSE) {
  # compute the non-linear predictor (eta) for brms models
  #
  # Args:
  #   x: a brmsfit object
  #   newdata: optional list as returned by amend_newdata.
  #            If NULL, the standata method will be called
  #   C: matrix containing values of the covariates
  #   re_formula: formula containing random effects 
  #               to be considered in the prediction
  #   subset: A numeric vector indicating the posterior samples to be used.
  #           If NULL, all samples are used.
  #
  # Returns:
  #   usually, an S x N matrix where S is the number of samples
  #   and N is the number of observations in the data.
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  # the nonlinear predictor will be based on an updated formula 
  # if re_formula is specified
  new_ranef <- check_re_formula(re_formula, old_ranef = x$ranef, 
                                data = x$data)
  new_formula <- update_re_terms(x$formula, re_formula = re_formula)
  family <- family(x)
  ee <- extract_effects(new_formula, family = family, nonlinear = x$nonlinear)
  args <- list(x = x, as.matrix = TRUE, subset = subset)
  nsamples <- if (!is.null(subset)) length(subset) else Nsamples(x)
  
  # prepare non-linear model of eta 
  nlmodel_list <- list()
  nlpars <- names(ee$nonlinear)
  for (i in seq_along(nlpars)) {
    nlfit <- x
    # call linear_predictor for each non-linear parameter
    nlfit$formula <- update(x$formula, rhs(x$nonlinear[[i]]))
    nlfit$nonlinear <- NULL
    nlfit$ranef <- gather_ranef(extract_effects(nlfit$formula), data = x$data)
    nlstandata <- amend_newdata(newdata, fit = nlfit, re_formula = re_formula, 
                                allow_new_levels = allow_new_levels)
    nlmodel_list[[nlpars[i]]] <- 
      linear_predictor(x = nlfit, standata = nlstandata, nlpar = nlpars[i],
                       re_formula = re_formula, subset = subset)
  }
  covars <- all.vars(rhs(ee$covars))
  if (!all(covars %in% colnames(C))) 
    stop("Covariate matrix is invalid. Please report a bug.")
  for (i in seq_along(covars)) {
    nlmodel_list[[covars[i]]] <- matrix(C[, covars[i]], nrow = nsamples, 
                                        ncol = nrow(C), byrow = TRUE)
  }
  # evaluate non-linear predictor
  out <- try(with(nlmodel_list, eval(ee$fixed[[3]])), silent = TRUE)
  if (is(out, "try-error")) {
    if (grepl("could not find function", out)) {
      out <- rename(out, "Error in eval(expr, envir, enclos) : ", "")
      stop(paste0(out, "Most likely this is because you used a Stan ",
                  "function in the non-linear model formula that ",
                  "is not defined in R. Currently, you have to write ",
                  "this function yourself making sure that it is ",
                  "vectorized. I apologize for the inconvenience."),
           call. = FALSE)
    } else {
      out <- rename(out, "^Error :", "", fixed = FALSE)
      stop(out, call. = FALSE)
    }
  }
  out
}

fixef_predictor <- function(X, b) {
  # compute eta for fixed effects
  #
  # Args:
  #   X: fixed effects design matrix
  #   b: fixed effects samples
  # 
  # Returns:
  #   linear predictor for fixed effects
  stopifnot(is.matrix(X))
  stopifnot(is.matrix(b))
  b %*% t(X)
}

monef_predictor <- function(Xm, bm, simplex) {
  stopifnot(is.vector(Xm))
  stopifnot(is.vector(bm))
  stopifnot(is.matrix(simplex))
  bm <- as.vector(bm)
  for (i in 2:ncol(simplex)) {
    # compute the cumulative representation of the simplex 
    simplex[, i] <- simplex[, i] + simplex[, i - 1]
  }
  simplex <- cbind(0, simplex)
  bm * simplex[, Xm + 1]
}

ranef_predictor <- function(Z, gf, r) {
  # compute eta for random effects
  #  
  # Args:
  #   Z: random effects design matrix
  #   gf: levels of grouping factor for each observation
  #   r: random effects samples
  #
  # Returns: 
  #   linear predictor for random effects
  stopifnot(is.matrix(Z))
  stopifnot(is.matrix(r))
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
  if (has_new_levels) max_levels <- max_levels + 1
  # compute RE part of eta
  Z <- expand_matrix(Z, gf)
  levels <- unique(gf)
  if (length(levels) < max_levels) {
    # if only a subset of levels is provided (only for newdata)
    take_levels <- ulapply(levels, function(l) 
      ((l - 1) * nranef + 1):(l * nranef))
    eta <- r[, take_levels, drop = FALSE] %*% 
      Matrix::t(Z[, take_levels, drop = FALSE])
  } else {
    eta <- r %*% Matrix::t(Z)
  }
  # Matrix should currently not be used outside of this function
  Matrix::as.matrix(eta)
}

arma_predictor <- function(standata, eta, ar = NULL, ma = NULL, 
                           link = "identity") {
  # compute eta for ARMA effects
  # ToDo: use C++ for this function
  #
  # Args:
  #   standata: the data initially passed to Stan
  #   eta: previous linear predictor samples
  #   ar: autoregressive samples (can be NULL)
  #   ma: moving average samples (can be NULL)
  #   link: the link function as character string
  #
  # Returns:
  #   new linear predictor samples updated by ARMA effects
  S <- nrow(eta)
  Kar <- ifelse(is.null(ar), 0, ncol(ar))
  Kma <- ifelse(is.null(ma), 0, ncol(ma))
  K <- max(Kar, Kma, 1)
  Ks <- 1:K
  Y <- link(standata$Y, link)
  N <- length(Y)
  tg <- c(rep(0, K), standata$tg)
  E <- array(0, dim = c(S, K, K + 1))
  e <- matrix(0, nrow = S, ncol = K)
  zero_mat <- e
  zero_vec <- rep(0, S)
  for (n in 1:N) {
    if (Kma) {
      # add MA effects
      eta[, n] <- eta[, n] + rowSums(ma * E[, 1:Kma, K])
    }
    e[, K] <- Y[n] - eta[, n]
    if (n < N) {
      I <- which(n < N & tg[n + 1 + K] == tg[n + 1 + K - Ks])
      E[, I, K + 1] <- e[, K + 1 - I]
    }
    if (Kar) {
      # add AR effects
      eta[, n] <- eta[, n] + rowSums(ar * E[, 1:Kar, K])
    }
    # allows to keep the object size of e and E small
    E <- abind(E[, , 2:(K + 1), drop = FALSE], zero_mat)
    if (K > 1) {
      e <- cbind(e[, 2:K, drop = FALSE], zero_vec)
    }
  }
  eta
}

cse_predictor <- function(Xp, p, eta, ncat) {
  # add category specific effects to eta
  # 
  # Args:
  #   Xp: category specific design matrix 
  #   p: category specific effects samples
  #   ncat: number of categories
  #   eta: linear predictor matrix
  #
  # Returns: 
  #   linear predictor including category specific effects as a 3D array
  stopifnot(is.matrix(Xp))
  stopifnot(is.matrix(p))
  ncat <- max(ncat)
  eta <- array(eta, dim = c(dim(eta), ncat - 1))
  indices <- seq(1, (ncat - 1) * ncol(Xp), ncat - 1) - 1
  Xp <- t(Xp)
  for (k in 1:(ncat - 1)) {
    eta[, , k] <- eta[, , k] + p[, indices + k, drop = FALSE] %*% Xp
  }
  eta
}

expand_matrix <- function(A, x) {
  # expand a matrix into a sparse matrix of higher dimension
  # 
  # Args:
  #   A: matrix to be expanded
  #   x: levels to expand the matrix
  # 
  # Notes:
  #   used in linear_predictor
  #
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