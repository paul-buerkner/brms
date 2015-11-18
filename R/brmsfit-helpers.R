array2list <- function(x) {
  # convert array to list of elements with reduced dimension
  #
  # Args: 
  #   x: an arrary of dimension d
  #
  # Returns: 
  #   A list of arrays of dimension d-1
  if (is.null(dim(x))) stop("Argument x has no dimension")
  n.dim <- length(dim(x))
  l <- list(length = dim(x)[n.dim])
  ind <- collapse(rep(",", n.dim-1))
  for (i in 1:dim(x)[n.dim])
    l[[i]] <- eval(parse(text = paste0("x[", ind, i,"]")))
  names(l) <- dimnames(x)[[n.dim]]
  l
}

first_greater <- function(A, target, i = 1) {
  # find the first element in A that is greater than target
  #
  # Args: 
  #   A: a matrix
  #   target: a vector of length nrow(A)
  #   i: column of A being checked first
  #
  # Returns: 
  #   A vector of the same length as target containing the column ids 
  #   where A[,i] was first greater than target
  ifelse(target <= A[,i] | ncol(A) == i, i, first_greater(A, target, i+1))
}

link <- function(x, link) {
  # apply a link function on x
  #
  # Args:
  #   x: An arrary of arbitrary dimension
  #   link: a character string defining the link
  #
  # Returns:
  #   an array of dimension dim(x) on which the link function was applied
  if (link == "identity") x
  else if (link == "log") log(x)
  else if (link == "inverse") 1/x
  else if (link == "sqrt") sqrt(x)
  else if (link == "1/mu^2") 1 / x^2
  else if (link == "logit") logit(x)
  else if (link == "probit") qnorm(x)
  else if (link == "probit_approx") qnorm(x)
  else if (link == "cloglog") log(-log(1-x))
  else if (link == "cauchit") qcauchy(x)
  else stop(paste("Link", link, "not supported"))
}

ilink <- function(x, link) {
  # apply the inverse link function on x
  #
  # Args:
  #   x: An arrary of arbitrary dimension
  #   link: a character string defining the link
  #
  # Returns:
  #   an array of dimension dim(x) on which the inverse link function was applied
  if (link == "identity") x
  else if (link == "log") exp(x)
  else if (link == "inverse") 1/x
  else if (link == "sqrt") x^2
  else if (link == "1/mu^2") 1 / sqrt(x)
  else if (link == "logit") ilogit(x)
  else if (link == "probit") pnorm(x)
  else if (link == "probit_approx") ilogit(0.07056*x^3 + 1.5976*x)
  else if (link == "cloglog") 1 - exp(-exp(x))
  else if (link == "cauchit") pcauchy(x)
  else stop(paste("Link", link, "not supported"))
}

get_cornames <- function(names, type = "cor", brackets = TRUE, subset = NULL, subtype = "") {
  # get correlation names as combinations of variable names
  #
  # Args:
  #  names: the variable names 
  #  type: of the correlation to be put in front of the returned strings
  #  brackets: should the correlation names contain brackets or underscores as seperators
  #  subset: subset of correlation parameters to be returned. Currently only used in summary.brmsfit (s3.methods.R)
  #  subtype: the subtype of the correlation (e.g., g1 in cor_g1_x_y). Only used when subset is not NULL
  #
  # Returns: 
  #  correlation names based on the variable names passed to the names argument
  cornames <- NULL
  if (is.null(subset) && length(names) > 1) {
    for (i in 2:length(names)) {
      for (j in 1:(i-1)) {
        if (brackets) {
          cornames <- c(cornames, paste0(type,"(",names[j],",",names[i],")"))
        } else {
          cornames <- c(cornames, paste0(type,"_",names[j],"_",names[i]))
        }
      }
    }
  } else if (!is.null(subset)) {
    possible_values <- get_cornames(names = names, type = "", brackets = FALSE)
    subset <- rename(subset, paste0("^",type, if (nchar(subtype)) paste0("_",subtype)),
                     "", fixed = FALSE)
    matches <- which(possible_values %in% subset)
    cornames <- get_cornames(names = names, type = type)[matches]
  }
  cornames
}

get_estimate <- function(coef, samples, margin = 2, to.array = FALSE, ...) {
  # calculate estimates over posterior samples 
  # 
  # Args:
  #   coef: coefficient to be applied on the samples (e.g., "mean")
  #   samples: the samples over which to apply coef
  #   margin: see apply
  #   to.array: logical; should the result be transformed into an array of increased dimension?
  #   ...: additional arguments passed to get(coef)
  #
  # Returns: 
  #   typically a matrix with colnames(samples) as colnames
  dots <- list(...)
  args <- list(X = samples, MARGIN = margin, FUN = coef)
  fun_args <- names(formals(coef))
  if (!"..." %in% fun_args)
    dots <- dots[names(dots) %in% fun_args]
  x <- do.call(apply, c(args, dots))
  if (is.null(dim(x))) 
    x <- matrix(x, dimnames = list(NULL, coef))
  else if (coef == "quantile") 
    x <- aperm(x, length(dim(x)):1)
  if (to.array && length(dim(x)) == 2) 
    x <- array(x, dim = c(dim(x), 1), dimnames = list(NULL, NULL, coef))
  x 
}

get_summary <- function(samples, probs = c(0.025, 0.975)) {
  # summarizes parameter samples based on mean, sd, and quantiles
  #
  # Args: 
  #  samples: a matrix or data.frame containing the samples to be summarized. 
  #    rows are samples, columns are parameters
  #  probs: quantiles to be computed
  #
  # Returns:
  #   a N x C matric where N is the number of observations and C is equal to \code{length(probs) + 2}.
  if (length(dim(samples)) == 2) {
    out <- do.call(cbind, lapply(c("mean", "sd", "quantile"), get_estimate, 
                                 samples = samples, probs = probs, na.rm = TRUE))
  } else if (length(dim(samples)) == 3) {
    out <- abind(lapply(1:dim(samples)[3], function(i)
      do.call(cbind, lapply(c("mean", "sd", "quantile"), get_estimate, 
                            samples = samples[, , i], probs = probs))), along = 3)
    dimnames(out) <- list(NULL, NULL, paste0("P(Y = ", 1:dim(out)[3], ")")) 
  } else { 
    stop("dimension of samples must be either 2 or 3") 
  }
  rownames(out) <- 1:nrow(out)
  colnames(out) <- c("Estimate", "Est.Error", paste0(probs * 100, "%ile"))
  out  
}

get_table <- function(samples, levels = sort(unique(samples))) {
  # compute absolute frequencies for each column
  # 
  # Args:
  #  samples: a S x N matrix
  #  levels: all possible values in \code{samples}
  # 
  # Returns:
  #   a N x \code{levels} matrix containing absolute frequencies in each column seperately
  if (!is.matrix(samples)) 
    stop("samples must be a matrix")
  out <- do.call(rbind, lapply(1:ncol(samples), function(n) 
    table(factor(samples[,n], levels = levels))))
  rownames(out) <- 1:nrow(out)
  colnames(out) <- paste0("N(Y = ", 1:ncol(out), ")")
  out
}

get_cov_matrix <- function(sd, cor = NULL) {
  # compute covariance and correlation matrices based on correlation and sd samples
  #
  # Args:
  #   sd: samples of standard deviations
  #   cor: samples of correlations
  #
  # Notes: 
  #   used in VarCorr.brmsfit
  # 
  # Returns: 
  #   samples of covariance and correlation matrices
  if (any(sd < 0)) stop("standard deviations must be non negative")
  if (!is.null(cor)) {
    if (ncol(cor) != ncol(sd) * (ncol(sd) - 1) / 2 || nrow(sd) != nrow(cor))
      stop("dimensions of standard deviations and corrrelations do not match")
    if (any(cor < -1 || cor > 1)) 
      stop("correlations must be between -1 and 1")  
  }
  nsamples <- nrow(sd)
  nranef <- ncol(sd)
  cor_matrix <- aperm(array(diag(1, nranef), dim = c(nranef, nranef, nsamples)), 
                      perm = c(3, 1, 2))
  cov_matrix <- cor_matrix
  for (i in 1:nranef) 
    cov_matrix[, i, i] <- sd[, i]^2 
  if (!is.null(cor)) {
    k <- 0 
    for (i in 2:nranef) {
      for (j in 1:(i-1)) {
        k = k + 1
        cor_matrix[, j, i] <- cor_matrix[, i, j] <- cor[, k]
        cov_matrix[, j, i] <- cov_matrix[, i, j] <- cor[, k] * sd[, i] * sd[, j]
      }
    }
  }
  list(cor = cor_matrix, cov = cov_matrix)
}

get_cov_matrix_ar1 <- function(ar, sigma, sq_se, nrows) {
  # compute the covariance matrix for an AR1 process
  # Args: 
  #   ar: AR1 autocorrelation samples
  #   sigma: standard deviation samples of the AR1 process
  #   sq_se: user defined standard errors (may be 0)
  #   nrows: number of rows of the covariance matrix
  # Returns:
  #   An nsamples x nrows x nrows AR1 covariance array (!)
  mat <- aperm(array(diag(sq_se, nrows), dim = c(nrows, nrows, nrow(ar))),
               perm = c(3, 1, 2))
  sigma2_adjusted <- sigma^2 / (1 - ar^2)
  pow_ar <- as.list(rep(1, nrows + 1))
  for (i in 1:nrows) {
    pow_ar[[i + 1]] <- ar^i
    mat[, i, i] <- mat[, i, i] + sigma2_adjusted
    if (i > 1) {
      for (j in 1:(i - 1)) { 
        mat[, i, j] <- sigma2_adjusted * pow_ar[[i - j + 1]]
        mat[, j, i] <- mat[, i, j]
      } 
    }
  } 
  mat 
}

get_cov_matrix_ma1 <- function(ma, sigma, sq_se, nrows) {
  # compute the covariance matrix for an MA1 process
  # Args: 
  #   ma: MA1 autocorrelation samples
  #   sigma: standard deviation samples of the AR1 process
  #   sq_se: user defined standard errors (may be 0)
  #   nrows: number of rows of the covariance matrix
  # Returns:
  #   An nsamples x nrows x nrows MA1 covariance array (!)
  mat <- aperm(array(diag(sq_se, nrows), dim = c(nrows, nrows, nrow(ma))),
               perm = c(3, 1, 2))
  sigma2 <- sigma^2
  sigma2_adjusted <- sigma2 * (1 + ma^2)
  sigma2_times_ma <- sigma2 * ma
  for (i in 1:nrows) { 
    mat[, i, i] <- mat[, i, i] + sigma2_adjusted
    if (i > 1) {
      mat[, i, i - 1] <- sigma2_times_ma
    }
    if (i < nrows) {
      mat[, i, i + 1] <- sigma2_times_ma
    }
  } 
  mat 
}

get_cov_matrix_arma1 <- function(ar, ma, sigma, sq_se, nrows) {
  # compute the covariance matrix for an AR1 process
  # Args: 
  #   ar: AR1 autocorrelation sample
  #   ma: MA1 autocorrelation sample
  #   sigma: standard deviation samples of the AR1 process
  #   sq_se: user defined standard errors (may be 0)
  #   nrows: number of rows of the covariance matrix
  # Returns:
  #   An nsamples x nrows x nrows ARMA1 covariance array (!)
  mat <- aperm(array(diag(sq_se, nrows), dim = c(nrows, nrows, nrow(ar))),
               perm = c(3, 1, 2))
  sigma2_adjusted <- sigma^2 / (1 - ar^2)
  gamma0 <- 1 + ma^2 + 2 * ar * ma
  gamma <- as.list(rep(NA, nrows))
  gamma[[1]] <- (1 + ar * ma) * (ar + ma)
  for (i in 1:nrows) {
    mat[, i, i] <- mat[, i, i] + sigma2_adjusted * gamma0
    gamma[[i]] <- gamma[[1]] * ar^(i - 1)
    if (i > 1) {
      for (j in 1:(i - 1)) { 
        mat[, i, j] <- sigma2_adjusted * gamma[[i - j]]
        mat[, j, i] <- mat[, i, j]
      } 
    }
  } 
  mat 
}

evidence_ratio <- function(x, cut = 0, wsign = c("equal", "less", "greater"), 
                           prior_samples = NULL, pow = 12, ...) {
  # calculate the evidence ratio between two disjunct hypotheses
  # 
  # Args:
  #   x: posterior samples 
  #   cut: the cut point between the two hypotheses
  #   wsign: direction of the hypothesis
  #   prior_samples: optional prior samples for undirected hypothesis
  #   pow influence: the accuracy of the density
  #   ...: optional arguments passed to density.default
  #
  # Returns:
  #   the evidence ratio of the two hypothesis
  wsign <- match.arg(wsign)
  if (wsign == "equal") {
    if (is.null(prior_samples)) {
      out <- NA
    } else {
      dots <- list(...)
      dots <- dots[names(dots) %in% names(formals("density.default"))]
      # compute prior and posterior densities
      prior_density <- do.call(density, c(list(x = prior_samples, n = 2^pow), dots))
      posterior_density <- do.call(density, c(list(x = x, n = 2^pow), dots))
      # evaluate densities at the cut point
      at_cut_prior <- match(min(abs(prior_density$x - cut)), 
                            abs(prior_density$x - cut))
      at_cut_posterior <- match(min(abs(posterior_density$x - cut)), 
                                abs(posterior_density$x - cut))
      out <- posterior_density$y[at_cut_posterior] / prior_density$y[at_cut_prior] 
    }
  } else if (wsign == "less") {
    out <- length(which(x < cut))
    out <- out / (length(x) - out)
  } else if (wsign == "greater") {
    out <- length(which(x > cut))
    out <- out / (length(x) - out)  
  }
  out  
}

get_sigma <- function(x, data, method, n) {
  # get residual standard devation of linear models
  # Args:
  #   x: a brmsfit object or posterior samples of sigma (can be NULL)
  #   data: data initially passed to Stan
  #   method: S3 method from which get_sigma is called
  #   n: meaning depends on the method argument:
  #      for predict and logLik this is the current observation number
  #      for fitted this is the number of samples
  if (!method %in% c("fitted", "predict", "logLik"))
    stop("Invalid method argument")
  if (is(x, "brmsfit")) {
    sigma <- posterior_samples(x, pars = "^sigma_")$sigma
  } else {
    sigma <- x
  }
  if (is.null(sigma)) {
    # user defined standard errors were applied
    sigma <- data$se
    if (is.null(sigma)) {
      # for backwards compatibility with brms <= 0.5.0
      sigma <- data$sigma
    }
    if (method %in% c("predict", "logLik")) {
      sigma <- sigma[n]
    } else {
      sigma <- matrix(rep(sigma, n), ncol = data$N, byrow = TRUE)
    }
  }
  sigma
}

linear_predictor <- function(x, newdata = NULL, re_formula = NULL) {
  # compute the linear predictor (eta) for brms models
  #
  # Args:
  #   x: a brmsfit object
  #   newdata: optional data.frame containing new data to make predictions for.
  #            If \code{NULL} (the default), the data used 
  #            to fit the model is applied.
  #   re_formula: formula containing random effects 
  #               to be considered in the prediction
  #
  # Returns:
  #   usually, an S x N matrix where S is the number of samples
  #   and N is the number of observations in the data.
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  if (is.null(newdata)) { 
    data <- standata(x, keep_intercept = TRUE)
  } else {
    data <- newdata
  }
  
  ee <- extract_effects(x$formula, family = x$family)
  n_samples <- nrow(posterior_samples(x, pars = "^lp__$"))
  eta <- matrix(0, nrow = n_samples, ncol = data$N)
  X <- data$X
  if (!is.null(X) && ncol(X) && x$family != "categorical") {
    b <- posterior_samples(x, pars = "^b_[^\\[]+$")
    eta <- eta + fixef_predictor(X = X, b = b)  
  }
  if (!is.null(data$offset)) {
    eta <- eta + matrix(rep(data$offset, n_samples), 
                        ncol = data$N, byrow = TRUE)
  }
  
  # incorporate random effects
  ee <- extract_effects(x$formula)
  group <- names(x$ranef)
  if (length(group) && is.null(re_formula)) {
    for (i in 1:length(group)) {
      if (any(grepl(paste0("^J_"), names(data)))) {  # implies brms > 0.4.1
        # create a single RE design matrix for every grouping factor
        Z <- lapply(which(ee$group == group[i]), 
                    function(k) get(paste0("Z_",k), data))
        Z <- do.call(cbind, Z)
        id <- match(group[i], ee$group)
        gf <- get(paste0("J_",id), data)
      } else {  # implies brms <= 0.4.1
        Z <- as.matrix(get(paste0("Z_",group[i]), data))
        gf <- get(group[i], data)
      }
      r <- posterior_samples(x, pars = paste0("^r_",group[i],"\\["),
                             as.matrix = TRUE)
      if (is.null(r)) {
        stop(paste("Random effects for each level of grouping factor",
                   group[i], "not found. Please set ranef = TRUE",
                   "when calling brm."))
      }
      eta <- eta + ranef_predictor(Z = Z, gf = gf, r = r) 
    }
  }
  # indicates if the model was fitted with brms <= 0.5.0
  old_autocor <- is.null(x$autocor$r)
  if (get_arr(x$autocor)) {
    # incorporate ARR effects
    if (old_autocor) {
      Yarr <- as.matrix(data$Yar)
      arr <- posterior_samples(x, pars = "^ar\\[")
    } else {
      # brms > 0.5.0
      Yarr <- as.matrix(data$Yarr)
      arr <- posterior_samples(x, pars = "^arr\\[")
    }
    eta <- eta + fixef_predictor(X = Yarr, b = arr)
  }
  if ((get_ar(x$autocor) || get_ma(x$autocor)) && !use_cov(x$autocor)) {
    # only run when ARMA effects were modeled as part of eta
    if (old_autocor) {
      ar <- NULL
    } else {
      ar <- posterior_samples(x, pars = "^ar\\[", as.matrix = TRUE)
    }
    ma <- posterior_samples(x, pars = "^ma\\[", as.matrix = TRUE)
    eta <- arma_predictor(data = data, ar = ar, ma = ma, 
                          eta = eta, link = x$link)
  }
  
  # transform eta to to etap for ordinal and categorical models
  if (is.ordinal(x$family)) {
    Intercept <- posterior_samples(x, "^b_Intercept\\[")
    if (!is.null(data$Xp) && ncol(data$Xp)) {
      p <- posterior_samples(x, paste0("^b_",colnames(data$Xp),"\\["))
      etap <- partial_predictor(Xp = data$Xp, p = p, ncat = data$max_obs)
    } else {
      etap <- array(0, dim = c(dim(eta), data$max_obs-1))
    } 
    for (k in 1:(data$max_obs-1)) {
      etap[, , k] <- etap[, , k] + eta
      if (x$family %in% c("cumulative", "sratio")) {
        etap[, , k] <-  Intercept[, k] - etap[, , k]
      } else {
        etap[, , k] <- etap[, , k] - Intercept[, k]
      }
    }
    eta <- etap
  } else if (x$family == "categorical") {
    if (!is.null(data$X)) {
      p <- posterior_samples(x, pars = "^b_")
      etap <- partial_predictor(data$X, p, data$max_obs)
    } else {
      etap <- array(0, dim = c(dim(eta), data$max_obs - 1))
    }
    for (k in 1:(data$max_obs-1)) {
      etap[, , k] <- etap[, , k] + eta
    }
    eta <- etap
  }
  eta
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
  as.matrix(b) %*% t(as.matrix(X))
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
  nranef <- ncol(Z)
  max_levels <- ncol(r) / nranef
  Z <- expand_matrix(Z, gf)
  levels <- unique(gf)
  # sort levels because we need row major instead of column major order
  sort_levels <- ulapply(1:max_levels, function(l) seq(l, ncol(r), max_levels))
  if (length(levels) < max_levels) {
    # if only a subset of levels is provided (only for newdata)
    take_levels <- ulapply(levels, function(l) ((l - 1) * nranef + 1):(l * nranef))
    eta <- as.matrix(r[, sort_levels])[, take_levels] %*% t(Z[, take_levels])
  } else {
    eta <- as.matrix(r[, sort_levels]) %*% t(Z)
  }
  eta
}

arma_predictor <- function(data, ar, ma, eta, link) {
  # compute eta for moving average effects
  #
  # Args:
  #   data: the data initially passed to stan
  #   ar: autoregressive samples (can be NULL)
  #   ma: moving average samples (can be NULL)
  #   eta: previous linear predictor samples
  #   link: the link function as character string
  #
  # Returns:
  #   new linear predictor samples updated by moving average effects
  S <- max(nrow(ar), nrow(ma))
  Kar <- ifelse(is.null(ar), 0, ncol(ar))
  Kma <- ifelse(is.null(ma), 0, ncol(ma))
  K <- max(Kar, Kma)
  Ks <- 1:K
  Y <- link(data$Y, link)
  N <- length(Y)
  tg <- c(rep(0, K), data$tgroup)
  E <- array(0, dim = c(S, K, N))
  e <- matrix(0, nrow = S, ncol = N)
  for (n in 1:N) {
    if (Kma) {
      # add MA effects
      eta[, n] <- eta[, n] + rowSums(ma * E[, 1:Kma, n])
    }
    e[, n] <- Y[n] - eta[, n]
    if (n < N) {
      I <- which(n < N & tg[n + 1 + K] == tg[n + 1 + K - Ks])
      E[, I, n + 1] <- e[, n + 1 - I]
    }
    if (Kar) {
      # add AR effects
      eta[, n] <- eta[, n] + rowSums(ar * E[, 1:Kar, n])
    }
  }
  eta
}

partial_predictor <- function(Xp, p, ncat) {
  # compute etap for partial and categorical effects
  # 
  # Args:
  #   Xp: partial design matrix 
  #   p: partial effects samples
  #   ncat: number of categories
  #
  # @return linear predictor of partial effects as a 3D array (not as a matrix)
  ncat <- max(ncat)
  etap <- array(0, dim = c(nrow(p), nrow(Xp), ncat - 1))
  indices <- seq(1, (ncat - 1) * ncol(Xp), ncat - 1) - 1
  for (k in 1:(ncat-1)) {
    etap[, , k] <- as.matrix(p[, indices + k]) %*% t(as.matrix(Xp))
  }
  etap
}

expand_matrix <- function(A, x) {
  # expand a matrix into a sparse matrix of higher dimension
  # 
  # Args:
  #   A: matrix to be expanded
  #   x: levels to expand the matrix
  # 
  # Notes:
  #   used in linear_predictor.brmsfit
  #
  # Returns:
  #   An expanded matrix of dimensions nrow(A) and ncol(A) * length(unique(x)) 
  A <- as.matrix(A)
  if (length(x) != nrow(A))
    stop("x must have nrow(A) elements")
  if (!all(is.wholenumber(x) & x > 0))
    stop("x must contain positive integers only")
  K <- ncol(A)
  v <- rep(0, K * max(x))
  do.call(rbind, lapply(1:nrow(A), function(n, v) {
    v[K * (x[n] - 1) + 1:K] <- A[n, ] 
    return(v)}, v = v))
}

compute_ic <- function(x, ic = c("waic", "loo"), ...) {
  # compute WAIC and LOO using the 'loo' package
  #
  # Args:
  #   x: an object of class brmsfit
  #   ic: the information criterion to be computed
  #
  # Returns:
  #   output of the loo package with amended class attribute
  ic <- match.arg(ic)
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples") 
  args <- list(x = logLik(x))
  if (ic == "loo") args <- c(args, ...)
  IC <- do.call(eval(parse(text = paste0("loo::",ic))), args)
  class(IC) <- c("ic", "loo")
  return(IC)
}

compare_ic <- function(x, ic = c("waic", "loo")) {
  # compare information criteria of different models
  #
  # Args:
  #   x: A list containing loo objects
  #   ic: the information criterion to be computed
  #
  # Returns:
  #   A matrix with differences in the ICs as well as corresponding standard errors
  ic <- match.arg(ic)
  n_models <- length(x)
  ic_diffs <- matrix(0, nrow = n_models * (n_models - 1) / 2, ncol = 2)
  rnames <- rep("", nrow(ic_diffs))
  # pairwise comparision to get differences in ICs and their SEs
  n <- 1
  for (i in 1:(n_models - 1)) {
    for (j in (i + 1):n_models) {
      temp <- loo::compare(x[[j]], x[[i]])
      ic_diffs[n, ] <- c(-2 * temp[["elpd_diff"]], 2 * temp[["se"]]) 
      rnames[n] <- paste(names(x)[i], "-", names(x)[j])
      n <- n + 1
    }
  }
  rownames(ic_diffs) <- rnames
  colnames(ic_diffs) <- c(toupper(ic), "SE")
  # compare all models at once to obtain weights
  all_compare <- do.call(loo::compare, x)
  if (n_models == 2) {
    # weights are named differently when comparing only 2 models
    weights <- unname(all_compare[c("weight1", "weight2")])
  } else {
    # weights must be resorted as loo::compare sorts models after weights
    get_input_names <- function(...) {
      # mimic the way loo::compare defines model names
      as.character(match.call())[-1L]
    }
    weights <- unname(all_compare[do.call(get_input_names, x), "weights"])
  }
  list(ic_diffs = ic_diffs, weights = weights)
}

find_names <- function(x) {
  # find all valid object names in a string (used in method hypothesis in s3.methods.R)
  # 
  # Args:
  #   x: a character string
  #
  # Notes:
  #   currently used in hypothesis.brmsfit
  #
  # Returns:
  #   all valid variable names within the string
  if (!is.character(x) || length(x) > 1) 
    stop("x must be a character string of length 1")
  x <- gsub(" ", "", x)
  reg_all <- paste0("([^([:digit:]|[:punct:])]|\\.)[[:alnum:]_\\.]*", 
                    "(\\[[[:digit:]]+(,[[:digit:]]+)*\\])?")
  pos_all <- gregexpr(reg_all, x)[[1]]
  reg_fun <- "([^([:digit:]|[:punct:])]|\\.)[[:alnum:]_\\.]*\\("
  pos_fun <- gregexpr(reg_fun, x)[[1]]
  pos_decnum <- gregexpr("\\.[[:digit:]]+", x)[[1]]
  pos_var <- list(rmMatch(pos_all, pos_fun, pos_decnum))
  unlist(regmatches(x, pos_var))
}

td_plot <- function(par, x) {
  # trace and density plots for one parameter
  #
  # Args:
  #   par: a single character string 
  #   x: a data.frame containing the samples
  #
  # Returns:
  #   a list containing the trace and density plot of one parameter
  if (!is.character(par) || length(par) != 1)
    stop("par must be a character string")
  if (!is.data.frame(x))
    stop("x must be a data.frame")
  names(x)[match(par, names(x))] <- "value" 
  trace <- ggplot(x, aes_string(x = "iter", y = "value", group = "chains", 
                                colour = "chains")) +
    geom_line(alpha = 0.7) + 
    xlab("") + ylab("") + ggtitle(paste("Trace of", par)) + 
    theme(legend.position = "none",
          plot.title = element_text(size = 15, vjust = 1),
          plot.margin = grid::unit(c(0.2, 0, -0.5, -0.5), "lines"))
  density <- ggplot(x, aes_string(x = "value")) + 
    geom_density(aes_string(fill = "chains"), alpha = 0.5) + 
    xlab("") + ylab("") + ggtitle(paste("Density of", par)) + 
    theme(plot.title = element_text(size = 15, vjust = 1),
          plot.margin = grid::unit(c(0.2, 0, -0.5, -0.5), "lines"))
  list(trace, density)
}

#' @export
print.brmssummary <- function(x, digits = 2, ...) {
  cat(paste0(" Family: ", x$family, " (", x$link, ") \n"))
  cat(paste("Formula:", gsub(" {1,}", " ", Reduce(paste, deparse(x$formula))), "\n"))
  cat(paste0("   Data: ", x$data.name, " (Number of observations: ",x$nobs,") \n"))
  if (x$sampler == "") {
    cat(paste("\nThe model does not contain posterior samples."))
  }
  else {
    final_samples <- (x$n.iter - x$n.warmup) / x$n.thin * x$n.chains
    waic <- ifelse(is.numeric(x$WAIC), round(x$WAIC, digits = digits), x$WAIC)
    cat(paste0("Samples: ", x$n.chains, " chains, each with n.iter = ", x$n.iter, 
               "; n.warmup = ", x$n.warmup, "; n.thin = ", x$n.thin, "; \n",
               "         total post-warmup samples = ", final_samples, "\n"))
    cat(paste0("   WAIC: ", waic, "\n \n"))
    
    if (length(x$group)) {
      cat("Random Effects: \n")
      for (i in 1:length(x$group)) {
        g <- x$group[i]
        cat(paste0("~",g," (Number of levels: ",x$ngrps[[g]],") \n"))
        x$random[[g]][, "Eff.Sample"] <- 
          round(x$random[[g]][, "Eff.Sample"], digits = 0)
        print(round(x$random[[g]], digits = digits))
        cat("\n")
      }
    }
    
    if (nrow(x$cor_pars)) {
      cat("Correlation Structure: ")
      print(x$autocor)
      cat("\n")
      x$cor_pars[, "Eff.Sample"] <- round(x$cor_pars[, "Eff.Sample"], digits = 0)
      print(round(x$cor_pars, digits = digits))
      cat("\n")
    }
    
    cat("Fixed Effects: \n")
    x$fixed[, "Eff.Sample"] <- round(x$fixed[, "Eff.Sample"], digits = 0)
    print(round(x$fixed, digits = digits)) 
    cat("\n")
    
    if (nrow(x$spec_pars)) {
      cat("Family Specific Parameters: \n")
      x$spec_pars[, "Eff.Sample"] <- 
        round(x$spec_pars[, "Eff.Sample"], digits = 0)
      print(round(x$spec_pars, digits = digits))
      cat("\n")
    }
    
    cat(paste0("Samples were drawn using ",x$sampler,". For each parameter, Eff.Sample is a \n",
               "crude measure of effective sample size, and Rhat is the potential scale \n",
               "reduction factor on split chains (at convergence, Rhat = 1)."))
  }
}

#' @export
as.data.frame.VarCorr_brmsfit <- function(x, ...) {
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
        # covariances and correlations present
        # add correlations
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
print.VarCorr_brmsfit <- function(x, digits = 2, ...) {
  dat <- as.data.frame(x)
  dat[, 4:ncol(dat)] <- round(as.matrix(dat[, 4:ncol(dat)]), digits = digits)
  dat[is.na(dat)] <- ""
  print(dat, row.names = FALSE, ...)
}

#' @export
print.brmshypothesis <- function(x, digits = 2, ...) {
  cat(paste0("Hypothesis Tests for class ", x$class, ":\n"))
  x$hypothesis[, 1:5] <- round(x$hypothesis[, 1:5], digits = digits)
  print(x$hypothesis, quote = FALSE)
  cat(paste0("---\n'*': The expected value under the hypothesis lies outside the ",
             (1 - x$alpha) * 100, "% CI."))
}

#' @export
print.brmsmodel <- function(x, ...) cat(x)

#' @export
print.ic <- function(x, digits = 2, ...) {
  # print the output of LOO(x) and WAIC(x)
  ic <- names(x)[3]
  mat <- matrix(c(x[[ic]], x[[paste0("se_",ic)]]), ncol = 2, 
                dimnames = list("", c(toupper(ic), "SE")))
  print(round(mat, digits = digits))
}

#' @export
print.iclist <- function(x, digits = 2, ...) {
  # print the output of LOO(x1, x2, ...) and WAIC(x1, x2, ...)
  ic <- names(x[[1]])[3]
  mat <- matrix(0, nrow = length(x), ncol = 2, 
                dimnames = list(names(x), c(toupper(ic), "SE")))
  for (i in 1:length(x)) { 
    mat[i, ] <- c(x[[i]][[ic]], x[[i]][[paste0("se_",ic)]])
  }
  if (is.matrix(attr(x, "compare"))) {
    # models were compared using the compare_ic function
    mat <- rbind(mat, attr(x, "compare"))
    weights <- c(attr(x, "weights"), rep(NA, nrow(attr(x, "compare")))) 
    mat <- cbind(mat, Weights = weights)
  }
  print(round(mat, digits = digits), na.print = "")
}