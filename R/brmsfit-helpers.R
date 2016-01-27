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
  ind <- collapse(rep(",", n.dim - 1))
  for (i in 1:dim(x)[n.dim])
    l[[i]] <- eval(parse(text = paste0("x[", ind, i,"]")))
  names(l) <- dimnames(x)[[n.dim]]
  l
}

Nsamples <- function(x) {
  # compute the number of posterior samples
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) {
    return(0)
  } 
  (x$fit@sim$iter - x$fit@sim$warmup) / x$fit@sim$thin * x$fit@sim$chains
}

algorithm <- function(x) {
  if (!is(x, "brmsfit")) {
    stop("x must be of class brmsfit")
  }
  if (is.null(x$algorithm)) "sampling"
  else x$algorithm
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
  ifelse(target <= A[, i] | ncol(A) == i, i, first_greater(A, target, i + 1))
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

get_cornames <- function(names, type = "cor", brackets = TRUE, 
                         subset = NULL, subtype = "") {
  # get correlation names as combinations of variable names
  #
  # Args:
  #  names: the variable names 
  #  type: of the correlation to be put in front of the returned strings
  #  brackets: should the correlation names contain brackets 
  #            or underscores as seperators
  #  subset: subset of correlation parameters to be returned. 
  #          Currently only used in summary.brmsfit (s3.methods.R)
  #  subtype: the subtype of the correlation (e.g., g1 in cor_g1_x_y). 
  #           Only used when subset is not NULL
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
  #   to.array: logical; should the result be transformed 
  #             into an array of increased dimension?
  #   ...: additional arguments passed to get(coef)
  #
  # Returns: 
  #   typically a matrix with colnames(samples) as colnames
  dots <- list(...)
  args <- list(X = samples, MARGIN = margin, FUN = coef)
  fun_args <- names(formals(coef))
  if (!"..." %in% fun_args) {
    dots <- dots[names(dots) %in% fun_args]
  }
  x <- do.call(apply, c(args, dots))
  if (is.null(dim(x))) {
    x <- matrix(x, dimnames = list(NULL, coef))
  } else if (coef == "quantile") {
    x <- aperm(x, length(dim(x)):1)
  }
  if (to.array && length(dim(x)) == 2) {
    x <- array(x, dim = c(dim(x), 1), 
               dimnames = list(NULL, NULL, coef))
  }
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

get_table <- function(samples, levels = sort(unique(as.numeric(samples)))) {
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
    table(factor(samples[, n], levels = levels))))
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

get_sigma <- function(x, data, n, method = c("fitted", "predict", "logLik")) {
  # get residual standard devation of linear models
  # Args:
  #   x: a brmsfit object or posterior samples of sigma (can be NULL)
  #   data: data initially passed to Stan
  #   method: S3 method from which get_sigma is called
  #   n: meaning depends on the method argument:
  #      for predict and logLik this is the current observation number
  #      for fitted this is the number of samples
  method <- match.arg(method)
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
    if (is.null(sigma)) {
      stop("No residual standard deviation(s) found")
    }
    if (method %in% c("predict", "logLik")) {
      sigma <- sigma[n]
    } else {
      sigma <- matrix(rep(sigma, n), ncol = data$N, byrow = TRUE)
    }
  }
  sigma
}

extract_pars <- function(pars, all_pars, exact_match = FALSE,
                         na_value = all_pars, ...) {
  # extract all valid parameter names that match pars
  # Args:
  #   pars: A character vector or regular expression
  #   all_pars: all parameter names of the fitted model
  #   exact_match: should parnames be matched exactly?
  #   na_value: what should be returned if pars is NA? 
  #   ...: Further arguments to be passed to grepl
  # Returns:
  #   A character vector of parameter names
  if (!(anyNA(pars) || is.character(pars))) 
    stop("Argument pars must be NA or a character vector")
  if (!anyNA(pars)) {
    if (exact_match) {
      pars <- all_pars[all_pars %in% pars]
    } else {
      pars <- all_pars[apply(sapply(pars, grepl, x = all_pars, ...), 1, any)]
    }
  } else {
    pars <- na_value
  }
  pars
}

linear_predictor <- function(x, newdata = NULL, re_formula = NULL,
                             allow_new_levels = FALSE, subset = NULL) {
  # compute the linear predictor (eta) for brms models
  #
  # Args:
  #   x: a brmsfit object
  #   newdata: optional list as returned by amend_newdata.
  #            If NULL, the standata method will be called
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
  # the linear predictor will be based on an updated formula 
  # if re_formula is specified
  new_ranef <- check_re_formula(re_formula, old_ranef = x$ranef, 
                                data = x$data)
  new_formula <- update_re_terms(x$formula, re_formula = re_formula)
  if (is.null(newdata)) { 
    data <- standata(x, re_formula = re_formula,
                     control = list(keep_intercept = TRUE))
  } else {
    data <- newdata
  }
  
  family <- family(x)
  ee <- extract_effects(new_formula, family = family)
  args <- list(x = x, as.matrix = TRUE, subset = subset)
  if (!is.null(subset)) {
    nsamples <- length(subset)
  } else {
    nsamples <- Nsamples(x)
  }
  eta <- matrix(0, nrow = nsamples, ncol = data$N)
  if (!is.null(data$X) && ncol(data$X) && !is.categorical(family)) {
    b <- do.call(posterior_samples, c(args, pars = "^b_[^\\[]+$"))
    eta <- eta + fixef_predictor(X = data$X, b = b)  
  }
  if (!is.null(data$offset)) {
    eta <- eta + matrix(rep(data$offset, nsamples), 
                        ncol = data$N, byrow = TRUE)
  }
  
  # incorporate random effects
  group <- names(new_ranef)
  for (i in seq_along(group)) {
    if (any(grepl(paste0("^J_"), names(data)))) {  # implies brms > 0.4.1
      # create a single RE design matrix for every grouping factor
      Z <- lapply(which(ee$random$group == group[i]), 
                  function(k) get(paste0("Z_",k), data))
      Z <- do.call(cbind, Z)
      id <- match(group[i], ee$random$group)
      gf <- get(paste0("J_",id), data)
    } else {  # implies brms <= 0.4.1
      Z <- as.matrix(get(paste0("Z_",group[i]), data))
      gf <- get(group[i], data)
    }
    r <- do.call(posterior_samples, 
                 c(args, pars = paste0("^r_",group[i],"\\[")))
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
      Yarr <- as.matrix(data$Yar)
      arr <- do.call(posterior_samples, c(args, pars = "^ar\\["))
    } else {
      # brms > 0.5.0
      Yarr <- as.matrix(data$Yarr)
      arr <- do.call(posterior_samples, c(args, pars = "^arr\\["))
    }
    eta <- eta + fixef_predictor(X = Yarr, b = arr)
  }
  if ((get_ar(x$autocor) || get_ma(x$autocor)) && !use_cov(x$autocor)) {
    # only run when ARMA effects were modeled as part of eta
    if (old_autocor) {
      ar <- NULL
    } else {
      ar <- do.call(posterior_samples, c(args, pars = "^ar\\["))
    }
    ma <- do.call(posterior_samples, c(args, pars = "^ma\\["))
    eta <- arma_predictor(data = data, ar = ar, ma = ma, 
                          eta = eta, link = x$link)
  }
  
  # transform eta to to etap for ordinal and categorical models
  if (is.ordinal(family)) {
    Intercept <- do.call(posterior_samples, c(args, pars = "^b_Intercept\\["))
    if (!is.null(data$Xp) && ncol(data$Xp)) {
      p <- do.call(posterior_samples, 
                   c(args, pars = paste0("^b_", colnames(data$Xp), "\\[")))
      eta <- cse_predictor(Xp = data$Xp, p = p, eta = eta, ncat = data$max_obs)
    } else {
      eta <- array(eta, dim = c(dim(eta), data$max_obs - 1))
    } 
    for (k in 1:(data$max_obs - 1)) {
      if (family$family %in% c("cumulative", "sratio")) {
        eta[, , k] <-  Intercept[, k] - eta[, , k]
      } else {
        eta[, , k] <- eta[, , k] - Intercept[, k]
      }
    }
  } else if (is.categorical(family)) {
    if (!is.null(data$Xp)) {
      p <- do.call(posterior_samples, c(args, pars = "^b_"))
      eta <- cse_predictor(Xp = data$Xp, p = p, eta = eta, ncat = data$max_obs)
    } else {
      eta <- array(eta, dim = c(dim(eta), data$max_obs - 1))
    }
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
  if (!is.matrix(X))
    stop("X must be a matrix")
  if (!is.matrix(b))
    stop("b must be a matrix")
  b %*% t(X)
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
  if (!is.matrix(Z))
    stop("Z must be a matrix")
  if (!is.matrix(r))
    stop("r must be a matrix")
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

arma_predictor <- function(data, eta, ar = NULL, ma = NULL, 
                           link = "identity") {
  # compute eta for ARMA effects
  # ToDo: use C++ for this function
  #
  # Args:
  #   data: the data initially passed to Stan
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
  Y <- link(data$Y, link)
  N <- length(Y)
  tg <- c(rep(0, K), data$tgroup)
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
  if (!is.matrix(Xp))
    stop("Xp must be a matrix")
  if (!is.matrix(p))
    stop("p must be a matrix")
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
  if (!is.matrix(A)) 
    stop("A must be a matrix")
  if (length(x) != nrow(A))
    stop("x must have nrow(A) elements")
  if (!all(is.wholenumber(x) & x > 0))
    stop("x must contain positive integers only")
  K <- ncol(A)
  i <- rep(seq_along(x), each = K)
  make_j <- function(n, K, x) K * (x[n] - 1) + 1:K
  j <- ulapply(seq_along(x), make_j, K = K, x = x)
  Matrix::sparseMatrix(i = i, j = j, x = as.vector(t(A)))
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
  IC <- do.call(eval(parse(text = paste0("loo::", ic))), args)
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
  #   A matrix with differences in the ICs 
  #   as well as corresponding standard errors
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
  nlist(ic_diffs, weights)
}

match_response <- function(models) {
  # compare the response parts of multiple brmsfit objects
  # Args:
  #  models: A list of brmsfit objects
  # Returns:
  #  TRUE if the response parts of all models match and FALSE else
  if (length(models) <= 1) return(TRUE)
  .match_fun <- function(x, y) {
    # checks if all relevant parts of the response are the same 
    # Args:
    #   x, y: named lists as returned by standata
    to_match <- c("Y", "se", "weights", "cens", "trunc")
    all(ulapply(to_match, function(v) 
      isTRUE(all.equal(as_matrix(x[[v]])[attr(x, "old_order"), ], 
                       as_matrix(y[[v]])[attr(y, "old_order"), ]))))
  } 
  standatas <- lapply(models, standata, control = list(save_order = TRUE))
  matches <- ulapply(standatas[-1], .match_fun, y = standatas[[1]]) 
  if (all(matches)) {
    out <- TRUE
  } else {
    out <- FALSE
    warning(paste("model comparisons are invalid as the response parts", 
                  "of at least two models do not match"), call. = FALSE)
  }
  out
}

find_names <- function(x) {
  # find all valid object names in a string 
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
                    "(\\[[^],]+(,[^],]+)*\\])?")
  pos_all <- gregexpr(reg_all, x)[[1]]
  reg_fun <- "([^([:digit:]|[:punct:])]|\\.)[[:alnum:]_\\.]*\\("
  pos_fun <- gregexpr(reg_fun, x)[[1]]
  pos_decnum <- gregexpr("\\.[[:digit:]]+", x)[[1]]
  pos_var <- list(rmMatch(pos_all, pos_fun, pos_decnum))
  unlist(regmatches(x, pos_var))
}

add_samples <- function(x, newpar, dim = numeric(0), dist = "norm", ...) {
  # add some random samples to a brmsfit object 
  # currently only used within tests
  # Args:
  #   x: a brmsfit object
  #   newpar: name of the new parameter to add; 
  #           a single character vector
  #   dim: dimension of the new parameter
  # Returns:
  #   a brmsfit object with new (standard normal) samples
  if (!is(x, "brmsfit")) {
    stop("x must be of class brmsfit")
  }
  if (!identical(dim, numeric(0))) {
    stop("currently dim must be numeric(0)")
  }
  for (i in seq_along(x$fit@sim$samples)) {
    x$fit@sim$samples[[i]][[newpar]] <- 
      do.call(paste0("r", dist), list(x$fit@sim$iter, ...))
  }
  x$fit@sim$fnames_oi <- c(x$fit@sim$fnames_oi, newpar) 
  x$fit@sim$dims_oi[[newpar]] <- dim
  x$fit@sim$pars_oi <- names(x$fit@sim$dims_oi)
  x
}