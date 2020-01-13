contains_samples <- function(x) {
  if (!(is.brmsfit(x) && length(x$fit@sim))) {
    stop2("The model does not contain posterior samples.")
  }
  invisible(TRUE)
}

is_mv <- function(x) {
  stopifnot(is.brmsfit(x))
  is.mvbrmsformula(x$formula)
}

stopifnot_resp <- function(x, resp = NULL) {
  if (is_mv(x) && length(resp) != 1L) {
    stop2("Argument 'resp' must be a single variable name ", 
          "when applying this method to a multivariate model.")
  }
  invisible(NULL)
}

# apply a link function 
# @param x an arrary of arbitrary dimension
# @param link character string defining the link
link <- function(x, link) {
  switch(link, 
    "identity" = x, 
    "log" = log(x), 
    "logm1" = logm1(x),
    "log1p" = log1p(x),
    "inverse" = 1 / x,
    "sqrt" = sqrt(x), 
    "1/mu^2" = 1 / x^2, 
    "tan_half" = tan(x / 2),
    "logit" = logit(x), 
    "probit" = qnorm(x), 
    "cauchit" = qcauchy(x),
    "cloglog" = cloglog(x), 
    "probit_approx" = qnorm(x),
    "softplus" = log_expm1(x),
    stop2("Link '", link, "' not supported.")
  )
}

# apply an inverse link function
# @param x an arrary of arbitrary dimension
# @param link a character string defining the link
ilink <- function(x, link) {
  switch(link, 
    "identity" = x, 
    "log" = exp(x),
    "logm1" = expp1(x),
    "log1p" = expm1(x),
    "inverse" = 1 / x,
    "sqrt" = x^2, 
    "1/mu^2" = 1 / sqrt(x), 
    "tan_half" = 2 * atan(x),
    "logit" = inv_logit(x), 
    "probit" = pnorm(x), 
    "cauchit" = pcauchy(x),
    "cloglog" = inv_cloglog(x), 
    "probit_approx" = pnorm(x),
    "softplus" = log1p_exp(x),
    stop2("Link '", link, "' not supported.")
  )
}

# generate integers indicating subsets of the posterior samples
subset_samples <- function(x, subset = NULL, nsamples = NULL) {
  stopifnot(is.brmsfit(x))
  if (is.null(subset) && !is.null(nsamples)) {
    subset <- sample(nsamples(x), nsamples)
  }
  subset
}

# get correlation names as combinations of variable names
# @param names the variable names 
# @param type character string to be put in front of the returned strings
# @param brackets should the correlation names contain brackets 
#   or underscores as seperators?
# @param sep character string to separate names; only used if !brackets
# @return a vector of character strings
get_cornames <- function(names, type = "cor", brackets = TRUE, sep = "__") {
  cornames <- NULL
  if (length(names) > 1) {
    for (i in seq_along(names)[-1]) {
      for (j in seq_len(i - 1)) {
        if (brackets) {
          c(cornames) <- paste0(type, "(", names[j], "," , names[i], ")")
        } else {
          c(cornames) <- paste0(type, sep, names[j], sep, names[i])
        }
      }
    }
  }
  cornames
}

# extract names of categorical variables in the model
get_cat_vars <- function(x) {
  stopifnot(is.brmsfit(x))
  like_factor <- sapply(model.frame(x), is_like_factor)
  valid_groups <- c(
    names(model.frame(x))[like_factor],
    get_group_vars(x)
  )
  unique(valid_groups[nzchar(valid_groups)])
}

# covariance matrices based on correlation and SD samples
# @param sd matrix of samples of standard deviations
# @param cor matrix of samples of correlations
get_cov_matrix <- function(sd, cor = NULL) {
  sd <- as.matrix(sd)
  stopifnot(all(sd >= 0))
  nsamples <- nrow(sd)
  size <- ncol(sd)
  out <- array(diag(1, size), dim = c(size, size, nsamples))
  out <- aperm(out, perm = c(3, 1, 2))
  for (i in seq_len(size)) { 
    out[, i, i] <- sd[, i]^2
  }
  if (length(cor)) {
    cor <- as.matrix(cor)
    stopifnot(nrow(sd) == nrow(cor))
    stopifnot(min(cor) >= -1, max(cor) <= 1)
    stopifnot(ncol(cor) == size * (size - 1) / 2)
    k <- 0 
    for (i in seq_len(size)[-1]) {
      for (j in seq_len(i - 1)) {
        k = k + 1
        out[, j, i] <- out[, i, j] <- cor[, k] * sd[, i] * sd[, j]
      }
    }
  }
  out
}

# correlation matrices based on correlation samples
# @param cor samples of correlations
# @param size optional size of the desired correlation matrix;
#   ignored is 'cor' is specified
# @param nsamples optional number of posterior samples;
#   ignored is 'cor' is specified
get_cor_matrix <- function(cor, size = NULL, nsamples = NULL) {
  if (length(cor)) {
    cor <- as.matrix(cor)
    size <- -1 / 2 + sqrt(1 / 4 + 2 * ncol(cor)) + 1
    nsamples <- nrow(cor)
  } 
  size <- as_one_numeric(size)
  nsamples <- as_one_numeric(nsamples)
  stopifnot(is_wholenumber(size) && size > 0)
  stopifnot(is_wholenumber(nsamples) && nsamples > 0)
  out <- array(diag(1, size), dim = c(size, size, nsamples))
  out <- aperm(out, perm = c(3, 1, 2))
  if (length(cor)) {
    k <- 0 
    for (i in seq_len(size)[-1]) {
      for (j in seq_len(i - 1)) {
        k = k + 1
        out[, j, i] <- out[, i, j] <- cor[, k]
      }
    }
  }
  out
}

# compute covariance matrices of autocor structures
# @param draws a brmsdraws object
# @param obs observations for which to compute the covariance matrix
# @param latent compute covariance matrix for latent residuals?
get_cov_matrix_autocor <- function(draws, obs, latent = FALSE) {
  nobs <- length(obs)
  nsamples <- draws$nsamples
  acef <- draws$ac$acef
  # prepare correlations
  if (has_ac_class(acef, "arma")) {
    ar <- as.numeric(draws$ac$ar)
    ma <- as.numeric(draws$ac$ma)
    if (length(ar) && !length(ma)) {
      cor <- get_cor_matrix_ar1(ar, nobs)
    } else if (!length(ar) && length(ma)) {
      cor <- get_cor_matrix_ma1(ma, nobs)
    } else if (length(ar) && length(ma)) {
      cor <- get_cor_matrix_arma1(ar, ma, nobs)
    } else {
      stop2("Neither 'ar' nor 'ma' were supplied. Please report a bug.")
    }
  } else if (has_ac_class(acef, "cosy")) {
    cosy <- as.numeric(draws$ac$cosy)
    cor <- get_cor_matrix_cosy(cosy, nobs)
  } else {
    cor <- get_cor_matrix_ident(nsamples, nobs)
  }
  # prepare known standard errors
  if (!is.null(draws$data$se)) {
    se2 <- draws$data$se[obs]^2
    se2 <- array(diag(se2, nobs), dim = c(nobs, nobs, nsamples))
    se2 <- aperm(se2, perm = c(3, 1, 2))
    # make sure not to add 'se' twice
    draws$data$se <- NULL
  } else {
    se2 <- rep(0, nobs)
  }
  # prepare residual standard deviations
  if (latent) {
    sigma2 <- as.numeric(draws$ac$sderr)^2
  } else {
    sigma <- get_dpar(draws, "sigma", i = obs)
    if (NCOL(sigma) > 1L) {
      # sigma varies across observations
      sigma2 <- array(dim = c(nsamples, nobs, nobs))
      for (s in seq_rows(sigma2)) {
        sigma2[s, , ] <- outer(sigma[s, ], sigma[s, ])
      }
    } else {
      sigma2 <- as.numeric(sigma)^2
    }
  }
  sigma2 * cor + se2
}

# compute AR1 correlation matrices
# @param ar AR1 autocorrelation samples
# @param nobs number of rows of the covariance matrix
# @return a numeric 'nsamples' x 'nobs' x 'nobs' array
get_cor_matrix_ar1 <- function(ar, nobs) {
  out <- array(0, dim = c(NROW(ar), nobs, nobs))
  fac <- 1 / (1 - ar^2)
  pow_ar <- as.list(rep(1, nobs + 1))
  for (i in seq_len(nobs)) {
    pow_ar[[i + 1]] <- ar^i
    out[, i, i] <- fac
    for (j in seq_len(i - 1)) { 
      out[, i, j] <- fac * pow_ar[[i - j + 1]]
      out[, j, i] <- out[, i, j]
    } 
  } 
  out
}

# compute MA1 correlation matrices
# @param ma MA1 autocorrelation samples
# @param nobs number of rows of the covariance matrix
# @return a numeric 'nsamples' x 'nobs' x 'nobs' array
get_cor_matrix_ma1 <- function(ma, nobs) {
  out <- array(0, dim = c(NROW(ma), nobs, nobs))
  gamma0 <- 1 + ma^2
  for (i in seq_len(nobs)) { 
    out[, i, i] <- gamma0
    if (i > 1) {
      out[, i, i - 1] <- ma
    }
    if (i < nobs) {
      out[, i, i + 1] <- ma
    }
  } 
  out 
}

# compute ARMA1 correlation matrices
# @param ar AR1 autocorrelation samples
# @param ma MA1 autocorrelation samples
# @param nobs number of rows of the covariance matrix
# @return a numeric 'nsamples' x 'nobs' x 'nobs' array
get_cor_matrix_arma1 <- function(ar, ma, nobs) {
  out <- array(0, dim = c(NROW(ar), nobs, nobs))
  fac <- 1 / (1 - ar^2)
  gamma0 <- 1 + ma^2 + 2 * ar * ma
  gamma <- as.list(rep(NA, nobs))
  gamma[[1]] <- (1 + ar * ma) * (ar + ma)
  for (i in seq_len(nobs)) {
    out[, i, i] <- fac * gamma0
    gamma[[i]] <- gamma[[1]] * ar^(i - 1)
    for (j in seq_len(i - 1)) { 
      out[, i, j] <- fac * gamma[[i - j]]
      out[, j, i] <- out[, i, j]
    } 
  } 
  out 
}

# compute compound symmetry correlation matrices
# @param cosy compund symmetry correlation samples
# @param nobs number of rows of the covariance matrix
# @return a numeric 'nsamples' x 'nobs' x 'nobs' array
get_cor_matrix_cosy <- function(cosy, nobs) {
  out <- array(0, dim = c(NROW(cosy), nobs, nobs))
  for (i in seq_len(nobs)) {
    out[, i, i] <- 1
    for (j in seq_len(i - 1)) { 
      out[, i, j] <- cosy
      out[, j, i] <- out[, i, j]
    } 
  } 
  out
}

# compute an identity correlation matrix
# @param nsamples number of posterior samples
# @param nobs number of rows of the covariance matrix
# @return a numeric 'nsamples' x 'nobs' x 'nobs' array
get_cor_matrix_ident <- function(nsamples, nobs) {
  out <- array(0, dim = c(nsamples, nobs, nobs))
  for (i in seq_len(nobs)) {
    out[, i, i] <- 1
  }
  out
}

# get samples of a distributional parameter
# @param x object to extract posterior samples from
# @param dpar name of the distributional parameter
# @param i the current observation number
# @param ilink should the inverse link function be applied?
#   if NULL the value is chosen internally
# @return 
#   If the parameter is predicted and i is NULL or 
#   length(i) > 1, an S x N matrix.
#   If the parameter it not predicted or length(i) == 1,
#   a vector of length S.
get_dpar <- function(draws, dpar, i = NULL, ilink = NULL) {
  stopifnot(is.brmsdraws(draws) || is.mvbrmsdraws(draws))
  x <- draws$dpars[[dpar]]
  stopifnot(!is.null(x))
  if (is.list(x)) {
    # compute samples of a predicted parameter
    out <- predictor(x, i = i, fdraws = draws)
    if (is.null(ilink)) {
      ilink <- apply_dpar_ilink(dpar, family = draws$family)
    }
    if (ilink) {
      out <- ilink(out, x$family$link)
    }
    if (length(i) == 1L) {
      out <- extract_col(out, 1)
    }
  } else if (!is.null(i) && !is.null(dim(x))) {
    out <- extract_col(x, i)
  } else {
    out <- x
  }
  if (dpar == "sigma") {
    out <- add_sigma_se(out, draws, i = i)
  }
  out
}

# get samples of a non-linear parameter
# @param x object to extract posterior samples from
# @param nlpar name of the non-linear parameter
# @param i the current observation number
# @return
#   If i is NULL or length(i) > 1: an S x N matrix
#   If length(i) == 1: a vector of length S
get_nlpar <- function(draws, nlpar, i = NULL) {
  stopifnot(is.brmsdraws(draws) || is.mvbrmsdraws(draws))
  x <- draws$nlpars[[nlpar]]
  stopifnot(!is.null(x))
  if (is.list(x)) {
    # compute samples of a predicted parameter
    out <- predictor(x, i = i, fdraws = draws)
    if (length(i) == 1L) {
      out <- extract_col(out, 1)
    }
  } else if (!is.null(i) && !is.null(dim(x))) {
    out <- extract_col(x, i)
  } else {
    out <- x
  }
  out
}

# get the mixing proportions of mixture models
get_theta <- function(draws, i = NULL) {
  stopifnot(is.brmsdraws(draws))
  if ("theta" %in% names(draws$dpars)) {
    # theta was not predicted; no need to call get_dpar
    theta <- draws$dpars$theta
  } else {
    # theta was predicted; apply softmax
    mix_family <- draws$family
    families <- family_names(mix_family)
    theta <- vector("list", length(families))
    for (j in seq_along(families)) {
      draws$family <- mix_family$mix[[j]]
      theta[[j]] <- as.matrix(get_dpar(draws, paste0("theta", j), i = i))
    }
    theta <- abind(theta, along = 3)
    for (n in seq_len(dim(theta)[2])) {
      theta[, n, ] <- softmax(theta[, n, ])
    }
    if (length(i) == 1L) {
      dim(theta) <- dim(theta)[c(1, 3)]
    }
  }
  theta
}

# get posterior samples of multivariate mean vectors
# only used in multivariate models with 'rescor'
get_Mu <- function(draws, i = NULL) {
  stopifnot(is.mvbrmsdraws(draws))
  Mu <- draws$mvpars$Mu
  if (is.null(Mu)) {
    Mu <- lapply(draws$resps, get_dpar, "mu", i = i)
    if (length(i) == 1L) {
      Mu <- do_call(cbind, Mu)
    } else {
      # keep correct dimension even if data has only 1 row
      Mu <- lapply(Mu, as.matrix)
      Mu <- do_call(abind, c(Mu, along = 3))
    }
  } else {
    stopifnot(!is.null(i))
    Mu <- extract_col(Mu, i)
  }
  Mu
}

# get posterior samples of residual covariance matrices
# only used in multivariate models with 'rescor'
get_Sigma <- function(draws, i = NULL) {
  stopifnot(is.mvbrmsdraws(draws))
  Sigma <- draws$mvpars$Sigma
  if (is.null(Sigma)) {
    stopifnot(!is.null(draws$mvpars$rescor))
    sigma <- lapply(draws$resps, get_dpar, "sigma", i = i)
    is_matrix <- ulapply(sigma, is.matrix)
    if (!any(is_matrix)) {
      # happens if length(i) == 1 or if no sigma was predicted
      sigma <- do_call(cbind, sigma)
      Sigma <- get_cov_matrix(sigma, draws$mvpars$rescor)
    } else {
      for (j in seq_along(sigma)) {
        # bring all sigmas to the same dimension
        if (!is_matrix[j]) {
          sigma[[j]] <- array(sigma[[j]], dim = dim_mu(draws))
        }
      }
      nsigma <- length(sigma)
      sigma <- abind(sigma, along = 3)
      Sigma <- array(dim = c(dim_mu(draws), nsigma, nsigma))
      for (n in seq_len(ncol(Sigma))) {
        Sigma[, n, , ] <- get_cov_matrix(sigma[, n, ], draws$mvpars$rescor)
      }
    }
  } else {
    stopifnot(!is.null(i))
    ldim <- length(dim(Sigma))
    stopifnot(ldim %in% 3:4)
    if (ldim == 4L) {
      Sigma <- extract_col(Sigma, i)
    }
  }
  Sigma
}

# extract user-defined standard errors
get_se <- function(draws, i = NULL) {
  stopifnot(is.brmsdraws(draws))
  se <- as.vector(draws$data[["se"]])
  if (!is.null(se)) {
    if (!is.null(i)) {
      se <- se[i]
    }
    if (length(se) > 1L) {
      dim <- c(draws$nsamples, length(se))
      se <- as_draws_matrix(se, dim = dim)
    }
  } else {
    se <- 0
  }
  se
}

# add user defined standard errors to 'sigma'
add_sigma_se <- function(sigma, draws, i = NULL) {
  needs_se <- "se" %in% names(draws$data) && 
    !isTRUE(attr(sigma, "se_added")) &&
    !isTRUE(grepl("_cov$", draws$family$fun))
  if (needs_se) {
    # 'se' will be incorporated directly into 'sigma'
    sigma <- sqrt(get_se(draws, i = i)^2 + sigma^2)
    # make sure not to add 'se' twice
    attr(sigma, "se_added") <- TRUE
  } 
  sigma
}

# return samples of ordinal thresholds for observation i
# @param draws a drawsl or drawsnl object
# @param i observation number
subset_thres <- function(draws, i) {
  thres <- draws$thres$thres
  Jthres <- draws$thres$Jthres
  if (!is.null(Jthres)) {
    thres <- thres[, Jthres[i, 1]:Jthres[i, 2], drop = FALSE]
  }
  thres
}

# helper function of 'get_dpar' to decide if
# the link function should be applied by default
apply_dpar_ilink <- function(dpar, family) {
  !((has_cat(family) || has_thres(family)) && dpar_class(dpar) == "mu") ||
    is.customfamily(family)
}

# insert zeros for the predictor term of the reference category
# in categorical-like models using the softmax response function
insert_refcat  <- function(eta, family) {
  stopifnot(is.matrix(eta), is.brmsfamily(family))
  if (!conv_cats_dpars(family) || isNA(family$refcat)) {
    return(eta)
  }
  # need to add zeros for the reference category
  zeros <- as.matrix(rep(0, nrow(eta)))
  if (is.null(family$refcat) || is.null(family$cats)) {
    # no information on the categories provided:
    # use the first category as the reference
    return(cbind(zeros, eta))
  }
  colnames(zeros) <- paste0("mu", family$refcat)
  iref <- match(family$refcat, family$cats)
  before <- seq_len(iref - 1)
  after <- setdiff(seq_cols(eta), before)
  cbind(eta[, before, drop = FALSE], zeros, eta[, after, drop = FALSE])
}

# validate the 'resp' argument of 'predict' and related methods
# @param resp response names to be validated
# @param x valid response names or brmsfit object to extract names from
# @param multiple allow multiple response variables?
# @return names of validated response variables
validate_resp <- function(resp, x, multiple = TRUE) {
  if (is.brmsfit(x)) {
    x <- parse_bf(x$formula)$responses
  }
  x <- as.character(x)
  if (!length(x)) {
    # resp is unused in univariate models
    return(NULL)
  }
  if (length(resp)) {
    resp <- as.character(resp)
    if (!all(resp %in% x)) {
      stop2("Invalid argument 'resp'. Valid response ",
            "variables are: ", collapse_comma(x))
    }
    if (!multiple) {
      resp <- as_one_character(resp)
    }
  } else {
    resp <- x
  }
  resp
}

# split '...' into a list of model objects and other arguments
# takes its argument names from parent.frame() 
# @param .... objects to split into model and non-model objects
# @param x object treated in the same way as '...'. Adding it is
#   necessary for substitute() to catch the name of the first 
#   argument passed to S3 methods.
# @param model_names optional names of the model objects  
# @param other: allow non-model arguments in '...'?
# @return
#   A list of arguments. All brmsfit objects are stored 
#   as a list in element 'models' unless 'other' is FALSE.
#   In the latter case just returns a list of models
split_dots <- function(x, ..., model_names = NULL, other = TRUE) {
  other <- as_one_logical(other)
  dots <- list(x, ...)
  names <- substitute(list(x, ...), env = parent.frame())[-1]
  names <- ulapply(names, deparse_combine)
  if (length(names)) {
    if (!length(names(dots))) {
      names(dots) <- names
    } else {
      has_no_name <- !nzchar(names(dots))
      names(dots)[has_no_name] <- names[has_no_name]
    }
  }
  is_brmsfit <- unlist(lapply(dots, is.brmsfit))
  models <- dots[is_brmsfit]
  models <- validate_models(models, model_names, names(models))
  out <- dots[!is_brmsfit]
  if (other) {
    out$models <- models
  } else {
    if (length(out)) {
      stop2("Only model objects can be passed to '...' for this method.")
    }
    out <- models
  }
  out
}

# reorder observations to be in the initial user-defined order
# currently only relevant for autocorrelation models 
# @param eta 'nsamples' x 'nobs' matrix
# @param old_order optional vector to retrieve the initial data order
# @param sort keep the new order as defined by the time-series?
# @return the 'eta' matrix with possibly reordered columns
reorder_obs <- function(eta, old_order = NULL, sort = FALSE) {
  stopifnot(length(dim(eta)) %in% c(2L, 3L))
  if (!is.null(old_order) && !sort) {
    if (isTRUE(length(old_order) == ncol(eta))) {
      if (length(dim(eta)) == 3L) {
        eta <- eta[, old_order, , drop = FALSE]   
      } else {
        eta <- eta[, old_order, drop = FALSE]   
      }
    } else {
      warning2("Cannot recover the original observation order.")
    }
  }
  eta
}

# extract argument names of a post-processing method
arg_names <- function(method) {
  opts <- c("posterior_predict", "pp_expect", "log_lik")
  method <- match.arg(method, opts)
  out <- names(formals(paste0(method, ".brmsfit")))
  c(out) <- names(formals(extract_draws.brmsfit))
  c(out) <- names(formals(validate_newdata))
  out <- unique(out)
  out <- setdiff(out, c("object", "x", "..."))
  out
}

# read a brmsfit object from a file
# @param file path to an rds file
# @return a brmsfit object or NULL
read_brmsfit <- function(file) {
  file <- check_brmsfit_file(file)
  x <- suppressWarnings(try(readRDS(file), silent = TRUE))
  if (!is(x, "try-error")) {
    if (!is.brmsfit(x)) {
      stop2("Object loaded via 'file' is not of class 'brmsfit'.")
    }
    x$file <- file
  } else {
    x <- NULL
  }
  x
}

# write a brmsfit object to a file
# @param x a brmsfit object
# @param file path to an rds file
# @return NULL
write_brmsfit <- function(x, file) {
  stopifnot(is.brmsfit(x))
  file <- check_brmsfit_file(file)
  x$file <- file
  saveRDS(x, file = file)
  invisible(NULL)
}

# check validity of file name to store a brmsfit object in
check_brmsfit_file <- function(file) {
  file <- as_one_character(file)
  file_ending <- tolower(get_matches("\\.[^\\.]+$", file))
  if (!isTRUE(file_ending == ".rds")) {
    file <- paste0(file, ".rds")
  }
  file
}

# add dummy samples to a brmsfit object for use in unit tests
# @param x a brmsfit object
# @param newpar name of the new parameter to add
# @param dim dimension of the new parameter
# @param dist name of the distribution from which to sample
# @param ... further arguments passed to r<dist>
# @return a brmsfit object including dummy samples of the new parameter
add_samples <- function(x, newpar, dim = numeric(0), dist = "norm", ...) {
  stopifnot(is.brmsfit(x))
  stopifnot(identical(dim, numeric(0)))
  newpar <- as_one_character(newpar)
  for (i in seq_along(x$fit@sim$samples)) {
    x$fit@sim$samples[[i]][[newpar]] <- 
      do_call(paste0("r", dist), list(x$fit@sim$iter, ...))
  }
  x$fit@sim$fnames_oi <- c(x$fit@sim$fnames_oi, newpar) 
  x$fit@sim$dims_oi[[newpar]] <- dim
  x$fit@sim$pars_oi <- names(x$fit@sim$dims_oi)
  x
}
