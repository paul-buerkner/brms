contains_draws <- function(x) {
  if (!(is.brmsfit(x) && length(x$fit@sim))) {
    stop2("The model does not contain posterior draws.")
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
# @param x an array of arbitrary dimension
# @param link character string defining the link
link <- function(x, link) {
  switch(link,
    identity = x,
    log = log(x),
    logm1 = logm1(x),
    log1p = log1p(x),
    inverse = 1 / x,
    sqrt = sqrt(x),
    "1/mu^2" = 1 / x^2,
    tan_half = tan(x / 2),
    logit = logit(x),
    probit = qnorm(x),
    cauchit = qcauchy(x),
    cloglog = cloglog(x),
    probit_approx = qnorm(x),
    softplus = log_expm1(x),
    squareplus = (x^2 - 1) / x,
    softit = softit(x),
    stop2("Link '", link, "' is not supported.")
  )
}

# apply an inverse link function
# @param x an array of arbitrary dimension
# @param link a character string defining the link
inv_link <- function(x, link) {
  switch(link,
    identity = x,
    log = exp(x),
    logm1 = expp1(x),
    log1p = expm1(x),
    inverse = 1 / x,
    sqrt = x^2,
    "1/mu^2" = 1 / sqrt(x),
    tan_half = 2 * atan(x),
    logit = inv_logit(x),
    probit = pnorm(x),
    cauchit = pcauchy(x),
    cloglog = inv_cloglog(x),
    probit_approx = pnorm(x),
    softplus = log1p_exp(x),
    squareplus = (x + sqrt(x^2 + 4)) / 2,
    softit = inv_softit(x),
    stop2("Link '", link, "' is not supported.")
  )
}

# log CDF for unit interval link functions
# @param x an array of arbitrary dimension
# @param link a character string defining the link
log_cdf <- function(x, link) {
  switch(link,
    logit = log_inv_logit(x),
    probit = pnorm(x, log.p = TRUE),
    cauchit = pcauchy(x, log.p = TRUE),
    cloglog = log1m_exp(-exp(x)),
    probit_approx = pnorm(x, log.p = TRUE),
    softit = log_inv_softit(x),
    stop2("Link '", link, "' is not supported.")
  )
}

# log CCDF for unit interval link functions
# @param x an array of arbitrary dimension
# @param link a character string defining the link
log_ccdf <- function(x, link) {
  switch(link,
    logit = log1m_inv_logit(x),
    probit = pnorm(x, log.p = TRUE, lower.tail = FALSE),
    cauchit = pcauchy(x, log.p = TRUE, lower.tail = FALSE),
    cloglog = -exp(x),
    probit_approx = pnorm(x, log.p = TRUE, lower.tail = FALSE),
    softit = log1m_inv_softit(x),
    stop2("Link '", link, "' is not supported.")
  )
}

# validate integers indicating which draws to subset
validate_draw_ids <- function(x, draw_ids = NULL, ndraws = NULL) {
  ndraws_total <- ndraws(x)
  if (is.null(draw_ids) && !is.null(ndraws)) {
    ndraws <- as_one_integer(ndraws)
    if (ndraws < 1 || ndraws > ndraws_total) {
      stop2("Argument 'ndraws' should be between 1 and ",
            "the maximum number of draws (", ndraws_total, ").")
    }
    draw_ids <- sample(seq_len(ndraws_total), ndraws)
  }
  if (!is.null(draw_ids)) {
    draw_ids <- as.integer(draw_ids)
    if (any(draw_ids < 1L) || any(draw_ids > ndraws_total)) {
      stop2("Some 'draw_ids' indices are out of range.")
    }
  }
  draw_ids
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

# covariance matrices based on correlation and SD draws
# @param sd matrix of draws of standard deviations
# @param cor matrix of draws of correlations
get_cov_matrix <- function(sd, cor = NULL) {
  sd <- as.matrix(sd)
  stopifnot(all(sd >= 0))
  ndraws <- nrow(sd)
  size <- ncol(sd)
  out <- array(diag(1, size), dim = c(size, size, ndraws))
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

# correlation matrices based on correlation draws
# @param cor draws of correlations
# @param size optional size of the desired correlation matrix;
#   ignored is 'cor' is specified
# @param ndraws optional number of posterior draws;
#   ignored is 'cor' is specified
get_cor_matrix <- function(cor, size = NULL, ndraws = NULL) {
  if (length(cor)) {
    cor <- as.matrix(cor)
    size <- -1 / 2 + sqrt(1 / 4 + 2 * ncol(cor)) + 1
    ndraws <- nrow(cor)
  }
  size <- as_one_numeric(size)
  ndraws <- as_one_numeric(ndraws)
  stopifnot(is_wholenumber(size) && size > 0)
  stopifnot(is_wholenumber(ndraws) && ndraws > 0)
  out <- array(diag(1, size), dim = c(size, size, ndraws))
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
# @param prep a brmsprep object
# @param obs observations for which to compute the covariance matrix
# @param latent compute covariance matrix for latent residuals?
get_cov_matrix_ac <- function(prep, obs = NULL, latent = FALSE) {
  if (is.null(obs)) {
    obs <- seq_len(prep$nobs)
  }
  nobs <- length(obs)
  ndraws <- prep$ndraws
  acef <- prep$ac$acef
  # prepare correlations
  if (has_ac_class(acef, "arma")) {
    ar <- as.numeric(prep$ac$ar)
    ma <- as.numeric(prep$ac$ma)
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
    cosy <- as.numeric(prep$ac$cosy)
    cor <- get_cor_matrix_cosy(cosy, nobs)
  } else if (has_ac_class(acef, "fcor")) {
    cor <- get_cor_matrix_fcor(prep$ac$Mfcor, ndraws)
  } else {
    cor <- get_cor_matrix_ident(ndraws, nobs)
  }
  # prepare known standard errors
  if (!is.null(prep$data$se)) {
    se2 <- prep$data$se[obs]^2
    se2 <- array(diag(se2, nobs), dim = c(nobs, nobs, ndraws))
    se2 <- aperm(se2, perm = c(3, 1, 2))
    # make sure not to add 'se' twice
    prep$data$se <- NULL
  } else {
    se2 <- rep(0, nobs)
  }
  # prepare residual standard deviations
  if (latent) {
    sigma2 <- as.numeric(prep$ac$sderr)^2
  } else {
    sigma <- get_dpar(prep, "sigma", i = obs)
    if (NCOL(sigma) > 1L) {
      # sigma varies across observations
      sigma2 <- array(dim = c(ndraws, nobs, nobs))
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
# @param ar AR1 autocorrelation draws
# @param nobs number of rows of the covariance matrix
# @return a numeric 'ndraws' x 'nobs' x 'nobs' array
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
# @param ma MA1 autocorrelation draws
# @param nobs number of rows of the covariance matrix
# @return a numeric 'ndraws' x 'nobs' x 'nobs' array
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
# @param ar AR1 autocorrelation draws
# @param ma MA1 autocorrelation draws
# @param nobs number of rows of the covariance matrix
# @return a numeric 'ndraws' x 'nobs' x 'nobs' array
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
# @param cosy compund symmetry correlation draws
# @param nobs number of rows of the covariance matrix
# @return a numeric 'ndraws' x 'nobs' x 'nobs' array
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

# prepare a fixed correlation matrix
# @param Mfcor correlation matrix to be prepared
# @param ndraws number of posterior draws
# @return a numeric 'ndraws' x 'nobs' x 'nobs' array
get_cor_matrix_fcor <- function(Mfcor, ndraws) {
  out <- array(Mfcor, dim = c(dim(Mfcor), ndraws))
  aperm(out, c(3, 1, 2))
}

# compute an identity correlation matrix
# @param ndraws number of posterior draws
# @param nobs number of rows of the covariance matrix
# @return a numeric 'ndraws' x 'nobs' x 'nobs' array
get_cor_matrix_ident <- function(ndraws, nobs) {
  out <- array(0, dim = c(ndraws, nobs, nobs))
  for (i in seq_len(nobs)) {
    out[, i, i] <- 1
  }
  out
}

#' Draws of a Distributional Parameter
#'
#' Get draws of a distributional parameter from a \code{brmsprep} or
#' \code{mvbrmsprep} object. This function is primarily useful when developing
#' custom families or packages depending on \pkg{brms}.
#' This function lets callers easily handle both the case when the
#' distributional parameter is predicted directly, via a (non-)linear
#' predictor or fixed to a constant. See the vignette
#' \code{vignette("brms_customfamilies")} for an example use case.
#'
#' @param prep A 'brmsprep' or 'mvbrmsprep' object created by
#'   \code{\link[brms:prepare_predictions.brmsfit]{prepare_predictions}}.
#' @param dpar Name of the distributional parameter.
#' @param i The observation numbers for which predictions shall be extracted.
#'   If \code{NULL} (the default), all observation will be extracted.
#'   Ignored if \code{dpar} is not predicted.
#' @param inv_link Should the inverse link function be applied?
#'   If \code{NULL} (the default), the value is chosen internally.
#'   In particular, \code{inv_link} is \code{TRUE} by default for custom
#'   families.
#' @return
#'   If the parameter is predicted and \code{i} is \code{NULL} or
#'   \code{length(i) > 1}, an \code{S x N} matrix. If the parameter it not
#'   predicted or \code{length(i) == 1}, a vector of length \code{S}. Here
#'   \code{S} is the number of draws and \code{N} is the number of
#'   observations or length of \code{i} if specified.
#'
#' @examples
#' \dontrun{
#' posterior_predict_my_dist <- function(i, prep, ...) {
#'   mu <- brms::get_dpar(prep, "mu", i = i)
#'   mypar <- brms::get_dpar(prep, "mypar", i = i)
#'   my_rng(mu, mypar)
#' }
#' }
#'
#' @export
get_dpar <- function(prep, dpar, i = NULL, inv_link = NULL) {
  stopifnot(is.brmsprep(prep) || is.mvbrmsprep(prep))
  dpar <- as_one_character(dpar)
  x <- prep$dpars[[dpar]]
  stopifnot(!is.null(x))
  if (is.list(x)) {
    # compute draws of a predicted parameter
    out <- predictor(x, i = i, fprep = prep)
    if (is.null(inv_link)) {
      inv_link <- apply_dpar_inv_link(dpar, family = prep$family)
    } else {
      inv_link <- as_one_logical(inv_link)
    }
    if (inv_link) {
      out <- inv_link(out, x$family$link)
    }
    if (length(i) == 1L) {
      out <- slice_col(out, 1)
    }
  } else if (!is.null(i) && !is.null(dim(x))) {
    out <- slice_col(x, i)
  } else {
    out <- x
  }
  out
}

# get draws of a non-linear parameter
# @param x object to extract posterior draws from
# @param nlpar name of the non-linear parameter
# @param i the current observation number
# @return
#   If i is NULL or length(i) > 1: an S x N matrix
#   If length(i) == 1: a vector of length S
get_nlpar <- function(prep, nlpar, i = NULL) {
  stopifnot(is.brmsprep(prep) || is.mvbrmsprep(prep))
  x <- prep$nlpars[[nlpar]]
  stopifnot(!is.null(x))
  if (is.list(x)) {
    # compute draws of a predicted parameter
    out <- predictor(x, i = i, fprep = prep)
    if (length(i) == 1L) {
      out <- slice_col(out, 1)
    }
  } else if (!is.null(i) && !is.null(dim(x))) {
    out <- slice_col(x, i)
  } else {
    out <- x
  }
  out
}

# get the mixing proportions of mixture models
get_theta <- function(prep, i = NULL) {
  stopifnot(is.brmsprep(prep))
  if ("theta" %in% names(prep$dpars)) {
    # theta was not predicted; no need to call get_dpar
    theta <- prep$dpars$theta
  } else {
    # theta was predicted; apply softmax
    mix_family <- prep$family
    families <- family_names(mix_family)
    theta <- vector("list", length(families))
    for (j in seq_along(families)) {
      prep$family <- mix_family$mix[[j]]
      theta[[j]] <- as.matrix(get_dpar(prep, paste0("theta", j), i = i))
    }
    theta <- abind(theta, along = 3)
    for (n in seq_len(dim(theta)[2])) {
      theta[, n, ] <- exp(log_softmax(slice(theta, 2, n)))
    }
    if (length(i) == 1L) {
      dim(theta) <- dim(theta)[c(1, 3)]
    }
  }
  theta
}

# get posterior draws of multivariate mean vectors
# only used in multivariate models with 'rescor'
# and in univariate models with multiple 'mu' pars such as logistic_normal
get_Mu <- function(prep, i = NULL) {
  is_mv <- is.mvbrmsprep(prep)
  if (is_mv) {
    Mu <- prep$mvpars$Mu
  } else {
    stopifnot(is.brmsprep(prep))
    Mu <- prep$dpars$Mu
  }
  if (!is.null(Mu)) {
    stopifnot(!is.null(i))
    Mu <- slice_col(Mu, i)
    return(Mu)
  }
  if (is_mv) {
    Mu <- lapply(prep$resps, get_dpar, "mu", i = i)
  } else {
    mu_dpars <- str_subset(names(prep$dpars), "^mu")
    Mu <- lapply(mu_dpars, get_dpar, prep = prep, i = i)
  }
  if (length(i) == 1L) {
    Mu <- do_call(cbind, Mu)
  } else {
    # keep correct dimension even if data has only 1 row
    Mu <- lapply(Mu, as.matrix)
    Mu <- abind::abind(Mu, along = 3)
  }
  Mu
}

# get posterior draws of residual covariance matrices
# only used in multivariate models with 'rescor'
# and in univariate models with multiple 'mu' pars such as logistic_normal
get_Sigma <- function(prep, i = NULL, cor_name = NULL) {
  is_mv <- is.mvbrmsprep(prep)
  if (is_mv) {
    cor_name <- "rescor"
    Sigma <- prep$mvpars$Sigma
  } else {
    stopifnot(is.brmsprep(prep))
    cor_name <- as_one_character(cor_name)
    Sigma <- prep$dpars$Sigma
  }
  if (!is.null(Sigma)) {
    # already computed before
    stopifnot(!is.null(i))
    ldim <- length(dim(Sigma))
    stopifnot(ldim %in% 3:4)
    if (ldim == 4L) {
      Sigma <- slice_col(Sigma, i)
    }
    return(Sigma)
  }
  if (is_mv) {
    cors <- prep$mvpars[[cor_name]]
    sigma <- named_list(names(prep$resps))
    for (j in seq_along(sigma)) {
      sigma[[j]] <- get_dpar(prep$resps[[j]], "sigma", i = i)
      sigma[[j]] <- add_sigma_se(sigma[[j]], prep$resps[[j]], i = i)
    }
  } else {
    cors <- prep$dpars[[cor_name]]
    sigma_names <- str_subset(names(prep$dpars), "^sigma")
    sigma <- named_list(sigma_names)
    for (j in seq_along(sigma)) {
      sigma[[j]] <- get_dpar(prep, sigma_names[j], i = i)
      sigma[[j]] <- add_sigma_se(sigma[[j]], prep, i = i)
    }
  }
  is_matrix <- ulapply(sigma, is.matrix)
  if (!any(is_matrix)) {
    # happens if length(i) == 1 or if no sigma was predicted
    sigma <- do_call(cbind, sigma)
    Sigma <- get_cov_matrix(sigma, cors)
  } else {
    for (j in seq_along(sigma)) {
      # bring all sigmas to the same dimension
      if (!is_matrix[j]) {
        sigma[[j]] <- array(sigma[[j]], dim = dim_mu(prep))
      }
    }
    nsigma <- length(sigma)
    sigma <- abind(sigma, along = 3)
    Sigma <- array(dim = c(dim_mu(prep), nsigma, nsigma))
    for (n in seq_len(ncol(Sigma))) {
      Sigma[, n, , ] <- get_cov_matrix(slice(sigma, 2, n), cors)
    }
  }
  Sigma
}

# extract user-defined standard errors
get_se <- function(prep, i = NULL) {
  stopifnot(is.brmsprep(prep))
  se <- as.vector(prep$data[["se"]])
  if (!is.null(se)) {
    if (!is.null(i)) {
      se <- se[i]
    }
    if (length(se) > 1L) {
      dim <- c(prep$ndraws, length(se))
      se <- data2draws(se, dim = dim)
    }
  } else {
    se <- 0
  }
  se
}

# add user defined standard errors to 'sigma'
# @param sigma draws of the 'sigma' parameter
add_sigma_se <- function(sigma, prep, i = NULL) {
  if ("se" %in% names(prep$data)) {
    se <- get_se(prep, i = i)
    sigma <- sqrt(se^2 + sigma^2)
  }
  sigma
}

# extract user-defined rate denominators
get_rate_denom <- function(prep, i = NULL) {
  stopifnot(is.brmsprep(prep))
  denom <- as.vector(prep$data[["denom"]])
  if (!is.null(denom)) {
    if (!is.null(i)) {
      denom <- denom[i]
    }
    if (length(denom) > 1L) {
      dim <- c(prep$ndraws, length(denom))
      denom <- data2draws(denom, dim = dim)
    }
  } else {
    denom <- 1
  }
  denom
}

# multiply a parameter with the 'rate' denominator
# @param dpar draws of the distributional parameter
multiply_dpar_rate_denom <- function(dpar, prep, i = NULL) {
  if ("denom" %in% names(prep$data)) {
    denom <- get_rate_denom(prep, i = i)
    dpar <- dpar * denom
  }
  dpar
}

# return draws of ordinal thresholds for observation i
# @param prep a bprepl or bprepnl object
# @param i observation number
subset_thres <- function(prep, i) {
  thres <- prep$thres$thres
  Jthres <- prep$thres$Jthres
  if (!is.null(Jthres)) {
    thres <- thres[, Jthres[i, 1]:Jthres[i, 2], drop = FALSE]
  }
  thres
}

# helper function of 'get_dpar' to decide if
# the link function should be applied directly
apply_dpar_inv_link <- function(dpar, family) {
  !(has_joint_link(family) && dpar_class(dpar, family) == "mu")
}

# insert zeros for the predictor term of the reference category
# in categorical-like models using the softmax response function
insert_refcat <- function(eta, refcat = 1) {
  stopifnot(is.array(eta))
  refcat <- as_one_integer(refcat)
  # need to add zeros for the reference category
  ndim <- length(dim(eta))
  dim_noncat <- dim(eta)[-ndim]
  zeros_arr <- array(0, dim = c(dim_noncat, 1))
  before <- seq_len(refcat - 1)
  after <- setdiff(seq_dim(eta, ndim), before)
  abind::abind(
    slice(eta, ndim, before, drop = FALSE),
    zeros_arr,
    slice(eta, ndim, after, drop = FALSE)
  )
}

# validate the 'resp' argument of 'predict' and related methods
# @param resp response names to be validated
# @param x valid response names or brmsfit object to extract names from
# @param multiple allow multiple response variables?
# @return names of validated response variables
validate_resp <- function(resp, x, multiple = TRUE) {
  if (is.brmsfit(x)) {
    x <- brmsterms(x$formula)$responses
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
# @param eta 'ndraws' x 'nobs' matrix or array
# @param old_order optional vector to retrieve the initial data order
# @param sort keep the new order as defined by the time-series?
# @return the 'eta' matrix with possibly reordered columns
reorder_obs <- function(eta, old_order = NULL, sort = FALSE) {
  stopifnot(length(dim(eta)) %in% c(2L, 3L))
  if (is.null(old_order) || sort) {
    return(eta)
  }
  stopifnot(length(old_order) == NCOL(eta))
  p(eta, old_order, row = FALSE)
}

# update .MISC environment of the stanfit object
# allows to call log_prob and other C++ using methods
# on objects not created in the current R session
# or objects created via another backend
update_misc_env <- function(x, only_windows = FALSE) {
  stopifnot(is.brmsfit(x))
  only_windows <- as_one_logical(only_windows)
  if (!has_rstan_model(x)) {
    x <- add_rstan_model(x)
  } else if (os_is_windows() || !only_windows) {
    # TODO: detect when updating .MISC is not required
    # TODO: find a more efficient way to update .MISC
    old_backend <- x$backend
    x$backend <- "rstan"
    x$fit@.MISC <- suppressMessages(brm(fit = x, chains = 0))$fit@.MISC
    x$backend <- old_backend
  }
  x
}

#' Add compiled \pkg{rstan} models to \code{brmsfit} objects
#'
#' Compile a \code{\link[rstan:stanmodel-class]{stanmodel}} and add
#' it to a \code{brmsfit} object. This enables some advanced functionality
#' of \pkg{rstan}, most notably \code{\link[rstan:log_prob]{log_prob}}
#' and friends, to be used with brms models fitted with other Stan backends.
#'
#' @param x A \code{brmsfit} object to be updated.
#' @param overwrite Logical. If \code{TRUE}, overwrite any existing
#' \code{\link[rstan:stanmodel-class]{stanmodel}}. Defaults to \code{FALSE}.
#'
#' @return A (possibly updated) \code{brmsfit} object.
#'
#' @export
add_rstan_model <- function(x, overwrite = FALSE) {
  stopifnot(is.brmsfit(x))
  overwrite <- as_one_logical(overwrite)
  if (!has_rstan_model(x) || overwrite) {
    message("Recompiling the model with 'rstan'")
    # threading is not yet supported by rstan and needs to be deactivated
    stanfit <- suppressMessages(rstan::stan(
      model_code = stancode(x, threads = threading(), backend = "rstan"),
      data = standata(x), chains = 0
    ))
    x$fit@stanmodel <- stanfit@stanmodel
    x$fit@.MISC <- stanfit@.MISC
    message("Recompilation done")
  }
  x
}

# does the model have a non-empty rstan 'stanmodel'
# that can be used for 'log_prob' and friends?
has_rstan_model <- function(x) {
  stopifnot(is.brmsfit(x))
  isTRUE(nzchar(x$fit@stanmodel@model_cpp$model_cppname)) &&
    length(ls(pos = x$fit@.MISC)) > 0
}

# extract argument names of a post-processing method
arg_names <- function(method) {
  opts <- c("posterior_predict", "posterior_epred", "log_lik")
  method <- match.arg(method, opts)
  out <- names(formals(paste0(method, ".brmsfit")))
  c(out) <- names(formals(prepare_predictions.brmsfit))
  c(out) <- names(formals(validate_newdata))
  out <- unique(out)
  out <- setdiff(out, c("object", "x", "..."))
  out
}

# validate 'cores' argument for use in post-processing functions
validate_cores_post_processing <- function(cores) {
  if (is.null(cores)) {
    if (os_is_windows()) {
      # multi cores often leads to a slowdown on windows
      # in post-processing functions as discussed in #1129
      cores <- 1L
    } else {
      cores <- getOption("mc.cores", 1L)
    }
  }
  cores <- as_one_integer(cores)
  if (cores < 1L) {
    cores <- 1L
  }
  cores
}

#' Check if cached fit can be used.
#'
#' Checks whether a given cached fit can be used without refitting when
#' \code{file_refit = "on_change"} is used.
#' This function is internal and exposed only to facilitate debugging problems
#' with cached fits. The function may change or be removed in future versions
#' and scripts should not use it.
#'
#' @param fit Old \code{brmsfit} object (e.g., loaded from file).
#' @param sdata New Stan data (result of a call to \code{\link{make_standata}}).
#'   Pass \code{NULL} to avoid this data check.
#' @param scode New Stan code (result of a call to \code{\link{make_stancode}}).
#'   Pass \code{NULL} to avoid this code check.
#' @param data New data to check consistency of factor level names.
#'   Pass \code{NULL} to avoid this data check.
#' @param algorithm New algorithm. Pass \code{NULL} to avoid algorithm check.
#' @param silent Logical. If \code{TRUE}, no messages will be given.
#' @param verbose Logical. If \code{TRUE} detailed report of the differences
#'   is printed to the console.
#' @return A boolean indicating whether a refit is needed.
#'
#' @details
#' Use with \code{verbose = TRUE} to get additional info on how the stored
#' fit differs from the given data and code.
#'
#' @export
#' @keywords internal
brmsfit_needs_refit <- function(fit, sdata = NULL, scode = NULL, data = NULL,
                                algorithm = NULL, silent = FALSE,
                                verbose = FALSE) {
  stopifnot(is.brmsfit(fit))
  silent <- as_one_logical(silent)
  verbose <- as_one_logical(verbose)
  if (!is.null(scode)) {
    scode <- as_one_character(scode)
    cached_scode <- stancode(fit)
  }
  if (!is.null(sdata)) {
    stopifnot(is.list(sdata))
    cached_sdata <- standata(fit)
  }
  if (!is.null(data)) {
    stopifnot(is.data.frame(data))
    cached_data <- fit$data
  }
  if (!is.null(algorithm)) {
    algorithm <- as_one_character(algorithm)
    stopifnot(!is.null(fit$algorithm))
  }

  refit <- FALSE
  if (!is.null(scode)) {
    if (normalize_stancode(scode) != normalize_stancode(cached_scode)) {
      if (!silent) {
        message("Stan code has changed beyond whitespace/comments.")
        if (verbose) {
          require_package("diffobj")
          print(diffobj::diffChr(scode, cached_scode, format = "ansi8"))
        }
      }
      refit <- TRUE
    }
  }
  if (!is.null(sdata)) {
    sdata_equality <- all.equal(sdata, cached_sdata, check.attributes = FALSE)
    if (!isTRUE(sdata_equality)) {
      if (!silent) {
        message("The processed data for Stan has changed.")
        if (verbose) {
          print(sdata_equality)
        }
      }
      refit <- TRUE
    }
  }
  if (!is.null(data)) {
    # check consistency of factor names
    # as they are only stored as attributes in sdata (#1128)
    factor_level_message <- FALSE
    for (var in names(cached_data)) {
      if (is_like_factor(cached_data[[var]])) {
        cached_levels <- levels(factor(cached_data[[var]]))
        new_levels <- levels(factor(data[[var]]))
        if (!is_equal(cached_levels, new_levels)) {
          if (!silent) {
            factor_level_message <- TRUE
            if (verbose) {
              cat(paste0(
                "Names of factor levels have changed for variable '", var, "' ",
                "with cached levels (", collapse_comma(cached_levels), ") ",
                "but new levels (", collapse_comma(new_levels), ").\n"
              ))
            }
          }
          refit <- TRUE
          if (!verbose) {
            # no need to check all variables if we trigger a refit anyway
            break
          }
        }
      }
    }
    if (factor_level_message) {
      message("Names of factor levels have changed.")
    }
  }
  if (!is.null(algorithm)) {
    if (algorithm != fit$algorithm) {
      if (!silent) {
        message("Algorithm has changed from '", fit$algorithm,
                "' to '", algorithm, "'.\n")
      }
      refit <- TRUE
    }
  }
  refit
}

# read a brmsfit object from a file
# @param file path to an rds file
# @return a brmsfit object or NULL
read_brmsfit <- function(file) {
  file <- check_brmsfit_file(file)
  dir <- dirname(file)
  if (!dir.exists(dir)) {
    stop2(
      "The directory '", dir, "' does not exist. Please choose an ",
      "existing directory where the model can be saved after fitting."
    )
  }
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
  invisible(x)
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

# check if a function requires an old default setting
# only used to ensure backwards compatibility
# @param version brms version in which the change to the default was made
# @return TRUE or FALSE
require_old_default <- function(version) {
  version <- as.package_version(version)
  brmsfit_version <- getOption(".brmsfit_version")
  isTRUE(brmsfit_version < version)
}

# add dummy draws to a brmsfit object for use in unit tests
# @param x a brmsfit object
# @param newpar name of the new parameter to add
# @param dim dimension of the new parameter
# @param dist name of the distribution from which to sample
# @param ... further arguments passed to r<dist>
# @return a brmsfit object including dummy draws of the new parameter
add_dummy_draws <- function(x, newpar, dim = numeric(0), dist = "norm", ...) {
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
