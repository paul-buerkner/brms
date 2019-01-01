contains_samples <- function(x) {
  if (!(is.brmsfit(x) && length(x$fit@sim))) {
    stop2("The model does not contain posterior samples.")
  }
  invisible(TRUE)
}

algorithm <- function(x) {
  stopifnot(is.brmsfit(x))
  if (is.null(x$algorithm)) "sampling"
  else x$algorithm
}

is_mv <- function(x) {
  stopifnot(is.brmsfit(x))
  is.mvbrmsformula(x$formula)
}

stopifnot_resp <- function(x, resp = NULL) {
  if (is_mv(x) && is.null(resp)) {
    stop2("Argument 'resp' is required when applying ", 
          "this method to a multivariate model.")
  }
  invisible(NULL)
}

link <- function(x, link) {
  # apply a link function on x
  # Args:
  #   x: An arrary of arbitrary dimension
  #   link: a character string defining the link
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
    stop2("Link '", link, "' not supported.")
  )
}

ilink <- function(x, link) {
  # apply the inverse link function on x
  # Args:
  #   x: An arrary of arbitrary dimension
  #   link: a character string defining the link
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
    stop2("Link '", link, "' not supported.")
  )
}

subset_samples <- function(x, subset = NULL, nsamples = NULL) {
  # generate integers indicating subsets of the posterior samples
  stopifnot(is.brmsfit(x))
  if (is.null(subset) && !is.null(nsamples)) {
    subset <- sample(nsamples(x), nsamples)
  }
  subset
}

get_rnames <- function(ranef, group = NULL, bylevels = NULL) {
  # extract names of group-level effects
  # Args:
  #  ranef: output of tidy_ranef()
  #  group: optinal name of a grouping factor for
  #    which to extract effect names
  #  bylevels: optional names of 'by' levels for 
  #    which to extract effect names
  stopifnot(is.data.frame(ranef))
  if (!is.null(group)) {
    group <- as_one_character(group)
    ranef <- subset2(ranef, group = group)
  }
  stopifnot(length(unique(ranef$group)) == 1L)
  out <- paste0(usc(combine_prefix(ranef), "suffix"), ranef$coef)
  if (isTRUE(nzchar(ranef$by[1]))) {
    if (!is.null(bylevels)) {
      stopifnot(all(bylevels %in% ranef$bylevels[[1]]))
    } else {
      bylevels <- ranef$bylevels[[1]]
    }
    bylabels <- paste0(ranef$by[1], bylevels)
    out <- outer(out, bylabels, paste, sep = ":")
  }
  out
}

get_cornames <- function(names, type = "cor", brackets = TRUE, sep = "__") {
  # get correlation names as combinations of variable names
  # Args:
  #   names: the variable names 
  #   type: character string to be put in front of the returned strings
  #   brackets: should the correlation names contain brackets 
  #            or underscores as seperators
  #   sep: character string to separate names;
  #        only used if brackets == FALSE
  # Returns: 
  #   correlation names based on the variable names 
  #   passed to the names argument
  cornames <- NULL
  if (length(names) > 1) {
    for (i in 2:length(names)) {
      for (j in 1:(i-1)) {
        if (brackets) {
          cornames <- c(cornames, 
            paste0(type, "(", names[j], "," , names[i], ")")
          )
        } else {
          cornames <- c(cornames, 
            paste0(type, sep, names[j], sep, names[i])
          )
        }
      }
    }
  }
  cornames
}

get_cat_vars <- function(x) {
  # extract names of categorical variables in the model
  stopifnot(is.brmsfit(x))
  like_factor <- sapply(model.frame(x), is_like_factor)
  valid_groups <- c(
    names(model.frame(x))[like_factor],
    get_group_vars(x)
  )
  unique(valid_groups[nzchar(valid_groups)])
}

get_levels <- function(...) {
  # extract list of levels with one element per grouping factor
  # Args:
  #   ...: object with a level attribute
  dots <- list(...)
  out <- vector("list", length(dots))
  for (i in seq_along(out)) {
    levels <- attr(dots[[i]], "levels", exact = TRUE)
    if (is.list(levels)) {
      stopifnot(!is.null(names(levels)))
      out[[i]] <- as.list(levels)
    } else if (!is.null(levels)) {
      stopifnot(isTRUE(nzchar(names(dots)[i])))
      out[[i]] <- setNames(list(levels), names(dots)[[i]])
    }
  }
  out <- unlist(out, recursive = FALSE)
  out[!duplicated(names(out))]
}

get_estimate <- function(coef, samples, margin = 2, ...) {
  # calculate estimates over posterior samples 
  # Args:
  #   coef: coefficient to be applied on the samples (e.g., "mean")
  #   samples: the samples over which to apply coef
  #   margin: see apply
  #   to.array: logical; should the result be transformed 
  #             into an array of increased dimension?
  #   ...: additional arguments passed to get(coef)
  # Returns: 
  #   typically a matrix with colnames(samples) as colnames
  dots <- list(...)
  args <- list(X = samples, MARGIN = margin, FUN = coef)
  fun_args <- names(formals(coef))
  if (!"..." %in% fun_args) {
    dots <- dots[names(dots) %in% fun_args]
  }
  x <- do_call(apply, c(args, dots))
  if (is.null(dim(x))) {
    x <- matrix(x, dimnames = list(NULL, coef))
  } else if (coef == "quantile") {
    x <- aperm(x, length(dim(x)):1)
  }
  x 
}

#' Table Creation for Posterior Samples
#' 
#' Create a table for unique values of posterior samples. 
#' This is usually only useful when summarizing predictions 
#' of ordinal models.
#' 
#' @param x A matrix of posterior samples where rows 
#'   indicate samples and columns indicate parameters. 
#' @param levels Optional values of possible posterior values.
#'   Defaults to all unique values in \code{x}.
#' 
#' @return A matrix where rows indicate parameters 
#'  and columns indicate the unique values of 
#'  posterior samples.
#'  
#' @examples 
#' \dontrun{
#' fit <- brm(rating ~ period + carry + treat, 
#'            data = inhaler, family = cumulative())
#' pr <- predict(fit, summary = FALSE)
#' posterior_table(pr)
#' }
#'  
#' @export
posterior_table <- function(x, levels = NULL) {
  x <- as.matrix(x)
  if (is.null(levels)) {
    levels <- sort(unique(as.vector(x)))
  }
  xlevels <- attr(x, "levels")
  if (length(xlevels) != length(levels)) {
    xlevels <- levels
  }
  out <- lapply(seq_len(ncol(x)), 
    function(n) table(factor(x[, n], levels = levels))
  )
  out <- do_call(rbind, out)
  # compute relative frequencies
  out <- out / sum(out[1, ])
  rownames(out) <- colnames(x)
  colnames(out) <- paste0("P(Y = ", xlevels, ")")
  out
}

get_cov_matrix <- function(sd, cor = NULL) {
  # compute covariance matrices based on correlation and sd samples
  # Args:
  #   sd: samples of standard deviations
  #   cor: samples of correlations
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

get_cor_matrix <- function(cor, size = NULL, nsamples = NULL) {
  # compute correlation matrices based on correlation samples
  # Args:
  #   cor: samples of correlations
  #   size: optional size of the desired correlation matrix
  #     ignored is cor is specified
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

get_cov_matrix_arma <- function(draws, obs) {
  # helper function to compute ARMA covariance matrices
  # currently, only ARMA1 processes are implemented
  # Args:
  #   draws: a brmsdraws object
  #   obs: observations for which to compute the covariance matrix
  nobs <- length(obs)
  # make sure not to add 'se' twice
  se <- draws$data$se[obs]
  draws$data$se <- NULL
  sigma <- get_dpar(draws, "sigma", i = obs)
  ar <- as.numeric(draws$ac$ar)
  ma <- as.numeric(draws$ac$ma)
  if (length(ar) && !length(ma)) {
    out <- get_cov_matrix_ar1(ar, sigma, nobs, se)
  } else if (!length(ar) && length(ma)) {
    out <- get_cov_matrix_ma1(ma, sigma, nobs, se)
  } else if (length(ar) && length(ma)) {
    out <- get_cov_matrix_arma1(ar, ma, sigma, nobs, se)
  } else {
    out <- get_cov_matrix_ident(sigma, nobs, se)
  }
  out
}

get_cov_matrix_ar1 <- function(ar, sigma, nrows, se = NULL) {
  # compute the covariance matrix for an AR1 process
  # Args: 
  #   ar: AR1 autocorrelation samples
  #   sigma: standard deviation samples of the AR1 process
  #   se: optional user defined standard errors
  #   nrows: number of rows of the covariance matrix
  # Returns:
  #   An nsamples x nrows x nrows AR1 covariance array (!)
  sigma <- as.matrix(sigma)
  se2 <- if (length(se)) se^2 else 0
  mat <- array(diag(se2, nrows), dim = c(nrows, nrows, nrow(sigma)))
  mat <- aperm(mat, perm = c(3, 1, 2))
  sigma2_adjusted <- sigma^2 / (1 - ar^2)
  pow_ar <- as.list(rep(1, nrows + 1))
  for (i in seq_len(nrows)) {
    pow_ar[[i + 1]] <- ar^i
    mat[, i, i] <- mat[, i, i] + sigma2_adjusted
    for (j in seq_len(i - 1)) { 
      mat[, i, j] <- sigma2_adjusted * pow_ar[[i - j + 1]]
      mat[, j, i] <- mat[, i, j]
    } 
  } 
  mat
}

get_cov_matrix_ma1 <- function(ma, sigma, nrows, se = NULL) {
  # compute the covariance matrix for an MA1 process
  # Args: 
  #   ma: MA1 autocorrelation samples
  #   sigma: standard deviation samples of the AR1 process
  #   se: optional user defined standard errors
  #   nrows: number of rows of the covariance matrix
  # Returns:
  #   An nsamples x nrows x nrows MA1 covariance array (!)
  sigma <- as.matrix(sigma)
  se2 <- if (length(se)) se^2 else 0
  mat <- array(diag(se2, nrows), dim = c(nrows, nrows, nrow(sigma)))
  mat <- aperm(mat, perm = c(3, 1, 2))
  sigma2 <- sigma^2
  sigma2_adjusted <- sigma2 * (1 + ma^2)
  sigma2_times_ma <- sigma2 * ma
  for (i in seq_len(nrows)) { 
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

get_cov_matrix_arma1 <- function(ar, ma, sigma, nrows, se = NULL) {
  # compute the covariance matrix for an AR1 process
  # Args: 
  #   ar: AR1 autocorrelation sample
  #   ma: MA1 autocorrelation sample
  #   sigma: standard deviation samples of the AR1 process
  #   se: optional user defined standard errors
  #   nrows: number of rows of the covariance matrix
  # Returns:
  #   An nsamples x nrows x nrows ARMA1 covariance array (!)
  sigma <- as.matrix(sigma)
  se2 <- if (length(se)) se^2 else 0
  mat <- array(diag(se2, nrows), dim = c(nrows, nrows, nrow(sigma)))
  mat <- aperm(mat, perm = c(3, 1, 2))
  sigma2_adjusted <- sigma^2 / (1 - ar^2)
  gamma0 <- 1 + ma^2 + 2 * ar * ma
  gamma <- as.list(rep(NA, nrows))
  gamma[[1]] <- (1 + ar * ma) * (ar + ma)
  for (i in seq_len(nrows)) {
    mat[, i, i] <- mat[, i, i] + sigma2_adjusted * gamma0
    gamma[[i]] <- gamma[[1]] * ar^(i - 1)
    for (j in seq_len(i - 1)) { 
      mat[, i, j] <- sigma2_adjusted * gamma[[i - j]]
      mat[, j, i] <- mat[, i, j]
    } 
  } 
  mat 
}

get_cov_matrix_ident <- function(sigma, nrows, se = 0) {
  # compute a variance matrix without including ARMA parameters
  # only used for ARMA covariance models when incl_autocor = FALSE
  # Args:
  #   sigma: standard deviation samples
  #   se: user defined standard errors (may be 0)
  #   nrows: number of rows of the covariance matrix
  # Returns:
  #   An nsamples x nrows x nrows sigma array
  sigma <- as.matrix(sigma)
  se2 <- if (length(se)) se^2 else 0
  mat <- array(diag(se2, nrows), dim = c(nrows, nrows, nrow(sigma)))
  mat <- aperm(mat, perm = c(3, 1, 2))
  sigma2 <- sigma^2
  for (i in seq_len(nrows)) {
    mat[, i, i] <- mat[, i, i] + sigma2
  }
  mat
}

get_dpar <- function(draws, dpar, i = NULL, ilink = NULL) {
  # get samples of a distributional parameter
  # Args:
  #   x: object to extract posterior samples from
  #   dpar: name of the distributional parameter
  #   i: the current observation number
  #   ilink: should the inverse link function be applied?
  #     if NULL the value is chosen internally
  # Returns:
  #   If the parameter is predicted and i is NULL or 
  #   length(i) > 1, an S x N matrix.
  #   If the parameter it not predicted or length(i) == 1,
  #   a vector of length S.
  stopifnot(is.brmsdraws(draws) || is.mvbrmsdraws(draws))
  x <- draws$dpars[[dpar]]
  stopifnot(!is.null(x))
  if (is.list(x)) {
    # compute samples of a predicted parameter
    out <- predictor(x, i = i, fdraws = draws)
    if (is.null(ilink)) {
      ilink <- apply_dpar_ilink(dpar, family = draws$f)
    }
    if (ilink) {
      out <- ilink(out, x$f$link)
    }
    if (length(i) == 1L) {
      out <- extract_col(out, 1)
    }
  } else if (!is.null(i) && !is.null(dim(x))) {
    out <- extract_col(x, i)
  } else {
    out <- x
  }
  if (dpar == "sigma" && !isTRUE(grepl("_cov$", draws$f$fun))) {
    # 'se' will be incorporated directly into 'sigma'
    if (!isTRUE(attr(x, "se_added")) && "se" %in% names(draws$data)) {
      out <- sqrt(get_se(draws, i = i)^2 + out^2)
      # make sure not to add 'se' twice
      attr(out, "se_added") <- TRUE
    }
  } else if (dpar == "disc" && is.null(i) && is.matrix(out)) {
    # 'disc' will be multiplied by a 3D array
    out <- array(out, dim = c(dim(out), draws$data$ncat - 1))
  }
  out
}

get_nlpar <- function(draws, nlpar, i = NULL) {
  # get samples of a non-linear parameter
  # Args:
  #   x: object to extract posterior samples from
  #   nlpar: name of the non-linear parameter
  #   i: the current observation number
  # Returns:
  #   If i is NULL or length(i) > 1: an S x N matrix
  #   If length(i) == 1: a vector of length S
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

get_theta <- function(draws, i = NULL) {
  # get the mixing proportions of mixture models
  stopifnot(is.brmsdraws(draws))
  if ("theta" %in% names(draws$dpars)) {
    # theta was not predicted; no need to call get_dpar
    theta <- draws$dpars$theta
  } else {
    # theta was predicted; apply softmax
    mix_family <- draws$f
    families <- family_names(mix_family)
    theta <- vector("list", length(families))
    for (j in seq_along(families)) {
      draws$f <- mix_family$mix[[j]]
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

get_se <- function(draws, i = NULL) {
  # extract user-defined standard errors
  # Args: see get_dpar
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

apply_dpar_ilink <- function(dpar, family) {
  # helper function of get_dpar to decide if
  # the link function should be applied by default
  # Returns:
  #   TRUE or FALSE
  if (is.mixfamily(family)) {
    dpar <- dpar_class(dpar) 
  }
  !(is_ordinal(family) && dpar == "mu" || is_categorical(family))
}

choose_N <- function(draws) {
  # choose N to be used in predict and log_lik
  stopifnot(is.brmsdraws(draws) || is.mvbrmsdraws(draws))
  if (!is.null(draws$ac$N_tg)) draws$ac$N_tg else draws$nobs
}

prepare_family <- function(x) {
  # prepare for calling family specific log_lik / predict functions
  stopifnot(is.brmsformula(x) || is.brmsterms(x))
  family <- x$family
  if (use_cov(x$autocor)) {
    family$fun <- paste0(family$family, "_cov")
  } else if (is.cor_sar(x$autocor)) {
    if (identical(x$autocor$type, "lag")) {
      family$fun <- paste0(family$family, "_lagsar")
    } else if (identical(x$autocor$type, "error")) {
      family$fun <- paste0(family$family, "_errorsar")
    }
  } else if (is.cor_fixed(x$autocor)) {
    family$fun <- paste0(family$family, "_fixed")
  } else {
    family$fun <- family$family
  }
  family
}

validate_resp <- function(resp, x, multiple = TRUE) {
  # validate the 'resp' argument of 'predict' and related methods
  # Args:
  #   resp: response names to be validated
  #   x: valid response names or brmsfit object to extract names from
  #   multiple: allow multiple response variables?
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

split_dots <- function(x, ..., model_names = NULL, other = TRUE) {
  # split '...' into a list of model objects and other arguments
  # takes its argument names from parent.frame() 
  # Args:
  #   ....: objects to split into model and non-model objects
  #   x: object treated in the same way as '...'. Adding it is
  #      necessary for substitute() to catch the name of the first 
  #      argument passed to S3 methods.
  #   model_names: optional names of the model objects  
  #   other: allow non-model arguments in '...'?
  # Returns
  #   A list of arguments. All brmsfit objects are stored 
  #   as a list in element 'models' unless 'other' is FALSE.
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
  if (is.null(out$subset) && !is.null(out$nsamples)) {
    out$subset <- sample(nsamples(models[[1]]), out$nsamples)
    out$nsamples <- NULL
  }
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

validate_weights <- function(weights, models, control = list()) {
  # validate weights passed to model averaging functions
  # Args: see pp_average.brmsfit
  if (!is.numeric(weights)) {
    weight_args <- c(unname(models), control)
    weight_args$weights <- weights
    weights <- do_call(model_weights, weight_args)
  } else {
    if (length(weights) != length(models)) {
      stop2("If numeric, 'weights' must have the same length ",
            "as the number of models.")
    }
    if (any(weights < 0)) {
      stop2("If numeric, 'weights' must be positive.")
    }
  }
  weights / sum(weights)
}

reorder_obs <- function(eta, old_order = NULL, sort = FALSE) {
  # reorder observations to be in the initial user-defined order
  # currently only relevant for autocorrelation models 
  # Args:
  #   eta: Nsamples x Nobs matrix
  #   old_order: optional vector to retrieve the initial data order
  #   sort: keep the new order as defined by the time-series?
  # Returns:
  #   eta with possibly reordered columns
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

fixef_pars <- function() {
  # regex to extract population-level coefficients
  types <- c("", "s", "cs", "sp", "mo", "me", "mi", "m")
  types <- paste0("(", types, ")", collapse = "|")
  paste0("^b(", types, ")_")
}

default_plot_pars <- function(family) {
  # list all parameter classes to be included in plots by default
  c(fixef_pars(), "^sd_", "^cor_", "^sigma_", "^rescor_", 
    paste0("^", valid_dpars(family), "$"), "^delta$",
    "^theta", "^ar", "^ma", "^arr", "^lagsar", "^errorsar", 
    "^car", "^sdcar", "^sigmaLL", "^sds_", "^sdgp_", "^lscale_")
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
  if (!(anyNA(pars) || is.character(pars))) {
    stop2("Argument 'pars' must be NA or a character vector.")
  }
  if (!anyNA(pars)) {
    if (exact_match) {
      out <- intersect(pars, all_pars)
    } else {
      out <- vector("list", length(pars))
      for (i in seq_along(pars)) {
        out[[i]] <- all_pars[grepl(pars[i], all_pars, ...)]
      }
      out <- unique(unlist(out))
    }
  } else {
    out <- na_value
  }
  out
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
  #   a brmsfit object with samples of a new parameter
  stopifnot(is.brmsfit(x))
  stopifnot(identical(dim, numeric(0)))
  for (i in seq_along(x$fit@sim$samples)) {
    x$fit@sim$samples[[i]][[newpar]] <- 
      do_call(paste0("r", dist), list(x$fit@sim$iter, ...))
  }
  x$fit@sim$fnames_oi <- c(x$fit@sim$fnames_oi, newpar) 
  x$fit@sim$dims_oi[[newpar]] <- dim
  x$fit@sim$pars_oi <- names(x$fit@sim$dims_oi)
  x
}
