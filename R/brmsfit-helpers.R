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

get_all_group_vars <- function(x) {
  # extract names of all grouping variables
  if (is.brmsfit(x)) {
    x <- x$ranef
  }
  unique(ulapply(x$gcall, "[[", "groups"))
}

name_model <- function(family) {
  # create the name of the fitted stan model
  # Args:
  #   family: A family object
  if (!is.family(family)) {
    out <- "brms-model"
  } else {
    out <- paste(summary(family), "brms-model")
  }
  out
}

link <- function(x, link) {
  # apply a link function on x
  # Args:
  #   x: An arrary of arbitrary dimension
  #   link: a character string defining the link
  # Returns:
  #   an array of dimension dim(x) on which the link function was applied
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
  # Returns:
  #   an array of dimension dim(x) on which the inverse link function was applied
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
            paste0(type, "(", names[j], "," , names[i], ")"))
        } else {
          cornames <- c(cornames, 
            paste0(type, sep, names[j], sep, names[i]))
        }
      }
    }
  }
  cornames
}

get_valid_groups <- function(x) {
  # extract names of valid grouping variables or factors
  stopifnot(is.brmsfit(x))
  like_factor <- sapply(model.frame(x), is_like_factor)
  valid_groups <- c(
    names(model.frame(x))[like_factor],
    parse_time(x$autocor$formula)$group,
    x$ranef$group
  )
  unique(valid_groups[nzchar(valid_groups)])
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
  x <- do.call(apply, c(args, dots))
  if (is.null(dim(x))) {
    x <- matrix(x, dimnames = list(NULL, coef))
  } else if (coef == "quantile") {
    x <- aperm(x, length(dim(x)):1)
  }
  x 
}

get_summary <- function(samples, probs = c(0.025, 0.975), robust = FALSE) {
  # summarizes parameter samples based on mean, sd, and quantiles
  # Args: 
  #   samples: a matrix or data.frame containing the samples to be summarized. 
  #            rows are samples, columns are parameters
  #   probs: quantiles to be computed
  #   robust: return median and mas instead of mean and sd?
  #   keep_names: keep dimnames of the samples?
  # Returns:
  #   a N x C matrix where N is the number of observations and C 
  #   is equal to \code{length(probs) + 2}.
  if (robust) {
    coefs <- c("median", "mad", "quantile")
  } else {
    coefs <- c("mean", "sd", "quantile")
  }
  .get_summary <- function(samples) {
    do.call(cbind, lapply(
      coefs, get_estimate, samples = samples, 
      probs = probs, na.rm = TRUE
    ))
  }
  stopifnot(length(dim(samples)) %in% 2:3)
  if (length(dim(samples)) == 2L) {
    out <- .get_summary(samples)
    rownames(out) <- colnames(samples)
  } else if (length(dim(samples)) == 3L) {
    out <- lapply(array2list(samples), .get_summary)
    out <- abind(out, along = 3)
    dimnames(out)[c(1, 3)] <- dimnames(samples)[c(2, 3)]
  }
  colnames(out) <- c("Estimate", "Est.Error", paste0(probs * 100, "%ile"))
  out  
}

get_table <- function(samples, levels = sort(unique(as.numeric(samples)))) {
  # compute absolute frequencies for each column
  # Args:
  #   samples: a S x N matrix
  #   levels: all possible values in samples
  # Returns:
  #    a N x levels matrix containing relative frequencies of each level
  stopifnot(is.matrix(samples))
  out <- do.call(rbind, lapply(seq_len(ncol(samples)), 
    function(n) table(factor(samples[, n], levels = levels))
  ))
  # compute relative frequencies
  out <- out / sum(out[1, ])
  rownames(out) <- colnames(samples)
  colnames(out) <- paste0("P(Y = ", seq_len(ncol(out)), ")")
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
    for (i in 2:size) {
      for (j in 1:(i-1)) {
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
    size <- - 1 / 2 + sqrt(1 / 4 + 2 * ncol(cor)) + 1
    nsamples <- nrow(cor)
  } 
  stopifnot(is_wholenumber(size) && size > 1, nsamples > 0)
  out <- array(diag(1, size), dim = c(size, size, nsamples))
  out <- aperm(out, perm = c(3, 1, 2))
  if (length(cor)) {
    k <- 0 
    for (i in 2:size) {
      for (j in 1:(i-1)) {
        k = k + 1
        out[, j, i] <- out[, i, j] <- cor[, k]
      }
    }
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
  # get samples of an distributional parameter
  # Args:
  #   x: object to extract postarior samples from
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
    if (!is.null(x$nlpars)) {
      out <- nonlinear_predictor(x, i = i)
    } else {
      out <- linear_predictor(x, i = i)
    }
    if (is.null(ilink)) {
      ilink <- apply_dpar_ilink(dpar, family = draws$f)
    }
    if (ilink) {
      out <- ilink(out, x$f$link)
    }
  } else if (!is.null(i) && isTRUE(ncol(x) > 1L)) {
    out <- index_col(x, i)
  } else {
    out <- x
  }
  if (is.matrix(out) && ncol(out) == 1L) {
    out <- as.vector(out)
  }
  if (dpar == "sigma" && !grepl("_cov$", draws$f$fun)) {
    if (!isTRUE(attr(x, "se_added"))) {
      out <- sqrt(get_se(draws, i = i)^2 + out^2)
      # make sure not to add 'se' twice
      attr(out, "se_added") <- TRUE
    }
  } else if (dpar == "disc" && is.matrix(out)) {
    # 'disc' will be multiplied by a 3D array
    out <- array(out, dim = c(dim(out), draws$data$ncat - 1))
  }
  out
}

get_theta <- function(draws, i = NULL) {
  # get the mixing proportions of mixture models
  stopifnot(is.brmsdraws(draws))
  if ("theta" %in% names(draws$dpars)) {
    theta <- get_dpar(draws, "theta", i = i)
  } else {
    # theta was predicted; apply softmax
    mix_family <- draws$f
    families <- family_names(mix_family)
    theta <- vector("list", length(families))
    for (j in seq_along(families)) {
      draws$f <- mix_family$mix[[j]]
      theta[[j]] <- get_dpar(draws, paste0("theta", j), i = i)
    }
    theta <- do.call(abind, c(theta, along = 3))
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
      Mu <- do.call(cbind, Mu)
    } else {
      Mu <- do.call(abind, c(Mu, along = 3))
    }
  } else {
    # at this point we now that i is not NULL
    Mu <- index_col(Mu, i)
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
      sigma <- do.call(cbind, sigma)
      Sigma <- get_cov_matrix(sigma, draws$mvpars$rescor)
    } else {
      for (j in seq_along(sigma)) {
        # bring all sigmas to the same dimension
        if (!is_matrix[j]) {
          sigma[[j]] <- array(sigma[[j]], dim = dim_mu(draws))
        }
      }
      nsigma <- length(sigma)
      sigma <- do.call(abind, c(sigma, along = 3))
      Sigma <- array(dim = c(dim_mu(draws), nsigma, nsigma))
      for (n in seq_len(ncol(Sigma))) {
        Sigma[, n, , ] <- get_cov_matrix(sigma[, n, ], draws$mvpars$rescor)
      }
    }
  } else {
    # at this point we know that i is not NULL
    ldim <- length(dim(Sigma))
    stopifnot(ldim %in% 3:4)
    if (ldim == 4L) {
      Sigma <- index_col(Sigma, i)
    }
  }
  Sigma
}

get_se <- function(draws, i = NULL) {
  # extract user-defined standard errors
  # Args: see get_dpar
  stopifnot(is.brmsdraws(draws))
  se <- draws$data[["se"]]
  if (!is.null(se)) {
    if (!is.null(i)) {
      se <- se[i]
    } else {
      dim <- c(draws$nsamples, draws$nobs)
      se <- as_draws_matrix(se, dim = dim)
    }
  } else {
    se <- 0
  }
  se
}

index_col <- function(x, i) {
  # savely index columns without dropping other dimensions
  # Args:
  #   x: an array
  #   i: colum index
  ldim <- length(dim(x))
  stopifnot(ldim %in% 2:4)
  if (ldim == 2L) {
    out <- x[, i]
  } else {
    if (ldim == 3L) {
      out <- x[, i, ]
    } else {
      out <- x[, i, , ]
    }
    if (length(i) == 1L && dim(out) != dim(x)[-2]) {
      # some non-column dims were unintentionally dropped
      dim(out) <- dim(x)[-2]
    }
  }
  out
}

apply_dpar_ilink <- function(dpar, family) {
  # helper function of get_dpar to decide if
  # the link function should be applied by default
  # Returns:
  #   TRUE or FALSE
  if (is.mixfamily(family)) {
    dpar <- dpar_class(dpar) 
  }
  no_link_family <- family$family %in% 
    c("weibull", "categorical", "cumulative", 
      "sratio", "cratio", "acat")
  !(no_link_family && dpar == "mu")
}

choose_N <- function(draws) {
  # choose N to be used in predict and log_lik
  stopifnot(is.brmsdraws(draws) || is.mvbrmsdraws(draws))
  if (!is.null(draws$data$N_tg)) draws$data$N_tg else draws$nobs
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

validate_resp <- function(resp, valid_resps) {
  # validate the 'resp' argument of 'predict' and related methods
  if (length(resp)) {
    if (!all(resp %in% valid_resps)) {
      stop2("Invalid argument 'resp'. Valid response ",
            "variables are: ", collapse_comma(valid_resps))
    }
  } else {
    resp <- valid_resps
  }
  resp
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
    N <- length(old_order)
    if (ncol(eta) %% N != 0) {
      # for compatibility with MV models fitted before brms 1.0.0
      stopifnot(N %% ncol(eta) == 0)
      old_order <- old_order[seq_len(ncol(eta))]
    }
    if (N < ncol(eta)) {
      # should occur for multivariate models only
      nresp <- ncol(eta) / N
      old_order <- rep(old_order, nresp)
      old_order <- old_order + rep(0:(nresp - 1) * N, each = N)
    }
    if (length(dim(eta)) == 3L) {
      eta <- eta[, old_order, , drop = FALSE]   
    } else {
      eta <- eta[, old_order, drop = FALSE]   
    }
  }
  colnames(eta) <- NULL
  eta
}

fixef_pars <- function() {
  # regex to extract population-level coefficients
  "^b(|cs|mo|me|m)_"
}

default_plot_pars <- function() {
  # list all parameter classes to be included in plots by default
  c(fixef_pars(), "^sd_", "^cor_", "^sigma_", "^rescor_", 
    paste0("^", dpars(), "[[:digit:]]*$"), "^delta$",
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
      pars <- all_pars[all_pars %in% pars]
    } else {
      pars <- all_pars[apply(sapply(pars, grepl, x = all_pars, ...), 1, any)]
    }
  } else {
    pars <- na_value
  }
  pars
}

find_names <- function(x) {
  # find all valid object names in a string 
  # Args:
  #   x: a character string
  # Notes:
  #   Does not use the R parser itself to allow for double points, 
  #   square brackets and kommas at the end of names.
  #   currently only used in 'hypothesis_internal'
  # Returns:
  #   all valid variable names within the string
  if (!is.character(x) || length(x) > 1) {
    stop2("Argument 'x' must be a character string of length one.")
  }
  x <- gsub("[[:space:]]", "", x)
  regex_all <- "([^([:digit:]|[:punct:])]|\\.)[[:alnum:]_\\.\\:]*"
  regex_all <- paste0(regex_all, "(\\[[^],]+(,[^],]+)*\\])?")
  pos_all <- gregexpr(regex_all, x)[[1]]
  regex_fun <- "([^([:digit:]|[:punct:])]|\\.)[[:alnum:]_\\.]*\\("
  pos_fun <- gregexpr(regex_fun, x)[[1]]
  pos_decnum <- gregexpr("\\.[[:digit:]]+", x)[[1]]
  pos_var <- list(rmMatch(pos_all, pos_fun, pos_decnum))
  unlist(regmatches(x, pos_var))
}

evidence_ratio <- function(x, cut = 0, wsign = c("equal", "less", "greater"), 
                           prior_samples = NULL, pow = 12, ...) {
  # compute the evidence ratio between two disjunct hypotheses
  # Args:
  #   x: posterior samples 
  #   cut: the cut point between the two hypotheses
  #   wsign: direction of the hypothesis
  #   prior_samples: optional prior samples for undirected hypothesis
  #   pow: influences the accuracy of the density
  #   ...: optional arguments passed to density.default
  # Returns:
  #   the evidence ratio of the two hypothesis
  wsign <- match.arg(wsign)
  if (wsign == "equal") {
    if (is.null(prior_samples)) {
      out <- NA
    } else {
      dots <- list(...)
      dots <- dots[names(dots) %in% names(formals("density.default"))]
      args <- c(list(n = 2^pow), dots)
      eval_dens <- function(x) {
        # evaluate density of x at cut
        from <- min(x)
        to <- max(x)
        if (from > cut) {
          from <- cut - sd(x) / 4
        } else if (to < cut) {
          to <- cut + sd(x) / 4
        }
        dens <- do.call(density, c(nlist(x, from, to), args))
        spline(dens$x, dens$y, xout = cut)$y
      }
      out <- eval_dens(x) / eval_dens(prior_samples)
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

hypothesis_internal <- function(x, hypothesis, class = "", alpha = 0.05, ...) {
  # internal function to evaluate hypotheses
  # Args:
  #   x: the primary object passed to the hypothesis method;
  #      Needs to be a brmsfit object or coercible to a data.frame
  #   hypothesis: Vector of character strings containing the hypotheses
  #   class: prefix of the parameters in the hypotheses
  #   alpha: alpha-level
  # Returns:
  #   an object of class 'brmshypothesis'
  if (!is.character(hypothesis)) {
    stop2("Argument 'hypothesis' must be a character vector.")
  }
  if (!is.character(class) || length(class) != 1L) {
    stop2("Argument 'class' must be a single character string.")
  }
  if (length(alpha) != 1L || alpha < 0 || alpha > 1) {
    stop2("Argument 'alpha' must be a single value in [0,1].")
  }
  
  .eval_hypothesis <- function(h) {
    stopifnot(length(h) == 1L && is.character(h))
    # parse hypothesis string
    h <- gsub("[ \t\r\n]", "", h)
    sign <- get_matches("=|<|>", h)
    lr <- get_matches("[^=<>]+", h)
    if (length(sign) != 1L || length(lr) != 2L) {
      stop2("Every hypothesis must be of the form 'left (= OR < OR >) right'.")
    }
    h <- paste0("(", lr[1], ")")
    h <- paste0(h, ifelse(lr[2] != "0", paste0("-(", lr[2], ")"), ""))
    varsH <- unique(find_names(h))
    parsH <- paste0(class, varsH)
    miss_pars <- setdiff(parsH, pars)
    if (length(miss_pars)) {
      miss_pars <- collapse_comma(miss_pars)
      stop2("Some parameters cannot be found in the model: \n", miss_pars)
    }
    # rename hypothesis for correct evaluation
    h_renamed <- rename(h, c(":", "[", "]", ","),  c("___", ".", ".", ".."))
    # get posterior and prior samples
    symbols <- c(paste0("^", class), ":", "\\[", "\\]", ",")
    subs <- c("", "___", ".", ".", "..")
    samples <- posterior_samples(x, pars = parsH, exact_match = TRUE)
    names(samples) <- rename(names(samples), symbols, subs, fixed = FALSE)
    samples <- as.matrix(eval2(h_renamed, samples))
    prior_samples <- prior_samples(x, pars = parsH, fixed = TRUE)
    if (!is.null(prior_samples) && ncol(prior_samples) == length(varsH)) {
      names(prior_samples) <- rename(
        names(prior_samples), symbols, subs, fixed = FALSE
      )
      prior_samples <- as.matrix(eval2(h_renamed, prior_samples))
    } else {
      prior_samples <- NULL
    }
    # summarize hypothesis
    wsign <- switch(sign, "=" = "equal", "<" = "less", ">" = "greater")
    probs <- switch(sign, 
      "=" = c(alpha / 2, 1 - alpha / 2), 
      "<" = c(0, 1 - alpha), ">" = c(alpha, 1)
    )
    sm <- lapply(
      c("mean", "sd", "quantile", "evidence_ratio"), 
      get_estimate, samples = samples, probs = probs, 
      wsign = wsign, prior_samples = prior_samples
    )
    sm <- as.data.frame(matrix(unlist(sm), nrow = 1))
    if (sign == "<") {
      sm[1, 3] <- -Inf
    } else if (sign == ">") {
      sm[1, 4] <- Inf
    }
    sm <- cbind(sm, ifelse(!(sm[1, 3] <= 0 && 0 <= sm[1, 4]), '*', ''))
    rownames(sm) <- paste(h, sign, "0")
    cl <- (1 - alpha) * 100
    colnames(sm) <- c(
      "Estimate", "Est.Error", paste0("l-", cl, "% CI"),
      paste0("u-", cl, "% CI"), "Evid.Ratio", "Star"
    )
    if (is.null(prior_samples)) {
      prior_samples <- as.matrix(rep(NA, nrow(samples)))
    }
    return(nlist(summary = sm, samples, prior_samples))
  }
  
  pars <- parnames(x)[grepl("^", class, parnames(x))]
  hlist <- lapply(hypothesis, .eval_hypothesis)
  hs <- do.call(rbind, lapply(hlist, function(h) h$summary))
  samples <- as.data.frame(
    do.call(cbind, lapply(hlist, function(h) h$samples))
  )
  prior_samples <- as.data.frame(
    do.call(cbind, lapply(hlist, function(h) h$prior_samples))
  )
  names(samples) <- names(prior_samples) <- paste0("H", seq_along(hlist))
  class <- sub("_+$", "", class)
  out <- nlist(hypothesis = hs, samples, prior_samples, class, alpha)
  class(out) <- "brmshypothesis"
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
      do.call(paste0("r", dist), list(x$fit@sim$iter, ...))
  }
  x$fit@sim$fnames_oi <- c(x$fit@sim$fnames_oi, newpar) 
  x$fit@sim$dims_oi[[newpar]] <- dim
  x$fit@sim$pars_oi <- names(x$fit@sim$dims_oi)
  x
}
