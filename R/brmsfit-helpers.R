array2list <- function(x) {
  # convert array to list of elements with reduced dimension
  # Args: 
  #   x: an arrary of dimension d
  # Returns: 
  #   A list of arrays of dimension d-1
  if (is.null(dim(x))) {
    stop("Argument 'x' has no dimension.")
  }
  ndim <- length(dim(x))
  l <- list(length = dim(x)[ndim])
  ind <- collapse(rep(",", ndim - 1))
  for (i in seq_len(dim(x)[ndim])) {
    l[[i]] <- eval(parse(text = paste0("x[", ind, i, "]"))) 
  }
  names(l) <- dimnames(x)[[ndim]]
  l
}

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

#' Retructure Old \code{brmsfit} Objects
#'
#' Restructure old \code{brmsfit} objects to work with 
#' the latest \pkg{brms} version. This function is called
#' internally when applying post-processing methods.
#' However, in order to avoid uncessary run time caused
#' by the restructuring, I recommend explicitely calling
#' \code{restructure} once per model after updating \pkg{brms}.
#' 
#' @param x An object of class \code{brmsfit}.
#' @param rstr_summary Logical; If \code{TRUE}, the cached summary
#'   stored by \pkg{rstan} is restructured as well.
#'   
#' @return A \code{brmsfit} object compatible with the latest version
#'   of \pkg{brms}.
#'   
#' @export
restructure <- function(x, rstr_summary = FALSE) {
  stopifnot(is.brmsfit(x))
  if (isTRUE(attr(x, "restructured"))) {
    return(x)  # already restructured
  }
  if (is.null(x$version)) {
    # this is the latest version without saving the version number
    x$version <- list(brms = package_version("0.9.1"))
  } else if (is.package_version(x$version)) {
    # also added the rstan version in brms 1.5.0
    x$version <- list(brms = x$version)
  }
  version <- x$version$brms
  if (version < utils::packageVersion("brms")) {
    # element 'nonlinear' deprecated as of brms > 0.9.1
    # element 'partial' deprecated as of brms > 0.8.0
    x$formula <- SW(amend_formula(
      formula(x), data = model.frame(x), family = family(x), 
      partial = x$partial, nonlinear = x$nonlinear
    ))
    x$nonlinear <- x$partial <- NULL
    x$formula[["old_mv"]] <- is_old_mv(x)
    bterms <- parse_bf(formula(x), family = family(x))
    x$ranef <- tidy_ranef(bterms, model.frame(x))
    if ("prior_frame" %in% class(x$prior)) {
      class(x$prior) <- c("brmsprior", "data.frame") 
    }
    if (is(x$autocor, "cov_fixed")) {
      # deprecated as of brms 1.4.0
      class(x$autocor) <- "cor_fixed"
    }
    if (version <= "0.9.1") {
      # update gaussian("log") to lognormal() family
      nresp <- length(bterms$response)
      if (is_old_lognormal(x$family, nresp = nresp, version = version)) {
        object$family <- object$formula$family <- lognormal()
      }
    }
    if (version <= "0.10.0.9000") {
      if (length(bterms$auxpars$mu$nlpars)) {
        # nlpar and group have changed positions
        change <- change_old_re(x$ranef, pars = parnames(x),
                                dims = x$fit@sim$dims_oi)
        x <- do_renaming(x, change)
      }
    }
    if (version < "1.0.0") {
      # double underscores were added to group-level parameters
      change <- change_old_re2(x$ranef, pars = parnames(x),
                               dims = x$fit@sim$dims_oi)
      x <- do_renaming(x, change)
    }
    if (version <= "1.0.1") {
      # names of spline parameters had to be changed after
      # allowing for multiple covariates in one spline term
      change <- change_old_sm(bterms, pars = parnames(x),
                              dims = x$fit@sim$dims_oi)
      x <- do_renaming(x, change)
    }
    if (version <= "1.2.0") {
      x$ranef$type[x$ranef$type == "mono"] <- "mo"
      x$ranef$type[x$ranef$type == "cse"] <- "cs"
    }
    stan_env <- attributes(x$fit)$.MISC
    if (rstr_summary && exists("summary", stan_env)) {
      stan_summary <- get("summary", stan_env)
      old_parnames <- rownames(stan_summary$msd)
      if (!identical(old_parnames, parnames(x))) {
        # do not rename parameters in the cached summary
        # just let rstan compute the summary again
        remove("summary", pos = stan_env)
      }
    }
  }
  structure(x, restructured = TRUE)
}

first_greater <- function(A, target, i = 1) {
  # find the first element in A that is greater than target
  # Args: 
  #   A: a matrix
  #   target: a vector of length nrow(A)
  #   i: column of A being checked first
  # Returns: 
  #   A vector of the same length as target containing the column ids 
  #   where A[,i] was first greater than target
  ifelse(target <= A[, i] | ncol(A) == i, i, first_greater(A, target, i + 1))
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

prepare_conditions <- function(x, conditions = NULL, effects = NULL,
                               re_formula = NA, rsv_vars = NULL) {
  # prepare marginal conditions
  # Args:
  #   x: an object of class 'brmsfit'
  #   conditions: optional data.frame containing user defined conditions
  #   effects: see marginal_effects
  #   re_formula: see marginal_effects
  #   rsv_vars: names of reserved variables
  # Returns:
  #   A data.frame with (possibly updated) conditions
  mf <- model.frame(x)
  new_formula <- update_re_terms(formula(x), re_formula = re_formula)
  bterms <- parse_bf(new_formula, family = family(x))
  re <- get_re(bterms)
  req_vars <- c(
    lapply(get_effect(bterms, "fe"), rhs), 
    lapply(get_effect(bterms, "mo"), rhs),
    lapply(get_effect(bterms, "me"), rhs),
    lapply(get_effect(bterms, "sm"), rhs),
    lapply(get_effect(bterms, "cs"), rhs),
    lapply(get_effect(bterms, "offset"), rhs),
    re$form, lapply(re$gcall, "[[", "weightvars"),
    bterms$adforms[c("se", "disp", "trials", "cat")],
    bterms$auxpars$mu$covars
  )
  req_vars <- unique(ulapply(req_vars, all.vars))
  req_vars <- setdiff(req_vars, rsv_vars)
  if (is.null(conditions)) {
    conditions <- as.data.frame(as.list(rep(NA, length(req_vars))))
    names(conditions) <- req_vars
  } else {
    conditions <- as.data.frame(conditions)
    if (!nrow(conditions)) {
      stop2("Argument 'conditions' must have a least one row.")
    }
    if (any(duplicated(rownames(conditions)))) {
      stop2("Row names of 'conditions' should be unique.")
    }
    conditions <- unique(conditions)
    req_vars <- setdiff(req_vars, names(conditions))
  }
  trial_vars <- all.vars(bterms$adforms$trials)
  if (length(trial_vars)) {
    write_msg <- any(ulapply(trial_vars, function(x) 
      !isTRUE(x %in% names(conditions)) || anyNA(conditions[[x]])
    ))
    if (write_msg) {
      message("Using the median number of trials by ", 
              "default if not specified otherwise.")
    }
  }
  # use default values for unspecified variables
  int_vars <- rmNULL(bterms$adforms[c("trials", "cat")])
  int_vars <- c(int_vars, get_effect(bterms, "mo"))
  int_vars <- unique(ulapply(int_vars, all.vars))
  for (v in req_vars) {
    if (is.numeric(mf[[v]])) {
      if (v %in% int_vars) {
        conditions[[v]] <- round(median(mf[[v]]))
      } else {
        conditions[[v]] <- mean(mf[[v]])
      }
    } else {
      # use reference category
      levels <- attr(as.factor(mf[[v]]), "levels")
      conditions[[v]] <- factor(
        levels[1], levels = levels, ordered = is.ordered(mf[[v]])
      )
    }
  }
  unused_vars <- setdiff(names(conditions), all.vars(bterms$allvars))
  if (length(unused_vars)) {
    warning2(
      "The following variables in 'conditions' are not ", 
      "part of the model:\n", collapse_comma(unused_vars)
    )
  }
  amend_newdata(
    conditions, fit = x, re_formula = re_formula,
    allow_new_levels = TRUE, incl_autocor = FALSE, 
    return_standata = FALSE
  )
}

prepare_marg_data <- function(data, conditions, int_vars = NULL,
                              surface = FALSE, resolution = 100) {
  # prepare data to be used in marginal_effects
  # Args:
  #  data: data.frame containing only data of the predictors of interest
  #  conditions: see argument 'conditions' of marginal_effects
  #  int_vars: names of variables being treated as integers
  #  surface: generate surface plots later on?
  #  resolution: number of distinct points at which to evaluate
  #              the predictors of interest
  effects <- names(data)
  stopifnot(length(effects) %in% c(1L, 2L))
  pred_types <- ifelse(ulapply(data, is.numeric), "numeric", "factor")
  # numeric effects should come first
  new_order <- order(pred_types, decreasing = TRUE)
  effects <- effects[new_order]
  pred_types <- pred_types[new_order]
  mono <- effects %in% int_vars
  if (pred_types[1] == "numeric") {
    min1 <- min(data[, effects[1]])
    max1 <- max(data[, effects[1]])
    if (mono[1]) {
      values <- seq(min1, max1, by = 1)
    } else {
      values <- seq(min1, max1, length.out = resolution)
    }
  }
  if (length(effects) == 2L) {
    if (pred_types[1] == "numeric") {
      values <- setNames(list(values, NA), effects)
      if (pred_types[2] == "numeric") {
        if (surface) {
          min2 <- min(data[, effects[2]])
          max2 <- max(data[, effects[2]])
          if (mono[2]) {
            values[[2]] <- seq(min2, max2, by = 1)
          } else {
            values[[2]] <- seq(min2, max2, length.out = resolution)
          }
        } else {
          if (mono[2]) {
            median2 <- median(data[, effects[2]])
            mad2 <- mad(data[, effects[2]])
            values[[2]] <- round((-1:1) * mad2 + median2)
          } else {
            mean2 <- mean(data[, effects[2]])
            sd2 <- sd(data[, effects[2]])
            values[[2]] <- (-1:1) * sd2 + mean2
          }
        }
      } else {
        values[[2]] <- unique(data[, effects[2]])
      }
      data <- do.call(expand.grid, values)
    }
  } else if (pred_types == "numeric") {
    # just a single numeric predictor
    data <- structure(data.frame(values), names = effects)
  }
  # no need to have the same value combination more than once
  data <- unique(data)
  data <- data[do.call(order, as.list(data)), , drop = FALSE]
  data <- replicate(nrow(conditions), data, simplify = FALSE)
  marg_vars <- setdiff(names(conditions), effects)
  for (j in seq_len(nrow(conditions))) {
    data[[j]][, marg_vars] <- conditions[j, marg_vars]
    data[[j]]$cond__ <- rownames(conditions)[j]
  }
  data <- do.call(rbind, data)
  data$cond__ <- factor(data$cond__, rownames(conditions))
  structure(data, effects = effects, types = pred_types, mono = mono)
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

get_estimate <- function(coef, samples, margin = 2, to.array = FALSE, ...) {
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
  if (to.array && length(dim(x)) == 2) {
    x <- array(x, dim = c(dim(x), 1), 
               dimnames = list(NULL, NULL, coef))
  }
  x 
}

get_summary <- function(samples, probs = c(0.025, 0.975),
                        robust = FALSE) {
  # summarizes parameter samples based on mean, sd, and quantiles
  # Args: 
  #   samples: a matrix or data.frame containing the samples to be summarized. 
  #            rows are samples, columns are parameters
  #   probs: quantiles to be computed
  # Returns:
  #   a N x C matric where N is the number of observations and C 
  #   is equal to \code{length(probs) + 2}.
  if (robust) {
    coefs <- c("median", "mad", "quantile")
  } else {
    coefs <- c("mean", "sd", "quantile")
  }
  if (length(dim(samples)) == 2L) {
    out <- do.call(cbind, lapply(coefs, get_estimate, samples = samples,
                                 probs = probs, na.rm = TRUE))
  } else if (length(dim(samples)) == 3L) {
    out <- abind(lapply(seq_len(dim(samples)[3]), function(i)
      do.call(cbind, lapply(coefs, get_estimate, samples = samples[, , i],
                            probs = probs))), along = 3)
    dimnames(out) <- list(NULL, NULL, paste0("P(Y = ", 1:dim(out)[3], ")")) 
  } else { 
    stop("Dimension of 'samples' must be either 2 or 3.") 
  }
  rownames(out) <- seq_len(nrow(out))
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
    function(n) table(factor(samples[, n], levels = levels))))
  # compute relative frequencies
  out <- out / sum(out[1, ])
  rownames(out) <- seq_len(nrow(out))
  colnames(out) <- paste0("P(Y = ", seq_len(ncol(out)), ")")
  out
}

get_cov_matrix <- function(sd, cor = NULL) {
  # compute covariance and correlation matrices based 
  # on correlation and sd samples
  # Args:
  #   sd: samples of standard deviations
  #   cor: samples of correlations
  # Notes: 
  #   used in VarCorr.brmsfit
  # Returns: 
  #   samples of covariance and correlation matrices
  sd <- as.matrix(sd)
  stopifnot(all(sd >= 0))
  if (!is.null(cor)) {
    cor <- as.matrix(cor)
    stopifnot(
      ncol(cor) == ncol(sd) * (ncol(sd) - 1) / 2, 
      nrow(sd) == nrow(cor), min(cor) >= -1 && max(cor) <= 1
    )
  }
  nsamples <- nrow(sd)
  nranef <- ncol(sd)
  cor_matrix <- array(diag(1, nranef), dim = c(nranef, nranef, nsamples))
  cor_matrix <- aperm(cor_matrix, perm = c(3, 1, 2))
  cov_matrix <- cor_matrix
  for (i in seq_len(nranef)) { 
    cov_matrix[, i, i] <- sd[, i]^2
  }
  if (!is.null(cor)) {
    k <- 0 
    for (i in 2:nranef) {
      for (j in 1:(i-1)) {
        k = k + 1
        cor_matrix[, j, i] <- cor_matrix[, i, j] <- cor[, k]
        cov_matrix[, j, i] <- 
          cov_matrix[, i, j] <- cor[, k] * sd[, i] * sd[, j]
      }
    }
  }
  list(cor = cor_matrix, cov = cov_matrix)
}

get_cov_matrix_ar1 <- function(ar, sigma, nrows, se2 = NULL) {
  # compute the covariance matrix for an AR1 process
  # Args: 
  #   ar: AR1 autocorrelation samples
  #   sigma: standard deviation samples of the AR1 process
  #   se2: optional square of user defined standard errors
  #   nrows: number of rows of the covariance matrix
  # Returns:
  #   An nsamples x nrows x nrows AR1 covariance array (!)
  sigma <- as.matrix(sigma)
  if (!length(se2)) se2 <- 0
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

get_cov_matrix_ma1 <- function(ma, sigma, nrows, se2 = NULL) {
  # compute the covariance matrix for an MA1 process
  # Args: 
  #   ma: MA1 autocorrelation samples
  #   sigma: standard deviation samples of the AR1 process
  #   se2: optional square of user defined standard errors
  #   nrows: number of rows of the covariance matrix
  # Returns:
  #   An nsamples x nrows x nrows MA1 covariance array (!)
  sigma <- as.matrix(sigma)
  if (!length(se2)) se2 <- 0
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

get_cov_matrix_arma1 <- function(ar, ma, sigma, nrows, se2 = NULL) {
  # compute the covariance matrix for an AR1 process
  # Args: 
  #   ar: AR1 autocorrelation sample
  #   ma: MA1 autocorrelation sample
  #   sigma: standard deviation samples of the AR1 process
  #   se2: optional square of user defined standard errors
  #   nrows: number of rows of the covariance matrix
  # Returns:
  #   An nsamples x nrows x nrows ARMA1 covariance array (!)
  sigma <- as.matrix(sigma)
  if (!length(se2)) se2 <- 0
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

get_cov_matrix_ident <- function(sigma, nrows, se2 = 0) {
  # compute a variance matrix without including ARMA parameters
  # only used for ARMA covariance models when incl_autocor = FALSE
  # Args:
  #   sigma: standard deviation samples
  #   se2: square of user defined standard errors (may be 0)
  #   nrows: number of rows of the covariance matrix
  # Returns:
  #   An nsamples x nrows x nrows sigma array
  sigma <- as.matrix(sigma)
  mat <- array(diag(se2, nrows), dim = c(nrows, nrows, nrow(sigma)))
  mat <- aperm(mat, perm = c(3, 1, 2))
  sigma2 <- sigma^2
  for (i in seq_len(nrows)) {
    mat[, i, i] <- mat[, i, i] + sigma2
  }
  mat
}

get_auxpar <- function(x, i = NULL) {
  # get samples of an auxiliary parameter
  # Args:
  #   x: object to extract postarior samples from
  #   i: the current observation number
  #      (used in predict and log_lik)
  if (is.list(x)) {
    # compute samples of a predicted parameter
    family <- x[["f"]]
    x <- get_eta(x, i = i)
    if (!nzchar(family$family)) {
      # apply links for auxiliary parameters only
      # the main family link is applied later on
      x <- ilink(x, family$link)
    }
  } else {
    if (!is.null(i) && is.matrix(x) && ncol(x) > 1L) {
      x <- x[, i, drop = FALSE]
    }
  }
  if (is.null(i) && is.matrix(x) && ncol(x) == 1L) {
    # for compatibility with fitted helper functions
    x <- as.vector(x)
  }
  x
}

get_sigma <- function(x, data, i = NULL, dim = NULL) {
  # get the residual standard devation of linear models
  # Args: 
  #    see get_auxpar
  #    dim: target dimension of output matrices (used in fitted)
  stopifnot(is.atomic(x) || is.list(x))
  out <- get_se(data = data, i = i, dim = dim)
  if (!is.null(x)) {
    out <- sqrt(out^2 + get_auxpar(x, i = i)^2)
  }
  mult_disp(out, data = data, i = i, dim = dim)
}

get_shape <- function(x, data, i = NULL, dim = NULL) {
  # get the shape parameter of gamma, weibull and negbinomial models
  # Args: see get_auxpar
  stopifnot(is.atomic(x) || is.list(x))
  x <- get_auxpar(x, i = i)
  mult_disp(x, data = data, i = i, dim = dim)
}

get_theta <- function(draws, i = NULL, par = c("zi", "hu")) {
  # convenience function to extract zi / hu parameters
  # also works with deprecated models fitted with brms < 1.0.0 
  # which were using multivariate syntax
  # Args:
  #   see get_auxpar
  #   par: parameter to extract; either 'zi' or 'hu'
  par <- match.arg(par)
  if (!is.null(draws$data$N_trait)) {
    j <- if (!is.null(i)) i else seq_len(draws$data$N_trait)
    theta <- ilink(get_eta(draws$mu, j + draws$data$N_trait), "logit")
  } else {
    theta <- get_auxpar(draws[[par]], i = i)
  }
  theta
}

get_disc <- function(draws, i = NULL, ncat = NULL) {
  # convenience function to extract discrimination parameters
  # Args:
  #   see get_auxpar 
  #   ncat: number of response categories
  if (!is.null(draws[["disc"]])) {
    disc <- get_auxpar(draws[["disc"]], i)
    if (!is.null(dim(disc))) {
      stopifnot(is.numeric(ncat))
      disc <- array(disc, dim = c(dim(disc), ncat - 1))
    }
  } else {
    disc <- 1
  }
  disc
}

get_se <- function(data, i = NULL, dim = NULL) {
  # extract user-defined standard errors
  # Args: see get_auxpar
  se <- data[["se"]]
  if (!is.null(se)) {
    if (!is.null(i)) {
      se <- se[i]
    } else {
      stopifnot(!is.null(dim))
      se <- matrix(se, nrow = dim[1], ncol = dim[2], byrow = TRUE)
    }
  } else {
    se <- 0
  }
  se
}

mult_disp <- function(x, data, i = NULL, dim = NULL) {
  # multiply existing samples by 'disp' data
  # Args: see get_auxpar
  if (!is.null(data$disp)) {
    if (!is.null(i)) {
      x <- x * data$disp[i]
    } else {
      # results in a nsamples x Nobs matrix
      if (is.matrix(x)) {
        stopifnot(!is.null(dim))
        disp <- matrix(disp, nrow = dim[1], ncol = dim[2], byrow = TRUE)
        x <- x * disp
      } else {
        disp <- matrix(data$disp, nrow = 1) 
        if (length(x) == 1L) {
          x <- x * disp
        } else {
          x <- x %*% disp 
        }
      }
    }
  }
  x
}

prepare_family <- function(x) {
  # prepare for calling family specific log_lik / predict functions
  family <- family(x)
  nresp <- length(parse_bf(x$formula, family = family)$response)
  if (is_old_lognormal(family, nresp = nresp, version = x$version$brms)) {
    family <- lognormal()
  } else if (is_linear(family) && nresp > 1L) {
    family$family <- paste0(family$family, "_mv")
  } else if (use_cov(x$autocor) && sum(x$autocor$p, x$autocor$q) > 0) {
    family$family <- paste0(family$family, "_cov")
  } else if (is.cor_fixed(x$autocor)) {
    family$family <- paste0(family$family, "_fixed")
  }
  family
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
    eta <- eta[, old_order, drop = FALSE]  
    colnames(eta) <- NULL
  }
  eta
}

fixef_pars <- function() {
  # regex to extract population-level coefficients
  "^b(|cs|mo|me|m)_"
}

default_plot_pars <- function() {
  # list all parameter classes to be included in plots by default
  c(fixef_pars(), "^sd_", "^cor_", "^sigma_", "^rescor_", 
    paste0("^", auxpars(), "[[:digit:]]*$"), "^delta$",
    "^theta", "^ar", "^ma", "^arr", "^sigmaLL", "^sds_")
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

compute_ic <- function(x, ic = c("waic", "loo", "psislw"), 
                       loo_args = list(), ...) {
  # compute WAIC and LOO using the 'loo' package
  # Args:
  #   x: an object of class brmsfit
  #   ic: the information criterion to be computed
  #   loo_args: passed to functions of the loo package
  #   ...: passed to log_lik.brmsfit
  # Returns:
  #   output of the loo functions with amended class attribute
  stopifnot(is.list(loo_args))
  ic <- match.arg(ic)
  dots <- list(...)
  contains_samples(x)
  loo_args$x <- do.call(log_lik, c(list(x), dots))
  pointwise <- is.function(loo_args$x)
  if (pointwise) {
    loo_args$args <- attr(loo_args$x, "args")
    attr(loo_args$x, "args") <- NULL
  }
  if (ic == "psislw") {
    if (pointwise) {
      loo_args[["llfun"]] <- loo_args[["x"]]
      loo_args[["llargs"]] <- loo_args[["args"]]
      loo_args[["x"]] <- loo_args[["args"]] <- NULL
    } else {
      loo_args[["lw"]] <- -loo_args[["x"]]
      loo_args[["x"]] <- NULL
    }
  }
  IC <- do.call(eval(parse(text = paste0("loo::", ic))), loo_args)
  class(IC) <- c("ic", "loo")
  IC
}

#' Compare Information Criteria of Different Models
#'
#' Compare information criteria of different models fitted
#' with \code{\link[brms:WAIC]{WAIC}} or \code{\link[brms:LOO]{LOO}}.
#' 
#' @param ... At least two objects returned by 
#'   \code{\link[brms:WAIC]{WAIC}} or \code{\link[brms:LOO]{LOO}}.
#' @param x A list of at least two objects returned by
#'   \code{\link[brms:WAIC]{WAIC}} or \code{\link[brms:LOO]{LOO}}.
#'   This argument can be used as an alternative to specifying the 
#'   models in \code{...}.
#'   
#' @return An object of class \code{iclist}.
#' 
#' @details For more details see \code{\link[loo:compare]{compare}}.
#' 
#' @seealso 
#'   \code{\link[brms:WAIC]{WAIC}}, 
#'   \code{\link[brms:LOO]{LOO}},
#'   \code{\link[loo:compare]{compare}}
#'   
#' @examples 
#' \dontrun{
#' # model with population-level effects only
#' fit1 <- brm(rating ~ treat + period + carry,
#'             data = inhaler, family = "gaussian")
#' w1 <- WAIC(fit1)
#' 
#' # model with an additional varying intercept for subjects
#' fit2 <- brm(rating ~ treat + period + carry + (1|subject),
#'             data = inhaler, family = "gaussian")
#' w2 <- WAIC(fit2)
#' 
#' # compare both models
#' compare_ic(w1, w2)
#' }
#' 
#' @export
compare_ic <- function(..., x = NULL) {
  # compare information criteria of different models
  # Args:
  #   x: A list containing loo objects
  #   ic: the information criterion to be computed
  # Returns:
  #   A matrix with differences in the ICs 
  #   as well as corresponding standard errors
  if (!(is.null(x) || is.list(x))) {
    stop2("Argument 'x' should be a list.")
  }
  x$ic_diffs__ <- NULL
  x <- c(list(...), x)
  if (!all(sapply(x, inherits, "ic"))) {
    stop2("All inputs should have class 'ic'.")
  }
  if (length(x) < 2L) {
    stop2("Expecting at least two objects.")
  }
  ics <- unname(sapply(x, function(y) names(y)[3]))
  if (!all(sapply(ics, identical, ics[1]))) {
    stop2("All inputs should be from the the same criterion.")
  }
  names(x) <- ulapply(x, "[[", "model_name")
  n_models <- length(x)
  ic_diffs <- matrix(0, nrow = n_models * (n_models - 1) / 2, ncol = 2)
  rnames <- rep("", nrow(ic_diffs))
  # pairwise comparision to get differences in ICs and their SEs
  n <- 1
  for (i in seq_len(n_models - 1)) {
    for (j in (i + 1):n_models) {
      temp <- loo::compare(x[[j]], x[[i]])
      ic_diffs[n, ] <- c(-2 * temp[["elpd_diff"]], 2 * temp[["se"]]) 
      rnames[n] <- paste(names(x)[i], "-", names(x)[j])
      n <- n + 1
    }
  }
  rownames(ic_diffs) <- rnames
  colnames(ic_diffs) <- c(toupper(ics[1]), "SE")
  x$ic_diffs__ <- ic_diffs
  class(x) <- c("iclist", "list")
  x
}

loo_weights <- function(object, lw = NULL, log = FALSE, 
                        loo_args = list(), ...) {
  # compute loo weights for use in loo_predict and related methods
  # Args:
  #   object: a brmsfit object
  #   lw: precomputed log weights matrix
  #   log: return log weights?
  #   loo_args: further arguments passed to functions of loo
  #   ...: further arguments passed to compute_ic
  # Returns:
  #   an S x N matrix
  if (!is.null(lw)) {
    stopifnot(is.matrix(lw))
  } else {
    message("Running PSIS to compute weights")
    psis <- compute_ic(object, ic = "psislw", loo_args = loo_args, ...)
    lw <- psis[["lw_smooth"]]
  }
  if (!log) {
    lw <- exp(lw) 
  } 
  lw
}

set_pointwise <- function(x, pointwise = NULL, newdata = NULL, 
                          subset = NULL, thres = 1e+08) {
  # set the pointwise argument based on the model size
  # Args:
  #   x: a brmsfit object
  #   newdata: optional data.frame containing new data
  #   subset: a vector to indicate a subset of the posterior samples
  #   thres: threshold above which pointwise is set to TRUE
  # Returns:
  #   TRUE or FALSE
  if (!is.null(pointwise)) {
    pointwise <- as.logical(pointwise)
    if (length(pointwise) != 1L || anyNA(pointwise)) {
      stop2("Argument 'pointwise' must be either TRUE or FALSE.")
    }
  } else {
    nsamples <- nsamples(x, subset = subset)
    if (is.data.frame(newdata)) {
      nobs <- nrow(newdata)
    } else {
      nobs <- nobs(x)
    }
    pointwise <- nsamples * nobs > thres
    if (pointwise) {
      message(
        "Switching to pointwise evaluation to reduce ",  
        "RAM requirements.\nThis will likely increase ",
        "computation time."
      )
    }
  }
  pointwise
}

match_response <- function(models) {
  # compare the response parts of multiple brmsfit objects
  # Args:
  #   models: A list of brmsfit objects
  # Returns:
  #   TRUE if the response parts of all models match and FALSE else
  if (length(models) <= 1L) {
    out <- TRUE  
  } else {
    add_funs <- lsp("brms", what = "exports", pattern = "^resp_")
    match_vars <- c("Y", sub("^resp_", "", add_funs))
    .match_fun <- function(x, y) {
      # checks if all relevant parts of the response are the same 
      # Args:
      #   x, y: named lists as returned by standata
      all(ulapply(match_vars, function(v) {
        a <- if (is.null(attr(x, "old_order"))) as.vector(x[[v]])
             else as.vector(x[[v]])[attr(x, "old_order")]
        b <- if (is.null(attr(y, "old_order"))) as.vector(y[[v]])
             else as.vector(y[[v]])[attr(y, "old_order")]
        is_equal(a, b)
      }))
    }
    sdatas <- lapply(models, standata, control = list(save_order = TRUE))
    matches <- ulapply(sdatas[-1], .match_fun, y = sdatas[[1]]) 
    if (all(matches)) {
      out <- TRUE
    } else {
      out <- FALSE
      warning2("Model comparisons are likely invalid as the response ", 
               "parts of at least two models do not match.")
    }
  }
  invisible(out)
}

find_names <- function(x) {
  # find all valid object names in a string 
  # Args:
  #   x: a character string
  # Notes:
  #   does not use the R parser itself to allow for 
  #   square brackets and kommas at the end of names
  #   currently only used in hypothesis.brmsfit
  # Returns:
  #   all valid variable names within the string
  if (!is.character(x) || length(x) > 1) {
    stop2("Argument 'x' must be a character string of length one.")
  }
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

make_point_frame <- function(mf, effects, conditions, groups, family) {
  # helper function for marginal_effects.brmsfit
  # allowing add data points to the marginal plots
  # Args:
  #   mf: the original model.frame
  #   effects: see argument 'effects' of marginal_effects
  #   conditions: see argument 'conditions' of marginal_effects
  #   groups: names of the grouping factors
  #   family: the model family
  # Returns:
  #   a data.frame containing the data points to be plotted
  points <- mf[, effects, drop = FALSE]
  points$resp__ <- model.response(mf)
  # get required variables i.e. (grouping) factors
  list_mf <- lapply(as.list(mf), function(x)
    if (is.numeric(x)) x else as.factor(x))
  req_vars <- names(mf)[sapply(list_mf, is.factor)]
  if (length(groups)) {
    req_vars <- c(req_vars, unlist(strsplit(groups, ":")))
  }
  req_vars <- unique(setdiff(req_vars, effects))
  req_vars <- intersect(req_vars, names(conditions))
  if (length(req_vars)) {
    # find out which data point is valid for which condition
    mf <- mf[, req_vars, drop = FALSE]
    conditions <- conditions[, req_vars, drop = FALSE]
    points$cond__ <- NA
    points <- replicate(nrow(conditions), points, simplify = FALSE)
    for (i in seq_along(points)) {
      cond <- conditions[i, , drop = FALSE]
      not_na <- which(!(is.na(cond) | cond == "zero__"))
      if (length(not_na)) {
        # do it like base::duplicated
        K <- do.call("paste", c(mf[, not_na, drop = FALSE], sep = "\r")) %in% 
             do.call("paste", c(cond[, not_na, drop = FALSE], sep = "\r"))
      } else {
        K <- seq_len(nrow(mf))
      }
      points[[i]]$cond__[K] <- rownames(conditions)[i] 
    }
    points <- do.call(rbind, points)
    # cond__ allows to assign points to conditions
    points$cond__ <- factor(points$cond__, rownames(conditions))
  }
  if (!is.numeric(points$resp__)) {
    points$resp__ <- as.numeric(as.factor(points$resp__))
    if (is_binary(family)) {
      points$resp__ <- points$resp__ - 1
    }
  }
  na.omit(points)
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
