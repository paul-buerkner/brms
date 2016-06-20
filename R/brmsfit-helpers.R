array2list <- function(x) {
  # convert array to list of elements with reduced dimension
  # Args: 
  #   x: an arrary of dimension d
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

Nsamples <- function(x, subset = NULL) {
  # compute the number of posterior samples
  # Args:
  #   x: a brmsfit object
  #   subset: a vector defining a subset of samples to be considered
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) {
    return(0)
  }
  ntsamples <- (x$fit@sim$iter - x$fit@sim$warmup) /
                 x$fit@sim$thin * x$fit@sim$chains
  if (length(subset)) {
    out <- length(subset)
    if (out > ntsamples || max(subset) > ntsamples) {
      stop("invalid 'subset' argument", call. = FALSE)
    }
  } else {
    out <- ntsamples
  }
  out
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
  if (link == "identity") x
  else if (link == "log") log(x)
  else if (link == "inverse") 1/x
  else if (link == "sqrt") sqrt(x)
  else if (link == "1/mu^2") 1 / x^2
  else if (link == "logit") logit(x)
  else if (link == "probit") qnorm(x)
  else if (link == "probit_approx") qnorm(x)
  else if (link == "cloglog") cloglog(x)
  else if (link == "cauchit") qcauchy(x)
  else stop(paste("Link", link, "not supported"))
}

ilink <- function(x, link) {
  # apply the inverse link function on x
  # Args:
  #   x: An arrary of arbitrary dimension
  #   link: a character string defining the link
  # Returns:
  #   an array of dimension dim(x) on which the inverse link function was applied
  if (link == "identity") x
  else if (link == "log") exp(x)
  else if (link == "inverse") 1/x
  else if (link == "sqrt") x^2
  else if (link == "1/mu^2") 1 / sqrt(x)
  else if (link == "logit") inv_logit(x)
  else if (link == "probit") pnorm(x)
  else if (link == "probit_approx") inv_logit(0.07056*x^3 + 1.5976*x)
  else if (link == "cloglog") inv_cloglog(x)
  else if (link == "cauchit") pcauchy(x)
  else stop(paste("Link", link, "not supported"))
}

get_cornames <- function(names, type = "cor", brackets = TRUE) {
  # get correlation names as combinations of variable names
  # Args:
  #   names: the variable names 
  #   type: of the correlation to be put in front of the returned strings
  #   brackets: should the correlation names contain brackets 
  #            or underscores as seperators
  # Returns: 
  #   correlation names based on the variable names passed to the names argument
  cornames <- NULL
  if (length(names) > 1) {
    for (i in 2:length(names)) {
      for (j in 1:(i-1)) {
        if (brackets) {
          cornames <- c(cornames, paste0(type,"(",names[j],",",names[i],")"))
        } else {
          cornames <- c(cornames, paste0(type,"_",names[j],"_",names[i]))
        }
      }
    }
  }
  cornames
}

get_nlpar <- function(x, suffix = "") {
  # extract name of a non-linear parameter
  nlpar <- attr(x, "nlpar")
  if (!is.null(nlpar)) paste0(nlpar, suffix) else ""
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
  if (length(dim(samples)) == 2) {
    out <- do.call(cbind, lapply(coefs, get_estimate, samples = samples,
                                 probs = probs, na.rm = TRUE))
  } else if (length(dim(samples)) == 3) {
    out <- abind(lapply(1:dim(samples)[3], function(i)
      do.call(cbind, lapply(coefs, get_estimate, samples = samples[, , i],
                            probs = probs))), along = 3)
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
  # Args:
  #   samples: a S x N matrix
  #   levels: all possible values in \code{samples}
  # Returns:
  #    a N x \code{levels} matrix containing absolute frequencies in each column seperately
  if (!is.matrix(samples)) 
    stop("samples must be a matrix")
  out <- do.call(rbind, lapply(1:ncol(samples), function(n) 
    table(factor(samples[, n], levels = levels))))
  rownames(out) <- 1:nrow(out)
  colnames(out) <- paste0("N(Y = ", 1:ncol(out), ")")
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
  if (any(sd < 0)) 
    stop("standard deviations must be non negative")
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

get_cov_matrix_ar1 <- function(ar, sigma, nrows, se2 = 0) {
  # compute the covariance matrix for an AR1 process
  # Args: 
  #   ar: AR1 autocorrelation samples
  #   sigma: standard deviation samples of the AR1 process
  #   se2: square of user defined standard errors (may be 0)
  #   nrows: number of rows of the covariance matrix
  # Returns:
  #   An nsamples x nrows x nrows AR1 covariance array (!)
  mat <- aperm(array(diag(se2, nrows), dim = c(nrows, nrows, nrow(ar))),
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

get_cov_matrix_ma1 <- function(ma, sigma, nrows, se2 = 0) {
  # compute the covariance matrix for an MA1 process
  # Args: 
  #   ma: MA1 autocorrelation samples
  #   sigma: standard deviation samples of the AR1 process
  #   se2: square of user defined standard errors (may be 0)
  #   nrows: number of rows of the covariance matrix
  # Returns:
  #   An nsamples x nrows x nrows MA1 covariance array (!)
  mat <- aperm(array(diag(se2, nrows), dim = c(nrows, nrows, nrow(ma))),
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

get_cov_matrix_arma1 <- function(ar, ma, sigma, nrows, se2 = 0) {
  # compute the covariance matrix for an AR1 process
  # Args: 
  #   ar: AR1 autocorrelation sample
  #   ma: MA1 autocorrelation sample
  #   sigma: standard deviation samples of the AR1 process
  #   se2: square of user defined standard errors (may be 0)
  #   nrows: number of rows of the covariance matrix
  # Returns:
  #   An nsamples x nrows x nrows ARMA1 covariance array (!)
  mat <- aperm(array(diag(se2, nrows), dim = c(nrows, nrows, nrow(ar))),
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

get_cov_matrix_ident <- function(sigma, nrows, se2 = 0) {
  # compute a variance matrix without including ARMA parameters
  # only used for ARMA covariance models when incl_autor = FALSE
  # Args:
  #   sigma: standard deviation samples of the AR1 process
  #   se2: square of user defined standard errors (may be 0)
  #   nrows: number of rows of the covariance matrix
  # Returns:
  #   An nsamples x nrows x nrows sigma array
  mat <- aperm(array(diag(se2, nrows), dim = c(nrows, nrows, nrow(sigma))),
               perm = c(3, 1, 2))
  sigma2 <- sigma^2
  for (i in 1:nrows) {
    mat[, i, i] <- mat[, i, i] + sigma2
  }
  mat
}

get_sigma <- function(x, data, i, method = c("fitted", "predict", "logLik")) {
  # get the residual standard devation of linear models
  # Args:
  #   x: a brmsfit object or posterior samples of sigma (can be NULL)
  #   data: data initially passed to Stan
  #   method: S3 method from which get_sigma is called
  #   i: meaning depends on the method argument:
  #      for predict and logLik this is the current observation number
  #      for fitted this is the number of samples
  method <- match.arg(method)
  if (is(x, "brmsfit")) {
    sigma <- posterior_samples(x, pars = "^sigma_")$sigma
  } else {
    sigma <- x
  }
  sigma <- as.vector(sigma)
  if (is.null(sigma)) {
    # user defined standard errors were applied
    sigma <- data$se
    if (is.null(sigma)) {
      # for backwards compatibility with brms <= 0.5.0
      sigma <- data$sigma
    }
    if (is.null(sigma)) {
      stop("no residual standard deviation(s) found")
    }
    if (method %in% c("predict", "logLik")) {
      sigma <- sigma[i]
    } else {
      sigma <- matrix(rep(sigma, i), ncol = data$N, byrow = TRUE)
    }
  } else if (!is.null(data$disp)) {
    if (method %in% c("predict", "logLik")) {
      sigma <- sigma * data$disp[i]
    } else {
      # results in a Nsamples x Nobs matrix
      sigma <- sigma %*% matrix(data$disp, nrow = 1)
    }
  }
  sigma
}

get_shape <- function(x, data, i = NULL,
                      method = c("fitted", "predict", "logLik")) {
  # get the shape parameter of gamma, weibull and negbinomial models
  # Args:
  #   x: a brmsfit object or posterior samples of shape (can be NULL)
  #   data: data initially passed to Stan
  #   method: S3 method from which get_sigma is called
  #   i: only used for "predict" and "logLik": 
  #      the current observation number
  method <- match.arg(method)
  if (is(x, "brmsfit")) {
    shape <- posterior_samples(x, pars = "^shape$")$shape
  } else {
    shape <- x
  }
  shape <- as.vector(shape)
  if (!is.null(data$disp)) {
    if (method %in% c("predict", "logLik")) {
      shape <- shape * data$disp[i]
    } else {
      # results in a Nsamples x Nobs matrix
      shape <- shape %*% matrix(data$disp, nrow = 1)
    }
  }
  shape
}

prepare_family <- function(x) {
  # prepare for calling family specific loglik / predict functions
  family <- family(x)
  nresp <- length(extract_effects(x$formula, family = family,
                                  nonlinear = x$nonlinear)$response)
  if (is.old_lognormal(family, nresp = nresp, version = x$version)) {
    family <- lognormal()
  } else if (is.linear(family) && nresp > 1L) {
    family$family <- paste0(family$family, "_multi")
  } else if (use_cov(x$autocor) && sum(x$autocor$p, x$autocor$q) > 0) {
    family$family <- paste0(family$family, "_cov")
  } else if (is(x$autocor, "cor_fixed")) {
    family$family <- paste0(family$family, "_fixed")
  }
  family
}

default_plot_pars <- function() {
  # list all parameter classes to be included in plots by default
  c("^b_", "^bm_", "^sd_", "^cor_", "^sigma", "^rescor", 
    "^nu$", "^shape$", "^delta$", "^phi$", "^ar", "^ma", 
    "^arr", "^simplex_", "^sds_")
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

compute_ic <- function(x, ic = c("waic", "loo"), ll_args = list(), ...) {
  # compute WAIC and LOO using the 'loo' package
  # Args:
  #   x: an object of class brmsfit
  #   ic: the information criterion to be computed
  #   ll_args: a list of additional arguments passed to logLik
  #   ...: passed to the loo package
  # Returns:
  #   output of the loo package with amended class attribute
  ic <- match.arg(ic)
  if (!is(x, "brmsfit")) 
    stop(paste("Cannot compute information criteria for", 
               "an object of class", class(x)), call. = FALSE)
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  args <- list(x = do.call(logLik, c(list(x), ll_args)))
  if (ll_args$pointwise) {
    args$args$draws <- attr(args$x, "draws")
    args$args$data <- data.frame()
    args$args$N <- attr(args$x, "N")
    args$args$S <- Nsamples(x, subset = ll_args$subset)
    attr(args$x, "draws") <- NULL
  }
  if (ic == "loo") {
    args <- c(args, ...)
  }
  IC <- do.call(eval(parse(text = paste0("loo::", ic))), args)
  class(IC) <- c("ic", "loo")
  return(IC)
}

compare_ic <- function(x, ic = c("waic", "loo")) {
  # compare information criteria of different models
  # Args:
  #   x: A list containing loo objects
  #   ic: the information criterion to be computed
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
    if ("weights" %in% colnames(all_compare)) {
      weights <- unname(all_compare[do.call(get_input_names, x), "weights"])
    } else {
      # weights have been temporarily removed in loo 0.1.5
      weights <- rep(NA, n_models)
    }
  }
  nlist(ic_diffs, weights)
}

set_pointwise <- function(x, newdata = NULL, subset = NULL, thres = 1e+07) {
  # set the pointwise argument based on the model size
  # Args:
  #   x: a brmsfit object
  #   newdata: optional data.frame containing new data
  #   subset: a vector to indicate a subset of the posterior samples
  #   thres: threshold above which pointwise is set to TRUE
  # Returns:
  #   TRUE or FALSE
  nsamples <- Nsamples(x, subset = subset)
  if (is.data.frame(newdata)) {
    nobs <- nrow(newdata)
  } else {
    nobs <- nobs(x)
  }
  pointwise <- nsamples * nobs > thres
  if (pointwise) {
    message(paste0("Switching to pointwise evaluation to reduce ",  
                   "RAM requirements.\nThis will likely increase ",
                   "computation time."))
  }
  pointwise
}

match_response <- function(models) {
  # compare the response parts of multiple brmsfit objects
  # Args:
  #   models: A list of brmsfit objects
  # Returns:
  #   TRUE if the response parts of all models match and FALSE else
  if (length(models) <= 1) return(TRUE)
  .match_fun <- function(x, y) {
    # checks if all relevant parts of the response are the same 
    # Args:
    #   x, y: named lists as returned by standata
    to_match <- c("Y", "se", "weights", "cens", "trunc")
    all(ulapply(to_match, function(v) {
      a <- if (is.null(attr(x, "old_order"))) as.vector(x[[v]])
           else as.vector(x[[v]])[attr(x, "old_order")]
      b <- if (is.null(attr(y, "old_order"))) as.vector(y[[v]])
           else as.vector(y[[v]])[attr(y, "old_order")]
      is_equal(a, b)
    }))
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
  # Args:
  #   x: a character string
  # Notes:
  #   currently used in hypothesis.brmsfit
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

evidence_ratio <- function(x, cut = 0, wsign = c("equal", "less", "greater"), 
                           prior_samples = NULL, pow = 12, ...) {
  # compute the evidence ratio between two disjunct hypotheses
  # Args:
  #   x: posterior samples 
  #   cut: the cut point between the two hypotheses
  #   wsign: direction of the hypothesis
  #   prior_samples: optional prior samples for undirected hypothesis
  #   pow influence: the accuracy of the density
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
  points <- mf[, effects[1], drop = FALSE]
  points$.RESP <- model.response(mf)
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
    points$MargCond <- NA
    points <- replicate(nrow(conditions), points, simplify = FALSE)
    for (i in seq_along(points)) {
      cond <- conditions[i, , drop = FALSE]
      not_na <- c(!is.na(cond))
      # do it like base::duplicated
      if (any(not_na)) {
        K <- do.call("paste", c(mf[, not_na, drop = FALSE], sep = "\r")) %in% 
             do.call("paste", c(cond[, not_na, drop = FALSE], sep = "\r"))
      } else {
        K <- 1:nrow(mf)
      }
      points[[i]]$MargCond[K] <- rownames(conditions)[i] 
    }
    points <- do.call(rbind, points)
    # MargCond allows to assign points to conditions
    points$MargCond <- factor(points$MargCond, rownames(conditions))
  }
  if (!is.numeric(points$.RESP)) {
    points$.RESP <- as.numeric(as.factor(points$.RESP))
    if (is.binary(family)) {
      points$.RESP <- points$.RESP - 1
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