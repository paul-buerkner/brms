# R helper functions for Gaussian Processes

#' Set up Gaussian process terms in \pkg{brms}
#'
#' Set up a Gaussian process (GP) term in \pkg{brms}. The function does not
#' evaluate its arguments -- it exists purely to help set up a model with
#' GP terms.
#'
#' @param ... One or more predictors for the GP.
#' @param by A numeric or factor variable of the same length as
#'   each predictor. In the numeric vector case, the elements multiply
#'   the values returned by the GP. In the factor variable
#'   case, a separate GP is fitted for each factor level.
#' @param k Optional number of basis functions for computing approximate
#'   GPs. If \code{NA} (the default), exact GPs are computed.
#' @param cov Name of the covariance kernel. By default,
#'   the exponentiated-quadratic kernel \code{"exp_quad"} is used.
#' @param iso A flag to indicate whether an isotropic (\code{TRUE}; the
#'   default) or a non-isotropic GP should be used.
#'   In the former case, the same amount of smoothing is applied to all
#'   predictors. In the latter case, predictors may have different smoothing.
#'   Ignored if only a single predictor is supplied.
#' @param gr Logical; Indicates if auto-grouping should be used (defaults
#'   to \code{TRUE}). If enabled, observations sharing the same
#'   predictor values will be represented by the same latent variable
#'   in the GP. This will improve sampling efficiency
#'   drastically if the number of unique predictor combinations is small
#'   relative to the number of observations.
#' @param cmc Logical; Only relevant if \code{by} is a factor. If \code{TRUE}
#'   (the default), cell-mean coding is used for the \code{by}-factor, that is
#'   one GP per level is estimated. If \code{FALSE}, contrast GPs are estimated
#'   according to the contrasts set for the \code{by}-factor.
#' @param scale Logical; If \code{TRUE} (the default), predictors are
#'   scaled so that the maximum Euclidean distance between two points
#'   is 1. This often improves sampling speed and convergence.
#'   Scaling also affects the estimated length-scale parameters
#'   in that they resemble those of scaled predictors (not of the original
#'   predictors) if \code{scale} is \code{TRUE}.
#' @param c Numeric value only used in approximate GPs. Defines the
#'   multiplicative constant of the predictors' range over which
#'   predictions should be computed. A good default could be \code{c = 5/4}
#'   but we are still working on providing better recommendations.
#'
#' @details A GP is a stochastic process, which
#'  describes the relation between one or more predictors
#'  \eqn{x = (x_1, ..., x_d)} and a response \eqn{f(x)}, where
#'  \eqn{d} is the number of predictors. A GP is the
#'  generalization of the multivariate normal distribution
#'  to an infinite number of dimensions. Thus, it can be
#'  interpreted as a prior over functions. The values of \eqn{f( )}
#'  at any finite set of locations are jointly multivariate
#'  normal, with a covariance matrix defined by the covariance
#'  kernel \eqn{k_p(x_i, x_j)}, where \eqn{p} is the vector of parameters
#'  of the GP:
#'  \deqn{(f(x_1), \ldots f(x_n) \sim MVN(0, (k_p(x_i, x_j))_{i,j=1}^n) .}
#'  The smoothness and general behavior of the function \eqn{f}
#'  depends only on the choice of covariance kernel.
#'  For a more detailed introduction to Gaussian processes,
#'  see \url{https://en.wikipedia.org/wiki/Gaussian_process}.
#'
#'  Below, we describe the currently supported covariance kernels:
#'  \itemize{
#'    \item "exp_quad": The exponentiated quadratic kernel is defined as
#'    \eqn{k(x_i, x_j) = sdgp^2 \exp(- || x_i - x_j ||^2 / (2 lscale^2))},
#'    where \eqn{|| . ||} is the Euclidean norm, \eqn{sdgp} is a
#'    standard deviation parameter, and \eqn{lscale} is characteristic
#'    length-scale parameter. The latter practically measures how close two
#'    points \eqn{x_i} and \eqn{x_j} have to be to influence each other
#'    substantially.
#'  }
#'
#'  In the current implementation, \code{"exp_quad"} is the only supported
#'  covariance kernel. More options will follow in the future.
#'
#' @return An object of class \code{'gp_term'}, which is a list
#'   of arguments to be interpreted by the formula
#'   parsing functions of \pkg{brms}.
#'
#' @examples
#' \dontrun{
#' # simulate data using the mgcv package
#' dat <- mgcv::gamSim(1, n = 30, scale = 2)
#'
#' # fit a simple GP model
#' fit1 <- brm(y ~ gp(x2), dat, chains = 2)
#' summary(fit1)
#' me1 <- conditional_effects(fit1, ndraws = 200, spaghetti = TRUE)
#' plot(me1, ask = FALSE, points = TRUE)
#'
#' # fit a more complicated GP model
#' fit2 <- brm(y ~ gp(x0) + x1 + gp(x2) + x3, dat, chains = 2)
#' summary(fit2)
#' me2 <- conditional_effects(fit2, ndraws = 200, spaghetti = TRUE)
#' plot(me2, ask = FALSE, points = TRUE)
#'
#' # fit a multivariate GP model
#' fit3 <- brm(y ~ gp(x1, x2), dat, chains = 2)
#' summary(fit3)
#' me3 <- conditional_effects(fit3, ndraws = 200, spaghetti = TRUE)
#' plot(me3, ask = FALSE, points = TRUE)
#'
#' # compare model fit
#' loo(fit1, fit2, fit3)
#'
#' # simulate data with a factor covariate
#' dat2 <- mgcv::gamSim(4, n = 90, scale = 2)
#'
#' # fit separate gaussian processes for different levels of 'fac'
#' fit4 <- brm(y ~ gp(x2, by = fac), dat2, chains = 2)
#' summary(fit4)
#' plot(conditional_effects(fit4), points = TRUE)
#' }
#'
#' @seealso \code{\link{brmsformula}}
#' @export
gp <- function(..., by = NA, k = NA, cov = "exp_quad", iso = TRUE,
               gr = TRUE, cmc = TRUE, scale = TRUE, c = 5/4) {
  cov_choices <- c("exp_quad", "matern52", "matern32", "exponential")
  cov <- match.arg(cov, choices = cov_choices)
  call <- match.call()
  label <- deparse0(call)
  vars <- as.list(substitute(list(...)))[-1]
  by <- deparse0(substitute(by))
  cmc <- as_one_logical(cmc)
  if (is.null(call[["gr"]]) && require_old_default("2.12.8")) {
    # the default of 'gr' has changed in version 2.12.8
    gr <- FALSE
  } else {
    gr <- as_one_logical(gr)
  }
  if (length(vars) > 1L) {
    iso <- as_one_logical(iso)
  } else {
    iso <- TRUE
  }
  if (!isNA(k)) {
    supported_hsgp_covs <- c("exp_quad")
    if (!cov %in% supported_hsgp_covs) {
      stop2("HSGPs with covariance kernel '", cov, "' are not yet supported.")
    }
    k <- as.integer(as_one_numeric(k))
    if (k < 1L) {
      stop2("'k' must be positive.")
    }
    c <- as.numeric(c)
    if (length(c) == 1L) {
      c <- rep(c, length(vars))
    }
    if (length(c) != length(vars)) {
      stop2("'c' must be of the same length as the number of covariates.")
    }
    if (any(c <= 0)) {
      stop2("'c' must be positive.")
    }
  } else {
    c <- NA
  }
  scale <- as_one_logical(scale)
  term <- ulapply(vars, deparse0, backtick = TRUE, width.cutoff = 500L)
  out <- nlist(term, label, by, cov, k, iso, gr, cmc, scale, c)
  structure(out, class = "gp_term")
}

# get labels of gaussian process terms
# @param x either a formula or a list containing an element "gp"
# @param data data frame containing the covariates
# @return a data.frame with one row per GP term
frame_gp <- function(x, data) {
  if (is.formula(x)) {
    x <- brmsterms(x, check_response = FALSE)$dpars$mu
  }
  form <- x[["gp"]]
  if (!is.formula(form)) {
    return(empty_data_frame())
  }
  out <- data.frame(
    term = all_terms(form),
    label = NA, cov = NA, k = NA, iso = NA, gr = NA, scale = NA,
    stringsAsFactors = FALSE
  )
  nterms <- nrow(out)
  out$cons <- out$byvars <- out$covars <-
    out$sfx1 <- out$sfx2 <- out$c <- vector("list", nterms)
  for (i in seq_len(nterms)) {
    gp <- eval2(out$term[i])
    out$label[i] <- paste0("gp", rename(collapse(gp$term)))
    out$cov[i] <- gp$cov
    out$k[i] <- gp$k
    out$c[[i]] <- gp$c
    out$iso[i] <- gp$iso
    out$cmc[i] <- gp$cmc
    out$gr[i] <- gp$gr
    out$scale[i] <- gp$scale
    out$covars[[i]] <- gp$term
    if (gp$by != "NA") {
      out$byvars[[i]] <- gp$by
      str_add(out$label[i]) <- rename(gp$by)
      byval <- get(gp$by, data)
      if (is_like_factor(byval)) {
        byval <- unique(as.factor(byval))
        byform <- str2formula(c(ifelse(gp$cmc, "0", "1"), "byval"))
        cons <- rename(colnames(model.matrix(byform)))
        out$cons[[i]] <- rm_wsp(sub("^byval", "", cons))
      }
    }
    # sfx1 is for sdgp and sfx2 is for lscale
    out$sfx1[[i]] <- paste0(out$label[i], out$cons[[i]])
    if (out$iso[i]) {
      out$sfx2[[i]] <- matrix(out$sfx1[[i]])
    } else {
      out$sfx2[[i]] <- outer(out$sfx1[[i]], out$covars[[i]], paste0)
    }
  }
  class(out) <- gpframe_class()
  out
}

gpframe_class <- function() {
  c("gpframe", "data.frame")
}

is.gpframe <- function(x) {
  inherits(x, "gpframe")
}

# covariance matrix of Gaussian processes
# not vectorized over parameter values
cov_gp <- function(x, x_new = NULL, sdgp = 1, lscale = 1, cov = "exp_quad") {
  sdgp <- as.numeric(sdgp)
  lscale <- as.numeric(lscale)
  Dls <- length(lscale)
  cov <- as_one_character(cov)
  kernel <- paste0("kernel_", cov)
  kernel <- get(kernel, asNamespace("brms"))
  if (Dls == 1L) {
    # one dimensional or isotropic GP
    diff_quad <- diff_quad(x = x, x_new = x_new)
    out <- kernel(diff_quad, sdgp = sdgp, lscale = lscale)
  } else {
    # multi-dimensional non-isotropic GP
    diff_quad <- diff_quad(x = x[, 1], x_new = x_new[, 1])
    out <- kernel(diff_abs, sdgp = sdgp, lscale = lscale[1])
    for (d in seq_len(Dls)[-1]) {
      diff_quad <- diff_quad(x = x[, d], x_new = x_new[, d])
      out <- out * kernel(diff_abs, sdgp = sdgp, lscale = lscale[d])
    }
  }
  out
}

# Squared exponential kernel
# @param diff_quad squared difference matrix
kernel_exp_quad <- function(diff_quad, sdgp, lscale) {
  sdgp^2 * exp(-diff_quad / (2 * lscale^2))
}

# Exponential kernel
kernel_exponential <- function(diff_quad, sdgp, lscale) {
  diff_abs <- sqrt(diff_quad)
  sdgp^2 * exp(-diff_abs / lscale)
}

# Matern 3/2 kernel
kernel_matern32 <- function(diff_quad, sdgp, lscale) {
  diff_abs <- sqrt(diff_quad)
  sdgp^2 * (1 + sqrt(3) * diff_abs / lscale) *
    exp(- sqrt(3) * diff_abs / lscale)
}

# Matern 5/2 kernel
kernel_matern52 <- function(diff_quad, sdgp, lscale) {
  diff_abs <- sqrt(diff_quad)
  sdgp^2 * (1 + sqrt(5) * diff_abs / lscale + 5 * diff_quad / (3 * lscale^2)) *
    exp(- sqrt(5) * diff_abs / lscale)
}

# compute squared differences
# @param x vector or matrix
# @param x_new optional vector of matrix with the same ncol as x
# @return an nrow(x) times nrow(x_new) matrix
# @details if matrices are passed results are summed over the columns
diff_quad <- function(x, x_new = NULL) {
  x <- as.matrix(x)
  if (is.null(x_new)) {
    x_new <- x
  } else {
    x_new <- as.matrix(x_new)
  }
  .diff_quad <- function(x1, x2) (x1 - x2)^2
  out <- 0
  for (i in seq_cols(x)) {
    out <- out + outer(x[, i], x_new[, i], .diff_quad)
  }
  out
}

# spectral density function
# vectorized over parameter values
spd_exp_quad <- function(x, sdgp = 1, lscale = 1) {
  NB <- NROW(x)
  D <- NCOL(x)
  Dls <- NCOL(lscale)
  out <- matrix(nrow = length(sdgp), ncol = NB)
  if (Dls == 1L) {
    # one dimensional or isotropic GP
    constant <- sdgp^2 * (sqrt(2 * pi) * lscale)^D
    neg_half_lscale2 <- -0.5 * lscale^2
    for (m in seq_len(NB)) {
      out[, m] <- constant * exp(neg_half_lscale2 * sum(x[m, ]^2))
    }
  } else {
    # multi-dimensional non-isotropic GP
    constant <- sdgp^2 * sqrt(2 * pi)^D * matrixStats::rowProds(lscale)
    neg_half_lscale2 = -0.5 * lscale^2
    for (m in seq_len(NB)) {
      x2 <- data2draws(x[m, ]^2, dim = dim(lscale))
      out[, m] <- constant * exp(rowSums(neg_half_lscale2 * x2))
    }
  }
  out
}

# compute the mth eigen value of an approximate GP
eigen_val_exp_quad <- function(m, L) {
  ((m * pi) / (2 * L))^2
}

# compute the mth eigen function of an approximate GP
eigen_fun_exp_quad <- function(x, m, L) {
  x <- as.matrix(x)
  D <- ncol(x)
  stopifnot(length(m) == D, length(L) == D)
  out <- vector("list", D)
  for (i in seq_cols(x)) {
    out[[i]] <- 1 / sqrt(L[i]) *
      sin((m[i] * pi) / (2 * L[i]) * (x[, i] + L[i]))
  }
  Reduce("*", out)
}

# extended range of input data for which predictions should be made
choose_L <- function(x, c) {
  if (!length(x)) {
    range <- 1
  } else {
    range <- max(1, max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }
  c * range
}

# try to evaluate a GP term and
# return an informative error message if it fails
try_nug <- function(expr, nug) {
  out <- try(expr, silent = TRUE)
  if (is_try_error(out)) {
    stop2("The Gaussian process covariance matrix is not positive ",
          "definite.\nThis occurs for numerical reasons. Setting ",
          "'nug' above ", nug, " may help.")
  }
  out
}
