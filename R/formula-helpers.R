#' Additional Response Information
#' 
#' Provide additional information on the response variable 
#' in \pkg{brms} models, such as censoring, truncation, or
#' known measurement error.
#' 
#' @name addition-terms
#' 
#' @param x A vector; usually a variable defined in the
#'  data. Allowed values depend on the function:
#'  \code{resp_se} and \code{resp_weights} require positive numeric values.
#'  \code{resp_trials} and \code{resp_cat} require positive integers.
#'  \code{resp_dec} requires \code{0} and \code{1}, or alternatively
#'  \code{'lower'} and \code{'upper'}; 
#'  \code{resp_subset} requires \code{0} and \code{1}, or alternatively
#'  \code{FALSE} and \code{TRUE};
#'  \code{resp_cens} requires \code{'left'}, \code{'none'}, \code{'right'},
#'  and \code{'interval'} (or equivalently \code{-1}, \code{0}, \code{1},
#'  and \code{2}) to indicate left, no, right, or interval censoring.
#' @param sigma Logical; Indicates whether the residual standard deviation
#'  parameter \code{sigma} should be included in addition to the known
#'  measurement error. Defaults to \code{FALSE} for backwards compatibility,
#'  but setting it to \code{TRUE} is usually the better choice.
#' @param scale Logical; Indicates whether weights should be scaled
#'  so that the average weight equals one. Defaults to \code{FALSE}.
#' @param y2 A vector specifying the upper bounds in interval censoring.
#' @param lb A numeric vector or single numeric value specifying 
#'   the lower truncation bound.
#' @param ub A numeric vector or single numeric value specifying 
#'   the upper truncation bound.
#' @param sdy Optional known measurement error of the response
#'   treated as standard deviation. If specified, handles
#'   measurement error and (completely) missing values
#'   at the same time using the plausible-values-technique.
#'
#' @return A vector containing additional information on the response
#'   variable in an appropriate format.
#'
#' @details 
#'   These functions are almost solely useful when
#'   called in formulas passed to the \pkg{brms} package.
#'   Within formulas, the \code{resp_} prefix may be omitted.
#'   More information is given in the 'Details' section
#'   of \code{\link{brmsformula}}.
#'   
#' @seealso 
#'   \code{\link{brm}}, 
#'   \code{\link{brmsformula}}   
#'  
#' @examples 
#' \dontrun{
#' ## Random effects meta-analysis
#' nstudies <- 20
#' true_effects <- rnorm(nstudies, 0.5, 0.2)
#' sei <- runif(nstudies, 0.05, 0.3)
#' outcomes <- rnorm(nstudies, true_effects, sei)
#' data1 <- data.frame(outcomes, sei)
#' fit1 <- brm(outcomes | se(sei, sigma = TRUE) ~ 1,
#'             data = data1)
#' summary(fit1)
#' 
#' ## Probit regression using the binomial family
#' n <- sample(1:10, 100, TRUE)  # number of trials
#' success <- rbinom(100, size = n, prob = 0.4)
#' x <- rnorm(100)
#' data2 <- data.frame(n, success, x)
#' fit2 <- brm(success | trials(n) ~ x, data = data2,
#'             family = binomial("probit"))
#' summary(fit2)
#' 
#' ## Survival regression modeling the time between the first 
#' ## and second recurrence of an infection in kidney patients.
#' fit3 <- brm(time | cens(censored) ~ age * sex + disease + (1|patient), 
#'             data = kidney, family = lognormal())
#' summary(fit3)
#' 
#' ## Poisson model with truncated counts  
#' fit4 <- brm(count | trunc(ub = 104) ~ zBase * Trt, 
#'             data = epilepsy, family = poisson())
#' summary(fit4)
#' }
#'   
NULL

#' @rdname addition-terms
#' @export
resp_se <- function(x, sigma = FALSE) {
  # standard errors for meta-analysis
  if (!is.numeric(x)) {
    stop2("Standard errors must be numeric.")
  }
  if (min(x) < 0) {
    stop2("Standard errors must be non-negative.")
  }
  sigma <- as_one_logical(sigma)
  structure(x, sigma = sigma)  
}

# only evaluate the sigma argument of 'resp_se'
resp_se_no_data <- function(x, sigma = FALSE) {
  resp_se(1, sigma = sigma)
}

#' @rdname addition-terms
#' @export
resp_weights <- function(x, scale = FALSE) {
  if (!is.numeric(x)) {
    stop2("Weights must be numeric.")
  }
  if (min(x) < 0) {
    stop2("Weights must be non-negative.")
  }
  scale <- as_one_logical(scale)
  if (scale) {
    x <- x / sum(x) * length(x)
  }
  x
}

#' @rdname addition-terms
#' @export
resp_trials <- function(x) {
  if (!is.numeric(x)) {
    stop2("Number of trials must be numeric.")
  }
  if (any(!is_wholenumber(x) | x < 1)) {
    stop2("Number of trials must be positive integers.")
  }
  x
}

#' @rdname addition-terms
#' @export
resp_cat <- function(x) {
  x <- as_one_numeric(x)
  if (!is_wholenumber(x) || x < 1) {
    stop2("Number of categories must be a positive integer.")
  }
  x
}

#' @rdname addition-terms
#' @export
resp_dec <- function(x) {
  if (is.character(x) || is.factor(x)) {
    if (!all(unique(x) %in% c("lower", "upper"))) {
      stop2("Decisions should be 'lower' or 'upper' ", 
            "when supplied as characters or factors.")
    }
    x <- ifelse(x == "lower", 0, 1)
  } else {
    x <- as.numeric(as.logical(x))
  }
  x
}

#' @rdname addition-terms
#' @export
resp_cens <- function(x, y2 = NULL) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  prepare_cens <- function(x) {
    stopifnot(length(x) == 1L)
    regx <- paste0("^", x)
    if (grepl(regx, "left")) {
      x <- -1
    } else if (grepl(regx, "none") || isFALSE(x)) {
      x <- 0
    } else if (grepl(regx, "right") || isTRUE(x)) {
      x <- 1
    } else if (grepl(regx, "interval")) {
      x <- 2
    }
    x
  }
  cens <- unname(ulapply(x, prepare_cens))
  if (!all(is_wholenumber(cens) & cens %in% -1:2)) {
    stop2(
      "Invalid censoring data. Accepted values are ", 
      "'left', 'none', 'right', and 'interval'\n",
      "(abbreviations are allowed) or -1, 0, 1, and 2.\n",
      "TRUE and FALSE are also accepted ",
      "and refer to 'right' and 'none' respectively."
    )
  }
  if (any(cens %in% 2)) {
    if (!length(y2)) {
      stop2("Argument 'y2' is required for interval censored data.")
    }
    attr(cens, "y2") <- unname(y2)
  }
  cens
}

#' @rdname addition-terms
#' @export
resp_trunc <- function(lb = -Inf, ub = Inf) {
  lb <- as.numeric(lb)
  ub <- as.numeric(ub)
  if (any(lb >= ub)) {
    stop2("Truncation bounds are invalid: lb >= ub")
  }
  nlist(lb, ub)
}

#' @rdname addition-terms
#' @export
resp_mi <- function(sdy = NULL) {
  if (!is.null(sdy) && !is.numeric(sdy)) {
    stop2("Measurement error should be numeric.")
  }
  sdy
}

#' @rdname addition-terms
#' @export
resp_subset <- function(x) {
  as.logical(x)
}

#' Defining smooths in \pkg{brms} formulas
#' 
#' Functions used in definition of smooth terms within a model formulas. 
#' The function does not evaluate a (spline) smooth - it exists purely 
#' to help set up a model using spline based smooths.
#' 
#' @param ... Arguments passed to \code{\link[mgcv:s]{mgcv::s}} or
#'  \code{\link[mgcv:t2]{mgcv::t2}}.
#'  
#' @details The function defined here are just simple wrappers
#'  of the respective functions of the \pkg{mgcv} package.
#'  
#' @seealso \code{\link{brmsformula}},
#'   \code{\link[mgcv:s]{mgcv::s}}, \code{\link[mgcv:t2]{mgcv::t2}}
#'  
#' @examples
#' \dontrun{
#' # simulate some data
#' dat <- mgcv::gamSim(1, n = 200, scale = 2)
#' 
#' # fit univariate smooths for all predictors
#' fit1 <- brm(y ~ s(x0) + s(x1) + s(x2) + s(x3), 
#'             data = dat, chains = 2)
#' summary(fit1)
#' plot(marginal_smooths(fit1), ask = FALSE)
#' 
#' # fit a more complicated smooth model
#' fit2 <- brm(y ~ t2(x0, x1) + s(x2, by = x3), 
#'             data = dat, chains = 2)
#' summary(fit2)
#' plot(marginal_smooths(fit2), ask = FALSE)
#' }
#' 
#' @export
s <- function(...) {
  mgcv::s(...)
}

#' @rdname s
#' @export
t2 <- function(...) {
  mgcv::t2(...)
}

#' Category Specific Predictors in \pkg{brms} Models
#' 
#' @aliases cse
#' 
#' @param expr Expression containing predictors,
#'  for which category specific effects should be estimated. 
#'  For evaluation, \R formula syntax is applied.
#'  
#' @details For detailed documentation see \code{help(brmsformula)}
#'   as well as \code{vignette("brms_overview")}.
#' 
#' This function is almost solely useful when
#' called in formulas passed to the \pkg{brms} package.
#' 
#' @seealso \code{\link{brmsformula}}
#'   
#' @examples   
#' \dontrun{
#' fit <- brm(rating ~ period + carry + cs(treat), 
#'            data = inhaler, family = sratio("cloglog"), 
#'            prior = set_prior("normal(0,5)"), chains = 2)
#' summary(fit)
#' plot(fit, ask = FALSE)
#' } 
#'  
#' @export
cs <- function(expr) {
  deparse_no_string(substitute(expr))
}

# alias of function 'cs' used in the JSS paper of brms
#' @export
cse <- function(expr) {
  deparse_no_string(substitute(expr))
}

#' Predictors with Measurement Error in \pkg{brms} Models
#' 
#' Specify predictors with measurement error. The function does not evaluate its
#' arguments -- it exists purely to help set up a model.
#' 
#' @param x The variable measured with error.
#' @param sdx Known measurement error of \code{x}
#'   treated as standard deviation.
#' @param gr Optional grouping factor to specify which
#'   values of \code{x} correspond to the same value of the
#'   latent variable. If \code{NULL} (the default) each
#'   observation will have its own value of the latent variable.
#' 
#' @details 
#' For detailed documentation see \code{help(brmsformula)}. 
#' 
#' By default, latent noise-free variables are assumed
#' to be correlated. To change that, add \code{set_mecor(FALSE)}
#' to your model formula object (see examples).
#' 
#' @seealso 
#' \code{\link{brmsformula}}, \code{\link{brmsformula-helpers}}
#'   
#' @examples 
#' \dontrun{
#' # sample some data
#' N <- 100
#' dat <- data.frame(
#'   y = rnorm(N), x1 = rnorm(N), 
#'   x2 = rnorm(N), sdx = abs(rnorm(N, 1))
#'  )
#' # fit a simple error-in-variables model 
#' fit1 <- brm(y ~ me(x1, sdx) + me(x2, sdx), data = dat, 
#'            save_mevars = TRUE)
#' summary(fit1)
#' 
#' # turn off modeling of correlations
#' bform <- bf(y ~ me(x1, sdx) + me(x2, sdx)) + set_mecor(FALSE)
#' fit2 <- brm(bform, data = dat, save_mevars = TRUE)
#' summary(fit2)
#' } 
#' 
#' @export
me <- function(x, sdx, gr = NULL) {
  # use 'term' for consistency with other special terms
  term <- deparse(substitute(x))
  sdx <- deparse(substitute(sdx))
  gr <- substitute(gr)
  if (!is.null(gr)) {
    gr <- deparse_combine(gr)
    stopif_illegal_group(gr)
  } else {
    gr <- ""
  }
  label <- deparse(match.call())
  out <- nlist(term, sdx, gr, label)
  class(out) <- c("me_term", "sp_term")
  out
}

#' Predictors with Missing Values in \pkg{brms} Models
#' 
#' Specify predictor term with missing values in \pkg{brms}. The function does
#' not evaluate its arguments -- it exists purely to help set up a model.
#' 
#' @param x The variable containing missings.
#' 
#' @details For detailed documentation see \code{help(brmsformula)}. 
#' 
#' @seealso \code{\link{brmsformula}}
#'   
#' @examples 
#' \dontrun{
#' data("nhanes", package = "mice")
#' bform <- bf(bmi | mi() ~ age * mi(chl)) +
#'   bf(chl | mi() ~ age) + set_rescor(FALSE)
#' fit <- brm(bform, data = nhanes)
#' summary(fit)
#' plot(marginal_effects(fit, resp = "bmi"), ask = FALSE)
#' LOO(fit, newdata = na.omit(fit$data))
#' } 
#' 
#' @export
mi <- function(x) {
  # use 'term' for consistency with other special terms
  term <- substitute(x)
  vars <- all.vars(term)
  term <- deparse(term)
  if (!is_equal(term, vars)) {
    stop2("'mi' only accepts single untransformed variables.")
  }
  label <- deparse(match.call())
  out <- nlist(term, label)
  class(out) <- c("mi_term", "sp_term")
  out
}

#' Monotonic Predictors in \pkg{brms} Models
#' 
#' Specify a monotonic predictor term in \pkg{brms}. The function does not
#' evaluate its arguments -- it exists purely to help set up a model.
#' 
#' @param x An integer variable or an ordered factor to be modeled as monotonic.
#'  
#' @details For detailed documentation see \code{help(brmsformula)}
#'   as well as \code{vignette("brms_monotonic")}.
#' 
#' @seealso \code{\link{brmsformula}}
#' 
#' @references 
#' Bürkner P. C. & Charpentier, E. (in review). Monotonic Effects: A Principled 
#' Approach for Including Ordinal Predictors in Regression Models. PsyArXiv 
#' preprint.
#'   
#' @examples   
#' \dontrun{
#' # generate some data
#' income_options <- c("below_20", "20_to_40", "40_to_100", "greater_100")
#' income <- factor(sample(income_options, 100, TRUE), 
#'                  levels = income_options, ordered = TRUE)
#' mean_ls <- c(30, 60, 70, 75)
#' ls <- mean_ls[income] + rnorm(100, sd = 7)
#' dat <- data.frame(income, ls)
#' 
#' # fit a simple monotonic model
#' fit1 <- brm(ls ~ mo(income), data = dat)
#' 
#' # summarise the model
#' summary(fit1)
#' plot(fit1, N = 6)
#' plot(marginal_effects(fit1), points = TRUE)
#' 
#' # model interaction with other variables
#' dat$x <- sample(c("a", "b", "c"), 100, TRUE)
#' fit2 <- brm(ls ~ mo(income)*x, data = dat)
#' 
#' # summarise the model
#' summary(fit2)
#' plot(marginal_effects(fit2), points = TRUE)
#' } 
#'  
#' @export
mo <- function(x) {
  # use 'term' for consistency with other special terms
  term <- deparse(substitute(x))
  label <- deparse(match.call())
  out <- nlist(term, label)
  class(out) <- c("mo_term", "sp_term")
  out
}

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
#'   default) of a non-isotropic GP should be used. 
#'   In the former case, the same amount of smoothing is applied to all
#'   predictors. In the latter case, predictors may have different smoothing.
#'   Ignored if only a single predictors is supplied.
#' @param gr Logical; Indicates if auto-grouping should be used (defaults 
#'   to \code{FALSE}). If enabled, observations sharing the same 
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
#'  interpreted as a prior over functions. Any finite sample 
#'  realized from this stochastic process is jointly multivariate 
#'  normal, with a covariance matrix defined by the covariance
#'  kernel \eqn{k_p(x)}, where \eqn{p} is the vector of parameters
#'  of the GP:
#'  \deqn{f(x) ~ MVN(0, k_p(x))}
#'  The smoothness and general behavior of the function \eqn{f} 
#'  depends only on the choice of covariance kernel. 
#'  For a more detailed introduction to Gaussian processes,
#'  see \url{https://en.wikipedia.org/wiki/Gaussian_process}.
#'  
#'  Below, we describe the currently supported covariance kernels:
#'  \itemize{
#'    \item{"exp_quad": }{The exponentiated-quadratic kernel is defined as
#'    \eqn{k(x_i, x_j) = sdgp^2 exp(- || x_i - x_j ||^2 / (2 lscale^2))},
#'    where \eqn{|| . ||} is the Euclidean norm, \eqn{sdgp} is a 
#'    standard deviation parameter, and \eqn{lscale} is characteristic 
#'    length-scale parameter. The latter practically measures how close two 
#'    points \eqn{x_i} and \eqn{x_j} have to be to influence each other 
#'    substantially.}
#'  }
#'
#'  In the current implementation, \code{"exp_quad"} is the only supported 
#'  covariance kernel. More options will follow in the future.
#'  
#' @return An object of class \code{'gpterm'}, which is a list 
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
#' me1 <- marginal_effects(fit1, nsamples = 200, spaghetti = TRUE)
#' plot(me1, ask = FALSE, points = TRUE)
#' 
#' # fit a more complicated GP model
#' fit2 <- brm(y ~ gp(x0) + x1 + gp(x2) + x3, dat, chains = 2)
#' summary(fit2)
#' me2 <- marginal_effects(fit2, nsamples = 200, spaghetti = TRUE)
#' plot(me2, ask = FALSE, points = TRUE)
#' 
#' # fit a multivariate GP model
#' fit3 <- brm(y ~ gp(x1, x2), dat, chains = 2)
#' summary(fit3)
#' me3 <- marginal_effects(fit3, nsamples = 200, spaghetti = TRUE)
#' plot(me3, ask = FALSE, points = TRUE)
#' 
#' # compare model fit
#' LOO(fit1, fit2, fit3)
#' 
#' # simulate data with a factor covariate
#' dat2 <- mgcv::gamSim(4, n = 90, scale = 2)
#' 
#' # fit separate gaussian processes for different levels of 'fac'
#' fit4 <- brm(y ~ gp(x2, by = fac), dat2, chains = 2)
#' summary(fit4)
#' plot(marginal_effects(fit4), points = TRUE)
#' }
#' 
#' @seealso \code{\link{brmsformula}}
#' @export
gp <- function(..., by = NA, k = NA, cov = "exp_quad", iso = TRUE, 
               gr = FALSE, cmc = TRUE, scale = TRUE, c = NULL) {
  cov <- match.arg(cov, choices = c("exp_quad"))
  label <- deparse(match.call())
  vars <- as.list(substitute(list(...)))[-1]
  by <- deparse(substitute(by))
  gr <- as_one_logical(gr)
  cmc <- as_one_logical(cmc)
  if (length(vars) > 1L) {
    iso <- as_one_logical(iso)
  } else {
    iso <- TRUE
  }
  if (!isNA(k)) {
    k <- as.integer(as_one_numeric(k))
    if (k < 1L) {
      stop2("'k' must be postive.")
    }
    if (is.null(c)) {
      stop2(
        "'c' must be specified for approximate GPs. ",
        "A good default could be c = 5/4 but we are still ",
        "working on providing better recommendations."
      )
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
  term <- ulapply(vars, deparse, backtick = TRUE, width.cutoff = 500)
  out <- nlist(term, label, by, cov, k, iso, gr, cmc, scale, c)
  structure(out, class = "gpterm")
}

#' Set up basic grouping terms in \pkg{brms}
#' 
#' Function used to set up a basic grouping term in \pkg{brms}.
#' The function does not evaluate its arguments --
#' it exists purely to help set up a model with grouping terms.
#' \code{gr} is called implicitly inside the package
#' and there is usually no need to call it directly.
#' 
#' @param ... One or more terms containing grouping factors.
#' @param by An optional factor variable, specifying sub-populations
#'   of the groups. For each level of the \code{by} variable, 
#'   a separate variance-covariance matrix will be fitted. 
#'   Levels of the grouping factor must be nested in levels 
#'   of the \code{by} variable.
#' @param dist Name of the distribution of the group-level effects.
#'   Currently \code{"gaussian"} is the only option.
#' 
#' @seealso \code{\link{brmsformula}}
#' 
#' @examples 
#' \dontrun{
#' # model using basic lme4-style formula
#' fit1 <- brm(count ~ Trt + (1|patient), data = epilepsy)
#' summary(fit1)
#' 
#' # equivalent model using 'gr' which is called anyway internally
#' fit2 <- brm(count ~ Trt + (1|gr(patient)), data = epilepsy)
#' summary(fit2)
#' 
#' # include Trt as a by variable
#' fit3 <- brm(count ~ Trt + (1|gr(patient, by = Trt)), data = epilepsy)
#' summary(fit3)
#' }
#' 
#' @export
gr <- function(..., by = NULL, dist = "gaussian") {
  label <- deparse(match.call())
  groups <- as.character(as.list(substitute(list(...)))[-1])
  if (length(groups) > 1L) {
    stop2("Grouping structure 'gr' expects only a single grouping term")
  }
  stopif_illegal_group(groups[1])
  by <- substitute(by)
  if (!is.null(by)) {
    by <- all.vars(by)
    if (length(by) != 1L) {
      stop2("Argument 'by' must contain exactly one variable.")
    }
  } else {
    by <- ""
  }
  dist <- match.arg(dist, c("gaussian", "student"))
  allvars <- str2formula(c(groups, by))
  nlist(groups, allvars, label, by, dist, type = "")
}

#' Set up multi-membership grouping terms in \pkg{brms}
#' 
#' Function to set up a multi-membership grouping term in \pkg{brms}.
#' The function does not evaluate its arguments --
#' it exists purely to help set up a model with grouping terms.
#'
#' @inheritParams gr
#' @param weights A matrix specifying the weights of each member.
#'  It should have as many columns as grouping terms specified in \code{...}.
#'  If \code{NULL} (the default), equally weights are used. 
#' @param scale Logical; if \code{TRUE} (the default), 
#'  weights are standardized in order to sum to one per row.
#'  If negative weights are specified, \code{scale} needs
#'  to be set to \code{FALSE}.
#'  
#' @seealso \code{\link{brmsformula}}, \code{\link{mmc}}
#'  
#' @examples 
#' \dontrun{
#' # simulate some data
#' dat <- data.frame(
#'  y = rnorm(100), x1 = rnorm(100), x2 = rnorm(100),
#'  g1 = sample(1:10, 100, TRUE), g2 = sample(1:10, 100, TRUE)
#' )
#' 
#' # multi-membership model with two members per group and equal weights
#' fit1 <- brm(y ~ x1 + (1|mm(g1, g2)), data = dat)
#' summary(fit1)
#' 
#' # weight the first member two times for than the second member
#' dat$w1 <- rep(2, 100)
#' dat$w2 <- rep(1, 100)
#' fit2 <- brm(y ~ x1 + (1|mm(g1, g2, weights = cbind(w1, w2))), data = dat)
#' summary(fit2)
#' 
#' # multi-membership model with level specific covariate values
#' dat$xc <- (dat$x1 + dat$x2) / 2
#' fit3 <- brm(y ~ xc + (1 + mmc(x1, x2) | mm(g1, g2)), data = dat)
#' summary(fit3)
#' }
#'   
#' @export
mm <- function(..., weights = NULL, scale = TRUE, dist = "gaussian") {
  label <- deparse(match.call())
  groups <- as.character(as.list(substitute(list(...)))[-1])
  if (length(groups) < 2) {
    stop2("Multi-membership terms require at least two grouping variables.")
  }
  for (i in seq_along(groups)) {
    stopif_illegal_group(groups[i])
  }
  dist <- match.arg(dist, c("gaussian", "student"))
  scale <- as_one_logical(scale)
  weights <- substitute(weights)
  weightvars <- all.vars(weights)
  allvars <- str2formula(c(groups, weightvars))
  if (!is.null(weights)) {
    weights <- str2formula(deparse_no_string(weights))
    attr(weights, "scale") <- scale
    weightvars <- str2formula(weightvars)
  }
  nlist(
    groups, weights, weightvars, allvars, 
    label, by = "", dist, type = "mm"
  )
}

#' Multi-Membership Covariates
#' 
#' Specify covarariates that vary over different levels 
#' of multi-membership grouping factors thus requiring
#' special treatment. This function is almost solely useful,
#' when called in combination with \code{\link{mm}}. 
#' Outside of multi-membership terms it will behave
#' very much like \code{\link{cbind}}.
#' 
#' @param ... One or more terms containing covariates 
#'   corresponding to the grouping levels specified in \code{\link{mm}}.
#' 
#' @return A matrix with covariates as columns.
#' 
#' @seealso \code{\link{mm}}
#' 
#' @examples 
#' \dontrun{
#' # simulate some data
#' dat <- data.frame(
#'   y = rnorm(100), x1 = rnorm(100), x2 = rnorm(100),
#'   g1 = sample(1:10, 100, TRUE), g2 = sample(1:10, 100, TRUE)
#' )
#' 
#' # multi-membership model with level specific covariate values
#' dat$xc <- (dat$x1 + dat$x2) / 2 
#' fit <- brm(y ~ xc + (1 + mmc(x1, x2) | mm(g1, g2)), data = dat)
#' summary(fit)
#' }
#' 
#' @export
mmc <- function(...) {
  dots <- list(...)
  if (any(ulapply(dots, is_like_factor))) {
    stop2("'mmc' requires numeric variables.")
  }
  out <- cbind(...)
  colnames(out) <- paste0("?", colnames(out))
  out
}

# return the right-hand side of a formula
rhs <- function(x) {
  attri <- attributes(x)
  x <- as.formula(x)
  x <- if (length(x) == 3) x[-2] else x
  do_call(structure, c(list(x), attri))
}

# return the left-hand side of a formula
lhs <- function(x) {
  x <- as.formula(x)
  if (length(x) == 3L) update(x, . ~ 1) else NULL
}

eval_rhs <- function(formula, data = NULL) {
  # computes data for addition arguments
  formula <- as.formula(formula)
  eval(rhs(formula)[[2]], data, environment(formula))
}

# convert a string to a formula
# @param x vector of strings to be converted
# @param ... passed to formula()
str2formula <- function(x, ..., collapse = "+") {
  has_chars <- nzchar(x)
  if (length(x) && any(has_chars)) {
    out <- paste(x[has_chars], collapse = collapse) 
  } else {
    out <- "1"
  }
  out <- formula(paste("~", out), ...)
  environment(out) <- parent.frame()
  out
}

# convert a formula to a character string
# @param formula a model formula
# @param rm a vector of to elements indicating how many characters 
#   should be removed at the beginning and end of the string respectively
# @param space how should whitespaces be treated?
formula2str <- function(formula, rm = c(0, 0), space = c("rm", "trim")) {
  space <- match.arg(space)
  if (!is.formula(formula)) {
    formula <- as.formula(formula)
  }
  if (is.na(rm[2])) rm[2] <- 0
  x <- Reduce(paste, deparse(formula))
  x <- gsub("[\t\r\n]+", "", x, perl = TRUE)
  if (space == "trim") {
    x <- gsub(" {1,}", " ", x, perl = TRUE)
  } else {
    x <- gsub(" ", "", x, perl = TRUE) 
  }
  substr(x, 1 + rm[1], nchar(x) - rm[2])
}

is.formula <- function(x) {
  inherits(x, "formula")
}

# expand the '.' variable in formula using stats::terms
expand_dot_formula <- function(formula, data = NULL) {
  if (isTRUE("." %in% all.vars(formula))) {
    att <- attributes(formula)
    try_terms <- try(
      stats::terms(formula, data = data), 
      silent = TRUE
    )
    if (!is(try_terms, "try-error")) {
      formula <- formula(try_terms)
    }
    attributes(formula) <- att
  }
  formula
}
