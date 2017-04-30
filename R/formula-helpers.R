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
#'  \code{resp_se}, \code{resp_weights}, and \code{resp_disp} 
#'  require positive numeric values;
#'  \code{resp_trials} and \code{resp_cat} require positive integers;
#'  \code{resp_dec} requires \code{0} and \code{1}, or alternatively
#'  \code{'lower'} and \code{'upper'}; 
#'  \code{resp_cens} requires \code{'left'}, \code{'none'}, \code{'right'},
#'  and \code{'interval'} (or equivalenty \code{-1}, \code{0}, \code{1},
#'  and \code{2}) to indicate left, no, right, or interval censoring.
#' @param sigma Logical; Indicates whether the residual standard deviation
#'  parameter \code{sigma} should be included in addition to the known
#'  measurement error. Defaults to \code{FALSE} for backwards compatibility,
#'  but setting it to \code{TRUE} is usually the better choice.
#' @param y2 A vector specifying the upper bounds in interval censoring.
#' @param lb A numeric vector or single numeric value specifying 
#'   the lower truncation bound.
#' @param ub A numeric vector or single numeric value specifying 
#'   the upper truncation bound.
#'
#' @return A vector containing additional information on the response
#'   variable in an appropriate format.
#'
#' @details 
#'   These functions are almost solely useful when
#'   called in formulas passed to the \pkg{brms} package.
#'   Within formulas, the \code{resp_} prefix may be omitted.
#'   More information is given in the 'Details' section
#'   of \code{\link[brms:brmsformula]{brmsformula}}.
#'   
#' @seealso 
#'   \code{\link[brms:brm]{brm}}, 
#'   \code{\link[brms:brmsformula]{brmsformula}}   
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
#' fit4 <- brm(count | trunc(ub = 104) ~ log_Base4_c * Trt_c, 
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
  sigma <- as.logical(sigma)
  if (length(sigma) != 1L) {
    stop2("Argument 'sigma' must be either TRUE or FALSE.")
  }
  structure(x, sigma = sigma)  
}

resp_se_no_data <- function(x, sigma = FALSE) {
  # only evaluate the sigma argument
  resp_se(1, sigma = sigma)
}

#' @rdname addition-terms
#' @export
resp_weights <- function(x) {
  # weights to be applied on any model
  if (!is.numeric(x)) {
    stop2("Weights must be numeric.")
  }
  if (min(x) < 0) {
    stop2("Weights must be non-negative.")
  }
  x
}

#' @rdname addition-terms
#' @export
resp_disp <- function(x) {
  # dispersion factors
  if (!is.numeric(x)) {
    stop2("Dispersion factors must be numeric.")
  }
  if (min(x) < 0) {
    stop2("Dispersion factors must be non-negative.")
  }
  x  
}

#' @rdname addition-terms
#' @export
resp_trials <- function(x) {
  # trials for binomial models
  if (!is.numeric(x)) {
    stop2("Number of trials must be numeric.")
  }
  if (any(!is_wholenumber(x) || x < 1)) {
    stop2("Number of trials must be positive integers.")
  }
  x
}

#' @rdname addition-terms
#' @export
resp_cat <- function(x) {
  # number of categories for ordinal models
  if (!is.numeric(x)) {
    stop2("Number of categories must be numeric.")
  }
  if (length(x) != 1L || !is_wholenumber(x) || x < 1) {
    stop2("Number of categories must be a positive integer.")
  }
  x
}

#' @rdname addition-terms
#' @export
resp_dec <- function(x) {
  # decisions for the wiener diffusion model
  if (is.character(x) || is.factor(x)) {
    x <- ifelse(x == "lower", 0, ifelse(x == "upper", 1, x))
    if (!is.numeric(x)) {
      stop2("Decisions should be 'lower' or 'upper' ", 
            "when supplied as characters or factors.")
    }
  } else {
    x <- as.numeric(as.logical(x))
  }
  x
}

#' @rdname addition-terms
#' @export
resp_cens <- function(x, y2 = NULL) {
  # indicator for censoring
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
    stop2("Invalid censoring data. Accepted values are ", 
          "'left', 'none', 'right', and 'interval'\n",
          "(abbreviations are allowed) or -1, 0, 1, and 2.\n",
          "TRUE and FALSE are also accepted ",
          "and refer to 'right' and 'none' respectively.")
  }
  if (any(cens %in% 2)) {
    if (length(y2) != length(cens)) {
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

#' Predictors with Measurement Error in \pkg{brms} Models
#' 
#' @param x The variable measured with error.
#' @param sdx Known measurement error of \code{x}
#'   treated as standard deviation.
#' 
#' @details For detailed documentation see \code{help(brmsformula)}. 
#' 
#' This function is almost solely useful when
#' called in formulas passed to the \pkg{brms} package.
#' 
#' @seealso \code{\link[brms:brmsformula]{brmsformula}}
#'   
#' @examples 
#' \dontrun{
#' # sample some data
#' N <- 100
#' dat <- data.frame(y = rnorm(N), x = rnorm(N), sdx = abs(rnorm(N, 1)))
#' # fit a simple error-in-variables model 
#' fit <- brm(y ~ me(x, sdx), data = dat, save_mevars = TRUE)
#' summary(fit)
#' } 
#' 
#' @export
me <- function(x, sdx = NULL) {
  x <- as.vector(x)
  sdx <- as.vector(sdx)
  if (length(sdx) == 0L) {
    stop2("Argument 'sdx' is missing in function 'me'.")
  } else if (length(sdx) == 1L) {
    sdx <- rep(sdx, length(x))
  }
  if (!is.numeric(x)) {
    stop2("Noisy variables should be numeric.")
  }
  if (!is.numeric(sdx)) {
    stop2("Measurement error should be numeric.")
  }
  if (any(sdx <= 0)) {
    stop2("Measurement error should be positive.")
  }
  out <- rep(1, length(x))
  structure(out, var = x, noise = sdx) 
}

#' Category Specific Predictors in \pkg{brms} Models
#' 
#' @aliases cse
#' 
#' @param expr Expression containing predictors,
#'  for which category specific effects should
#'  be estimated. For evaluation, \R formula syntax is applied.
#'  
#' @details For detailed documentation see \code{help(brmsformula)}
#'   as well as \code{vignette("brms_overview")}.
#' 
#' This function is almost solely useful when
#' called in formulas passed to the \pkg{brms} package.
#' 
#' @seealso \code{\link[brms:brmsformula]{brmsformula}}
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

#' @export
cse <- function(expr) {
  # alias of function 'cs'
  deparse_no_string(substitute(expr))
}

#' Monotonic Predictors in \pkg{brms} Models
#' 
#' @aliases mono monotonic
#' 
#' @param expr Expression containing predictors,
#'  for which monotonic effects should
#'  be estimated. For evaluation, \R formula syntax is applied.
#'  
#' @details For detailed documentation see \code{help(brmsformula)}
#'   as well as \code{vignette("brms_monotonic")}.
#' 
#' This function is almost solely useful when
#' called in formulas passed to the \pkg{brms} package.
#' 
#' @seealso \code{\link[brms:brmsformula]{brmsformula}}
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
#' fit <- brm(ls ~ mo(income), data = dat)
#' 
#' # summarise the model
#' summary(fit)
#' plot(fit, N = 6)
#' plot(marginal_effects(fit), points = TRUE)
#' } 
#'  
#' @export
mo <- function(expr) {
  deparse_no_string(substitute(expr))
}

#' @export
mono <- function(expr) {
  # alias of function 'mo'
  deparse_no_string(substitute(expr))
}

#' @export
monotonic <- function(expr) {
  # alias of function 'mo'
  deparse_no_string(substitute(expr))
}

#' Set up Gaussian process terms in \pkg{brms}
#' 
#' Function used to set up a Gaussian process term in \pkg{brms}.
#' The function does not evaluate its arguments --
#' it exists purely to help set up a model with Gaussian process terms.
#' 
#' @param ... One or more predictors for the Gaussian process.
#' @param cov Name of the covariance kernel. By default, 
#'   the exponentiated-quadratic kernel \code{"exp_quad"} is used.
#' @param scale Logical; If \code{TRUE} (the default), predictors are
#'   scaled so that the maximum Euclidean distance between two points
#'   is 1. Since the default prior on \code{lscale} expects scaled
#'   predictors, it is recommended to manually specify priors
#'   on \code{lscale}, if \code{scale} is set to \code{FALSE}.
#'   
#' @details A Gaussian process is a stochastic process, whichs
#'  describes the relation between one or more predictors 
#'  \eqn{x = (x_1, ..., x_d)} and a response \eqn{f(x)}, where 
#'  \eqn{d} is the number of predictors. A Gaussian process is the
#'  generalization of the multivariate normal distribution
#'  to an infinite number of dimensions. Thus, it can be
#'  interpreted as a prior over functions. Any finite sample 
#'  realized from this stochastic process is jointly multivariate 
#'  normal, with a covariance matrix defined by the covariance
#'  kernel \eqn{k_p(x)}, where \eqn{p} is the vector of parameters
#'  of the Gaussian process:
#'  \deqn{f(x) ~ MVN(0, k_p(x))}
#'  The smoothness and general behavior of the function \eqn{f} 
#'  depends only on the choice of covariance kernel. 
#'  For a more detailed introduction to Gaussian processes,
#'  see \url{https://en.wikipedia.org/wiki/Gaussian_process}.
#'  
#'  Below, we describe the currently supported covariance kernels:
#'  \itemize{
#'    \item{"exp_quad": }{The exponentiated-quadratic kernel is defined as
#'    \eqn{k(x_i, x_j) = sdgp^2 exp(- || x_i - x_j || / (2 lscale^2)},
#'    where \eqn{|| . ||} is the Euclidean norm, \eqn{sdgp} is a 
#'    standard deviation parameter, and \eqn{lscale} is characteristic 
#'    length-scale parameter. The latter practically measures how close two 
#'    points \eqn{x_i} and \eqn{x_j} have to be to influence each other 
#'    substantially.}
#'  }
#'
#'  In the current implementation, only a single predictor can be 
#'  passed, and \code{"exp_quad"} is the only supported 
#'  covariance kernel. More options will follow in the future.  
#'  
#' @return An object of class \code{'gpterm'}, which is a list 
#'   of arguments to be interpreted by the formula 
#'   parsing functions of \code{brms}.
#'   
#' @examples
#' \dontrun{
#' # simulate data using the mgcv package
#' dat <- mgcv::gamSim(1, n = 30, scale = 2)
#' 
#' # fit a simple gaussian process model
#' fit1 <- brm(y ~ gp(x2), dat)
#' summary(fit1)
#' me1 <- marginal_effects(fit1, nsamples = 200, spaghetti = TRUE)
#' plot(me1, ask = FALSE, points = TRUE)
#' 
#' # fit a more complicated gaussian process model
#' fit2 <- brm(y ~ gp(x0) + x1 + gp(x2) + x3, dat)
#' summary(fit2)
#' me2 <- marginal_effects(fit2, nsamples = 200, spaghetti = TRUE)
#' plot(me2, ask = FALSE points = TRUE)
#' 
#' # compare model fit
#' LOO(fit1, fit2)
#' }
#' 
#' @seealso \code{\link[brms:brmsformula]{brmsformula}}
#' @export
gp <- function(..., cov = "exp_quad", scale = TRUE) {
  cov <- match.arg(cov)
  label <- deparse(match.call())
  vars <- as.list(substitute(list(...)))[-1]
  if (length(vars) != 1L) {
    stop2("'gp' currently allows only one predictor.")
  }
  scale <- as.logical(scale)
  if (anyNA(scale) || length(scale) != 1L) {
    stop2("'scale' should be either TRUE or FALSE.")
  }
  term <- ulapply(vars, deparse, backtick = TRUE, width.cutoff = 500)
  structure(nlist(term, label, cov, scale), class = "gpterm")
}

#' Set up basic grouping terms in \pkg{brms}
#' 
#' Function used to set up a basic grouping term in \pkg{brms}.
#' The function does not evaluate its arguments --
#' it exists purely to help set up a model with grouping terms.
#' \code{gr} is called implicitely inside the package
#' and there is usually no need to call it directly.
#' 
#' @param ... One or more terms containing grouping factors.
#' 
#' @seealso \code{\link[brms:brmsformula]{brmsformula}}
#' 
#' @examples 
#' \dontrun{
#' # model using basic lme4-style formula
#' fit1 <- brm(count ~ Trt_c + (1|patient), data = epilepsy)
#' summary(fit1)
#' # equivalent model using 'gr' which is called anyway internally
#' fit2 <- brm(count ~ Trt_c + (1|gr(patient)), data = epilepsy)
#' summary(fit2)
#' WAIC(fit1, fit2)
#' }
#' 
#' @export
gr <- function(...) {
  groups <- as.character(as.list(substitute(list(...)))[-1])
  if (length(groups) > 1L) {
    stop2("Grouping structure 'gr' expects only a single grouping term")
  }
  if (illegal_group_expr(groups[1])) {
    stop2("Illegal grouping term: ", groups[1], "\nIt may contain ",
          "only variable names combined by the symbol ':'")
  }
  allvars <- str2formula(groups)
  nlist(groups, allvars, type = "")
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
#'  Weights are standardized in order to sum to one per row.
#'  If \code{NULL} (the default), equally weights are used. 
#'  
#' @seealso \code{\link[brms:brmsformula]{brmsformula}}
#'  
#' @examples 
#' \dontrun{
#' # simulate some data
#' dat <- data.frame(y = rnorm(100), x = rnorm(100), 
#'                   g1 = sample(1:10, 100, TRUE),
#'                   g2 = sample(1:10, 100, TRUE))
#' # multi-membership model with two members per group and equal weights
#' fit1 <- brm(y ~ x + (1|mm(g1, g2)), data = dat)
#' summary(fit1)
#' 
#' # weight the first member two times for than the second member
#' dat$w1 <- rep(2, 100)
#' dat$w2 <- rep(1, 100)
#' fit2 <- brm(y ~ x + (1|mm(g1, g2, weights = cbind(w1, w2))), data = dat)
#' summary(fit2)
#' }
#'   
#' @export
mm <- function(..., weights = NULL) {
  groups <- as.character(as.list(substitute(list(...)))[-1])
  for (i in seq_along(groups)) {
    if (illegal_group_expr(groups[i])) {
      stop2("Illegal grouping term: ", groups[i], "\nIt may contain ",
            "only variable names combined by the symbol ':'")
    }
  }
  weights <- substitute(weights)
  weightvars <- all.vars(weights)
  allvars <- str2formula(c(groups, weightvars))
  if (!is.null(weights)) {
    weights <- str2formula(deparse_no_string(weights))
    weightvars <- str2formula(weightvars)
  }
  nlist(groups, weights, weightvars, allvars, type = "mm")
}

rhs <- function(x) {
  # return the righthand side of a formula
  attri <- attributes(x)
  x <- as.formula(x)
  x <- if (length(x) == 3) x[-2] else x
  do.call(structure, c(list(x), attri))
}

lhs <- function(x) {
  # return the lefthand side of a formula
  x <- as.formula(x)
  if (length(x) == 3L) update(x, . ~ 1) else NULL
}

eval_rhs <- function(formula, data = NULL) {
  # computes data for addition arguments
  formula <- as.formula(formula)
  eval(rhs(formula)[[2]], data, environment(formula))
}

str2formula <- function(x, ...) {
  # converts a string to a formula
  # Args:
  #   x: vector of strings to be converted
  #   ...: passed to formula(.)
  # Returns:
  #   a formula
  if (length(x)) {
    x <- paste(x, collapse = "+") 
  } else {
    x <- "1"
  }
  formula(paste("~", x), ...)
}

formula2str <- function(formula, rm = c(0, 0), space = c("rm", "trim")) {
  # converts a formula to a string
  # Args:
  #   formula: a model formula
  #   rm: a vector of to elements indicating how many characters 
  #       should be removed at the beginning
  #       and end of the string respectively
  #   space: how should whitespaces be treated?
  # Returns:
  #   a string
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
