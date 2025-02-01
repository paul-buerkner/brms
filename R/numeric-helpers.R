# Most of the functions below have equivalents in Stan. Defining them in R is
# necessary to evaluate non-linear formulas containing these functions.

logit <- function(p) {
  log(p) - log1p(-p)
}

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

cloglog <- function(x) {
  log(-log1p(-x))
}

inv_cloglog <- function(x) {
  1 - exp(-exp(x))
}

Phi <- function(x) {
  pnorm(x)
}

# incomplete gamma funcion
incgamma <- function(a, x) {
  pgamma(x, shape = a) * gamma(a)
}

square <- function(x) {
  x^2
}

cbrt <- function(x) {
  x^(1/3)
}

exp2 <- function(x) {
  2^x
}

pow <- function(x, y) {
  x^y
}

inv <- function(x) {
  1/x
}

inv_sqrt <- function(x) {
  1/sqrt(x)
}

inv_square <- function(x) {
  1/x^2
}

hypot <- function(x, y) {
  stopifnot(all(x >= 0))
  stopifnot(all(y >= 0))
  sqrt(x^2 + y^2)
}

log1m <- function(x) {
  log(1 - x)
}

step <- function(x) {
  if (is.logical(x)) {
    warning(
      "Consider using int_step() instead of step(), ",
      "because step(FALSE) = step(TRUE) = 1 in Stan"
    )
  }
  ifelse(x >= 0, 1, 0)
}

int_step <- function(x) {
  ifelse(x > 0, 1, 0)
}

#' Logarithm with a minus one offset.
#'
#' Computes \code{log(x - 1)}.
#'
#' @param x A numeric or complex vector.
#' @param base A positive or complex number: the base with respect to which
#'   logarithms are computed. Defaults to \emph{e} = \code{exp(1)}.
#'
#' @export
logm1 <- function(x, base = exp(1)) {
  log(x - 1, base = base)
}

#' Exponential function plus one.
#'
#' Computes \code{exp(x) + 1}.
#'
#' @param x A numeric or complex vector.
#'
#' @export
expp1 <- function(x) {
  exp(x) + 1
}

#' Scaled logit-link
#'
#' Computes \code{logit((x - lb) / (ub - lb))}
#'
#' @param x A numeric or complex vector.
#' @param lb Lower bound defaulting to \code{0}.
#' @param ub Upper bound defaulting to \code{1}.
#'
#' @return A numeric or complex vector.
#'
#' @export
logit_scaled <- function(x, lb = 0, ub = 1) {
  logit((x - lb) / (ub - lb))
}

#' Scaled inverse logit-link
#'
#' Computes \code{inv_logit(x) * (ub - lb) + lb}
#'
#' @param x A numeric or complex vector.
#' @param lb Lower bound defaulting to \code{0}.
#' @param ub Upper bound defaulting to \code{1}.
#'
#' @return A numeric or complex vector between \code{lb} and \code{ub}.
#'
#' @export
inv_logit_scaled <- function(x, lb = 0, ub = 1) {
  inv_logit(x) * (ub - lb) + lb
}

multiply_log <- function(x, y) {
  ifelse(x == y & x == 0, 0, x * log(y))
}

log1p_exp <- function(x) {
  # approaches identity(x) for x -> Inf
  out <- log1p(exp(x))
  ifelse(out < Inf, out, x)
}

log1m_exp <- function(x) {
  ifelse(x < 0, log1p(-exp(x)), NaN)
}

log_diff_exp <- function(x, y) {
  stopifnot(length(x) == length(y))
  ifelse(x > y, log(exp(x) - exp(y)), NaN)
}

log_sum_exp <- function(x, y) {
  max <- pmax(x, y)
  max + log(exp(x - max) + exp(y - max))
}

log_mean_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x))) - log(length(x))
}

log_expm1 <- function(x) {
  # approaches identity(x) for x -> Inf
  out <- log(expm1(x))
  ifelse(out < Inf, out, x)
}

log_inv_logit <- function(x) {
  log(inv_logit(x))
}

log1m_inv_logit <- function(x) {
  log(1 - inv_logit(x))
}

scale_unit <- function(x, lb = min(x), ub = max(x)) {
  (x - lb) / (ub - lb)
}

fabs <- function(x) {
  abs(x)
}

log_softmax <- function(x) {
  ndim <- length(dim(x))
  if (ndim <= 1) {
    x <- matrix(x, nrow = 1)
    ndim <- length(dim(x))
  }
  dim_noncat <- dim(x)[-ndim]
  marg_noncat <- seq_along(dim(x))[-ndim]
  catsum <- log(array(apply(exp(x), marg_noncat, sum), dim = dim_noncat))
  sweep(x, marg_noncat, catsum, "-")
}

softmax <- function(x) {
  # log_softmax is more numerically stable #1401
  exp(log_softmax(x))
}

inv_odds <- function(x) {
  x / (1 + x)
}

# inspired by logit but with softplus instead of log
softit <- function(x) {
  log_expm1(x / (1 - x))
}

# inspired by inv_logit but with softplus instead of exp
inv_softit <- function(x) {
  y <- log1p_exp(x)
  y / (1 + y)
}

# inspired by inv_logit but with softplus instead of exp
log_inv_softit <- function(x) {
  y <- log1p_exp(x)
  log(y) - log1p(y)
}

# inspired by inv_logit but with softplus instead of exp
log1m_inv_softit <- function(x) {
  y <- log1p_exp(x)
  -log1p(y)
}

# names of built-in stan functons reimplemented in R within brms
names_stan_functions <- function() {
  c("logit", "inv_logit", "cloglog", "inv_cloglog", "Phi", "incgamma",
    "square", "cbrt", "exp2", "pow", "inv", "inv_sqrt", "inv_square",
    "hypot", "log1m", "step", "logm1", "expp1", "logit_scaled",
    "inv_logit_scaled", "multiply_log", "log1p_exp", "log1m_exp",
    "log_diff_exp", "log_sum_exp", "log_mean_exp", "log_expm1",
    "log_inv_logit", "log1m_inv_logit", "scale_unit", "fabs", "log_softmax",
    "softmax", "inv_odds", "softit", "inv_softit", "log_inv_softit",
    "log1m_inv_softit")
}

# create an environement with all the reimplemented stan functions in it
# see issue #1635 for discussion of this approach
env_stan_functions <- function(...) {
  env <- new.env(...)
  brms_env <- asNamespace("brms")
  for (f in names_stan_functions()) {
    env[[f]] <- get(f, brms_env)
  }
  env
}
