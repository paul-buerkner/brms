# Most of the functions below have equivalents in Stan. Defining them in R is
# necessary to evaluate non-linear formulas containing these functions.

logit <- function(p) {
  # logit link
  log(p / (1 - p))
}

inv_logit <- function(x) { 
  # inverse of logit link
  1 / (1 + exp(-x))
}

cloglog <- function(x) {
  # cloglog link
  log(-log(1 - x))
}

inv_cloglog <- function(x) {
  # inverse of the cloglog link
  1 - exp(-exp(x))
}

Phi <- function(x) {
  pnorm(x)
}

incgamma <- function(x, a) {
  # incomplete gamma funcion
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
  log(1 + exp(x))
}

log1m_exp <- function(x) {
  ifelse(x < 0, log(1 - exp(x)), NaN)
}

log_diff_exp <- function(x, y) {
  stopifnot(length(x) == length(y))
  ifelse(x > y, log(exp(x) - exp(y)), NaN)
}

log_sum_exp <- function(x, y) {
  max <- max(x, y)
  max + log(exp(x - max) + exp(y - max))
}

log_mean_exp <- function(x) {
  # just log_sum_exp(x) - log(length(x))
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x))) - log(length(x))
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

softmax <- function(x) {
  if (!is.matrix(x)) {
    x <- matrix(x, nrow = 1)
  }
  x <- exp(x) 
  x / rowSums(x)
}

log_softmax <- function(x) {
  if (!is.matrix(x)) {
    x <- matrix(x, nrow = 1)
  }
  x - log(rowSums(exp(x)))
}
