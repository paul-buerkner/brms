# R helper functions for Gaussian Processes 

cov_exp_quad <- function(x, x_new = NULL, sdgp = 1, lscale = 1) {
  # exponential-quadratic covariance matrix
  diff_quad <- diff_quad(x = x, x_new = x_new)
  sdgp^2 * exp(-diff_quad / (2 * lscale^2))
}

diff_quad <- function(x, x_new = NULL) {
  # compute squared differences
  # Args:
  #   x: vector or matrix
  #   x_new: optional vector of matrix with the same ncol as x
  # Returns:
  #   An nrow(x) times nrow(x_new) matrix
  # Details:
  #   If matrices are passed results are summed over the columns
  x <- as.matrix(x)
  if (is.null(x_new)) {
    x_new <- x
  } else {
    x_new <- as.matrix(x_new)
  }
  .diff_quad <- function(x1, x2) (x1 - x2)^2
  out <- 0
  for (i in seq_len(ncol(x))) {
    out <- out + outer(x[, i], x_new[, i], .diff_quad)
  }
  out
}

spd_cov_exp_quad <- function(x, sdgp = 1, lscale = 1) {
  # spectral density function
  # vectorized over parameter values
  out <- matrix(NA, nrow = length(sdgp), ncol = length(x))
  for (m in seq_along(x)) {
    out[, m] <- sdgp^2 * sqrt(2 * pi) * 
      lscale * exp(-0.5 * lscale^2 * x[m]^2);
  }
  out
}

eigen_val_cov_exp_quad <- function(m, L) {
  # compute the mth eigen value of an approximate GP
  ((m * pi) / (2 * L))^2
}

eigen_fun_cov_exp_quad <- function(x, m, L) {
  # compute the mth eigen function of an approximate GP
  1 / sqrt(L) * sin((m * pi) / (2 * L) * (x + L))
}

choose_L <- function(x, L) {
  # extended range of input data for which predictions should be made
  if (!length(x)) {
    range <- 1
  } else {
    range <- max(1, max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }
  L * range
}

try_nug <- function(expr, nug) {
  # try to evaluate a GP term and 
  # return an informative error message if it fails
  out <- try(expr, silent = TRUE)
  if (is(out, "try-error")) {
    stop2("The Gaussian process covariance matrix is not positive ", 
          "definite.\nThis occurs for numerical reasons. Setting ",
          "'nug' above ", nug, " may help.")
  }
  out
}
