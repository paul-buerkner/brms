# R helper functions for Gaussian Processes 

# get labels of gaussian process terms
# @param x either a formula or a list containing an element "gp"
# @param data data frame containing the covariates
# @return a data.frame with one row per GP term
tidy_gpef <- function(x, data) {
  if (is.formula(x)) {
    x <- parse_bf(x, check_response = FALSE)$dpars$mu
  }
  form <- x[["gp"]]
  if (!is.formula(form)) {
    return(empty_data_frame())
  }
  out <- data.frame(term = all_terms(form), stringsAsFactors = FALSE)
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
  out
}

# exponential-quadratic covariance matrix
# not vectorized over parameter values
cov_exp_quad <- function(x, x_new = NULL, sdgp = 1, lscale = 1) {
  Dls <- length(lscale)
  if (Dls == 1L) {
    # one dimensional or isotropic GP
    diff_quad <- diff_quad(x = x, x_new = x_new)
    out <- sdgp^2 * exp(-diff_quad / (2 * lscale^2))
  } else {
    # multi-dimensional non-isotropic GP
    diff_quad <- diff_quad(x = x[, 1], x_new = x_new[, 1])
    out <- sdgp^2 * exp(-diff_quad / (2 * lscale[1]^2))
    for (d in seq_len(Dls)[-1]) {
      diff_quad <- diff_quad(x = x[, d], x_new = x_new[, d])
      out <- out * exp(-diff_quad / (2 * lscale[d]^2))
    }
  }
  out
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
spd_cov_exp_quad <- function(x, sdgp = 1, lscale = 1) {
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
      x2 <- as_draws_matrix(x[m, ]^2, dim = dim(lscale))
      out[, m] <- constant * exp(rowSums(neg_half_lscale2 * x2))
    }
  }
  out
}

# compute the mth eigen value of an approximate GP
eigen_val_cov_exp_quad <- function(m, L) {
  ((m * pi) / (2 * L))^2
}

# compute the mth eigen function of an approximate GP
eigen_fun_cov_exp_quad <- function(x, m, L) {
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
  if (is(out, "try-error")) {
    stop2("The Gaussian process covariance matrix is not positive ", 
          "definite.\nThis occurs for numerical reasons. Setting ",
          "'nug' above ", nug, " may help.")
  }
  out
}
