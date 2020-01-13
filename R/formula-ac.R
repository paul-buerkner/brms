#' @export
arma <- function(time = NA, gr = NA, p = 1, q = 1, cov = FALSE) {
  label <- deparse(match.call())
  time <- deparse(substitute(time))
  gr <- deparse(substitute(gr))
  if (gr != "NA") {
    stopif_illegal_group(gr)
  }
  p <- as_one_numeric(p)
  q <- as_one_numeric(q)
  if (!(p >= 0 && is_wholenumber(p))) {
    stop2("Autoregressive order must be a non-negative integer.")
  }
  if (!(q >= 0 && is_wholenumber(q))) {
    stop2("Moving-average order must be a non-negative integer.")
  }
  if (!sum(p, q)) {
    stop2("At least one of 'p' and 'q' should be greater zero.")
  }
  cov <- as_one_logical(cov)
  if (cov && (p > 1 || q > 1)) {
    stop2("Covariance formulation of ARMA structures is ", 
          "only possible for effects of maximal order one.")
  }
  out <- nlist(time, gr, p, q, cov, label)
  class(out) <- c("arma_term", "ac_term")
  out
}

#' @export
cosy <- function(time = NA, gr = NA) {
  label <- deparse(match.call())
  time <- deparse(substitute(time))
  gr <- deparse(substitute(gr))
  if (gr != "NA") {
    stopif_illegal_group(gr)
  }
  out <- nlist(time, gr, label)
  class(out) <- c("cosy_term", "ac_term")
  out
}

#' @export
sar <- function(W, type = "lag") {
  label <- deparse(match.call())
  W <- deparse(substitute(W))
  if (!nzchar(W)) {
    stop2("Argument 'W' is missing in 'sar'.")
  }
  options <- c("lag", "error")
  type <- match.arg(type, options)
  out <- nlist(W, type, label)
  class(out) <- c("sar_term", "ac_term")
  out
}

#' @export
car <- function(W, gr = NA, type = "escar") {
  label <- deparse(match.call())
  W <- deparse(substitute(W))
  if (!nzchar(W)) {
    stop2("Argument 'W' is missing in 'car'.")
  }
  gr <- deparse(substitute(gr))
  if (gr != "NA") {
    stopif_illegal_group(gr)
  }
  options <- c("escar", "esicar", "icar", "bym2")
  type <- match.arg(type, options)
  out <- nlist(W, gr, type, label)
  class(out) <- c("car_term", "ac_term")
  out
}

#' @export
fcor <- function(V) {
  # TODO: support estimating sigma additionally
  label <- deparse(match.call())
  V <- deparse(substitute(V))
  if (!nzchar(V)) {
    stop2("Argument 'V' is missing in 'fcor'.")
  }
  out <- nlist(V, label)
  class(out) <- c("fcor_term", "ac_term")
  out
}

# gather information on autocor terms
# @return a data.frame with one row per autocor term
tidy_acef <- function(x, ...) {
  UseMethod("tidy_acef")
}

#' @export
tidy_acef.default <- function(x, ...) {
  x <- parse_bf(x, check_response = FALSE)
  stopifnot(is.brmsterms(x))
  tidy_acef(x, ...)
}

tidy_acef.brmsterms <- function(x, ...) {
  tidy_acef(x$dpars$mu, ...)
}

#' @export
tidy_acef.acef <- function(x, ...) {
  x
}

#' @export
tidy_acef.NULL <- function(x, ...) {
  empty_acef()
}

#' @export
tidy_acef.btl <- function(x, data = NULL, ...) {
  form <- x[["ac"]]
  if (!is.formula(form)) {
    return(empty_acef())
  }
  out <- data.frame(term = all_terms(form), stringsAsFactors = FALSE)
  nterms <- NROW(out)
  cnames <- c("class", "dim", "type", "time", "gr", "p", "q", "W", "V")
  out[cnames] <- list(NA)
  out$cov <- FALSE
  for (i in seq_len(nterms)) {
    ac <- eval2(out$term[i])
    if (is.arma_term(ac)) {
      out$class[i] <- "arma"
      out$dim[i] <- "time"
      out$time[i] <- ac$time
      out$gr[i] <- ac$gr
      out$p[i] <- ac$p
      out$q[i] <- ac$q
      out$cov[i] <- ac$cov
    }
    if (is.cosy_term(ac)) {
      out$class[i] <- "cosy"
      out$dim[i] <- "time"
      out$time[i] <- ac$time
      out$gr[i] <- ac$gr
      out$cov[i] <- TRUE
    }
    if (is.sar_term(ac)) {
      out$class[i] <- "sar"
      out$dim[i] <- "space"
      out$type[i] <- ac$type
      out$W[i] <- ac$W
      out$cov[i] <- TRUE
    }
    if (is.car_term(ac)) {
      out$class[i] <- "car"
      out$dim[i] <- "space"
      out$type[i] <- ac$type
      out$gr[i] <- ac$gr
      out$W[i] <- ac$W
    }
    if (is.fcor_term(ac)) {
      out$class[i] <- "fcor"
      out$V[i] <- ac$V
      out$cov[i] <- TRUE
    }
  }
  if (sum(out$class %in% "arma") > 1L) {
    stop2("Formulas may not contain more than one ARMA term.")
  }
  if (sum(out$class %in% "cosy") > 1L) {
    stop2("Formulas may not contain more than one COSY term.")
  }
  acef_time_cov <- subset2(out, dim = "time", cov = TRUE)
  if (NROW(acef_time_cov) > 1) {
    stop2("Can only model one time-related covariance structure at a time.")
  }
  # TODO: add more validity checks
  class(out) <- c("acef", "data.frame")
  out
}

#' @export
tidy_acef.btnl <- function(x, ... ) {
  tidy_acef.btl(x, ...)
}

empty_acef <- function() {
  structure(empty_data_frame(), class = c("acef", "data.frame"))
}

# is certain subset of autocor terms is present?
has_ac_subset <- function(x, ...) {
  NROW(subset2(tidy_acef(x), ...)) > 0L
}

# is a certain autocorrelation class present?
has_ac_class <- function(x, class) {
  has_ac_subset(x, class = class)
}

# use explicit residual covariance structure?
use_ac_cov <- function(x) {
  has_ac_subset(x, cov = TRUE)
}

# use explicit residual covariance structure for time-series?
use_ac_cov_time <- function(x) {
  has_ac_subset(x, cov = TRUE, dim = "time")
}

# should natural residuals be modeled as correlated?
has_cor_natural_residuals <- function(bterms) {
  has_natural_residuals(bterms) && use_ac_cov(bterms)
}

# has the model latent residuals to be used in autocor structures
# TODO: rename to has_cor_latent_residuals?
has_latent_residuals <- function(bterms) {
  !has_natural_residuals(bterms) && use_ac_cov(bterms)
}

# helper function to prepare spatial weights matrices
sar_weights <- function(W) {
  if (is(W, "listw")) {
    require_package("spdep")
    W <- spdep::listw2mat(W)
  } else if (is(W, "nb")) {
    require_package("spdep")
    W <- spdep::nb2mat(W)
  }
  if (!is.matrix(W)) {
    stop2("'W' must be of class 'matrix', 'listw', or 'nb'.")
  }
  W
}

is.ac_term <- function(x) {
  inherits(x, "ac_term")
}

is.arma_term <- function(x) {
  inherits(x, "arma_term")
}

is.cosy_term <- function(x) {
  inherits(x, "cosy_term")
}

is.sar_term <- function(x) {
  inherits(x, "sar_term")
}

is.car_term <- function(x) {
  inherits(x, "car_term")
}

is.fcor_term <- function(x) {
  inherits(x, "fcor_term")
}
