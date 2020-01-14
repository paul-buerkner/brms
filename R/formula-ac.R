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
sar <- function(M, type = "lag") {
  label <- deparse(match.call())
  M <- deparse(substitute(M))
  if (!nzchar(M)) {
    stop2("Argument 'M' is missing in 'sar'.")
  }
  options <- c("lag", "error")
  type <- match.arg(type, options)
  out <- nlist(M, type, label)
  class(out) <- c("sar_term", "ac_term")
  out
}

#' @export
car <- function(M, gr = NA, type = "escar") {
  label <- deparse(match.call())
  M <- deparse(substitute(M))
  if (!nzchar(M)) {
    stop2("Argument 'M' is missing in 'car'.")
  }
  gr <- deparse(substitute(gr))
  if (gr != "NA") {
    stopif_illegal_group(gr)
  }
  options <- c("escar", "esicar", "icar", "bym2")
  type <- match.arg(type, options)
  out <- nlist(M, gr, type, label)
  class(out) <- c("car_term", "ac_term")
  out
}

#' @export
fcor <- function(M) {
  label <- deparse(match.call())
  M <- deparse(substitute(M))
  if (!nzchar(M)) {
    stop2("Argument 'M' is missing in 'fcor'.")
  }
  out <- nlist(M, label)
  class(out) <- c("fcor_term", "ac_term")
  out
}

# validate 'autocor' argument
validate_autocor <- function(autocor) {
  if (is.null(autocor) || is.cor_empty(autocor)) {
    return(NULL)
  }
  if (is.cor_brms(autocor)) {
    warning2("Using 'cor_brms' objects for 'autocor' is deprecated. ",
             "Please see ?cor_brms for details.")
    autocor <- as_formula_cor_brms(autocor)
  }
  if (is.null(autocor)) {
    return(NULL)
  }
  autocor <- as.formula(autocor)
  att <- attributes(autocor)
  autocor <- parse_ac(autocor)
  if (!is.null(autocor) && !is.formula(autocor)) {
    stop2("Argument 'autocor' must be coercible to a formula.")
  }
  attributes(autocor)[names(att)] <- att
  autocor
}

# gather information on autocor terms
# @return a data.frame with one row per autocor term
tidy_acef <- function(x, ...) {
  UseMethod("tidy_acef")
}

#' @export
tidy_acef.default <- function(x, ...) {
  x <- parse_bf(x, check_response = FALSE)
  tidy_acef(x, ...)
}

#' @export
tidy_acef.mvbrmsterms <- function(x, ...) {
  out <- lapply(x$terms, tidy_acef, ...)
  out <- do_call(rbind, out)
  structure(out, class = acef_class()) 
}

#' @export
tidy_acef.brmsterms <- function(x, ...) {
  out <- lapply(x$dpars, tidy_acef, ...)
  out <- do_call(rbind, out) 
  structure(out, class = acef_class())
}

#' @export
tidy_acef.btl <- function(x, data = NULL, ...) {
  form <- x[["ac"]]
  if (!is.formula(form)) {
    return(empty_acef())
  }
  px <- check_prefix(x)
  out <- data.frame(term = all_terms(form), stringsAsFactors = FALSE)
  nterms <- NROW(out)
  cnames <- c("class", "dim", "type", "time", "gr", "p", "q", "M")
  out[cnames] <- list(NA)
  out$cov <- FALSE
  out[names(px)] <- px
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
      out$M[i] <- ac$M
      out$cov[i] <- TRUE
    }
    if (is.car_term(ac)) {
      out$class[i] <- "car"
      out$dim[i] <- "space"
      out$type[i] <- ac$type
      out$gr[i] <- ac$gr
      out$M[i] <- ac$M
    }
    if (is.fcor_term(ac)) {
      out$class[i] <- "fcor"
      out$M[i] <- ac$M
      out$cov[i] <- TRUE
    }
  }
  if (any(duplicated(out$class))) {
    stop2("Can only model one term per autocorrelation class.")
  }
  if (NROW(subset2(out, dim = "time")) > 1) {
    stop2("Can only model one time-series term.")
  }
  if (NROW(subset2(out, dim = "space")) > 1) {
    stop2("Can only model one spatial term.")
  }
  if (NROW(subset2(out, cov = TRUE)) > 1) {
    stop2("Can only model one explicit covariance term.")
  }
  if (NROW(subset2(out, cov = TRUE))) {
    if (!dpar_class(px$dpar) %in% c("", "mu") || nzchar(px$nlpar)) {
      stop2("Explicit covariance terms can only be specified on 'mu'.")
    }
  }
  structure(out, class = acef_class())
}

#' @export
tidy_acef.btnl <- function(x, ... ) {
  tidy_acef.btl(x, ...)
}

#' @export
tidy_acef.acef <- function(x, ...) {
  x
}

#' @export
tidy_acef.NULL <- function(x, ...) {
  empty_acef()
}

empty_acef <- function() {
  structure(empty_data_frame(), class = acef_class())
}

acef_class <- function() {
  c("acef", "data.frame")
}

# get names of certain autocor variables
get_ac_vars <- function(x, var, ...) {
  var <- match.arg(var, c("time", "gr", "M"))
  acef <- subset2(tidy_acef(x), ...)
  out <- unique(acef[[var]])
  setdiff(na.omit(out), "NA")
}

# get names of autocor grouping variables
get_ac_groups <- function(x, ...) {
  get_ac_vars(x, "gr", ...)
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

# validate SAR matrices
validate_sar_matrix <- function(M) {
  if (is(M, "listw")) {
    require_package("spdep")
    M <- spdep::listw2mat(M)
  } else if (is(M, "nb")) {
    require_package("spdep")
    M <- spdep::nb2mat(M)
  }
  if (length(dim(M)) != 2L) {
    stop2("'M' for SAR terms must be of class 'matrix', 'listw', or 'nb'.")
  }
  M <- Matrix::Matrix(M, sparse = TRUE)
  if (!Matrix::isSymmetric(M, check.attributes = FALSE)) {
    stop2("'M' for CAR terms must be symmetric.")
  }
  M
}

# validate CAR matrices
validate_car_matrix <- function(M) {
  if (length(dim(M)) != 2L) {
    stop2("'M' for CAR terms must be a matrix.")
  }
  M <- Matrix::Matrix(M, sparse = TRUE)
  if (!Matrix::isSymmetric(M, check.attributes = FALSE)) {
    stop2("'M' for CAR terms must be symmetric.")
  }
  if (is.null(rownames(M))) {
    stop2("Row names are required for 'M' for CAR terms.")
  }
  colnames(M) <- rownames(M)
  not_binary <- !(M == 0 | M == 1)
  if (any(not_binary)) {
    message("Converting all non-zero values in 'M' to 1")
    M[not_binary] <- 1
  }
  M
}

# validate FCOR matrices
validate_fcor_matrix <- function(M) {
  if (length(dim(M)) <= 1L) {
    M <- diag(as.vector(M), length(M))
  }
  if (length(dim(M)) != 2L) {
    stop2("'M' for FCOR terms must be a matrix.")
  }
  M <- as.matrix(M)
  if (!isSymmetric(M, check.attributes = FALSE)) {
    stop2("'M' for FCOR terms must be symmetric.")
  }
  if (min(eigen(M)$values <= 0)) {
    stop2("'M' for FCOR terms must be positive definite.")
  }
  M
}

# regex to extract all parameter names of autocorrelation structures
regex_autocor_pars <- function() {
  p <- c("ar", "ma", "sderr", "cosy", "lagsar", "errorsar", "car", "sdcar")
  p <- paste0("(", p, ")", collapse = "|")
  paste0("^(", p, ")(\\[|_|$)")
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
