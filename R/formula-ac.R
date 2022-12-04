#' Autocorrelation structures
#'
#' Specify autocorrelation terms in \pkg{brms} models. Currently supported terms
#' are \code{\link{arma}}, \code{\link{ar}}, \code{\link{ma}},
#' \code{\link{cosy}}, \code{\link{sar}}, \code{\link{car}}, and
#' \code{\link{fcor}}. Terms can be directly specified within the formula, or
#' passed to the \code{autocor} argument of \code{\link{brmsformula}} in the
#' form of a one-sided formula. For deprecated ways of specifying
#' autocorrelation terms, see \code{\link{cor_brms}}.
#'
#' @name autocor-terms
#'
#' @details The autocor term functions are almost solely useful when called in
#' formulas passed to the \pkg{brms} package. They do not evaluate its
#' arguments -- but exist purely to help set up a model with autocorrelation
#' terms.
#'
#' @seealso \code{\link{brmsformula}}, \code{\link{acformula}},
#'   \code{\link{arma}}, \code{\link{ar}}, \code{\link{ma}},
#'   \code{\link{cosy}}, \code{\link{sar}}, \code{\link{car}},
#'   \code{\link{fcor}}
#'
#' @examples
#' # specify autocor terms within the formula
#' y ~ x + arma(p = 1, q = 1) + car(M)
#'
#' # specify autocor terms in the 'autocor' argument
#' bf(y ~ x, autocor = ~ arma(p = 1, q = 1) + car(M))
#'
#' # specify autocor terms via 'acformula'
#' bf(y ~ x) + acformula(~ arma(p = 1, q = 1) + car(M))
NULL

#' Set up ARMA(p,q) correlation structures
#'
#' Set up an autoregressive moving average (ARMA) term of order (p, q) in
#' \pkg{brms}. The function does not evaluate its arguments -- it exists purely
#' to help set up a model with ARMA terms.
#'
#' @param time An optional time variable specifying the time ordering
#'   of the observations. By default, the existing order of the observations
#'   in the data is used.
#' @param gr An optional grouping variable. If specified, the correlation
#'   structure is assumed to apply only to observations within the same grouping
#'   level.
#' @param p A non-negative integer specifying the autoregressive (AR)
#'   order of the ARMA structure. Default is \code{1}.
#' @param q A non-negative integer specifying the moving average (MA)
#'   order of the ARMA structure. Default is \code{1}.
#' @param cov A flag indicating whether ARMA effects should be estimated by
#'   means of residual covariance matrices. This is currently only possible for
#'   stationary ARMA effects of order 1. If the model family does not have
#'   natural residuals, latent residuals are added automatically. If
#'   \code{FALSE} (the default), a regression formulation is used that is
#'   considerably faster and allows for ARMA effects of order higher than 1 but
#'   is only available for \code{gaussian} models and some of its
#'   generalizations.
#'
#' @return An object of class \code{'arma_term'}, which is a list
#'   of arguments to be interpreted by the formula
#'   parsing functions of \pkg{brms}.
#'
#' @seealso \code{\link{autocor-terms}}, \code{\link{ar}}, \code{\link{ma}},
#'
#' @examples
#' \dontrun{
#' data("LakeHuron")
#' LakeHuron <- as.data.frame(LakeHuron)
#' fit <- brm(x ~ arma(p = 2, q = 1), data = LakeHuron)
#' summary(fit)
#' }
#'
#' @export
arma <- function(time = NA, gr = NA, p = 1, q = 1, cov = FALSE) {
  label <- deparse(match.call())
  time <- deparse(substitute(time))
  gr <- deparse(substitute(gr))
  .arma(time = time, gr = gr, p = p, q = q, cov = cov, label = label)
}

#' Set up AR(p) correlation structures
#'
#' Set up an autoregressive (AR) term of order p in \pkg{brms}. The function
#' does not evaluate its arguments -- it exists purely to help set up a model
#' with AR terms.
#'
#' @inheritParams arma
#'
#' @return An object of class \code{'arma_term'}, which is a list
#'   of arguments to be interpreted by the formula
#'   parsing functions of \pkg{brms}.
#'
#' @seealso \code{\link{autocor-terms}}, \code{\link{arma}}, \code{\link{ma}}
#'
#' @examples
#' \dontrun{
#' data("LakeHuron")
#' LakeHuron <- as.data.frame(LakeHuron)
#' fit <- brm(x ~ ar(p = 2), data = LakeHuron)
#' summary(fit)
#' }
#'
#' @export
ar <- function(time = NA, gr = NA, p = 1, cov = FALSE) {
  label <- deparse(match.call())
  time <- deparse(substitute(time))
  gr <- deparse(substitute(gr))
  .arma(time = time, gr = gr, p = p, q = 0, cov = cov, label = label)
}

#' Set up MA(q) correlation structures
#'
#' Set up a moving average (MA) term of order q in \pkg{brms}. The function does
#' not evaluate its arguments -- it exists purely to help set up a model with
#' MA terms.
#'
#' @inheritParams arma
#'
#' @return An object of class \code{'arma_term'}, which is a list
#'   of arguments to be interpreted by the formula
#'   parsing functions of \pkg{brms}.
#'
#' @seealso \code{\link{autocor-terms}}, \code{\link{arma}}, \code{\link{ar}}
#'
#' @examples
#' \dontrun{
#' data("LakeHuron")
#' LakeHuron <- as.data.frame(LakeHuron)
#' fit <- brm(x ~ ma(p = 2), data = LakeHuron)
#' summary(fit)
#' }
#'
#' @export
ma <- function(time = NA, gr = NA, q = 1, cov = FALSE) {
  label <- deparse(match.call())
  time <- deparse(substitute(time))
  gr <- deparse(substitute(gr))
  .arma(time = time, gr = gr, p = 0, q = q, cov = cov, label = label)
}

# helper function to validate input to arma()
.arma <- function(time, gr, p, q, cov, label) {
  time <- as_one_variable(time)
  gr <- as_one_character(gr)
  stopif_illegal_group(gr)
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
  label <- as_one_character(label)
  out <- nlist(time, gr, p, q, cov, label)
  class(out) <- c("arma_term", "ac_term")
  out
}

#' Set up COSY correlation structures
#'
#' Set up a compounds symmetry (COSY) term in \pkg{brms}. The function does
#' not evaluate its arguments -- it exists purely to help set up a model with
#' COSY terms.
#'
#' @inheritParams arma
#'
#' @return An object of class \code{'cosy_term'}, which is a list
#'   of arguments to be interpreted by the formula
#'   parsing functions of \pkg{brms}.
#'
#' @seealso \code{\link{autocor-terms}}
#'
#' @examples
#' \dontrun{
#' data("lh")
#' lh <- as.data.frame(lh)
#' fit <- brm(x ~ cosy(), data = lh)
#' summary(fit)
#' }
#'
#' @export
cosy <- function(time = NA, gr = NA) {
  label <- deparse(match.call())
  time <- deparse(substitute(time))
  time <- as_one_variable(time)
  gr <- deparse(substitute(gr))
  stopif_illegal_group(gr)
  out <- nlist(time, gr, label)
  class(out) <- c("cosy_term", "ac_term")
  out
}

#' @export
unstr <- function(time = NA, gr = NA) {
  label <- deparse(match.call())
  time <- deparse(substitute(time))
  time <- as_one_variable(time)
  gr <- deparse(substitute(gr))
  stopif_illegal_group(gr)
  out <- nlist(time, gr, label)
  class(out) <- c("unstr_term", "ac_term")
  out
}

#' Spatial simultaneous autoregressive (SAR) structures
#'
#' Set up an spatial simultaneous autoregressive (SAR) term in \pkg{brms}. The
#' function does not evaluate its arguments -- it exists purely to help set up a
#' model with SAR terms.
#'
#' @param M An object specifying the spatial weighting matrix.
#'   Can be either the spatial weight matrix itself or an
#'   object of class \code{listw} or \code{nb}, from which
#'   the spatial weighting matrix can be computed.
#' @param type Type of the SAR structure. Either \code{"lag"}
#'   (for SAR of the response values) or \code{"error"}
#'   (for SAR of the residuals). More information is
#'   provided in the 'Details' section.
#'
#' @details The \code{lagsar} structure implements SAR of the response values:
#'   \deqn{y = \rho W y + \eta + e}
#'   The \code{errorsar} structure implements SAR of the residuals:
#'   \deqn{y = \eta + u, u = \rho W u + e}
#'   In the above equations, \eqn{\eta} is the predictor term and \eqn{e} are
#'   independent normally or t-distributed residuals. Currently, only families
#'   \code{gaussian} and \code{student} support SAR structures.
#'
#' @return An object of class \code{'sar_term'}, which is a list
#'   of arguments to be interpreted by the formula
#'   parsing functions of \pkg{brms}.
#'
#' @seealso \code{\link{autocor-terms}}
#'
#' @examples
#' \dontrun{
#' data(oldcol, package = "spdep")
#' fit1 <- brm(CRIME ~ INC + HOVAL + sar(COL.nb, type = "lag"),
#'             data = COL.OLD, data2 = list(COL.nb = COL.nb),
#'             chains = 2, cores = 2)
#' summary(fit1)
#' plot(fit1)
#'
#' fit2 <- brm(CRIME ~ INC + HOVAL + sar(COL.nb, type = "error"),
#'             data = COL.OLD, data2 = list(COL.nb = COL.nb),
#'             chains = 2, cores = 2)
#' summary(fit2)
#' plot(fit2)
#' }
#'
#' @export
sar <- function(M, type = "lag") {
  label <- deparse(match.call())
  if (missing(M)) {
    stop2("Argument 'M' is missing in sar().")
  }
  M <- deparse(substitute(M))
  M <- as_one_variable(M)
  options <- c("lag", "error")
  type <- match.arg(type, options)
  out <- nlist(M, type, label)
  class(out) <- c("sar_term", "ac_term")
  out
}

#' Spatial conditional autoregressive (CAR) structures
#'
#' Set up an spatial conditional autoregressive (CAR) term in \pkg{brms}. The
#' function does not evaluate its arguments -- it exists purely to help set up a
#' model with CAR terms.
#'
#' @param M Adjacency matrix of locations. All non-zero entries are treated as
#'   if the two locations are adjacent. If \code{gr} is specified, the row names
#'   of \code{M} have to match the levels of the grouping factor.
#' @param gr An optional grouping factor mapping observations to spatial
#'   locations. If not specified, each observation is treated as a separate
#'   location. It is recommended to always specify a grouping factor to allow
#'   for handling of new data in post-processing methods.
#' @param type Type of the CAR structure. Currently implemented are
#'   \code{"escar"} (exact sparse CAR), \code{"esicar"} (exact sparse intrinsic
#'   CAR), \code{"icar"} (intrinsic CAR), and \code{"bym2"}. More information is
#'   provided in the 'Details' section.
#'
#' @return An object of class \code{'car_term'}, which is a list
#'   of arguments to be interpreted by the formula
#'   parsing functions of \pkg{brms}.
#'
#' @seealso \code{\link{autocor-terms}}
#'
#' @details The \code{escar} and \code{esicar} types are
#'   implemented based on the case study of Max Joseph
#'   (\url{https://github.com/mbjoseph/CARstan}). The \code{icar} and
#'   \code{bym2} type is implemented based on the case study of Mitzi Morris
#'   (\url{https://mc-stan.org/users/documentation/case-studies/icar_stan.html}).
#'
#' @examples
#' \dontrun{
#' # generate some spatial data
#' east <- north <- 1:10
#' Grid <- expand.grid(east, north)
#' K <- nrow(Grid)
#'
#' # set up distance and neighbourhood matrices
#' distance <- as.matrix(dist(Grid))
#' W <- array(0, c(K, K))
#' W[distance == 1] <- 1
#'
#' # generate the covariates and response data
#' x1 <- rnorm(K)
#' x2 <- rnorm(K)
#' theta <- rnorm(K, sd = 0.05)
#' phi <- rmulti_normal(
#'   1, mu = rep(0, K), Sigma = 0.4 * exp(-0.1 * distance)
#' )
#' eta <- x1 + x2 + phi
#' prob <- exp(eta) / (1 + exp(eta))
#' size <- rep(50, K)
#' y <- rbinom(n = K, size = size, prob = prob)
#' dat <- data.frame(y, size, x1, x2)
#'
#' # fit a CAR model
#' fit <- brm(y | trials(size) ~ x1 + x2 + car(W),
#'            data = dat, data2 = list(W = W),
#'            family = binomial())
#' summary(fit)
#' }
#'
#' @export
car <- function(M, gr = NA, type = "escar") {
  label <- deparse(match.call())
  if (missing(M)) {
    stop2("Argument 'M' is missing in car().")
  }
  M <- deparse(substitute(M))
  M <- as_one_variable(M)
  gr <- deparse(substitute(gr))
  stopif_illegal_group(gr)
  options <- c("escar", "esicar", "icar", "bym2")
  type <- match.arg(type, options)
  out <- nlist(M, gr, type, label)
  class(out) <- c("car_term", "ac_term")
  out
}

#' Fixed residual correlation (FCOR) structures
#'
#' Set up a fixed residual correlation (FCOR) term in \pkg{brms}. The function
#' does not evaluate its arguments -- it exists purely to help set up a model
#' with FCOR terms.
#'
#' @param M Known correlation/covariance matrix of the response variable.
#'   If a vector is passed, it will be used as diagonal entries
#'   (variances) and correlations/covariances will be set to zero.
#'   The actual covariance matrix used in the likelihood is obtained
#'   by multiplying \code{M} by the square of the residual standard
#'   deviation parameter \code{sigma} estimated as part of the model.
#'
#' @return An object of class \code{'fcor_term'}, which is a list
#'   of arguments to be interpreted by the formula
#'   parsing functions of \pkg{brms}.
#'
#' @seealso \code{\link{autocor-terms}}
#'
#' @examples
#' \dontrun{
#' dat <- data.frame(y = rnorm(3))
#' V <- cbind(c(0.5, 0.3, 0.2), c(0.3, 1, 0.1), c(0.2, 0.1, 0.2))
#' fit <- brm(y ~ 1 + fcor(V), data = dat, data2 = list(V = V))
#' }
#'
#' @export
fcor <- function(M) {
  label <- deparse(match.call())
  if (missing(M)) {
    stop2("Argument 'M' is missing in fcor().")
  }
  M <- deparse(substitute(M))
  M <- as_one_variable(M)
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
  autocor <- terms_ac(autocor)
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
  x <- brmsterms(x, check_response = FALSE)
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
  if (!NROW(out)) {
    return(empty_acef())
  }
  out <- structure(out, class = acef_class())
  if (has_ac_class(out, "sar")) {
    if (any(c("sigma", "nu") %in% names(x$dpars))) {
      stop2("SAR models are not implemented when predicting 'sigma' or 'nu'.")
    }
  }
  if (use_ac_cov(out)) {
    if (isTRUE(x$rescor)) {
      stop2("Explicit covariance terms cannot be modeled ",
            "when 'rescor' is estimated at the same time.")
    }
  }
  out
}

#' @export
tidy_acef.btl <- function(x, data = NULL, ...) {
  form <- x[["ac"]]
  if (!is.formula(form)) {
    return(empty_acef())
  }
  if (is.mixfamily(x$family)) {
    stop2("Autocorrelation terms cannot be applied in mixture models.")
  }
  px <- check_prefix(x)
  out <- data.frame(term = all_terms(form), stringsAsFactors = FALSE)
  nterms <- NROW(out)
  cnames <- c("class", "dim", "type", "time", "gr", "p", "q", "M")
  out[cnames] <- list(NA)
  out$cov <- out$nat_cov <- FALSE
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
    if (is.unstr_term(ac)) {
      out$class[i] <- "unstr"
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
  # covariance matrices of natural residuals will be handled
  # directly in the likelihood function while latent residuals will
  # be added to the linear predictor of the main parameter 'mu'
  out$nat_cov <- out$cov & has_natural_residuals(x)
  class(out) <- acef_class()
  # validate specified autocor terms
  if (any(duplicated(out$class))) {
    stop2("Can only model one term per autocorrelation class.")
  }
  if (NROW(subset2(out, dim = "time")) > 1) {
    stop2("Can only model one time-series term.")
  }
  if (NROW(subset2(out, dim = "space")) > 1) {
    stop2("Can only model one spatial term.")
  }
  if (NROW(subset2(out, nat_cov = TRUE)) > 1) {
    stop2("Can only model one covariance matrix of natural residuals.")
  }
  if (use_ac_cov(out) || has_ac_class(out, "arma")) {
    if (any(!out$dpar %in% c("", "mu") | nzchar(out$nlpar))) {
      stop2("Explicit covariance terms can only be specified on 'mu'.")
    }
  }
  out
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

# does the model need latent residuals for autocor structures?
has_ac_latent_residuals <- function(bterms) {
  !has_natural_residuals(bterms) &&
    (use_ac_cov(bterms) || has_ac_class(bterms, "arma"))
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
  colnames(M) <- rownames(M)
  not_binary <- M@x != 1
  if (any(not_binary)) {
    message("Converting all non-zero values in 'M' to 1.")
    M@x[not_binary] <- 1
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
  p <- c("ar", "ma", "sderr", "cosy", "lagsar", "errorsar",
         "car", "sdcar", "rhocar")
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

is.unstr_term <- function(x) {
  inherits(x, "unstr_term")
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
