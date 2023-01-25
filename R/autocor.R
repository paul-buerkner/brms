# All functions in this file belong to the deprecated 'cor_brms' class
# for specifying autocorrelation structures. They will be removed in brms 3.

#' (Deprecated) Correlation structure classes for the \pkg{brms} package
#'
#' Classes of correlation structures available in the \pkg{brms} package.
#' \code{cor_brms} is not a correlation structure itself,
#' but the class common to all correlation structures implemented in \pkg{brms}.
#'
#' @name cor_brms
#' @aliases cor_brms-class
#'
#' @section Available correlation structures:
#' \describe{
#'   \item{cor_arma}{autoregressive-moving average (ARMA) structure,
#'   with arbitrary orders for the autoregressive and moving
#'   average components}
#'   \item{cor_ar}{autoregressive (AR) structure of arbitrary order}
#'   \item{cor_ma}{moving average (MA) structure of arbitrary order}
#'   \item{cor_car}{Spatial conditional autoregressive (CAR) structure}
#'   \item{cor_sar}{Spatial simultaneous autoregressive (SAR) structure}
#'   \item{cor_fixed}{fixed user-defined covariance structure}
#' }
#'
#' @seealso
#' \code{\link{cor_arma}, \link{cor_ar}, \link{cor_ma},
#'       \link{cor_car}, \link{cor_sar}, \link{cor_fixed}}
#'
NULL

#' (Deprecated) ARMA(p,q) correlation structure
#'
#' This function is deprecated. Please see \code{\link{arma}} for the new syntax.
#' This functions is a constructor for the \code{cor_arma} class, representing
#' an autoregression-moving average correlation structure of order (p, q).
#'
#' @aliases cor_arma-class
#'
#' @param formula A one sided formula of the form \code{~ t}, or \code{~ t | g},
#'   specifying a time covariate \code{t} and, optionally, a grouping factor
#'   \code{g}. A covariate for this correlation structure must be integer
#'   valued. When a grouping factor is present in \code{formula}, the
#'   correlation structure is assumed to apply only to observations within the
#'   same grouping level; observations with different grouping levels are
#'   assumed to be uncorrelated. Defaults to \code{~ 1}, which corresponds to
#'   using the order of the observations in the data as a covariate, and no
#'   groups.
#' @param p A non-negative integer specifying the autoregressive (AR)
#'   order of the ARMA structure. Default is 0.
#' @param q A non-negative integer specifying the moving average (MA)
#'   order of the ARMA structure. Default is 0.
#' @param r No longer supported.
#' @param cov A flag indicating whether ARMA effects should be estimated by
#'   means of residual covariance matrices. This is currently only possible for
#'   stationary ARMA effects of order 1. If the model family does not have
#'   natural residuals, latent residuals are added automatically. If
#'   \code{FALSE} (the default) a regression formulation is used that is
#'   considerably faster and allows for ARMA effects of order higher than 1 but
#'   is only available for \code{gaussian} models and some of its
#'   generalizations.
#'
#' @return An object of class \code{cor_arma}, representing an
#'   autoregression-moving-average correlation structure.
#'
#' @seealso \code{\link{cor_ar}}, \code{\link{cor_ma}}
#'
#' @examples
#' cor_arma(~ visit | patient, p = 2, q = 2)
#'
#' @export
cor_arma <- function(formula = ~1, p = 0, q = 0, r = 0, cov = FALSE) {
  formula <- as.formula(formula)
  p <- as_one_numeric(p)
  q <- as_one_numeric(q)
  cov <- as_one_logical(cov)
  if ("r" %in% names(match.call())) {
    warning2("The ARR structure is no longer supported and ignored.")
  }
  if (!(p >= 0 && p == round(p))) {
    stop2("Autoregressive order must be a non-negative integer.")
  }
  if (!(q >= 0 && q == round(q))) {
    stop2("Moving-average order must be a non-negative integer.")
  }
  if (!sum(p, q)) {
    stop2("At least one of 'p' and 'q' should be greater zero.")
  }
  if (cov && (p > 1 || q > 1)) {
    stop2("Covariance formulation of ARMA structures is ",
          "only possible for effects of maximal order one.")
  }
  x <- nlist(formula, p, q, cov)
  class(x) <- c("cor_arma", "cor_brms")
  x
}

#' (Deprecated) AR(p) correlation structure
#'
#' This function is deprecated. Please see \code{\link{ar}} for the new syntax.
#' This function is a constructor for the \code{cor_arma} class,
#' allowing for autoregression terms only.
#'
#' @inheritParams cor_arma
#' @param p A non-negative integer specifying the autoregressive (AR)
#'   order of the ARMA structure. Default is 1.
#'
#' @return An object of class \code{cor_arma} containing solely autoregression terms.
#'
#' @details AR refers to autoregressive effects of residuals, which
#'   is what is typically understood as autoregressive effects.
#'   However, one may also model autoregressive effects of the response
#'   variable, which is called ARR in \pkg{brms}.
#'
#' @seealso \code{\link{cor_arma}}
#'
#' @examples
#' cor_ar(~visit|patient, p = 2)
#'
#' @export
cor_ar <- function(formula = ~1, p = 1, cov = FALSE) {
  cor_arma(formula = formula, p = p, q = 0, cov = cov)
}

#' (Deprecated) MA(q) correlation structure
#'
#' This function is deprecated. Please see \code{\link{ma}} for the new syntax.
#' This function is a constructor for the \code{cor_arma} class,
#' allowing for moving average terms only.
#'
#' @inheritParams cor_arma
#' @param q A non-negative integer specifying the moving average (MA)
#'   order of the ARMA structure. Default is 1.
#'
#' @return An object of class \code{cor_arma} containing solely moving
#' average terms.
#'
#' @seealso \code{\link{cor_arma}}
#'
#' @examples
#' cor_ma(~visit|patient, q = 2)
#'
#' @export
cor_ma <- function(formula = ~1, q = 1, cov = FALSE) {
  cor_arma(formula = formula, p = 0, q = q, cov = cov)
}

#' (Defunct) ARR correlation structure
#'
#' The ARR correlation structure is no longer supported.
#'
#' @inheritParams cor_arma
#'
#' @keywords internal
#' @export
cor_arr <- function(formula = ~1, r = 1) {
  cor_arma(formula = formula, p = 0, q = 0, r = r)
}

#' (Deprecated) Compound Symmetry (COSY) Correlation Structure
#'
#' This function is deprecated. Please see \code{\link{cosy}} for the new syntax.
#' This functions is a constructor for the \code{cor_cosy} class, representing
#' a compound symmetry structure corresponding to uniform correlation.
#'
#' @aliases cor_cosy-class
#'
#' @inheritParams cor_arma
#'
#' @return An object of class \code{cor_cosy}, representing a compound symmetry
#'   correlation structure.
#'
#' @examples
#' cor_cosy(~ visit | patient)
#'
#' @export
cor_cosy <- function(formula = ~1) {
  formula <- as.formula(formula)
  x <- nlist(formula)
  class(x) <- c("cor_cosy", "cor_brms")
  x
}

#' (Deprecated) Spatial simultaneous autoregressive (SAR) structures
#'
#' Thse functions are deprecated. Please see \code{\link{sar}} for the new
#' syntax. These functions are constructors for the \code{cor_sar} class
#' implementing spatial simultaneous autoregressive structures.
#' The \code{lagsar} structure implements SAR of the response values:
#' \deqn{y = \rho W y + \eta + e}
#' The \code{errorsar} structure implements SAR of the residuals:
#' \deqn{y = \eta + u, u = \rho W u + e}
#' In the above equations, \eqn{\eta} is the predictor term and
#' \eqn{e} are independent normally or t-distributed residuals.
#'
#' @param W An object specifying the spatial weighting matrix.
#'   Can be either the spatial weight matrix itself or an
#'   object of class \code{listw} or \code{nb}, from which
#'   the spatial weighting matrix can be computed.
#' @param type Type of the SAR structure. Either \code{"lag"}
#'   (for SAR of the response values) or \code{"error"}
#'   (for SAR of the residuals).
#'
#' @details Currently, only families \code{gaussian} and \code{student}
#'   support SAR structures.
#'
#' @return An object of class \code{cor_sar} to be used in calls to
#'   \code{\link{brm}}.
#'
#' @examples
#' \dontrun{
#' data(oldcol, package = "spdep")
#' fit1 <- brm(CRIME ~ INC + HOVAL, data = COL.OLD,
#'             autocor = cor_lagsar(COL.nb),
#'             chains = 2, cores = 2)
#' summary(fit1)
#' plot(fit1)
#'
#' fit2 <- brm(CRIME ~ INC + HOVAL, data = COL.OLD,
#'             autocor = cor_errorsar(COL.nb),
#'             chains = 2, cores = 2)
#' summary(fit2)
#' plot(fit2)
#' }
#'
#' @export
cor_sar <- function(W, type = c("lag", "error")) {
  type <- match.arg(type)
  W_name <- deparse0(substitute(W))
  W <- validate_sar_matrix(W)
  structure(
    nlist(W, W_name, type),
    class = c("cor_sar", "cor_brms")
  )
}

#' @rdname cor_sar
#' @export
cor_lagsar <- function(W) {
  out <- cor_sar(W, type = "lag")
  out$W_name <- deparse0(substitute(W))
  out
}

#' @rdname cor_sar
#' @export
cor_errorsar <- function(W) {
  out <- cor_sar(W, type = "error")
  out$W_name <- deparse0(substitute(W))
  out
}

#' (Deprecated) Spatial conditional autoregressive (CAR) structures
#'
#' These function are deprecated. Please see \code{\link{car}} for the new
#' syntax. These functions are constructors for the \code{cor_car} class
#' implementing spatial conditional autoregressive structures.
#'
#' @param W Adjacency matrix of locations.
#'   All non-zero entries are treated as if the two locations
#'   are adjacent. If \code{formula} contains a grouping factor,
#'   the row names of \code{W} have to match the levels
#'   of the grouping factor.
#' @param formula An optional one-sided formula of the form
#'   \code{~ 1 | g}, where \code{g} is a grouping factor mapping
#'   observations to spatial locations. If not specified,
#'   each observation is treated as a separate location.
#'   It is recommended to always specify a grouping factor
#'   to allow for handling of new data in post-processing methods.
#' @param type Type of the CAR structure. Currently implemented
#'   are \code{"escar"} (exact sparse CAR), \code{"esicar"}
#'   (exact sparse intrinsic CAR), \code{"icar"} (intrinsic CAR),
#'   and \code{"bym2"}. More information is provided in the 'Details' section.
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
#' fit <- brm(y | trials(size) ~ x1 + x2, data = dat,
#'            family = binomial(), autocor = cor_car(W))
#' summary(fit)
#' }
#'
#' @export
cor_car <- function(W, formula = ~1, type = "escar") {
  options <- c("escar", "esicar", "icar", "bym2")
  type <- match.arg(type, options)
  W_name <- deparse0(substitute(W))
  W <- validate_car_matrix(W)
  formula <- as.formula(formula)
  if (!is.null(lhs(formula))) {
    stop2("'formula' should be a one-sided formula.")
  }
  if (length(attr(terms(formula), "term.labels")) > 1L) {
    stop2("'formula' should not contain more than one term.")
  }
  structure(
    nlist(W, W_name, formula, type),
    class = c("cor_car", "cor_brms")
  )
}

#' @rdname cor_car
#' @export
cor_icar <- function(W, formula = ~1) {
  out <- cor_car(W, formula, type = "icar")
  out$W_name <- deparse0(substitute(W))
  out
}

#' (Deprecated) Fixed user-defined covariance matrices
#'
#' This function is deprecated. Please see \code{\link{fcor}} for the new
#' syntax. Define a fixed covariance matrix of the response variable for
#' instance to model multivariate effect sizes in meta-analysis.
#'
#' @aliases cov_fixed
#'
#' @param V Known covariance matrix of the response variable.
#'   If a vector is passed, it will be used as diagonal entries
#'   (variances) and covariances will be set to zero.
#'
#' @return An object of class \code{cor_fixed}.
#'
#' @examples
#' \dontrun{
#' dat <- data.frame(y = rnorm(3))
#' V <- cbind(c(0.5, 0.3, 0.2), c(0.3, 1, 0.1), c(0.2, 0.1, 0.2))
#' fit <- brm(y~1, data = dat, autocor = cor_fixed(V))
#' }
#'
#' @export
cor_fixed <- function(V) {
  V_name <- deparse0(substitute(V))
  if (is.vector(V)) {
    V <- diag(V)
  } else {
    V <- as.matrix(V)
  }
  if (!isSymmetric(unname(V))) {
    stop2("'V' must be symmetric")
  }
  structure(nlist(V, V_name), class = c("cor_fixed", "cor_brms"))
}

#' (Defunct) Basic Bayesian Structural Time Series
#'
#' The BSTS correlation structure is no longer supported.
#'
#' @inheritParams cor_arma
#'
#' @keywords internal
#' @export
cor_bsts <- function(formula = ~1) {
  stop2("The BSTS structure is no longer supported.")
}

#' Check if argument is a correlation structure
#'
#' Check if argument is one of the correlation structures
#' used in \pkg{brms}.
#'
#' @param x An \R object.
#'
#' @export
is.cor_brms <- function(x) {
  inherits(x, "cor_brms")
}

#' @rdname is.cor_brms
#' @export
is.cor_arma <- function(x) {
  inherits(x, "cor_arma")
}

#' @rdname is.cor_brms
#' @export
is.cor_cosy <- function(x) {
  inherits(x, "cor_cosy")
}

#' @rdname is.cor_brms
#' @export
is.cor_sar <- function(x) {
  inherits(x, "cor_sar")
}

#' @rdname is.cor_brms
#' @export
is.cor_car <- function(x) {
  inherits(x, "cor_car")
}

#' @rdname is.cor_brms
#' @export
is.cor_fixed <- function(x) {
  inherits(x, "cor_fixed")
}

#' @export
print.cor_empty <- function(x, ...) {
  cat("empty()\n")
}

#' @export
print.cor_arma <- function(x, ...) {
  cat(paste0("arma(", formula2str(x$formula), ", ", x$p, ", ", x$q, ")\n"))
  invisible(x)
}

#' @export
print.cor_cosy <- function(x, ...) {
  cat(paste0("cosy(", formula2str(x$formula), ")\n"))
  invisible(x)
}

#' @export
print.cor_sar <- function(x, ...) {
  cat(paste0("sar(", x$W_name, ", '", x$type, "')\n"))
  invisible(x)
}

#' @export
print.cor_car <- function(x, ...) {
  form <- formula2str(x$formula)
  cat(paste0("car(", x$W_name, ", ", form, ", '", x$type, "')\n"))
  invisible(x)
}

#' @export
print.cor_fixed <- function(x, ...) {
  cat("Fixed covariance matrix: \n")
  print(x$V)
  invisible(x)
}

#' @export
print.cov_fixed <- function(x, ...) {
  class(x) <- "cor_fixed"
  print.cor_fixed(x)
}

stop_not_cor_brms <- function(x) {
  if (!(is.null(x) || is.cor_brms(x))) {
    stop2("Argument 'autocor' must be of class 'cor_brms'.")
  }
  TRUE
}

# empty 'cor_brms' object
cor_empty <- function() {
  structure(list(), class = c("cor_empty", "cor_brms"))
}

is.cor_empty <- function(x) {
  inherits(x, "cor_empty")
}

#' (Deprecated) Extract Autocorrelation Objects
#'
#' @inheritParams posterior_predict.brmsfit
#' @param ... Currently unused.
#'
#' @return A \code{cor_brms} object or a list of such objects for multivariate
#'   models. Not supported for models fitted with brms 2.11.1 or higher.
#'
#' @export
autocor.brmsfit <- function(object, resp = NULL, ...) {
  warning2("Method 'autocor' is deprecated and will be removed in the future.")
  object <- restructure(object)
  resp <- validate_resp(resp, object)
  if (!is.null(resp)) {
    # multivariate model
    autocor <- object$autocor[resp]
    if (length(resp) == 1L) {
      autocor <- autocor[[1]]
    }
  } else {
    # univariate model
    autocor <- object$autocor
  }
  autocor
}

#' @rdname autocor.brmsfit
#' @export
autocor <- function(object, ...) {
  UseMethod("autocor")
}

# extract variables for autocorrelation structures
# @param autocor object of class 'cor_brms'
# @return a list with elements 'time', and 'group'
terms_autocor <- function(autocor) {
  out <- list()
  formula <- autocor$formula
  if (is.null(formula)) {
    formula <- ~1
  }
  if (!is.null(lhs(formula))) {
    stop2("Autocorrelation formulas must be one-sided.")
  }
  formula <- formula2str(formula)
  time <- as.formula(paste("~", gsub("~|\\|[[:print:]]*", "", formula)))
  time_vars <- all_vars(time)
  if (is.cor_car(autocor) && length(time_vars) > 0L) {
    stop2("The CAR structure should not contain a 'time' variable.")
  }
  if (length(time_vars) > 1L) {
    stop2("Autocorrelation structures may only contain 1 time variable.")
  }
  if (length(time_vars)) {
    out$time <- time_vars
  } else {
    out$time <- NA
  }
  group <- sub("^\\|*", "", sub("~[^\\|]*", "", formula))
  stopif_illegal_group(group)
  group_vars <- all_vars(group)
  if (length(group_vars)) {
    out$group <- paste0(group_vars, collapse = ":")
  } else {
    out$group <- NA
  }
  out
}

# transform a 'cor_brms' object into a formula
# this ensure compatibility with brms <= 2.11
as_formula_cor_brms <- function(x) {
  stop_not_cor_brms(x)
  if (is.cor_empty(x))  {
    return(NULL)
  }
  args <- data2 <- list()
  pac <- terms_autocor(x)
  if (is.cor_arma(x)) {
    fun <- "arma"
    args$time <- pac$time
    args$gr <- pac$group
    args$p <- x$p
    args$q <- x$q
    args$cov <- x$cov
    out <- paste0(names(args), " = ", args, collapse = ", ")
    out <- paste0("arma(", out, ")")
  } else if (is.cor_cosy(x)) {
    fun <- "cosy"
    args$time <- pac$time
    args$gr <- pac$group
  } else if (is.cor_sar(x)) {
    fun <- "sar"
    args$M <- make_M_names(x$W_name)
    args$type <- paste0("'", x$type, "'")
    data2[[args$M]] <- x$W
  } else if (is.cor_car(x)) {
    fun <- "car"
    args$M <- make_M_names(x$W_name)
    args$gr <- pac$group
    args$type <- paste0("'", x$type, "'")
    data2[[args$M]] <- x$W
  } else if (is.cor_fixed(x)) {
    fun <- "fcor"
    args$M <- make_M_names(x$V_name)
    data2[[args$M]] <- x$V
  }
  out <- paste0(names(args), " = ", args, collapse = ", ")
  out <- paste0(fun, "(", out, ")")
  out <- str2formula(out)
  attr(out, "data2") <- data2
  class(out) <- c("cor_brms_formula", "formula")
  out
}

# ensures covariance matrix inputs are named reasonably
make_M_names <- function(x) {
  out <- make.names(x)
  if (!length(out)) {
    # likely unique random name for the matrix argument
    out <- paste0("M", collapse(sample(0:9, 5, TRUE)))
  }
  out
}

# get data objects from 'autocor' for use in 'data2'
# for backwards compatibility with brms <= 2.11
get_data2_autocor <- function(x, ...) {
  UseMethod("get_data2_autocor")
}

#' @export
get_data2_autocor.brmsformula <- function(x, ...) {
  attr(attr(x$formula, "autocor"), "data2")
}

#' @export
get_data2_autocor.mvbrmsformula <- function(x, ...) {
  ulapply(x$forms, get_data2_autocor, recursive = FALSE)
}

#' @export
print.cor_brms_formula <- function(x, ...) {
  y <- x
  attr(y, "data2") <- NULL
  class(y) <- "formula"
  print(y)
  invisible(x)
}
