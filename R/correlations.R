#' Correlation structure classes for the \pkg{brms} package
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
#'   \item{cor_arr}{response autoregressive (ARR) structure}
#'   \item{cor_car}{Spatial conditional autoregressive (CAR) structure}
#'   \item{cor_sar}{Spatial simultaneous autoregressive (SAR) structure}
#'   \item{cor_bsts}{Bayesian structural time series (BSTS) structure}
#'   \item{cor_fixed}{fixed user-defined covariance structure}
#' }
#' 
#' @seealso 
#' \code{\link{cor_arma}, \link{cor_ar}, \link{cor_ma}, \link{cor_arr}, 
#'       \link{cor_car}, \link{cor_sar}, \link{cor_bsts}, \link{cor_fixed}}
#' 
NULL

#' ARMA(p,q) correlation structure
#' 
#' This functions is a constructor for the \code{cor_arma} class, representing 
#' an autoregression-moving average correlation structure of order (p, q).
#' 
#' @aliases cor_arma-class
#' 
#' @param formula A one sided formula of the form \code{~ t}, or \code{~ t | g}, 
#'   specifying a time covariate \code{t} and, optionally, a grouping factor \code{g}. 
#'   A covariate for this correlation structure must be integer valued. 
#'   When a grouping factor is present in \code{formula}, the correlation structure 
#'   is assumed to apply only to observations within the same grouping level; 
#'   observations with different grouping levels are assumed to be uncorrelated. 
#'   Defaults to \code{~ 1}, which corresponds to using the order of the observations 
#'   in the data as a covariate, and no groups.
#' @param p A non-negative integer specifying the autoregressive (AR) 
#'   order of the ARMA structure. Default is 0.  
#' @param q A non-negative integer specifying the moving average (MA) 
#'   order of the ARMA structure. Default is 0. 
#' @param r A non-negative integer specifying the autoregressive 
#'   response (ARR) order. See 'Details' for differences of AR and ARR 
#'   effects. Default is 0. 
#' @param cov A flag indicating whether ARMA effects should be estimated 
#'   by means of residual covariance matrices
#'   (currently only possible for stationary ARMA effects of order 1). 
#'   If \code{FALSE} (the default) a regression formulation
#'   is used that is considerably faster and allows for ARMA effects 
#'   of order higher than 1 but cannot handle user defined standard errors.
#'   
#' @return An object of class \code{cor_arma}, representing an 
#'   autoregression-moving-average correlation structure.
#' 
#' @details AR refers to autoregressive effects of residuals, which
#'   is what is typically understood as autoregressive effects.
#'   However, one may also model autoregressive effects of the response
#'   variable, which is called ARR in \pkg{brms}.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @seealso \code{\link{cor_ar}, \link{cor_ma}, \link{cor_arr}}
#' 
#' @examples
#' cor_arma(~visit|patient, p = 2, q = 2)
#' 
#' @export
cor_arma <- function(formula = ~ 1, p = 0, q = 0, r = 0, cov = FALSE) {
  formula <- as.formula(formula)
  if (!(p >= 0 && (p == round(p)))) {
    stop2("Autoregressive order must be a non-negative integer.")
  }
  if (!(q >= 0 && (q == round(q)))) {
    stop2("Moving-average order must be a non-negative integer.")
  }
  if (!(r >= 0 && (r == round(r)))) {
    stop2("Response autoregressive order must be a non-negative integer.")
  }
  if (!sum(p, q, r)) {
    stop2("At least one of 'p', 'q', and 'r' should be greater zero.")
  }
  if (cov && (p > 1 || q > 1)) {
    stop2("Covariance formulation of ARMA structures is ", 
          "only possible for effects of maximal order one.")
  }
  x <- nlist(formula, p, q, r, cov = as.logical(cov))
  class(x) <- c("cor_arma", "cor_brms")
  x
}

#' AR(p) correlation structure
#' 
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
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @seealso \code{\link{cor_arma}}
#' 
#' @examples
#' cor_ar(~visit|patient, p = 2)
#' 
#' @export
cor_ar <- function(formula = ~ 1, p = 1, cov = FALSE) {
  cor_arma(formula = formula, p = p, q = 0, r = 0, cov = cov)
}
  
#' MA(q) correlation structure
#' 
#' This function is a constructor for the \code{cor_arma} class, 
#' allowing for moving average terms only.
#' 
#' @inheritParams cor_arma
#' @param q A non-negative integer specifying the moving average (MA) 
#'   order of the ARMA structure. Default is 1.  
#' 
#' @return An object of class \code{cor_arma} containing solely moving average terms.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @seealso \code{\link{cor_arma}}
#' 
#' @examples
#' cor_ma(~visit|patient, q = 2)
#' 
#' @export
cor_ma <- function(formula = ~ 1, q = 1, cov = FALSE) {
  cor_arma(formula = formula, p = 0, q = q, r = 0, cov = cov)
}

#' ARR(r) correlation structure
#' 
#' This function is a constructor for the \code{cor_arma} class 
#' allowing for autoregressive effects of the response only.
#' 
#' @inheritParams cor_arma
#' @param r A non-negative integer specifying the autoregressive
#'   response (ARR) order of the ARMA structure. Default is 1.  
#' 
#' @return An object of class \code{cor_arma} containing solely
#'   autoregressive response terms.
#'   
#' @details AR refers to autoregressive effects of residuals, which
#'   is what is typically understood as autoregressive effects.
#'   However, one may also model autoregressive effects of the response
#'   variable, which is called ARR in \pkg{brms}.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @seealso \code{\link{cor_arma}}
#' 
#' @examples
#' cor_arr(~visit|patient, r = 2)
#' 
#' @export
cor_arr <- function(formula = ~ 1, r = 1) {
  cor_arma(formula = formula, p = 0, q = 0, r = r)
}

#' Spatial simultaneous autoregressive (SAR) structures
#' 
#' These functions are constructors for the \code{cor_sar} class
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
  W_name <- deparse(substitute(W))
  W <- sar_weights(W)
  structure(
    nlist(W, W_name, type), 
    class = c("cor_sar", "cor_brms")
  )
}

#' @rdname cor_sar
#' @export
cor_lagsar <- function(W) {
  out <- cor_sar(W, type = "lag")
  out$W_name <- deparse(substitute(W))
  out
}

#' @rdname cor_sar
#' @export
cor_errorsar <- function(W) {
  out <- cor_sar(W, type = "error")
  out$W_name <- deparse(substitute(W))
  out
}

sar_weights <- function(W) {
  # helper function to prepare spatial weights matrices
  require_package("spdep")
  if (is(W, "listw")) {
    W <- spdep::listw2mat(W)
  } else if (is(W, "nb")) {
    W <- spdep::nb2mat(W)
  } else if (!is.matrix(W)) {
    stop2("'W' must be of class 'matrix', 'listw', or 'nb'.")
  }
  W
}

#' Spatial conditional autoregressive (CAR) structures
#' 
#' These functions are constructors for the \code{cor_car} class
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
#'   are \code{"escar"} (exact sparse CAR) and \code{"esicar"}
#'   (exact sparse intrinsic CAR) and \code{"icar"} (intrinsic CAR). 
#'   More information is provided in the 'Details' section.
#' 
#' @details The \code{escar} and \code{esicar} types are 
#'   implemented based on the case study of Max Joseph
#'   (\url{https://github.com/mbjoseph/CARstan}). The \code{icar}
#'   type is implemented based on the case study of Mitzi Morris
#'   (\url{http://mc-stan.org/users/documentation/case-studies/icar_stan.html}).
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
cor_car <- function(W, formula = ~1, type = c("escar", "esicar", "icar")) {
  type <- match.arg(type)
  W_name <- deparse(substitute(W))
  W <- Matrix::Matrix(W, sparse = TRUE)
  if (!Matrix::isSymmetric(W, check.attributes = FALSE)) {
    stop2("'W' must be symmetric.")
  }
  not_binary <- !(W == 0 | W == 1)
  if (any(not_binary)) {
    message("Converting all non-zero values in 'W' to 1")
    W[not_binary] <- 1
  }
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
  out$W_name <- deparse(substitute(W))
  out
}

#' Fixed user-defined covariance matrices 
#' 
#' Define a fixed covariance matrix of the response variable
#' for instance to model multivariate effect sizes in meta-analysis.
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
  if (is.vector(V)) {
    V <- diag(V)
  } else {
    V <- as.matrix(V)
  }
  if (!isSymmetric(unname(V))) {
    stop2("'V' must be symmetric")
  }
  structure(list(V = V), class = c("cor_fixed", "cor_brms"))
}

#' Basic Bayesian Structural Time Series
#' 
#' Add a basic Bayesian structural time series component to a brms model
#' 
#' @aliases cor_bsts-class
#' 
#' @inheritParams cor_arma
#' 
#' @return An object of class \code{cor_bsts}.
#' 
#' @details Bayesian structural time series models offer an alternative 
#'   to classical AR(I)MA models (they are in fact a generalization).
#'   The basic version currently implemented in \pkg{brms} introduces
#'   local level terms for each observation, whereas each local level term 
#'   depends on the former local level term:
#'   \deqn{LL_t ~ N(LL_{t-1}, sigmaLL)}
#'   A simple introduction can be found in this blogpost: 
#'   \url{http://multithreaded.stitchfix.com/blog/2016/04/21/forget-arima/}.
#'   More complicated Bayesian structural time series models may follow
#'   in the future.
#'   
#' @examples 
#' \dontrun{
#' dat <- data.frame(y = rnorm(100), x = rnorm(100))
#' fit <- brm(y~x, data = dat, autocor = cor_bsts())
#' }   
#' 
#' @export
cor_bsts <- function(formula = ~1) {
  x <- list(formula = as.formula(formula))
  class(x) <- c("cor_bsts", "cor_brms")
  x
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

#' @rdname is.cor_brms
#' @export
is.cor_bsts <- function(x) {
  inherits(x, "cor_bsts")
}

#' @export
print.cor_empty <- function(x, ...) {
  cat("empty()")
}

#' @export
print.cor_arma <- function(x, ...) {
  cat(paste0("arma(", formula2str(x$formula), ", ", 
             get_ar(x), ", ", get_ma(x), ", ", get_arr(x),")"))
  invisible(x)
}

#' @export
print.cor_sar <- function(x, ...) {
  cat(paste0("sar(", x$W_name, ", '", x$type, "')"))
  invisible(x)
}

#' @export
print.cor_car <- function(x, ...) {
  cat(paste0(
    "car(", x$W_name, ", ", formula2str(x$formula), ", '", x$type, "')"
  ))
  invisible(x)
}

#' @export
print.cor_bsts <- function(x, ...) {
  cat(paste0("bsts(", formula2str(x$formula), ")"))
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

has_arma <- function(x) {
  # checks if any autocorrelation effects are present
  stop_not_cor_brms(x)
  isTRUE(sum(x$p, x$q, x$r) > 0)
}

get_ar <- function(x) {
  # get AR (autoregressive effects of residuals) order 
  # for which AR effects were not implemented in the present form
  stop_not_cor_brms(x)
  ifelse(is.null(x$p), 0, x$p)
}

get_ma <- function(x) {
  # get MA (moving-average) order
  stop_not_cor_brms(x)
  ifelse(is.null(x$q), 0, x$q)
}

get_arr <- function(x) {
  # get ARR (autoregressive effects of the response) order 
  # for which ARR was labled as AR
  stop_not_cor_brms(x)
  ifelse(is.null(x$r), 0, x$r)
}

use_cov <- function(x) {
  stop_not_cor_brms(x)
  if (!is.null(x$cov) && isTRUE(sum(x$p, x$q) > 0)) {
    x$cov
  } else {
    FALSE
  }
}

stop_not_cor_brms <- function(x) {
  if (!(is.null(x) || is.cor_brms(x))) {
    stop2("Argument 'autocor' must be of class 'cor_brms'.")
  }
  TRUE
}

cor_empty <- function() {
  # empty 'cor_brms' object
  structure(list(), class = c("cor_empty", "cor_brms"))
}

is.cor_empty <- function(x) {
  inherits(x, "cor_empty")
}

check_autocor <- function(autocor) {
  # check validity of autocor argument
  if (is.null(autocor))  {
    autocor <- cor_empty()
  }
  stop_not_cor_brms(autocor)
  autocor
}

remove_autocor <- function(x) {
  # convenience function to ignore autocorrelation terms
  # currently excludes ARMA, SAR, and CAR structures
  if (is_mv(x)) {
    for (r in names(x$formula$forms)) {
      ac <- x$formula$forms[[r]]$autocor
      excl_cor <- is.cor_arma(ac) || is.cor_sar(ac) || is.cor_car(ac)
      if (excl_cor) {
        x$autocor[[r]] <- x$formula$forms[[r]]$autocor <- cor_empty()
      }
    }
  } else {
    ac <- x$formula$autocor
    excl_cor <- is.cor_arma(ac) || is.cor_sar(ac) || is.cor_car(ac)
    if (excl_cor) {
      x$autocor <- x$formula$autocor <- cor_empty()
    }
  }
  x
}

subset_autocor <- function(x, subset, autocor = NULL) {
  # subset matrices stored in cor_brms objects
  # Args:
  #   x: a brmsfit object to be updated
  #   subset: indices of observations to keep
  #   autocor: optional (list of) cor_brms objects
  #     from which to take matrices
  # Returns:
  #   an updated brmsfit object
  .subset_autocor <- function(autocor) {
    if (is.cor_sar(autocor)) {
      autocor$W <- autocor$W[subset, subset, drop = FALSE]
    } else if (is.cor_fixed(autocor)) {
      autocor$V <- autocor$V[subset, subset, drop = FALSE]
    }
    return(autocor)
  }
  if (is.null(autocor)) {
    autocor <- autocor(x)
  }
  if (is_mv(x)) {
    for (i in seq_along(x$formula$forms)) {
      new_autocor <- .subset_autocor(autocor[[i]])
      x$formula$forms[[i]]$autocor <- x$autocor[[i]] <- new_autocor
    }
  } else {
    x$formula$autocor <- x$autocor <- .subset_autocor(autocor)
  }
  # prevents double updating in add_new_objects()
  structure(x, autocor_updated = TRUE)
}

regex_cor_pars <- function() {
  # regex to extract all pars of cor structures
  # used in summary.brmsfit
  p <- c("ar", "ma", "arr", "lagsar", "errorsar", "car", "sdcar", "sigmaLL")
  p <- paste0("(", p, ")", collapse = "|")
  paste0("^(", p, ")(\\[|_|$)")
}
