#' ARMA(p,q) correlation structure
#' 
#' This functions is a constructor for the \code{cor_arma} class, representing 
#' an autoregression-moving average correlation structure of order (p, q).
#' 
#' @aliases cor.arma
#' @aliases cor_arma-class
#' 
#' @param formula A one sided formula of the form ~ t, or ~ t | g, 
#'   specifying a time covariate t and, optionally, a grouping factor g. 
#'   A covariate for this correlation structure must be integer valued. 
#'   When a grouping factor is present in \code{formula}, the correlation structure 
#'   is assumed to apply only to observations within the same grouping level; 
#'   observations with different grouping levels are assumed to be uncorrelated. 
#'   Defaults to ~ 1, which corresponds to using the order of the observations 
#'   in the data as a covariate, and no groups.
#' @param p A non-negative integer specifying the autoregressive (AR) order of the ARMA structure. 
#'   Default is 0.  
#' @param q A non-negative integer specifying the moving average (MA) order of the ARMA structure. 
#'   Default is 0. 
#' @param r A non-negative integer specifying the autoregressive response (ARR) order. 
#'   See 'Details' for differences of AR and ARR effects. Default is 0. 
#' @param cov A flag indicating whether ARMA effects should be estimated by means
#'   of residual covariance matrices
#'   (currently only possible for stationary ARMA effects of order 1). 
#'   If \code{FALSE} (the default) a regression formulation
#'   is used that is considerably faster and allows for ARMA effects 
#'   of order higher than 1 but cannot handle user defined standard errors.
#'   
#' @return An object of class \code{cor_arma}, representing an 
#'   autoregression-moving-average correlation structure.
#' 
#' @details As of \pkg{brms} version 0.6.0, the AR structure refers to autoregressive effects of residuals
#'   to match the naming and implementation in other packages such as nlme. 
#'   Previously, the AR term in \pkg{brms} referred to autoregressive effects of the response.
#'   The latter are now named ARR effects and can be modeled using argument \code{r} in the
#'   \code{cor_arma} and \code{cor_arr} functions.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @seealso \code{\link[brms:cor_ar]{cor_ar}} \code{\link[brms:cor_ma]{cor_ma}}
#'   \code{\link[brms:cor_arr]{cor_arr}}
#' 
#' @examples
#' cor_arma(~visit|patient, p = 2, q = 2)
#' 
#' @export
cor_arma <- function(formula = ~ 1, p = 0, q = 0, r = 0, cov = FALSE) {
  if (!is.formula(formula)) {
    stop("formula must be of class formula")
  }
  if (!(p >= 0 && (p == round(p)))) {
    stop("autoregressive order must be a non-negative integer")
  }
  if (!(q >= 0 && (q == round(q)))) {
    stop("moving-average order must be a non-negative integer")
  }
  if (!(r >= 0 && (r == round(r)))) {
    stop("response autoregressive order must be a non-negative integer")
  }
  if (cov && (p > 1 || q > 1)) {
    stop(paste("covariance formulation of ARMA structures", 
               "is only possible for effects of maximal order 1"))
  }
  x <- nlist(formula, p, q, r, cov = as.logical(cov))
  class(x) <- c("cor_arma", "cor_brms")
  x
}

#' @export
cor.arma <- function(formula = ~ 1, p = 0, q = 0, r = 0, cov = FALSE) {
  # deprecated alias of cor_carma
  cor_arma(formula = formula, p = p, q = q, r = r, cov = cov)
}

#' AR(p) correlation structure
#' 
#' This function is a constructor for the \code{cor_arma} class, allowing for autoregression terms only.
#' 
#' @aliases cor.ar
#' 
#' @inheritParams cor_arma
#' 
#' @return An object of class \code{cor_arma} containing solely autoregression terms.
#' 
#' @details As of \pkg{brms} version 0.6.0, the AR structure refers to autoregressive effects of residuals
#'   to match the naming and implementation in other packages such as nlme. 
#'   Previously, the AR term in \pkg{brms} referred to autoregressive effects of the response.
#'   The latter are now named ARR effects and can be modeled using argument \code{r} in the
#'   \code{cor_arma} and \code{cor_arr} functions.
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

#' @export
cor.ar <- function(formula = ~ 1, p = 1, cov = FALSE) {
  # deprecated alias of cor_ar
  cor_ar(formula = formula, p = p, cov = cov)
}
  
#' MA(q) correlation structure
#' 
#' This function is a constructor for the \code{cor_arma} class, allowing for moving average terms only.
#' 
#' @aliases cor.ma
#' 
#' @inheritParams cor_arma
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

#' @export
cor.ma <- function(formula = ~ 1, q = 1, cov = FALSE) {
  # deprecated alias of cor_ma
  cor_ma(formula = formula, q = q, cov = cov)
}

#' ARR(r) correlation structure
#' 
#' This function is a constructor for the \code{cor_arma} class 
#' allowing for autoregressive effects of the response only.
#' 
#' @inheritParams cor_arma
#' 
#' @return An object of class \code{cor_arma} containing solely
#'   autoregressive response terms.
#'   
#' @details In most packages, AR effects refer to autocorrelation of residuals.
#'   \pkg{brms} also implements autocorrelation of the response, which can be specified 
#'   using argument \code{r} in the \code{cor_arma} and \code{cor_arr} functions.
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

#' @export
print.cor_arma <- function(x, ...) {
  cat(paste0("arma(", gsub(" ", "", Reduce(paste, deparse(x$formula))),
             ", ",get_ar(x),", ",get_ma(x),", ",get_arr(x),")"))
  invisible(x)
}

has_arma <- function(x) {
  # checks if any autocorrelation effects are present
  if (!is(x, "cor_brms")) {
    stop("x must be of class cor_brms")
  }
  sum(x$p, x$q, x$r) > 0
}

get_ar <- function(x) {
  # get AR (autoregressive effects of residuals) order 
  # and ensure backwards compatibility with old models (brms <= 0.5.0),
  # for which AR effects were not implemented in the present form
  if (!(is(x, "cor_arma") || is(x, "cor.arma"))) {
    stop("x must be of class cor_arma")
  }
  if (is.null(x$r)) {
    # for models fitted with brms <= 0.5.0
    0
  } else {
    x$p
  }
}

get_ma <- function(x) {
  # get MA (moving-average) order
  if (!(is(x, "cor_arma") || is(x, "cor.arma"))) {
    stop("x must be of class cor_arma")
  }
  x$q
}

get_arr <- function(x) {
  # get ARR (autoregressive effects of the response) order 
  # and ensure backwards compatibility with old models (brms <= 0.5.0),
  # for which ARR was labled as AR
  if (!(is(x, "cor_arma") || is(x, "cor.arma"))) {
    stop("x must be of class cor_arma")
  }
  if (is.null(x$r)) {
    # for models fitted with brms <= 0.5.0
    x$p
  } else {
    x$r
  }
}

use_cov <- function(x) {
  if (!(is(x, "cor_arma") || is(x, "cor.arma"))) {
    stop("x must be of class cor_arma")
  }
  if (!is.null(x$cov) && sum(x$p, x$q) > 0) {
    x$cov
  } else {
    FALSE
  }
}

check_autocor <- function(autocor) {
  # check validity of autocor argument
  if (is.null(autocor)) 
    autocor <- cor_arma()
  if (!is(autocor, "cor_brms")) { 
    stop("autocor must be of class cor_brms")
  }
  autocor
}