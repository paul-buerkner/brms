#' ARMA(p,q) correlation structure
#' 
#' This functions is a constructor for the \code{cor_arma} class, representing 
#' an autoregression-moving average correlation structure of order (p, q).
#' 
#' @aliases cor.arma
#' @aliases cor_arma-class
#' 
#' @param formula A one sided formula of the form \code{~ t}, or \code{~ t | g}, 
#'   specifying a time covariate \code{t} and, optionally, 
#'   a grouping factor \code{g}. 
#'   A covariate for this correlation structure must be integer valued. 
#'   When a grouping factor is present in \code{formula}, the correlation structure 
#'   is assumed to apply only to observations within the same grouping level; 
#'   observations with different grouping levels are assumed to be uncorrelated. 
#'   Defaults to \code{~ 1}, which corresponds to using the order of the observations 
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
#' @details AR refers to autoregressive effects of residuals, which
#'   is what is typcially understood as autoregressive effects.
#'   However, one may also model autoregressive effects of the response
#'   variable, which is called ARR in \pkg{brms}.
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
  if (cov && (p > 1 || q > 1)) {
    stop2("Covariance formulation of ARMA structures is ", 
          "only possible for effects of maximal order one.")
  }
  x <- nlist(formula, p, q, r, cov = as.logical(cov))
  class(x) <- c("cor_arma", "cor_brms")
  x
}

#' @export
cor.arma <- function(formula = ~ 1, p = 0, q = 0, r = 0, cov = FALSE) {
  warn_deprecated("cor_arma")
  cor_arma(formula = formula, p = p, q = q, r = r, cov = cov)
}

#' AR(p) correlation structure
#' 
#' This function is a constructor for the \code{cor_arma} class, 
#' allowing for autoregression terms only.
#' 
#' @aliases cor.ar
#' 
#' @inheritParams cor_arma
#' 
#' @return An object of class \code{cor_arma} containing solely autoregression terms.
#' 
#' @details AR refers to autoregressive effects of residuals, which
#'   is what is typcially understood as autoregressive effects.
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

#' @export
cor.ar <- function(formula = ~ 1, p = 1, cov = FALSE) {
  warn_deprecated("cor_ar")
  cor_ar(formula = formula, p = p, cov = cov)
}
  
#' MA(q) correlation structure
#' 
#' This function is a constructor for the \code{cor_arma} class, 
#' allowing for moving average terms only.
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
  warn_deprecated("cor_ma")
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
#' @details AR refers to autoregressive effects of residuals, which
#'   is what is typcially understood as autoregressive effects.
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

#' @export
cov_fixed <- function(V) {
  warn_deprecated("cor_fixed")
  cor_fixed(V)
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
is.cor_fixed <- function(x) {
  inherits(x, "cor_fixed")
}

#' @rdname is.cor_brms
#' @export
is.cor_bsts <- function(x) {
  inherits(x, "cor_bsts")
}

#' @export
print.cor_arma <- function(x, ...) {
  cat(paste0("arma(", formula2str(x$formula), ", ", 
             get_ar(x), ", ", get_ma(x), ", ", get_arr(x),")"))
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
  if (!(is.null(x) || is.cor_brms(x) || is(x, "cor.brms"))) {
    stop2("Argument 'autocor' must be of class 'cor_brms'.")
  }
  TRUE
}

check_autocor <- function(autocor) {
  # check validity of autocor argument
  if (is.null(autocor))  {
    autocor <- cor_arma()
  }
  stop_not_cor_brms(autocor)
  autocor
}

remove_autocor <- function(x, keep = FALSE) {
  # convenience function to ignore autocorrelation terms
  if (!keep && !is.cor_bsts(x$autocor)) {
    x$autocor <- cor_arma()
  }
  x
}
