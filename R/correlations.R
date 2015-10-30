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
  x <- list(formula = formula, p = p, q = q, r = r, cov = as.logical(cov))
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
  if (!is.null(x$cov)) {
    x$cov
  } else {
    FALSE
  }
}

#' Autocorrelation Function Estimation based on MCMC-Samples (experimental)
#' 
#' Compute (and plot by default) the autocorrelation function based on MCMC samples
#' 
#' @param x A matrix, data.frame or numeric vector representing the time series. 
#'   Each row is interpreted as a sample from an MCMC procedure, whereas the columns are taken to be points of a time series.
#'   The colnames of x are taken to be the groups and it is assumed that observations can only correlate within each group.
#' @param lag.max  Maximum lag at which to calculate the autocorrelation. 
#'   Default is \code{10*log10(N/g)} where \code{N} is the number of observations and \code{g} the number of groups. 
#'   Will be automatically limited to one less than the number of maximum observations within one group.
#' @param plot logical. If \code{TRUE} (the default) results are plotted directly
#' @param ... Further arguments to be passed to \code{plot.acf}.
#'                
#' @return An object of class \code{\link[stats:acf]{acf}}.
#' 
#' @examples 
#' \dontrun{
#' ## investigate the autorcorrelation in the residuals of a fitted model
#' ## simulate data
#' set.seed(123)
#' phi <- c(0.4, 0.7, -0.3)
#' y <- 0
#' y[2] <- phi[1] * y[1] + rnorm(1)
#' y[3] <- phi[1] * y[2] + phi[2] * y[1] + rnorm(1)
#' for (i in 4:300) y[i] <- sum(phi * y[(i-1):(i-3)]) + rnorm(1)
#' 
#' ## fit the model
#' fit1 <- brm(y ~ 1)
#' summary(fit1)
#' 
#' ## investigate the residuals (autocorrelation clearly visible)
#' bacf(residuals(fit1, summary = FALSE), lag.max = 10)
#' 
#' ## fit the model again with autoregressive coefficients
#' fit2 <- brm(y ~ 1, autocor = cor_ar(p = 3))
#' summary(fit2)
#' 
#' ## investigate the residuals again (autocorrelation is gone)
#' bacf(residuals(fit2, summary = FALSE), lag.max = 10)
#' }
#'              
#' @export
macf <- function(x, lag.max = NULL, plot = TRUE, ...) {
  series <- Reduce(paste, deparse(substitute(x)))
  if (!is.matrix(x)) {
    if (is.data.frame(x)) {
      x <- as.matrix(x)
    } else if (is.numeric(x)) {
      x <- matrix(x, nrow = 1)
    } else { 
      stop("x must be a matrix, data.frame or numeric vector")
    }
  }
  if (is.null(colnames(x))) {
    group <- rep(1, ncol(x))
  } else {
    group <- colnames(x)
  }
  if (is.null(lag.max)) {
    lag.max <- min(sort(table(group), decreasing = TRUE)[1],
                   10*log(ncol(x)/length(unique(group)), base = 10))
  }
  ac_names <- paste0("ac",1:lag.max) 
  lm_call <- parse(text = paste0("lm(y ~ ",paste(ac_names, collapse = "+"),", data = D)"))
  coefs <- do.call(rbind, lapply(1:nrow(x), function(i) {
    D <- as.data.frame(arr_design_matrix(x[i, ], r = lag.max, group = group))
    names(D) <- paste0("ac", 1:lag.max) 
    D <- cbind(D, y = x[i, ])
    fit <- eval(lm_call)
    fit$coefficients[2:length(fit$coefficients)]
  }))
  out <- list(lag = array(0:lag.max, dim = c(lag.max + 1, 1, 1)),
              acf = array(c(1, colMeans(coefs)), dim = c(lag.max + 1, 1, 1)), 
              type = "correlation", 
              n.used = ncol(x) / length(unique(group)), 
              series = series, 
              snames = NULL)
  class(out) <- "acf"
  if (plot) {
    plot(out, ...)
    invisible(out)
  } else {
    out
  }
}