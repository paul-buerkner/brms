#' ARMA(p,q) correlation structure
#' 
#' This functions is a constructor for the \code{cor.arma} class, representing an autoregression-moving average correlation structure of order (p, q).
#' 
#' @aliases cor.arma-class
#' 
#' @param formula A one sided formula of the form ~ t, or ~ t | g, specifying a time covariate t and, optionally, a grouping factor g. 
#'   A covariate for this correlation structure must be integer valued. When a grouping factor is present in \code{formula}, the correlation structure 
#'   is assumed to apply only to observations within the same grouping level; observations with different grouping levels are assumed to be uncorrelated. 
#'   Defaults to ~ 1, which corresponds to using the order of the observations in the data as a covariate, and no groups.
#' @param p A non-negative integer specifying the autoregressive order of the ARMA structure. Default is 0.  
#' @param q A non-negative integer specifying the moving average order of the ARMA structure. Default is 0. 
#' 
#' @return An object of class \code{cor.arma}, representing an autoregression-moving average correlation structure.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' cor.arma(~visit|patient, p = 2, q = 2)
#' 
#' @export
cor.arma <- function(formula = ~ 1, p = 0, q = 0) {
  if (!(p >= 0 && (p == round(p)))) {
    stop("autoregressive order must be a non-negative integer")
  }
  if (!(q >= 0 && (q == round(q)))) {
    stop("moving average order must be a non-negative integer")
  }
  x <- list(formula = formula, p = p, q = q)
  class(x) <- c("cor.arma", "cor.brms")
  x
}

#' AR(p) correlation structure
#' 
#' This function is a constructor for the \code{cor.arma} class, allowing for autoregression terms only.
#' 
#' @param formula A one sided formula of the form ~ t, or ~ t | g, specifying a time covariate t and, optionally, a grouping factor g. 
#'   A covariate for this correlation structure must be integer valued. When a grouping factor is present in \code{formula}, the correlation structure 
#'   is assumed to apply only to observations within the same grouping level; observations with different grouping levels are assumed to be uncorrelated. 
#'   Defaults to ~ 1, which corresponds to using the order of the observations in the data as a covariate, and no groups.
#' @param p A non-negative integer specifying the autoregressive order of the ARMA structure. Default is 1. 
#' 
#' @return An object of class \code{cor.arma} containing solely autoregression terms.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @seealso \code{\link{cor.arma}}
#' 
#' @examples
#' cor.ar(~visit|patient, p = 2)
#' 
#' @export
cor.ar <- function(formula = ~ 1, p = 1) {
  cor.arma(formula = formula, p = p, q = 0)
}

#' MA(q) correlation structure
#' 
#' This function is a constructor for the \code{cor.arma} class, allowing for moving average terms only.
#' 
#' @param formula A one sided formula of the form ~ t, or ~ t | g, specifying a time covariate t and, optionally, a grouping factor g. 
#'   A covariate for this correlation structure must be integer valued. When a grouping factor is present in \code{formula}, the correlation structure 
#'   is assumed to apply only to observations within the same grouping level; observations with different grouping levels are assumed to be uncorrelated. 
#'   Defaults to ~ 1, which corresponds to using the order of the observations in the data as a covariate, and no groups.
#' @param q A non-negative integer specifying the moving average order of the ARMA structure. Default is 1.
#' 
#' @return An object of class \code{cor.arma} containing solely moving average terms.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @seealso \code{\link{cor.arma}}
#' 
#' @examples
#' cor.ma(~visit|patient, q = 2)
#' 
#' @export
cor.ma <- function(formula = ~ 1, q = 1) {
  cor.arma(formula = formula, p = 0, q = q)
}

#' @export
print.cor.arma <- function(x, ...) {
  cat(paste0("arma(", gsub(" ", "", Reduce(paste, deparse(x$formula))),
             ", ",x$p,", ",x$q,")"))
}

#' Autocorrelation Function Estimation based on MCMC samples (experimental)
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
#' @return An object of class \code{link[stats:acf]{acf}}.
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
#' bacf(residuals(fit1), lag.max = 10)
#' 
#' ## fit the model again with autoregressive coefficients
#' fit2 <- brm(y ~ 1, autocor = cor.ar(p = 3))
#' summary(fit2)
#' 
#' ## investigate the residuals again (autocorrelation is gone)
#' bacf(residuals(fit2), lag.max = 10)
#' }
#'              
#' @export
bacf <- function(x, lag.max = NULL, plot = TRUE, ...) {
  series <- Reduce(paste, deparse(substitute(x)))
  if (!is.matrix(x)) {
    if (is.data.frame(x)) x <- as.matrix(x)
    else if (is.numeric(x)) x <- matrix(x, nrow = 1)
    else stop("x must be a matrix, data.frame or numeric vector")
  }
  if (is.null(colnames(x))) group <- rep(1, ncol(x))
  else group <- colnames(x)
  if (is.null(lag.max))
    lag.max <- min(sort(table(group), decreasing = TRUE)[1],
                   10*log(ncol(x)/length(unique(group)), base = 10))
  ac_names <- paste0("ac",1:lag.max) 
  lm_call <- parse(text = paste0("lm(y ~ ",paste(ac_names, collapse = "+"),", data = D)"))
  coefs <- do.call(rbind, lapply(1:nrow(x), function(i) {
    D <- as.data.frame(ar_design_matrix(x[i,], p = lag.max, group = group))
    names(D) <- paste0("ac",1:lag.max) 
    D <- cbind(D, y = x[i,])
    fit <- eval(lm_call)
    fit$coefficients[2:length(fit$coefficients)]
  }))
  out <- list(lag = array(0:lag.max, dim = c(lag.max+1, 1, 1)),
              acf = array(c(1, colMeans(coefs)), dim = c(lag.max+1, 1, 1)), 
              type = "correlation", n.used = ncol(x)/length(unique(group)), 
              series = series, snames = NULL)
  class(out) <- "acf"
  if (plot) {
    plot(out, ...)
    invisible(out)
  } else out
}