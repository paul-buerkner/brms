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

#'@export
print.cor.arma <- function(x, ...) {
  cat(paste0("arma(", gsub(" ", "", Reduce(paste, deparse(x$formula))),
             ", ",x$p,", ",x$q,")"))
}