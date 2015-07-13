# brmsfit class
brmsfit <- function(formula = NULL, family = "", link = "", data.name = "", data = data.frame(), 
                    model = "", pars = NULL, autocor = NULL, partial = NULL, fit = NA) {
  x <- list(formula = formula, family = family, link = link, data.name = data.name,
            data = data, model = model, pars = pars, autocor = autocor, partial = partial,
            fit = fit)
  class(x) <- "brmsfit"
  return(x)
}

# brmssummary class
brmssummary <- function(formula = NULL, family = "", link = "", data.name = "", group = NULL,
                 nobs = NULL, ngrps = NULL, n.chain = 1, n.iter = 2000, n.warmup = 500, n.thin = 1,
                 sampler = "", fixed = NULL, random = list(), cor.pars = NULL, autocor = NULL, 
                 spec.pars = NULL) {
  x <- list(formula = formula, family = family, link = link, data.name = data.name, group = group, 
            nobs = nobs, ngrps = ngrps, n.chain = n.chain, n.iter = n.iter,  n.warmup = n.warmup, 
            n.thin = n.thin, sampler = sampler, fixed = fixed, random = random, 
            cor.pars = cor.pars, autocor = autocor, spec.pars = spec.pars)
  class(x) <- "brmssummary"
  x
}

#' Extract Fixed Effects for \code{brmsfit} objects
#' 
#' A generic function to extract the fixed effects from a fitted model object. 
#' 
#' @aliases fixef.brmsfit
#' 
#' @usage ## S3 method for class 'brmsfit'
#' fixef(x, estimate = "mean", ...) 
#' 
#' @param x An object of class \code{brmsfit}
#' @param estimate A character vector specifying which coefficients (e.g., "mean", "median", "sd", or "quantile") 
#' should be calculated for the fixed effects.
#' @param ... Further arguments to be passed to the functions specified in \code{estimate}
#' 
#' @return A matrix with one row per fixed effect and one column per calculated estimate.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' fixef(brm(time | cens ~ age + sex + disease, data=kidney, family="exponential"),
#'       estimate = c("mean", "sd"))
#' }
#' 
#' @export
fixef <- function(x, estimate = "mean", ...) 
  UseMethod("fixef")

#' Extract Random Effects for \code{brmsfit} objects
#' 
#' A generic function to extract the random effects of each level from a fitted model object. 
#' 
#' @aliases ranef.brmsfit
#' @usage ## S3 method for class 'brmsfit'
#' ranef(x, estimate = "mean", var = FALSE, center.zero = TRUE, ...)
#' 
#' @param x An object of a class of fitted models with random effects, typically a \code{brmsfit} object.
#' @param estimate The point estimate to be calculated for the random effects, either "mean" or "median".
#' @param var logical; indicating if the covariance matrix for each random effects should be computed.
#' @param center.zero logical; indicating if the random effects are centered around the corresponding
#'   fixed effect (if present) or around zero (the default).
#' @param ... Further arguments to be passed to the function specified in \code{estimate}
#'
#' @return A list of matrices (one per grouping factor), each with one row per level
#'  and one single column for the estimate (either mean or median).
#'     
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}   
#'   
#' @examples
#' \dontrun{
#' fit_e <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1+Trt_c|visit), 
#'              data = epilepsy, family = "poisson", n.chains = 1)
#' ## random effects means centered around zero with corresponding covariances
#' rf <- ranef(fit_e, var = TRUE)
#' attr(rf, "var")
#' ## random effects medians centered around the corresponding fixed effect
#' ranef(fit_e, estimate = "median", center.zero = FALSE)                                                        
#' }
#' 
#' @export
ranef <- function(x, estimate = "mean", var = FALSE, center.zero = TRUE, ...) 
  UseMethod("ranef")

#' Extract variance and correlation components
#' 
#' This function calculates the estimated standard deviations, correlations and covariances of the
#' random-effects terms in a mixed-effects model of class \code{brmsfit}. For linear models, the residual
#' standard deviations, correlations and covariances are also returned. 
#' 
#' @aliases VarCorr.brmsfit
#' 
#' @usage ## S3 method for class 'brmsfit'
#' VarCorr(x, estimate = "mean", as.list = TRUE, ...) 
#' 
#' @param x An object of class \code{brmsfit}.
#' @param estimate A character vector specifying which coefficients (e.g., "mean", "median", "sd", or "quantile")
#'  should be calculated for the random effects.
#' @param as.list logical; Indicates if covariance and correlation matrices should be returned as 
#'   lists of matrices (the default), or as 3-dimensional arrays.
#' @param ... Further arguments to be passed to the functions specified in \code{estimate}
#' 
#' @return A list of lists (one per grouping factor), each containing 3 elements:
#'  a matrix containing the standard deviations, a list of correlation matrices, and a list of covariance matrices.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' fit_e <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1+Trt_c|visit), 
#'              data = epilepsy, family = "poisson", n.chains = 1)
#' ## return the means of random effects covariances
#' VarCorr(fit_e)
#' ## return 2.5% and 97.5% quantiles of random effects covariances
#' VarCorr(fit_e, estimate = "quantile", probs = c(0.025, 0.975))
#' }
#' 
#' @import abind
#' @export
VarCorr <- function(x, estimate = "mean", as.list = TRUE, ...) 
  UseMethod("VarCorr")

#' Number of levels
#' 
#' Number of levels of one or more grouping factor
#' 
#' @aliases ngrps.brmsfit
#' 
#' @param object An \code{R} object typically of class \code{brmsfit}.
#' @param ... Currently ignored.
#' 
#' @details Currently there are methods for \code{brmsfit} objects.
#' @return Number(s) of levels
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @export
ngrps <- function(object, ...) 
  UseMethod("ngrps")

#' Non-linear hypothesis testing
#' 
#' Perform non-linear hypothesis testing of fixed effects parameters
#' 
#' @aliases hypothesis.brmsfit
#' 
#' @param x An \code{R} object typically of class \code{brmsfit}
#' @param hypothesis A character vector specifying one or more non-linear hypothesis concerning fixed effects
#' @param ... Currently ignored
#' 
#' @details Currently there are methods for \code{brmsfit} objects.
#' @return Summary statistics of the posterior distributions related to the hypotheses
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' fit_i <- brm(rating ~ treat + period + carry, data = inhaler, family = "cumulative")
#' 
#' hypothesis(fit_i, "treat = period + carry")
#' hypothesis(fit_i, "exp(treat) - 3 = 0")
#' 
#' ## test both of the above hypotheses with the same call 
#' hypothesis(fit_i, c("treat = period + carry", "exp(treat) - 3 = 0"))
#' }
#' 
#' @export
hypothesis <- function(x, hypothesis, ...)
  UseMethod("hypothesis")

#' Extract posterior samples
#' 
#' Extract posterior samples of specified parameters 
#' 
#' @aliases posterior.samples.brmsfit
#' 
#' @param x An \code{R} object typically of class \code{brmsfit}
#' @param parameters Name of parameters for which posterior samples should be returned, as given by a character vector or regular expressions.
#'   By default, all posterior samples of all parameters are extracted
#' @param ... Currently ignored
#'   
#' @details Currently there are methods for \code{brmsfit} objects.
#' @return A data frame containing the posterior samples, with one column per parameter.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' fit_i <- brm(rating ~ treat + period + carry + (1|subject), data = inhaler, family = "cumulative")
#' 
#' #extract posterior samples of fixed effects 
#' samples1 <- posterior.samples(fit_i, "b_")
#' head(samples1)
#' 
#' #extract posterior samples of standard deviations of random effects
#' samples2 <- posterior.samples(fit_i, "sd_")
#' head(samples2)
#' }
#' 
#' @export 
posterior.samples <- function(x, parameters = NA, ...)
  UseMethod("posterior.samples")

#' Extract parameter names
#' 
#' Extract all parameter names of a given model
#' 
#' @param x An \code{R} object typically of class \code{brmsfit}
#' @param ... Currently ignored
#' 
#' @details Currently there are methods for \code{brmsfit} objects.
#' @return A character vector containing the parameter names of the model
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @export
par.names <- function(x, ...)
  UseMethod("par.names")