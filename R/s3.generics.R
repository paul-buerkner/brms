# brmsfit class
brmsfit <- function(formula = NULL, family = "", link = "", data.name = "", data = data.frame(), 
                    model = "", exclude = NULL, prior = list(), ranef = NULL, autocor = NULL,
                    partial = NULL, fit = NA) {
  x <- list(formula = formula, family = family, link = link, data.name = data.name,
            data = data, model = model, exclude = exclude, prior = prior, 
            ranef = ranef, autocor = autocor, partial = partial, fit = fit)
  class(x) <- "brmsfit"
  return(x)
}

# brmssummary class
brmssummary <- function(formula = NULL, family = "", link = "", data.name = "", group = NULL,
                 nobs = NULL, ngrps = NULL, n.chains = 1, n.iter = 2000, n.warmup = 500, n.thin = 1,
                 sampler = "", fixed = NULL, random = list(), cor.pars = NULL, autocor = NULL, 
                 spec.pars = NULL, WAIC = "Not computed") {
  x <- list(formula = formula, family = family, link = link, data.name = data.name, group = group, 
            nobs = nobs, ngrps = ngrps, n.chains = n.chains, n.iter = n.iter,  n.warmup = n.warmup, 
            n.thin = n.thin, sampler = sampler, fixed = fixed, random = random, WAIC = WAIC,
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
#' ranef(x, estimate = "mean", var = FALSE, ...)
#' 
#' @param x An object of a class of fitted models with random effects, typically a \code{brmsfit} object.
#' @param estimate The point estimate to be calculated for the random effects, either "mean" or "median".
#' @param var logical; indicating if the covariance matrix for each random effects should be computed.
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
#' ## random effects means with corresponding covariances
#' rf <- ranef(fit_e, var = TRUE)
#' attr(rf, "var")
#' ## random effects medians
#' ranef(fit_e, estimate = "median")                                                        
#' }
#' 
#' @export
ranef <- function(x, estimate = "mean", var = FALSE, ...) 
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
#' @param hypothesis A character vector specifying one or more non-linear hypothesis concerning parameters of the model
#' @param class A string specifying the class of parameters being tested. Default is "b" for fixed effects. 
#'        Other typical options are "sd" or "cor". If \code{class = NULL}, all parameters can be tested
#'        against each other, but have to be specified with their full name (see also \code{\link[brms:par.names]{par.names}}) 
#' @param alpha the alpha-level of the tests (default is 0.05)        
#' @param ... Currently ignored
#' 
#' @details Among others, \code{hypothesis} calculates an evidence ratio for each hypothesis. 
#'   For a directed hypothesis, this is just the posterior probability under the hypothesis against its alternative.
#'   For an undirected (i.e. point) hypothesis the evidence ratio is a Bayes factor between the hypothesis and its alternative.
#'   In order to calculate this Bayes factor, all parameters related to the hypothesis must have proper priors
#'   and argument \code{sample.priors} of function \code{brm} must be set to \code{TRUE}.
#' 
#' @return Summary statistics of the posterior distributions related to the hypotheses. 
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' fit_i <- brm(rating ~ treat + period + carry + (1+treat|subject),
#'              data = inhaler, family = "gaussian", sample.prior = TRUE,
#'              prior = list(b = "normal(0,2)"), n.cluster = 2)
#' 
#' hypothesis(fit_i, "treat = period + carry")
#' hypothesis(fit_i, "exp(treat) - 3 = 0")
#' 
#' ## perform one-sided hypothesis testing
#' hypothesis(fit_i, "period + carry - 3 < 0")
#' 
#' ## compare random effects standard deviations
#' hypothesis(fit_i, "treat < Intercept", class = "sd_subject")
#' 
#' ## test the amount of random intercept variance on all variance
#' h <- paste("sd_subject_Intercept^2 / (sd_subject_Intercept^2 +",
#'            "sd_subject_treat^2 + sigma_rating^2) = 0")
#' hypothesis(fit_i, h, class = NULL)
#' 
#' ## test more than one hypothesis at once
#' hypothesis(fit_i, c("treat = period + carry", "exp(treat) - 3 = 0"))
#' }
#' 
#' @export
hypothesis <- function(x, hypothesis, class = "b", alpha = 0.05, ...)
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
#' @param add.chains A flag indicating if the returned data.frame should contain information on the chains
#' @param ... Currently ignored
#'   
#' @details Currently there are methods for \code{brmsfit} objects.
#' @return A data frame containing the posterior samples, with one column per parameter.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' fit_i <- brm(rating ~ treat + period + carry + (1|subject), 
#'              data = inhaler, family = "cumulative")
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
posterior.samples <- function(x, parameters = NA, add.chains = FALSE,...)
  UseMethod("posterior.samples")


#' Extract prior samples
#' 
#' Extract prior samples of specified parameters 
#' 
#' @aliases prior.samples.brmsfit
#' 
#' @param x An \code{R} object typically of class \code{brmsfit}
#' @param parameters Name of parameters for which posterior samples should be returned, as given by a character vector or regular expressions.
#'   By default, all prior samples are extracted
#' @param ... Currently ignored
#'   
#' @details To make use of this function, the model must contain samples of prior distributions.
#'  This can be ensured by setting \code{prior.samples = TRUE} in function \code{brm}.
#'  Currently there are methods for \code{brmsfit} objects.
#' @return A data frame containing the prior samples.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' fit_i <- brm(rating ~ treat + period + carry + (1|subject), 
#'              data = inhaler, family = "cumulative", 
#'              prior = list(b = "normal(0,2)"), sample.prior = TRUE)
#' 
#' #extract all prior samples
#' samples1 <- prior.samples(fit_i)
#' head(samples1)
#' 
#' #extract prior samples for the fixed effect of \code{treat}.
#' samples2 <- posterior.samples(fit_i, "b_treat")
#' head(samples2)
#' }
#' 
#' @export 
prior.samples <- function(x, parameters = NA, ...)
  UseMethod("prior.samples")

#' Extract parameter names
#' 
#' Extract all parameter names of a given model or formula. This help page describes the functionality for
#'  an object of class \code{brmsfit}. See \code{\link[brms:par.names.formula]{par.names.formula}} for
#'  help on the \code{formula} method.
#' 
#' @param x An \code{R} object
#' @param ... Further arguments passed to or from other methods
#' 
#' @details Currently there are methods for \code{brmsfit} and \code{formula} objects.
#' @return A character vector containing the parameter names of the model
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @export
par.names <- function(x, ...)
  UseMethod("par.names")


#' Compute the WAIC
#' 
#' Compute the Watanabe-Akaike Information Criterion based on the posterior likelihood
#' by using the \pkg{loo} package
#' 
#' @param x A fitted model object typically of class \code{brmsfit}. 
#' @param ... Optionally more fitted model objects.
#' @param compare A flag indicating if the WAICs of the models should be compared to each other
#' 
#' @details When comparing models fitted to the same data, the smaller the WAIC, the better the fit.
#' @return If just one object is provided, an object of class \code{ic}. 
#' If multiple objects are provided, an object of class \code{iclist}.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' #model with fixed effects only
#' fit_i1 <- brm(rating ~ treat + period + carry,
#'               data = inhaler, family = "gaussian")
#' WAIC(fit_i1)
#' 
#' #model with an additional random intercept for subjects
#' fit_i2 <- brm(rating ~ treat + period + carry + (1|subject),
#'              data = inhaler, family = "gaussian")
#' #compare both models
#' WAIC(fit_i1, fit_i2)                          
#' }
#' 
#' @references 
#' Vehtari, A., Gelman, A., and Gabry, J. (2015). Efficient implementation of leave-one-out cross-validation and WAIC for evaluating fitted Bayesian models.
#' 
#' Gelman, A., Hwang, J., & Vehtari, A. (2014). Understanding predictive information criteria for Bayesian models. 
#' Statistics and Computing, 24, 997-1016.
#' 
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. 
#' The Journal of Machine Learning Research, 11, 3571-3594.
#' 
#' @export
WAIC <- function(x, ..., compare = TRUE)
  UseMethod("WAIC")

#' Compute the LOO
#' 
#' Compute the Leave-one-out cross-validation based on the posterior likelihood
#' by using the \pkg{loo} package
#' 
#' @inheritParams WAIC
#' @param compare A flag indicating if the LOOs of the models should be compared to each other
#' 
#' @details When comparing models fitted to the same data, the smaller the LOO, the better the fit.
#' @return If just one object is provided, an object of class \code{ic}. 
#' If multiple objects are provided, an object of class \code{iclist}.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' #model with fixed effects only
#' fit_i1 <- brm(rating ~ treat + period + carry,
#'               data = inhaler, family = "gaussian")
#' LOO(fit_i1)
#' 
#' #model with an additional random intercept for subjects
#' fit_i2 <- brm(rating ~ treat + period + carry + (1|subject),
#'              data = inhaler, family = "gaussian")
#' #compare both models
#' LOO(fit_i1, fit_i2)                          
#' }
#' 
#' @references 
#' Vehtari, A., Gelman, A., and Gabry, J. (2015). Efficient implementation of leave-one-out cross-validation and WAIC for evaluating fitted Bayesian models.
#' 
#' Gelman, A., Hwang, J., & Vehtari, A. (2014). Understanding predictive information criteria for Bayesian models. 
#' Statistics and Computing, 24, 997-1016.
#' 
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. 
#' The Journal of Machine Learning Research, 11, 3571-3594.
#' 
#' @export
LOO <- function(x, ..., compare = TRUE)
  UseMethod("LOO")

#' Compute the pointwise log-likelihood
#' 
#' @param x A fitted model object typically of class \code{brmsfit}. 
#' @param ... Currently ignored
#' 
#' @return Usually, an S x N matrix containing the pointwise log-likelihood samples, where S is the number of samples
#'   and N is the number of observations in the data. 
#' 
#' @export
loglik <- function(x, ...)
  UseMethod("loglik")

#' Compute linear predictors term
#' 
#' @param x A fitted model object typically of class \code{brmsfit}. 
#' @param ... Currently ignored
#' 
#' @return Usually, an S x N matrix containing the linear predictor samples, where S is the number of samples
#'   and N is the number of observations in the data. 
#'   
#'  @examples
#'  \dontrun{
#'  fit_i2 <- brm(rating ~ treat + period + carry + (1|subject),
#'              data = inhaler, family = "gaussian")
#'  eta <- linear.predictor(fit_i2)                         
#'  }
#' 
#' @export
linear.predictor <- function(x, ...)
  UseMethod("linear.predictor")

#' Interface to \pkg{shinystan}
#' 
#' Provide an interface to \pkg{shinystan} for models fitted with \pkg{brms}
#' 
#' @param x A fitted model object typically of class \code{brmsfit}. 
#' @param rstudio Only relevant for RStudio users. The default (\code{rstudio=FALSE}) is to launch the app 
#' in the default web browser rather than RStudio's pop-up Viewer. Users can change the default to \code{TRUE} 
#' by setting the global option \code{options(shinystan.rstudio = TRUE)}.
#' @param ... Optional arguments to pass to \code{\link[shiny:runApp]{runApp}}
#' 
#' @return An S4 shinystan object
#' 
#' @examples
#' \dontrun{
#' fit_i2 <- brm(rating ~ treat + period + carry + (1|subject),
#'              data = inhaler, family = "gaussian")
#' launch.shiny(fit_i2)                         
#' }
#' 
#' @seealso \code{\link[shinystan:launch_shinystan]{launch_shinystan}}
#' 
#' @export
launch.shiny <- function(x, rstudio = getOption("shinystan.rstudio"), ...)
  UseMethod("launch.shiny")