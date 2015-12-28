brmsfit <- function(formula = NULL, family = "", link = "", data.name = "", 
                    data = data.frame(), model = "", exclude = NULL,
                    prior = list(), ranef = NULL, autocor = NULL,
                    multiply = NULL, partial = NULL, 
                    cov.ranef = NULL, fit = NA) {
  # brmsfit class
  x <- nlist(formula, family, link, data.name, data, model, exclude, prior, 
             ranef, autocor, multiply, partial, cov.ranef, fit)
  class(x) <- "brmsfit"
  x
}

brmssummary <- function(formula = NULL, family = "", link = "", 
                        data.name = "", group = NULL, nobs = NULL, 
                        ngrps = NULL, n.chains = 1, n.iter = 2000, 
                        n.warmup = 500, n.thin = 1, sampler = "", 
                        autocor = NULL, multiply = NULL,
                        fixed = NULL, random = list(), 
                        cor_pars = NULL, spec_pars = NULL, 
                        mult_pars = NULL, WAIC = "Not computed") {
  # brmssummary class
  x <- nlist(formula, family, link, data.name, group, nobs, ngrps, n.chains, 
             n.iter,  n.warmup, n.thin, sampler, autocor, multiply, 
             fixed, random, cor_pars, spec_pars, mult_pars, WAIC)
  class(x) <- "brmssummary"
  x
}

#' Extract Fixed Effects Estimates
#' 
#' A generic function to extract the fixed effects from a fitted model object. 
#' 
#' @aliases fixef.brmsfit
#' 
#' @param x An object of class \code{brmsfit}
#' @param estimate A character vector specifying which coefficients 
#'  (e.g., "mean", "median", "sd", or "quantile") 
#'  should be calculated for the fixed effects.
#' @param ... Further arguments to be passed to the functions 
#'  specified in \code{estimate}
#' 
#' @return A matrix with one row per fixed effect 
#'   and one column per calculated estimate.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' fit <- brm(time | cens(censored) ~ age + sex + disease, 
#'            data = kidney, family = "exponential")
#' fixef(fit, estimate = c("mean", "sd"))
#' }
#' 
#' @export
fixef <- function(x, ...) 
  UseMethod("fixef")

#' Extract Random Effects for \code{brmsfit} objects
#' 
#' A generic function to extract the random effects 
#' of each level from a fitted model object. 
#' 
#' @aliases ranef.brmsfit
#' 
#' @param x An object of a class of fitted models with random effects, 
#'  typically a \code{brmsfit} object.
#' @param estimate The point estimate to be calculated 
#'  for the random effects, either "mean" or "median".
#' @param var logical; indicating if the covariance matrix 
#'  for each random effects should be computed.
#' @param ... Further arguments to be passed to the function 
#'  specified in \code{estimate}
#'
#' @return A list of matrices (one per grouping factor), 
#'  each with one row per level and one single column 
#'  for the estimate (either mean or median).
#'     
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}   
#'   
#' @examples
#' \dontrun{
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1+Trt_c|visit), 
#'              data = epilepsy, family = "poisson", n.chains = 1)
#' ## random effects means with corresponding covariances
#' rf <- ranef(fit, var = TRUE)
#' attr(rf, "var")
#' ## random effects medians
#' ranef(fit, estimate = "median")                                                        
#' }
#' 
#' @export
ranef <- function(x, ...) 
  UseMethod("ranef")

#' Extract variance and correlation components
#' 
#' This function calculates the estimated standard deviations, 
#' correlations and covariances of therandom-effects terms 
#' in a mixed-effects model of class \code{brmsfit}. 
#' For linear models, the residual standard deviations, 
#' correlations and covariances are also returned. 
#' 
#' @aliases VarCorr.brmsfit
#' 
#' @param x A fitted model object usually of class \code{brmsift}
#' @param estimate A character vector specifying which coefficients 
#'  (e.g., "mean", "median", "sd", or "quantile")
#'  should be calculated for the random effects.
#' @param as.list logical; Indicates if covariance 
#'  and correlation matrices should be returned as 
#'  lists of matrices (the default), or as 3-dimensional arrays.
#' @param ... Further arguments to be passed to the functions 
#'  specified in \code{estimate}
#' 
#' @return An object of class \code{VarCorr_brmsfit}, 
#' which is a list of lists (one per grouping factor), 
#' each containing 3 elements: a matrix containing the standard deviations, 
#' a list of correlation matrices, and a list of covariance matrices. 
#' Can be coerced to a \code{data.frame} 
#' by using the \code{as.data.frame} method.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1+Trt_c|visit), 
#'              data = epilepsy, family = "poisson", n.chains = 1)
#' ## return the means of random effects covariances
#' (vc <- VarCorr(fit))
#' as.data.frame(vc)
#' 
#' ## return 2.5% and 97.5% quantiles of random effects covariances
#' VarCorr(fit, estimate = "quantile", probs = c(0.025, 0.975))
#' }
#' 
#' @export
VarCorr <- function(x, ...) 
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
#' @param hypothesis A character vector specifying one or more 
#'  non-linear hypothesis concerning parameters of the model
#' @param class A string specifying the class of parameters being tested. 
#'  Default is "b" for fixed effects. 
#'  Other typical options are "sd" or "cor". 
#'  If \code{class = NULL}, all parameters can be tested
#'  against each other, but have to be specified with their full name 
#'  (see also \code{\link[brms:parnames]{parnames}}) 
#' @param group Name of a grouping factor to evaluate only 
#'  random effects parameters related to this grouping factor.
#'  Ignored if \code{class} is not \code{"sd"} or \code{"cor"}.
#' @param alpha the alpha-level of the tests (default is 0.05)        
#' @param ... Currently ignored
#' 
#' @details Among others, \code{hypothesis} computes an 
#'  evidence ratio for each hypothesis. 
#'  For a directed hypothesis, this is just the posterior probability 
#'  under the hypothesis against its alternative.
#'  For an undirected (i.e. point) hypothesis the evidence ratio 
#'  is a Bayes factor between the hypothesis and its alternative.
#'  In order to calculate this Bayes factor, all parameters related 
#'  to the hypothesis must have proper priors
#'  and argument \code{sample.prior} of function \code{brm} 
#'  must be set to \code{TRUE}. 
#'  When interpreting Bayes factors, make sure 
#'  that your priors are reasonable and carefully chosen,
#'  as the result will depend heavily on the priors. 
#'  It particular, avoid using default priors.
#' 
#' @return Summary statistics of the posterior distributions 
#'  related to the hypotheses. 
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' fit <- brm(rating ~ treat + period + carry + (1+treat|subject),
#'              data = inhaler, family = "gaussian", sample.prior = TRUE,
#'              prior = set_prior("normal(0,2)", class = "b"), n.cluster = 2)
#' 
#' hypothesis(fit, "treat = period + carry")
#' hypothesis(fit, "exp(treat) - 3 = 0")
#' 
#' ## perform one-sided hypothesis testing
#' hypothesis(fit, "period + carry - 3 < 0")
#' 
#' ## compare random effects standard deviations
#' hypothesis(fit, "treat < Intercept", class = "sd", group  = "subject")
#' 
#' ## test the amount of random intercept variance on all variance
#' h <- paste("sd_subject_Intercept^2 / (sd_subject_Intercept^2 +",
#'            "sd_subject_treat^2 + sigma_rating^2) = 0")
#' hypothesis(fit, h, class = NULL)
#' 
#' ## test more than one hypothesis at once
#' hypothesis(fit, c("treat = period + carry", "exp(treat) - 3 = 0"))
#' }
#' 
#' @export
hypothesis <- function(x, hypothesis, ...)
  UseMethod("hypothesis")

#' Extract posterior samples
#' 
#' Extract posterior samples of specified parameters 
#' 
#' @aliases posterior.samples posterior_samples.brmsfit posterior.samples.brmsfit
#' 
#' @param x An \code{R} object typically of class \code{brmsfit}
#' @param pars Names of parameters for which posterior samples 
#'   should be returned, as given by a character vector or regular expressions.
#'   By default, all posterior samples of all parameters are extracted
#' @param parameters A deprecated alias of \code{pars}   
#' @param exact_match Indicates whether parameter names 
#'   should be matched exactly or treated as regular expression. 
#'   Default is \code{FALSE}.
#' @param add_chains A flag indicating if the returned data.frame 
#'   should contain information on the chains
#' @param subset A numeric vector indicating the rows 
#'        (i.e., posterior samples) to be returned. 
#'        If \code{NULL} (the default), all  posterior samples are returned.
#' @param as.matrix Should the output be a \code{matrix} 
#'   instead of a \code{data.frame}? Defaults to \code{FALSE}
#' @param ... additional arguments
#'   
#' @details Currently there are methods for \code{brmsfit} objects.
#' @return A data frame containing the posterior samples, 
#'   with one column per parameter.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' fit <- brm(rating ~ treat + period + carry + (1|subject), 
#'            data = inhaler, family = "cumulative")
#' 
#' #extract posterior samples of fixed effects 
#' samples1 <- posterior_samples(fit, "^b")
#' head(samples1)
#' 
#' #extract posterior samples of standard deviations of random effects
#' samples2 <- posterior_samples(fit, "^sd")
#' head(samples2)
#' }
#' 
#' @export 
posterior_samples <- function(x, pars = NA, ...)
  UseMethod("posterior_samples")

# deprecated alias of posterior_samples
#' @export 
posterior.samples <- function(x, pars = NA, ...)
  UseMethod("posterior_samples")

#' Extract prior samples
#' 
#' Extract prior samples of specified parameters 
#' 
#' @aliases prior_samples.brmsfit
#' 
#' @param x An \code{R} object typically of class \code{brmsfit}
#' @param pars Names of parameters for which prior samples should be returned, 
#'   as given by a character vector or regular expressions.
#'   By default, all prior samples are extracted
#' @param parameters A deprecated alias of \code{pars}       
#' @param ... Currently ignored
#'   
#' @details To make use of this function, 
#'  the model must contain samples of prior distributions.
#'  This can be ensured by setting \code{sample.prior = TRUE} 
#'  in function \code{brm}.
#'  Currently there are methods for \code{brmsfit} objects.
#' @return A data frame containing the prior samples.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' fit <- brm(rating ~ treat + period + carry + (1|subject), 
#'            data = inhaler, family = "cumulative", 
#'            prior = set_prior("normal(0,2)", class = "b"), 
#'            sample.prior = TRUE)
#' 
#' #extract all prior samples
#' samples1 <- prior_samples(fit)
#' head(samples1)
#' 
#' #extract prior samples for the fixed effect of \code{treat}.
#' samples2 <- posterior_samples(fit, "b_treat")
#' head(samples2)
#' }
#' 
#' @export 
prior_samples <- function(x, pars = NA, ...)
  UseMethod("prior_samples")

#' Extract Parameter Names
#' 
#' Extract all parameter names of a given model.
#'  
#' @aliases par.names parnames.brmsfit par.names.brmsfit
#' 
#' @param x An \code{R} object
#' @param ... Further arguments passed to or from other methods
#' 
#' @details Currently there are methods for \code{brmsfit} 
#'   and \code{formula} objects.
#' @return A character vector containing the parameter names of the model.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @export
parnames <- function(x, ...)
  UseMethod("parnames")

# deprecated alias of parnames
#' @export
par.names <- function(x, ...)
  UseMethod("parnames")

#' Compute the WAIC
#' 
#' Compute the Watanabe-Akaike Information Criterion 
#' based on the posterior likelihood by using the \pkg{loo} package
#' 
#' @aliases WAIC.brmsfit
#' 
#' @param x A fitted model object typically of class \code{brmsfit}. 
#' @param ... Optionally more fitted model objects.
#' @param compare A flag indicating if the WAICs 
#'  of the models should be compared to each other.
#' 
#' @details When comparing models fitted to the same data, 
#'  the smaller the WAIC, the better the fit.
#' @return If just one object is provided, an object of class \code{ic}. 
#'  If multiple objects are provided, an object of class \code{iclist}.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' #model with fixed effects only
#' fit1 <- brm(rating ~ treat + period + carry,
#'             data = inhaler, family = "gaussian")
#' WAIC(fit1)
#' 
#' #model with an additional random intercept for subjects
#' fit2 <- brm(rating ~ treat + period + carry + (1|subject),
#'             data = inhaler, family = "gaussian")
#' #compare both models
#' WAIC(fit1, fit2)                          
#' }
#' 
#' @references 
#' Vehtari, A., Gelman, A., and Gabry, J. (2015). 
#' Efficient implementation of leave-one-out cross-validation 
#' and WAIC for evaluating fitted Bayesian models.
#' 
#' Gelman, A., Hwang, J., & Vehtari, A. (2014). 
#' Understanding predictive information criteria for Bayesian models. 
#' Statistics and Computing, 24, 997-1016.
#' 
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation 
#' and widely applicable information criterion in singular learning theory. 
#' The Journal of Machine Learning Research, 11, 3571-3594.
#' 
#' @export
WAIC <- function(x, ..., compare = TRUE)
  UseMethod("WAIC")

#' Compute LOO
#' 
#' Compute Leave-one-out cross-validation based on the posterior likelihood
#' by using the \pkg{loo} package
#' 
#' @aliases LOO.brmsfit
#' 
#' @inheritParams WAIC
#' @param cores The number of cores to use for parallelization. 
#'  This can be set for an entire R session 
#'  by \code{options(loo.cores = NUMBER)}. 
#'  The default is \code{\link[parallel:detectCores]{detectCores()}}.
#' @param wcp,wtrunc Parameters used for 
#'  the Pareto smoothed importance sampling. 
#'  See \code{\link[loo:loo]{loo}} for details.
#' 
#' @details When comparing models fitted to the same data, 
#'  the smaller the LOO, the better the fit.
#' @return If just one object is provided, an object of class \code{ic}. 
#'  If multiple objects are provided, an object of class \code{iclist}.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' #model with fixed effects only
#' fit1 <- brm(rating ~ treat + period + carry,
#'             data = inhaler, family = "gaussian")
#' LOO(fit1)
#' 
#' #model with an additional random intercept for subjects
#' fit2 <- brm(rating ~ treat + period + carry + (1|subject),
#'             data = inhaler, family = "gaussian")
#' #compare both models
#' LOO(fit1, fit2)                          
#' }
#' 
#' @references 
#' Vehtari, A., Gelman, A., and Gabry, J. (2015). 
#' Efficient implementation of leave-one-out cross-validation 
#' and WAIC for evaluating fitted Bayesian models.
#' 
#' Gelman, A., Hwang, J., & Vehtari, A. (2014). 
#' Understanding predictive information criteria for Bayesian models. 
#' Statistics and Computing, 24, 997-1016.
#' 
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation 
#' and widely applicable information criterion in singular learning theory. 
#' The Journal of Machine Learning Research, 11, 3571-3594.
#' 
#' @export
LOO <- function(x, ..., compare = TRUE)
  UseMethod("LOO")

#' Interface to \pkg{shinystan}
#' 
#' Provide an interface to \pkg{shinystan} for models fitted with \pkg{brms}
#' 
#' @aliases launch_shiny.brmsfit
#' 
#' @param x A fitted model object typically of class \code{brmsfit}. 
#' @param rstudio Only relevant for RStudio users. 
#' The default (\code{rstudio=FALSE}) is to launch the app 
#' in the default web browser rather than RStudio's pop-up Viewer. 
#' Users can change the default to \code{TRUE} 
#' by setting the global option \cr \code{options(shinystan.rstudio = TRUE)}.
#' @param ... Optional arguments to pass to \code{\link[shiny:runApp]{runApp}}
#' 
#' @return An S4 shinystan object
#' 
#' @examples
#' \dontrun{
#' fit <- brm(rating ~ treat + period + carry + (1|subject),
#'            data = inhaler, family = "gaussian")
#' launch_shiny(fit)                         
#' }
#' 
#' @seealso \code{\link[shinystan:launch_shinystan]{launch_shinystan}}
#' 
#' @export
launch_shiny <- function(x, rstudio = getOption("shinystan.rstudio"), ...)
  UseMethod("launch_shiny")

#' Extract Stan Model Code
#' 
#' Extract the model code in Stan language
#' 
#' @aliases stancode.brmsfit
#' 
#' @param object An object of class \code{brmsfit}
#' @param ... Currently ignored
#' 
#' @return model code in stan language for further processing.
#' 
#' @export
stancode <- function(object, ...)
  UseMethod("stancode")

#' Extract Data passed to Stan
#' 
#' Extract all data that was used by Stan to fit the model
#' 
#' @aliases standata.brmsfit
#' 
#' @param object An object of class \code{brmsfit}
#' @param ... Currently ignored
#' 
#' @return A named list containing the data passed to Stan
#' 
#' @export
standata <- function(object, ...)
  UseMethod("standata")

#' Various Plotting Functions implemented in \pkg{rstan} 
#' 
#' Conveniant way to call plotting functions 
#' implemented in the \pkg{rstan} package. 
#' 
#' @inheritParams posterior_samples
#' @param object An R object typically of class \code{brmsfit}
#' @param pars Names of parameters to be plotted, 
#'   as given by a character vector or regular expressions. 
#'   By default, the first 10 parameters are plotted.
#' @param type The type of the plot. 
#'   Supported types are (as names) \code{plot},
#'   \code{trace}, \code{hist}, \code{dens}, \code{scat}, 
#'   \code{diag}, \code{rhat}, \code{ess}, \code{mcse}, \code{ac}. 
#'   For an overview on the various plot types see
#'   \code{\link[rstan:plotting-functions]{plotting-functions}}.
#' @param quiet A flag indicating whether messages 
#'   produced by \pkg{ggplot2} during the plotting process 
#'   should be silenced. Default is \code{FALSE}.
#' @param ... Additional arguments passed to the plotting functions.
#' 
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object 
#'   that can be further customized using the \pkg{ggplot2} package.
#' 
#' @details Instead of using \code{stanplot(<brmsfit-object>)}, 
#'   the plotting functions can be called directly 
#'   via \code{stan_<plot-type>(<brmsfit-object>$fit)}. 
#'   For more details on the plotting functions see 
#'   \code{\link[rstan:stan_plot]{Plots}} as well as 
#'   \code{\link[rstan:stan_diag]{Diagnostic plots}}.
#'   Note that the plotting functions themselves 
#'   only accept full parameter names,
#'   while \code{stanplot} allows for partial matching 
#'   and regular expressions.
#'   You should also consider using 
#'   the \pkg{shinystan} package available via method 
#'   \code{\link[brms:launch_shiny]{launch_shiny}} 
#'   in \pkg{brms} for flexible and interactive visual analysis. 
#' 
#' @examples
#' \dontrun{
#' model <- brm(count ~ log_Age_c + log_Base4_c * Trt_c 
#'              + (1|patient) + (1|visit),
#'              data = epilepsy, family = "poisson")
#' # plot 95% CIs
#' stanplot(model, type = "plot", ci_level = 0.95)
#' # equivalent to
#' stan_plot(model$fit, ci_level = 0.95)
#' 
#' # only show fixed effects in the plots
#' # this will not work when calling stan_plot directly
#' stanplot(model, pars = "^b", type = "plot", ci_level = 0.95)
#' 
#' # plot some diagnostics on the sampler
#' stanplot(model, type = "diag")
#' # equivalent to 
#' stan_diag(model$fit)                           
#' }
#' 
#' @export
stanplot <- function(object, pars, ...)
  UseMethod("stanplot")