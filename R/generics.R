#' Class \code{brmsfit} of models fitted with the \pkg{brms} package
#' 
#' Models fitted with the \code{\link[brms:brms]{brms}} package are 
#' represented as a \code{brmsfit} object, which contains the posterior 
#' samples, model formula, Stan code, relevant data, and other information.
#' 
#' @name brmsfit-class
#' @aliases brmsfit
#' @docType class
#' 
#' @details 
#' See \code{methods(class = "brmsfit")} for an overview of available methods.
#' 
#' @slot formula A \code{\link{brmsformula}} object
#' @slot family A \code{\link{brmsfamily}} object
#' @slot data A \code{data.frame} containing all variables used in the model
#' @slot data.name The name of \code{data} as specified by the user 
#' @slot model The model code in \pkg{Stan} language
#' @slot prior A \code{\link{brmsprior}} object containing
#'   information on the priors used in the model
#' @slot autocor An \code{\link{cor_brms}} object containing 
#'   the autocorrelation structure if specified
#' @slot ranef A \code{data.frame} containing the group-level structure
#' @slot cov_ranef A \code{list} of customized group-level covariance matrices
#' @slot stanvars A \code{\link{stanvars}} object or \code{NULL}
#' @slot stan_funs A character string of length one or \code{NULL}
#' @slot loo An empty slot for adding the \code{\link{loo}} 
#'   criterion after model fitting
#' @slot waic An empty slot for adding the \code{\link{waic}} 
#'   criterion after model fitting
#' @slot kfold An empty slot for adding the \code{\link{kfold}} 
#'   criterion after model fitting  
#' @slot R2 An empty slot for adding the \code{\link{bayes_R2}}
#'   (Bayesian R-squared) value after model fitting 
#' @slot marglik An empty slot for adding a \code{bridge} object 
#'   after model fitting containing the log marginal likelihood 
#'   (see \code{\link{bridge_sampler}} for details)
#' @slot fit An object of class \code{\link[rstan:stanfit]{stanfit}}
#'   among others containing the posterior samples
#' @slot exclude The names of the parameters for which samples are not saved
#' @slot algorithm The name of the algorithm used to fit the model
#' @slot version The versions of \pkg{brms} and \pkg{rstan} with 
#'   which the model was fitted
#' 
#' @seealso 
#'   \code{\link{brms}}, 
#'   \code{\link{brm}}, 
#'   \code{\link{brmsformula}}, 
#'   \code{\link{brmsfamily}}
#' 
NULL

brmsfit <- function(formula = NULL, family = NULL, data = data.frame(), 
                    data.name = "", model = "", prior = empty_brmsprior(), 
                    autocor = NULL, ranef = empty_ranef(), cov_ranef = NULL, 
                    loo = NULL, waic = NULL, kfold = NULL, R2 = NULL,
                    marglik = NULL, stanvars = NULL, stan_funs = NULL, 
                    fit = NA, exclude = NULL, algorithm = "sampling") {
  # brmsfit class
  version <- list(
    brms = utils::packageVersion("brms"),
    rstan = utils::packageVersion("rstan")
  )
  x <- nlist(
    formula, family, data, data.name, model, prior,
    autocor, ranef, cov_ranef, loo, waic, kfold, R2, marglik,
    stanvars, stan_funs, fit, exclude, algorithm, version
  )
  class(x) <- "brmsfit"
  x
}

#' Checks if argument is a \code{brmsfit} object
#' 
#' @param x An \R object
#' 
#' @export
is.brmsfit <- function(x) {
  inherits(x, "brmsfit")
}

#' Checks if argument is a \code{brmsfit_multiple} object
#' 
#' @param x An \R object
#' 
#' @export
is.brmsfit_multiple <- function(x) {
  inherits(x, "brmsfit_multiple")
}

#' Descriptions of \code{brmshypothesis} Objects
#' 
#' A \code{brmshypothesis} object contains posterior samples
#' as well as summary statistics of non-linear hypotheses as 
#' returned by \code{\link[brms:hypothesis]{hypothesis}}.
#' 
#' @name brmshypothesis
#' 
#' @param ignore_prior A flag indicating if prior distributions 
#'  should also be plotted. Only used if priors were specified on
#'  the relevant parameters.
#' @param digits Minimal number of significant digits, 
#'   see \code{\link[base:print.default]{print.default}}.
#' @param chars Maximum number of characters of each hypothesis
#'  to print or plot. If \code{NULL}, print the full hypotheses.
#'  Defaults to \code{20}.
#' @param colors Two values specifying the colors of the posterior
#'  and prior density respectively. If \code{NULL} (the default)
#'  colors are taken from the current color scheme of 
#'  the \pkg{bayesplot} package.
#' @param ... Currently ignored.
#' @inheritParams plot.brmsfit
#' 
#' @details 
#' The two most important elements of a \code{brmshypothesis} object are
#' \code{hypothesis}, which is a data.frame containing the summary estimates
#' of the hypotheses, and \code{samples}, which is a data.frame containing 
#' the corresponding posterior samples.
#' 
#' @seealso \code{\link[brms:hypothesis]{hypothesis}}
NULL

#' Extract posterior samples
#' 
#' Extract posterior samples of specified parameters 
#' 
#' @aliases posterior.samples posterior_samples.brmsfit posterior.samples.brmsfit
#' 
#' @param x An \code{R} object typically of class \code{brmsfit}
#' @param pars Names of parameters for which posterior samples 
#'   should be returned, as given by a character vector or regular expressions.
#'   By default, all posterior samples of all parameters are extracted.
#' @param exact_match Indicates whether parameter names 
#'   should be matched exactly or treated as regular expression. 
#'   Default is \code{FALSE}.
#' @param add_chain A flag indicating if the returned \code{data.frame} 
#'   should contain two additional columns. The \code{chain} column 
#'   indicates the chain in which each sample was generated, the \code{iter} 
#'   column indicates the iteration number within each chain.
#' @param subset A numeric vector indicating the rows 
#'   (i.e., posterior samples) to be returned. 
#'   If \code{NULL} (the default), all  posterior samples are returned.
#' @param as.matrix Should the output be a \code{matrix} 
#'   instead of a \code{data.frame}? Defaults to \code{FALSE}.
#' @param as.array Should the output be an \code{array} 
#'   instead of a \code{data.frame}? Defaults to \code{FALSE}.
#' @param row.names,optional See \code{\link[base:as.data.frame]{as.data.frame}}.
#' @param ... For \code{as.data.frame}, \code{as.matrix}, and \code{as.array}:
#'   Further arguments to be passed to \code{posterior_samples}.
#'   
#' @details Currently there are methods for \code{brmsfit} objects.
#'   \code{as.data.frame.brmsfit}, \code{as.matrix.brmsfit}, and
#'   \code{as.array.brmsfit} are basically aliases of 
#'   \code{posterior_samples.brmsfit} and differ from
#'   each other only in type of the returned object.
#'   
#' @return A data frame (matrix or array) containing the posterior samples, 
#'   with one column per parameter. In case an array is returned,
#'   it contains one additional dimension for the chains.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' fit <- brm(rating ~ treat + period + carry + (1|subject), 
#'            data = inhaler, family = "cumulative")
#' 
#' # extract posterior samples of population-level effects 
#' samples1 <- posterior_samples(fit, "^b")
#' head(samples1)
#' 
#' # extract posterior samples of group-level standard deviations
#' samples2 <- posterior_samples(fit, "^sd_")
#' head(samples2)
#' }
#' 
#' @export 
posterior_samples <- function(x, pars = NA, ...) {
  UseMethod("posterior_samples")
}

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
#' @param ... Currently ignored
#'   
#' @details To make use of this function, 
#'  the model must contain samples of prior distributions.
#'  This can be ensured by setting \code{sample_prior = TRUE} 
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
#'            sample_prior = TRUE)
#' 
#' # extract all prior samples
#' samples1 <- prior_samples(fit)
#' head(samples1)
#' 
#' # extract prior samples for the population-level effects of 'treat'
#' samples2 <- prior_samples(fit, "b_treat")
#' head(samples2)
#' }
#' 
#' @export 
prior_samples <- function(x, pars = NA, ...) {
  UseMethod("prior_samples")
}

#' Extract Parameter Names
#' 
#' Extract all parameter names of a given model.
#'  
#' @aliases par.names parnames.brmsfit par.names.brmsfit
#' 
#' @param x An \R object
#' @param ... Further arguments passed to or from other methods.
#' 
#' @details Currently there are methods for \code{brmsfit} objects.
#' 
#' @return A character vector containing the parameter names of the model.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @export
parnames <- function(x, ...) {
  UseMethod("parnames")
}

#' Number of Posterior Samples
#' 
#' Extract the number of posterior samples 
#' stored in a fitted Bayesian model.
#' 
#' @param x An \R object
#' @param ... Further arguments passed to or from other methods.
#' @param subset An optional integer vector defining a 
#'   subset of samples to be considered.
#' @param incl_warmup A flag indicating whether to also 
#'   count warmup / burn-in samples.
#' 
#' @details Currently there are methods for \code{brmsfit} objects.
#' 
#' @export
nsamples <- function(x, ...) {
  UseMethod("nsamples")
}

#' Number of levels
#' 
#' Extract the number of levels of one or more grouping factors.
#' 
#' @aliases ngrps.brmsfit
#' 
#' @param object An \R object.
#' @param ... Currently ignored.
#' 
#' @return A named list containing the number of levels per
#'   grouping factor.
#' 
#' @export
ngrps <- function(object, ...) {
  UseMethod("ngrps")
}

#' Extract Stan model code
#' 
#' @aliases stancode.brmsfit
#' 
#' @param object An object of class \code{brmsfit}
#' @param version Logical; indicates if the first line containing
#'   the \pkg{brms} version number should be included.
#'   Defaults to \code{TRUE}.
#' @param ... Currently ignored
#' 
#' @return Stan model code for further processing.
#' 
#' @export
stancode <- function(object, ...) {
  UseMethod("stancode")
}

#' Extract Data passed to Stan
#' 
#' Extract all data that was used by Stan to fit the model
#' 
#' @aliases standata.brmsfit
#' 
#' @param object An object of class \code{brmsfit}.
#' @param internal Logical, indicates if the data should be prepared 
#'   for internal use in other post-processing methods.
#' @param control A named list currently for internal usage only.
#' @param ... More arguments passed to \code{\link{make_standata}}.
#' @inheritParams predict.brmsfit
#' 
#' @return A named list containing the data originally passed to Stan.
#' 
#' @export
standata <- function(object, ...) {
  UseMethod("standata")
}

#' Extract Autocorrelation Structures
#' 
#' @param object An \R object.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return What exactly is returned depends on the specific method.
#' 
#' @export
autocor <- function(object, ...) {
  UseMethod("autocor")
}

#' MCMC Plots Implemented in \pkg{bayesplot} 
#' 
#' Convenient way to call MCMC plotting functions 
#' implemented in the \pkg{bayesplot} package. 
#' 
#' @inheritParams posterior_samples
#' @param object An \R object typically of class \code{brmsfit}
#' @param pars Names of parameters to be plotted, 
#'   as given by a character vector or regular expressions. 
#'   By default, all parameters except for group-level and 
#'   smooth effects are plotted. May be ignored for some plots.
#' @param type The type of the plot. 
#'   Supported types are (as names) \code{hist}, \code{dens}, 
#'   \code{hist_by_chain}, \code{dens_overlay}, 
#'   \code{violin}, \code{intervals}, \code{areas}, \code{acf}, 
#'   \code{acf_bar},\code{trace}, \code{trace_highlight}, \code{scatter},
#'   \code{rhat}, \code{rhat_hist}, \code{neff}, \code{neff_hist}
#'   \code{nuts_acceptance}, \code{nuts_divergence},
#'   \code{nuts_stepsize}, \code{nuts_treedepth}, and \code{nuts_energy}. 
#'   For an overview on the various plot types see
#'   \code{\link[bayesplot:MCMC-overview]{MCMC-overview}}.
#' @param ... Additional arguments passed to the plotting functions.
#'   See \code{\link[bayesplot:MCMC-overview]{MCMC-overview}} for
#'   more details.
#' 
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object 
#'   that can be further customized using the \pkg{ggplot2} package.
#' 
#' @details 
#'   Also consider using the \pkg{shinystan} package available via 
#'   method \code{\link{launch_shinystan}} in \pkg{brms} for flexible 
#'   and interactive visual analysis. 
#' 
#' @examples
#' \dontrun{
#' model <- brm(count ~ log_Age_c + log_Base4_c * Trt 
#'              + (1|patient) + (1|visit),
#'              data = epilepsy, family = "poisson")
#'              
#' # plot posterior intervals
#' stanplot(model)
#' 
#' # only show population-level effects in the plots
#' stanplot(model, pars = "^b_")
#' 
#' # show histograms of the posterior distributions
#' stanplot(model, type = "hist")
#' 
#' # plot some diagnostics of the sampler
#' stanplot(model, type = "neff")
#' stanplot(model, type = "rhat")
#' 
#' # plot some diagnostics specific to the NUTS sampler
#' stanplot(model, type = "nuts_acceptance")
#' stanplot(model, type = "nuts_divergence")
#' }
#' 
#' @export
stanplot <- function(object, ...) {
  UseMethod("stanplot")
}

#' Model Weighting Methods
#' 
#' Compute model weights in various ways, for instance via
#' stacking of predictive distributions, Akaike weights, or
#' marginal likelihoods.
#' 
#' @inheritParams loo.brmsfit
#' @param weights Name of the criterion to compute weights from. 
#'   Should be one of \code{"loo"}, \code{"waic"}, \code{"kfold"}, 
#'   \code{"loo2"} (current default), or \code{"marglik"}. 
#'   For the former three options, Akaike weights will be computed
#'   based on the information criterion values returned by
#'   the respective methods. For \code{"loo2"}, method
#'   \code{\link{loo_model_weights}} will be used to obtain weights. 
#'   For \code{"marglik"}, method \code{\link{post_prob}} 
#'   will be used to compute weights based on log marginal 
#'   likelihood values (make sure to specify reasonable priors in 
#'   this case). Alternatively, \code{weights} can be a numeric vector 
#'   of pre-specified weights.
#'   
#' @return A numeric vector of weights for the models.
#'   
#' @examples 
#' \dontrun{
#' # model with 'treat' as predictor
#' fit1 <- brm(rating ~ treat + period + carry, data = inhaler)
#' summary(fit1)
#' 
#' # model without 'treat' as predictor
#' fit2 <- brm(rating ~ period + carry, data = inhaler)
#' summary(fit2)
#' 
#' # obtain Akaike weights based on the WAIC
#' model_weights(fit1, fit2, weights = "waic")
#' }
#' 
#' @export
model_weights <- function(x, ...) {
  UseMethod("model_weights")
}

#' Posterior predictive samples averaged across models
#' 
#' Compute posterior predictive samples averaged across models.
#' Weighting can be done in various ways, for instance using
#' Akaike weights based on information criteria or 
#' marginal likelihoods.
#' 
#' @inheritParams model_weights
#' @param method Type of predictions to average. Should be one of 
#'   \code{"predict"} (default), \code{"fitted"}, or \code{"residuals"}. 
#' @param control Optional \code{list} of further arguments 
#'   passed to the function specified in \code{weights}.
#' @param nsamples Total number of posterior samples to use.
#' @param seed A single numeric value passed to \code{\link{set.seed}}
#'   to make results reproducible.
#' @param summary Should summary statistics 
#'   (i.e. means, sds, and 95\% intervals) be returned
#'  instead of the raw values? Default is \code{TRUE}.
#' @param robust If \code{FALSE} (the default) the mean is used as 
#'  the measure of central tendency and the standard deviation as 
#'  the measure of variability. If \code{TRUE}, the median and the 
#'  median absolute deviation (MAD) are applied instead.
#'  Only used if \code{summary} is \code{TRUE}.
#' @param probs  The percentiles to be computed by the \code{quantile} 
#'  function. Only used if \code{summary} is \code{TRUE}. 
#' 
#' @return Same as the output of the method specified 
#'   in argument \code{method}.
#'   
#' @details Weights are computed with the \code{\link{model_weights}} method.
#'   
#' @seealso \code{\link{model_weights}}, \code{\link{posterior_average}}
#'   
#' @examples 
#' \dontrun{
#' # model with 'treat' as predictor
#' fit1 <- brm(rating ~ treat + period + carry, data = inhaler)
#' summary(fit1)
#' 
#' # model without 'treat' as predictor
#' fit2 <- brm(rating ~ period + carry, data = inhaler)
#' summary(fit2)
#' 
#' # compute model-averaged predicted values
#' (df <- unique(inhaler[, c("treat", "period", "carry")]))
#' pp_average(fit1, fit2, newdata = df)
#' 
#' # compute model-averaged fitted values
#' pp_average(fit1, fit2, method = "fitted", newdata = df)
#' }
#' 
#' @export
pp_average <- function(x, ...) {
  UseMethod("pp_average")
}

#' Posterior samples of parameters averaged across models
#' 
#' Extract posterior samples of parameters averaged across models.
#' Weighting can be done in various ways, for instance using
#' Akaike weights based on information criteria or 
#' marginal likelihoods.
#' 
#' @inheritParams pp_average
#' @param pars Names of parameters for which to average across models.
#'   Only those parameters can be averaged that appear in every model.
#'   Defaults to all overlapping parameters.
#' @param missing An optional numeric value or a named list of numeric values 
#'   to use if a model does not contain a parameter for which posterior samples 
#'   should be averaged. Defaults to \code{NULL}, in which case only those
#'   parameters can be averaged that are present in all of the models.
#' 
#' @return A \code{data.frame} of posterior samples. Samples are rows
#'   and parameters are columns.
#' 
#' @details Weights are computed with the \code{\link{model_weights}} method.
#' 
#' @seealso \code{\link{model_weights}}, \code{\link{pp_average}}
#'   
#' @examples 
#' \dontrun{
#' # model with 'treat' as predictor
#' fit1 <- brm(rating ~ treat + period + carry, data = inhaler)
#' summary(fit1)
#' 
#' # model without 'treat' as predictor
#' fit2 <- brm(rating ~ period + carry, data = inhaler)
#' summary(fit2)
#' 
#' # compute model-averaged posteriors of overlapping parameters
#' posterior_average(fit1, fit2, weights = "waic")
#' }
#' 
#' @export
posterior_average <- function(x, ...) {
  UseMethod("posterior_average")
}

#' Posterior Probabilities of Mixture Component Memberships
#' 
#' Compute the posterior probabilities of mixture component 
#' memberships for each observation including uncertainty
#' estimates.
#' 
#' @inheritParams predict.brmsfit
#' @param x An \R object usually of class \code{brmsfit}.
#' @param log Logical; Indicates whether to return 
#'   probabilities on the log-scale.
#' 
#' @return 
#' If \code{summary = TRUE}, an N x E x K array,
#' where N is the number of observations, K is the number
#' of mixture components, and E is equal to \code{length(probs) + 2}.
#' If \code{summary = FALSE}, an S x N x K array, where
#' S is the number of posterior samples.
#' 
#' @details 
#' The returned probabilities can be written as
#' \eqn{P(Kn = k | Yn)}, that is the posterior probability 
#' that observation n originates from component k. 
#' They are computed using Bayes' Theorem
#' \deqn{P(Kn = k | Yn) = P(Yn | Kn = k) P(Kn = k) / P(Yn),}
#' where \eqn{P(Yn | Kn = k)} is the (posterior) likelihood
#' of observation n for component k, \eqn{P(Kn = k)} is 
#' the (posterior) mixing probability of component k 
#' (i.e. parameter \code{theta<k>}), and 
#' \deqn{P(Yn) = \sum (k=1,...,K) P(Yn | Kn = k) P(Kn = k)}
#' is a normalizing constant.
#' 
#' @examples 
#' \dontrun{
#' ## simulate some data
#' set.seed(1234)
#' dat <- data.frame(
#'   y = c(rnorm(100), rnorm(50, 2)), 
#'   x = rnorm(150)
#' )
#' ## fit a simple normal mixture model
#' mix <- mixture(gaussian, nmix = 2)
#' prior <- c(
#'   prior(normal(0, 5), Intercept, nlpar = mu1),
#'   prior(normal(0, 5), Intercept, nlpar = mu2),
#'   prior(dirichlet(2, 2), theta)
#' )
#' fit1 <- brm(bf(y ~ x), dat, family = mix,
#'             prior = prior, chains = 2, inits = 0)
#' summary(fit1)
#'    
#' ## compute the membership probabilities         
#' ppm <- pp_mixture(fit1)
#' str(ppm)
#' 
#' ## extract point estimates for each observation
#' head(ppm[, 1, ])
#' 
#' ## classify every observation according to 
#' ## the most likely component
#' apply(ppm[, 1, ], 1, which.max)
#' }
#' 
#' @export
pp_mixture <- function(x, ...) {
  UseMethod("pp_mixture")
}

#' Expose user-defined \pkg{Stan} functions
#' 
#' Export user-defined \pkg{Stan} function and
#' optionally vectorize them. For more details see 
#' \code{\link[rstan:expose_stan_functions]{expose_stan_functions}}.
#' 
#' @param x An \R object
#' @param vectorize Logical; Indicates if the exposed functions
#'   should be vectorized via \code{\link{Vectorize}}. 
#'   Defaults to \code{FALSE}.
#' @param env Environment where the functions should be made
#'   available. Defaults to the global environment.
#' @param ... Further arguments passed to 
#'   \code{\link[rstan:expose_stan_functions]{expose_stan_functions}}.
#' 
#' @export
expose_functions <- function(x, ...) {
  UseMethod("expose_functions")
}

#' Extract Control Parameters of the NUTS Sampler
#' 
#' Extract control parameters of the NUTS sampler such as 
#' \code{adapt_delta} or \code{max_treedepth}.
#' 
#' @param x An \R object
#' @param pars Optional names of the control parameters to be returned.
#'  If \code{NULL} (the default) all control parameters are returned.
#'  See \code{\link[rstan:stan]{stan}} for more details.
#' @param ... Currently ignored.
#' 
#' @return A named \code{list} with control parameter values.
#' 
#' @export
control_params <- function(x, ...) {
  UseMethod("control_params")
}

#' Summarize Posterior Samples
#' 
#' Summarizes posterior samples based on point estimates (mean or median),
#' estimation errors (SD or MAD) and quantiles.
#' 
#' @param x An \R object.
#' @param probs The percentiles to be computed by the 
#'   \code{quantile} function.
#' @param robust If \code{FALSE} (the default) the mean is used as 
#'  the measure of central tendency and the standard deviation as 
#'  the measure of variability. If \code{TRUE}, the median and the 
#'  median absolute deviation (MAD) are applied instead.
#' @param ... More arguments passed to or from other methods.
#' @inheritParams posterior_samples
#' 
#' @return A matrix where rows indicate parameters 
#'  and columns indicate the summary estimates.
#'  
#' @examples 
#' \dontrun{
#' fit <- brm(time ~ age * sex, data = kidney)
#' posterior_summary(fit)
#' }
#' 
#' @export
posterior_summary <- function(x, ...) {
  UseMethod("posterior_summary")
}

#' Extract Diagnostic Quantities of \pkg{brms} Models
#' 
#' Extract quantities that can be used to diagnose sampling behavior
#' of the algorithms applied by \pkg{Stan} at the back-end of \pkg{brms}.
#' 
#' @name diagnostic-quantities
#' @aliases log_posterior nuts_params rhat neff_ratio
#'     
#' @param object A \code{brmsfit} object.
#' @param pars An optional character vector of parameter names. 
#'   For \code{nuts_params} these will be NUTS sampler parameter 
#'   names rather than model parameters. If pars is omitted 
#'   all parameters are included.
#' @param ... Arguments passed to individual methods.
#' 
#' @return The exact form of the output depends on the method.
#' 
#' @details For more details see 
#'   \code{\link[bayesplot:bayesplot-extractors]{bayesplot-extractors}}.
#'   
#' @examples 
#' \dontrun{
#' fit <- brm(time ~ age * sex, data = kidney)
#' 
#' lp <- log_posterior(fit)
#' head(lp)
#' 
#' np <- nuts_params(fit)
#' str(np)
#' # extract the number of divergence transitions
#' sum(subset(np, Parameter == "divergent__")$Value)
#' 
#' head(rhat(fit))
#' head(neff_ratio(fit))
#' }
NULL
