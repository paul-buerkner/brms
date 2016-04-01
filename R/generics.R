brmsfit <- function(formula = NULL, family = "", link = "", data.name = "", 
                    data = data.frame(), model = "", exclude = NULL,
                    prior = prior_frame(), ranef = TRUE, autocor = NULL,
                    nonlinear = NULL, partial = NULL, threshold = "", 
                    cov_ranef = NULL, fit = NA, algorithm = "sampling") {
  # brmsfit class
  x <- nlist(formula, family, link, data.name, data, model, exclude, prior, 
             ranef, autocor, nonlinear, partial, threshold, cov_ranef, fit, 
             algorithm)
  class(x) <- "brmsfit"
  x
}

brmssummary <- function(formula = NULL, family = "", link = "", 
                        data.name = "", group = NULL, nobs = NULL, 
                        ngrps = NULL, chains = 1, iter = 2000, 
                        warmup = 500, thin = 1, sampler = "", 
                        nonlinear = NULL, autocor = NULL, fixed = NULL, 
                        random = list(), cor_pars = NULL, spec_pars = NULL, 
                        mult_pars = NULL, WAIC = "Not computed",
                        algorithm = "sampling") {
  # brmssummary class
  x <- nlist(formula, family, link, data.name, group, nobs, ngrps, chains, 
             iter,  warmup, thin, sampler, nonlinear, autocor, fixed, 
             random, cor_pars, spec_pars, mult_pars, WAIC, algorithm)
  class(x) <- "brmssummary"
  x
}

#' Non-linear hypothesis testing
#' 
#' Perform non-linear hypothesis testing for all model parameters. 
#' 
#' @aliases hypothesis.brmsfit
#' 
#' @param x An \code{R} object typically of class \code{brmsfit}.
#' @param hypothesis A character vector specifying one or more 
#'  non-linear hypothesis concerning parameters of the model.
#' @param class A string specifying the class of parameters being tested. 
#'  Default is "b" for fixed effects. 
#'  Other typical options are "sd" or "cor". 
#'  If \code{class = NULL}, all parameters can be tested
#'  against each other, but have to be specified with their full name 
#'  (see also \code{\link[brms:parnames]{parnames}}) 
#' @param group Name of a grouping factor to evaluate only 
#'  random effects parameters related to this grouping factor.
#'  Ignored if \code{class} is not \code{"sd"} or \code{"cor"}.
#' @param alpha The alpha-level of the tests (default is 0.05).
#' @param ignore_prior A flag indicating if prior distributions 
#'  should also be plotted. Only used if priors were specified on
#'  the relevant parameters.
#' @param digits Minimal number of significant digits, 
#'   see \code{\link[base:print.default]{print.default}}.
#' @param chars Maximum number of characters of each hypothesis
#'  to print or plot. If \code{NULL}, print the full hypotheses.
#'  Defaults to \code{20}.
#' @inheritParams plot.brmsfit
#' @param ... Currently ignored.
#' 
#' @details Among others, \code{hypothesis} computes an 
#'  evidence ratio for each hypothesis. 
#'  For a directed hypothesis, this is just the posterior probability 
#'  under the hypothesis against its alternative.
#'  For an undirected (i.e. point) hypothesis the evidence ratio 
#'  is a Bayes factor between the hypothesis and its alternative.
#'  In order to calculate this Bayes factor, all parameters related 
#'  to the hypothesis must have proper priors
#'  and argument \code{sample_prior} of function \code{brm} 
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
#' ## define priors
#' prior <- c(set_prior("normal(0,2)", class = "b"),
#'            set_prior("student_t(10,0,1)", class = "sigma"),
#'            set_prior("student_t(10,0,1)", class = "sd"))
#' 
#' ## fit a linear mixed effects models
#' fit <- brm(time ~ age + sex + disease + (1 + age|patient),
#'            data = kidney, family = gaussian("log"),
#'            prior = prior, sample_prior = TRUE, 
#'            control = list(adapt_delta = 0.95))
#' 
#' ## perform two-sided hypothesis testing
#' (hyp1 <- hypothesis(fit, "sexfemale = age + diseasePKD"))
#' plot(hyp1)
#' hypothesis(fit, "exp(age) - 3 = 0", alpha = 0.01)
#' 
#' ## perform one-sided hypothesis testing
#' hypothesis(fit, "diseasePKD + diseaseGN - 3 < 0")
#' 
#' hypothesis(fit, "age < Intercept", 
#'            class = "sd", group  = "patient")
#' 
#' ## test the amount of random intercept variance on all variance
#' h <- paste("sd_patient_Intercept^2 / (sd_patient_Intercept^2 +",
#'            "sd_patient_age^2 + sigma_time^2) = 0")
#' (hyp2 <- hypothesis(fit, h, class = NULL))
#' plot(hyp2)
#' 
#' ## test more than one hypothesis at once
#' (hyp3 <- hypothesis(fit, c("diseaseGN = diseaseAN", 
#'                            "2 * diseaseGN - diseasePKD = 0")))
#' plot(hyp3, ignore_prior = TRUE)
#' }
#' 
#' @export
hypothesis <- function(x, hypothesis, ...) {
  UseMethod("hypothesis")
}

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
#' @param add_chain A flag indicating if the returned \code{data.frame} 
#'   should contain two additional columns. The \code{chain} column 
#'   indicates the chain in which each sample was generated, the \code{iter} 
#'   column indicates the iteration number within each chain.
#' @param add_chains A deprecated alias of \code{add_chain}.
#'   Note that the \code{chain} column will be named \code{chains} instead.
#' @param subset A numeric vector indicating the rows 
#'   (i.e., posterior samples) to be returned. 
#'   If \code{NULL} (the default), all  posterior samples are returned.
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
posterior_samples <- function(x, pars = NA, ...) {
  UseMethod("posterior_samples")
}

# deprecated alias of posterior_samples
#' @export 
posterior.samples <- function(x, pars = NA, ...) {
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
#' @param parameters A deprecated alias of \code{pars}       
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
prior_samples <- function(x, pars = NA, ...) {
  UseMethod("prior_samples")
}

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
parnames <- function(x, ...) {
  UseMethod("parnames")
}

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
#' @inheritParams predict.brmsfit
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
WAIC <- function(x, ..., compare = TRUE) {
  UseMethod("WAIC")
}

#' Compute LOO
#' 
#' Compute Leave-one-out cross-validation based on the posterior likelihood
#' by using the \pkg{loo} package
#' 
#' @aliases LOO.brmsfit
#' 
#' @inheritParams WAIC
#' @param cores The number of cores to use for parallelization. 
#'  Default is \code{1}.
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
LOO <- function(x, ..., compare = TRUE) {
  UseMethod("LOO")
}

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
launch_shiny <- function(x, rstudio = getOption("shinystan.rstudio"), ...) {
  UseMethod("launch_shiny")
}

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
stancode <- function(object, ...) {
  UseMethod("stancode")
}

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
standata <- function(object, ...) {
  UseMethod("standata")
}

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
stanplot <- function(object, pars, ...) {
  UseMethod("stanplot")
}

#' Display marginal effects of predictors
#' 
#' Display marginal effects of one or more numeric and/or categorical 
#' predictors including interaction effects of order 2.
#' 
#' @param x An object usually of class \code{brmsfit}
#' @param effects An optional character vector naming effects
#'   (main effects or interactions) for which to compute marginal plots.
#'   If \code{NULL} (the default), plots for all effects are generated.
#' @param conditions An optional \code{data.frame} containing variable values
#'   to marginalize on. Each effect defined in \code{effects} will
#'   be plotted separately for each row of \code{data}. 
#'   The row names of \code{data} will be treated as titles of the subplots. 
#'   It is recommended to only define a few rows in order to keep the plots clear.
#'   If \code{NULL} (the default), numeric variables will be marginalized
#'   by using their means and factors will get their reference level assigned.   
#' @param re_formula A formula containing random effects to be considered 
#'   in the marginal predictions. If \code{NULL}, include all random effects; 
#'   if \code{NA} (default), include no random effects.
#' @param probs The quantiles to be used in the computation of credible
#'   intervals (defaults to 2.5 and 97.5 percent quantiles)
#' @param method Either \code{"fitted"} or \code{"predict"}. 
#'   If \code{"fitted"}, plot marginal predictions of the regression curve. 
#'   If \code{"predict"}, plot marginal predictions of the responses.
#' @param ncol Number of plots to display per column for each effect.
#'   If \code{NULL} (default), \code{ncol} is computed internally based
#'   on the number of rows of \code{data}.
#' @param points Logical; indicating whether the original data points
#'   should be added via \code{\link[ggplot2:geom_point]{geom_point}}.
#'   Default is \code{FALSE}. Note that only those data points will be added
#'   that match the specified conditions defined in \code{conditions}.
#' @param rug Logical; indicating whether a rug representation of predictor
#'   values should be added via \code{\link[ggplot2:geom_rug]{geom_rug}}.
#'   Default is \code{FALSE}.
#' @inheritParams plot.brmsfit
#' @param ... Currently ignored.
#' 
#' @return An object of class \code{brmsMarginalEffects}, which is a named list
#'   with one element per effect containing all information required to generate
#'   marginal effects plots. The corresponding \code{plot} method returns a named 
#'   list of \code{\link[ggplot2:ggplot]{ggplot}} objects, which can be further 
#'   customized using the \pkg{ggplot2} package.
#' 
#' @examples 
#' \dontrun{
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1 | patient),
#'            data = epilepsy, family = poisson()) 
#' ## plot all marginal effects
#' plot(marginal_effects(fit), ask = FALSE)
#' ## only plot the marginal interaction effect of 'log_Base4_c:Trt_c'
#' ## for different values for 'log_Age_c'
#' mdata <- data.frame(log_Age_c = c(-0.3, 0, 0.3))
#' plot(marginal_effects(fit, effects = "log_Base4_c:Trt_c", 
#'                       data = mdata))
#' ## also incorporate random effects variance over patients
#' ## and add a rug representation of predictor values
#' plot(marginal_effects(fit, effects = "log_Base4_c:Trt_c", 
#'                       data = mdata, re_formula = NULL), rug = TRUE)
#' }
#' 
#' @export
marginal_effects <- function(x, ...) {
  UseMethod("marginal_effects")
}

#' Expose user-defined \pkg{Stan} functions
#' 
#' Export user-defined \pkg{Stan} function to the 
#' \code{\link[base:environment]{.GlobalEnv}}.
#' For more details see 
#' \code{\link[rstan:expose_stan_functions]{expose_stan_functions}}.
#' 
#' @param x An \code{R} object
#' @param ... Further arguments
#' 
#' @export
expose_functions <- function(x, ...) {
  UseMethod("expose_functions")
}
