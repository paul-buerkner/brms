brmsfit <- function(formula = NULL, family = "", link = "", data.name = "", 
                    data = data.frame(), model = "", exclude = NULL,
                    prior = brmsprior(), ranef = TRUE, autocor = NULL,
                    threshold = "", cov_ranef = NULL, fit = NA, 
                    algorithm = "sampling") {
  # brmsfit class
  version <- list(
    brms = utils::packageVersion("brms"),
    rstan = utils::packageVersion("rstan")
  )
  x <- nlist(formula, family, link, data.name, data, model, exclude, prior, 
             ranef, autocor, threshold, cov_ranef, fit, algorithm, version)
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

brmssummary <- function(formula = NULL, family = "", link = "", 
                        data.name = "", group = NULL, nobs = NULL, 
                        ngrps = NULL, chains = 1, iter = 2000, 
                        warmup = 500, thin = 1, sampler = "", 
                        autocor = NULL, fixed = NULL, random = list(), 
                        cor_pars = NULL, spec_pars = NULL, 
                        mult_pars = NULL, prior = empty_brmsprior(),
                        WAIC = "Not computed", algorithm = "sampling") {
  # brmssummary class
  x <- nlist(formula, family, link, data.name, group, nobs, ngrps, chains, 
             iter,  warmup, thin, sampler, autocor, fixed, random, cor_pars, 
             spec_pars, mult_pars, prior, WAIC, algorithm)
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
#'  Default is "b" for population-level effects. 
#'  Other typical options are "sd" or "cor". 
#'  If \code{class = NULL}, all parameters can be tested
#'  against each other, but have to be specified with their full name 
#'  (see also \code{\link[brms:parnames]{parnames}}) 
#' @param group Name of a grouping factor to evaluate only 
#'  group-level effects parameters related to this grouping factor.
#'  Ignored if \code{class} is not \code{"sd"} or \code{"cor"}.
#' @param alpha The alpha-level of the tests (default is 0.05;
#'  see 'Details' for more information).
#' @param seed A single numeric value passed to \code{set.seed} 
#'  to make results reproducible.
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
#' @inheritParams plot.brmsfit
#' @param ... Currently ignored.
#' 
#' @details Among others, \code{hypothesis} computes an 
#'  evidence ratio (\code{Evid.Ratio}) for each hypothesis. 
#'  For a directed hypothesis, this is just the posterior probability 
#'  under the hypothesis against its alternative.
#'  That is, when the hypothesis if of the form \code{a > b}, 
#'  the evidence ratio is the ratio of the posterior probability 
#'  of \code{a > b} and the posterior probability of \code{a < b}.
#'  In this example, values greater than one indicate that the evidence in
#'  favour of \code{a > b} is larger than evidence in favour of \code{a < b}.
#'  For an undirected (point) hypothesis, the evidence ratio 
#'  is a Bayes factor between the hypothesis and its alternative
#'  computed via the Savage-Dickey density ratio method.
#'  That is the posterior density at the point of interest divided
#'  by the prior density at that point.
#'  Values greater than one indicate that evidence in favour of the point
#'  hypothesis has increased after seeing the data.
#'  In order to calculate this Bayes factor, all parameters related 
#'  to the hypothesis must have proper priors
#'  and argument \code{sample_prior} of function \code{brm} 
#'  must be set to \code{TRUE}. 
#'  When interpreting Bayes factors, make sure 
#'  that your priors are reasonable and carefully chosen,
#'  as the result will depend heavily on the priors. 
#'  In particular, avoid using default priors.
#'  
#'  The argument \code{alpha} specifies the size of the credible interval
#'  (i.e., Bayesian confidence interval).
#'  For instance, if \code{alpha = 0.05} (5\%), the credible interval
#'  will contain \code{1 - alpha = 0.95} (95\%) of the posterior values.
#'  Hence, \code{alpha * 100}\% of the posterior values will lie
#'  outside of the credible interval. Although this allows testing of
#'  hypotheses in a similar manner as in the frequentist null-hypothesis
#'  testing framework, we strongly argue against using arbitrary cutoffs 
#'  (e.g., \code{p < .05}) to determine the 'existence' of an effect.
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
#'            data = kidney, family = lognormal(),
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
#'            "sd_patient_age^2 + sigma^2) = 0")
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
#'   By default, all posterior samples of all parameters are extracted.
#' @param parameters A deprecated alias of \code{pars}.  
#' @param exact_match Indicates whether parameter names 
#'   should be matched exactly or treated as regular expression. 
#'   Default is \code{FALSE}.
#' @param add_chain A flag indicating if the returned \code{data.frame} 
#'   should contain two additional columns. The \code{chain} column 
#'   indicates the chain in which each sample was generated, the \code{iter} 
#'   column indicates the iteration number within each chain.
#' @param add_chains A deprecated alias of \code{add_chain}.
#' @param subset A numeric vector indicating the rows 
#'   (i.e., posterior samples) to be returned. 
#'   If \code{NULL} (the default), all  posterior samples are returned.
#' @param as.matrix Should the output be a \code{matrix} 
#'   instead of a \code{data.frame}? Defaults to \code{FALSE}.
#' @param row.names,optional See \code{\link[base:as.data.frame]{as.data.frame}}.
#' @param ... For \code{as.data.frame} and \code{as.matrix}:
#'   Further arguments to be passed to \code{posterior_samples}.
#'   
#' @details Currently there are methods for \code{brmsfit} objects.
#'   \code{as.data.frame.brmsfit} and \code{as.matrix.brmsfit} are basically 
#'   just aliases of \code{posterior_samples.brmsfit} and differ from
#'   each other only in type of the returend object.
#' @return A data frame (or matrix) containing the posterior samples, 
#'   with one column per parameter.
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

#' @export 
posterior.samples <- function(x, pars = NA, ...) {
  # deprecated alias of posterior_samples
  warning("Method 'posterior.samples' is deprecated. ", 
          "Please use method 'posterior_samples' instead.", 
          call. = FALSE)
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
#' # extract all prior samples
#' samples1 <- prior_samples(fit)
#' head(samples1)
#' 
#' # extract prior samples for the population-level effects of 'treat'
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
#' @param x An \R object
#' @param ... Further arguments passed to or from other methods
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
#' @param ... Further arguments passed to or from other methods
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

#' Compute the WAIC
#' 
#' Compute the widely applicable information criterion (WAIC)
#' based on the posterior likelihood using the \pkg{loo} package.
#' 
#' @aliases WAIC.brmsfit waic.brmsfit waic
#' 
#' @param x A fitted model object typically of class \code{brmsfit}. 
#' @param ... Optionally more fitted model objects.
#' @param compare A flag indicating if the information criteria
#'  of the models should be compared to each other
#'  via \code{\link[brms:compare_ic]{compare_ic}}.
#' @param pointwise A flag indicating whether to compute the full
#'  log-likelihood matrix at once or separately for each observation. 
#'  The latter approach is usually considerably slower but 
#'  requires much less working memory. Accordingly, if one runs 
#'  into memory issues, \code{pointwise = TRUE} is the way to go.
#'  By default, \code{pointwise} is automatically chosen based on 
#'  the size of the model.
#' @inheritParams predict.brmsfit
#' 
#' @details When comparing models fitted to the same data, 
#'  the smaller the WAIC, the better the fit.
#'  For \code{brmsfit} objects, \code{waic} is an alias of \code{WAIC}.
#' @return If just one object is provided, an object of class \code{ic}. 
#'  If multiple objects are provided, an object of class \code{iclist}.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' # model with population-level effects only
#' fit1 <- brm(rating ~ treat + period + carry,
#'             data = inhaler, family = "gaussian")
#' WAIC(fit1)
#' 
#' # model with an additional varying intercept for subjects
#' fit2 <- brm(rating ~ treat + period + carry + (1|subject),
#'             data = inhaler, family = "gaussian")
#' # compare both models
#' WAIC(fit1, fit2)                          
#' }
#' 
#' @references 
#' Vehtari, A., Gelman, A., & Gabry J. (2016). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. In Statistics 
#' and Computing, doi:10.1007/s11222-016-9696-4. arXiv preprint arXiv:1507.04544.
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
WAIC <- function(x, ...) {
  UseMethod("WAIC")
}

#' Compute the LOO information criterion
#' 
#' Perform Leave-one-out cross-validation based on the posterior likelihood
#' using the \pkg{loo} package.
#' 
#' @aliases LOO.brmsfit loo.brmsfit loo
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
#'  For \code{brmsfit} objects, \code{loo} is an alias of \code{LOO}.
#' @return If just one object is provided, an object of class \code{ic}. 
#'  If multiple objects are provided, an object of class \code{iclist}.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' # model with population-level effects only
#' fit1 <- brm(rating ~ treat + period + carry,
#'             data = inhaler, family = "gaussian")
#' LOO(fit1)
#' 
#' # model with an additional varying intercept for subjects
#' fit2 <- brm(rating ~ treat + period + carry + (1|subject),
#'             data = inhaler, family = "gaussian")
#' # compare both models
#' LOO(fit1, fit2)                          
#' }
#' 
#' @references 
#' Vehtari, A., Gelman, A., & Gabry J. (2016). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. In Statistics 
#' and Computing, doi:10.1007/s11222-016-9696-4. arXiv preprint arXiv:1507.04544.
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
LOO <- function(x, ...) {
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

#' MCMC Plots Implemented in \pkg{bayesplot} 
#' 
#' Conveniant way to call MCMC plotting functions 
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
#'   method \code{\link[brms:launch_shiny]{launch_shiny}} 
#'   in \pkg{brms} for flexible and interactive visual analysis. 
#' 
#' @examples
#' \dontrun{
#' model <- brm(count ~ log_Age_c + log_Base4_c * Trt_c 
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

#' Display marginal effects of predictors
#' 
#' Display marginal effects of one or more numeric and/or categorical 
#' predictors including two-way interaction effects.
#' 
#' @param x An \R object usually of class \code{brmsfit}.
#' @param effects An optional character vector naming effects
#'   (main effects or interactions) for which to compute marginal plots.
#'   Interactions are specified by a \code{:} between variable names.
#'   If \code{NULL} (the default), plots are generated for all main effects
#'   and two-way interactions estimated in the model. When specifying
#'   \code{effects} manually, \emph{all} two-way interactions may be plotted
#'   even if not orginally modeled.
#' @param conditions An optional \code{data.frame} containing variable values
#'   to condition on. Each effect defined in \code{effects} will
#'   be plotted separately for each row of \code{data}. 
#'   The row names of \code{data} will be treated as titles of the subplots. 
#'   It is recommended to only define a few rows in order to keep the plots clear.
#'   If \code{NULL} (the default), numeric variables will be marginalized
#'   by using their means and factors will get their reference level assigned.
#' @param int_conditions An optional named \code{list} whose elements are numeric
#'   vectors of values of the second variables in two-way interactions. 
#'   At these values, predictions are evaluated. The names of 
#'   \code{int_conditions} have to match the variable names exactly.
#'   Additionally, the elements of the numeric vectors may be named themselves,
#'   in which case their names appear as labels for the conditions in the plots.
#'   Instead of vectors, functions returning vectors may be passed and are
#'   applied on the original values of the corresponding variable.
#'   If \code{NULL} (the default), predictions are evaluated at the 
#'   \eqn{mean} and at \eqn{mean +/- sd}. 
#' @param re_formula A formula containing random effects to be considered 
#'   in the marginal predictions. If \code{NULL}, include all random effects; 
#'   if \code{NA} (default), include no random effects.
#' @param robust If \code{TRUE} (the default) the median is used as the 
#'   measure of central tendency. If \code{FALSE} the mean is used instead.
#' @param probs The quantiles to be used in the computation of credible
#'   intervals (defaults to 2.5 and 97.5 percent quantiles)
#' @param method Either \code{"fitted"} or \code{"predict"}. 
#'   If \code{"fitted"}, plot marginal predictions of the regression curve. 
#'   If \code{"predict"}, plot marginal predictions of the responses.
#' @param surface Logical; Indicates whether interactions or 
#'   two-dimensional smooths should be visualized as a surface. 
#'   Defaults to \code{FALSE}. The surface type can be controlled 
#'   via argument \code{stype} of the related plotting method.
#' @param resolution Number of support points used to generate 
#'   the plots. Higher resolution leads to smoother plots. 
#'   Defaults to \code{100}. If \code{surface} is \code{TRUE},
#'   this implies \code{10000} support points for interaction terms,
#'   so it might be necessary to reduce \code{resolution} 
#'   when only few RAM is available.
#' @param too_far Positive number. 
#'   For surface plots only: Grid points that are too 
#'   far away from the actual data points can be excluded from the plot. 
#'   \code{too_far} determines what is too far. The grid is scaled into 
#'   the unit square and then grid points more than \code{too_far} 
#'   from the predictor variables are excluded. By default, all
#'   grid points are used. Ignored for non-surface plots.
#' @param select_points Positive number. 
#'   Only relevant if \code{points} or \code{rug} are set to \code{TRUE}: 
#'   Actual data points of numeric variables that 
#'   are too far away from the values specified in \code{conditions} 
#'   can be excluded from the plot. Values are scaled into 
#'   the unit interval and then points more than \code{select_points} 
#'   from the values in \code{conditions} are excluded. 
#'   By default, all points are used.
#' @param ncol Number of plots to display per column for each effect.
#'   If \code{NULL} (default), \code{ncol} is computed internally based
#'   on the number of rows of \code{data}.
#' @param points Logical; indicating whether the original data points
#'   should be added via \code{\link[ggplot2:geom_point]{geom_point}}.
#'   Default is \code{FALSE}. Note that only those data points will be added
#'   that match the specified conditions defined in \code{conditions}.
#'   For categorical predictors, the conditions have to match exactly. 
#'   For numeric predictors, argument \code{select_points} is used to
#'   determine, which points do match a condition.
#' @param rug Logical; indicating whether a rug representation of predictor
#'   values should be added via \code{\link[ggplot2:geom_rug]{geom_rug}}.
#'   Default is \code{FALSE}. Depends on \code{select_points} in the same
#'   way as \code{points} does.
#' @param stype Indicates how surface plots should be displayed.
#'   Either \code{"contour"} or \code{"raster"}.
#' @inheritParams plot.brmsfit
#' @param ... Further arguments such as \code{subset} or \code{nsamples}
#'   passed to \code{\link[brms:predict.brmsfit]{predict}} or 
#'   \code{\link[brms:fitted.brmsfit]{fitted}}.
#' 
#' @return An object of class \code{brmsMarginalEffects}, which is a named list
#'   with one data.frame per effect containing all information required 
#'   to generate marginal effects plots. Among others, these data.frames
#'   contain some special variables, namely \code{estimate__} (predicted values
#'   of the response), \code{se__} (standard error of the predicted response),
#'   \code{lower__} and \code{upper__} (lower and upper bounds of the uncertainty
#'   interval of the response), as well as \code{cond__} (used in faceting when 
#'   \code{conditions} contains multiple rows).
#'   
#'   The corresponding \code{plot} method returns a named 
#'   list of \code{\link[ggplot2:ggplot]{ggplot}} objects, which can be further 
#'   customized using the \pkg{ggplot2} package.
#'   
#' @details When creating \code{marginal_effects} for a particular predictor 
#'   (or interaction of two predictors), one has to choose the values of all 
#'   other predictors to condition on. 
#'   By default, the mean is used for continuous variables
#'   and the reference category is used for factors, but you may change these
#'   values via argument \code{conditions}. 
#'   This also has an implication for the \code{points} argument: 
#'   In the created plots, only those points will be shown that correspond 
#'   to the factor levels actually used in the conditioning, in order not 
#'   to create the false impressivion of bad model fit, where it is just 
#'   due to conditioning on certain factor levels.
#'   Since we condition on rather than actually marginalizing variables, 
#'   the name  \code{marginal_effects} is possibly not ideally chosen in 
#'   retrospect. 
#' 
#'   \code{NA} values within factors in \code{conditions}, 
#'   are interpreted as if all dummy variables of this factor are 
#'   zero. This allows, for instance, to make predictions of the grand mean 
#'   when using sum coding. 
#'   
#'   To fully change colours of the created plots, 
#'   one has to amend both \code{scale_colour} and \code{scale_fill}.
#'   See \code{\link[ggplot2:scale_colour_grey]{scale_colour_grey}} or
#'   \code{\link[ggplot2:scale_colour_gradient]{scale_colour_gradient}}
#'   for more details.
#' 
#' @examples 
#' \dontrun{
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1 | patient),
#'            data = epilepsy, family = poisson()) 
#'            
#' ## plot all marginal effects
#' plot(marginal_effects(fit), ask = FALSE)
#' 
#' ## change colours to grey scale
#' me <- marginal_effects(fit, "log_Base4_c:Trt_c")
#' plot(me, plot = FALSE)[[1]] + 
#'   scale_color_grey() +
#'   scale_fill_grey()
#' 
#' ## only plot the marginal interaction effect of 'log_Base4_c:Trt_c'
#' ## for different values for 'log_Age_c'
#' conditions <- data.frame(log_Age_c = c(-0.3, 0, 0.3))
#' plot(marginal_effects(fit, effects = "log_Base4_c:Trt_c", 
#'                       conditions = conditions))
#'                       
#' ## also incorporate random effects variance over patients
#' ## also add data points and a rug representation of predictor values
#' plot(marginal_effects(fit, effects = "log_Base4_c:Trt_c", 
#'                       conditions = conditions, re_formula = NULL), 
#'      points = TRUE, rug = TRUE)
#'  
#' ## change handling of two-way interactions
#' int_conditions <- list(
#'   log_Base4_c = setNames(c(-2, 1, 0), c("b", "c", "a"))
#' )
#' marginal_effects(fit, effects = "Trt_c:log_Base4_c",
#'                  int_conditions = int_conditions)
#' marginal_effects(fit, effects = "Trt_c:log_Base4_c",
#'                  int_conditions = list(log_Base4_c = quantile))        
#'      
#' ## fit a model to illustrate how to plot 3-way interactions
#' fit3way <- brm(count ~ log_Age_c * log_Base4_c * Trt_c, data = epilepsy)
#' conditions <- data.frame(log_Age_c = c(-0.3, 0, 0.3))
#' rownames(conditions) <- paste("log_Age_c =", conditions$log_Age_c)
#' marginal_effects(
#'   fit3way, "log_Base4_c:Trt_c", conditions = conditions
#' )
#' ## only include points close to the specified values of log_Age_c
#' me <- marginal_effects(
#'  fit3way, "log_Base4_c:Trt_c", conditions = conditions, 
#'  select_points = 0.1
#' )
#' plot(me, points = TRUE)
#' }
#' 
#' @export
marginal_effects <- function(x, ...) {
  UseMethod("marginal_effects")
}

#' Display Smooth Terms
#' 
#' Display smooth \code{s} and \code{t2} terms of models
#' fitted with \pkg{brms}.
#' 
#' @inheritParams marginal_effects
#' @param smooths Optional character vector of smooth terms
#'   to display. If \code{NULL} (the default) all smooth terms
#'   are shown.
#' @param subset A numeric vector specifying
#'  the posterior samples to be used. 
#'  If \code{NULL} (the default), all samples are used.
#' @param nsamples Positive integer indicating how many 
#'  posterior samples should be used. 
#'  If \code{NULL} (the default) all samples are used.
#'  Ignored if \code{subset} is not \code{NULL}.
#' @param ... Currently ignored.
#'   
#' @return For the \code{brmsfit} method, 
#' an object of class \code{brmsMarginalEffects}. See
#' \code{\link[brms:marginal_effects]{marginal_effects}} for 
#' more details and documentation of the related plotting function.
#' 
#' @details Two-dimensional smooth terms will be visualized using
#'   either contour or raster plots.
#'   
#' @examples 
#' \dontrun{
#' set.seed(0) 
#' dat <- mgcv::gamSim(1, n = 200,scale = 2)
#' fit <- brm(y ~ s(x0) + s(x1) + s(x2) + s(x3), data = dat)
#' # show all smooth terms
#' plot(marginal_smooths(fit), rug = TRUE, ask = FALSE)
#' # show only the smooth term s(x2)
#' plot(marginal_smooths(fit, smooths = "s(x2)"), ask = FALSE)
#' 
#' # fit and plot a two-dimensional smooth term
#' fit2 <- brm(y ~ t2(x0, x2), data = dat)
#' ms <- marginal_smooths(fit2)
#' plot(ms, stype = "contour")
#' plot(ms, stype = "raster")
#' }
#' 
#' @export
marginal_smooths <- function(x, ...) {
  UseMethod("marginal_smooths")
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
#' If \code{summary = FALSE}, an S x N x K arrary, where
#' S is the number of posterior samples.
#' 
#' @details 
#' The returned probabilities can be written as
#' \eqn{P(K_n = k | Y_n)}, that is the posterior probability 
#' that observation n orginiates from component k. 
#' They are computed using Bayes' Theorem
#' \deqn{P(K_n = k | Y_n) = P(Y_n | K_n = k) P(K_n = k) / P(Y_n),}
#' where \eqn{P(Y_n | K_n = k)} is the (posterior) likelihood
#' of observation n for component k, \eqn{P(K_n = k)} is 
#' the (posterior) mixing probability of component k 
#' (i.e. parameter \code{theta<k>}), and 
#' \deqn{P(Y_n) = \sum {k=1}^K P(Y_n | K_n = k) P(K_n = k)}
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
#' Export user-defined \pkg{Stan} function to the 
#' \code{\link[base:environment]{.GlobalEnv}}.
#' For more details see 
#' \code{\link[rstan:expose_stan_functions]{expose_stan_functions}}.
#' 
#' @param x An \R object
#' @param ... Further arguments
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
#' @param x an \R object
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


# ----- internal generics -----
get_re <- function(x, ...) {
  # extract group-level terms
  UseMethod("get_re")
}

get_effect <- function(x, ...) {
  # extract various kind of effects
  UseMethod("get_effect")
}

get_all_effects <- function(x, ...) {
  # extract combinations of predictor variables
  UseMethod("get_all_effects")
}

prior_effects <- function(x, ...) {
  # generate priors various kind of effects 
  UseMethod("prior_effects")
}

data_effects <- function(x, ...) {
  # generate data for various kind of effects 
  UseMethod("data_effects")
}

stan_effects <- function(x, ...) {
  # generate stan code various kind of effects 
  UseMethod("stan_effects")
}

change_effects <- function(x, ...) {
  # helps in renaming parameters after model fitting
  UseMethod("change_effects")
}

extract_draws <- function(x, ...) {
  # extract data and posterior draws
  UseMethod("extract_draws")
}

check_prior_special <- function(x, ...) {
  # prepare special priors for use in Stan
  UseMethod("check_prior_special")
}

auxpar_family <- function(family, auxpar, ...) {
  # generate a family object of an auxiliary parameter
  UseMethod("auxpar_family")
}

family_names <- function(family, ...) {
  # extract family names
  UseMethod("family_names")
}

valid_auxpars <- function(family, ...) {
  # get valid auxiliary parameters for a family
  UseMethod("valid_auxpars")
}

stan_llh <- function(family, ...) {
  UseMethod("stan_llh")
}
