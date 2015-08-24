#' Fit Bayesian Generalized Linear and Ordinal Mixed Models
#' 
#' Fit a Bayesian generalized linear or ordinal mixed model using Stan
#' 
#' @param formula An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. 
#'   The details of model specification are given under 'Details'.
#' @param data An optional data frame, list or environment  (or object coercible by \code{as.data.frame} to a data frame) containing 
#'  the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, 
#'  typically the environment from which \code{brm} is called.
#' @param family A vector of one or two character strings. The first string indicates the distribution of the dependent variable (the 'family'). Currently, the following families are supported:
#'  \code{"gaussian"}, \code{"student"}, \code{"cauchy"}, \code{"binomial"}, \code{"bernoulli"}, \code{"categorical"}, \code{"poisson"}, \code{"negbinomial"}, \cr
#'  \code{"geometric"}, \code{"gamma"}, \code{"exponential"}, \code{"weibull"}, \code{"cumulative"}, \code{"cratio"}, \code{"sratio"}, and \code{"acat"}.
#'  The second string indicates the link function, which must be supported by the distribution of the dependent variable. 
#'  If not specified, default link functions are used. Further information is provided under 'Details'.
#' @param prior A named list of character strings specifing the prior distributions of the parameters. Further information
#'  is provided under 'Details'.
#' @param addition A named list of one sided formulas each containing additional information on the response variable. The following names are allowed:
#'  \code{se} for specifying standard errors for meta-analysis, \code{weights} to fit weighted regression models, 
#'  \code{trials} to specify the number of trials per observation in binomial models, \code{cat} to specify the number of categories in 
#'  categorical or ordinal models, and \code{cens} to indicate censoring. Alternatively, the \code{addition} arguments can be incorporated directly into \code{formula}.
#'  See 'Formula Syntax' under 'Details' for further information.
#' @param autocor An optional \code{\link{cor.brms}} object describing the correlation structure within the response variable (i.e. the 'autocorrelation'). 
#'   See the documentation of \code{\link{cor.brms}} for a description of the available correlation structures. 
#'   Defaults to NULL, corresponding to no correlations.
#' @param partial A one sided formula of the form \code{~partial.effects} specifing the predictors that can vary between categories in non-cumulative ordinal models
#'  (i.e. in families \code{"cratio"}, \code{"sratio"}, or \code{"acat"}).
#' @param threshold A character string indicating the type of thresholds (i.e. intercepts) used in an ordinal model.
#'  \code{"flexible"} provides the standard unstructured thresholds and \code{"equidistant"} restricts the distance between consecutive thresholds to the same value.
#' @param cov.ranef A list of matrices that are proportional to the (within) covariance structure of the random effects. 
#'   The names of the matrices should correspond to columns in \code{data} that are used as grouping factors. 
#'   All levels of the grouping factor should appear as rownames of the corresponding matrix. 
#' @param ranef A flag to indicate if random effects for each level of the grouping factor(s) should be saved (default is \code{TRUE}). 
#'   Set to \code{FALSE} to save memory. The argument has no impact on the model fitting itself.
#' @param sample.prior A flag to indicate if samples from all specified proper priors should be additionally drawn. 
#'   Among others, these samples can be used to calculate Bayes factors for point hypotheses. Default is \code{FALSE}. 
#' @param fit An instance of S3 class \code{brmsfit} derived from a previous fit; defaults to \code{NA}. If \code{fit} is of class \code{brmsfit}, the compiled model associated 
#'   with the fitted result is re-used and the arguments \code{formula}, \code{data}, \code{family}, \code{prior}, \code{addition}, \code{autocor}, \code{partial}, \code{threshold},
#'  \code{cov.ranef}, and \code{ranef}, are ignored.
#' @param n.chains Number of Markov chains (default: 2)
#' @param n.iter Number of total iterations per chain (including burnin; default: 2000)
#' @param n.warmup A positive integer specifying number of warmup (aka burnin) iterations. This also specifies the number of iterations used for stepsize adaptation, 
#'   so warmup samples should not be used for inference. The number of warmup should not be larger than \code{n.iter} and the default is 500.
#' @param n.thin Thinning rate. Must be a positive integer. Set \code{n.thin > 1} to save memory and computation time if \code{n.iter} is large. Default is 1, that is no thinning.
#' @param n.cluster	Number of clusters to use to run parallel chains. Default is 1.   
#' @param inits Either \code{"random"} or \code{"0"}. If inits is \code{"random"} (the default), Stan will randomly generate initial values for parameters. 
#'   If it is \code{"0"}, all parameters are initiliazed to zero (recommended for \code{weibull} models). 
#'   Alternatively, \code{inits} can be a list of lists containing the initial values, or a function (or function name) generating initial values. 
#'   The latter options are mainly implemented for internal testing.
#' @param save.model Either \code{NULL} or a character string. In the latter case, the model code is
#'   saved in a file named after the string supplied in \code{save.model}, which may also contain the full path where to save the file.
#'   If only a name is given, the file is save in the current working directory. 
#' @param silent logical; If \code{TRUE}, most intermediate output from Stan is suppressed.
#' @param seed Positive integer. Used by \code{set.seed} to make results reproducable.  
#' @param ... Further arguments to be passed to Stan.
#' 
#' @return An object of class \code{brmsfit}, which contains the posterior samples along with many other useful information about the model.
#'  If rstan is not installed, \code{brmsfit} will not contain posterior samples.
#'  
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @details Fit a generalized linear mixed model, which incorporates both fixed-effects parameters and random effects in a linear predictor  
#'   via full bayesian inference using Stan. During warmup aka burnin phase, Stan may print out quite a few informational
#'   messages that \cr
#'   \code{"the current Metropolis proposal is about the be rejected ..."}. \cr 
#'   These messages can be ignored in nearly all cases. 
#'   Use \code{silent = TRUE} to stop these messages from being printed out. \cr
#'   
#'   \bold{Formula syntax}
#'   
#'   The \code{formula} argument accepts formulas of the following syntax: 
#'   
#'   \code{response | addition ~ fixed + (random | group)} 
#'   
#'   Multiple grouping factors each with multiple random effects are possible. Instead of | you may use || in random effects terms
#'   to prevent random effects correlations from being modeled.
#'   With the exception of \code{addition}, this is basically \code{lme4} syntax. 
#'   The optional \code{addition} term may contain multiple terms of the form \code{fun(variable)} seperated by \code{|} each providing
#'   special information on the response variable. \code{fun} can be replaced with either \code{se}, \code{weights}, \code{trials},
#'   \code{cat}, or \code{cens} (their meanings are explained below). Using the \code{addition} term in \code{formula} is equivalent
#'   to using argument \code{addition}: Instead of writing \code{fun(variable)} in \code{formula}, we may use \code{addition = list(fun = ~variable)}.
#'   
#'   For families \code{gaussian}, \code{student}, and \code{cauchy} it is possible to specify
#'   standard errors of the observation, thus allowing to perform meta-analysis. 
#'   Suppose that the variable \code{yi} contains the effect sizes from the studies and \code{sei} the 
#'   corresponding standard errors. Then, fixed and random effects meta-analyses can be conducted
#'   using the formulae \code{yi | se(sei) ~ 1} and \code{yi | se(sei) ~ 1 + (1|study)}, respectively, where 
#'   \code{study} is a variable uniquely identifying every study.
#'   If desired, meta-regressen can be performed via \code{yi | se(sei) ~ 1 + mod1 + mod2 + (1|study)} 
#'   or \cr \code{yi | se(sei) ~ 1 + mod1 + mod2 + (1 + mod1 + mod2|study)}, where
#'   \code{mod1} and \code{mod2} represent moderator variables. 
#'   
#'   For all families, weighted regression may be performed using
#'   \code{weights} in the addition part. Suppose that variable \code{wei} contains the weights and that \code{yi} is the response variable.
#'   Then, formula \code{yi | weights(wei) ~ predictors} implements a weighted regression. 
#'   
#'   For family \code{binomial}, addition may contain a variable indicating the number of trials 
#'   underlying each observation. In \code{lme4} syntax, we may write for instance 
#'   \code{cbind(success, n - success)}, which is equivalent
#'   to \code{success | trials(n)} in \code{brms} syntax. If the number of trials
#'   is constant across all observation (say \code{10}), we may also write \code{success | trials(10)}. 
#'   
#'   For family \code{categorical} and all ordinal families, \code{addition} may contain a term \code{cat(categories)} to
#'   specify the number categories for each observation, either with a variable name (e.g, \code{categories} in this example) or a single number.
#'   
#'   With the expection of \code{categorical} and ordinal families, left and right censoring can be modeled through \code{yi | cens(censored) ~ predictors}.
#'   The censoring variable (named \code{censored} in this example) should contain the values \code{'left'}, \code{'none'}, and \code{'right'}  
#'   (or equivalenty -1, 0, and 1) to indicate that the corresponding observation is left censored, not censored, or right censored. 
#' 
#'   Mutiple \code{addition} terms may be specified at the same time, for instance \cr 
#'   \code{formula = yi | se(sei) | cens(censored) ~ 1} for a censored meta-analytic model, equivalent to 
#'   \code{formula = yi ~ 1} and \code{addition = list(se = ~sei, cens = ~censored)} when using argument \code{addition}. \cr
#'   
#'   Family \code{gaussian} allows to perform multivariate (normal) regression using \code{cbind} notation. Suppose that
#'   \code{y1} and \code{y2} are response variables and \code{x} is a predictor, then \code{cbind(y1,y2) ~ x} speficies a multivariate model, where
#'   \code{x} has the same effect on \code{y1} and \code{y2}.
#'   To indicate different effects on each response variable, the word \code{trait} (which is reserved in models with mutiple responses) can be used
#'   as an additional categorical predictor. For instance, \code{cbind(y1,y2) ~ -1 + x:trait} leads to seperate effects
#'   of \code{x} on \code{y1} and \code{y2}. In this case, \code{trait} has two levels, namely \code{"y1"} and \code{"y2"}. 
#'   By default, \code{trait} is dummy-coded. 
#'   It may also be used within random effects terms, both as grouping factor or as random effect within a grouping factor. 
#'   Note that variable \code{trait} is generated internally and may not be specified in the data passed to \code{brm}. \cr
#'   
#'   
#'   \bold{Families and link functions}
#'   
#'   Family \code{gaussian} with \code{identity} link leads to linear regression. Families \code{student}, and \code{cauchy}
#'   with \code{identity} link leads to robust linear regression that is less influenced by outliers. 
#'   Families \code{poisson}, \code{negbinomial}, and \code{geometric} with \code{log} link lead to regression models for count data. 
#'   Families \code{binomial} and \code{bernoulli} with \code{logit} link leads to logistic regression and family \code{categorical} to
#'   multi-logistic regression when there are more than two possible outcomes.
#'   Families \code{cumulative}, \code{cratio} ('contiuation ratio'), \code{sratio} ('stopping ratio'), 
#'   and \code{acat} ('adjacent category') leads to ordinal regression. Families \code{gamma}, \code{weibull}, and \code{exponential}
#'   can be used (among others) for survival regression when combined with the \code{log} link.
#'   
#'   In the following, we list all possible links for each family.
#'   The families \code{gaussian}, \code{student}, and \code{cauchy} accept the links (as names) \code{identity}, \code{log}, and \code{inverse};
#'   families \code{poisson}, \code{negbinomial}, and \code{geometric} the links \code{log}, \code{identity}, and \code{sqrt}; 
#'   families \code{binomial}, \code{bernoulli}, \code{cumulative}, \code{cratio}, \code{sratio}, and \code{acat} the links \code{logit}, \code{probit}, \code{probit_approx}, and \code{cloglog};
#'   family  \code{categorical} the link \code{logit}; families \code{gamma}, \code{weibull}, and \code{exponential} the links \code{log}, \code{identity}, and \code{inverse}. 
#'   The first link mentioned for each family is the default. \cr    
#'   
#'   
#'   
#'   \bold{Prior distributions}
#'   
#'   Below, we describe the usage of the \code{prior} argument and list some common prior distributions 
#'   for parameters in \code{brms} models. 
#'   A complete overview on possible prior distributions is given in the Stan Reference Manual available at 
#'   \url{http://mc-stan.org/}.
#'   
#'   \code{brm} performs no checks if the priors are written in correct Stan language.
#'   Instead, Stan will check their correctness when the model is parsed to C++ and returns an error if they are not.
#'   Currently, there are five types of parameters in \code{brms} models, 
#'   for which the user can specify prior distributions. \cr
#'   
#'   1. Fixed effects 
#'   
#'   Every fixed (and partial) effect has its corresponding regression parameter. These parameters are named as
#'   \code{b_(fixed)}, where \code{(fixed)} represents the name of the corresponding fixed effect. 
#'   Suppose, for instance, that \code{y} is predicted by \code{x1} and \code{x2} 
#'   (i.e. \code{y ~ x1+x2} in formula syntax). 
#'   Then, \code{x1} and \code{x2} have regression parameters \code{b_x1} and \code{b_x2} respectively. 
#'   The default prior for fixed effects parameters is an improper flat prior over the reals. 
#'   Other common options are normal priors or uniform priors over a finite interval.
#'   If we want to have a normal prior with mean 0 and standard deviation 5 for \code{b_x1}, 
#'   and a uniform prior between -10 and 10 for \code{b_x2},
#'   we can specify this via \cr
#'   \code{prior = list(b_x1 = "normal(0,5)", b_x2 = "uniform(-10,10)")}. 
#'   To put the same prior (e.g. a normal prior) on all fixed effects at once, 
#'   we may write as a shortcut \code{prior = } \cr \code{list(b = "normal(0,5)")}. In addition, this
#'   leads to faster sampling in Stan, because priors can be vectorized. \cr
#'   
#'   2. Autocorrelation parameters
#'   
#'   The autocorrelation parameters currently implemented are named \code{ar} (autoregression) and \code{ma} (moving average).
#'   The default prior for autocorrelation parameters is an improper flat prior over the reals. It should be noted that \code{ar} will
#'   only take one values between -1 and 1 if the response variable is wide-sence stationay, i.e. if there is no drift in the responses. \cr
#'   
#'   3. Standard deviations of random effects
#'   
#'   Each random effect of each grouping factor has a standard deviation named
#'   \code{sd_(group)_(random)}. Consider, for instance, the formula \code{y ~ x1+x2+(1+x1|z)}.
#'   We see that the intercept as well as \code{x1} are random effects nested in the grouping factor \code{z}. 
#'   The corresponding standard deviation parameters are named as \code{sd_z_Intercept} and \code{sd_z_x1} respectively. 
#'   These parameters are restriced to be non-negative and, by default, 
#'   have a half cauchy prior with 'mean' 0 and 'standard deviation' 5. 
#'   We could make this explicit by writing \code{prior = list(sd = "cauchy(0,5)")}. 
#'   One common alternative is a uniform prior over a positive interval. \cr
#'   
#'   4. Correlations of random effects 
#'   
#'   If there is more than one random effect per grouping factor, the correlations between those random
#'   effects have to be estimated. 
#'   However, in \code{brms} models, the corresponding correlation matrix \eqn{C} does not have prior itself. 
#'   Instead, a prior is defined for the cholesky factor \eqn{L} of \eqn{C}. They are related through the equation
#'     \deqn{L * L' = C.} 
#'   The prior \code{"lkj_corr_cholesky(eta)"} with \code{eta > 0} is essentially the only prior for 
#'   cholesky factors of correlation matrices.
#'   If \code{eta = 1} (the default) all correlations matrices are equally likely a priori. If \code{eta > 1}, 
#'   extreme correlations become less likely, 
#'   whereas \code{0 < eta < 1} results in higher probabilities for extreme correlations. 
#'   The cholesky factors in \code{brms} models are named as 
#'   \code{L_(group)}, (e.g., \code{L_z} if \code{z} is the grouping factor). \cr
#'   
#'   5. Parameters for specific families 
#'   
#'   Some families need additional parameters to be estimated. 
#'   Families \code{gaussian}, \code{student}, and \code{cauchy} need the parameter \code{sigma} 
#'   to account for the standard deviation of the response variable around the regression line
#'   (not to be confused with the standard deviations of random effects). 
#'   By default, \code{sigma} has a half cauchy prior with 'mean' 0 and 'standard deviation' 5. 
#'   Furthermore, family \code{student} needs the parameter \code{nu} representing 
#'   the degrees of freedom of students t distribution. 
#'   By default, \code{nu} has prior \code{"uniform(1,100)"}. 
#'   Families \code{gamma} and \code{weibull} need the parameter \code{shape} 
#'   that has a \code{"gamma(0.01,0.01)"} prior by default. For families \code{cumulative}, \code{cratio}, \code{sratio}, 
#'   and \code{acat}, and only if \code{threshold = "equidistant"}, the parameter \code{delta} is used to model the distance
#'   between to adjacent thresholds. By default, \code{delta} has an improper flat prior over the reals. \cr
#' 
#' @examples
#' \dontrun{ 
#' ## Poisson Regression for the number of seizures in epileptic patients
#' ## using half cauchy priors for standard deviations of random effects 
#' fit_e <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|visit), 
#'            data = epilepsy, family = "poisson", prior = list(sd = "cauchy(0,2.5)"))  
#' ## generate a summary of the results
#' summary(fit_e)
#' ## plot the MCMC chains as well as the posterior distributions
#' plot(fit_e)
#' ## extract random effects standard devations, correlation and covariance matrices
#' VarCorr(fit_e)
#' ## extract random effects for each level
#' ranef(fit_e)
#'  
#' ## Ordinal regression (with family 'sratio') modeling patient's rating 
#' ## of inhaler instructions using normal priors for fixed effects parameters
#' fit_i <- brm(rating ~ treat + period + carry, data = inhaler, 
#'               family = "sratio", prior = list(b = "normal(0,5)"))
#' summary(fit_i)
#' plot(fit_i)    
#' 
#' ## Surivival Regression (with family 'weibull') modeling time between 
#' ## first and second recurrence of an infection in kidney patients
#' ## time | cens indicates which values in variable time are right censored
#' fit_k <- brm(time | cens(censored) ~ age + sex + disease, data = kidney, 
#'              family = "weibull", silent = TRUE, inits = "0")
#' summary(fit_k) 
#' plot(fit_k)    
#' 
#' ## Logisitic Regression (with family 'binomial') to illustrate model specification
#' ## variable n contains the number of trials, success contains the number of successes
#' n <- sample(1:10, 100, TRUE)
#' success <- rbinom(100, size = n, prob = 0.4)
#' x <- rnorm(100)
#' fit_b <- brm(success | trials(n) ~ x, family = "binomial")
#' summary(fit_b)
#'                                           
#' }
#' 
#' @import rstan
#' @import parallel
#' @import methods
#' @import stats   
#' @export 
brm <- function(formula, data = NULL, family = c("gaussian", "identity"), prior = list(),
                addition = NULL, autocor = NULL, partial = NULL, threshold = "flexible", cov.ranef = NULL, 
                ranef = TRUE, sample.prior = FALSE, fit = NA, 
                n.chains = 2, n.iter = 2000, n.warmup = 500, n.thin = 1, n.cluster = 1, inits = "random", 
                silent = FALSE, seed = 12345, save.model = NULL, ...) {
  
  if (n.chains %% n.cluster != 0) stop("n.chains must be a multiple of n.cluster")
  if (is.null(autocor)) autocor <- cor.arma()
  if (!is(autocor, "cor.brms")) stop("cor must be of class cor.brms")
  if (!threshold %in% c("flexible","equidistant")) 
    stop("threshold must be either flexible or equidistant")
  dots <- list(...) 
  if ("WAIC" %in% names(dots))
    warning("Argument WAIC is depricated. Just use method WAIC on the fitted model.")
  if ("predict" %in% names(dots)) 
    warning("Argument predict is depricated. Just use method predict on the fitted model.")
  dots[c("WAIC", "predict")] <- NULL
  set.seed(seed)
  
  if (is(fit, "brmsfit")) {
    x <- fit
    x$fit <- rstan::get_stanmodel(x$fit)
  }
  else {
    link <- link4family(family)
    family <- family[1]
    formula <- brm.update.formula(formula, addition = addition)
    et <- extract.time(autocor$formula)
    ee <- extract.effects(formula, family = family, partial, et$all)
    data.name <- Reduce(paste, deparse(substitute(data)))
    data <- updateData(data, family = family, effects = ee, et$group)
    x <- brmsfit(formula = formula, family = family, link = link, partial = partial,
                 data.name = data.name, autocor = autocor, prior = prior)
    x$ranef <- setNames(lapply(lapply(ee$random, brm.model.matrix, data = data), colnames), 
                        gsub("__", ":", ee$group))
    x$exclude <- exclude_pars(formula, ranef = ranef)
    x$data <- brm.data(formula, data = data, family = family, prior = prior, cov.ranef = cov.ranef,
                       autocor = autocor, partial = partial) 
    x$model <- stan.model(formula = formula, data = data, family = family, link = link, prior = prior, 
                          autocor = autocor, partial = partial, threshold = threshold, 
                          cov.ranef = names(cov.ranef), sample.prior = sample.prior, 
                          save.model = save.model)
    x$fit <- get_stanmodel(suppressMessages(stan(model_code = x$model, data = x$data,
               model_name = paste0(family,"(",link,") brms-model"), chains = 0)))
  }
  
  if (is.character(inits) && !inits %in% c("random", "0")) 
    inits <- get(inits, mode = "function", envir = parent.frame())
  args <- list(object = x$fit, data = x$data, pars = x$exclude, init = inits,
            iter = n.iter, warmup = n.warmup, thin = n.thin, chains = n.chains, 
            include = FALSE)
  args[names(dots)] <- dots 
  if (n.cluster > 1 || silent && n.chains > 0) {
    if (is.character(inits) || is.numeric(inits)) inits <- rep(inits, n.chains)
    cl <- makeCluster(n.cluster)
    clusterEvalQ(cl, require(rstan))
    clusterExport(cl = cl, varlist = "args", envir = environment())
    sflist <- parLapply(cl, 1:n.chains, fun = function(i) { 
      args[c("chains", "chain_id", "init")] <- list(chains = 1, chain_id = i, init = args$init[i])
      do.call(rstan::sampling, args)})
    x$fit <- rstan::sflist2stanfit(rmNULL(lapply(1:length(sflist), function(i)
      if (!is(sflist[[i]], "stanfit") || length(sflist[[i]]@sim$samples) == 0) {
        warning(paste("chain", i, "did not contain samples and was removed from the fitted model"))
        NULL
      } else sflist[[i]])))
    stopCluster(cl)
  } 
  else x$fit <- do.call(rstan::sampling, args)
  return(rename.pars(x))
}