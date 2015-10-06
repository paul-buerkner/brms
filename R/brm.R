#' Fit Bayesian Generalized Linear and Ordinal Mixed Models
#' 
#' Fit a Bayesian generalized linear or ordinal mixed model using Stan
#' 
#' @param formula An object of class "formula" (or one that can be coerced to that class): 
#'   a symbolic description of the model to be fitted. 
#'   The details of model specification are given under 'Details'.
#' @param data An optional data frame, list or environment  (or object coercible by 
#'   \code{as.data.frame} to a data frame) containing the variables in the model. If not found in data, 
#'   the variables are taken from \code{environment(formula)}, 
#'   typically the environment from which \code{brm} is called. 
#'   Although it is optional, we strongly recommend to supply a data.frame. 
#' @param family A description of the error distribution and link function to be used in the model. 
#'   This can be a vector of one or two character strings. 
#'   The first string indicates the distribution of the dependent variable (the 'family'). 
#'   Currently, the following families are supported:
#'   \code{"gaussian"}, \code{"student"}, \code{"cauchy"}, \code{"binomial"}, \code{"bernoulli"}, 
#'   \code{"categorical"}, \code{"poisson"}, \code{"negbinomial"}, \code{"geometric"}, \code{"gamma"}, 
#'   \code{"inverse.gaussian"}, \code{"exponential"}, \code{"weibull"}, 
#'   \code{hurdle_poisson}, \code{hurdle_negbinomial}, \code{hurdle_gamma}, 
#'   \code{zero_inflated_poisson}, \code{zero_inflated_negbinomial},
#'   \code{"cumulative"}, \code{"cratio"}, \code{"sratio"}, and \code{"acat"}.
#'   The second string indicates the link function, which must be supported by the family 
#'   (if not specified, default links are used). 
#'   Alternatively, a family function or the result of a call to a family function are also accepted
#'   (see \code{\link[stats:family]{family}} for help on standard family functions and 
#'   \code{\link[brms:brmsfamily]{brmsfamily}} for brms specific family functions). 
#'   Further information is provided under 'Details'.
#' @param prior One or more \code{brmsprior} objects created by function 
#'   \code{\link[brms:set_prior]{set_prior}} and combined using the \code{c} method. 
#'   A single \code{brmsprior} object may be passed without \code{c()} surrounding it. 
#'   See also  \code{\link[brms:get_prior]{get_prior}} for more help.
#' @param addition A named list of one sided formulas each containing additional information 
#'   on the response variable. The following names are allowed:
#'   \code{se} for specifying standard errors for meta-analysis, 
#'   \code{weights} to fit weighted regression models, 
#'   \code{trials} to specify the number of trials per observation in binomial models, 
#'   \code{cat} to specify the number of categories in 
#'   categorical or ordinal models, and \code{cens} to indicate censoring. 
#'   Alternatively, the \code{addition} arguments can be incorporated directly into \code{formula}.
#'   See 'Formula Syntax' under 'Details' for further information.
#' @param autocor An optional \code{\link{cor_brms}} object describing the correlation structure 
#'   within the response variable (i.e. the 'autocorrelation'). 
#'   See the documentation of \code{\link{cor_brms}} for a description 
#'   of the available correlation structures. Defaults to NULL, corresponding to no correlations.
#' @param partial A one sided formula of the form \code{~expression} specifying the predictors with 
#'   category specific effects in non-cumulative ordinal models
#'   (i.e. in families \code{"cratio"}, \code{"sratio"}, or \code{"acat"}).
#' @param threshold A character string indicating the type of thresholds (i.e. intercepts) 
#'   used in an ordinal model. \code{"flexible"} provides the standard unstructured thresholds and 
#'   \code{"equidistant"} restricts the distance between consecutive thresholds to the same value.
#' @param cov.ranef A list of matrices that are proportional to the (within) covariance structure of the random effects. 
#'   The names of the matrices should correspond to columns in \code{data} that are used as grouping factors. 
#'   All levels of the grouping factor should appear as rownames of the corresponding matrix. 
#' @param ranef A flag to indicate if random effects for each level of the grouping factor(s) 
#'   should be saved (default is \code{TRUE}). Set to \code{FALSE} to save memory. 
#'   The argument has no impact on the model fitting itself.
#' @param sample.prior A flag to indicate if samples from all specified proper priors 
#'   should be additionally drawn. Among others, these samples can be used to calculate 
#'   Bayes factors for point hypotheses. Default is \code{FALSE}. 
#' @param fit An instance of S3 class \code{brmsfit} derived from a previous fit; defaults to \code{NA}. 
#'   If \code{fit} is of class \code{brmsfit}, the compiled model associated 
#'   with the fitted result is re-used and the arguments \code{formula}, \code{data}, 
#'   \code{family}, \code{prior}, \code{addition}, \code{autocor}, \code{partial}, \code{threshold},
#'  \code{cov.ranef}, and \code{ranef}, are ignored.
#' @param inits Either \code{"random"} or \code{"0"}. If inits is \code{"random"} (the default), 
#'   Stan will randomly generate initial values for parameters. 
#'   If it is \code{"0"}, all parameters are initiliazed to zero. 
#'   This option is recommended for \code{exponential} and \code{weibull} models, as it
#'   happens that default (\code{"random"}) inits cause samples to be essentially constant. 
#'   Generally, setting \code{inits = "0"} is worth a try, if chains do not behave well.
#'   Alternatively, \code{inits} can be a list of lists containing the initial values, 
#'   or a function (or function name) generating initial values. 
#'   The latter options are mainly implemented for internal testing.
#' @param n.chains Number of Markov chains (default: 2)
#' @param n.iter Number of total iterations per chain (including burnin; default: 2000)
#' @param n.warmup A positive integer specifying number of warmup (aka burnin) iterations. 
#'   This also specifies the number of iterations used for stepsize adaptation, 
#'   so warmup samples should not be used for inference. The number of warmup should not 
#'   be larger than \code{n.iter} and the default is 500.
#' @param n.thin Thinning rate. Must be a positive integer. 
#'   Set \code{n.thin > 1} to save memory and computation time if \code{n.iter} is large. 
#'   Default is 1, that is no thinning.
#' @param n.cluster	Number of clusters to use to run parallel chains. Default is 1.   
#' @param cluster_type A character string specifying the type of cluster created by 
#'   \code{\link[parallel:makeCluster]{makeCluster}} when sampling in parallel 
#'   (i.e. when \code{n.cluster} is greater \code{1}). Default is \code{"PSOCK"} working on all platforms. 
#'   For OS X and Linux, \code{"FORK"} may be a faster and more stable option, 
#'   but it does not work on Windows.
#' @param save.model Either \code{NULL} or a character string. In the latter case, the model code is
#'   saved in a file named after the string supplied in \code{save.model}, 
#'   which may also contain the full path where to save the file.
#'   If only a name is given, the file is save in the current working directory. 
#' @param silent logical; If \code{TRUE}, most intermediate output from Stan is suppressed.
#' @param seed Positive integer. Used by \code{set.seed} to make results reproducable.  
#' @param ... Further arguments to be passed to Stan.
#' 
#' @return An object of class \code{brmsfit}, which contains the posterior samples along 
#'   with many other useful information about the model.
#'   Use \code{methods(class = "brmsfit")} for an overview on available methods.
#'  
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @details Fit a generalized linear mixed model, which incorporates both fixed-effects parameters 
#'   and random effects in a linear predictor via full bayesian inference using Stan. 
#'   During warmup aka burnin phase, Stan may print out quite a few informational messages that \cr
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
#'   Multiple grouping factors each with multiple random effects are possible. 
#'   Instead of | you may use || in random effects terms
#'   to prevent random effects correlations from being modeled.
#'   With the exception of \code{addition}, this is basically \code{lme4} syntax. 
#'   The optional \code{addition} term may contain multiple terms of the form \code{fun(variable)} 
#'   seperated by \code{|} each providing special information on the response variable. 
#'   \code{fun} can be replaced with either \code{se}, \code{weights}, \code{trials},
#'   \code{cat}, \code{cens}, or \code{trunc} (their meanings are explained below). 
#'   Using the \code{addition} term in \code{formula} is equivalent
#'   to using argument \code{addition}: Instead of writing \code{fun(variable)} in \code{formula}, 
#'   we may use \code{addition = list(fun = ~variable)}.
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
#'   \code{weights} in the addition part. Suppose that variable \code{wei} contains the weights 
#'   and that \code{yi} is the response variable. Then, formula \code{yi | weights(wei) ~ predictors} 
#'   implements a weighted regression. 
#'   
#'   For family \code{binomial}, addition may contain a variable indicating the number of trials 
#'   underlying each observation. In \code{lme4} syntax, we may write for instance 
#'   \code{cbind(success, n - success)}, which is equivalent
#'   to \code{success | trials(n)} in \code{brms} syntax. If the number of trials
#'   is constant across all observation (say \code{10}), we may also write \code{success | trials(10)}. 
#'   
#'   For family \code{categorical} and all ordinal families, \code{addition} may contain a term \code{cat(number)} to
#'   specify the number categories (e.g, \code{cat(7)}. 
#'   If not given, the number of categories is calculated from the data.
#'   
#'   With the expection of \code{categorical} and ordinal families, left and right censoring 
#'   can be modeled through \code{yi | cens(censored) ~ predictors}.
#'   The censoring variable (named \code{censored} in this example) should 
#'   contain the values \code{'left'}, \code{'none'}, and \code{'right'}  
#'   (or equivalenty -1, 0, and 1) to indicate that the corresponding observation is 
#'   left censored, not censored, or right censored. 
#'   
#'   With the expection of \code{categorical} and ordinal families, the response 
#'   distribution can be truncated using the \code{trunc} function in the addition part.
#'   If the response variable is truncated between, say, 0 and 100, we can specify this via
#'   \code{yi | trunc(lb = 0, ub = 100) ~ predictors}. Defining only one of the two arguments 
#'   in \code{trunc} leads to one-sided truncation.
#' 
#'   Mutiple \code{addition} terms may be specified at the same time, for instance \cr 
#'   \code{formula = yi | se(sei) | cens(censored) ~ 1} for a censored meta-analytic model, equivalent to 
#'   \code{formula = yi ~ 1} and \code{addition = list(se = ~sei, cens = ~censored)} 
#'   when using argument \code{addition}. \cr
#'   
#'   Family \code{gaussian} allows to perform multivariate (normal) regression using \code{cbind} notation. 
#'   Suppose that \code{y1} and \code{y2} are response variables and \code{x} is a predictor, 
#'   then \code{cbind(y1,y2) ~ x} speficies a multivariate model, 
#'   where \code{x} has the same effect on \code{y1} and \code{y2}.
#'   To indicate different effects on each response variable, the word \code{trait} 
#'   (which is reserved in multivariate models) can be used as an additional categorical predictor. 
#'   For instance, \code{cbind(y1,y2) ~ 0 + x:trait} leads to seperate effects
#'   of \code{x} on \code{y1} and \code{y2}. 
#'   In this case, \code{trait} has two levels, namely \code{"y1"} and \code{"y2"}. 
#'   By default, \code{trait} is dummy-coded. 
#'   It may also be used within random effects terms, both as grouping factor or 
#'   as random effect within a grouping factor. Note that variable \code{trait} is generated 
#'   internally and may not be specified in the data passed to \code{brm}. \cr
#'   
#'   Zero-inflated and hurdle families are bivariate and also make use of the special internal
#'   variable \code{trait} having two levels in this case. 
#'   However, only the actual response must be specified in \code{formula}, 
#'   as the second response variable used for the zero-inflation / hurdle part 
#'   is internally generated.
#'   A \code{formula} for this type of models may, for instance, look like this: 
#'   \code{y ~ 0 + trait * (x1 + x2) + (0 + trait | g)}. In this example, the fixed effects
#'   \code{x1} and \code{x1} influence the zero-inflation / hurdle part differently
#'   than the actual response part as indicated by their interaction with \code{trait}.
#'   In addition, a random effect of \code{trait} was added while the random intercept 
#'   was removed leading to the estimation of two random effects, 
#'   one for the zero-inflation / hurdle part and one for the actual response. 
#'   In the example above, the correlation between the two random effects will also be estimated.
#'   
#'   \bold{Families and link functions}
#'   
#'   Family \code{gaussian} with \code{identity} link leads to linear regression. 
#'   Families \code{student}, and \code{cauchy} with \code{identity} link leads to 
#'   robust linear regression that is less influenced by outliers. 
#'   Families \code{poisson}, \code{negbinomial}, and \code{geometric} with \code{log} link lead to 
#'   regression models for count data. 
#'   Families \code{binomial} and \code{bernoulli} with \code{logit} link leads to 
#'   logistic regression and family \code{categorical} to multi-logistic regression 
#'   when there are more than two possible outcomes.
#'   Families \code{cumulative}, \code{cratio} ('contiuation ratio'), \code{sratio} ('stopping ratio'), 
#'   and \code{acat} ('adjacent category') leads to ordinal regression. Families \code{gamma}, 
#'   \code{weibull}, \code{exponential}, and \code{inverse.gaussian} can be used (among others) 
#'   for survival regression when combined with the \code{log} link. 
#'   Families \code{hurdle_poisson}, \code{hurdle_negbinomial}, \code{hurdle_gamma}, 
#'   \code{zero_inflated_poisson}, and \code{zero_inflated_negbinomial} combined with the 
#'   \code{log} link allow to estimate zero-inflated and hurdle models. These models 
#'   can be very helpful when there are many zeros in the data that cannot be explained 
#'   by the primary distribution of the response. Family \code{hurdle_gamma} is 
#'   especially useful, as a traditional \code{gamma} model cannot be reasonably fitted for
#'   data containing zeros in the response. 
#'   
#'   In the following, we list all possible links for each family.
#'   The families \code{gaussian}, \code{student}, and \code{cauchy} accept the links (as names) 
#'   \code{identity}, \code{log}, and \code{inverse};
#'   families \code{poisson}, \code{negbinomial}, and \code{geometric} the links 
#'   \code{log}, \code{identity}, and \code{sqrt}; 
#'   families \code{binomial}, \code{bernoulli}, \code{cumulative}, \code{cratio}, \code{sratio}, 
#'   and \code{acat} the links \code{logit}, \code{probit}, \code{probit_approx}, 
#'   \code{cloglog}, and \code{cauchit}; 
#'   family \code{categorical} the link \code{logit}; families \code{gamma}, \code{weibull}, 
#'   and \code{exponential} the links \code{log}, \code{identity}, and \code{inverse};
#'   family \code{inverse.gaussian} the links \code{1/mu^2}, \code{inverse}, \code{identity} 
#'   and \code{log}; families \code{hurdle_poisson}, \code{hurdle_negbinomial},
#'   \code{hurdle_gamma}, \code{zero_inflated_poisson}, and
#'   \code{zero_inflated_negbinomial} the link \code{log}. 
#'   The first link mentioned for each family is the default.     
#'   
#'   Please note that when calling the \code{\link[stats:family]{Gamma}} family function, 
#'   the default link will be \code{inverse} not \code{log}. 
#'   Also, the \code{probit_approx} link cannot be used when calling the
#'   \code{\link[stats:family]{binomial}} family function. 
#'   
#'   The current implementation of \code{inverse.gaussian} models has some 
#'   convergence problems and requires carefully chosen prior distributions 
#'   to work efficiently. For this reason, we currently do not recommend
#'   to use the \code{inverse.gaussian} family, unless you really feel
#'   that your data requires exactly this type of model. \cr
#'   
#'   \bold{Prior distributions}
#'   
#'   As of \pkg{brms} 0.5.0, priors should be specified using the 
#'   \code{\link[brms:set_prior]{set_prior}} function. 
#'   Its documentation contains detailed information on how to correctly specify priors. 
#'   To find out on which parameters or parameter classes priors can be defined, 
#'   use \code{\link[brms:get_prior]{get_prior}}.
#'   
#' @examples
#' \dontrun{ 
#' ## Poisson Regression for the number of seizures in epileptic patients
#' ## using student_t priors for fixed effects 
#' ## and half cauchy priors for standard deviations of random effects 
#' fit_e <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|visit), 
#'            data = epilepsy, family = "poisson", 
#'            prior = c(set_prior("student_t(5,0,10)", class = "b"),
#'                      set_prior("cauchy(0,2)", class = "sd")))
#' ## generate a summary of the results
#' summary(fit_e)
#' ## plot the MCMC chains as well as the posterior distributions
#' plot(fit_e)
#' ## extract random effects standard devations, correlation and covariance matrices
#' VarCorr(fit_e)
#' ## extract random effects for each level
#' ranef(fit_e)
#'  
#' ## Ordinal regression (with family \code{sratio} and \code{cloglog} link) modeling patient's rating 
#' ## of inhaler instructions using normal priors for fixed effects parameters
#' fit_i <- brm(rating ~ treat + period + carry, data = inhaler, 
#'              family = sratio("cloglog"), prior = set_prior("normal(0,5)"))
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
#' fit_b <- brm(success | trials(n) ~ x, family = binomial("probit"))
#' summary(fit_b)
#' }
#' 
#' @import rstan
#' @import parallel
#' @import methods
#' @import stats   
#' @export 
brm <- function(formula, data = NULL, family = "gaussian", 
                prior = NULL, addition = NULL, autocor = NULL, partial = NULL, 
                threshold = c("flexible", "equidistant"), cov.ranef = NULL, 
                ranef = TRUE, sample.prior = FALSE, fit = NA, inits = "random", 
                n.chains = 2, n.iter = 2000, n.warmup = 500, n.thin = 1, n.cluster = 1, 
                cluster_type = "PSOCK", silent = FALSE, seed = 12345, 
                save.model = NULL, ...) {
  
  if (n.chains %% n.cluster != 0) 
    stop("n.chains must be a multiple of n.cluster")
  if (is.null(autocor)) 
    autocor <- cor_arma()
  if (!is(autocor, "cor_brms")) 
    stop("cor must be of class cor_brms")
  threshold <- match.arg(threshold)
  
  dots <- list(...) 
  if ("WAIC" %in% names(dots))
    warning("Argument WAIC is depricated. Just use method WAIC on the fitted model.")
  if ("predict" %in% names(dots)) 
    warning("Argument predict is depricated. Just use method predict on the fitted model.")
  dots[c("WAIC", "predict")] <- NULL
  
  set.seed(seed)
  if (is(fit, "brmsfit")) {  
    x <- fit  # re-use existing model
    x$fit <- rstan::get_stanmodel(x$fit)  # extract the compiled model
    standata <- standata(x)  # compute data to be passed to Stan
  } else {  # build new model
    # see validate.R for function definitions
    family <- check_family(family) 
    link <- family$link
    family <- family$family
    formula <- update_formula(formula, addition = addition) 
    prior <- check_prior(prior, formula = formula, data = data, 
                         family = family, link = link, 
                         autocor = autocor, partial = partial, 
                         threshold = threshold) 
    et <- extract_time(autocor$formula)  
    ee <- extract_effects(formula, family = family, partial, et$all)
    data.name <- Reduce(paste, deparse(substitute(data)))
    
    # initialize S3 object
    x <- brmsfit(formula = formula, family = family, link = link, 
                 partial = partial, data.name = data.name, 
                 autocor = autocor, prior = prior, 
                 cov.ranef = cov.ranef)  
    x$data <- update_data(data, family = family, effects = ee, et$group)  # see data.R
    x$ranef <- gather_ranef(effects = ee, data = x$data)  # see validate.R
    x$exclude <- exclude_pars(formula, ranef = ranef)  # see validate.R
    x$model <- stan_model(formula = formula, data = x$data, 
                          family = family, link = link, 
                          prior = prior,  autocor = autocor, 
                          partial = partial, threshold = threshold, 
                          cov.ranef = cov.ranef, sample.prior = sample.prior, 
                          save.model = save.model)  # see stan.R
    standata <- standata(x)  # compute data to be passed to Stan
    message("Compiling the C++ model")
    x$fit <- rstan::stanc(model_code = x$model,
                          model_name = paste0(family,"(",link,") brms-model"))
    x$fit <- rstan::stan_model(stanc_ret = x$fit) 
  }
  
  if (is.character(inits) && !inits %in% c("random", "0")) 
    inits <- get(inits, mode = "function", envir = parent.frame()) 
  if (x$family %in% c("exponential", "weibull") && inits == "random") {
    warning(paste("Families exponential and weibull may not work well",
                   "with default initial values. \n",
                   " It is thus recommended to set inits = '0'"))
  }
  if (x$family == "inverse.gaussian") {
    warning(paste("inverse gaussian models require carefully chosen prior distributions",
                  "to ensure convergence of the chains"))
  }
  if (x$link == "sqrt") {
    warning(paste(x$family, "model with sqrt link may not be uniquely identified"))
  }
  
  # arguments to be passed to stan
  args <- list(object = x$fit, data = standata, pars = x$exclude, 
               init = inits,  iter = n.iter, warmup = n.warmup, 
               thin = n.thin, chains = n.chains, include = FALSE)  
  args[names(dots)] <- dots 
  
  if (n.cluster > 1 || silent && n.chains > 0) {  # sample in parallel
    message("Start sampling")
    if (is.character(args$init) || is.numeric(args$init)) 
      args$init <- rep(args$init, n.chains)
    cl <- makeCluster(n.cluster, type = cluster_type)
    on.exit(stopCluster(cl))  # close all clusters when exiting brm
    clusterExport(cl = cl, varlist = "args", envir = environment())
    clusterEvalQ(cl, require(rstan))
    run_chain <- function(i) {
      args$chains <- 1L
      args$chain_id <- i
      args$init <- args$init[i]
      Sys.sleep(0.5 * i)
      do.call(rstan::sampling, args = args)
    }
    sflist <- parLapply(cl, X = 1:n.chains, run_chain)
    # remove chains that failed to run correctly; see validate.R
    x$fit <- lapply(1:length(sflist), remove_chains, sflist = sflist)  
    x$fit <- rstan::sflist2stanfit(rmNULL(x$fit))
  } else {  # do not sample in parallel
    x$fit <- do.call(rstan::sampling, args = args)
  }
  return(rename_pars(x))  # see rename.R
}