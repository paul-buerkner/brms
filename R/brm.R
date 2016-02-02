#' Fit Bayesian Generalized Linear and Ordinal Mixed Models
#' 
#' Fit a Bayesian generalized linear or ordinal mixed model using Stan
#' 
#' @param formula An object of class "formula" (or one that can be coerced to that class): 
#'   a symbolic description of the model to be fitted. 
#'   The details of model specification are given under 'Details'.
#' @param data An optional data frame, list or environment  (or object coercible by 
#'   \code{as.data.frame} to a data frame) containing the variables in the model. 
#'   If not found in data, the variables are taken from \code{environment(formula)}, 
#'   typically the environment from which \code{brm} is called. 
#'   Although it is optional, we strongly recommend to supply a data.frame. 
#' @param family A description of the error distribution and link function 
#'   to be used in the model. This can be a family function, 
#'   a call to a family function or a character string naming the family.
#'   Currently, the following families are supported:
#'   \code{gaussian}, \code{student}, \code{cauchy} (deprecated), \code{binomial}, 
#'   \code{bernoulli}, \code{Beta}, \code{poisson}, \code{negbinomial}, 
#'   \code{geometric}, \code{Gamma}, \code{inverse.gaussian}, 
#'   \code{exponential}, \code{weibull}, \code{categorical}, \code{cumulative}, 
#'   \code{cratio}, \code{sratio}, \code{acat}, \code{hurdle_poisson}, 
#'   \code{hurdle_negbinomial}, \code{hurdle_gamma}, \code{zero_inflated_binomial},
#'   \code{zero_inflated_beta}, \code{zero_inflated_negbinomial}, 
#'   and \code{zero_inflated_poisson}.
#'   Every family function has a \code{link} argument allowing to specify
#'   the link function to be applied on the response variable.
#'   If not specified, default links are used.
#'   See \code{\link[stats:family]{family}} for help on standard family functions 
#'   and \code{\link[brms:brmsfamily]{brmsfamily}} for family functions
#'   specific to the \pkg{brms} package. 
#'   For backwards compatibility, \code{family} may also be a vector of two
#'   character strings, the first naming the family and the second naming the link.
#'   Further information is provided under 'Details'.
#' @param prior One or more \code{brmsprior} objects created by function 
#'   \code{\link[brms:set_prior]{set_prior}} and combined using the \code{c} method. 
#'   A single \code{brmsprior} object may be passed without \code{c()} surrounding it. 
#'   See also  \code{\link[brms:get_prior]{get_prior}} for more help.
#' @param addition Deprecated.
#'   All additional information on the response variable should be incorporated 
#'   directly into \code{formula}. 
#'   See 'Formula Syntax' under 'Details' for further information.
#' @param autocor An optional \code{\link{cor_brms}} object describing 
#'   the correlation structure 
#'   within the response variable (i.e. the 'autocorrelation'). 
#'   See the documentation of \code{\link{cor_brms}} for a description 
#'   of the available correlation structures. Defaults to NULL, 
#'   corresponding to no correlations.
#' @param partial A one sided formula of the form \code{~expression} 
#'   allowing to specify predictors with category specific effects 
#'   in non-cumulative ordinal models 
#'   (i.e. in families \code{cratio}, \code{sratio}, or \code{acat}).
#' @param threshold A character string indicating the type of thresholds 
#'   (i.e. intercepts) used in an ordinal model. 
#'   \code{"flexible"} provides the standard unstructured thresholds and 
#'   \code{"equidistant"} restricts the distance between 
#'   consecutive thresholds to the same value.
#' @param cov_ranef A list of matrices that are proportional to the 
#'   (within) covariance structure of the random effects. 
#'   The names of the matrices should correspond to columns 
#'   in \code{data} that are used as grouping factors. 
#'   All levels of the grouping factor should appear as rownames 
#'   of the corresponding matrix. 
#' @param ranef A flag to indicate if random effects 
#'   for each level of the grouping factor(s) 
#'   should be saved (default is \code{TRUE}). 
#'   Set to \code{FALSE} to save memory. 
#'   The argument has no impact on the model fitting itself.
#' @param sample_prior A flag to indicate if samples from all specified proper priors 
#'   should be additionally drawn. Among others, these samples can be used to calculate 
#'   Bayes factors for point hypotheses. Default is \code{FALSE}. 
#' @param fit An instance of S3 class \code{brmsfit} derived from a previous fit; 
#'   defaults to \code{NA}. 
#'   If \code{fit} is of class \code{brmsfit}, the compiled model associated 
#'   with the fitted result is re-used and all arguments 
#'   modifying the model code or data are ignored.
#' @param inits Either \code{"random"} or \code{"0"}. 
#'   If inits is \code{"random"} (the default), 
#'   Stan will randomly generate initial values for parameters. 
#'   If it is \code{"0"}, all parameters are initiliazed to zero. 
#'   This option is recommended for \code{exponential} and \code{weibull} models, 
#'   as it happens that default (\code{"random"}) inits cause samples 
#'   to be essentially constant. 
#'   Generally, setting \code{inits = "0"} is worth a try, 
#'   if chains do not behave well.
#'   Alternatively, \code{inits} can be a list of lists containing 
#'   the initial values, or a function (or function name) generating initial values. 
#'   The latter options are mainly implemented for internal testing.
#' @param chains Number of Markov chains (defaults to 2). 
#'   A deprecated alias is \code{n.chains}.
#' @param iter Number of total iterations per chain (including warmup; defaults to 2000).
#'   A deprecated alias is \code{n.iter}.
#' @param warmup A positive integer specifying number of warmup (aka burnin) iterations. 
#'   This also specifies the number of iterations used for stepsize adaptation, 
#'   so warmup samples should not be used for inference. The number of warmup should not 
#'   be larger than \code{iter} and the default is 500.
#'   A deprecated alias is \code{n.warmup}.
#' @param thin Thinning rate. Must be a positive integer. 
#'   Set \code{thin > 1} to save memory and computation time if \code{iter} is large. 
#'   Default is 1, that is no thinning. A deprecated alias is \code{n.thin}.
#' @param cluster	Number of clusters to use to run parallel chains. Default is 1.  
#'   A deprecated alias is \code{n.cluster}. To use the built-in parallel execution
#'   of \pkg{rstan}, specify argument \code{cores} instead of \code{cluster}. 
#' @param cluster_type A character string specifying the type of cluster created by 
#'   \code{\link[parallel:makeCluster]{makeCluster}} when sampling in parallel 
#'   (i.e. when \code{cluster} is greater \code{1}). 
#'   Default is \code{"PSOCK"} working on all platforms. 
#'   For OS X and Linux, \code{"FORK"} may be a faster and more stable option, 
#'   but it does not work on Windows.
#' @param save_model Either \code{NULL} or a character string. 
#'   In the latter case, the model code is
#'   saved in a file named after the string supplied in \code{save_model}, 
#'   which may also contain the full path where to save the file.
#'   If only a name is given, the file is save in the current working directory. 
#' @param algorithm Character string indicating the estimation approach to use. 
#'   Can be \code{"sampling"} for MCMC (the default), \code{"meanfield"} for
#'   variational inference with independent normal distributions, or
#'   \code{"fullrank"} for variational inference with a multivariate normal
#'   distribution.
#' @param silent logical; If \code{TRUE}, warning messages of the sampler are suppressed.
#' @param seed Positive integer. Used by \code{set.seed} to make results reproducable.  
#' @param ... Further arguments to be passed to Stan.
#' 
#' @return An object of class \code{brmsfit}, which contains the posterior samples along 
#'   with many other useful information about the model.
#'   Use \code{methods(class = "brmsfit")} for an overview on available methods.
#'  
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @details Fit a generalized linear mixed model, 
#'   which incorporates both fixed-effects parameters 
#'   and random effects in a linear predictor 
#'   via full bayesian inference using Stan. 
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
#'   \code{cat}, \code{cens}, or \code{trunc}. Their meanings are explained below. 
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
#'   For families \code{binomial} and \code{zero_inflated_binomial}, 
#'   addition should contain a variable indicating the number of trials 
#'   underlying each observation. In \code{lme4} syntax, we may write for instance 
#'   \code{cbind(success, n - success)}, which is equivalent
#'   to \code{success | trials(n)} in \code{brms} syntax. If the number of trials
#'   is constant across all observation (say \code{10}), we may also write \code{success | trials(10)}. 
#'   
#'   For family \code{categorical} and all ordinal families, 
#'   \code{addition} may contain a term \code{cat(number)} to
#'   specify the number categories (e.g, \code{cat(7)}). 
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
#'   \code{formula = yi | se(sei) | cens(censored) ~ 1} for a censored meta-analytic model. \cr
#'   
#'   For families \code{gaussian}, \code{student}, and \code{cauchy} 
#'   multivariate models may be specified using \code{cbind} notation. 
#'   Suppose that \code{y1} and \code{y2} are response variables 
#'   and \code{x} is a predictor.
#'   Then \code{cbind(y1,y2) ~ x} specifies a multivariate model, 
#'   where \code{x} has the same effect on \code{y1} and \code{y2}.
#'   To indicate different effects on each response variable, 
#'   the variable \code{trait} (which is reserved in multivariate models) 
#'   can be used as an additional categorical predictor. 
#'   For instance, \code{cbind(y1,y2) ~ 0 + x:trait} leads to seperate effects
#'   of \code{x} on \code{y1} and \code{y2}. 
#'   In this case, \code{trait} has two levels, namely \code{"y1"} and \code{"y2"}. 
#'   By default, \code{trait} is dummy-coded. 
#'   It may also be used within random effects terms, both as grouping factor or 
#'   as random effect within a grouping factor. Note that variable \code{trait} is generated 
#'   internally and may not be specified in the data passed to \code{brm}. \cr
#'   
#'   Zero-inflated and hurdle families are bivariate and also make use 
#'   of the special internal variable \code{trait} having two levels in this case. 
#'   However, only the actual response must be specified in \code{formula}, 
#'   as the second response variable used for the zero-inflation / hurdle
#'   (ZIH) part is internally generated.
#'   A \code{formula} for this type of models may, for instance, look like this: \cr
#'   \code{y ~ 0 + trait * (x1 + x2) + (0 + trait | g)}. In this example, the fixed effects
#'   \code{x1} and \code{x1} influence the ZIH part differently
#'   than the actual response part as indicated by their interaction with \code{trait}.
#'   In addition, a random effect of \code{trait} was added while the random intercept 
#'   was removed leading to the estimation of two random effects, 
#'   one for the ZIH part and one for the actual response. 
#'   In the example above, the correlation between the two random effects 
#'   will also be estimated.
#'   Sometimes, predictors should only influence the ZIH part
#'   but not the actual response (or vice versa). As this cannot be modeled
#'   with the \code{trait} variable, two other internally generated and 
#'   reserved (numeric) variables, namely \code{main} and \code{spec}, are supported.
#'   \code{main} is \code{1} for the response part and \code{0} for the
#'   ZIH part of the model. For \code{spec} it is the other way round.
#'   Suppose that \code{x1} should only influence the actual response,
#'   and \code{x2} only the ZIH process. We can write this as follows:
#'   \code{formula = y ~ 0 + main + spec + main:x1 + spec:x2}. 
#'   The main effects of \code{main} or \code{spec} serve as intercepts,
#'   while the interaction terms \code{main:x1} and \code{spec:x2} ensure
#'   that \code{x1} and \code{x2} only predict one part of the model, respectively.
#'   
#'   Using the same syntax as for zero-inflated and hurdle models, it is
#'   possible to specify multiplicative effects in family \code{bernoulli}
#'   (make sure to set argument \code{type} to \code{"2PL"}; 
#'   see \code{\link[brms:brmsfamily]{brmsfamily}} for more details).
#'   In Item Response Theory (IRT), these models are known as 2PL models.
#'   Suppose that we have the variables \code{item} and \code{person} and
#'   want to model fixed effects for items and random effects for persons.
#'   The discriminality (multiplicative effect) should depend only on the items. 
#'   We can specify this by setting
#'   \code{formula = response ~ 0 + (main + spec):item + (0 + main|person)}. 
#'   The random term \code{0 + main} ensures that \code{person} does 
#'   not influence discriminalities. Of course it is possible
#'   to predict only discriminalities by using
#'   variable \code{spec} in the model formulation. 
#'   To identify the model, multiplicative effects
#'   are estimated on the log scale. 
#'   In addition, we strongly recommend setting proper priors 
#'   on fixed effects in this case to increase sampling efficiency 
#'   (for details on priors see \code{\link[brms:set_prior]{set_prior}}).     
#'   
#'   \bold{Parameterization of the fixed effects intercept}
#'   
#'   The fixed effects intercept (if incorporated) is estimated separately 
#'   and not as part of the fixed effects parameter vector \code{b}. 
#'   This has the side effect that priors on the intercept 
#'   also have to be specified separately
#'   (see \code{\link[brms:set_prior]{set_prior}} for more details).
#'   Furthermore, to increase sampling efficiency, the fixed effects 
#'   design matrix \code{X} is centered around its column means 
#'   \code{X_means} if the intercept is incorporated. 
#'   This leads to a temporary bias in the intercept equal to 
#'   \code{<X_means, b>}, where \code{<,>} is the scalar product. 
#'   The bias is corrected after fitting the model, but be aware 
#'   that you are effectively defining a prior on the temporary
#'   intercept of the centered design matrix not on the real intercept.
#'   
#'   This behavior can be avoided by using the reserved 
#'   (and internally generated) variable \code{intercept}. 
#'   Instead of \code{y ~ x}, you may write
#'   \code{y ~ -1 + intercept + x}. This way, priors can be
#'   defined on the real intercept, directly. In addition,
#'   the intercept is just treated as an ordinary fixed effect
#'   and thus priors defined on \code{b} will also apply to it. 
#'   Note that this parameterization may be a bit less efficient
#'   than the default parameterization discussed above.  
#'   
#'   \bold{Families and link functions}
#'   
#'   Family \code{gaussian} with \code{identity} link leads to linear regression. 
#'   Families \code{student}, and \code{cauchy} with \code{identity} link leads to 
#'   robust linear regression that is less influenced by outliers. 
#'   Families \code{poisson}, \code{negbinomial}, and \code{geometric} 
#'   with \code{log} link lead to regression models for count data. 
#'   Families \code{binomial} and \code{bernoulli} with \code{logit} link leads to 
#'   logistic regression and family \code{categorical} to multi-logistic regression 
#'   when there are more than two possible outcomes.
#'   Families \code{cumulative}, \code{cratio} ('contiuation ratio'), 
#'   \code{sratio} ('stopping ratio'), and \code{acat} ('adjacent category') 
#'   leads to ordinal regression. Families \code{Gamma}, \code{weibull}, 
#'   \code{exponential}, and \code{inverse.gaussian} can be used (among others) 
#'   for survival regression when combined with the \code{log} link. 
#'   Families \code{hurdle_poisson}, \code{hurdle_negbinomial}, \code{hurdle_gamma}, 
#'   \code{zero_inflated_poisson}, and \cr
#'   \code{zero_inflated_negbinomial} combined with the 
#'   \code{log} link, and  \code{zero_inflated_binomial} with the \code{logit} link, 
#'   allow to estimate zero-inflated and hurdle models. These models 
#'   can be very helpful when there are many zeros in the data that cannot be explained 
#'   by the primary distribution of the response. Family \code{hurdle_gamma} is 
#'   especially useful, as a traditional \code{Gamma} model cannot be reasonably 
#'   fitted for data containing zeros in the response. 
#'   
#'   In the following, we list all possible links for each family.
#'   The families \code{gaussian}, \code{student}, and \code{cauchy} 
#'   accept the links (as names) \code{identity}, \code{log}, and \code{inverse};
#'   families \code{poisson}, \code{negbinomial}, and \code{geometric} the links 
#'   \code{log}, \code{identity}, and \code{sqrt}; 
#'   families \code{binomial}, \code{bernoulli}, \code{Beta},
#'   \code{cumulative}, \code{cratio}, \code{sratio}, and \code{acat} 
#'   the links \code{logit}, \code{probit}, \code{probit_approx}, 
#'   \code{cloglog}, and \code{cauchit}; 
#'   family \code{categorical} the link \code{logit}; 
#'   families \code{Gamma}, \code{weibull}, and \code{exponential} 
#'   the links \code{log}, \code{identity}, and \code{inverse};
#'   family \code{inverse.gaussian} the links \code{1/mu^2}, 
#'   \code{inverse}, \code{identity} and \code{log}; 
#'   families \code{hurdle_poisson}, \code{hurdle_negbinomial},
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
#'   Its documentation contains detailed information 
#'   on how to correctly specify priors. To find out on 
#'   which parameters or parameter classes priors can be defined, 
#'   use \code{\link[brms:get_prior]{get_prior}}. \cr
#'   
#'   \bold{Adjusting the sampling behavior of \pkg{Stan}}
#'   
#'   Despite choosing the number of iterations, chains, etc., 
#'   users can directly change the sampling behavior of the NUTS sampler, 
#'   by using the \code{control} argument (a named list), 
#'   which is passed directly to \pkg{Stan} when specified in \code{brm}. 
#'   The most important reason to use \code{control} is to decrease 
#'   (or eliminate at best) the number of divergent transitions
#'   that cause a bias in the obtained posterior samples. 
#'   Whenever you see the warning
#'   "There were x divergent transitions after warmup. 
#'   Increasing adapt_delta may help." 
#'   you should really think about increasing \code{adapt_delta}.
#'   To do this, write \code{control = list(adapt_delta = <x>)}, where \code{<x>}
#'   should usually be value between \code{0.8} (default) and \code{1}.
#'   Increasing \code{adapt_delta} will slow down the sampler but will 
#'   decrease the number of divergent transitions threatening
#'   the validity of your posterior samples. 
#'   For more details on the \code{control} argument see 
#'   \code{\link[rstan:stan]{stan}}.
#'   
#' @examples
#' \dontrun{ 
#' ## Poisson regression for the number of seizures in epileptic patients
#' ## using student_t priors for fixed effects 
#' ## and half cauchy priors for standard deviations of random effects 
#' fit1 <- brm(count ~ log_Age_c + log_Base4_c * Trt_c 
#'             + (1|patient) + (1|visit) + (1|obs), 
#'             data = epilepsy, family = poisson(), 
#'             prior = c(set_prior("student_t(5,0,10)", class = "b"),
#'                       set_prior("cauchy(0,2)", class = "sd")))
#' ## generate a summary of the results
#' summary(fit1)
#' ## plot the MCMC chains as well as the posterior distributions
#' plot(fit1)
#' ## extract random effects standard devations and covariance matrices
#' VarCorr(fit1)
#' ## extract random effects for each level
#' ranef(fit1)
#' ## predict responses based on the fitted model
#' head(predict(fit1))  
#'  
#' ## Ordinal regression modeling patient's rating 
#' ## of inhaler instructions with normal priors on fixed effects
#' fit2 <- brm(rating ~ treat + period + carry, 
#'             data = inhaler, family = sratio("cloglog"), 
#'             prior = set_prior("normal(0,5)"))
#' summary(fit2)
#' plot(fit2)    
#' 
#' ## Surivival regression (with family 'weibull') modeling time between 
#' ## first and second recurrence of an infection in kidney patients.
#' fit3 <- brm(time | cens(censored) ~ age + sex + disease + (1|patient), 
#'             data = kidney, family = weibull(), inits = "0")
#' summary(fit3) 
#' plot(fit3)    
#' 
#' ## Probit regression using the binomial family
#' n <- sample(1:10, 100, TRUE)  # number of trials
#' success <- rbinom(100, size = n, prob = 0.4)
#' x <- rnorm(100)
#' fit4 <- brm(success | trials(n) ~ x, 
#'             family = binomial("probit"))
#' summary(fit4)
#' }
#' 
#' @import rstan
#' @import parallel
#' @import methods
#' @import stats   
#' @export 
brm <- function(formula, data = NULL, family = gaussian(), 
                prior = NULL, addition = NULL, autocor = NULL, 
                partial = NULL, threshold = c("flexible", "equidistant"), 
                cov_ranef = NULL, ranef = TRUE, sample_prior = FALSE, 
                fit = NA, inits = "random", chains = 2, iter = 2000, 
                warmup = 500, thin = 1, cluster = 1, cluster_type = "PSOCK", 
                algorithm = c("sampling", "meanfield", "fullrank"),
                silent = TRUE, seed = 12345, save_model = NULL, ...) {
  
  dots <- list(...) 
  # use deprecated arguments if specified
  iter <- use_alias(iter, dots$n.iter)
  warmup <- use_alias(warmup, dots$n.warmup)
  thin <- use_alias(thin, dots$n.thin)
  chains <- use_alias(chains, dots$n.chains)
  cluster <- use_alias(cluster, dots$n.cluster)
  cov_ranef <- use_alias(cov_ranef, dots$cov.ranef)
  sample_prior <- use_alias(sample_prior, dots$sample.prior)
  save_model <- use_alias(save_model, dots$save.model)
  dots[c("n.iter", "n.warmup", "n.thin", "n.chains", "n.cluster",
         "cov.ranef", "sample.prior", "save.model")] <- NULL
  # some input checks 
  check_brm_input(nlist(family, chains, cluster, inits))
  autocor <- check_autocor(autocor)
  threshold <- match.arg(threshold)
  algorithm <- match.arg(algorithm)
  
  testmode <- dots$testmode
  dots$testmode <- NULL
  if (is(fit, "brmsfit")) {  
    x <- fit  # re-use existing model
    x$fit <- rstan::get_stanmodel(x$fit)  # extract the compiled model
    # compute data to be passed to Stan
    standata <- standata(x, is_newdata = dots$is_newdata)
    dots$is_newdata <- NULL
  } else {  # build new model
    # see validate.R and priors.R for function definitions
    family <- check_family(family)
    formula <- update_formula(formula, addition = addition, data = data) 
    prior <- check_prior(prior, formula = formula, data = data, 
                         family = family, autocor = autocor,
                         partial = partial, threshold = threshold) 
    et <- extract_time(autocor$formula)  
    ee <- extract_effects(formula, family = family, partial, et$all)
    data.name <- Reduce(paste, deparse(substitute(data)))
    
    # initialize S3 object
    x <- brmsfit(formula = formula, family = family, link = family$link, 
                 partial = partial, data.name = data.name, autocor = autocor, 
                 prior = prior, cov_ranef = cov_ranef,
                 algorithm = algorithm)  
    # see data.R
    x$data <- update_data(data, family = family, effects = ee, et$group) 
    # see validate.R
    x$ranef <- gather_ranef(ee$random, data = x$data, 
                            is_forked = is.forked(family))  
    x$exclude <- exclude_pars(formula, ranef = ranef)
    # see stan.R
    x$model <- make_stancode(formula = formula, data = data, 
                             family = family, prior = prior,  
                             autocor = autocor, partial = partial, 
                             threshold = threshold, 
                             cov_ranef = cov_ranef, 
                             sample_prior = sample_prior, 
                             save_model = save_model)
    # generate standata before compiling the model to avoid
    # unnecessary compilations in case that the data is invalid
    standata <- standata(x, newdata = dots$is_newdata)
    message("Compiling the C++ model")
    x$fit <- rstan::stanc(model_code = x$model, 
                          model_name = model_name(family))
    x$fit <- rstan::stan_model(stanc_ret = x$fit) 
  }
  
  # arguments to be passed to stan
  if (is.character(inits) && !inits %in% c("random", "0")) {
    inits <- get(inits, mode = "function", envir = parent.frame())
  }
  args <- list(object = x$fit, data = standata, pars = x$exclude, 
               include = FALSE, algorithm = algorithm)
  args[names(dots)] <- dots 
  if (algorithm == "sampling") {
    args <- c(args, init = inits, iter = iter, warmup = warmup, 
              thin = thin, chains = chains, show_messages = !silent)
  }
  
  set.seed(seed)
  if (cluster > 1) {  # sample in parallel
    message("Start sampling")
    if (is.character(args$init) || is.numeric(args$init)) 
      args$init <- rep(args$init, chains)
    cl <- makeCluster(cluster, type = cluster_type)
    on.exit(stopCluster(cl))  # close all clusters when exiting brm
    clusterExport(cl = cl, varlist = "args", envir = environment())
    clusterEvalQ(cl, require(rstan))
    run_chain <- function(i) {
      args$chains <- 1L
      args$chain_id <- i
      args$init <- args$init[i]
      Sys.sleep(0.5 * i)
      if (args$algorithm == "sampling") {
        args$algorithm <- NULL
        do.call(rstan::sampling, args = args)
      } else {
        do.call(rstan::vb, args = args)
      } 
    }
    sflist <- parLapply(cl, X = 1:chains, run_chain)
    # remove chains that failed to run correctly; see validate.R
    sflist <- rmNULL(lapply(seq_along(sflist), remove_chains, sflist = sflist))  
    if (length(sflist) == 0) {
      stop(paste("All chains failed to run correctly." ,
                 "For more detailed error reporting",
                 "fit the model in non-parallel mode."), 
           call. = FALSE)
    }
    x$fit <- rstan::sflist2stanfit(sflist)
  } else {  # do not sample in parallel
    if (args$algorithm == "sampling") {
      args$algorithm <- NULL
      x$fit <- do.call(rstan::sampling, args = args)
    } else {
      x$fit <- do.call(rstan::vb, args = args)
    } 
  }
  if (!isTRUE(testmode)) x <- rename_pars(x) # see rename.R
  x
}