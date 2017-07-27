#' Fit Bayesian Generalized (Non-)Linear Multilevel Models
#' 
#' Fit Bayesian generalized (non-)linear multilevel models 
#' using Stan for full Bayesian inference. A wide range of distributions 
#' and link functions are supported, allowing users to fit -- among others -- 
#' linear, robust linear, count data, survival, response times, ordinal, 
#' zero-inflated, hurdle, and even self-defined mixture models all in a 
#' multilevel context. Further modeling options include non-linear and 
#' smooth terms, auto-correlation structures, censored data, meta-analytic 
#' standard errors, and quite a few more. In addition, all parameters of the 
#' response distribution can be predicted in order to perform distributional 
#' regression. Prior specifications are flexible and explicitly encourage 
#' users to apply prior distributions that actually reflect their beliefs.
#' In addition, model fit can easily be assessed and compared with
#' posterior predictive checks and leave-one-out cross-validation.
#' 
#' @param formula An object of class 
#'   \code{\link[stats:formula]{formula}} or
#'   \code{\link[brms:brmsformula]{brmsformula}}
#'   (or one that can be coerced to that classes): 
#'   A symbolic description of the model to be fitted. 
#'   The details of model specification are explained in 
#'   \code{\link[brms:brmsformula]{brmsformula}}.
#' @param data An object of class \code{data.frame} 
#'   (or one that can be coerced to that class) 
#'   containing data of all variables used in the model.
#' @param family A description of the response distribution and link function 
#'   to be used in the model. This can be a family function, 
#'   a call to a family function or a character string naming the family.
#'   Every family function has a \code{link} argument allowing to specify
#'   the link function to be applied on the response variable.
#'   If not specified, default links are used.
#'   For details of supported families see 
#'   \code{\link[brms:brmsfamily]{brmsfamily}}.
#'   By default, a linear \code{gaussian} model is applied.
#' @param prior One or more \code{brmsprior} objects created by
#'   \code{\link[brms:set_prior]{set_prior}} or related functions 
#'   and combined using the \code{c} method. A single \code{brmsprior} 
#'   object may be passed without \code{c()} surrounding it. 
#'   See also  \code{\link[brms:get_prior]{get_prior}} for more help.
#' @param autocor An optional \code{\link{cor_brms}} object describing 
#'   the correlation structure within the response variable 
#'   (i.e., the 'autocorrelation'). 
#'   See the documentation of \code{\link{cor_brms}} for a description 
#'   of the available correlation structures. Defaults to \code{NULL}, 
#'   corresponding to no correlations.
#' @param nonlinear (Deprecated) An optional list of formulas, specifying 
#'   linear models for non-linear parameters. If \code{NULL} (the default)
#'   \code{formula} is treated as an ordinary formula. 
#'   If not \code{NULL}, \code{formula} is treated as a non-linear model
#'   and \code{nonlinear} should contain a formula for each non-linear 
#'   parameter, which has the parameter on the left hand side and its
#'   linear predictor on the right hand side.
#'   Alternatively, it can be a single formula with all non-linear
#'   parameters on the left hand side (separated by a \code{+}) and a
#'   common linear predictor on the right hand side.
#'   As of \pkg{brms} 1.4.0, we recommend specifying non-linear
#'   parameters directly within \code{formula}.
#' @param threshold (Deprecated) A character string indicating the type 
#'   of thresholds (i.e. intercepts) used in an ordinal model. 
#'   \code{"flexible"} provides the standard unstructured thresholds and 
#'   \code{"equidistant"} restricts the distance between 
#'   consecutive thresholds to the same value.
#'   As of \pkg{brms} 1.8.0, we recommend specifying threshold
#'   directly within the ordinal family functions.
#' @param sparse Logical; indicates whether the population-level 
#'   design matrix should be treated as sparse (defaults to \code{FALSE}). 
#'   For design matrices with many zeros, this can considerably 
#'   reduce required memory. Sampling speed is currently not 
#'   improved or even slightly decreased.
#' @param cov_ranef A list of matrices that are proportional to the 
#'   (within) covariance structure of the group-level effects. 
#'   The names of the matrices should correspond to columns 
#'   in \code{data} that are used as grouping factors. 
#'   All levels of the grouping factor should appear as rownames 
#'   of the corresponding matrix. This argument can be used,
#'   among others to model pedigrees and phylogenetic effects.
#'   See \code{vignette("brms_phylogenetics")} for more details.
#' @param save_ranef A flag to indicate if group-level effects 
#'   for each level of the grouping factor(s) 
#'   should be saved (default is \code{TRUE}). 
#'   Set to \code{FALSE} to save memory. 
#'   The argument has no impact on the model fitting itself.
#'   A deprecated alias is \code{ranef}.
#' @param save_mevars A flag to indicate if samples
#'   of noise-free variables obtained by using \code{me} terms
#'   should be saved (default is \code{FALSE}).
#'   Saving these samples allows to use methods such as
#'   \code{predict} with the noise-free variables but 
#'   leads to very large \R objects even for models
#'   of moderate size and complexity.
#' @param sample_prior Indicate if samples from all specified 
#'   proper priors should be drawn additionally to the posterior samples
#'   (defaults to \code{"no"}). Among others, these samples can be used 
#'   to calculate Bayes factors for point hypotheses. 
#'   If set to \code{"only"}, samples are drawn solely from
#'   the priors ignoring the likelihood. In this case, 
#'   all parameters must have proper priors.
#' @param knots Optional list containing user specified knot values to be 
#'   used for basis construction of smoothing terms. 
#'   See \code{\link[mgcv:gamm]{gamm}} for more details.
#' @param stan_funs An optional character string containing self-defined 
#'   \pkg{Stan} functions, which will be included in the functions block 
#'   of the generated \pkg{Stan} code. 
#'   Note that these functions must additionally be defined 
#'   as \emph{vectorized} \R functions in the global environment for 
#'   various post-processing methods to work on the returned model object.
#' @param fit An instance of S3 class \code{brmsfit} derived from a previous fit; 
#'   defaults to \code{NA}. 
#'   If \code{fit} is of class \code{brmsfit}, the compiled model associated 
#'   with the fitted result is re-used and all arguments 
#'   modifying the model code or data are ignored.
#'   It is not recommended to use this argument directly, but to call 
#'   the \code{\link[brms:update.brmsfit]{update}} method, instead.
#' @param inits Either \code{"random"} or \code{"0"}. 
#'   If inits is \code{"random"} (the default), 
#'   Stan will randomly generate initial values for parameters. 
#'   If it is \code{"0"}, all parameters are initiliazed to zero. 
#'   This option is recommended for \code{exponential} and \code{weibull} models, 
#'   as it happens that default (\code{"random"}) inits cause samples 
#'   to be essentially constant. 
#'   Generally, setting \code{inits = "0"} is worth a try, if chains do not behave well.
#'   Alternatively, \code{inits} can be a list of lists containing 
#'   the initial values, or a function (or function name) generating initial values. 
#'   The latter options are mainly implemented for internal testing.
#' @param chains Number of Markov chains (defaults to 4). 
#' @param iter Number of total iterations per chain (including warmup; defaults to 2000).
#' @param warmup A positive integer specifying number of warmup (aka burnin) iterations. 
#'   This also specifies the number of iterations used for stepsize adaptation, 
#'   so warmup samples should not be used for inference. The number of warmup should not 
#'   be larger than \code{iter} and the default is \code{iter/2}.
#' @param thin Thinning rate. Must be a positive integer. 
#'   Set \code{thin > 1} to save memory and computation time if \code{iter} is large. 
#' @param cores	Number of cores to use when executing the chains in parallel, 
#'   which defaults to 1 but we recommend setting the \code{mc.cores} option 
#'   to be as many processors as the hardware and RAM allow (up to the number of chains).
#'   For non-Windows OS in non-interactive \R sessions, forking is used
#'   instead of PSOCK clusters. A deprecated alias is \code{cluster}.
#' @param algorithm Character string indicating the estimation approach to use. 
#'   Can be \code{"sampling"} for MCMC (the default), \code{"meanfield"} for
#'   variational inference with independent normal distributions, or
#'   \code{"fullrank"} for variational inference with a multivariate normal
#'   distribution.
#' @param control A named \code{list} of parameters to control the sampler's behavior. 
#'   It defaults to \code{NULL} so all the default values are used. 
#'   The most important control parameters are discussed in the 'Details'
#'   section below. For a comprehensive overview see \code{\link[rstan:stan]{stan}}.
#' @param future Logical; If \code{TRUE}, the \pkg{\link[future:future]{future}}
#'   package is used for parallel execution of the chains and argument \code{cores}
#'   will be ignored. Can be set globally for the current \R session via the 
#'   \code{future} option. The execution type is controlled via 
#'   \code{\link[future:plan]{plan}} (see the examples section below).
#' @param silent logical; If \code{TRUE} (the default), most of the
#'   informational messages of compiler and sampler are suppressed.
#'   The actual sampling progress is still printed. 
#'   Set \code{refresh = 0} to turn this off as well.
#' @param seed Used by \code{set.seed} to make results reproducable.
#'   Be aware that \code{brm} resets the seed to the value specified
#'   in \code{seed} (default: \code{12345}) every time it is run.
#'   If you want to use different seeds per run, use, for instance,
#'   \code{seed = sample(1e+7, size = 1)}. Be aware that generally, 
#'   the seed also affects subsequently called functions (such as 
#'   \code{predict}), which make use of the random number generator of \R.
#' @param save_model Either \code{NULL} or a character string. 
#'   In the latter case, the model code is
#'   saved in a file named after the string supplied in \code{save_model}, 
#'   which may also contain the full path where to save the file.
#'   If only a name is given, the file is saved in the current working directory.
#' @param save_dso Logical, defaulting to \code{TRUE}, indicating whether 
#'   the dynamic shared object (DSO) compiled from the C++ code for the model 
#'   will be saved or not. If \code{TRUE}, we can draw samples from the same 
#'   model in another \R session using the saved DSO 
#'   (i.e., without compiling the C++ code again).
#' @param ... Further arguments to be passed to Stan.
#' 
#' @return An object of class \code{brmsfit}, which contains the posterior samples along 
#'   with many other useful information about the model.
#'   Use \code{methods(class = "brmsfit")} for an overview on available methods.
#'  
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @details Fit a generalized (non-)linear multilevel model
#'   via full Bayesian inference using Stan. A general overview is provided 
#'   in the vignettes \code{vignette("brms_overview")} and 
#'   \code{vignette("brms_multilevel")}. For a full list of available 
#'   vignettes see \code{vignette(package = "brms")}.
#'  
#'   \bold{Formula syntax of brms models}
#'   
#'   Details of the formula syntax applied in \pkg{brms} 
#'   can be found in \code{\link[brms:brmsformula]{brmsformula}}.
#'   
#'   \bold{Families and link functions}
#'   
#'   Details of families supported by \pkg{brms} 
#'   can be found in \code{\link[brms:brmsfamily]{brmsfamily}}.
#'   
#'   \bold{Prior distributions}
#'   
#'   Priors should be specified using the 
#'   \code{\link[brms:set_prior]{set_prior}} function. 
#'   Its documentation contains detailed information 
#'   on how to correctly specify priors. To find out on 
#'   which parameters or parameter classes priors can be defined, 
#'   use \code{\link[brms:get_prior]{get_prior}}.
#'   Default priors are chosen to be non or very weakly informative 
#'   so that their influence on the results will be negligable and
#'   you don't have to worry about them.
#'   However, after getting more familiar with Bayesian statistics, 
#'   I recommend you to start thinking about reasonable informative
#'   priors for your model parameters: Nearly always, there is at least some
#'   prior information available that can be used to improve your inference.
#'   
#'   \bold{Adjusting the sampling behavior of \pkg{Stan}}
#'   
#'   In addition to choosing the number of iterations, warmup samples, 
#'   and chains, users can control the behavior of the NUTS sampler, 
#'   by using the \code{control} argument.
#'   The most important reason to use \code{control} is to decrease 
#'   (or eliminate at best) the number of divergent transitions
#'   that cause a bias in the obtained posterior samples. 
#'   Whenever you see the warning
#'   "There were x divergent transitions after warmup." 
#'   you should really think about increasing \code{adapt_delta}.
#'   To do this, write \code{control = list(adapt_delta = <x>)}, 
#'   where \code{<x>} should usually be value between \code{0.8} 
#'   (current default) and \code{1}. Increasing \code{adapt_delta} 
#'   will slow down the sampler but will decrease the number of 
#'   divergent transitions threatening the validity of your 
#'   posterior samples.
#'   
#'   Another problem arises when the depth of the tree being evaluated 
#'   in each iteration is exceeded. This is less common than having 
#'   divergent transitions, but may also bias the posterior samples.
#'   When it happens, \pkg{Stan} will throw out a warning suggesting 
#'   to increase \code{max_treedepth}, which can be accomplished by 
#'   writing \code{control = list(max_treedepth = <x>)} with a positive 
#'   integer \code{<x>} that should usually be larger than the current 
#'   default of \code{10}. For more details on the \code{control} argument 
#'   see \code{\link[rstan:stan]{stan}}.
#'   
#' @seealso
#'   \code{\link[brms:brms]{brms}}, 
#'   \code{\link[brms:brmsformula]{brmsformula}}, 
#'   \code{\link[brms:brmsfamily]{brmsfamily}},
#'   \code{\link[brms:brmsfit-class]{brmsfit}}
#'   
#' @examples
#' \dontrun{ 
#' ## Poisson regression for the number of seizures in epileptic patients
#' ## using student_t priors for population-level effects 
#' ## and half cauchy priors for standard deviations of group-level effects 
#' fit1 <- brm(count ~ log_Age_c + log_Base4_c * Trt_c  
#'               + (1|patient) + (1|obs), 
#'             data = epilepsy, family = poisson(), 
#'             prior = c(prior(student_t(5,0,10), class = b),
#'                       prior(cauchy(0,2), class = sd)))
#' ## generate a summary of the results
#' summary(fit1)
#' ## plot the MCMC chains as well as the posterior distributions
#' plot(fit1, ask = FALSE)
#' ## extract random effects standard devations and covariance matrices
#' VarCorr(fit1)
#' ## extract group specific effects of each level
#' ranef(fit1)
#' ## predict responses based on the fitted model
#' head(predict(fit1))  
#' ## plot marginal effects of each predictor
#' plot(marginal_effects(fit1), ask = FALSE)
#' ## investigate model fit
#' WAIC(fit1)
#' pp_check(fit1)
#'  
#' ## Ordinal regression modeling patient's rating of inhaler instructions 
#' ## category specific effects are estimated for variable 'treat'
#' fit2 <- brm(rating ~ period + carry + cs(treat), 
#'             data = inhaler, family = sratio("cloglog"), 
#'             prior = set_prior("normal(0,5)"), chains = 2)
#' summary(fit2)
#' plot(fit2, ask = FALSE) 
#' WAIC(fit2)   
#' head(predict(fit2))
#' 
#' ## Survival regression modeling the time between the first 
#' ## and second recurrence of an infection in kidney patients.
#' fit3 <- brm(time | cens(censored) ~ age * sex + disease + (1|patient), 
#'             data = kidney, family = lognormal())
#' summary(fit3) 
#' plot(fit3, ask = FALSE)
#' plot(marginal_effects(fit3), ask = FALSE)   
#' 
#' ## Probit regression using the binomial family
#' n <- sample(1:10, 100, TRUE)  # number of trials
#' success <- rbinom(100, size = n, prob = 0.4)
#' x <- rnorm(100)
#' data4 <- data.frame(n, success, x)
#' fit4 <- brm(success | trials(n) ~ x, data = data4,
#'             family = binomial("probit"))
#' summary(fit4)
#' 
#' ## Simple non-linear gaussian model
#' x <- rnorm(100)
#' y <- rnorm(100, mean = 2 - 1.5^x, sd = 1)
#' data5 <- data.frame(x, y)
#' fit5 <- brm(bf(y ~ a1 - a2^x, a1 + a2 ~ 1, nl = TRUE),  
#'             data = data5,
#'             prior = c(prior(normal(0, 2), nlpar = a1),
#'                       prior(normal(0, 2), nlpar = a2)))
#' summary(fit5)
#' plot(marginal_effects(fit5), ask = FALSE)
#' 
#' ## Normal model with heterogeneous variances
#' data_het <- data.frame(y = c(rnorm(50), rnorm(50, 1, 2)),
#'                        x = factor(rep(c("a", "b"), each = 50)))
#' fit6 <- brm(bf(y ~ x, sigma ~ 0 + x), data = data_het)
#' summary(fit6)
#' plot(fit6)
#' marginal_effects(fit6)
#' # extract estimated residual SDs of both groups
#' sigmas <- exp(posterior_samples(fit6, "^b_sigma_"))
#' ggplot(stack(sigmas), aes(values)) + 
#'   geom_density(aes(fill = ind))
#'   
#' ## Quantile regression predicting the 25%-quantile
#' fit7 <- brm(bf(y ~ x, quantile = 0.25), data = data_het, 
#'             family = asym_laplace())
#' summary(fit7)
#' marginal_effects(fit7)
#' 
#' ## use the future package for parallelization
#' library(future)
#' plan(multiprocess)
#' fit7 <- update(fit7, future = TRUE)
#' }
#' 
#' @import rstan
#' @import Rcpp
#' @import parallel
#' @import methods
#' @import stats   
#' @export 
brm <- function(formula, data, family = gaussian(), prior = NULL, 
                autocor = NULL, nonlinear = NULL, 
                threshold = c("flexible", "equidistant"), 
                cov_ranef = NULL, save_ranef = TRUE, save_mevars = FALSE, 
                sparse = FALSE, sample_prior = c("no", "yes", "only"), 
                knots = NULL, stan_funs = NULL, fit = NA, inits = "random", 
                chains = 4, iter = 2000, warmup = floor(iter / 2),
                thin = 1, cores = getOption("mc.cores", 1L), control = NULL,
                algorithm = c("sampling", "meanfield", "fullrank"),
                future = getOption("future", FALSE), silent = TRUE, 
                seed = 12345, save_model = NULL, save_dso = TRUE, ...) {
  
  dots <- list(...) 
  # use deprecated arguments if specified
  iter <- use_alias(iter, dots[["n.iter"]])
  warmup <- use_alias(warmup, dots[["n.warmup"]])
  thin <- use_alias(thin, dots[["n.thin"]])
  chains <- use_alias(chains, dots[["n.chains"]])
  cores <- use_alias(cores, dots[["cluster"]])
  cov_ranef <- use_alias(cov_ranef, dots[["cov.ranef"]])
  save_ranef <- use_alias(save_ranef, dots[["ranef"]])
  sample_prior <- use_alias(sample_prior, dots[["sample.prior"]])
  save_model <- use_alias(save_model, dots[["save.model"]])
  if (!is.null(dots[["cluster_type"]])) {
    warning2("Argument 'cluster_type' is deprecated and unused.\n",
             "Forking is now automatically applied when appropriate.")
  }
  dots[deprecated_brm_args()] <- NULL
  autocor <- check_autocor(autocor)
  algorithm <- match.arg(algorithm)
  
  testmode <- dots$testmode
  dots$testmode <- NULL
  if (is(fit, "brmsfit")) {
    # re-use existing model
    x <- fit
    # compute data to be passed to Stan
    sdata <- standata(x, is_newdata = dots$is_newdata)
    dots$is_newdata <- NULL
    # extract the compiled model
    x$fit <- rstan::get_stanmodel(x$fit)
  } else {  
    # build new model
    family <- check_family(family, threshold = threshold)
    formula <- amend_formula(
      formula, data = data, family = family, nonlinear = nonlinear
    )
    family <- formula$family
    check_brm_input(nlist(family))
    bterms <- parse_bf(formula, autocor = autocor)
    if (is.null(dots$data.name)) {
      data.name <- substr(Reduce(paste, deparse(substitute(data))), 1, 50)
    } else {
      data.name <- dots$data.name
      dots$data.name <- NULL
    }
    data <- update_data(data, bterms = bterms)
    prior <- check_prior(
      prior, formula = formula, data = data, family = family, 
      sample_prior = sample_prior, autocor = autocor, warn = TRUE
    )
    # initialize S3 object
    x <- brmsfit(
      formula = formula, family = family, data = data, 
      data.name = data.name, prior = prior, 
      autocor = autocor, cov_ranef = cov_ranef, 
      algorithm = algorithm
    )
    x$ranef <- tidy_ranef(bterms, data = x$data)  
    x$exclude <- exclude_pars(
      bterms, data = x$data, ranef = x$ranef, 
      save_ranef = save_ranef, save_mevars = save_mevars
    )
    x$model <- make_stancode(
      formula = formula, data = data, family = family, 
      prior = prior, autocor = autocor, 
      sparse = sparse, cov_ranef = cov_ranef,
      sample_prior = sample_prior, knots = knots, 
      stan_funs = stan_funs, save_model = save_model, 
      brm_call = TRUE, silent = silent
    )
    # generate Stan data before compiling the model to avoid
    # unnecessary compilations in case of invalid data
    sdata <- standata(x, newdata = dots$is_newdata)
    message("Compiling the C++ model")
    x$fit <- eval_silent(
      rstan::stan_model(stanc_ret = x$model, save_dso = save_dso)
    )
    x$model <- x$model$model_code
  }
  
  # arguments to be passed to Stan
  if (is.character(inits) && !inits %in% c("random", "0")) {
    inits <- get(inits, mode = "function", envir = parent.frame())
  }
  args <- nlist(
    object = x$fit, data = sdata, pars = x$exclude, 
    include = FALSE, algorithm, iter
  )
  args[names(dots)] <- dots
  
  set.seed(seed)
  message("Start sampling")
  if (args$algorithm == "sampling") {
    args$algorithm <- NULL
    args <- c(args,
      nlist(init = inits, warmup, thin, control, show_messages = !silent)
    )
    if (future) {
      require_package("future")
      if (cores > 1L) {
        warning("Argument 'cores' is ignored when using 'future'.")
      }
      args$chains <- 1L
      futures <- fits <- vector("list", chains)
      for (i in seq_len(chains)) {
        args$chain_id <- i
        if (is.list(inits)) {
          args$init <- inits[i]
        }
        futures[[i]] <- future::future(
          do.call(rstan::sampling, args), packages = "rstan"
        )
      }
      for (i in seq_len(chains)) {
        fits[[i]] <- future::value(futures[[i]]) 
      }
      x$fit <- rstan::sflist2stanfit(fits)
      rm(futures, fits)
    } else {
      args <- c(args, nlist(chains, cores))
      x$fit <- do.call(rstan::sampling, args) 
    }
  } else {
    # vb does not support parallel execution
    x$fit <- do.call(rstan::vb, args)
  }
  if (!isTRUE(testmode)) {
    x <- rename_pars(x)
  }
  x
}

check_brm_input <- function(x) {
  # misc checks on brm arguments 
  # Args:
  #   x: A named list
  family <- check_family(x$family) 
  if (family$family == "inverse.gaussian") {
    warning2("Inverse gaussian models require carefully chosen ", 
             "prior distributions to ensure convergence of the chains.")
  }
  if (family$family == "geometric") {
    warning2("Family 'geometric' is deprecated. Use 'negbinomial' ", 
             "instead and fix the 'shape' parameter to 1.")
  }
  if (family$link == "sqrt") {
    warning2(family$family, " model with sqrt link may not be ", 
             "uniquely identified")
  }
  invisible(NULL)
}

deprecated_brm_args <- function() {
  # list all deprecated arguments of the brm function
  c("n.iter", "n.warmup", "n.thin", "n.chains", "cluster", "cov.ranef",
    "ranef", "sample.prior", "save.model", "cluster_type")
}
