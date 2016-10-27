#' Prior Definitions for \pkg{brms} Models
#'
#' Define priors for specific parameters or classes of parameters
#'
#' @param prior A character string defining a distribution in \pkg{Stan} language
#' @param class The parameter class. Defaults to \code{"b"} (fixed effects). 
#'   See 'Details' for other valid parameter classes. 
#' @param coef Name of the (population- or group-level) parameter  
#' @param group Grouping factor of group-level parameters.
#' @param nlpar Name of a non-linear / auxiliary parameter. 
#'   Only used in non-linear / distributional models.
#' @param resp Name of the response variable / category.
#'   Only used in multivariate and categorical models.
#'   Is internally handled as an alias of \code{nlpar}.
#' @param lb Lower bound for parameter restriction. Currently only allowed
#'   for classes \code{"b"}, \code{"ar"}, \code{"ma"}, and \code{"arr"}.
#'   Defaults to \code{NULL}, that is no restriction.
#' @param ub Upper bound for parameter restriction. Currently only allowed
#'   for classes \code{"b"}, \code{"ar"}, \code{"ma"}, and \code{"arr"}.
#'   Defaults to \code{NULL}, that is no restriction.
#' @param ... Arguments passed to \code{set_prior}.
#' 
#' @return An object of class \code{brmsprior} to be used in the \code{prior}
#'   argument of \code{\link[brms:brm]{brm}}.
#' 
#' @details 
#'   \code{set_prior} is used to define prior distributions for parameters 
#'   in \pkg{brms} models. The functions \code{prior} and \code{prior_string} 
#'   are both aliases of \code{set_prior}, the former allowing to specify 
#'   arguments without quotes \code{""} using non-standard evaluation.
#'   
#'   Below, we explain its usage and list some common 
#'   prior distributions for parameters. 
#'   A complete overview on possible prior distributions is given 
#'   in the Stan Reference Manual available at \url{http://mc-stan.org/}.
#'   
#'   To combine multiple priors, use \code{c(...)}, 
#'   e.g., \code{c(set_prior(...), set_prior(...))}.
#'   \pkg{brms} does not check if the priors are written in correct \pkg{Stan} language. 
#'   Instead, \pkg{Stan} will check their syntactical correctness when the model 
#'   is parsed to \code{C++} and returns an error if they are not. 
#'   This, however, does not imply that priors are always meaningful if they are 
#'   accepted by \pkg{Stan}. Although \pkg{brms} trys to find common problems 
#'   (e.g., setting bounded priors on unbounded parameters), there is no guarantee 
#'   that the defined priors are reasonable for the model.
#'   Currently, there are seven types of parameters in \pkg{brms} models, 
#'   for which the user can specify prior distributions. \cr
#'   
#'   1. Population-level ('fixed') effects
#'   
#'   Every Population-level effect has its own regression parameter 
#    These parameters are internally named as \code{b_<fixed>}, where \code{<fixed>} 
#'   represents the name of the corresponding population-level effect. 
#'   Suppose, for instance, that \code{y} is predicted by \code{x1} and \code{x2} 
#'   (i.e. \code{y ~ x1+x2} in formula syntax). 
#'   Then, \code{x1} and \code{x2} have regression parameters 
#'   \code{b_x1} and \code{b_x2} respectively. 
#'   The default prior for population-level effects (including monotonic and 
#'   category specific effects) is an improper flat prior over the reals. 
#'   Other common options are normal priors or student-t priors. 
#'   If we want to have a normal prior with mean 0 and 
#'   standard deviation 5 for \code{x1}, and a unit student-t prior with 10 
#'   degrees of freedom for \code{x2}, we can specify this via
#'   \code{set_prior("normal(0,5)", class = "b", coef = "x1")} and \cr
#'   \code{set_prior("student_t(10,0,1)", class = "b", coef = "x2")}.
#'   To put the same prior on all fixed effects at once, 
#'   we may write as a shortcut \code{set_prior("<prior>", class = "b")}. 
#'   This also leads to faster sampling, because priors can be vectorized in this case. 
#'   Both ways of defining priors can be combined using for instance 
#'   \code{set_prior("normal(0,2)", class = "b")} and \cr
#'   \code{set_prior("normal(0,10)", class = "b", coef = "x1")}
#'   at the same time. This will set a \code{normal(0,10)} prior on 
#'   the fixed effect of \code{x1} and a \code{normal(0,2)} prior 
#'   on all other fixed effects. However, this will break vectorization and
#'   may slow down the sampling procedure a bit.
#'   
#'   In case of the default intercept parameterization 
#'   (discussed in the 'Details' section of \code{\link[brms:brm]{brm}}),
#'   the fixed effects intercept has its own parameter class 
#'   named \code{"Intercept"} and priors can thus be 
#'   specified via \code{set_prior("<prior>", class = "Intercept")}.
#'   Setting a prior on the intercept will not break vectorization
#'   of the other population-level effects.
#'   
#'   A special shrinkage prior to be applied on population-level effects 
#'   is the horseshoe prior.
#'   It is symmetric around zero with fat tails and an infinitely large spike
#'   at zero. This makes it ideal for sparse models that have 
#'   many regression coefficients, although only a minority of them is non-zero. 
#'   For more details see Carvalho et al. (2009).
#'   The horseshoe prior can be applied on all population-level effects at once 
#'   (excluding the intercept) by using \code{set_prior("horseshoe(1)")}.
#'   The \code{1} implies that the student-t prior of the local shrinkage 
#'   parameters has 1 degrees of freedom. This may, however, lead to an 
#'   increased number of divergent transition in \pkg{Stan}.
#'   Accordingly, increasing the degrees of freedom to slightly higher values 
#'   (e.g., \code{3}) may often be a better option, although the prior 
#'   no longer resembles a horseshoe in this case. 
#'   Further, the scale of the global shrinkage parameter plays an important role
#'   in amount of shrinkage applied. It defaults to \code{1},
#'   but this may result in too few shrinkage (Piironen & Vehtari, 2016).
#'   It is thus possible to change the scale using argument \code{scale_global}
#'   of the horseshoe prior, for instance \code{horseshoe(1, scale_global = 0.5)}.
#'   For recommendations how to properly set the global scale see 
#'   Piironen and Vehtari (2016).
#'   Generally, models with horseshoe priors a more likely than other models
#'   to have divergent transitions so that increasing \code{adapt_delta} 
#'   from \code{0.8} to values closer to \code{1} will often be necessary.
#'   See the documentation of \code{\link[brms:brm]{brm}} for instructions
#'   on how to increase \code{adapt_delta}. \cr
#'   
#'   In non-linear models, population-level effects are defined separately 
#'   for each non-linear parameter. Accordingly, it is necessary to specify
#'   the non-linear parameter in \code{set_prior} so that priors
#'   we can be assigned correctly. 
#'   If, for instance, \code{alpha} is the parameter and \code{x} the predictor
#'   for which we want to define the prior, we can write
#'   \code{set_prior("<prior>", coef = "x", nlpar = "alpha")}. 
#'   As a shortcut we can use \code{set_prior("<prior>", nlpar = "alpha")}
#'   to set the same prior on all population-level effects of \code{alpha} at once.
#'   
#'   If desired, population-level effects can be restricted to fall only 
#'   within a certain interval using the \code{lb} and \code{ub} arguments
#'   of \code{set_prior}. This is often required when defining priors
#'   that are not defined everywhere on the real line, such as uniform
#'   or gamma priors. When defining a \code{uniform(2,4)} prior, 
#'   you should write \code{set_prior("uniform(2,4)", lb = 2, ub = 4)}. 
#'   When using a prior that is defined on the postive reals only 
#'   (such as a gamma prior) set \code{lb = 0}. 
#'   In most situations, it is not useful to restrict population-level
#'   parameters through bounded priors 
#'   (non-linear models are an important exception), 
#'   but if you really want to this is the way to go.
#'   
#'   2. Standard deviations of group-level ('random') effects
#'   
#'   Each group-level effect of each grouping factor has a standard deviation named
#'   \code{sd_<group>_<random>}. Consider, for instance, the formula 
#'   \code{y ~ x1+x2+(1+x1|g)}.
#'   We see that the intercept as well as \code{x1} are group-level effects
#'   nested in the grouping factor \code{g}. 
#'   The corresponding standard deviation parameters are named as 
#'   \code{sd_g_Intercept} and \code{sd_g_x1} respectively. 
#'   These parameters are restriced to be non-negative and, by default, 
#'   have a half student-t prior with 3 degrees of freedom and a 
#'   scale parameter that depends on the standard deviation of the response 
#'   after applying the link function. Minimally, the scale parameter is 10. 
#'   This prior is used (a) to be only very weakly informative in order to influence
#'   results as few as possible, while (b) providing at least some regularization
#'   to considerably improve convergence and sampling efficiency.
#'   To define a prior distribution only for standard deviations 
#'   of a specific grouping factor,
#'   use \cr \code{set_prior("<prior>", class = "sd", group = "<group>")}. 
#'   To define a prior distribution only for a specific standard deviation 
#'   of a specific grouping factor, you may write \cr
#'   \code{set_prior("<prior>", class = "sd", group = "<group>", coef = "<coef>")}. 
#'   Recommendations on useful prior distributions for 
#'   standard deviations are given in Gelman (2006), but note that he
#'   is no longer recommending uniform priors, anymore. \cr
#'   
#'   When defining priors on group-level effects parameters in non-linear models, 
#'   please make sure to specify the corresponding non-linear parameter 
#'   through the \code{nlpar} argument in the same way as 
#'   for population-level effects.
#'   
#'   3. Correlations of group-level ('random') effects 
#'   
#'   If there is more than one group-level effect per grouping factor, 
#'   the correlations between those effects have to be estimated. 
#'   The prior \code{"lkj_corr_cholesky(eta)"} or in short 
#'   \code{"lkj(eta)"} with \code{eta > 0} 
#'   is essentially the only prior for (choelsky factors) of correlation matrices. 
#'   If \code{eta = 1} (the default) all correlations matrices 
#'   are equally likely a priori. If \code{eta > 1}, extreme correlations 
#'   become less likely, whereas \code{0 < eta < 1} results in 
#'   higher probabilities for extreme correlations. 
#'   Correlation matrix parameters in \code{brms} models are named as 
#'   \code{cor_(group)}, (e.g., \code{cor_g} if \code{g} is the grouping factor).
#'   To set the same prior on every correlation matrix, 
#'   use for instance \code{set_prior("lkj(2)", class = "cor")}.
#'   
#'   4. Standard deviations of smoothing terms
#'   
#'   GAMMs are implemented in \pkg{brms} using the 'random effects' 
#'   formulation of smoothing terms (for details see 
#'   \code{\link[mgcv:gamm]{gamm}}). Thus, each smoothing term
#'   has its corresponding standard deviation modeling
#'   the variability within this term. In \pkg{brms}, this 
#'   parameter class is called \code{sds} and priors can
#'   be specified via \code{set_prior("<prior>", class = "sds", 
#'   coef = "<term label>")}. The default prior is the same as
#'   for standard deviations of group-level effects.
#'   
#'   5. Autocorrelation parameters
#'   
#'   The autocorrelation parameters currently implemented are named 
#'   \code{ar} (autoregression), \code{ma} (moving average),
#'   and \code{arr} (autoregression of the response).
#'   
#'   Priors can be defined by \code{set_prior("<prior>", class = "ar")} 
#'   for \code{ar} and similar for \code{ma} and \code{arr} effects.
#'   By default, \code{ar} and \code{ma} are bounded between \code{-1} 
#'   and \code{1} and \code{arr} is unbounded (you may change this 
#'   by using the arguments \code{lb} and \code{ub}). The default
#'   prior is flat over the definition area.
#'   
#'   6. Distance parameters of monotonic effects
#'   
#'   As explained in the details section of \code{\link[brms:brm]{brm}},
#'   monotonic effects make use of a special parameter vector to
#'   estimate the 'normalized distances' between consecutive predictor 
#'   categories. This is realized in \pkg{Stan} using the \code{simplex}
#'   parameter type and thus this class is also named \code{"simplex"} in
#'   \pkg{brms}. The only valid prior for simplex parameters is the
#'   dirichlet prior, which accepts a vector of length \code{K - 1}
#'   (K = number of predictor categories) as input defining the
#'   'concentration' of the distribution. Explaining the dirichlet prior 
#'   is beyond the scope of this documentation, but we want to describe
#'   how to define this prior syntactically correct.
#'   If a predictor \code{x} with \code{K} categories is modeled as monotonic, 
#'   we can define a prior on its corresponding simplex via \cr
#'   \code{set_prior("dirichlet(<vector>)", class = "simplex", coef = "x")}.
#'   For \code{<vector>}, we can put in any \code{R} expression
#'   defining a vector of length \code{K - 1}. The default is a uniform 
#'   prior (i.e. \code{<vector> = rep(1, K-1)}) over all simplexes
#'   of the respective dimension.   
#'   
#'   7. Parameters for specific families 
#'   
#'   Some families need additional parameters to be estimated. 
#'   Families \code{gaussian}, \code{student}, and \code{cauchy} 
#'   need the parameter \code{sigma} 
#'   to account for the residual standard deviation.
#'   By default, \code{sigma} has a half student-t prior that scales 
#'   in the same way as the random effects standard deviations.
#'   Furthermore, family \code{student} needs the parameter 
#'   \code{nu} representing the degrees of freedom of students t distribution. 
#'   By default, \code{nu} has prior \code{"gamma(2,0.1)"}
#'   and a fixed lower bound of \code{0}.
#'   Families \code{gamma}, \code{weibull}, \code{inverse.gaussian}, and
#'   \code{negbinomial} need a \code{shape} parameter that has a 
#'   \code{"gamma(0.01,0.01)"} prior by default. 
#'   For families \code{cumulative}, \code{cratio}, \code{sratio}, 
#'   and \code{acat}, and only if \code{threshold = "equidistant"}, 
#'   the parameter \code{delta} is used to model the distance between 
#'   two adjacent thresholds. 
#'   By default, \code{delta} has an improper flat prior over the reals.
#'   The \code{von_mises} family needs the parameter \code{kappa}, representing
#'   the concentration parameter. By default, \code{kappa} has prior 
#'   \code{"gamma(2, 0.01)"}. \cr
#'   Every family specific parameter has its own prior class, so that
#'   \code{set_prior("<prior>", class = "<parameter>")} is the right way to go.
#'   All of these priors are chosen to be weakly informative,
#'   having only minimal influence on the estimations,
#'   while improving convergence and sampling efficiency.
#' 
#'   Often, it may not be immediately clear, 
#'   which parameters are present in the model.
#'   To get a full list of parameters and parameter classes for which 
#'   priors can be specified (depending on the model) 
#'   use function \code{\link[brms:get_prior]{get_prior}}.
#'
#' @seealso \code{\link[brms:get_prior]{get_prior}}
#' 
#' @references 
#' Carvalho, C. M., Polson, N. G., & Scott, J. G. (2009). 
#'   Handling sparsity via the horseshoe. 
#'   In International Conference on Artificial Intelligence and Statistics (pp. 73-80).
#' 
#' Gelman A. (2006). Prior distributions for variance parameters in hierarchical models.
#'    Bayesian analysis, 1(3), 515 -- 534.
#'    
#' Piironen J. & Vehtari A. (2016). On the Hyperprior Choice for the Global 
#'    Shrinkage Parameter in the Horseshoe Prior. 
#'    \url{https://arxiv.org/pdf/1610.05559v1.pdf}
#'   
#' 
#' @examples
#' ## check which parameters can have priors
#' get_prior(rating ~ treat + period + carry + (1|subject),
#'           data = inhaler, family = sratio(), 
#'           threshold = "equidistant")
#'          
#' ## define some priors          
#' prior <- c(set_prior("normal(0,10)", class = "b"),
#'            set_prior("normal(1,2)", class = "b", coef = "treat"),
#'            set_prior("cauchy(0,2)", class = "sd", 
#'                      group = "subject", coef = "Intercept"),
#'            set_prior("uniform(-5,5)", class = "delta"))
#'               
#' ## verify that the priors indeed found their way into Stan's model code
#' make_stancode(rating ~ period + carry + cse(treat) + (1|subject),
#'               data = inhaler, family = sratio(), 
#'               threshold = "equidistant",
#'               prior = prior)
#'               
#' ## use horseshoe priors to model sparsity in population-level effects parameters
#' make_stancode(count ~ log_Age_c + log_Base4_c * Trt_c,
#'               data = epilepsy, family = poisson(),
#'               prior = set_prior("horseshoe(3)"))
#'               
#' ## use alias functions
#' (prior1 <- prior_string("cauchy(0, 1)", class = "sd"))
#' (prior2 <- prior(cauchy(0, 1), class = sd))
#' identical(prior1, prior2)
#'
#' @export
set_prior <- function(prior, class = "b", coef = "", group = "",
                      nlpar = "", resp = NULL, lb = NULL, ub = NULL) {
  prior <- as.character(prior)
  class <- as.character(class)
  group <- as.character(group)
  coef <- as.character(coef)
  nlpar <- as.character(use_alias(nlpar, resp, warn = FALSE))
  lb <- as.numeric(lb)
  ub <- as.numeric(ub)
  if (length(prior) != 1 || length(class) != 1 || length(coef) != 1 || 
      length(group) != 1 || length(nlpar) != 1 || length(lb) > 1 || 
      length(ub) > 1)
    stop("All arguments of set_prior must be of length 1.", call. = FALSE)
  valid_classes <- c("Intercept", "b", "sd", "sds", "simplex", "cor", "L", 
                     "ar", "ma", "arr", "sigma", "sigmaLL", "rescor", 
                     "Lrescor", "nu", "shape", "delta", "phi", "kappa")
  if (!class %in% valid_classes) {
    stop(paste(class, "is not a valid parameter class"), call. = FALSE)
  }
  if (nchar(group) && !class %in% c("sd", "cor", "L")) {
    stop(paste("argument 'group' not meaningful for class", class), 
         call. = FALSE)
  }
  coef_classes <- c("Intercept", "b", "sd", "sds", "sigma", "simplex")
  if (nchar(coef) && !class %in% coef_classes) {
    stop(paste("argument 'coef' not meaningful for class", class),
         call. = FALSE)
  }
  if (nchar(nlpar) && !class %in% valid_classes[1:5]) {
    stop(paste("argument 'nlpar' not meaningful for class", class),
         call. = FALSE)
  }
  is_arma <- class %in% c("ar", "ma")
  if (length(lb) || length(ub) || is_arma) {
    if (!(class %in% c("b", "arr") || is_arma))
      stop("Currently boundaries are only allowed for ", 
           "population-level and ARMA effects.", call. = FALSE)
    if (nchar(coef)) {
      stop("'coef' may not be specified when using boundaries")
    }
    if (is_arma) {
      lb <- ifelse(length(lb), lb, -1)
      ub <- ifelse(length(ub), ub, 1) 
      if (is.na(lb) || is.na(ub) || abs(lb) > 1 || abs(ub) > 1) {
        warning("Setting boundaries of ARMA parameters outside of ", 
                "[-1,1] may not be appropriate.", call. = FALSE)
      }
    }
    # don't put spaces in boundary declarations
    lb <- if (length(lb) && !is.na(lb)) paste0("lower=", lb)
    ub <- if (length(ub) && !is.na(ub)) paste0("upper=", ub)
    if (!is.null(lb) || !is.null(ub)) {
      bound <- paste0("<", paste(c(lb, ub), collapse = ","), ">")
    } else {
      bound <- ""
    }
  } else {
    bound <- ""
  }
  if (grepl("^increment_log_prob\\(", prior)) {
    # increment_log_prob can be used to directly add a term 
    # to the log posterior
    class <- coef <- group <- nlpar <- ""
  }
  do.call(brmsprior, nlist(prior, class, coef, group, nlpar, bound))
}

#' @describeIn set_prior Alias of \code{set_prior}.
#' @export
prior_string <- function(prior, ...) {
  set_prior(prior, ...)
}

#' @describeIn set_prior Alias of \code{set_prior} allowing to specify 
#'   arguments without quotes \code{""} using non-standard evaluation.
#' @export
prior <- function(prior, ...) {
  call <- as.list(match.call()[-1])
  call <- lapply(call, deparse_no_string)
  do.call(set_prior, call)
}

#' Overview on Priors for \pkg{brms} Models
#' 
#' Get information on all parameters (and parameter classes) for which priors 
#' may be specified including default priors.
#' 
#' @inheritParams brm
#' @param internal A flag indicating if the names of additional internal parameters should be displayed. 
#'   Setting priors on these parameters is not recommended
#' 
#' @return A data.frame with columns \code{prior}, \code{class}, \code{coef}, and \code{group}
#'   and several rows, each providing information on a parameter (or parameter class) on which
#'   priors can be specified. The prior column is empty except for internal default priors.
#'   
#' @seealso \code{\link[brms:set_prior]{set_prior}}
#' 
#' @examples 
#' ## get all parameters and parameters classes to define priors on
#' (prior <- get_prior(count ~ log_Age_c + log_Base4_c * Trt_c
#'                     + (1|patient) + (1|visit),
#'                     data = epilepsy, family = poisson()))   
#'          
#' ## define a prior on all population-level effects a once
#' prior$prior[1] <- "normal(0,10)"
#' 
#' ## define a specific prior on the population-level effect of Trt_c
#' prior$prior[5] <- "student_t(10, 0, 5)"       
#' 
#' ## verify that the priors indeed found their way into Stan's model code
#' make_stancode(count ~ log_Age_c + log_Base4_c * Trt_c 
#'               + (1|patient) + (1|visit),
#'               data = epilepsy, family = poisson(), 
#'               prior = prior)
#' 
#' @export
get_prior <- function(formula, data = NULL, family = gaussian(),
                      autocor = NULL, nonlinear = NULL, partial = NULL, 
                      threshold = c("flexible", "equidistant"), 
                      internal = FALSE) {
  # note that default priors are stored in this function
  family <- check_family(family) 
  link <- family$link
  formula <- update_formula(formula, data = data, family = family, 
                            partial = partial, nonlinear = nonlinear)
  threshold <- match.arg(threshold)
  autocor <- check_autocor(autocor)
  ee <- extract_effects(formula, family = family)
  data <- update_data(data, family = family, effects = ee)
  
  # ensure that RE and residual SDs only have a weakly informative prior by default
  Y <- unname(model.response(data))
  prior_scale <- 10
  if (is.lognormal(family)) link <- "log"
  if (link %in% c("identity", "log", "inverse", "sqrt", "1/mu^2")) {
    if (link %in% c("log", "inverse", "1/mu^2")) {
      Y <- ifelse(Y == 0, Y + 0.1, Y)  # avoid Inf in link(Y)
    }
    suggested_scale <- SW(round(link(sd(Y), link = link)))
    if (!is.nan(suggested_scale)) {
      prior_scale <- max(prior_scale, suggested_scale, na.rm = TRUE)
    } 
  }
  def_scale_prior <- paste0("student_t(3, 0, ", prior_scale, ")")
  
  # initialize output
  prior <- empty_brmsprior()
  # priors for primary regression effects
  if (length(ee$nonlinear)) {
    nlpars <- names(ee$nonlinear)
    for (i in seq_along(nlpars)) {
      prior_eff <- get_prior_effects(ee$nonlinear[[i]], data = data, 
                                     autocor = autocor, nlpar = nlpars[i],
                                     spec_intercept = FALSE,
                                     def_scale_prior = def_scale_prior,
                                     internal = internal)
      prior <- rbind(prior, prior_eff)
    }
  } else {
    if (length(ee$response) > 1L) {
      # priors for effects in multivariate models
      for (r in c("", ee$response)) {
        # r = "" adds global priors affecting parameters of all responses
        prior_eff <- get_prior_effects(ee, data = data, autocor = autocor,
                                       def_scale_prior = def_scale_prior,
                                       internal = internal, nlpar = r)
        prior <- rbind(prior, prior_eff)
      }
    } else {
      # priors for effects in univariate models
      prior_eff <- get_prior_effects(ee, data = data, autocor = autocor,
                                     def_scale_prior = def_scale_prior,
                                     internal = internal)
      prior <- rbind(prior, prior_eff)
    }
  }
  # priors for auxiliary parameters
  def_auxprior <- c(sigma = def_scale_prior, shape = "gamma(0.01, 0.01)",
                    nu = "gamma(2, 0.1)", phi = "gamma(0.01, 0.01)",
                    kappa = "gamma(2, 0.01)", zi = "beta(1, 1)", 
                    hu = "beta(1, 1)")
  valid_auxpars <- valid_auxpars(family, effects = ee, autocor = autocor)
  for (ap in valid_auxpars) {
    if (!is.null(ee[[ap]])) {
      auxprior <- get_prior_effects(ee[[ap]], data = data,
                                    autocor = autocor, nlpar = ap,
                                    spec_intercept = FALSE,
                                    def_scale_prior = def_scale_prior,
                                    internal = internal)
    } else {
      auxprior <- brmsprior(class = ap, prior = def_auxprior[ap])
    }
    prior <- rbind(prior, auxprior)
  }
  # priors of group-level parameters
  ranef <- tidy_ranef(ee, data)
  prior_ranef <- get_prior_ranef(ranef, def_scale_prior = def_scale_prior,
                                 global_sd = length(ee$response) > 1L,
                                 internal = internal)
  prior <- rbind(prior, prior_ranef)
  
  # prior for the delta parameter for equidistant thresholds
  if (is.ordinal(family) && threshold == "equidistant") {
    prior <- rbind(prior, brmsprior(class = "delta"))
  }
  # priors for auxiliary parameters of multivariate models
  if (is.linear(family) && length(ee$response) > 1L) {
    sigma_coef <- c("", ee$response)
    sigma_prior <- c(def_scale_prior, rep("", length(ee$response)))
    sigma_prior <- brmsprior(class = "sigma", coef = sigma_coef,
                             prior = sigma_prior)
    prior <- rbind(prior, sigma_prior)
    if (internal) {
      prior <- rbind(prior, brmsprior(class = "Lrescor", 
                                      prior = "lkj_corr_cholesky(1)"))
    } else {
      prior <- rbind(prior, brmsprior(class = "rescor", prior = "lkj(1)"))
    }
  }
  # priors for autocor parameters
  cbound <- "<lower=-1,upper=1>"
  if (get_ar(autocor)) {
    prior <- rbind(prior, brmsprior(class = "ar", bound = cbound))
  }
  if (get_ma(autocor)) {
    prior <- rbind(prior, brmsprior(class = "ma", bound = cbound))
  }
  if (get_arr(autocor)) {
    prior <- rbind(prior, brmsprior(class = "arr"))
  }
  if (is(autocor, "cor_bsts")) {
    prior <- rbind(prior, brmsprior(class = "sigmaLL", 
                                    prior = def_scale_prior))
  }
  # do not remove unique(.)
  prior <- unique(prior[with(prior, order(nlpar, class, group, coef)), ])
  rownames(prior) <- seq_len(nrow(prior))
  structure(prior, class = c("brmsprior", "data.frame"))
}

get_prior_effects <- function(effects, data, autocor = cor_arma(), 
                              nlpar = "", spec_intercept = TRUE,
                              def_scale_prior = "", internal = FALSE) {
  # wrapper function to get priors for various kinds of effects
  # don't use the family argument here to avoid
  # removal of the intercept for ordinal models
  # Args:
  #   spec_intercept: special parameter class for the FE Intercept? 
  fixef <- colnames(data_fixef(effects, data, autocor = autocor)$X)
  spec_intercept <- has_intercept(effects$fixed) && spec_intercept
  prior_fixef <- get_prior_fixef(fixef, spec_intercept = spec_intercept,
                                 nlpar = nlpar, internal = internal)
  monef <- colnames(get_model_matrix(effects$mono, data))
  prior_monef <- get_prior_monef(monef, fixef = fixef, nlpar = nlpar)
  # group-level priors are prepared separately
  splines <- get_spline_labels(effects)
  prior_splines <- get_prior_splines(splines, def_scale_prior, nlpar = nlpar)
  csef <- colnames(get_model_matrix(effects$cse, data = data))
  prior_csef <- get_prior_csef(csef, fixef = fixef)
  rbind(prior_fixef, prior_monef, prior_splines, prior_csef)
}

get_prior_fixef <- function(fixef, spec_intercept = TRUE, 
                            nlpar = "", internal = FALSE) {
  # priors for fixed effects parameters
  # Args:
  #   fixef: names of the fixed effects
  #   spec_intercept: special parameter class for the Intercept? 
  #   internal: see get_prior
  # Returns:
  #   an object of class brmsprior
  prior <- empty_brmsprior()
  if (length(fixef)) {
    prior <- rbind(prior, brmsprior(class = "b", coef = c("", fixef),
                                    nlpar = nlpar)) 
  }
  if (spec_intercept) {
    prior <- rbind(prior, brmsprior(class = "Intercept", coef = "",
                                    nlpar = nlpar))
    if (internal) {
      prior <- rbind(prior, brmsprior(class = "temp_Intercept",
                                      coef = "", nlpar = nlpar))
    }
  }
  prior
}

get_prior_monef <- function(monef, fixef = NULL, nlpar = "") {
  # priors for monotonic effects parameters
  # Args:
  #   monef: names of the monotonic effects
  #   fixef: names of the fixed effects
  #   nlpar: optional name of a non-linear parameter
  # Returns:
  #   an object of class brmsprior
  prior <- empty_brmsprior()
  if (length(monef)) {
    invalid <- intersect(fixef, monef)
    if (length(invalid)) {
      stop(paste("Variables cannot be modeled as fixed and", 
                 "monotonic effects at the same time.", 
                 "\nError occured for variables:", 
                 paste(invalid, collapse = ", ")), call. = FALSE)
    }
    prior <- rbind(brmsprior(class = "b", coef = c("", monef), 
                             nlpar = nlpar),
                   brmsprior(class = "simplex", coef = monef, 
                             nlpar = nlpar))
  }
  prior
}

get_prior_csef <- function(csef, fixef = NULL) {
  # priors for category spcific effects parameters
  # Args:
  #   csef: names of the category specific effects
  #   fixef: names of the fixed effects
  # Returns:
  #   an object of class brmsprior
  prior <- empty_brmsprior()
  if (length(csef)) {
    invalid <- intersect(fixef, csef)
    if (length(invalid)) {
      stop(paste("Variables cannot be modeled as fixed and", 
                 "category specific effects at the same time.", 
                 "\nError occured for variables:", 
                 paste(invalid, collapse = ", ")), call. = FALSE)
    }
    prior <- brmsprior(class = "b", coef = c("", csef))
  }
  prior
}

get_prior_ranef <- function(ranef, def_scale_prior, 
                            global_sd = FALSE, internal = FALSE) {
  # priors for random effects parameters
  # Args:
  #   ranef: a list returned by tidy_ranef
  #   def_scale_prior: a character string defining the default
  #                    prior for random effects SDs
  #   global_sd: allow to set a global SD prior
  #              affecting all non-linear parameters?
  #   internal: see get_prior
  # Returns:
  #   an object of class brmsprior
  prior <- empty_brmsprior()
  if (nrow(ranef)) {
    # global sd class
    nlpars <- unique(ranef$nlpar)
    if (global_sd) {
      global_sd_prior <- rep("", length(setdiff(nlpars, "")))
      global_sd_prior <- c(def_scale_prior, global_sd_prior)
      global_sd_prior <- brmsprior(class = "sd", 
                                   prior = global_sd_prior,
                                   nlpar = union("", nlpars))
    } else {
      global_sd_prior <- brmsprior(class = "sd", 
                                   prior = def_scale_prior,
                                   nlpar = nlpars)
    }
    prior <- rbind(prior, global_sd_prior)
    for (id in unique(ranef$id)) {
      r <- ranef[ranef$id == id, ]
      group <- r$group[1]
      # include group-level standard deviations
      prior <- rbind(prior,
        brmsprior(class = "sd", group = group, 
                  nlpar = unique(r$nlpar)),
        brmsprior(class = "sd", coef = r$coef, 
                  group = group, nlpar = r$nlpar))
      # detect duplicated group-level effects
      J <- with(prior, class == "sd" & nchar(coef))
      dupli <- duplicated(prior[J, ])
      if (any(dupli)) {
        stop("Duplicated group-level effects detected for group ", 
             group, call. = FALSE)
      }
      # include correlation parameters
      if (isTRUE(r$cor[1]) && nrow(r) > 1L) {
        if (internal) {
          prior <- rbind(prior, 
            brmsprior(class = "L", group = c("", group),
                      prior = c("lkj_corr_cholesky(1)", "")))
        } else {
          prior <- rbind(prior, 
            brmsprior(class = "cor", group = c("", group),
                      prior = c("lkj(1)", "")))
        }
      }
    }
  } 
  prior
}

get_prior_splines <- function(splines, def_scale_prior, nlpar = "") {
  # priors for GAMM models
  # Args:
  #   splines: names of the spline terms
  #   def_scale_prior: a character string defining the default
  #                    prior for spline SDs
  #   nlpar: optional name of a non-linear parameter
  if (length(splines)) {
    prior_strings <- c(def_scale_prior, rep("", length(splines)))
    prior <- brmsprior(class = "sds", coef = c("", splines), 
                       prior = prior_strings, nlpar = nlpar)
  } else {
    prior <- empty_brmsprior()
  }
  prior
}

check_prior <- function(prior, formula, data = NULL, family = gaussian(), 
                        sample_prior = FALSE, autocor = NULL, 
                        threshold = "flexible", check_rows = NULL, 
                        warn = FALSE) {
  # check prior input and amend it if needed
  # Args:
  #   same as the respective parameters in brm
  #   (nonlinear is expected to be an attribute of formula)
  #   check_rows: if not NULL, check only the rows given in check_rows
  #   warn: passed to check_prior_content
  # Returns:
  #   a data.frame of prior specifications to be used in stan_prior (see stan.R)
  if (isTRUE(attr(prior, "checked"))) {
    # prior has already been checked; no need to do it twice
    return(prior)
  }
  stopifnot(is(formula, "brmsformula"))
  ee <- extract_effects(formula, family = family)  
  all_priors <- get_prior(formula = formula, data = data, 
                          family = family, autocor = autocor, 
                          threshold = threshold, internal = TRUE)
  if (is.null(prior)) {
    prior <- all_priors  
  }
  # exclude priors using increment_log_prob to readd them at the end
  has_incr_lp <- grepl("^increment_log_prob\\(", prior$prior)
  prior_incr_lp <- prior[has_incr_lp, ]
  prior <- prior[!has_incr_lp, ]
  # check for duplicated priors
  prior$class <- rename(prior$class, symbols = c("^cor$", "^rescor$"), 
                        subs = c("L", "Lrescor"), fixed = FALSE)
  duplicated_input <- duplicated(prior[, 2:5])
  if (any(duplicated_input)) {
    stop("Duplicated prior specifications are not allowed.", call. = FALSE)
  }
  # handle special priors that are not explictly coded as functions in Stan
  has_specef <- is.formula(ee[["mono"]]) || is.formula(ee[["cse"]])
  prior <- handle_special_priors(prior, has_specef = has_specef)  
  # check if parameters in prior are valid
  if (nrow(prior)) {
    valid <- which(duplicated(rbind(all_priors[, 2:5], prior[, 2:5])))
    invalid <- which(!1:nrow(prior) %in% (valid - nrow(all_priors)))
    if (length(invalid)) {
      msg_priors <- .print_prior(prior[invalid, ])
      message("The following priors don't correspond to any ", 
              "model parameter \nand will thus not affect the results: \n",
              collapse(.print_prior(prior[invalid, ]), "\n"))
      prior <- prior[-invalid, ]
    }
  }
  check_prior_content(prior, family = family, warn = warn)
  # merge prior with all_priors
  prior <- rbind(prior, all_priors)
  prior <- prior[!duplicated(prior[, 2:5]), ]
  rows2remove <- NULL
  # special treatment of fixed effects Intercepts
  int_index <- which(prior$class == "Intercept")
  if (length(int_index)) {
    int_prior <- prior[int_index, ]
    if (length(int_index) > 1L) {
      intercepts <- prior$coef[int_index]
      intercepts <- intercepts[nchar(intercepts) > 0]
    } else intercepts <- "Intercept"
    bint_index <- which(prior$class == "b" & prior$coef %in% intercepts)
    bint_prior <- prior[bint_index, ]
    for (t in which(prior$class %in% "temp_Intercept")) {
      ti <- int_prior$coef == prior$coef[t]
      tb <- bint_prior$coef %in% c(prior$coef[t], "Intercept")
      if (sum(ti) && nchar(int_prior$prior[ti]) > 0) {
        # take 'Intercept' priors first if specified
        prior$prior[t] <- int_prior$prior[ti]
      } else if (sum(tb) && nchar(bint_prior$prior[tb]) > 0) {
        # fall back to 'b' (fixed effects) priors
        prior$prior[t] <- bint_prior$prior[tb]
      }
    }
    rows2remove <- c(rows2remove, int_index, bint_index)
  }
  # prepare priors of monotonic effects
  mono_forms <- get_effect(ee, "mono")
  for (k in seq_along(mono_forms)) {
    monef <- colnames(get_model_matrix(mono_forms[[k]], data = data))
    for (i in seq_along(monef)) {
      take <- with(prior, class == "simplex" & coef == monef[i] &
                          nlpar == names(mono_forms)[k])
      simplex_prior <- paste0(".", prior$prior[take])
      if (nchar(simplex_prior) > 1L) {
        simplex_prior <- paste(eval(parse(text = simplex_prior)),
                               collapse = ",")
        prior$prior[take] <- paste0("dirichlet(c(", simplex_prior, "))")
      }
    }
  }
  # check if priors for non-linear parameters are defined
  if (length(ee$nonlinear)) {
    nlpars <- names(ee$nonlinear)
    for (nlp in nlpars) {
      nlp_prior <- prior$prior[with(prior, nlpar == nlp & class == "b")]
      if (!any(as.logical(nchar(nlp_prior)))) {
        stop(paste0("Priors on fixed effects are required in non-linear ", 
                    "models, but none were found for parameter '", nlp, 
                    "'. \nSee help(set_prior) for more details."), 
             call. = FALSE)
      }
    }
  }
  if (length(rows2remove)) {   
    prior <- prior[-rows2remove, ]
  }
  prior <- prior[with(prior, order(nlpar, class, group, coef)), ]
  prior <- rbind(prior, prior_incr_lp)
  rownames(prior) <- seq_len(nrow(prior))
  attr(prior, "prior_only") <- identical(sample_prior, "only")
  attr(prior, "checked") <- TRUE
  prior
}

check_prior_content <- function(prior, family = gaussian(), warn = TRUE) {
  # try to check if prior distributions are reasonable
  # Args:
  #  prior: A brmsprior object
  #  family: the model family
  #  warn: logical; print boundary warnings?
  if (!is(prior, "brmsprior")) {
    return(invisible(NULL))
  }
  stopifnot(is(family, "family"))
  family <- family$family
  if (nrow(prior)) {
    lb_priors <- c("lognormal", "chi_square", "inv_chi_square",
                   "scaled_inv_chi_square", "exponential", "gamma",
                   "inv_gamma", "weibull", "frechet", "rayleigh",
                   "pareto", "pareto_type_2")
    lb_priors_reg <- paste0("^(", paste0(lb_priors, collapse = "|"), ")")
    ulb_priors <- c("beta", "uniform", "von_mises")
    ulb_priors_reg <- paste0("^(", paste0(ulb_priors, collapse = "|"), ")")
    nb_pars <- c("b", "Intercept", if (!family %in% "cumulative") "delta")
    lb_pars <- c("sd", "sigma", "nu", "shape", "phi", "kappa",
                 if (family %in% "cumulative") "delta")
    cor_pars <- c("cor", "L", "rescor", "Lrescor")
    autocor_pars <- c("ar", "ma")
    lb_warning <- ub_warning <- ""
    autocor_warning <- FALSE
    for (i in seq_len(nrow(prior))) {
      msg_prior <- .print_prior(prior[i, , drop = FALSE])
      has_lb_prior <- grepl(lb_priors_reg, prior$prior[i])
      has_ulb_prior <- grepl(ulb_priors_reg, prior$prior[i])
      # priors with nchar(coef) inherit their boundaries 
      j <- with(prior, which(class == class[i] & group == group[i] & 
                               nlpar == nlpar[i] & !nchar(coef)))
      bound <- if (length(j)) prior$bound[j] else ""
      has_lb <- grepl("lower", bound)
      has_ub <- grepl("upper", bound)
      if (prior$class[i] %in% nb_pars) {
        if ((has_lb_prior || has_ulb_prior) && !has_lb) {
          lb_warning <- paste0(lb_warning, msg_prior, "\n")
        } 
        if (has_ulb_prior && !has_ub) {
          ub_warning <- paste0(ub_warning, msg_prior, "\n")
        }
      } else if (prior$class[i] %in% lb_pars) {
        if (has_ulb_prior && !has_ub) {
          ub_warning <- paste0(ub_warning, msg_prior, "\n")
        }
      } else if (prior$class[i] %in% cor_pars) {
        if (nchar(prior$prior[i]) && !grepl("^lkj", prior$prior[i])) {
          stop(paste("Currently 'lkj' is the only valid prior",
                     "for group-level correlations. See help(set_prior)",
                     "for more details."), call. = FALSE)
        }
      } else if (prior$class[i] %in% autocor_pars) {
        if (prior$bound[i] != "<lower=-1,upper=1>") {
          autocor_warning <- TRUE
        } 
      } else if (prior$class[i] == "simplex") {
        if (nchar(prior$prior[i]) && !grepl("^dirichlet\\(", prior$prior[i])) {
          stop(paste("Currently 'dirichlet' is the only valid prior",
                     "for simplex parameters. See help(set_prior)",
                     "for more details."), call. = FALSE)
        }
      }
    }  # end for  
    if (nchar(lb_warning) && warn) {
      warning2("It appears that you have specified a lower bounded ", 
               "prior on a parameter that has no natural lower bound.",
               "\nIf this is really what you want, please specify ",
               "argument 'lb' of 'set_prior' appropriately.",
               "\nWarning occurred for prior \n", lb_warning)
    }
    if (nchar(ub_warning) && warn) {
      warning2("It appears that you have specified an upper bounded ", 
               "prior on a parameter that has no natural upper bound.",
               "\nIf this is really what you want, please specify ",
               "argument 'ub' of 'set_prior' appropriately.",
               "\nWarning occurred for prior \n", ub_warning)
    }
    if (autocor_warning && warn) {
      warning2("Changing the boundaries of autocorrelation ", 
               "parameters is not recommended.")
    }
  }
  invisible(NULL)
}

handle_special_priors <- function(prior, has_specef = FALSE) {
  # look for special priors such as horseshoe and process them appropriately
  # Args:
  #   prior: an object of class brmsprior
  #   has_specef: are monotonic or category specific effects present?
  # Returns:
  #   a possibly amended prior.frame with additional attributes
  prior_attr <- list()
  b_index <- which(prior$class == "b" & !nchar(prior$coef))
  if (length(b_index) && grepl("^horseshoe\\(", prior$prior[b_index])) {
    # horseshoe prior for fixed effects parameters
    if (any(nchar(prior$nlpar))) {
      stop2("Horseshoe priors are not yet allowed in non-linear models.")
    }
    if (has_specef) {
      stop2("Horseshoe priors are not yet allowed in models with ", 
            "monotonic or category specific effects.")
    }
    b_coef_indices <- which(prior$class == "b" & nchar(prior$coef) &
                            prior$coef != "Intercept")
    if (any(nchar(prior$prior[b_coef_indices]))) {
      stop2("Defining priors for single population-level parameters",
            "is not allowed when using the horseshoe prior",
            "(except for the Intercept).")
    }
    hs <- eval2(prior$prior[b_index])
    prior_attr[c("hs_df", "hs_scale_global")] <- hs[c("df", "scale_global")]
    prior$prior[b_index] <- hs$prior
  }
  # expand lkj correlation prior to full name
  prior$prior <- sub("^(lkj\\(|lkj_corr\\()", "lkj_corr_cholesky(", prior$prior)
  do.call(structure, c(list(prior), prior_attr))
}

get_bound <- function(prior, class = "b", coef = "", 
                      group = "", nlpar = "") {
  # extract the boundaries of a parameter described by class etc.
  # Args:
  #   prior: object of class brmsprior
  #   class, coef, group, nlpar: strings of length 1
  stopifnot(length(class) == 1L)
  if (!length(coef)) coef <- ""
  if (!length(group)) group <- ""
  if (!length(nlpar)) nlpar <- ""
  take <- prior$class == class & prior$coef == coef & 
          prior$group == group & prior$nlpar == nlpar
  if (sum(take) > 1L) {
    stop("extracted more than one boundary at once")
  }
  prior$bound[take]
}

brmsprior <- function(prior = "", class = "", coef = "", group = "", 
                      nlpar = "", bound = "") {
  # helper function to create data.frames containing prior information 
  out <- data.frame(prior = prior, class = class, coef = coef, 
                    group = group, nlpar = nlpar, bound = bound, 
                    stringsAsFactors = FALSE)
  class(out) <- c("brmsprior", "data.frame")
  out
}

empty_brmsprior <- function() {
  # define a brmsprior object with zero rows
  brmsprior(prior = character(0), class = character(0), 
            coef = character(0), group = character(0),
            nlpar = character(0), bound = character(0))
}

#' @export
print.brmsprior <- function(x, show_df, ...) {
  if (missing(show_df)) {
    show_df <- nrow(x) > 1L
  }
  if (show_df) {
    NextMethod()
  } else {
    cat(paste(.print_prior(x), collapse = "\n")) 
  }
  invisible(x)
}

.print_prior <- function(x) {
  # prepare text for print.brmsprior
  group <-  usc(x$group, "prefix")
  coef <- usc(x$coef, "prefix")
  nlpar <- usc(x$nlpar, "prefix")
  bound <- ifelse(nchar(x$bound), paste0(x$bound, " "), "")
  tilde <- ifelse(nchar(x$class) + nchar(group) + nchar(coef), " ~ ", "")
  prior <- ifelse(nchar(x$prior), x$prior, "(no prior)")
  paste0(bound, x$class, nlpar, group, coef, tilde, prior)
}

#' @export
c.brmsprior <- function(x, ...) {
  # combines multiple brmsprior objects into one brmsprior
  if (all(sapply(list(...), is, class2 = "brmsprior"))) {
    out <- do.call(rbind, list(x, ...)) 
  } else {
    out <- c(as.data.frame(x), ...)
  }
  out
}

dirichlet <- function(...) {
  # dirichlet prior of simplex parameters
  out <- as.numeric(c(...))
  if (anyNA(out) || any(out <= 0)) {
    stop2("The dirichlet prior expects positive values.")
  }
  out
}

horseshoe <- function(df = 1, scale_global = 1) {
  # validate input for the horseshoe prior 
  # Args:
  #   df: degrees of freedom of the local parameters
  #   scale_global: scale of the global cauchy prior
  df <- as.numeric(df)
  scale_global <- as.numeric(scale_global)
  if (!isTRUE(df > 0)) {
    stop2("Invalid horseshoe prior: Degrees of freedom of ", 
          "the local priors must be a single positive number.")
  }
  if (!isTRUE(scale_global > 0)) {
    stop2("Invalid horseshoe prior: Scale of the global ", 
          "prior must be a single positive number.")
  }
  nlist(prior = "normal(0, hs_local * hs_global)", df, scale_global)
}
