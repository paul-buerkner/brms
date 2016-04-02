#' Prior Definitions for \pkg{brms} Models
#'
#' Define priors for specific parameters or classes of parameters
#'
#' @param prior A character string defining a distribution in \pkg{Stan} language
#' @param class The parameter class. Defaults to \code{"b"} (fixed effects). 
#'   See 'Details' for other valid parameter classes. 
#' @param coef Name of the (population- or group-level) parameter  
#' @param group Grouping factor of group-level parameters.
#' @param nlpar Name of a non-linear parameter. Only used in non-linear models.
#' @param lb Lower bound for parameter restriction. Currently only allowed
#'   if \code{class %in% c("b", "ar", "ma", "arr")}. 
#'   Defaults to \code{NULL}, that is no restriction.
#' @param ub Upper bound for parameter restriction. Currently only allowed
#'   if  \code{class %in% c("b", "ar", "ma", "arr")}. 
#'   Defaults to \code{NULL}, that is no restriction.
#' 
#' @return An object of class \code{brmsprior} to be used in the \code{prior}
#'   argument of \code{\link[brms:brm]{brm}}.
#' 
#' @details 
#'   \code{set_prior} is used to define prior distributions for parameters 
#'   in \pkg{brms} models. Below, we explain its usage and list some common 
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
#'   Currently, there are six types of parameters in \pkg{brms} models, 
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
#'   The default prior for population-level effects (including monotonous and 
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
#'   many regression coefficients,although only a minority of them is non-zero. 
#'   For more details see Carvalho et al. (2009).
#'   The horseshoe prior can be applied on all population-level effects at once 
#'   (excluding the intercept) by using \code{set_prior("horseshoe(1)")}.
#'   The \code{1} implies that the student-t prior of the local shrinkage 
#'   parameters has 1 degrees of freedom. This may, however, lead to an 
#'   increased number of divergent transition in \pkg{Stan}.
#'   Accordingly, increasing the degrees of freedom to slightly higher values 
#'   (e.g., \code{3}) may often be a better option, although the prior 
#'   no longer resembles a horseshoe in this case. 
#'   Generally, models with horseshoe priors a more likely than other models
#'   to have divergent transitions so that increasing \code{adapt_delta} 
#'   from \code{0.8} to values closer to \code{1} will often be necessary.
#'   See the documentation of \code{\link[brms:brm]{brm}} for instructions
#'   on how to increase \code{adapt_delta}. \cr
#'   
#'   In non-linear models, population-level effects are defined separately 
#'   for each non-linear parameter. Accordingly, it is necessary to specify
#'   the corresponding non-linear parameter in \code{set_prior} so that priors
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
#'   To define a prior distribution only for standard deviations 
#'   of a specific grouping factor,
#'   use \cr \code{set_prior("<prior>", class = "sd", group = "<group>")}. 
#'   To define a prior distribution only for a specific standard deviation 
#'   of a specific grouping factor, you may write \cr
#'   \code{set_prior("<prior>", class = "sd", group = "<group>", coef = "<coef>")}. 
#'   Recommendations on useful prior distributions for 
#'   standard deviations are given in Gelman (2006). \cr
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
#'   4. Autocorrelation parameters
#'   
#'   The autocorrelation parameters currently implemented are named 
#'   \code{ar} (autoregression), \code{ma} (moving average),
#'   and \code{arr} (autoregression of the response).
#'   Priors can be defined by \code{set_prior("<prior>", class = "ar")} 
#'   for \code{ar} and similar for \code{ma} and \code{arr} effects.
#'   By default, \code{ar} and \code{ma} are bounded between \code{-1} 
#'   and \code{1} and \code{arr} is unbounded (you may change this 
#'   by using the arguments \code{lb} and \code{ub}). The default
#'   prior is flat over the definition area.
#'   
#'   5. Distance parameters of monotonous effects
#'   
#'   As explained in the details section of \code{\link[brms:brm]{brm}},
#'   monotonous effects make use of a special parameter vector to
#'   estimate the 'normalized distances' between consecutive predictor 
#'   categories. This is realized in \pkg{Stan} using the \code{simplex}
#'   parameter type and thus this class is also named \code{"simplex"} in
#'   \pkg{brms}. The only valid prior for simplex parameters is the
#'   dirichlet prior, which accepts a vector of length \code{K - 1}
#'   (K = number of predictor categories) as input defining the
#'   'concentration' of the distribution. Explaining the dirichlet prior 
#'   is beyond the scope of this documentation, but we want to describe
#'   how to define this prior syntactically correct.
#'   If a predictor \code{x} with \code{K} categories is modeled as monotonous, 
#'   we can defined a prior on its corresponding simplex via
#'   \code{set_prior("dirichlet(<vector>)", class = "simplex", coef = "x")}.
#'   For \code{<vector>}, we can put in any \code{R} expression
#'   defining a vector of length \code{K - 1}. The default is a uniform 
#'   prior (i.e. \code{<vector> = rep(1, K-1)}) over all simplexes
#'   of respective dimension.   
#'   
#'   6. Parameters for specific families 
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
#'   and a fixed lower bound of \code{1}.
#'   Families \code{gamma}, \code{weibull}, \code{inverse.gaussian}, and
#'   \code{negbinomial} need a \code{shape} parameter that has a 
#'   \code{"gamma(0.01,0.01)"} prior by default. 
#'   For families \code{cumulative}, \code{cratio}, \code{sratio}, 
#'   and \code{acat}, and only if \code{threshold = "equidistant"}, 
#'   the parameter \code{delta} is used to model the distance between 
#'   two adjacent thresholds. 
#'   By default, \code{delta} has an improper flat prior over the reals. \cr
#'   Every family specific parameter has its own prior class, so that \cr
#'   \code{set_prior("<prior>", class = "<parameter>")} it the right way to go.
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
#' Gelman A (2006). Prior distributions for variance parameters in hierarchical models.
#'    Bayesian analysis, 1(3), 515 -- 534.
#'    
#' Carvalho, C. M., Polson, N. G., & Scott, J. G. (2009). 
#'   Handling sparsity via the horseshoe. 
#'   In International Conference on Artificial Intelligence and Statistics (pp. 73-80).
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
#' @export
set_prior <- function(prior, class = "b", coef = "", group = "",
                      nlpar = "", lb = NULL, ub = NULL) {
  prior <- as.character(prior)
  class <- as.character(class)
  group <- as.character(group)
  coef <- as.character(coef)
  nlpar <- as.character(nlpar)
  lb <- as.numeric(lb)
  ub <- as.numeric(ub)
  if (length(prior) != 1 || length(class) != 1 || length(coef) != 1 || 
      length(group) != 1 || length(nlpar) != 1 || length(lb) > 1 || 
      length(ub) > 1)
    stop("All arguments of set_prior must be of length 1.", call. = FALSE)
  valid_classes <- c("Intercept", "b", "sd", "cor", "L", "ar", "ma", "arr", 
                     "simplex", "sigma", "rescor", "Lrescor", "nu", "shape", 
                     "delta", "phi")
  if (!class %in% valid_classes)
    stop(paste(class, "is not a valid parameter class"), call. = FALSE)
  if (nchar(group) && !class %in% c("sd", "cor", "L"))
    stop(paste("argument 'group' not meaningful for class", class), 
         call. = FALSE)
  if (nchar(coef) && !class %in% c("Intercept", "b", "sd", "sigma", "simplex"))
    stop(paste("argument 'coef' not meaningful for class", class),
         call. = FALSE)
  if (nchar(nlpar) && !class %in% valid_classes[1:5])
    stop(paste("argument 'nlpar' not meaningful for class", class),
         call. = FALSE)
  is_arma <- class %in% c("ar", "ma")
  if (length(lb) || length(ub) || is_arma) {
    if (!(class %in% c("b", "arr") || is_arma))
      stop(paste("Currently boundaries are only allowed", 
                 "for population-level and ARMA effects."), call. = FALSE)
    if (coef != "")
      stop("'coef' may not be specified when using boundaries")
    if (is_arma) {
      lb <- ifelse(length(lb), lb, -1)
      ub <- ifelse(length(ub), ub, 1) 
      if (abs(lb) > 1 || abs(ub) > 1) {
        warning(paste("Setting boundaries of ARMA parameters outside of", 
                      "[-1,1] may not be appropriate."), call. = FALSE)
      }
    }
    # don't put spaces in boundary declarations
    lb <- if (length(lb)) paste0("lower=", lb)
    ub <- if (length(ub)) paste0("upper=", ub)
    bound <- paste0("<", paste(c(lb, ub), collapse = ","), ">")
  } else {
    bound <- ""
  }
  if (grepl("^increment_log_prob\\(", prior)) {
    # increment_log_prob can be used to directly add a term 
    # to the log posterior
    class <- coef <- group <- nlpar <- ""
  }
  out <- nlist(prior, class, coef, group, nlpar, bound)
  class(out) <- c("brmsprior", "list")
  out
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
  if (!(is.null(data) || is.list(data)))
    stop("argument 'data' must be a data.frame or list", call. = FALSE)
  family <- check_family(family) 
  link <- family$link
  nonlinear <- nonlinear2list(nonlinear) 
  formula <- update_formula(formula, data = data, family = family, 
                            partial = partial, nonlinear = nonlinear)
  threshold <- match.arg(threshold)
  autocor <- check_autocor(autocor)
  ee <- extract_effects(formula, family = family,
                        nonlinear = nonlinear)
  data <- update_data(data, family = family, effects = ee)
  ranef <- gather_ranef(ee, data = data, forked = is.forked(family))  
  
  # ensure that RE and residual SDs only have a weakly informative prior by default
  Y <- unname(model.response(data))
  prior_scale <- 10
  if (link %in% c("identity", "log", "inverse", "sqrt", "1/mu^2")) {
    if (link %in% c("log", "inverse", "1/mu^2")) {
      Y <- ifelse(Y == 0, Y + 0.1, Y)  # avoid Inf in link(Y)
    }
    suggested_scale <- round(sd(link(Y, link = link)))
    if (!is.nan(suggested_scale)) {
      prior_scale <- max(prior_scale, suggested_scale, na.rm = TRUE)
    } 
  }
  def_scale_prior <- paste0("student_t(3, 0, ", prior_scale, ")")
  
  # initialize output
  prior <- empty_prior_frame()
  if (length(nonlinear)) {
    nlpars <- names(ee$nonlinear)
    for (i in seq_along(nlpars)) {
      fixef <- colnames(get_model_matrix(ee$nonlinear[[i]]$fixed, data))
      prior_fixef <- get_prior_fixef(fixef, nlpar = nlpars[i],
                                     internal = internal)
      monef <- colnames(get_model_matrix(ee$nonlinear[[i]]$mono, data))
      prior_monef <- get_prior_monef(monef, fixef = fixef, nlpar = nlpars[i])
      prior_ranef <- get_prior_ranef(ranef, def_scale_prior, nlpar = nlpars[i], 
                                     internal = internal)
      prior <- rbind(prior, prior_fixef, prior_monef, prior_ranef)
    }
  } else {
    # don't remove the intercept columns here!
    fixef <- colnames(get_model_matrix(rhs(ee$fixed), data = data,
                                       forked = is.forked(family)))
    intercepts <- names(get_intercepts(ee, data = data, family = family))
    prior_fixef <- get_prior_fixef(fixef, intercepts = intercepts,
                                   internal = internal)
    monef <- colnames(get_model_matrix(ee$mono, data = data))
    prior_monef <- get_prior_monef(monef, fixef = fixef)
    csef <- colnames(get_model_matrix(ee$cse, data = data))
    prior_csef <- get_prior_csef(csef, fixef = fixef)
    prior_ranef <- get_prior_ranef(ranef, def_scale_prior, 
                                   internal = internal)
    prior <- rbind(prior, prior_fixef, prior_monef, prior_csef, prior_ranef)
  }
  # handle additional parameters
  is_ordinal <- is.ordinal(family)
  is_linear <- is.linear(family)
  nresp <- length(ee$response)
  cbound <- "<lower=-1,upper=1>"
  if (get_ar(autocor)) 
    prior <- rbind(prior, prior_frame(class = "ar", bound = cbound))
  if (get_ma(autocor)) 
    prior <- rbind(prior, prior_frame(class = "ma", bound = cbound))
  if (get_arr(autocor)) 
    prior <- rbind(prior, prior_frame(class = "arr"))
  if (has_sigma(family, se = is.formula(ee$se), autocor = autocor)) {
    sigma_prior <- prior_frame(class = "sigma", coef = c("", ee$response),
                               prior = c(def_scale_prior, rep("", nresp)))
    prior <- rbind(prior, sigma_prior)
  }
  if (is_linear && nresp > 1L) {
    if (internal) {
      prior <- rbind(prior, prior_frame(class = "Lrescor", 
                                        prior = "lkj_corr_cholesky(1)"))
    } else {
      prior <- rbind(prior, prior_frame(class = "rescor", prior = "lkj(1)"))
    }
  }
  if (family$family == "student") {
    prior <- rbind(prior, prior_frame(class = "nu", prior = "gamma(2, 0.1)"))
  }
  if (family$family == "beta") {
    prior <- rbind(prior, prior_frame(class = "phi", 
                                      prior = "gamma(0.01, 0.01)"))
  }
  if (has_shape(family)) {
    prior <- rbind(prior, prior_frame(class = "shape", 
                                      prior = "gamma(0.01, 0.01)"))
  }
  if (is_ordinal && threshold == "equidistant") {
    prior <- rbind(prior, prior_frame(class = "delta"))
  }
  # do not remove unique(.)
  prior <- unique(prior[with(prior, order(nlpar, class, group, coef)), ])
  rownames(prior) <- 1:nrow(prior)
  prior
}

get_prior_fixef <- function(fixef, intercepts = "Intercept", 
                            nlpar = "", internal = FALSE) {
  # priors for fixed effects parameters
  # Args:
  #   fixef: names of the fixed effects
  #   intercepts: names of the fixed effects Intercept(s)
  #   nlpar: optional name of a non-linear parameter
  #   internal: see get_prior
  # Returns:
  #   an object of class prior_frame
  prior <- empty_prior_frame()
  if (length(fixef)) {
    prior <- rbind(prior, prior_frame(class = "b", coef = c("", fixef),
                                      nlpar = nlpar)) 
  }
  if (length(intercepts)) {
    int_coefs <- "" 
    if (!is_equal(intercepts, "Intercept")) {
      int_coefs <- c(int_coefs, intercepts)
    }
    prior <- rbind(prior, prior_frame(class = "Intercept", coef = int_coefs,
                                      nlpar = nlpar))
    if (internal) {
      prior <- rbind(prior, prior_frame(class = "temp_Intercept",
                                        coef = int_coefs, nlpar = nlpar))
    }
  }
  prior
}

get_prior_monef <- function(monef, fixef = NULL, nlpar = "") {
  # priors for monotonous effects parameters
  # Args:
  #   monef: names of the monotonous effects
  #   fixef: names of the fixed effects
  #   nlpar: optional name of a non-linear parameter
  # Returns:
  #   an object of class prior_frame
  prior <- empty_prior_frame()
  if (length(monef)) {
    invalid <- intersect(fixef, monef)
    if (length(invalid)) {
      stop(paste("Variables cannot be modeled as fixed and", 
                 "monotonous effects at the same time.", 
                 "\nError occured for variables:", 
                 paste(invalid, collapse = ", ")), call. = FALSE)
    }
    prior <- rbind(prior_frame(class = "b", coef = c("", monef), nlpar = nlpar),
                   prior_frame(class = "simplex", coef = monef, nlpar = nlpar))
  }
  prior
}

get_prior_csef <- function(csef, fixef = NULL) {
  # priors for category spcific effects parameters
  # Args:
  #   csef: names of the category specific effects
  #   fixef: names of the fixed effects
  # Returns:
  #   an object of class prior_frame
  prior <- empty_prior_frame()
  if (length(csef)) {
    invalid <- intersect(fixef, csef)
    if (length(invalid)) {
      stop(paste("Variables cannot be modeled as fixed and", 
                 "category specific effects at the same time.", 
                 "\nError occured for variables:", 
                 paste(invalid, collapse = ", ")), call. = FALSE)
    }
    prior <- prior_frame(class = "b", coef = c("", csef))
  }
  prior
}

get_prior_ranef <- function(ranef, def_scale_prior, nlpar = "", 
                            internal = FALSE) {
  # priors for random effects parameters
  # Args:
  #   ranef: a list returned by gather_ranef
  #   def_scale_prior: a character string defining the default
  #                    prior for random effects SDs
  #   nlpar: optional name of a non-linear parameter
  #   internal: see get_prior
  # Returns:
  #   an object of class prior_frame
  if (nchar(nlpar)) {
    # extract only the relevant random effects
    ranef <- rmNULL(lapply(ranef, function(y) 
      if (identical(attr(y, "nlpar"), nlpar)) y else NULL))
  }
  prior <- empty_prior_frame()
  if (length(ranef)) {
    # global sd class
    prior <- rbind(prior, prior_frame(class = "sd", nlpar = nlpar,
                                      prior = def_scale_prior))
    gs <- names(ranef)
    for (i in seq_along(ranef)) {
      # include random effects standard deviations
      prior <- rbind(prior, prior_frame(class = "sd", coef = c("", ranef[[i]]),
                                        group = gs[i], nlpar = nlpar))
      # detect duplicated random effects
      J <- with(prior, class == "sd" & group == gs[i] & 
                       nlpar == nlpar & nchar(coef))
      dupli <- duplicated(prior[J, ])
      if (any(dupli)) {
        stop(paste("Duplicated random effects detected for group", gs[i]),
             call. = FALSE)
      }
      # include correlation parameters
      if (attr(ranef[[i]], "cor") && length(ranef[[i]]) > 1L) {
        if (internal) {
          prior <- rbind(prior, 
            prior_frame(class = "L", group = c("", gs[i]), nlpar = nlpar,
                        prior = c("lkj_corr_cholesky(1)", "")))
        } else {
          prior <- rbind(prior, 
            prior_frame(class = "cor", group = c("", gs[i]),
                        nlpar = nlpar, prior = c("lkj(1)", "")))
        }
      }
    }
  } 
  prior
}

check_prior <- function(prior, formula, data = NULL, family = gaussian(), 
                        sample_prior = FALSE, autocor = NULL, nonlinear = NULL, 
                        threshold = "flexible", check_rows = NULL) {
  # check prior input and amend it if needed
  #
  # Args:
  #   same as the respective parameters in brm
  #   check_rows: if not NULL, check only the rows given in check_rows
  #
  # Returns:
  #   a data.frame of prior specifications to be used in stan_prior (see stan.R)
  if (isTRUE(attr(prior, "checked"))) {
    return(prior)  # prior has already been checked; no need to do it twice
  }
  ee <- extract_effects(formula, family = family, nonlinear = nonlinear)  
  all_priors <- get_prior(formula = formula, data = data, 
                          family = family, autocor = autocor, 
                          threshold = threshold, nonlinear = nonlinear, 
                          internal = TRUE)
  if (is.null(prior)) {
    prior <- all_priors  
  } else {
    prior <- as.prior_frame(prior)
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
  has_specef <- is.formula(ee[c("mono", "cse")])
  temp <- handle_special_priors(prior, has_specef = has_specef)  
  prior <- temp$prior
  attrib <- temp$attrib 
  # check if parameters in prior are valid
  if (nrow(prior)) {
    valid <- which(duplicated(rbind(all_priors[, 2:5], prior[, 2:5])))
    invalid <- which(!1:nrow(prior) %in% (valid - nrow(all_priors)))
    if (length(invalid)) {
      msg_priors <- lapply(as.brmsprior(prior[invalid, ]), .print_prior)
      message(paste("The following priors don't correspond to any", 
                    "model parameter \nand will thus not affect the results:",
                    collapse("  \n", msg_priors)), "\n")
      prior <- prior[-invalid, ]
    }
  }
  check_prior_content(prior, family = family)
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
  if (is.formula(ee$mono)) {
    monef <- colnames(get_model_matrix(ee$mono, data = data))
    for (i in seq_along(monef)) {
      take <- with(prior, class == "simplex" & coef == monef[i])
      simplex_prior <- paste0(".", prior$prior[take])
      if (nchar(simplex_prior) > 1L) {
        simplex_prior <- paste(eval(parse(text = simplex_prior)),
                               collapse = ",")
        prior$prior[take] <- paste0("dirichlet(c(", simplex_prior, "))")
      }
    }
  }
  # check if priors for non-linear parameters are defined
  if (length(nonlinear)) {
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
  prior <- prior[with(prior, order(class, group, coef)), ]
  prior <- rbind(prior, prior_incr_lp)
  rownames(prior) <- 1:nrow(prior)
  # add attributes to prior generated in handle_special_priors
  for (i in seq_along(attrib)) {
    attr(prior, names(attrib)[i]) <- attrib[[i]]
  }
  attr(prior, "prior_only") <- identical(sample_prior, "only")
  attr(prior, "checked") <- TRUE
  prior
}

check_prior_content <- function(prior, family = gaussian()) {
  # try to check if prior distributions are reasonable
  # Args:
  #  prior: A prior_frame
  #  family: the model family
  stopifnot(is(prior, "prior_frame"))
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
    lb_pars <- c("sd", "sigma", "nu", "shape", "phi",
                 if (family %in% "cumulative") "delta")
    cor_pars <- c("cor", "L", "rescor", "Lrescor")
    autocor_pars <- c("ar", "ma")
    lb_warning <- ub_warning <- ""
    autocor_warning <- FALSE
    for (i in 1:nrow(prior)) {
      msg_prior <- .print_prior(as.brmsprior(prior[i, , drop = FALSE])[[1]])
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
    if (nchar(lb_warning)) {
      warning(paste0("It appears that you have specified a lower bounded ", 
                     "prior on a parameter that has no natural lower bound.",
                     "\nIf this is really what you want, please specify ",
                     "argument 'lb' of 'set_prior' appropriately.",
                     "\nWarning occurred for prior \n", lb_warning), 
              call. = FALSE)
    }
    if (nchar(ub_warning)) {
      warning(paste0("It appears that you have specified an upper bounded ", 
                     "prior on a parameter that has no natural upper bound.",
                     "\nIf this is really what you want, please specify ",
                     "argument 'ub' of 'set_prior' appropriately.",
                     "\nWarning occurred for prior \n", ub_warning), 
              call. = FALSE)
    }
    if (autocor_warning) {
      warning(paste("Changing the boundaries of autocorrelation", 
                    "parameters is not recommended."), call. = FALSE)
    }
  }
  invisible(NULL)
}

handle_special_priors <- function(prior, has_specef = FALSE) {
  # look for special priors such as horseshoe and process them appropriately
  #
  # Args:
  #   prior: an object of class prior_frame
  #   has_specef: are monotonous or category specific effects present?
  #
  # Returns:
  #   an named list of two objects: 
  #   prior: an updated version of prior
  #   attrib: a named list containing future attributes of prior
  attrib <- list()
  b_index <- which(prior$class == "b" & !nchar(prior$coef))
  if (length(b_index) && grepl("^horseshoe\\(.+\\)$", prior$prior[b_index])) {
    # horseshoe prior for fixed effects parameters
    if (any(nchar(prior$nlpar))) {
      stop("Horseshoe priors are not yet allowed in non-linear models.",
           call. = FALSE)
    }
    if (has_specef) {
      stop(paste("Horseshoe priors are not yet allowed in models with", 
                 "monotonous or category specific effects."), 
           call. = FALSE)
    }
    hs_df <- gsub("^horseshoe\\(|\\)$", "", prior$prior[b_index])
    hs_df <- suppressWarnings(as.numeric(hs_df))
    if (!is.na(hs_df) && hs_df > 0) {
      b_coef_indices <- which(prior$class == "b" & nchar(prior$coef)
                              & prior$coef != "Intercept")
      if (any(nchar(prior$prior[b_coef_indices]))) {
        stop(paste("Defining priors for single fixed effects parameters",
                   "is not allowed when using horseshoe priors",
                   "(except for the Intercept)"), call. = FALSE)
      }
      attrib$hs_df <- hs_df
      prior$prior[b_index] <- "normal(0, hs_local * hs_global)"
    } else {
      stop("degrees of freedom of horseshoe prior must be a positive number",
           call. = FALSE)
    }
  }
  # expand lkj correlation prior to full name
  prior$prior <- sub("^(lkj\\(|lkj_corr\\()", "lkj_corr_cholesky(", prior$prior)
  list(prior = prior, attrib = attrib)
}

get_bound <- function(prior, class = "b", coef = "", 
                      group = "", nlpar = "") {
  # extract the boundaries of a parameter described by class etc.
  # Args:
  #   prior: object of class prior_frame5
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

prior_frame <- function(prior = "", class = "", coef = "", 
                        group = "", nlpar = "", bound = "") {
  # helper function to create data.frames containing prior information 
  out <- data.frame(prior = prior, class = class, coef = coef, 
                    group = group, nlpar = nlpar, bound = bound,
                    stringsAsFactors = FALSE)
  class(out) <- c("prior_frame", "data.frame")
  out
}

empty_prior_frame <- function() {
  # define a prior_frame with zero rows
  prior_frame(prior = character(0), class = character(0), 
              coef = character(0), group = character(0),
              nlpar = character(0), bound = character(0))
}

update_prior_frame <- function(object, ranef = list(), ...) {
  # update prior_frames of models fitted with brms <= 0.8.0
  # Args:
  #   object: an object of class prior_frame
  #   ranef: a named list of group specific terms
  #   ...: currently ignored
  has_group <- nchar(object$group) > 0
  num_group <- suppressWarnings(as.numeric(object$group[has_group]))
  if (length(num_group) && !anyNA(num_group)) {
    if (max(num_group) != length(ranef)) {
      warning(paste("Priors for standard deviation and correlation", 
                    "parameters of group specific terms cannot be upated.", 
                    "\nReturning to default priors."), call. = FALSE)
    } else {
      object$group[has_group] <- names(ranef)[num_group]
    }
  }
  attr(object, "checked") <- NULL
  object
}

#' @export
print.brmsprior <- function(x, ...) {
  cat(.print_prior(x))
  invisible(x)
}

.print_prior <- function(x) {
  # prepare text for print.brmsprior
  group <- ifelse(nchar(x$group), paste0("_", x$group), "")
  coef <- ifelse(nchar(x$coef), paste0("_", x$coef), "")
  nlpar <- ifelse(nchar(x$nlpar), paste0("_", x$nlpar), "")
  bound <- ifelse(nchar(x$bound), paste0(x$bound, " "), "")
  tilde <- ifelse(nchar(x$class) + nchar(group) + nchar(coef), " ~ ", "")
  prior <- ifelse(nchar(x$prior), x$prior, "(no prior)")
  paste0(bound, x$class, nlpar, group, coef, tilde, prior)
}

#' @export
c.brmsprior <- function(x, ...) {
  # combines multiple brmsprior objects into one prior_frame
  if(any(!sapply(list(...), is, class2 = "brmsprior")))
    stop("All arguments must be of class brmsprior")
  prior <- data.frame(matrix(unlist(list(x, ...)), ncol = 6, byrow = TRUE),
                      stringsAsFactors = FALSE)
  names(prior) <- c("prior", "class", "coef", "group", "nlpar", "bound") 
  class(prior) <- c("prior_frame", "data.frame")
  prior
}

as.brmsprior <- function(prior) {
  # convert a prior_frame into a list of brmsprior objects
  # Args:
  #   prior: an object of class 'prior_frame' or 'brmsprior'
  stopifnot(is(prior, "prior_frame") || is(prior, "brmsprior"))
  if (is(prior, "prior_frame")) {
    .convert <- function(x) {
      structure(as.list(x), class = c("brmsprior", "list"))
    } 
    prior <- unname(apply(prior, MARGIN = 1, FUN = .convert))
  }
  prior
}

as.prior_frame <- function(prior) {
  # convert a brmsprior object into a prior_frame object
  # Args:
  #   prior: an object of class 'prior_frame' or 'brmsprior'
  if (is.null(prior)) {
    prior <- prior_frame()
  } else if (is(prior, "brmsprior")) {
    prior <- c(prior)
  } else if (!is(prior, "prior_frame")) {
    stop(paste("Invalid 'prior' argument. See help(set_prior)", 
               "for further information."), call. = FALSE)
  }
  prior
}

.dirichlet <- function(...) {
  # helper function for dirichlet priors of simplex parameters
  out <- as.numeric(c(...))
  if (any(out <= 0)) {
    stop("The dirichlet prior expects positive values.", call. = FALSE)
  }
  out
}