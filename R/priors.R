#' Prior Definitions for \pkg{brms} Models
#'
#' Define priors for specific parameters or classes of parameters
#' 
#' @aliases brmsprior brmsprior-class
#'
#' @param prior A character string defining a distribution in \pkg{Stan} language
#' @param class The parameter class. Defaults to \code{"b"} 
#'   (i.e. population-level effects). 
#'   See 'Details' for other valid parameter classes. 
#' @param coef Name of the (population- or group-level) parameter.  
#' @param group Grouping factor of group-level parameters.
#' @param resp Name of the response variable / category.
#'   Only used in multivariate and categorical models.
#' @param dpar Name of a distributional parameter.
#'   Only used in distributional models.
#' @param nlpar Name of a non-linear parameter. 
#'   Only used in non-linear models.
#' @param lb Lower bound for parameter restriction. Currently only allowed
#'   for classes \code{"b"}, \code{"ar"}, \code{"ma"}, and \code{"arr"}.
#'   Defaults to \code{NULL}, that is no restriction.
#' @param ub Upper bound for parameter restriction. Currently only allowed
#'   for classes \code{"b"}, \code{"ar"}, \code{"ma"}, and \code{"arr"}.
#'   Defaults to \code{NULL}, that is no restriction.
#' @param check Logical; Indicates whether priors
#'   should be checked for validity (as far as possible).
#'   Defaults to \code{TRUE}. If \code{FALSE}, \code{prior} is passed
#'   to the Stan code as is, and all other arguments are ignored.
#' @param ... Arguments passed to \code{set_prior}.
#' 
#' @return An object of class \code{brmsprior} to be used in the \code{prior}
#'   argument of \code{\link[brms:brm]{brm}}.
#' 
#' @details 
#'   \code{set_prior} is used to define prior distributions for parameters 
#'   in \pkg{brms} models. The functions \code{prior}, \code{prior_}, and
#'   \code{prior_string} are aliases of \code{set_prior} each allowing
#'   for a differnt kind of argument specification. 
#'   \code{prior} allows specifying arguments as expression without
#'   quotation marks using non-standard evaluation. 
#'   \code{prior_} allows specifying arguments as one-sided formulas
#'   or wrapped in \code{quote}.
#'   \code{prior_string} allows specifying arguments as strings just
#'   as \code{set_prior} itself.
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
#    These parameters are internally named as \code{b_<coef>}, where \code{<coef>} 
#'   represents the name of the corresponding population-level effect. 
#'   Suppose, for instance, that \code{y} is predicted by \code{x1} and \code{x2} 
#'   (i.e., \code{y ~ x1 + x2} in formula syntax). 
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
#'   To put the same prior on all population-level effects at once, 
#'   we may write as a shortcut \code{set_prior("<prior>", class = "b")}. 
#'   This also leads to faster sampling, because priors can be vectorized in this case. 
#'   Both ways of defining priors can be combined using for instance 
#'   \code{set_prior("normal(0,2)", class = "b")} and \cr
#'   \code{set_prior("normal(0,10)", class = "b", coef = "x1")}
#'   at the same time. This will set a \code{normal(0,10)} prior on 
#'   the effect of \code{x1} and a \code{normal(0,2)} prior 
#'   on all other population-level effects. 
#'   However, this will break vectorization and
#'   may slow down the sampling procedure a bit.
#'   
#'   In case of the default intercept parameterization 
#'   (discussed in the 'Details' section of 
#'   \code{\link[brms:brmsformula]{brmsformula}}),
#'   general priors on class \code{"b"} will \emph{not} affect 
#'   the intercept. Instead, the intercept has its own parameter class 
#'   named \code{"Intercept"} and priors can thus be 
#'   specified via \code{set_prior("<prior>", class = "Intercept")}.
#'   Setting a prior on the intercept will not break vectorization
#'   of the other population-level effects.
#'   Note that technially, this prior is set on an intercept that
#'   results when internally centering all population-level predictors 
#'   around zero to improve sampling efficiency. On this centered 
#'   intercept, specifying a prior is actually much easier and 
#'   intuitive than on the original intercept, since the former 
#'   represents the expected response value when all predictors 
#'   are at their means. To treat the intercept as an ordinary 
#'   population-level effect and avoid the centering parameterization, 
#'   use \code{0 + intercept} on the right-hand side of the model formula.
#'   
#'   A special shrinkage prior to be applied on population-level effects 
#'   is the horseshoe prior. See \code{\link[brms:horseshoe]{horseshoe}}
#'   for details. Another shrinkage prior is the so-called lasso prior.
#'   See \code{\link[brms:lasso]{lasso}} for details.
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
#'   \code{sd_<group>_<coef>}. Consider, for instance, the formula 
#'   \code{y ~ x1 + x2 + (1 + x1 | g)}.
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
#'   When defining priors on group-level parameters in non-linear models, 
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
#'   is essentially the only prior for (Cholesky factors) of correlation matrices. 
#'   If \code{eta = 1} (the default) all correlations matrices 
#'   are equally likely a priori. If \code{eta > 1}, extreme correlations 
#'   become less likely, whereas \code{0 < eta < 1} results in 
#'   higher probabilities for extreme correlations. 
#'   Correlation matrix parameters in \code{brms} models are named as 
#'   \code{cor_<group>}, (e.g., \code{cor_g} if \code{g} is the grouping factor).
#'   To set the same prior on every correlation matrix, 
#'   use for instance \code{set_prior("lkj(2)", class = "cor")}.
#'   Internally, the priors are transformed to be put on the Cholesky factors
#'   of the correlation matrices to improve efficiency and numerical stability.
#'   The corresponding parameter class of the Cholesky factors is \code{L},
#'   but it is not recommended to specify priors for this parameter class directly.
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
#'   in the same way as the group-level standard deviations.
#'   Furthermore, family \code{student} needs the parameter 
#'   \code{nu} representing the degrees of freedom of students-t distribution. 
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
#' Gelman A. (2006). Prior distributions for variance parameters in hierarchical models.
#'    Bayesian analysis, 1(3), 515 -- 534.
#' 
#' @examples
#' ## use alias functions
#' (prior1 <- prior(cauchy(0, 1), class = sd))
#' (prior2 <- prior_(~cauchy(0, 1), class = ~sd))
#' (prior3 <- prior_string("cauchy(0, 1)", class = "sd"))
#' identical(prior1, prior2)
#' identical(prior1, prior3)
#' 
#' ## check which parameters can have priors
#' get_prior(rating ~ treat + period + carry + (1|subject),
#'           data = inhaler, family = cumulative())
#'          
#' ## define some priors          
#' prior <- c(prior_string("normal(0,10)", class = "b"),
#'            prior(normal(1,2), class = b, coef = treat),
#'            prior_(~cauchy(0,2), class = ~sd, 
#'                   group = ~subject, coef = ~Intercept))
#'               
#' ## verify that the priors indeed found their way into Stan's model code
#' make_stancode(rating ~ treat + period + carry + (1|subject),
#'               data = inhaler, family = cumulative(),
#'               prior = prior)
#'               
#' ## use the horseshoe prior to model sparsity in population-level effects
#' make_stancode(count ~ log_Age_c + log_Base4_c * Trt_c,
#'               data = epilepsy, family = poisson(),
#'               prior = set_prior("horseshoe(3)"))
#'               
#' ## alternatively use the lasso prior
#' make_stancode(count ~ log_Age_c + log_Base4_c * Trt_c,
#'               data = epilepsy, family = poisson(),
#'               prior = set_prior("lasso(1)"))
#' 
#' ## pass priors to Stan without checking
#' prior <- prior_string("target += normal_lpdf(b[1] | 0, 1)", check = FALSE)
#' make_stancode(count ~ Trt_c, data = epilepsy, prior = prior)
#'
#' @export
set_prior <- function(prior, class = "b", coef = "", group = "",
                      resp = "", dpar = "", nlpar = "", 
                      lb = NULL, ub = NULL, check = TRUE) {
  prior <- as.character(prior)
  class <- as.character(class)
  group <- as.character(group)
  coef <- as.character(coef)
  resp <- as.character(resp)
  dpar <- as.character(dpar)
  nlpar <- as.character(nlpar)
  lb <- as.numeric(lb)
  ub <- as.numeric(ub)
  check <- as_one_logical(check)
  if (length(prior) != 1 || length(class) != 1 || length(coef) != 1 || 
      length(group) != 1 || length(resp) != 1 || length(dpar) != 1 ||
      length(nlpar) != 1 || length(lb) > 1 || length(ub) > 1) {
    stop2("All arguments of set_prior must be of length 1.")
  }
  # validate boundaries
  is_arma <- class %in% c("ar", "ma")
  if (length(lb) || length(ub) || is_arma) {
    if (!(class %in% c("b", "ar", "ma", "arr"))) {
      stop2(
        "Currently boundaries are only allowed for ", 
        "population-level and autocorrelation parameters."
      ) 
    }
    if (nzchar(coef)) {
      stop2("Argument 'coef' may not be specified when using boundaries.")
    }
    if (is_arma) {
      lb <- ifelse(length(lb), lb, -1)
      ub <- ifelse(length(ub), ub, 1) 
      if (is.na(lb) || is.na(ub) || abs(lb) > 1 || abs(ub) > 1) {
        warning2(
          "Setting boundaries of autocorrelation parameters ", 
          "outside of [-1,1] may not be appropriate."
        )
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
  if (!check) {
    # prior will be added to the log-posterior as is
    class <- coef <- group <- resp <- dpar <- nlpar <- bound <- ""
  }
  do.call(brmsprior, 
    nlist(prior, class, coef, group, resp, dpar, nlpar, bound)
  )
}

#' @describeIn set_prior Alias of \code{set_prior} allowing to 
#'   specify arguments as expressions without quotation marks.
#' @export
prior <- function(prior, ...) {
  call <- as.list(match.call()[-1])
  seval <- rmNULL(call[prior_seval_args()])
  call[prior_seval_args()] <- NULL
  call <- lapply(call, deparse_no_string)
  do.call(set_prior, c(call, seval))
}

#' @describeIn set_prior Alias of \code{set_prior} allowing to specify 
#'   arguments as as one-sided formulas or wrapped in \code{quote}.
#' @export
prior_ <- function(prior, ...) {
  call <- nlist(prior, ...)
  seval <- rmNULL(call[prior_seval_args()])
  call[prior_seval_args()] <- NULL
  as_string <- function(x) {
    if (is.formula(x) && length(x) == 2) {
      deparse_no_string(x[[2]])
    } else if (is.call(x) || is.name(x) || is.atomic(x)) {
      deparse_no_string(x)
    } else {
      stop2("Arguments must be one-sided formula, call, name, or constant.") 
    }
  }
  call <- lapply(call, as_string)
  do.call(set_prior, c(call, seval))
}

prior_seval_args <- function() {
  # arguments for which to use standard evaluation
  c("check")
}

#' @describeIn set_prior Alias of \code{set_prior} allowing to
#'   specify arguments as strings.
#' @export
prior_string <- function(prior, ...) {
  set_prior(prior, ...)
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
get_prior <- function(formula, data, family = gaussian(),
                      autocor = NULL, nonlinear = NULL,
                      threshold = c("flexible", "equidistant"), 
                      internal = FALSE) {
  # note that default priors are stored in this function
  formula <- amend_formula(
    formula, data = data, family = family, 
    autocor = autocor, threshold = threshold,
    nonlinear = nonlinear
  )
  family <- formula$family
  autocor <- formula$autocor
  bterms <- parse_bf(formula)
  data <- update_data(data, bterms = bterms)
  ranef <- tidy_ranef(bterms, data)
  
  # ensure that RE and residual SDs only have 
  # a weakly informative prior by default
  Y <- unname(model.response(data))
  prior_scale <- 10
  link <- family$link
  if (is_lognormal(family)) {
    link <- "log"
  }
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
  if (length(bterms$response) > 1L) {
    # priors for effects in multivariate models
    for (r in bterms$response) {
      bterms$dpars[["mu"]]$resp <- r
      prior_eff <- prior_effects(
        bterms$dpars[["mu"]], data = data,
        def_scale_prior = def_scale_prior,
        internal = internal
      )
      prior <- prior + prior_eff
    }
    bterms$dpars[["mu"]] <- NULL
    # add "global" priors for population-level effects
    # in 1.8.0 as users keep asking about this
    for (cl in c("b", "Intercept")) {
      if (any(with(prior, class == cl & coef == ""))) {
        prior <- prior + brmsprior(class = cl) 
      }
    }
  }
  # priors for distributional parameters
  def_auxprior <- c(
    sigma = def_scale_prior, 
    shape = "gamma(0.01, 0.01)",
    nu = "gamma(2, 0.1)", 
    phi = "gamma(0.01, 0.01)",
    kappa = "gamma(2, 0.01)", 
    beta = "gamma(1, 0.1)", 
    zi = "beta(1, 1)", 
    hu = "beta(1, 1)", 
    zoi = "beta(1, 1)",
    coi = "beta(1, 1)",
    bs = "gamma(1, 1)", 
    ndt = "uniform(0, min_Y)", 
    bias = "beta(1, 1)", 
    quantile = "beta(1, 1)",
    xi = "normal(0, 2.5)",
    alpha = "normal(0, 4)",
    disc = NA,
    mu = NA
  )
  valid_dpars <- valid_dpars(family, bterms = bterms)
  for (dp in valid_dpars) {
    dp_class <- dpar_class(dp)
    if (!is.null(bterms$dpars[[dp]])) {
      dp_prior <- prior_effects(
        bterms$dpars[[dp]], data = data,
        def_scale_prior = def_scale_prior
      )
    } else if (!is.na(def_auxprior[dp_class])) {
      dp_prior <- brmsprior(class = dp, prior = def_auxprior[dp_class])
    } else {
      dp_prior <- empty_brmsprior()
    }
    prior <- prior + dp_prior
  }
  # priors of group-level parameters
  prior_re <- prior_re(
    ranef, def_scale_prior = def_scale_prior,
    global_sd = length(bterms$response) > 1L,
    internal = internal
  )
  prior <- prior + prior_re
  # prior for the delta parameter for equidistant thresholds
  if (is_ordinal(family) && is_equal(family$threshold, "equidistant")) {
    bound <- ifelse(family$family == "cumulative", "<lower=0>", "")
    prior <- prior + brmsprior(class = "delta", bound = bound)
  }
  # priors for mixture models
  ap_classes <- dpar_class(names(c(bterms$dpars, bterms$fdpars)))
  if (is.mixfamily(family) && !any(ap_classes == "theta")) {
    prior <- prior + brmsprior(class = "theta")
  }
  # priors for distributional parameters of multivariate models
  if (is_linear(family) && length(bterms$response) > 1L) {
    sigma_coef <- c("", bterms$response)
    sigma_prior <- c(def_scale_prior, rep("", length(bterms$response)))
    sigma_prior <- brmsprior(
      class = "sigma", coef = sigma_coef, prior = sigma_prior
    )
    prior <- prior + sigma_prior
    if (internal) {
      prior <- prior + 
        brmsprior(class = "Lrescor", prior = "lkj_corr_cholesky(1)")
    } else {
      prior <- prior + brmsprior(class = "rescor", prior = "lkj(1)")
    }
  }
  # priors for autocor parameters
  cbound <- "<lower=-1,upper=1>"
  if (get_ar(autocor)) {
    prior <- prior + brmsprior(class = "ar", bound = cbound)
  }
  if (get_ma(autocor)) {
    prior <- prior + brmsprior(class = "ma", bound = cbound)
  }
  if (get_arr(autocor)) {
    prior <- prior + brmsprior(class = "arr")
  }
  if (is.cor_sar(autocor)) {
    if (identical(autocor$type, "lag")) {
      prior <- prior + brmsprior(class = "lagsar")
    }
    if (identical(autocor$type, "error")) {
      prior <- prior + brmsprior(class = "errorsar")
    }
  }
  if (is.cor_car(autocor)) {
    prior <- prior +  
      brmsprior(def_scale_prior, class = "sdcar")
    if (identical(autocor$type, "escar")) {
      prior <- prior + brmsprior(class = "car")
    }
  }
  if (is(autocor, "cor_bsts")) {
    prior <- prior +
      brmsprior(class = "sigmaLL", prior = def_scale_prior)
  }
  # do not remove unique(.)
  prior <- unique(prior[with(prior, 
    order(resp, dpar, nlpar, class, group, coef)
  ), ])
  rownames(prior) <- NULL
  structure(prior, class = c("brmsprior", "data.frame"))
}

#' @export
prior_effects.btl <- function(x, data, spec_intercept = TRUE,
                              def_scale_prior = "", ...) {
  # collect default priors for various kinds of effects
  # Args:
  #   spec_intercept: special parameter class for the Intercept?
  #   def_scale_prior: default prior for SD parameters
  # Return:
  #   An object of class brmsprior 
  prior_fe(x, data, spec_intercept = spec_intercept) +
    prior_cs(x, data) +
    prior_mo(x, data) +
    prior_sm(x, data, def_scale_prior = def_scale_prior) + 
    prior_me(x, data) + 
    prior_gp(x, data, def_scale_prior = def_scale_prior)
}

#' @export
prior_effects.btnl <- function(x, data, def_scale_prior = "", ...) {
  # collect default priors for non-linear parameters
  # Args:
  #   see prior_effects.btl
  nlpars <- names(x$nlpars)
  prior <- empty_brmsprior()
  for (i in seq_along(nlpars)) {
    prior_eff <- prior_effects(
      x$nlpars[[i]], data = data, 
      def_scale_prior = def_scale_prior,
      spec_intercept = FALSE
    )
    prior <- prior + prior_eff
  }
  prior
}

prior_fe <- function(bterms, data, spec_intercept = TRUE) {
  # priors for population-level parameters
  # Args:
  #   spec_intercept: special parameter class for the Intercept? 
  # Returns:
  #   an object of class brmsprior
  prior <- empty_brmsprior()
  fixef <- colnames(data_fe(bterms, data)$X)
  px <- check_prefix(bterms)
  if (has_intercept(bterms$fe) && spec_intercept) {
    prior <- prior +
      brmsprior(class = "Intercept", coef = "", ls = px) +
      brmsprior(class = "b", coef = "Intercept", ls = px)
    fixef <- setdiff(fixef, "Intercept")
  }
  if (length(fixef)) {
    prior <- prior + 
      brmsprior(class = "b", coef = c("", fixef), ls = px)
  }
  prior
}

prior_mo <- function(bterms, data) {
  # priors for monotonic effects parameters
  # Returns:
  #   an object of class brmsprior
  prior <- empty_brmsprior()
  monef <- all_terms(bterms$mo)
  if (length(monef)) {
    px <- check_prefix(bterms)
    prior <- prior + 
      brmsprior(class = "b", coef = c("", monef), ls = px) + 
      brmsprior(class = "simplex", coef = monef, ls = px)
  }
  prior
}

prior_cs <- function(bterms, data) {
  # priors for category spcific effects parameters
  # Returns:
  #   an object of class brmsprior
  prior <- empty_brmsprior()
  csef <- colnames(get_model_matrix(bterms$cs, data = data))
  if (length(csef)) {
    px <- check_prefix(bterms)
    prior <- prior + 
      brmsprior(class = "b", coef = c("", csef), ls = px)
  }
  prior
}

prior_me <- function(bterms, data) {
  # default priors of coefficients of noisy terms
  # Returns:
  #   an object of class brmsprior
  prior <- empty_brmsprior()
  meef <- get_me_labels(bterms, data)
  if (length(meef)) {
    px <- check_prefix(bterms)
    prior <- prior + 
      brmsprior(class = "b", coef = c("", rename(meef)), ls = px)
  }
  prior
}

prior_gp <- function(bterms, data, def_scale_prior) {
  # default priors of gaussian processes
  # Returns:
  #   an object of class brmsprior
  #   def_scale_prior: a character string defining 
  #     the default prior for random effects SDs
  prior <- empty_brmsprior()
  gpef <- get_gp_labels(bterms)
  if (length(gpef)) {
    px <- check_prefix(bterms)
    prior <- prior +
      brmsprior(class = "sdgp", prior = def_scale_prior, ls = px) +
      brmsprior(class = "sdgp", coef = gpef, ls = px) +
      brmsprior(class = "lscale", prior = "normal(0, 0.5)", ls = px) +
      brmsprior(class = "lscale", coef = gpef, ls = px)
  }
  prior
}

prior_re <- function(ranef, def_scale_prior, global_sd = FALSE,
                     internal = FALSE) {
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
  if (!nrow(ranef)) {
    return(prior)
  }
  # global sd class
  px <- check_prefix(ranef)
  upx <- unique(px)
  if (global_sd) {
    global_sd_prior <- rep("", nrow(upx))
    global_sd_prior <- c(def_scale_prior, global_sd_prior)
    upx <- lapply(upx, function(x) union("", x))
    global_sd_prior <- brmsprior(
      class = "sd", prior = global_sd_prior, ls = upx
    )
  } else {
    global_sd_prior <- brmsprior(
      class = "sd", prior = def_scale_prior, ls = px
    )
  }
  prior <- prior + global_sd_prior
  for (id in unique(ranef$id)) {
    r <- subset2(ranef, id = id)
    group <- r$group[1]
    rpx <- check_prefix(r)
    urpx <- unique(rpx)
    # include group-level standard deviations
    prior <- prior + 
      brmsprior(class = "sd", group = group, ls = urpx) +
      brmsprior(class = "sd", coef = r$coef, group = group, ls = rpx)
    # detect duplicated group-level effects
    J <- with(prior, class == "sd" & nzchar(coef))
    dupli <- duplicated(prior[J, ])
    if (any(dupli)) {
      stop2("Duplicated group-level effects detected for group ", group)
    }
    # include correlation parameters
    if (isTRUE(r$cor[1]) && nrow(r) > 1L) {
      if (internal) {
        prior <- prior +
          brmsprior(
            class = "L", group = c("", group),
            prior = c("lkj_corr_cholesky(1)", "")
          )
      } else {
        prior <- prior +
          brmsprior(
            class = "cor", group = c("", group),
            prior = c("lkj(1)", "")
          )
      }
    }
  }
  prior
}

prior_sm <- function(bterms, data, def_scale_prior) {
  # priors for smooth terms
  # Args:
  #   def_scale_prior: a character string defining 
  #     the default prior for smooth SDs
  prior <- empty_brmsprior()
  smooths <- get_sm_labels(bterms)
  if (length(smooths)) {
    px <- check_prefix(bterms)
    prior_strings <- c(def_scale_prior, rep("", length(smooths)))
    prior <- prior + brmsprior(
      class = "sds", coef = c("", smooths), 
      prior = prior_strings, ls = px
    )
  }
  prior
}

check_prior <- function(prior, formula, data = NULL, family = NULL, 
                        autocor = NULL, sample_prior = c("no", "yes", "only"),
                        check_rows = NULL, warn = FALSE) {
  # check prior input and amend it if needed
  # Args:
  #   same as the respective parameters in brm
  #   check_rows: if not NULL, check only the rows given in check_rows
  #   warn: passed to check_prior_content
  # Returns:
  #   a data.frame of prior specifications to be used in stan_prior (see stan.R)
  sample_prior <- check_sample_prior(sample_prior)
  if (isTRUE(attr(prior, "checked"))) {
    # prior has already been checked; no need to do it twice
    # attributes may still need to be updated
    attr(prior, "sample_prior") <- sample_prior
    return(prior)
  }
  formula <- bf(formula, family = family, autocor = autocor)
  family <- formula$family
  autocor <- formula$autocor
  bterms <- parse_bf(formula)
  all_priors <- get_prior(formula = formula, data = data, internal = TRUE)
  if (is.null(prior)) {
    prior <- all_priors  
  }
  # temporarily exclude priors that should not be checked
  no_checks <- !nzchar(prior$class)
  prior_no_checks <- prior[no_checks, ]
  prior <- prior[!no_checks, ]
  # check for duplicated priors
  prior$class <- rename(
    prior$class, symbols = c("^cor$", "^rescor$"), 
    subs = c("L", "Lrescor"), fixed = FALSE
  )
  duplicated_input <- duplicated(prior[, 2:7])
  if (any(duplicated_input)) {
    stop2("Duplicated prior specifications are not allowed.")
  }
  valid_dpars <- valid_dpars(family, bterms)
  nlpars_in_dpars <- prior$nlpar %in% valid_dpars
  if (any(nlpars_in_dpars)) {
    warning2(
      "Specifying priors of distributional parameters via ",
      "'nlpar' is deprecated. Please use 'dpar' instead."
    )
    prior$dpar[nlpars_in_dpars] <- prior$nlpar[nlpars_in_dpars]
    prior$nlpar[nlpars_in_dpars] <- ""
  }
  if (length(bterms$response) > 1L) {
    nlpars_in_resp <- prior$nlpar %in% bterms$response
    if (any(nlpars_in_resp)) {
      warning2(
        "Specifying priors in multivariate models via ",
        "'nlpar' is deprecated. Please use 'resp' instead."
      )
      prior$resp[nlpars_in_resp] <- prior$nlpar[nlpars_in_resp]
      prior$nlpar[nlpars_in_resp] <- ""
    }
  }
  # check if parameters in prior are valid
  if (nrow(prior)) {
    valid <- which(duplicated(rbind(all_priors[, 2:7], prior[, 2:7])))
    invalid <- which(!1:nrow(prior) %in% (valid - nrow(all_priors)))
    if (length(invalid)) {
      msg_priors <- .print_prior(prior[invalid, ])
      stop2(
        "The following priors do not correspond ", 
        "to any model parameter: \n",
        collapse(.print_prior(prior[invalid, ]), "\n")
      )
    }
  }
  prior$prior <- sub("^(lkj|lkj_corr)\\(", "lkj_corr_cholesky(", prior$prior)
  check_prior_content(prior, family = family, warn = warn)
  # check if priors for non-linear parameters are defined
  for (dp in names(bterms$dpars)) {
    nlpars <- names(bterms$dpars[[dp]]$nlpars)
    dp <- ifelse(dp == "mu", "", dp)
    for (nlp in nlpars) {
      nlp_prior <- subset2(prior, dpar = dp, nlpar = nlp, class = "b")
      if (!any(nzchar(nlp_prior$prior))) {
        stop2(
          "Priors on population-level effects are required in ",
          "non-linear models,but none were found for parameter ", 
          "'", nlp, "'. See help(set_prior) for more details."
        )
      }
    }
  }
  # prepare special priors for use in Stan
  prior <- check_prior_special(bterms, prior)
  # merge user-specified priors with default priors
  prior <- prior + all_priors
  prior <- prior[!duplicated(prior[, 2:7]), ]
  rows2remove <- NULL
  # copy over the global population-level prior in MV models
  if (length(bterms$response) > 1L) {
    for (cl in c("b", "Intercept")) {
      g_index <- find_rows(prior, class = cl, coef = "", resp = "")
      for (r in bterms$response) {
        r_index <- find_rows(prior, class = cl, coef = "", resp = r)
        if (isTRUE(!nzchar(prior$prior[r_index]))) {
          prior$prior[r_index] <- prior$prior[g_index]
        }
      }
      rows2remove <- c(rows2remove, which(g_index))
    }
  }
  # special treatment of population-level intercepts
  int_index <- which(prior$class == "Intercept")
  if (length(int_index)) {
    int_prior <- prior[int_index, ]
    bint_index <- with(prior, class == "b" & coef %in% "Intercept")
    bint_prior <- prior[bint_index, ]
    for (t in int_index) {
      tb <- match(prior$nlpar[t], bint_prior$nlpar) 
      if (!is.na(tb) && nzchar(bint_prior$prior[tb])) {
        # fall back to 'b' priors deprecated as of brms 1.5.0
        if (nzchar(prior$prior[t])) {
          stop2("Duplicated prior definitions detected ", 
                "for the population-level intercept.")
        } else {
          warning2(
            "Setting a prior on the population-level intercept",
            "\nvia (class = 'b', coef = 'Intercept') is deprecated.",
            "\nPlease use (class = 'Intercept', coef = '') instead."
          )
          prior$prior[t] <- bint_prior$prior[tb]
        }
      }
    }
    rows2remove <- c(rows2remove, which(bint_index))
  }
  # prepare priors of monotonic effects
  mo_forms <- get_effect(bterms, "mo")
  for (k in seq_along(mo_forms)) {
    monef <- colnames(get_model_matrix(mo_forms[[k]], data = data))
    for (i in seq_along(monef)) {
      take <- with(prior,
        class == "simplex" & coef == monef[i] & nlpar == names(mo_forms)[k]
      )
      simplex_prior <- prior$prior[take]
      if (isTRUE(nzchar(simplex_prior))) {
        # hard code prior concentration 
        # in order not to depend on external objects
        simplex_prior <- paste0(eval2(simplex_prior), collapse = ", ")
        prior$prior[take] <- paste0("dirichlet(c(", simplex_prior, "))")
      }
    }
  }
  # prepare priors for mixture probabilities
  if (is.mixfamily(family)) {
    take <- prior$class == "theta"
    theta_prior <- prior$prior[take]
    if (isTRUE(nzchar(theta_prior))) {
      # hard code prior concentration
      theta_prior <- paste0(eval2(theta_prior), collapse = ", ")
      prior$prior[take] <- paste0("dirichlet(c(", theta_prior, "))")
    }
  }
  if (length(rows2remove)) {   
    prior <- prior[-rows2remove, ]
  }
  prior <- prior[with(prior, order(resp, dpar, nlpar, class, group, coef)), ]
  prior <- prior + prior_no_checks
  rownames(prior) <- NULL
  attr(prior, "sample_prior") <- sample_prior
  attr(prior, "checked") <- TRUE
  prior
}

check_prior_content <- function(prior, family = gaussian(), warn = TRUE) {
  # try to check if prior distributions are reasonable
  # Args:
  #  prior: A brmsprior object
  #  family: the model family
  #  warn: logical; print boundary warnings?
  if (!is.brmsprior(prior)) {
    return(invisible(TRUE))
  }
  stopifnot(is.family(family))
  family <- family$family
  if (nrow(prior)) {
    lb_priors <- c(
      "lognormal", "chi_square", "inv_chi_square",
      "scaled_inv_chi_square", "exponential", "gamma",
      "inv_gamma", "weibull", "frechet", "rayleigh",
      "pareto", "pareto_type_2"
    )
    lb_priors_reg <- paste0("^(", paste0(lb_priors, collapse = "|"), ")")
    ulb_priors <- c("beta", "uniform", "von_mises")
    ulb_priors_reg <- paste0("^(", paste0(ulb_priors, collapse = "|"), ")")
    nb_pars <- c(
      "b", "Intercept", "alpha", "xi",
      if (!family %in% "cumulative") "delta"
    )
    lb_pars <- c(
      "sigma", "shape", "nu", "phi", "kappa", "beta", "bs", 
      "disc", "sdcar", "sigmaLL", "sd", "sds", "sdgp", "lscale", 
      if (family %in% "cumulative") "delta"
    )
    cor_pars <- c("cor", "L", "rescor", "Lrescor")
    autocor_pars <- c("ar", "ma")
    lb_warning <- ub_warning <- ""
    autocor_warning <- FALSE
    for (i in seq_len(nrow(prior))) {
      msg_prior <- .print_prior(prior[i, , drop = FALSE])
      has_lb_prior <- grepl(lb_priors_reg, prior$prior[i])
      has_ulb_prior <- grepl(ulb_priors_reg, prior$prior[i])
      # priors with nchar(coef) inherit their boundaries 
      j <- which(with(prior, 
        class == class[i] & group == group[i] & 
        nlpar == nlpar[i] & !nchar(coef)
      ))
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
          stop2(
            "Currently 'lkj' is the only valid prior ",
            "for group-level correlations. See help(set_prior) ",
            "for more details."
          )
        }
      } else if (prior$class[i] %in% autocor_pars) {
        if (prior$bound[i] != "<lower=-1,upper=1>") {
          autocor_warning <- TRUE
        }
      } else if (prior$class[i] %in% c("simplex", "theta")) {
        if (nchar(prior$prior[i]) && !grepl("^dirichlet\\(", prior$prior[i])) {
          stop2(
            "Currently 'dirichlet' is the only valid prior ",
            "for simplex parameters. See help(set_prior) ",
            "for more details."
          )
        }
      }
    } 
    if (nchar(lb_warning) && warn) {
      warning2(
        "It appears as if you have specified a lower bounded ", 
        "prior on a parameter that has no natural lower bound.",
        "\nIf this is really what you want, please specify ",
        "argument 'lb' of 'set_prior' appropriately.",
        "\nWarning occurred for prior \n", lb_warning
      )
    }
    if (nchar(ub_warning) && warn) {
      warning2(
        "It appears as if you have specified an upper bounded ", 
        "prior on a parameter that has no natural upper bound.",
        "\nIf this is really what you want, please specify ",
        "argument 'ub' of 'set_prior' appropriately.",
        "\nWarning occurred for prior \n", ub_warning
      )
    }
    if (autocor_warning && warn) {
      warning2(
        "Changing the boundaries of autocorrelation ", 
        "parameters is not recommended."
      )
    }
  }
  invisible(TRUE)
}

#' @export
check_prior_special.brmsterms <- function(x, prior = NULL, ...) {
  if (isTRUE(attr(prior, "checked"))) {
    return(prior) 
  }
  if (is.null(prior)) {
    prior <- empty_brmsprior()
  }
  simple_sigma <- has_sigma(x$family, x) && is.null(x$dpars$sigma)
  for (dp in names(x$dpars)) {
    allow_autoscale <- simple_sigma && identical(dp, "mu") 
    prior <- check_prior_special(
      x$dpars[[dp]], prior, allow_autoscale = allow_autoscale, ...
    )
  }
  prior
}

#' @export
check_prior_special.btnl <- function(x, prior, ...) {
  stopifnot(is.brmsprior(prior))
  for (nlp in names(x$nlpars)) {
    prior <- check_prior_special(x$nlpars[[nlp]], prior, ...)
  }
  prior
}

#' @export
check_prior_special.btl <- function(x, prior, allow_autoscale = TRUE, ...) {
  # prepare special priors such as horseshoe or lasso
  # Args:
  #   prior: an object of class brmsprior
  #   allow_autoscale: allow autoscaling using sigma?
  # Returns:
  #   a possibly amended brmsprior object with additional attributes
  px <- check_prefix(x)
  prior_special <- list()
  b_index <- which(find_rows(prior, class = "b", coef = "", ls = px))
  stopifnot(length(b_index) <= 1L)
  if (length(b_index)) {
    b_prior <- prior$prior[b_index]
    if (any(grepl("^(horseshoe|lasso)\\(", b_prior))) {
      # horseshoe prior for population-level parameters
      if (any(nzchar(prior[b_index, "bound"]))) {
        stop2("Boundaries for population-level effects are not", 
              "allowed when using the horseshoe or lasso priors.")
      }
      if (any(ulapply(x[c("me", "mo", "cs")], is.formula))) {
        stop2("Horseshoe or lasso priors are not yet allowed ",
              "in models with special population-level effects.")
      }
      b_coef_indices <- which(
        find_rows(prior, class = "b", ls = px) &
          !find_rows(prior, coef = c("", "Intercept"))
      )
      if (any(nchar(prior$prior[b_coef_indices]))) {
        stop2(
          "Defining priors for single population-level parameters",
          "is not allowed when using horseshoe or lasso priors",
          "(except for the Intercept)."
        )
      }
      if (grepl("^horseshoe\\(", b_prior)) {
        hs <- eval2(b_prior)
        prior$prior[b_index] <- ""
        hs_obj_names <- c(
          "df", "df_global", "df_slab", "scale_global", 
          "scale_slab", "par_ratio", "autoscale"
        )
        hs_att <- attributes(hs)[hs_obj_names]
        names(hs_att) <- paste0("hs_", names(hs_att))
        prior_special <- c(prior_special, hs_att)
        prior_special$hs_autoscale <- 
          isTRUE(prior_special$hs_autoscale) && allow_autoscale
      } else if (grepl("^lasso\\(", b_prior)) {
        lasso <- eval2(b_prior)
        # the parameterization via double_exponential appears to be more
        # efficient than an indirect parameterization via normal and 
        # exponential distributions; tested on 2017-06-09
        p <- usc(combine_prefix(px))
        lasso_scale <- paste0(
          "lasso_scale", p, " * lasso_inv_lambda", p
        )
        lasso_prior <- paste0(
          "double_exponential(0, ", lasso_scale, ")"
        )
        prior$prior[b_index] <- lasso_prior
        lasso_att <- attributes(lasso)
        prior_special$lasso_df <- lasso_att[["df"]]
        prior_special$lasso_scale <- lasso_att[["scale"]]
      }
    }
  }
  prefix <- combine_prefix(px, keep_mu = TRUE)
  attributes(prior)$special[[prefix]] <- prior_special
  prior
}

check_sample_prior <- function(sample_prior) {
  # validate argument 'sample_prior'
  options <- c("no", "yes", "only")
  if (!is.character(sample_prior)) {
    sample_prior <- as_one_logical(sample_prior)
    sample_prior <- if (sample_prior) "yes" else "no"
  }
  match.arg(sample_prior, options)
}

get_bound <- function(prior, class = "b", coef = "", 
                      group = "", px = list()) {
  # extract the boundaries of a parameter described by class etc.
  # Args:
  #   prior: object of class brmsprior
  #   class, coef, group, nlpar: strings of length 1
  stopifnot(length(class) == 1L)
  if (!length(coef)) coef <- ""
  if (!length(group)) group <- ""
  bound <- subset2(prior, ls = c(nlist(class, coef, group), px))$bound
  if (length(bound) > 1L) {
    stop("extracted more than one boundary at once")
  }
  bound
}

brmsprior <- function(prior = "", class = "", coef = "", group = "", 
                      resp = "", dpar = "", nlpar = "", bound = "",
                      ls = list()) {
  # helper function to create data.frames containing prior information
  if (length(ls)) {
    if (is.null(names(ls))) {
      stop("Argument 'ls' must be named.")
    }
    names <- c(
      "prior", "class", "coef", "group", 
      "resp", "dpar", "nlpar", "bound"
    )
    if (!all(names(ls) %in% names)) {
      stop("Names of 'ls' must some of ", collapse_comma(names))
    }
    for (v in names(ls)) {
      assign(v, ls[[v]])
    }
  }
  out <- data.frame(
    prior, class, coef, group, resp, dpar, nlpar, bound, 
    stringsAsFactors = FALSE
  )
  class(out) <- c("brmsprior", "data.frame")
  out
}

empty_brmsprior <- function() {
  # define a brmsprior object with zero rows
  char0 <- character(0)
  brmsprior(
    prior = char0, class = char0, coef = char0, 
    group = char0, resp = char0, dpar = char0,
    nlpar = char0, bound = char0
  )
}

prior_bounds <- function(prior) {
  # natural upper and lower bounds for priors
  # Returns:
  #   A named list with elements 'lb and 'ub'
  switch(prior,
    lognormal = list(lb = 0, ub = Inf),
    chi_square = list(lb = 0, ub = Inf),
    inv_chi_square = list(lb = 0, ub = Inf),
    scaled_inv_chi_square = list(lb = 0, ub = Inf),
    exponential = list(lb = 0, ub = Inf),
    gamma = list(lb = 0, ub = Inf),
    inv_gamma = list(lb = 0, ub = Inf),
    weibull = list(lb = 0, ub = Inf),
    frechet = list(lb = 0, ub = Inf),
    rayleigh = list(lb = 0, ub = Inf),
    pareto = list(lb = 0, ub = Inf),
    pareto_type_2 = list(lb = 0, ub = Inf),
    beta = list(lb = 0, ub = 1),
    von_mises = list(lb = -pi, ub = pi),
    list(lb = -Inf, ub = Inf)
  )
}

par_bounds <- function(par, bound = "") {
  # upper and lower bounds for parameter classes
  # Returns:
  #   A named list with elements 'lb and 'ub'
  out <- switch(par,
    sigma = list(lb = 0, ub = Inf),
    shape = list(lb = 0, ub = Inf),
    nu = list(lb = 1, ub = Inf),
    phi = list(lb = 0, ub = Inf),
    kappa = list(lb = 0, ub = Inf), 
    beta = list(lb = 0, ub = Inf),
    zi = list(lb = 0, ub = 1),
    hu = list(lb = 0, ub = 1),
    zoi = list(lb = 0, ub = 1),
    coi = list(lb = 0, ub = 1),
    bs = list(lb = 0, ub = Inf),
    ndt = list(lb = 0, ub = "min_Y"), 
    bias = list(lb = 0, ub = 1), 
    disc = list(lb = 0, ub = Inf),
    quantile = list(lb = 0, ub = 1),
    ar = list(lb = -1, ub = 1),
    ma = list(lb = -1, ub = 1),
    lagsar = list(lb = 0, ub = 1),
    errorsar = list(lb = 0, ub = 1),
    car = list(lb = 0, ub = 1),
    sdcar = list(lb = 0, ub = Inf),
    sigmaLL = list(lb = 0, ub = Inf),
    sd = list(lb = 0, ub = Inf),
    sds = list(lb = 0, ub = Inf),
    sdgp = list(lb = 0, ub = Inf),
    lscale = list(lb = 0, ub = Inf),
    list(lb = -Inf, ub = Inf)
  )
  if (isTRUE(nzchar(bound))) {
    opt_lb <- get_matches("(<|,)lower=[^,>]+", bound)
    if (isTRUE(nzchar(opt_lb))) {
      out$lb <- substr(opt_lb, 8, nchar(opt_lb))
    } 
    opt_ub <- get_matches("(<|,)upper=[^,>]+", bound)
    if (isTRUE(nzchar(opt_ub))) {
      out$ub <- substr(opt_ub, 8, nchar(opt_ub)) 
    } 
  }
  out
}

#' Checks if argument is a \code{brmsprior} object
#' 
#' @param x An \R object
#' 
#' @export
is.brmsprior <- function(x) {
  inherits(x, "brmsprior")
}

#' Print method for \code{brmsprior} objects
#' 
#' @param x An object of class \code{brmsprior}.
#' @param show_df Logical; Print priors as a single 
#'   \code{data.frame} (\code{TRUE}) or as a sequence of 
#'   sampling statements (\code{FALSE})?
#' @param ... Currently ignored.
#' 
#' @export
print.brmsprior <- function(x, show_df, ...) {
  if (missing(show_df)) {
    show_df <- nrow(x) > 1L
  }
  if (show_df) {
    NextMethod()
  } else {
    cat(collapse(.print_prior(x), "\n"))
  }
  invisible(x)
}

.print_prior <- function(x) {
  # prepare text for print.brmsprior
  group <-  usc(x$group)
  resp <- usc(x$resp)
  dpar <- usc(x$dpar)
  nlpar <- usc(x$nlpar)
  coef <- usc(x$coef)
  if (any(nzchar(c(resp, dpar, nlpar, coef)))) {
    group <- usc(group, "suffix")
  }
  bound <- ifelse(nzchar(x$bound), paste0(x$bound, " "), "")
  tilde <- ifelse(nzchar(x$class) | nzchar(group) | nzchar(coef), " ~ ", "")
  prior <- ifelse(nzchar(x$prior), x$prior, "(no prior)")
  paste0(bound, x$class, group, resp, dpar, nlpar, coef, tilde, prior)
}

#' @export
c.brmsprior <- function(x, ...) {
  # combine multiple brmsprior objects into one brmsprior
  if (all(sapply(list(...), is.brmsprior))) {
    out <- do.call(rbind, list(x, ...)) 
  } else {
    out <- c(as.data.frame(x), ...)
  }
  out
}

#' @export
"+.brmsprior" <- function(e1, e2) {
  c(e1, e2)
}

dirichlet <- function(...) {
  # dirichlet prior of simplex parameters
  out <- as.numeric(c(...))
  if (anyNA(out) || any(out <= 0)) {
    stop2("The dirichlet prior expects positive values.")
  }
  out
}

#' Set up a horseshoe prior in \pkg{brms}
#' 
#' Function used to set up a horseshoe prior for population-level effects 
#' in \pkg{brms}. The function does not evaluate its arguments --
#' it exists purely to help set up the model.
#' 
#' @param df Degrees of freedom of student-t prior of the 
#'   local shrinkage parameters. Defaults to \code{1}.
#' @param scale_global Scale of the student-t prior of the global shrinkage 
#'   parameter. Defaults to \code{1}. 
#'   In linear models, \code{scale_global} will internally be 
#'   multiplied by the residual standard deviation parameter \code{sigma}.
#' @param df_global Degrees of freedom of student-t prior of the 
#'   global shrinkage parameter. Defaults to \code{1}.
#' @param scale_slab Scale of the student-t prior of the regularization
#'   parameter. Defaults to \code{2}. 
#' @param df_slab Degrees of freedom of the student-t prior of 
#'   the regularization parameter. Defaults to \code{4}. 
#' @param par_ratio Ratio of the expected number of non-zero coefficients 
#'   to the expected number of zero coefficients. If specified,
#'   \code{scale_global} is ignored and internally computed as
#'   \code{par_ratio / sqrt(N)}, where \code{N} is the total number 
#'   of observations in the data.
#' @param autoscale Logical; indicating whether the horseshoe
#'   prior should be scaled using the residual standard deviation
#'   \code{sigma} if possible and sensible (defaults to \code{TRUE}).
#'   Autoscaling is not applied for distributional parameters or 
#'   when the model does not contain the parameter \code{sigma}.
#'   
#' @return A character string obtained by \code{match.call()} with
#'   additional arguments.
#'   
#' @details 
#'   The horseshoe prior is a special shrinkage prior initially proposed by
#'   Carvalho et al. (2009).
#'   It is symmetric around zero with fat tails and an infinitely large spike
#'   at zero. This makes it ideal for sparse models that have 
#'   many regression coefficients, although only a minority of them is non-zero. 
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
#'   In linear models, \code{scale_global} will internally be multiplied by the 
#'   residual standard deviation parameter \code{sigma}. See Piironen and 
#'   Vehtari (2016) for recommendations how to properly set the global scale.
#'   The degrees of freedom of the global shrinkage prior may also be 
#'   adjusted via argument \code{df_global}. 
#'   Piironen and Vehtari (2017) recommend to specifying the ratio of the 
#'   expected number of non-zero coefficients to the expected number of zero 
#'   coefficients \code{par_ratio} rather than \code{scale_global} directly.
#'   As proposed by Piironen and Vehtari (2017), an additional regularization 
#'   is applied that only affects non-zero coefficients. The amount of
#'   regularization can be controlled via \code{scale_slab} and \code{df_slab}.
#'   To make sure that shrinkage can equally affect all coefficients, 
#'   predictors should be one the same scale. 
#'   Generally, models with horseshoe priors a more likely than other models
#'   to have divergent transitions so that increasing \code{adapt_delta} 
#'   from \code{0.8} to values closer to \code{1} will often be necessary.
#'   See the documentation of \code{\link[brms:brm]{brm}} for instructions
#'   on how to increase \code{adapt_delta}. 
#'   
#' @references 
#' Carvalho, C. M., Polson, N. G., & Scott, J. G. (2009). 
#'   Handling sparsity via the horseshoe. 
#'   In International Conference on Artificial Intelligence and Statistics (pp. 73-80).
#'    
#' Piironen J. & Vehtari A. (2016). On the Hyperprior Choice for the Global 
#'    Shrinkage Parameter in the Horseshoe Prior. 
#'    \url{https://arxiv.org/pdf/1610.05559v1.pdf}
#'    
#' Piironen, J., and Vehtari, A. (2017). Sparsity information and regularization
#'    in the horseshoe and other shrinkage priors. 
#'    \url{https://arxiv.org/abs/1707.01694}    
#'   
#' @seealso \code{\link[brms:set_prior]{set_prior}}
#'   
#' @examples 
#' set_prior(horseshoe(df = 3, par_ratio = 0.1))
#' 
#' @export
horseshoe <- function(df = 1, scale_global = 1, df_global = 1, 
                      scale_slab = 2, df_slab = 4, par_ratio = NULL,
                      autoscale = TRUE) {
  out <- deparse(match.call(), width.cutoff = 500L)
  df <- as.numeric(df)
  df_global <- as.numeric(df_global)
  df_slab <- as.numeric(df_slab)
  scale_global <- as.numeric(scale_global)
  scale_slab <- as.numeric(scale_slab)
  if (!isTRUE(df > 0)) {
    stop2("Invalid horseshoe prior: Degrees of freedom of ", 
          "the local priors must be a single positive number.")
  }
  if (!isTRUE(df_global > 0)) {
    stop2("Invalid horseshoe prior: Degrees of freedom of ", 
          "the global prior must be a single positive number.")
  }
  if (!isTRUE(scale_global > 0)) {
    stop2("Invalid horseshoe prior: Scale of the global ", 
          "prior must be a single positive number.")
  }
  if (!isTRUE(df_slab > 0)) {
    stop2("Invalid horseshoe prior: Degrees of freedom of ", 
          "the slab part must be a single positive number.")
  }
  if (!isTRUE(scale_slab > 0)) {
    stop2("Invalid horseshoe prior: Scale of the slab ", 
          "part must be a single positive number.")
  }
  if (!is.null(par_ratio)) {
    par_ratio <- as.numeric(par_ratio)
    if (!isTRUE(par_ratio > 0 && par_ratio <= 1)) {
      stop2("Argument 'par_ratio' must be within [0, 1].")
    }
  }
  autoscale <- as_one_logical(autoscale)
  att <- nlist(
    df, df_global, df_slab, scale_global, 
    scale_slab, par_ratio, autoscale
  )
  attributes(out)[names(att)] <- att
  out
}

#' Set up a lasso prior in \pkg{brms}
#' 
#' Function used to set up a lasso prior for population-level effects 
#' in \pkg{brms}. The function does not evaluate its arguments --
#' it exists purely to help set up the model.
#' 
#' @param df Degrees of freedom of the chi-sqaure prior of the inverse tuning
#'   parameter. Defaults to \code{1}.
#' @param scale Scale of the lasso prior. Defaults to \code{1}.
#'   
#' @return A character string obtained by \code{match.call()} with
#'   additional arguments.
#'   
#' @details  
#'   The lasso prior is the Bayesian equivalent to the LASSO method for performing
#'   variable selection (Park & Casella, 2008).
#'   With this prior, independent Laplace (i.e. double exponential) priors 
#'   are placed on the population-level effects. 
#'   The scale of the Laplace priors depends on a tuning parameter
#'   that controls the amount of shrinkage. In \pkg{brms}, the inverse
#'   of the tuning parameter is used so that smaller values imply
#'   more shrinkage. The inverse tuning parameter has a chi-square distribution
#'   and with degrees of freedom controlled via argument \code{df}
#'   of function \code{lasso} (defaults to \code{1}). For instance,
#'   one can specify a lasso prior using \code{set_prior("lasso(1)")}.
#'   To make sure that shrinkage can equally affect all coefficients, 
#'   predictors should be one the same scale.
#'   If you do not want to standarized all variables,
#'   you can adjust the general scale of the lasso prior via argument
#'   \code{scale}, for instance, \code{lasso(1, scale = 10)}.
#' 
#' @references
#' Park, T., & Casella, G. (2008). The Bayesian Lasso. Journal of the American 
#'    Statistical Association, 103(482), 681-686.
#'    
#' @seealso \code{\link[brms:set_prior]{set_prior}}
#'   
#' @examples 
#' set_prior(lasso(df = 1, scale = 10))
#' 
#' @export
lasso <- function(df = 1, scale = 1) {
  out <- deparse(match.call(), width.cutoff = 500L)
  df <- as.numeric(df)
  scale <- as.numeric(scale)
  if (!isTRUE(df > 0)) {
    stop2("Invalid lasso prior: Degrees of freedom of the shrinkage ", 
          "parameter prior must be a single positive number.")
  }
  if (!isTRUE(scale > 0)) {
    stop2("Invalid lasso prior: Scale of the Laplace ", 
          "priors must be a single positive number.")
  }
  att <- nlist(df, scale)
  attributes(out)[names(att)] <- att
  out
}
