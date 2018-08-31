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
#'   Only used in multivariate models.
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
#'   argument of \code{\link{brm}}.
#' 
#' @details 
#'   \code{set_prior} is used to define prior distributions for parameters 
#'   in \pkg{brms} models. The functions \code{prior}, \code{prior_}, and
#'   \code{prior_string} are aliases of \code{set_prior} each allowing
#'   for a different kind of argument specification. 
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
#'   To combine multiple priors, use \code{c(...)} or the \code{+} operator 
#'   (see 'Examples'). \pkg{brms} does not check if the priors are written 
#'   in correct \pkg{Stan} language. Instead, \pkg{Stan} will check their 
#'   syntactical correctness when the model is parsed to \code{C++} and 
#'   returns an error if they are not. 
#'   This, however, does not imply that priors are always meaningful if they are 
#'   accepted by \pkg{Stan}. Although \pkg{brms} trys to find common problems 
#'   (e.g., setting bounded priors on unbounded parameters), there is no guarantee 
#'   that the defined priors are reasonable for the model.
#'   Below, we list the types of parameters in \pkg{brms} models, 
#'   for which the user can specify prior distributions.
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
#'   (discussed in the 'Details' section of \code{\link{brmsformula}}), 
#'   general priors on class \code{"b"} will \emph{not} affect 
#'   the intercept. Instead, the intercept has its own parameter class 
#'   named \code{"Intercept"} and priors can thus be 
#'   specified via \code{set_prior("<prior>", class = "Intercept")}.
#'   Setting a prior on the intercept will not break vectorization
#'   of the other population-level effects.
#'   Note that technically, this prior is set on an intercept that
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
#'   is the horseshoe prior. See \code{\link{horseshoe}}
#'   for details. Another shrinkage prior is the so-called lasso prior.
#'   See \code{\link{lasso}} for details.
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
#'   When using a prior that is defined on the positive reals only 
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
#'   These parameters are restricted to be non-negative and, by default, 
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
#'   4. Splines
#'   
#'   Splines are implemented in \pkg{brms} using the 'random effects' 
#'   formulation as explained in \code{\link[mgcv:gamm]{gamm}}). 
#'   Thus, each spline has its corresponding standard deviations 
#'   modeling the variability within this term. In \pkg{brms}, this 
#'   parameter class is called \code{sds} and priors can
#'   be specified via \code{set_prior("<prior>", class = "sds", 
#'   coef = "<term label>")}. The default prior is the same as
#'   for standard deviations of group-level effects.
#'   
#'   5. Gaussian processes
#'   
#'   Gaussian processes as currently implemented in \pkg{brms} have
#'   two parameters, the standard deviation parameter \code{sdgp}, 
#'   and characteristic length-scale parameter \code{lscale} 
#'   (see \code{\link{gp}} for more details). The default prior 
#'   of \code{sdgp} is the same as for standard deviations of 
#'   group-level effects. The default prior of \code{lscale}
#'   is an informative inverse-gamma prior specifically tuned 
#'   to the covariates of the Gaussian process (for more details see
#'   \url{https://betanalpha.github.io/assets/case_studies/gp_part3/part3.html}).
#'   This tuned prior may be overly informative in some cases, so please 
#'   consider other priors as well to make sure inference is
#'   robust to the prior specification. If tuning fails, a half-normal prior 
#'   is used instead.
#'   
#'   6. Autocorrelation parameters
#'   
#'   The autocorrelation parameters currently implemented are named 
#'   \code{ar} (autoregression), \code{ma} (moving average),
#'   \code{arr} (autoregression of the response), \code{car} 
#'   (spatial conditional autoregression), as well as \code{lagsar} 
#'   and \code{errorsar} (Spatial simultaneous autoregression).
#'   
#'   Priors can be defined by \code{set_prior("<prior>", class = "ar")} 
#'   for \code{ar} and similar for other autocorrelation parameters.
#'   By default, \code{ar} and \code{ma} are bounded between \code{-1} 
#'   and \code{1}, \code{car}, \code{lagsar}, and \code{errorsar} are 
#'   bounded between \code{0}, and \code{1}, and \code{arr} is unbounded 
#'   (you may change this by using the arguments \code{lb} and \code{ub}). 
#'   The default prior is flat over the definition area.
#'   
#'   7. Distance parameters of monotonic effects
#'   
#'   As explained in the details section of \code{\link{brm}},
#'   monotonic effects make use of a special parameter vector to
#'   estimate the 'normalized distances' between consecutive predictor 
#'   categories. This is realized in \pkg{Stan} using the \code{simplex}
#'   parameter type. This class is named \code{"simo"} (short for 
#'   simplex monotonic) in \pkg{brms}. 
#'   The only valid prior for simplex parameters is the
#'   dirichlet prior, which accepts a vector of length \code{K - 1}
#'   (K = number of predictor categories) as input defining the
#'   'concentration' of the distribution. Explaining the dirichlet prior 
#'   is beyond the scope of this documentation, but we want to describe
#'   how to define this prior syntactically correct.
#'   If a predictor \code{x} with \code{K} categories is modeled as monotonic, 
#'   we can define a prior on its corresponding simplex via \cr
#'   \code{prior(dirichlet(<vector>), class = simo, coef = mox1)}.
#'   The \code{1} in the end of \code{coef} indicates that this is the first
#'   simplex in this term. If interactions between multiple monotonic
#'   variables are modeled, multiple simplexes per term are required.  
#'   For \code{<vector>}, we can put in any \code{R} expression
#'   defining a vector of length \code{K - 1}. The default is a uniform 
#'   prior (i.e. \code{<vector> = rep(1, K-1)}) over all simplexes
#'   of the respective dimension.   
#'   
#'   8. Parameters for specific families 
#'   
#'   Some families need additional parameters to be estimated. 
#'   Families \code{gaussian}, \code{student}, \code{skew_normal},
#'   \code{lognormal}, and \code{gen_extreme_value} need the parameter 
#'   \code{sigma} to account for the residual standard deviation.
#'   By default, \code{sigma} has a half student-t prior that scales 
#'   in the same way as the group-level standard deviations.
#'   Further, family \code{student} needs the parameter 
#'   \code{nu} representing the degrees of freedom of students-t distribution. 
#'   By default, \code{nu} has prior \code{"gamma(2, 0.1)"}
#'   and a fixed lower bound of \code{1}.
#'   Families \code{gamma}, \code{weibull}, \code{inverse.gaussian}, and
#'   \code{negbinomial} need a \code{shape} parameter that has a 
#'   \code{"gamma(0.01, 0.01)"} prior by default. 
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
#'   use function \code{\link{get_prior}}.
#'
#' @seealso \code{\link{get_prior}}
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
                      lb = NA, ub = NA, check = TRUE) {
  input <- nlist(prior, class, coef, group, resp, dpar, nlpar, lb, ub, check)
  input <- try(as.data.frame(input), silent = TRUE)
  if (is(input, "try-error")) {
    stop2("Processing arguments of 'set_prior' has failed:\n", input)
  }
  out <- vector("list", nrow(input))
  for (i in seq_along(out)) {
    out[[i]] <- do.call(.set_prior, input[i, ])
  }
  Reduce("+", out)
}
  
.set_prior <- function(prior, class, coef, group, resp, 
                       dpar, nlpar, lb, ub, check) {
  # validate arguments passed to 'set_prior'
  prior <- as_one_character(prior)
  class <- as_one_character(class)
  group <- as_one_character(group)
  coef <- as_one_character(coef)
  resp <- as_one_character(resp)
  dpar <- as_one_character(dpar)
  nlpar <- as_one_character(nlpar)
  lb <- as_one_character(lb, allow_na = TRUE)
  ub <- as_one_character(ub, allow_na = TRUE)
  check <- as_one_logical(check)
  # validate boundaries
  bound <- ""
  is_arma <- class %in% c("ar", "ma")
  if (!is.na(lb) || !is.na(ub) || is_arma) {
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
      lb <- ifelse(!is.na(lb), lb, -1)
      ub <- ifelse(!is.na(ub), ub, 1)
    }
    # don't put spaces in boundary declarations
    lb <- if (!is.na(lb)) paste0("lower=", lb)
    ub <- if (!is.na(ub)) paste0("upper=", ub)
    if (!is.null(lb) || !is.null(ub)) {
      bound <- paste0("<", paste(c(lb, ub), collapse = ","), ">")
    }
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
#' @param internal A flag indicating if the names of additional internal parameters 
#'   should be displayed. Setting priors on these parameters is not recommended
#' 
#' @return A data.frame with columns \code{prior}, \code{class}, \code{coef}, and \code{group}
#'   and several rows, each providing information on a parameter (or parameter class) on which
#'   priors can be specified. The prior column is empty except for internal default priors.
#'   
#' @seealso \code{\link{set_prior}}
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
                      autocor = NULL, internal = FALSE) {
  # note that default priors are stored in this function
  if (is.brmsfit(formula)) {
    stop2("Use 'prior_summary' to extract priors from 'brmsfit' objects.")
  }
  formula <- validate_formula(
    formula, data = data, family = family, autocor = autocor
  )
  bterms <- parse_bf(formula)
  data <- update_data(data, bterms = bterms)
  ranef <- tidy_ranef(bterms, data)
  meef <- tidy_meef(bterms, data)
  # initialize output
  prior <- empty_brmsprior()
  # priors for distributional parameters
  prior <- prior + prior_effects(
    bterms, data = data, internal = internal
  )
  # priors of group-level parameters
  def_scale_prior <- def_scale_prior(bterms, data)
  prior <- prior + prior_re(
    ranef, def_scale_prior = def_scale_prior,
    internal = internal
  )
  # priors for noise-free variables
  prior <- prior + prior_Xme(meef, internal = internal)
  # do not remove unique(.)
  to_order <- with(prior, order(resp, dpar, nlpar, class, group, coef))
  prior <- unique(prior[to_order, , drop = FALSE])
  rownames(prior) <- NULL
  structure(prior, class = c("brmsprior", "data.frame"))
}

prior_effects <- function(x, ...) {
  # generate priors various kind of effects 
  UseMethod("prior_effects")
}

#' @export
prior_effects.default <- function(x, ...) {
  empty_brmsprior()
}

prior_effects.mvbrmsterms <- function(x, internal = FALSE, ...) {
  prior <- empty_brmsprior()
  for (i in seq_along(x$terms)) {
    prior <- prior + prior_effects(x$terms[[i]], ...) 
  }
  # add "global" priors for population-level effects
  # in 1.8.0 as users keep asking about this
  for (cl in c("b", "Intercept")) {
    if (any(with(prior, class == cl & coef == ""))) {
      prior <- prior + brmsprior(class = cl)
    }
  }
  if (x$rescor) {
    if (internal) {
      prior <- prior +
        brmsprior(class = "Lrescor", prior = "lkj_corr_cholesky(1)")
    } else {
      prior <- prior + brmsprior(class = "rescor", prior = "lkj(1)")
    }
    if (family_names(x)[1] %in% "student") {
      prior <- prior + brmsprior(class = "nu", prior = "gamma(2, 0.1)")
    }
  }
  prior
}

prior_effects.brmsterms <- function(x, data, ...) {
  def_scale_prior <- def_scale_prior(x, data)
  valid_dpars <- valid_dpars(x$family, bterms = x)
  prior <- empty_brmsprior()
  for (dp in valid_dpars) {
    def_dprior <- def_dprior(x, dp, data = data)
    if (!is.null(x$dpars[[dp]])) {
      # parameter is predicted
      dp_prior <- prior_effects(
        x$dpars[[dp]], data = data,
        def_scale_prior = def_scale_prior,
        def_dprior = def_dprior
      )
    } else if (!is.null(x$fdpars[[dp]])){
      # parameter is fixed
      dp_prior <- empty_brmsprior()
    } else {
      # parameter is estimated
      dp_prior <- brmsprior(def_dprior, class = dp, resp = x$resp)
    }
    prior <- prior + dp_prior
  }
  for (nlp in names(x$nlpars)) {
    nlp_prior <- prior_effects(
      x$nlpars[[nlp]], data = data,
      def_scale_prior = def_scale_prior,
      def_dprior = def_dprior, 
      spec_intercept = FALSE
    )
    prior <- prior + nlp_prior
  }
  # global population-level priors for categorical models
  if (is_categorical(x$family)) {
    for (cl in c("b", "Intercept")) {
      if (any(find_rows(prior, class = cl, coef = "", resp = x$resp))) {
        prior <- prior + brmsprior(class = cl, resp  = x$resp)
      }
    }
  }
  # prior for the delta parameter for equidistant thresholds
  if (is_ordinal(x$family) && is_equal(x$family$threshold, "equidistant")) {
    bound <- ifelse(x$family$family == "cumulative", "<lower=0>", "")
    prior <- prior + brmsprior(class = "delta", bound = bound, resp = x$resp)
  }
  # priors for mixture models
  dp_classes <- dpar_class(names(c(x$dpars, x$fdpars)))
  if (is.mixfamily(x$family) && !any(dp_classes == "theta")) {
    prior <- prior + brmsprior(class = "theta", resp = x$resp)
  }
  # priors for noise-free response variables
  sdy <- get_sdy(x, data)
  if (!is.null(sdy)) {
    prior <- prior + 
      brmsprior(class = "meanme", resp = x$resp) +
      brmsprior(class = "sdme", resp = x$resp)
  }
  # priors for autocorrelation parameters
  prior <- prior + prior_autocor(x, def_scale_prior = def_scale_prior)
  prior
}

#' @export
prior_effects.btl <- function(x, data, spec_intercept = TRUE,
                              def_scale_prior = "", def_dprior = "", 
                              ...) {
  # collect default priors for various kinds of effects
  # Args:
  #   spec_intercept: special parameter class for the Intercept?
  #   def_scale_prior: default prior for SD parameters
  # Return:
  #   An object of class brmsprior
  prior_fe(x, data, def_dprior = def_dprior, spec_intercept = spec_intercept) +
    prior_sp(x, data) +
    prior_cs(x, data) +
    prior_sm(x, data, def_scale_prior = def_scale_prior) + 
    prior_gp(x, data, def_scale_prior = def_scale_prior)
}

#' @export
prior_effects.btnl <- function(x, data, ...) {
  empty_brmsprior()
}

prior_fe <- function(bterms, data, spec_intercept = TRUE, def_dprior = "") {
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
      brmsprior(def_dprior, class = "Intercept", coef = "", ls = px)
    fixef <- setdiff(fixef, "Intercept")
  }
  if (length(fixef)) {
    prior <- prior + 
      brmsprior(class = "b", coef = c("", fixef), ls = px)
  }
  prior
}

prior_sp <- function(bterms, data) {
  # priors for special effects parameters
  # Returns:
  #   an object of class brmsprior
  prior <- empty_brmsprior()
  spef <- tidy_spef(bterms, data)
  if (nrow(spef)) {
    px <- check_prefix(bterms)
    prior <- prior + brmsprior(
      class = "b", coef = c("", spef$coef), ls = px
    )
    simo_coef <- get_simo_labels(spef)
    if (length(simo_coef)) {
      prior <- prior + brmsprior(
        class = "simo", coef = simo_coef, ls = px
      ) 
    }
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

prior_Xme <- function(meef, internal = FALSE) {
  # default priors for hyper-parameters of noise-free variables
  # Returns:
  #   an object of class brmsprior
  stopifnot(is.meef_frame(meef))
  prior <- empty_brmsprior()
  if (nrow(meef)) {
    prior <- prior + 
      brmsprior(class = "meanme", coef = c("", meef$coef)) +
      brmsprior(class = "sdme", coef = c("", meef$coef))
    # priors for correlation parameters
    groups <- unique(meef$grname)
    for (i in seq_along(groups)) {
      g <- groups[i]
      K <- which(meef$grname %in% g)
      if (meef$cor[K[1]] && length(K) > 1L) {
        if (internal) {
          prior <- prior + brmsprior("lkj_corr_cholesky(1)", class = "Lme")
          if (nzchar(g)) {
            prior <- prior + brmsprior(class = "Lme", group = g)
          }
        } else {
          prior <- prior + brmsprior("lkj(1)", class = "corme")
          if (nzchar(g)) {
            prior <- prior + brmsprior(class = "corme", group = g)
          }
        }
      }
    }
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
  gpef <- tidy_gpef(bterms, data)
  if (nrow(gpef)) {
    px <- check_prefix(bterms)
    lscale_prior <- def_lscale_prior(bterms, data)
    # GPs of each 'by' level get their own 'lscale' prior
    all_gpterms <- ulapply(seq_rows(gpef), 
      function(i) paste0(gpef$term[i], gpef$bylevels[[i]])
    )
    prior <- prior +
      brmsprior(class = "sdgp", prior = def_scale_prior, ls = px) +
      brmsprior(class = "sdgp", coef = gpef$term, ls = px) +
      brmsprior(class = "lscale", prior = "normal(0, 0.5)", ls = px) +
      brmsprior(class = "lscale", prior = lscale_prior, 
                coef = all_gpterms, ls = px)
  }
  prior
}

def_lscale_prior <- function(bterms, data, plb = 0.01, pub = 0.01) {
  # default priors for length-scale parameters of GPs
  # see https://betanalpha.github.io/assets/case_studies/gp_part3/part3.html
  # Args:
  #   plb: prior probability of being lower than minimum length-scale
  #   pub: prior probability of being higher than maximum length-scale
  .opt_fun <- function(x, lb, ub) {
    # optimize parameters on the log-scale to make them positive only
    x <- exp(x)
    y1 <- pinvgamma(lb, x[1], x[2], log.p = TRUE)
    y2 <- pinvgamma(ub, x[1], x[2], lower.tail = FALSE, log.p = TRUE)
    c(y1 - log(plb), y2 - log(pub))
  }
  gp_dat <- data_gp(bterms, data)
  Xgp_names <- names(gp_dat)[grepl("Xgp_", names(gp_dat))]
  out <- rep("", length(Xgp_names))
  for (i in seq_along(out)) {
    dq <- diff_quad(gp_dat[[Xgp_names[i]]])
    lb <- sqrt(min(dq[dq > 0]))
    ub <- sqrt(max(dq))
    opt_res <- nleqslv::nleqslv(
      c(0, 0), .opt_fun, lb = lb, ub = ub,
      control = list(allowSingular = TRUE)
    )
    if (opt_res$termcd %in% 1:2) {
      # use the inverse-gamma prior only in case of convergence
      pars <- exp(opt_res$x)
      out[i] <- paste0("inv_gamma(", sargs(round(pars, 6)), ")") 
    }
  }
  out
}

prior_re <- function(ranef, def_scale_prior, internal = FALSE) {
  # priors for random effects parameters
  # Args:
  #   ranef: a list returned by tidy_ranef
  #   def_scale_prior: a character string defining the default
  #                    prior for random effects SDs
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
  if (length(def_scale_prior) > 1L) {
    def_scale_prior <- def_scale_prior[px$resp] 
  }
  global_sd_prior <- brmsprior(
    class = "sd", prior = def_scale_prior, ls = px
  )
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
  tranef <- get_dist_groups(ranef, "student")
  if (isTRUE(nrow(tranef) > 0L)) {
    prior <- prior + 
      brmsprior("gamma(2, 0.1)", class = "df", group = tranef$group)
  }
  prior
}

prior_sm <- function(bterms, data, def_scale_prior) {
  # priors for smooth terms
  # Args:
  #   def_scale_prior: a character string defining 
  #     the default prior for smooth SDs
  prior <- empty_brmsprior()
  smterms <- all_terms(bterms[["sm"]])
  if (length(smterms)) {
    px <- check_prefix(bterms)
    prior_strings <- c(def_scale_prior, rep("", length(smterms)))
    prior <- prior + brmsprior(
      class = "sds", coef = c("", smterms), 
      prior = prior_strings, ls = px
    )
  }
  prior
}

prior_autocor <- function(bterms, def_scale_prior) {
  # priors for autocor parameters
  stopifnot(is.brmsterms(bterms))
  autocor <- bterms$autocor
  resp <- bterms$resp
  cbound <- "<lower=-1,upper=1>"
  prior <- empty_brmsprior()
  if (get_ar(autocor)) {
    prior <- prior + brmsprior(class = "ar", resp = resp, bound = cbound)
  }
  if (get_ma(autocor)) {
    prior <- prior + brmsprior(class = "ma", resp = resp, bound = cbound)
  }
  if (get_arr(autocor)) {
    prior <- prior + brmsprior(class = "arr", resp = resp)
  }
  if (is.cor_sar(autocor)) {
    if (identical(autocor$type, "lag")) {
      prior <- prior + brmsprior(class = "lagsar", resp = resp)
    }
    if (identical(autocor$type, "error")) {
      prior <- prior + brmsprior(class = "errorsar", resp = resp)
    }
  }
  if (is.cor_car(autocor)) {
    prior <- prior +  
      brmsprior(def_scale_prior, class = "sdcar", resp = resp)
    if (identical(autocor$type, "escar")) {
      prior <- prior + brmsprior(class = "car", resp = resp)
    }
  }
  if (is.cor_bsts(autocor)) {
    prior <- prior +
      brmsprior(class = "sigmaLL", prior = def_scale_prior, resp = resp)
  }
  prior
}

def_dprior <- function(x, dpar, data = NULL) {
  # default priors for distributional parameters
  stopifnot(is.brmsterms(x))
  dpar <- as_one_character(dpar)
  dpar_class <- dpar_class(dpar)
  link <- x$dpars[[dpar]]$family$link
  if (is.null(link)) {
    link <- "identity"
  }
  # ensures reasonable scaling in def_scale_prior
  x$family$link <- link
  if (link == "identity") {
    # dpar is estimated or predicted on the linear scale
    out <- switch(dpar_class, "",
      mu = def_scale_prior(x, data, center = FALSE),
      sigma = def_scale_prior(x, data), 
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
      disc = "lognormal(0, 1)",
      theta = "logistic(0, 1)"
    )
  } else {
    # except for 'mu' all parameters only support one link other than identity
    out <- switch(dpar_class, "",
      mu = def_scale_prior(x, data, center = FALSE),
      sigma = def_scale_prior(x, data),
      shape = "student_t(3, 0, 10)",
      nu = "normal(2.7, 0.8)", 
      phi = "student_t(3, 0, 10)",
      kappa = "normal(5.0, 0.8)", 
      beta = "normal(1.7, 1.3)", 
      zi = "logistic(0, 1)", 
      hu = "logistic(0, 1)", 
      zoi = "logistic(0, 1)",
      coi = "logistic(0, 1)",
      bs = "normal(-0.6, 1.3)", 
      bias = "logistic(0, 1)", 
      quantile = "logistic(0, 1)",
      xi = "normal(0, 4)",
      alpha = "normal(0, 4)",
      disc = "normal(0, 1)"
    )
  }
  out
}

def_scale_prior <- function(x, data, ...) {
  # ensure that SD parameters have a weakly informative prior by default
  UseMethod("def_scale_prior")
}

#' @export
def_scale_prior.mvbrmsterms <- function(x, data, ...) {
  out <- ulapply(x$terms, def_scale_prior, data = data, ...)
  names(out) <- x$responses
  out
}

#' @export
def_scale_prior.brmsterms <- function(x, data, center = TRUE, ...) {
  # Args:
  #   center: Should the prior be centererd around zero?
  #     If FALSE, the prior location is computed based on Y.
  Y <- unname(model.response(model.frame(x$respform, data)))
  prior_location <- 0
  prior_scale <- 10
  link <- x$family$link
  if (has_logscale(x$family)) {
    link <- "log"
  }
  tlinks <- c("identity", "log", "inverse", "sqrt", "1/mu^2")
  if (link %in% tlinks && !is_like_factor(Y)) {
    if (link %in% c("log", "inverse", "1/mu^2")) {
      # avoid Inf in link(Y)
      Y <- ifelse(Y == 0, Y + 0.1, Y) 
    }
    sgst_scale <- SW(round(mad(link(Y, link = link))))
    if (is.finite(sgst_scale)) {
      prior_scale <- max(prior_scale, sgst_scale)
    } 
    if (!center) {
      sgst_location <- SW(round(median(link(Y, link = link))))
      if (is.finite(sgst_location)) {
        prior_location <- sgst_location
      }
    }
  }
  paste0("student_t(", sargs("3", prior_location, prior_scale), ")")
}

check_prior <- function(prior, formula, data = NULL, 
                        sample_prior = c("no", "yes", "only"),
                        check_rows = NULL, warn = FALSE) {
  # check prior input and amend it if needed
  # Args:
  #   same as the respective parameters in brm
  #   check_rows: if not NULL, check only the rows given in check_rows
  #   warn: passed to check_prior_content
  # Returns:
  #   a data.frame of prior specifications to be used in stan_prior
  sample_prior <- check_sample_prior(sample_prior)
  if (isTRUE(attr(prior, "checked"))) {
    # prior has already been checked; no need to do it twice
    # attributes may still need to be updated
    attr(prior, "sample_prior") <- sample_prior
    return(prior)
  }
  bterms <- parse_bf(formula)
  all_priors <- get_prior(formula = formula, data = data, internal = TRUE)
  if (is.null(prior)) {
    prior <- all_priors  
  } else if (!is.brmsprior(prior)) {
    stop2("Argument 'prior' must be a 'brmsprior' object.")
  }
  # temporarily exclude priors that should not be checked
  no_checks <- !nzchar(prior$class)
  prior_no_checks <- prior[no_checks, ]
  prior <- prior[!no_checks, ]
  # check for duplicated priors
  prior$class <- rename(
    prior$class, c("^cor$", "^rescor$", "^corme$"), 
    c("L", "Lrescor", "Lme"), fixed = FALSE
  )
  rcols <- rcols_prior()
  duplicated_input <- duplicated(prior[, rcols])
  if (any(duplicated_input)) {
    stop2("Duplicated prior specifications are not allowed.")
  }
  # check if parameters in prior are valid
  if (nrow(prior)) {
    valid <- which(duplicated(rbind(all_priors[, rcols], prior[, rcols])))
    invalid <- which(!seq_rows(prior) %in% (valid - nrow(all_priors)))
    if (length(invalid)) {
      msg_priors <- .print_prior(prior[invalid, ])
      stop2(
        "The following priors do not correspond ", 
        "to any model parameter: \n",
        collapse(.print_prior(prior[invalid, ]), "\n"),
        "Function 'get_prior' might be helpful to you."
      )
    }
  }
  prior$prior <- sub("^(lkj|lkj_corr)\\(", "lkj_corr_cholesky(", prior$prior)
  check_prior_content(prior, warn = warn)
  # merge user-specified priors with default priors
  prior$new <- rep(TRUE, nrow(prior))
  all_priors$new <- rep(FALSE, nrow(all_priors))
  prior <- prior + all_priors
  prior <- prior[!duplicated(prior[, rcols]), ]
  prior <- check_prior_special(prior, bterms = bterms, data = data)
  prior <- prior[with(prior, order(class, group, resp, dpar, nlpar, coef)), ]
  prior <- prior + prior_no_checks
  rownames(prior) <- NULL
  attr(prior, "sample_prior") <- sample_prior
  attr(prior, "checked") <- TRUE
  prior
}

check_prior_content <- function(prior, warn = TRUE) {
  # try to check if prior distributions are reasonable
  # Args:
  #  prior: A brmsprior object
  #  warn: logical; print boundary warnings?
  if (!is.brmsprior(prior)) {
    return(invisible(TRUE))
  }
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
    nb_pars <- c("b", "alpha", "xi")
    lb_pars <- c(
      "sigma", "shape", "nu", "phi", "kappa", "beta", "bs", 
      "disc", "sdcar", "sigmaLL", "sd", "sds", "sdgp", "lscale" 
    )
    cor_pars <- c("cor", "rescor", "corme", "L", "Lrescor", "Lme")
    autocor_pars <- c("ar", "ma")
    lb_warning <- ub_warning <- ""
    autocor_warning <- FALSE
    for (i in seq_rows(prior)) {
      msg_prior <- .print_prior(prior[i, , drop = FALSE])
      has_lb_prior <- grepl(lb_priors_reg, prior$prior[i])
      has_ulb_prior <- grepl(ulb_priors_reg, prior$prior[i])
      # priors with nchar(coef) inherit their boundaries 
      j <- which(with(prior, 
        class == class[i] & group == group[i] & 
        nlpar == nlpar[i] & !nzchar(coef)
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
        if (nzchar(prior$prior[i]) && !grepl("^lkj", prior$prior[i])) {
          stop2(
            "The only supported prior for correlation matrices is ", 
            "the 'lkj' prior. See help(set_prior) for more details."
          )
        }
      } else if (prior$class[i] %in% autocor_pars) {
        if (prior$bound[i] != "<lower=-1,upper=1>") {
          autocor_warning <- TRUE
        }
      } else if (prior$class[i] %in% c("simo", "theta")) {
        if (nchar(prior$prior[i]) && !grepl("^dirichlet\\(", prior$prior[i])) {
          stop2(
            "Currently 'dirichlet' is the only valid prior for ", 
            "simplex parameters. See help(set_prior) for more details."
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

check_prior_special <- function(x, ...) {
  # prepare special priors for use in Stan
  UseMethod("check_prior_special")
}

#' @export
check_prior_special.default <- function(x, prior = empty_brmsprior(), ...) {
  prior
}

#' @export
check_prior_special.brmsprior <- function(x, bterms, ...) {
  if (!nrow(x) || isTRUE(attr(x, "checked"))) {
    return(x) 
  }
  if (is.null(x$new)) x$new <- TRUE
  x$remove <- FALSE
  x <- check_prior_special(bterms, prior = x, ...)
  x <- x[!x$remove, ]
  x$new <- x$remove <- NULL
  x
}

#' @export
check_prior_special.mvbrmsterms <- function(x, prior = NULL, ...) {
  for (cl in c("b", "Intercept")) {
    # copy over the global population-level prior in MV models
    gi <- find_rows(prior, class = cl, coef = "", resp = "")
    prior$remove[gi] <- TRUE
    for (r in x$responses) {
      ri <- find_rows(prior, class = cl, coef = "", resp = r)
      if (isTRUE(!prior$new[ri] || !nzchar(prior$prior[ri]))) {
        prior$prior[ri] <- prior$prior[gi]
      }
    }
  }
  for (i in seq_along(x$terms)) {
    prior <- check_prior_special(x$terms[[i]], prior = prior, ...)
  }
  prior
}

#' @export
check_prior_special.brmsterms <- function(x, prior = NULL, ...) {
  if (is.null(prior)) {
    prior <- empty_brmsprior()
  }
  simple_sigma <- simple_sigma(x)
  for (dp in names(x$dpars)) {
    allow_autoscale <- simple_sigma && identical(dp, "mu") 
    prior <- check_prior_special(
      x$dpars[[dp]], prior, allow_autoscale = allow_autoscale, ...
    )
  }
  for (nlp in names(x$nlpars)) {
    prior <- check_prior_special(
      x$nlpars[[nlp]], prior, is_nlpar = TRUE, ...
    )
  }
  # copy over the global population-level prior in categorical models
  if (is_categorical(x$family)) {
    for (cl in c("b", "Intercept")) {
      gi <- find_rows(
        prior, class = cl, coef = "", dpar = "", resp = x$resp
      )
      prior$remove[gi] <- TRUE
      for (dp in names(x$dpars)) {
        dpi <- find_rows(
          prior, class = cl, coef = "", dpar = dp, resp = x$resp
        )
        if (isTRUE(!prior$new[dpi] || !nzchar(prior$prior[dpi]))) {
          prior$prior[dpi] <- prior$prior[gi]
        }
      }
    }
  }
  # prepare priors for mixture probabilities
  if (is.mixfamily(x$family)) {
    take <- find_rows(prior, class = "theta", resp = x$resp)
    theta_prior <- prior$prior[take]
    if (isTRUE(nzchar(theta_prior))) {
      # hard code prior concentration
      theta_prior <- paste0(eval2(theta_prior), collapse = ", ")
      prior$prior[take] <- paste0("dirichlet(c(", theta_prior, "))")
    }
  }
  prior
}

#' @export
check_prior_special.btnl <- function(x, prior, ...) {
  prior
}

#' @export
check_prior_special.btl <- function(x, prior, data,
                                    check_nlpar_prior = TRUE,
                                    allow_autoscale = TRUE, ...) {
  # prepare special priors that cannot be passed to Stan as is
  # Args:
  #   prior: an object of class brmsprior
  #   allow_autoscale: allow autoscaling using sigma?
  #   check_nlpar_prior: check for priors on non-linear parameters?
  # Returns:
  #   a possibly amended brmsprior object with additional attributes
  px <- check_prefix(x)
  if (is_nlpar(x) && check_nlpar_prior) {
    nlp_prior <- subset2(prior, ls = px)
    if (!any(nzchar(nlp_prior$prior))) {
      stop2(
        "Priors on population-level effects are required in ",
        "non-linear models, but none were found for parameter ", 
        "'", px$nlpar, "'. See help(set_prior) for more details."
      )
    }
  }
  # prepare priors of monotonic effects
  spef <- tidy_spef(x, data)
  monef <- spef[lengths(spef$call_mo) > 0, "coef"]
  for (mo in monef) {
    take <- find_rows(prior, class = "simo", coef = mo, ls = px)
    simo_prior <- prior$prior[take]
    if (isTRUE(nzchar(simo_prior))) {
      # hard code prior concentration 
      # in order not to depend on external objects
      simo_prior <- paste0(eval2(simo_prior), collapse = ", ")
      prior$prior[take] <- paste0("dirichlet(c(", simo_prior, "))")
    }
  }
  # prepare special priors such as horseshoe or lasso
  prior_special <- list()
  b_index <- which(find_rows(prior, class = "b", coef = "", ls = px))
  stopifnot(length(b_index) <= 1L)
  if (length(b_index)) {
    b_prior <- prior$prior[b_index]
    if (any(grepl("^(horseshoe|lasso)\\(", b_prior))) {
      # horseshoe prior for population-level parameters
      if (any(nzchar(prior[b_index, "bound"]))) {
        stop2("Boundaries for population-level effects are not ", 
              "allowed when using the horseshoe or lasso priors.")
      }
      if (is.formula(x[["cs"]])) {
        stop2("Horseshoe or lasso priors are not yet allowed ",
              "in models with category-specific effects.")
      }
      b_coef_indices <- which(
        find_rows(prior, class = "b", ls = px) &
          !find_rows(prior, coef = c("", "Intercept"))
      )
      if (any(nchar(prior$prior[b_coef_indices]))) {
        stop2(
          "Defining priors for single population-level parameters ",
          "is not allowed when using horseshoe or lasso priors ",
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
        names(hs_att) <- paste0("hs_", hs_obj_names)
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

rcols_prior <- function() {
  # relevant columns for duplication checks in brmsprior objects
  c("class", "coef", "group", "resp", "dpar", "nlpar")
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
#'   See the documentation of \code{\link{brm}} for instructions
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
#' @seealso \code{\link{set_prior}}
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
#' @param df Degrees of freedom of the chi-square prior of the inverse tuning
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
#'   If you do not want to standardized all variables,
#'   you can adjust the general scale of the lasso prior via argument
#'   \code{scale}, for instance, \code{lasso(1, scale = 10)}.
#' 
#' @references
#' Park, T., & Casella, G. (2008). The Bayesian Lasso. Journal of the American 
#'    Statistical Association, 103(482), 681-686.
#'    
#' @seealso \code{\link{set_prior}}
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
