#' Prior Definitions for \pkg{brms} Models
#'
#' Define priors for specific parameters or classes of parameters.
#'
#' @aliases brmsprior brmsprior-class
#'
#' @param prior A character string defining a distribution in \pkg{Stan} language
#' @param class The parameter class. Defaults to \code{"b"}
#'   (i.e. population-level effects).
#'   See 'Details' for other valid parameter classes.
#' @param coef Name of the coefficient within the parameter class.
#' @param group Grouping factor for group-level parameters.
#' @param resp Name of the response variable.
#'   Only used in multivariate models.
#' @param dpar Name of a distributional parameter.
#'   Only used in distributional models.
#' @param nlpar Name of a non-linear parameter.
#'   Only used in non-linear models.
#' @param lb Lower bound for parameter restriction. Currently only allowed
#'   for classes \code{"b"}. Defaults to \code{NULL}, that is no restriction.
#' @param ub Upper bound for parameter restriction. Currently only allowed
#'   for classes \code{"b"}. Defaults to \code{NULL}, that is no restriction.
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
#'   in the Stan Reference Manual available at \url{https://mc-stan.org/}.
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
#'   Below, we provide details for the individual parameter classes that you can
#'   set priors on. Often, it may not be immediately clear, which parameters are
#'   present in the model. To get a full list of parameters and parameter
#'   classes for which priors can be specified (depending on the model) use
#'   function \code{\link[brms:default_prior.default]{default_prior}}.
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
#'   \code{set_prior("student_t(10, 0, 1)", class = "b", coef = "x2")}.
#'   To put the same prior on all population-level effects at once,
#'   we may write as a shortcut \code{set_prior("<prior>", class = "b")}.
#'   This also leads to faster sampling, because priors can be vectorized in this case.
#'   Both ways of defining priors can be combined using for instance
#'   \code{set_prior("normal(0, 2)", class = "b")} and \cr
#'   \code{set_prior("normal(0, 10)", class = "b", coef = "x1")}
#'   at the same time. This will set a \code{normal(0, 10)} prior on
#'   the effect of \code{x1} and a \code{normal(0, 2)} prior
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
#'   use \code{0 + Intercept} on the right-hand side of the model formula.
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
#'   The same goes for specifying priors for specific distributional
#'   parameters in the context of distributional regression, for example,
#'   \code{set_prior("<prior>", coef = "x", dpar = "sigma")}.
#'   For most other parameter classes (see below), you need to indicate
#'   non-linear and distributional parameters in the same way as shown here.
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
#'   2. Group-level ('random') effects
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
#'   after applying the link function. Minimally, the scale parameter is 2.5.
#'   This prior is used (a) to be only weakly informative in order to influence
#'   results as few as possible, while (b) providing at least some regularization
#'   to considerably improve convergence and sampling efficiency.
#'   To define a prior distribution only for standard deviations
#'   of a specific grouping factor,
#'   use \cr \code{set_prior("<prior>", class = "sd", group = "<group>")}.
#'   To define a prior distribution only for a specific standard deviation
#'   of a specific grouping factor, you may write \cr
#'   \code{set_prior("<prior>", class = "sd", group = "<group>", coef = "<coef>")}.
#'
#'   If there is more than one group-level effect per grouping factor,
#'   the correlations between those effects have to be estimated.
#'   The prior \code{lkj_corr_cholesky(eta)} or in short
#'   \code{lkj(eta)} with \code{eta > 0}
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
#'   4. Smoothing Splines
#'
#'   Smoothing splines are implemented in \pkg{brms} using the 'random effects'
#'   formulation as explained in \code{\link[mgcv:gamm]{gamm}}). Thus, each
#'   spline has its corresponding standard deviations modeling the variability
#'   within this term. In \pkg{brms}, this parameter class is called \code{sds}
#'   and priors can be specified via
#'   \code{set_prior("<prior>", class = "sds", coef = "<term label>")}.
#'   The default prior is the same as for standard deviations of group-level effects.
#'
#'   5. Gaussian processes
#'
#'   Gaussian processes as currently implemented in \pkg{brms} have two
#'   parameters, the standard deviation parameter \code{sdgp}, and
#'   characteristic length-scale parameter \code{lscale} (see \code{\link{gp}}
#'   for more details). The default prior of \code{sdgp} is the same as for
#'   standard deviations of group-level effects. The default prior of
#'   \code{lscale} is an informative inverse-gamma prior specifically tuned to
#'   the covariates of the Gaussian process (for more details see
#'   \url{https://betanalpha.github.io/assets/case_studies/gp_part3/part3.html}).
#'   This tuned prior may be overly informative in some cases, so please
#'   consider other priors as well to make sure inference is robust to the prior
#'   specification. If tuning fails, a half-normal prior is used instead.
#'
#'   6. Autocorrelation parameters
#'
#'   The autocorrelation parameters currently implemented are named \code{ar}
#'   (autoregression), \code{ma} (moving average), \code{sderr} (standard
#'   deviation of latent residuals in latent ARMA models), \code{cosy} (compound
#'   symmetry correlation), \code{car} (spatial conditional autoregression), as
#'   well as \code{lagsar} and \code{errorsar} (spatial simultaneous
#'   autoregression).
#'
#'   Priors can be defined by \code{set_prior("<prior>", class = "ar")} for
#'   \code{ar} and similar for other autocorrelation parameters. By default,
#'   \code{ar} and \code{ma} are bounded between \code{-1} and \code{1};
#'   \code{cosy}, \code{car}, \code{lagsar}, and \code{errorsar} are bounded
#'   between \code{0} and \code{1}. The default priors are flat over the
#'   respective definition areas.
#'
#'   7. Parameters of measurement error terms
#'
#'   Latent variables induced via measurement error \code{\link{me}} terms
#'   require both mean and standard deviation parameters, whose prior classes
#'   are named \code{"meanme"} and \code{"sdme"}, respectively. If multiple
#'   latent variables are induced this way, their correlation matrix will
#'   be modeled as well and corresponding priors can be specified via the
#'   \code{"corme"} class. All of the above parameters have flat priors over
#'   their respective definition spaces by default.
#'
#'   8. Distance parameters of monotonic effects
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
#'   9. Parameters for specific families
#'
#'   Some families need additional parameters to be estimated.
#'   Families \code{gaussian}, \code{student}, \code{skew_normal},
#'   \code{lognormal}, and \code{gen_extreme_value} need the parameter
#'   \code{sigma} to account for the residual standard deviation.
#'   By default, \code{sigma} has a half student-t prior that scales
#'   in the same way as the group-level standard deviations.
#'   Further, family \code{student} needs the parameter
#'   \code{nu} representing the degrees of freedom of students-t distribution.
#'   By default, \code{nu} has prior \code{gamma(2, 0.1)}
#'   and a fixed lower bound of \code{1}.
#'   Families \code{gamma}, \code{weibull}, \code{inverse.gaussian}, and
#'   \code{negbinomial} need a \code{shape} parameter that has a
#'   \code{gamma(0.01, 0.01)} prior by default.
#'   For families \code{cumulative}, \code{cratio}, \code{sratio},
#'   and \code{acat}, and only if \code{threshold = "equidistant"},
#'   the parameter \code{delta} is used to model the distance between
#'   two adjacent thresholds.
#'   By default, \code{delta} has an improper flat prior over the reals.
#'   The \code{von_mises} family needs the parameter \code{kappa}, representing
#'   the concentration parameter. By default, \code{kappa} has prior
#'   \code{gamma(2, 0.01)}.
#'
#'   Every family specific parameter has its own prior class, so that
#'   \code{set_prior("<prior>", class = "<parameter>")} is the right way to go.
#'   All of these priors are chosen to be weakly informative,
#'   having only minimal influence on the estimations,
#'   while improving convergence and sampling efficiency.
#'
#'   10. Shrinkage priors
#'
#'   To reduce the danger of overfitting in models with many predictor terms fit
#'   on comparably sparse data, brms supports special shrinkage priors, namely
#'   the (regularized) \code{\link{horseshoe}} and the \code{\link{R2D2}} prior.
#'   These priors can be applied on many parameter classes, either directly on
#'   the coefficient classes (e.g., class \code{b}), if directly setting priors
#'   on them is supported, or on the corresponding standard deviation
#'   hyperparameters (e.g., class \code{sd}) otherwise. Currently, the following
#'   classes support shrinkage priors: \code{b} (overall regression
#'   coefficients), \code{sds} (SDs of smoothing splines), \code{sdgp} (SDs of
#'   Gaussian processes), \code{ar} (autoregressive coefficients), \code{ma}
#'   (moving average coefficients), \code{sderr} (SD of latent residuals),
#'   \code{sdcar} (SD of spatial CAR structures), \code{sd} (SD of varying
#'   coefficients).
#'
#'   11. Fixing parameters to constants
#'
#'   Fixing parameters to constants is possible by using the \code{constant}
#'   function, for example, \code{constant(1)} to fix a parameter to 1.
#'   Broadcasting to vectors and matrices is done automatically.
#'
#' @seealso \code{\link[brms:default_prior.default]{default_prior}}
#'
#' @examples
#' ## use alias functions
#' (prior1 <- prior(cauchy(0, 1), class = sd))
#' (prior2 <- prior_(~cauchy(0, 1), class = ~sd))
#' (prior3 <- prior_string("cauchy(0, 1)", class = "sd"))
#' identical(prior1, prior2)
#' identical(prior1, prior3)
#'
#' # check which parameters can have priors
#' default_prior(rating ~ treat + period + carry + (1|subject),
#'              data = inhaler, family = cumulative())
#'
#' # define some priors
#' bprior <- c(prior_string("normal(0,10)", class = "b"),
#'             prior(normal(1,2), class = b, coef = treat),
#'             prior_(~cauchy(0,2), class = ~sd,
#'                    group = ~subject, coef = ~Intercept))
#'
#' # verify that the priors indeed found their way into Stan's model code
#' stancode(rating ~ treat + period + carry + (1|subject),
#'          data = inhaler, family = cumulative(),
#'          prior = bprior)
#'
#' # use the horseshoe prior to model sparsity in regression coefficients
#' stancode(count ~ zAge + zBase * Trt,
#'          data = epilepsy, family = poisson(),
#'          prior = set_prior("horseshoe(3)"))
#'
#' # fix certain priors to constants
#' bprior <- prior(constant(1), class = "b") +
#'   prior(constant(2), class = "b", coef = "zBase") +
#'   prior(constant(0.5), class = "sd")
#' stancode(count ~ zAge + zBase + (1 | patient),
#'               data = epilepsy, prior = bprior)
#'
#' # pass priors to Stan without checking
#' prior <- prior_string("target += normal_lpdf(b[1] | 0, 1)", check = FALSE)
#' stancode(count ~ Trt, data = epilepsy, prior = prior)
#'
#' # define priors in a vectorized manner
#' # useful in particular for categorical or multivariate models
#' set_prior("normal(0, 2)", dpar = c("muX", "muY", "muZ"))
#'
#' @export
set_prior <- function(prior, class = "b", coef = "", group = "",
                      resp = "", dpar = "", nlpar = "",
                      lb = NA, ub = NA, check = TRUE) {
  input <- nlist(prior, class, coef, group, resp, dpar, nlpar, lb, ub, check)
  input <- try(as.data.frame(input), silent = TRUE)
  if (is_try_error(input)) {
    stop2("Processing arguments of 'set_prior' has failed:\n", input)
  }
  out <- vector("list", nrow(input))
  for (i in seq_along(out)) {
    out[[i]] <- do_call(.set_prior, input[i, ])
  }
  Reduce("+", out)
}

# validate arguments passed to 'set_prior'
.set_prior <- function(prior, class, coef, group, resp,
                       dpar, nlpar, lb, ub, check) {
  prior <- as_one_character(prior)
  class <- as_one_character(class)
  group <- as_one_character(group)
  coef <- as_one_character(coef)
  resp <- as_one_character(resp)
  dpar <- as_one_character(dpar)
  nlpar <- as_one_character(nlpar)
  check <- as_one_logical(check)
  lb <- as_one_character(lb, allow_na = TRUE)
  ub <- as_one_character(ub, allow_na = TRUE)
  if (dpar == "mu") {
    # distributional parameter 'mu' is currently implicit #1368
    dpar <- ""
  }
  if (!check) {
    # prior will be added to the log-posterior as is
    class <- coef <- group <- resp <- dpar <- nlpar <- lb <- ub <- ""
  }
  source <- "user"
  out <- nlist(prior, source, class, coef, group, resp, dpar, nlpar, lb, ub)
  do_call(brmsprior, out)
}

#' @describeIn set_prior Alias of \code{set_prior} allowing to
#'   specify arguments as expressions without quotation marks.
#' @export
prior <- function(prior, ...) {
  call <- as.list(match.call()[-1])
  seval <- rmNULL(call[prior_seval_args()])
  call[prior_seval_args()] <- NULL
  call <- lapply(call, deparse_no_string)
  do_call(set_prior, c(call, seval))
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
  do_call(set_prior, c(call, seval))
}

# arguments for which to use standard evaluation
prior_seval_args <- function() {
  c("check")
}

#' @describeIn set_prior Alias of \code{set_prior} allowing to
#'   specify arguments as strings.
#' @export
prior_string <- function(prior, ...) {
  set_prior(prior, ...)
}

#' Default priors for Bayesian models
#'
#' @description \code{default_prior} is a generic function that can be used to
#'   get default priors for Bayesian models. It's original use is
#'   within the \pkg{brms} package, but new methods for use
#'   with objects from other packages can be registered to the same generic.
#'
#' @param object An object whose class will determine which method will
#'   be used. A symbolic description of the model to be fitted.
#' @param formula Synonym of \code{object} for use in \code{get_prior}.
#' @param ... Further arguments passed to the specific method.
#'
#' @return Usually, a \code{brmsprior} object. See
#'   \code{\link{default_prior.default}} for more details.
#'
#' @details
#' See \code{\link{default_prior.default}} for the default method applied for
#' \pkg{brms} models. You can view the available methods by typing
#' \code{methods(default_prior)}.
#'
#' @seealso
#'   \code{\link{set_prior}}, \code{\link{default_prior.default}}
#'
#' @examples
#' ## get all parameters and parameters classes to define priors on
#' (prior <- default_prior(count ~ zAge + zBase * Trt + (1|patient) + (1|obs),
#'                         data = epilepsy, family = poisson()))
#'
#' @export
default_prior <- function(object, ...) {
  UseMethod("default_prior")
}

#' @rdname default_prior
#' @export
get_prior <- function(formula, ...) {
  # became an alias of default_prior in 2.20.14.
  default_prior(formula, ...)
}

#' Default Priors for \pkg{brms} Models
#'
#' Get information on all parameters (and parameter classes) for which priors
#' may be specified including default priors.
#'
#' @inheritParams brm
#' @param object An object of class \code{\link[stats:formula]{formula}},
#'   \code{\link{brmsformula}}, or \code{\link{mvbrmsformula}} (or one that can
#'   be coerced to that classes): A symbolic description of the model to be
#'   fitted. The details of model specification are explained in
#'   \code{\link{brmsformula}}.
#' @param ... Other arguments for internal usage only.
#'
#' @return A \code{brmsprior} object. That is, a data.frame with specific
#'   columns including \code{prior}, \code{class}, \code{coef}, and \code{group}
#'   and several rows, each providing information on a parameter (or parameter
#'   class) on which priors can be specified. The prior column is empty except
#'   for internal default priors.
#'
#' @seealso \code{\link{default_prior}}, \code{\link{set_prior}}
#'
#' @examples
#' # get all parameters and parameters classes to define priors on
#' (prior <- default_prior(count ~ zAge + zBase * Trt + (1|patient) + (1|obs),
#'                         data = epilepsy, family = poisson()))
#'
#' # define a prior on all population-level effects a once
#' prior$prior[1] <- "normal(0,10)"
#'
#' # define a specific prior on the population-level effect of Trt
#' prior$prior[5] <- "student_t(10, 0, 5)"
#'
#' # verify that the priors indeed found their way into Stan's model code
#' stancode(count ~ zAge + zBase * Trt + (1|patient) + (1|obs),
#'          data = epilepsy, family = poisson(),
#'          prior = prior)
#'
#' @export
default_prior.default <- function(object, data, family = gaussian(), autocor = NULL,
                                  data2 = NULL, knots = NULL, drop_unused_levels = TRUE,
                                  sparse = NULL, ...) {

  object <- validate_formula(
    object, data = data, family = family,
    autocor = autocor, sparse = sparse
  )
  bterms <- brmsterms(object)
  data2 <- validate_data2(
    data2, bterms = bterms,
    get_data2_autocor(object)
  )
  data <- validate_data(
    data, bterms = bterms,
    data2 = data2, knots = knots,
    drop_unused_levels = drop_unused_levels
  )
  .default_prior(bterms, data, ...)
}

# internal work function of 'default_prior'
# @param internal return priors for internal use?
# @return a brmsprior object
.default_prior <- function(bterms, data, internal = FALSE, ...) {
  ranef <- tidy_ranef(bterms, data)
  meef <- tidy_meef(bterms, data)
  # initialize output
  prior <- empty_prior()
  # priors for distributional parameters
  prior <- prior + prior_predictor(
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
  # explicitly label default priors as such
  prior$source <- "default"
  # apply 'unique' as the same prior may have been included multiple times
  to_order <- with(prior, order(resp, dpar, nlpar, class, group, coef))
  prior <- unique(prior[to_order, , drop = FALSE])
  rownames(prior) <- NULL
  class(prior) <- c("brmsprior", "data.frame")
  prior
}

# generate priors for predictor terms
# @return a 'brmsprior' object
prior_predictor <- function(x, ...) {
  UseMethod("prior_predictor")
}

#' @export
prior_predictor.default <- function(x, ...) {
  empty_prior()
}

prior_predictor.mvbrmsterms <- function(x, internal = FALSE, ...) {
  prior <- empty_prior()
  for (i in seq_along(x$terms)) {
    prior <- prior + prior_predictor(x$terms[[i]], internal = internal, ...)
  }
  for (cl in c("b", "Intercept")) {
    # deprecated; see warning in 'validate_special_prior'
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
      prior <- prior +
        brmsprior(class = "nu", prior = "gamma(2, 0.1)", lb = "1")
    }
  }
  prior
}

prior_predictor.brmsterms <- function(x, data, internal = FALSE, ...) {
  data <- subset_data(data, x)
  def_scale_prior <- def_scale_prior(x, data)
  valid_dpars <- valid_dpars(x)
  prior <- empty_prior()
  # priors for mixture models
  if (is.mixfamily(x$family)) {
    if (has_joint_theta(x)) {
      # individual theta parameters should not have a prior in this case
      theta_dpars <- str_subset(valid_dpars, "^theta[[:digit:]]+")
      valid_dpars <- setdiff(valid_dpars, theta_dpars)
      prior <- prior +
        brmsprior(prior = "dirichlet(1)", class = "theta", resp = x$resp)
    }
    if (fix_intercepts(x)) {
      # fixing thresholds across mixture components
      # requires a single set of priors at the top level
      stopifnot(is_ordinal(x))
      prior <- prior + prior_thres(x, def_scale_prior = def_scale_prior)
    }
  }
  # priors for distributional parameters
  for (dp in valid_dpars) {
    def_dpar_prior <- def_dpar_prior(x, dp, data = data)
    if (!is.null(x$dpars[[dp]])) {
      # parameter is predicted
      dp_prior <- prior_predictor(
        x$dpars[[dp]], data = data,
        def_scale_prior = def_scale_prior,
        def_dpar_prior = def_dpar_prior,
        internal = internal
      )
    } else if (!is.null(x$fdpars[[dp]])) {
      # parameter is fixed
      dp_prior <- empty_prior()
    } else {
      # parameter is estimated
      dp_bound <- dpar_bounds(dp, suffix = x$resp, family = x$family)
      dp_prior <- brmsprior(
        def_dpar_prior, class = dp, resp = x$resp,
        lb = dp_bound$lb, ub = dp_bound$ub
      )
    }
    prior <- prior + dp_prior
  }
  # priors for non-linear parameters
  for (nlp in names(x$nlpars)) {
    nlp_prior <- prior_predictor(
      x$nlpars[[nlp]], data = data,
      def_scale_prior = def_scale_prior,
      internal = internal
    )
    prior <- prior + nlp_prior
  }
  if (is_logistic_normal(x$family)) {
    if (internal) {
      prior <- prior +
        brmsprior("lkj_corr_cholesky(1)", class = "Llncor", resp = x$resp)
    } else {
      prior <- prior +
        brmsprior("lkj(1)", class = "lncor", resp = x$resp)
    }
  }
  prior
}

# prior for linear predictor termss
#' @export
prior_predictor.btl <- function(x, ...) {
  prior_fe(x, ...) +
    prior_thres(x, ...) +
    prior_sp(x, ...) +
    prior_cs(x, ...) +
    prior_sm(x, ...) +
    prior_gp(x, ...) +
    prior_ac(x, ...) +
    prior_bhaz(x, ...)
}

# priors for non-linear predictor terms
#' @export
prior_predictor.btnl <- function(x, ...) {
  # thresholds are required even in non-linear ordinal models
  prior_thres(x, ...) +
    prior_ac(x, ...) +
    prior_bhaz(x, ...)
}

# priors for population-level parameters
prior_fe <- function(bterms, data, def_dpar_prior = "", ...) {
  prior <- empty_prior()
  fixef <- colnames(data_fe(bterms, data)$X)
  px <- check_prefix(bterms)
  center_X <- stan_center_X(bterms)
  if (center_X && !is_ordinal(bterms)) {
    # priors for ordinal thresholds are provided in 'prior_thres'
    prior <- prior + brmsprior(def_dpar_prior, class = "Intercept", ls = px)
    fixef <- setdiff(fixef, "Intercept")
  }
  if (length(fixef)) {
    prior <- prior + brmsprior(class = "b", coef = c("", fixef), ls = px)
  }
  prior
}

# priors for thresholds of ordinal models
prior_thres <- function(bterms, def_scale_prior = "", ...) {
  prior <- empty_prior()
  if (!is_ordinal(bterms)) {
    # thresholds only exist in ordinal models
    return(prior)
  }
  if (fix_intercepts(bterms) && !is.mixfamily(bterms$family)) {
    # fixed thresholds cannot have separate priors
    return(prior)
  }

  # create priors for threshold per group
  .prior_thres <- function(thres, thres_prior = "", group = "") {
    prior <- empty_prior()
    if (has_equidistant_thres(bterms)) {
      # prior for the delta parameter for equidistant thresholds
      thres <- character(0)
      lb <- str_if(has_ordered_thres(bterms), "0")
      prior <- prior + brmsprior(
        class = "delta", group = group, lb = lb, ls = px
      )
    }
    prior <- prior + brmsprior(
      prior = c(thres_prior, rep("", length(thres))),
      class = "Intercept", coef = c("", thres),
      group = group, ls = px
    )
  }

  px <- check_prefix(bterms)
  groups <- get_thres_groups(bterms)
  if (any(nzchar(groups))) {
    # for models with multiple threshold vectors
    prior <- prior + .prior_thres(character(0), def_scale_prior)
    for (g in groups) {
      prior <- prior + .prior_thres(get_thres(bterms, group = g), group = g)
    }
  } else {
    # for models with a single threshold vector
    prior <- prior + .prior_thres(get_thres(bterms), def_scale_prior)
  }
  prior
}

# priors for coefficients of baseline hazards in the Cox model
prior_bhaz <- function(bterms, ...) {
  prior <- empty_prior()
  if (!is_cox(bterms$family)) {
    return(prior)
  }
  px <- check_prefix(bterms)
  # the scale of sbhaz is not identified when an intercept is part of mu
  # thus a sum-to-one constraint ensures identification
  prior <- prior + brmsprior("dirichlet(1)", class = "sbhaz", ls = px)
  prior
}

# priors for special effects parameters
prior_sp <- function(bterms, data, ...) {
  prior <- empty_prior()
  spef <- tidy_spef(bterms, data)
  if (nrow(spef)) {
    px <- check_prefix(bterms)
    prior <- prior + brmsprior(
      class = "b", coef = c("", spef$coef), ls = px
    )
    simo_coef <- get_simo_labels(spef, use_id = TRUE)
    if (length(simo_coef)) {
      prior <- prior + brmsprior(
        prior  = "dirichlet(1)", class = "simo",
        coef = simo_coef, ls = px
      )
    }
  }
  prior
}

# priors for category spcific effects parameters
prior_cs <- function(bterms, data, ...) {
  prior <- empty_prior()
  csef <- colnames(get_model_matrix(bterms$cs, data = data))
  if (length(csef)) {
    px <- check_prefix(bterms)
    prior <- prior +
      brmsprior(class = "b", coef = c("", csef), ls = px)
  }
  prior
}

# default priors for hyper-parameters of noise-free variables
prior_Xme <- function(meef, internal = FALSE, ...) {
  stopifnot(is.meef_frame(meef))
  prior <- empty_prior()
  if (nrow(meef)) {
    prior <- prior +
      brmsprior(class = "meanme") +
      brmsprior(class = "meanme", coef = meef$coef) +
      brmsprior(class = "sdme", lb = "0") +
      brmsprior(class = "sdme", coef = meef$coef)
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

# default priors of gaussian processes
# @param def_scale_prior: a character string defining
#   the default prior SD parameters
prior_gp <- function(bterms, data, def_scale_prior, ...) {
  prior <- empty_prior()
  gpef <- tidy_gpef(bterms, data)
  if (nrow(gpef)) {
    px <- check_prefix(bterms)
    lscale_prior <- def_lscale_prior(bterms, data)
    prior <- prior +
      brmsprior(class = "sdgp", prior = def_scale_prior, ls = px,
                lb = "0") +
      brmsprior(class = "sdgp", coef = unlist(gpef$sfx1), ls = px) +
      brmsprior(class = "lscale", ls = px, lb = "0") +
      brmsprior(class = "lscale", prior = lscale_prior,
                coef = names(lscale_prior), ls = px)
  }
  prior
}

# default priors for length-scale parameters of GPs
# see https://betanalpha.github.io/assets/case_studies/gp_part3/part3.html
# @param plb prior probability of being lower than minimum length-scale
# @param pub prior probability of being higher than maximum length-scale
def_lscale_prior <- function(bterms, data, plb = 0.01, pub = 0.01) {
  .opt_fun <- function(x, lb, ub) {
    # optimize parameters on the log-scale to make them positive only
    x <- exp(x)
    y1 <- pinvgamma(lb, x[1], x[2], log.p = TRUE)
    y2 <- pinvgamma(ub, x[1], x[2], lower.tail = FALSE, log.p = TRUE)
    c(y1 - log(plb), y2 - log(pub))
  }
  .def_lscale_prior <- function(X) {
    dq <- diff_quad(X)
    ub <- sqrt(max(dq))
    lb <- sqrt(min(dq[dq > 0]))
    # prevent extreme priors
    lb <- max(lb, 0.01 * ub)
    opt_res <- nleqslv::nleqslv(
      c(0, 0), .opt_fun, lb = lb, ub = ub,
      control = list(allowSingular = TRUE)
    )
    prior <- "normal(0, 0.5)"
    if (opt_res$termcd %in% 1:2) {
      # use the inverse-gamma prior only in case of convergence
      pars <- exp(opt_res$x)
      prior <- paste0("inv_gamma(", sargs(round(pars, 6)), ")")
    }
    return(prior)
  }
  p <- usc(combine_prefix(bterms))
  gpef <- tidy_gpef(bterms, data)
  data_gp <- data_gp(bterms, data, internal = TRUE)
  out <- vector("list", NROW(gpef))
  for (i in seq_along(out)) {
    pi <- paste0(p, "_", i)
    iso <- gpef$iso[i]
    cons <- gpef$cons[[i]]
    if (length(cons) > 0L) {
      for (j in seq_along(cons)) {
        Xgp <- data_gp[[paste0("Xgp_prior", pi, "_", j)]]
        if (iso) {
          c(out[[i]]) <- .def_lscale_prior(Xgp)
        } else {
          c(out[[i]]) <- apply(Xgp, 2, .def_lscale_prior)
        }
      }
    } else {
      Xgp <- data_gp[[paste0("Xgp_prior", pi)]]
      if (iso) {
        out[[i]] <- .def_lscale_prior(Xgp)
      } else {
        out[[i]] <- apply(Xgp, 2, .def_lscale_prior)
      }
    }
    # transpose so that by-levels vary last
    names(out[[i]]) <- as.vector(t(gpef$sfx2[[i]]))
  }
  unlist(out)
}

# priors for varying effects parameters
# @param ranef: a list returned by tidy_ranef
# @param def_scale_prior a character string defining
#   the default prior for SD parameters
# @param internal: see 'default_prior'
prior_re <- function(ranef, def_scale_prior, internal = FALSE, ...) {
  prior <- empty_prior()
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
    class = "sd", prior = def_scale_prior,
    lb = "0", ls = px
  )
  prior <- prior + global_sd_prior
  for (id in unique(ranef$id)) {
    r <- subset2(ranef, id = id)
    group <- r$group[1]
    rpx <- check_prefix(r)
    urpx <- unique(rpx)
    # include group-level standard deviations
    prior <- prior +
      # don't specify lb as we already have it above
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
      brmsprior("gamma(2, 0.1)", class = "df", group = tranef$group,
                lb = "1")
  }
  prior
}

# priors for smooth terms
prior_sm <- function(bterms, data, def_scale_prior, ...) {
  prior <- empty_prior()
  smef <- tidy_smef(bterms, data)
  if (NROW(smef)) {
    px <- check_prefix(bterms)
    # prior for the FE coefficients
    Xs_names <- attr(smef, "Xs_names")
    if (length(Xs_names)) {
      prior <- prior + brmsprior(
        class = "b", coef = c("", Xs_names), ls = px
      )
    }
    # prior for SD parameters of the RE coefficients
    smterms <- unique(smef$term)
    prior <- prior +
      brmsprior(prior = def_scale_prior, class = "sds",
                lb = "0", ls = px) +
      brmsprior(class = "sds", coef = smterms, ls = px)
  }
  prior
}

# priors for autocor parameters
prior_ac <- function(bterms, def_scale_prior, internal = FALSE, ...) {
  prior <- empty_prior()
  acef <- tidy_acef(bterms)
  if (!NROW(acef)) {
    return(prior)
  }
  px <- check_prefix(bterms)
  p <- combine_prefix(px)
  has_ac_latent_residuals <- has_ac_latent_residuals(bterms)
  if (has_ac_class(acef, "arma")) {
    acef_arma <- subset2(acef, class = "arma")
    # no boundaries are required in the conditional formulation
    # when natural residuals automatically define the scale
    need_arma_bound <- acef_arma$cov || has_ac_latent_residuals
    arma_lb <- str_if(need_arma_bound, "-1")
    arma_ub <- str_if(need_arma_bound, "1")
    if (acef_arma$p > 0) {
      prior <- prior +
        brmsprior(class = "ar", ls = px, lb = arma_lb, ub = arma_ub)
    }
    if (acef_arma$q > 0) {
      prior <- prior +
        brmsprior(class = "ma", ls = px, lb = arma_lb, ub = arma_ub)
    }
  }
  if (has_ac_class(acef, "cosy")) {
    # cosy correlations may be negative in theory but
    # this causes problems with divergent transitions (#878)
    prior <- prior +
      brmsprior(class = "cosy", ls = px, lb = "0", ub = "1")
  }
  if (has_ac_class(acef, "unstr")) {
    if (internal) {
      prior <- prior +
        brmsprior("lkj_corr_cholesky(1)", class = "Lcortime", ls = px)
    } else {
      prior <- prior +
        brmsprior("lkj(1)", class = "cortime", ls = px)
    }
  }
  if (has_ac_latent_residuals(bterms)) {
    prior <- prior +
      brmsprior(def_scale_prior, class = "sderr", ls = px, lb = "0")
  }
  if (has_ac_class(acef, "sar")) {
    acef_sar <- subset2(acef, class = "sar")
    sar_lb <- glue("min_eigenMsar{p}")
    sar_ub <- glue("max_eigenMsar{p}")
    if (acef_sar$type == "lag") {
      prior <- prior +
        brmsprior(class = "lagsar", lb = sar_lb, ub = sar_ub, ls = px)
    }
    if (acef_sar$type == "error") {
      prior <- prior +
        brmsprior(class = "errorsar", lb = sar_lb, ub = sar_ub, ls = px)
    }
  }
  if (has_ac_class(acef, "car")) {
    acef_car <- subset2(acef, class = "car")
    prior <- prior +
      brmsprior(def_scale_prior, class = "sdcar", lb = "0", ls = px)
    if (acef_car$type %in% "escar") {
      prior <- prior +
        brmsprior(class = "car", lb = "0", ub = "1", ls = px)
    } else if (acef_car$type %in% "bym2") {
      prior <- prior +
        brmsprior("beta(1, 1)", class = "rhocar", lb = "0", ub = "1", ls = px)
    }
  }
  prior
}

# default priors for distributional parameters
def_dpar_prior <- function(x, dpar, data = NULL) {
  stopifnot(is.brmsterms(x))
  dpar <- as_one_character(dpar)
  resp <- usc(x$resp)
  dpar_class <- dpar_class(dpar, family = x)
  link <- x$dpars[[dpar]]$family$link
  if (is.null(link)) {
    link <- "identity"
  }
  # ensures reasonable scaling in def_scale_prior
  x$family$link <- link
  if (link == "identity") {
    # dpar is estimated or predicted on the linear scale
    out <- switch(dpar_class, "",
      mu = def_scale_prior(x, data, center = FALSE, dpar = dpar),
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
      ndt = glue("uniform(0, min_Y{resp})"),
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
      mu = def_scale_prior(x, data, center = FALSE, dpar = dpar),
      sigma = def_scale_prior(x, data),
      shape = "student_t(3, 0, 2.5)",
      nu = "normal(2.7, 0.8)",
      phi = "student_t(3, 0, 2.5)",
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

# default priors for scale/SD parameters
def_scale_prior <- function(x, data, ...) {
  UseMethod("def_scale_prior")
}

#' @export
def_scale_prior.mvbrmsterms <- function(x, data, ...) {
  out <- ulapply(x$terms, def_scale_prior, data = data, ...)
  names(out) <- x$responses
  out
}

# @param center Should the prior be centered around zero?
#   If FALSE, the prior location is computed based on Y.
#' @export
def_scale_prior.brmsterms <- function(x, data, center = TRUE, df = 3,
                                      location = 0, scale = 2.5,
                                      dpar = NULL, ...) {
  y <- unname(model.response(model.frame(x$respform, data)))
  link <- x$family$link
  if (has_logscale(x$family)) {
    link <- "log"
  }
  tlinks <- c("identity", "log", "inverse", "sqrt", "1/mu^2")
  if (link %in% tlinks && !is_like_factor(y) && !conv_cats_dpars(x)) {
    if (link %in% c("log", "inverse", "1/mu^2")) {
      # avoid Inf in link(y)
      y <- ifelse(y == 0, y + 0.1, y)
    }
    y_link <- SW(link(y, link = link))
    scale_y <- round(mad(y_link), 1)
    if (is.finite(scale_y)) {
      scale <- max(scale, scale_y)
    }
    if (!center) {
      location_y <- round(median(y_link), 1)
      if (is.finite(location_y)) {
        location <- location_y
      }
      # offsets may render default intercept priors not sensible
      dpar <- as_one_character(dpar)
      offset <- unname(unlist(data_offset(x$dpars[[dpar]], data)))
      if (length(offset)) {
        mean_offset <- mean(offset)
        if (is.finite(mean_offset)) {
          location <- location - mean_offset
        }
      }
    }
  }
  paste0("student_t(", sargs(df, location, scale), ")")
}

#' Validate Prior for \pkg{brms} Models
#'
#' Validate priors supplied by the user. Return a complete
#' set of priors for the given model, including default priors.
#'
#' @inheritParams default_prior.default
#' @inheritParams brm
#'
#' @return An object of class \code{brmsprior}.
#'
#' @seealso \code{\link[brms:default_prior.default]{default_prior}}, \code{\link{set_prior}}.
#'
#' @examples
#' prior1 <- prior(normal(0,10), class = b) +
#'   prior(cauchy(0,2), class = sd)
#' validate_prior(prior1, count ~ zAge + zBase * Trt + (1|patient),
#'                data = epilepsy, family = poisson())
#'
#' @export
validate_prior <- function(prior, formula, data, family = gaussian(),
                           sample_prior = "no", data2 = NULL, knots = NULL,
                           drop_unused_levels = TRUE, ...) {
  formula <- validate_formula(formula, data = data, family = family)
  bterms <- brmsterms(formula)
  data2 <- validate_data2(data2, bterms = bterms)
  data <- validate_data(
    data, bterms = bterms,
    data2 = data2, knots = knots,
    drop_unused_levels = drop_unused_levels
  )
  .validate_prior(
    prior, bterms = bterms, data = data,
    sample_prior = sample_prior, ...
  )
}

# internal work function of 'validate_prior'
.validate_prior <- function(prior, bterms, data, sample_prior, ...) {
  sample_prior <- validate_sample_prior(sample_prior)
  all_priors <- .default_prior(bterms, data, internal = TRUE)
  if (is.null(prior)) {
    prior <- all_priors
  } else if (!is.brmsprior(prior)) {
    stop2("Argument 'prior' must be a 'brmsprior' object.")
  }
  # when updating existing priors, invalid priors should be allowed
  allow_invalid_prior <- isTRUE(attr(prior, "allow_invalid_prior"))
  # temporarily exclude priors that should not be checked
  no_checks <- !nzchar(prior$class)
  prior_no_checks <- prior[no_checks, ]
  prior <- prior[!no_checks, ]
  # check for duplicated priors
  prior$class <- rename(
    prior$class,
    c("^cor$", "^rescor$", "^corme$", "^lncor$", "^cortime$"),
    c("L", "Lrescor", "Lme", "Llncor", "Lcortime"),
    fixed = FALSE
  )
  if (any(duplicated(prior))) {
    stop2("Duplicated prior specifications are not allowed.")
  }
  # check for invalid priors
  # it is good to let the user know beforehand that some of their priors
  # were invalid in the model to avoid unnecessary refits
  if (nrow(prior)) {
    valid_ids <- which(duplicated(rbind(all_priors, prior)))
    invalid <- !seq_rows(prior) %in% (valid_ids - nrow(all_priors))
    if (any(invalid) && !allow_invalid_prior) {
      stop2(
        "The following priors do not correspond ",
        "to any model parameter: \n",
        collapse(.print_prior(prior[invalid, ]), "\n"),
        "Function 'default_prior' might be helpful to you."
      )
    }
    prior <- prior[!invalid, ]
  }
  prior$prior <- sub("^(lkj|lkj_corr)\\(", "lkj_corr_cholesky(", prior$prior)

  # include default parameter bounds; only new priors need bounds
  which_needs_lb <- which(is.na(prior$lb) & !nzchar(prior$coef))
  for (i in which_needs_lb) {
    if (!is.na(prior$ub[i]) && nzchar(prior$ub[i])) {
      # if ub is specified lb should be specified in the same line as well
      prior$lb[i] <- stan_base_prior(all_priors, "lb", sel_prior = prior[i, ])
    } else {
      # take the corresponding lb from the default prior
      prior_sub_i <- rbind(prior[i, ], all_priors)
      prior_sub_i <- prior_sub_i[duplicated(prior_sub_i), ]
      stopifnot(NROW(prior_sub_i) == 1L)
      prior$lb[i] <- prior_sub_i$lb
    }
  }
  which_needs_ub <- which(is.na(prior$ub) & !nzchar(prior$coef))
  for (i in which_needs_ub) {
    if (!is.na(prior$lb[i]) && nzchar(prior$lb[i])) {
      # if lb is specified ub should be specified in the same line as well
      prior$ub[i] <- stan_base_prior(all_priors, "ub", sel_prior = prior[i, ])
    } else {
      # take the corresponding lb from the default prior
      prior_sub_i <- rbind(prior[i, ], all_priors)
      prior_sub_i <- prior_sub_i[duplicated(prior_sub_i), ]
      stopifnot(NROW(prior_sub_i) == 1L)
      prior$ub[i] <- prior_sub_i$ub
    }
  }
  # the remaining NAs are in coef priors which cannot have bounds yet
  prior$lb[is.na(prior$lb)] <- prior$ub[is.na(prior$ub)] <- ""

  # boundaries on individual coefficients are not yet supported
  # TODO: enable bounds for coefficients as well?
  if (any((nzchar(prior$lb) | nzchar(prior$ub)) & nzchar(prior$coef))) {
    stop2("Prior argument 'coef' may not be specified when using boundaries.")
  }

  # merge user-specified priors with default priors
  prior$new <- rep(TRUE, nrow(prior))
  all_priors$new <- rep(FALSE, nrow(all_priors))
  prior <- c(all_priors, prior, replace = TRUE)
  check_prior_content(prior)

  prior <- validate_special_prior(prior, bterms = bterms, data = data, ...)
  prior <- prior[with(prior, order(class, group, resp, dpar, nlpar, coef)), ]
  # check and warn valid but unused priors
  for (i in which(nzchar(prior$prior) & !nzchar(prior$coef))) {
    ls <- prior[i, c("class", "group", "resp", "dpar", "nlpar")]
    class(ls) <- "data.frame"
    prior_sub_coef <- subset2(prior, ls = ls)
    prior_sub_coef <- prior_sub_coef[nzchar(prior_sub_coef$coef), ]
    if (nrow(prior_sub_coef) && all(nzchar(prior_sub_coef$prior))) {
      warning2(
        "The global prior '", prior$prior[i], "' of class '", prior$class[i],
        "' will not be used in the model as all related coefficients have ",
        "individual priors already. If you did not set those ",
        "priors yourself, then maybe brms has assigned default priors. ",
        "See ?set_prior and ?default_prior for more details."
      )
    }
  }
  prior <- prior + prior_no_checks
  rownames(prior) <- NULL
  attr(prior, "sample_prior") <- sample_prior
  if (is_verbose()) {
    # show remaining default priors added to the model
    def_prior <- prepare_print_prior(prior)
    def_prior <- subset2(def_prior, source = "default")
    if (nrow(def_prior)) {
      message("The following priors were automatically added to the model:")
      print(def_prior, show_df = TRUE)
    }
  }
  prior
}

# try to check if prior distributions are reasonable
# @param prior A brmsprior object
check_prior_content <- function(prior) {
  if (!is.brmsprior(prior) || !NROW(prior)) {
    return(invisible(TRUE))
  }
  lb_priors <- c(
    "lognormal", "chi_square", "inv_chi_square", "scaled_inv_chi_square",
    "exponential", "gamma", "inv_gamma", "weibull", "frechet", "rayleigh",
    "pareto", "pareto_type_2"
  )
  lb_priors_regex <- paste0("^(", paste0(lb_priors, collapse = "|"), ")")
  ulb_priors <- c("beta", "uniform", "von_mises", "beta_proportion")
  ulb_priors_regex <- paste0("^(", paste0(ulb_priors, collapse = "|"), ")")
  cormat_pars <- c(
    "cor", "L", "rescor", "Lrescor", "corme", "Lme",
    "lncor", "Llncor", "cortime", "Lcortime"
  )
  cormat_regex <- "^((lkj)|(constant))"
  simplex_pars <- c("simo", "theta", "sbhaz")
  simplex_regex <- "^((dirichlet)|(constant))\\("

  lb_warning <- ub_warning <- ""
  for (i in seq_rows(prior)) {
    if (!nzchar(prior$prior[i]) || !prior$new[i]) {
      next
    }
    msg_prior <- .print_prior(prior[i, ])
    has_lb_prior <- grepl(lb_priors_regex, prior$prior[i])
    has_ulb_prior <- grepl(ulb_priors_regex, prior$prior[i])
    base_bounds <- stan_base_prior(prior, c("lb", "ub"), sel_prior = prior[i, ])
    has_lb <- nzchar(base_bounds[, "lb"])
    has_ub <- nzchar(base_bounds[, "ub"])
    if ((has_lb_prior || has_ulb_prior) && !has_lb) {
      lb_warning <- paste0(lb_warning, msg_prior, "\n")
    }
    if (has_ulb_prior && !has_ub) {
      ub_warning <- paste0(ub_warning, msg_prior, "\n")
    }
    if (prior$class[i] %in% cormat_pars &&
        !grepl(cormat_regex, prior$prior[i])) {
      stop2(
        "The only supported prior for correlation matrices is ",
        "the 'lkj' prior. See help(set_prior) for more details."
      )
    }
    if (prior$class[i] %in% simplex_pars &&
        !grepl(simplex_regex, prior$prior[i])) {
      stop2(
        "Currently 'dirichlet' is the only valid prior for ",
        "simplex parameters. See help(set_prior) for more details."
      )
    }
  }
  if (nzchar(lb_warning)) {
    warning2(
      "It appears as if you have specified a lower bounded ",
      "prior on a parameter that has no natural lower bound.",
      "\nIf this is really what you want, please specify ",
      "argument 'lb' of 'set_prior' appropriately.",
      "\nWarning occurred for prior \n", lb_warning
    )
  }
  if (nzchar(ub_warning)) {
    warning2(
      "It appears as if you have specified an upper bounded ",
      "prior on a parameter that has no natural upper bound.",
      "\nIf this is really what you want, please specify ",
      "argument 'ub' of 'set_prior' appropriately.",
      "\nWarning occurred for prior \n", ub_warning
    )
  }
  invisible(TRUE)
}

# prepare special priors for use in Stan
# required for priors that are not natively supported by Stan
validate_special_prior <- function(x, ...) {
  UseMethod("validate_special_prior")
}

#' @export
validate_special_prior.default <- function(x, prior = empty_prior(), ...) {
  prior
}

#' @export
validate_special_prior.brmsprior <- function(x, bterms, ...) {
  if (!NROW(x)) {
    return(x)
  }
  if (is.null(x$new)) {
    x$new <- TRUE
  }
  x$remove <- FALSE
  x <- validate_special_prior(bterms, prior = x, ...)
  x <- x[!x$remove, ]
  x$new <- x$remove <- NULL
  x
}

#' @export
validate_special_prior.mvbrmsterms <- function(x, prior = NULL, ...) {
  for (i in seq_along(x$terms)) {
    prior <- validate_special_prior(x$terms[[i]], prior = prior, ...)
  }
  prior
}

#' @export
validate_special_prior.brmsterms <- function(x, data, prior = NULL, ...) {
  data <- subset_data(data, x)
  if (is.null(prior)) {
    prior <- empty_prior()
  }
  simple_sigma <- simple_sigma(x)
  for (dp in names(x$dpars)) {
    allow_autoscale <- dp == "mu" && simple_sigma
    prior <- validate_special_prior(
      x$dpars[[dp]], prior = prior, data = data,
      allow_autoscale = allow_autoscale, ...
    )
  }
  for (nlp in names(x$nlpars)) {
    prior <- validate_special_prior(
      x$nlpars[[nlp]], prior = prior, data = data,
      allow_autoscale = simple_sigma, ...
    )
  }
  prior
}

#' @export
validate_special_prior.btnl <- function(x, prior, ...) {
  prior
}

# prepare special priors that cannot be passed to Stan as is
# @param allow_autoscale allow autoscaling by parameter sigma?
# @return a possibly updated brmsprior object with additional attributes
#' @export
validate_special_prior.btl <- function(x, prior, data, allow_autoscale = TRUE,
                                       ...) {
  allow_autoscale <- as_one_logical(allow_autoscale)
  px <- check_prefix(x)
  # prepare special priors such as horseshoe
  special <- list()
  # the order of the classes doesn't matter but for consistency
  # it is still the same as the order in the Stan code
  special_classes <- c("b", "sds", "sdgp", "ar", "ma", "sderr", "sdcar", "sd")
  for (sc in special_classes) {
    index <- which(find_rows(prior, class = sc, coef = "", group = "", ls = px))
    if (!length(index)) {
      next
    }
    stopifnot(length(index) <= 1L)
    sub_prior <- prior$prior[index]
    if (any(is_special_prior(sub_prior))) {
      # shrinkage priors have been specified
      if (sc %in% c("b", "ar", "ma")) {
        if (any(nzchar(prior[index, "lb"]) | nzchar(prior[index, "ub"]))) {
          stop2("Setting boundaries on coefficients is not ",
                "allowed when using special priors.")
        }
        # TODO: allow special priors also for 'cs' coefficients
        if (is.formula(x[["cs"]])) {
          stop2("Special priors are not yet allowed ",
                "in models with category-specific effects.")
        }
      }
      if (sc %in% c("sds", "sdgp", "sderr", "sdcar", "sd")) {
        if (any(prior[index, "lb"] != "0" | nzchar(prior[index, "ub"]))) {
          stop2("Setting custom boundaries on SD parameters is not ",
                "allowed when using special priors.")
        }
      }
      coef_indices <- which(
        find_rows(prior, class = sc, ls = px) &
          !find_rows(prior, class = sc, group = "", coef = "")
      )
      if (any(nzchar(prior$prior[coef_indices]))) {
        stop2(
          "Defining separate priors for single coefficients or groups is not ",
          "allowed when using special priors for the whole class."
        )
      }
      tmp <- attributes(eval2(sub_prior))
      tmp$autoscale <- isTRUE(tmp$autoscale) && allow_autoscale
      special[[sc]] <- tmp
    }
  }
  special_names <- unique(ufrom_list(special, "name"))
  if (length(special_names) > 1L) {
    stop2("Currently only one special prior per formula is allowed.")
  }
  prefix <- combine_prefix(px, keep_mu = TRUE)
  attributes(prior)$special[[prefix]] <- special
  prior
}

# validate argument 'sample_prior'
validate_sample_prior <- function(sample_prior) {
  options <- c("no", "yes", "only")
  if (is.null(sample_prior)) {
    sample_prior <- "no"
  }
  if (!is.character(sample_prior)) {
    sample_prior <- as_one_logical(sample_prior)
    sample_prior <- if (sample_prior) "yes" else "no"
  }
  match.arg(sample_prior, options)
}

# get stored 'sample_prior' argument
get_sample_prior <- function(prior) {
  validate_sample_prior(attr(prior, "sample_prior", TRUE))
}

# create data.frames containing prior information
brmsprior <- function(prior = "", class = "", coef = "", group = "",
                      resp = "", dpar = "", nlpar = "", lb = "", ub = "",
                      source = "", ls = list()) {
  if (length(ls)) {
    if (is.null(names(ls))) {
      stop("Argument 'ls' must be named.")
    }
    names <- all_cols_prior()
    if (!all(names(ls) %in% names)) {
      stop("Names of 'ls' must some of ", collapse_comma(names))
    }
    for (v in names(ls)) {
      assign(v, ls[[v]])
    }
  }
  out <- data.frame(
    prior, class, coef, group,
    resp, dpar, nlpar, lb, ub, source,
    stringsAsFactors = FALSE
  )
  class(out) <- c("brmsprior", "data.frame")
  out
}

#' @describeIn set_prior Create an empty \code{brmsprior} object.
#' @export
empty_prior <- function() {
  char0 <- character(0)
  brmsprior(
    prior = char0, source = char0, class = char0,
    coef = char0, group = char0, resp = char0,
    dpar = char0, nlpar = char0, lb = char0, ub = char0
  )
}

# natural upper and lower bounds for priors
# @param a named list with elements 'lb' and 'ub'
prior_bounds <- function(prior) {
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

# all columns of brmsprior objects
all_cols_prior <- function() {
  c("prior", "class", "coef", "group", "resp",
    "dpar", "nlpar", "lb", "ub", "source")
}

# relevant columns for duplication checks in brmsprior objects
rcols_prior <- function() {
  c("class", "coef", "group", "resp", "dpar", "nlpar")
}

# default Stan definitions for distributional parameters
# @param dpar name of a distributional parameter
# @param suffix optional suffix of the parameter name
# @param family optional brmsfamily object
# @return a named list with numeric elements 'lb' and 'ub'
dpar_bounds <- function(dpar, suffix = "", family = NULL) {
  dpar <- as_one_character(dpar)
  suffix <- usc(as_one_character(suffix))
  if (is.mixfamily(family)) {
    if (dpar_class(dpar) == "theta") {
      return(list(lb = -Inf, ub = Inf))
    }
    family <- family$mix[[as.numeric(dpar_id(dpar))]]
  }
  dpar_class <- dpar_class(dpar, family)
  if (is.customfamily(family)) {
    lb <- family$lb[[dpar_class]]
    ub <- family$ub[[dpar_class]]
    return(nlist(lb, ub))
  }
  min_Y <- glue("min_Y{suffix}")
  out <- switch(dpar_class,
    sigma = list(lb = "0", ub = ""),
    shape = list(lb = "0", ub = ""),
    nu = list(lb = "1", ub = ""),
    phi = list(lb = "0", ub = ""),
    kappa = list(lb = "0", ub = ""),
    beta = list(lb = "0", ub = ""),
    zi = list(lb = "0", ub = "1"),
    hu = list(lb = "0", ub = "1"),
    zoi = list(lb = "0", ub = "1"),
    coi = list(lb = "0", ub = "1"),
    bs = list(lb = "0", ub = ""),
    ndt = list(lb = "0", ub = min_Y),
    bias = list(lb = "0", ub = "1"),
    disc = list(lb = "0", ub = ""),
    quantile = list(lb = "0", ub = "1"),
    xi = list(lb = "", ub = ""),
    alpha = list(lb = "", ub = "")
  )
  out
}

# convert parameter bounds to Stan syntax
# vectorized over both 'lb' and 'ub' vectors
# @param bounds a named list with elements 'lb' and 'ub'
# @param default default output if no proper bounds are specified
convert_bounds2stan <- function(bounds, default = "") {
  lb <- as.character(bounds$lb)
  ub <- as.character(bounds$ub)
  stopifnot(length(lb) == length(ub))
  default <- as_one_character(default, allow_na = TRUE)
  if (any(lb %in% "Inf")) {
    stop2("Lower boundaries cannot be positive infinite.")
  }
  if (any(ub %in% "-Inf")) {
    stop2("Upper boundaries cannot be negative infinite.")
  }
  lb <- ifelse(
    !is.na(lb) & !lb %in% c("NA", "-Inf", ""),
    paste0("lower=", lb), ""
  )
  ub <- ifelse(
    !is.na(ub) & !ub %in% c("NA", "Inf", ""),
    paste0("upper=", ub), ""
  )
  out <- ifelse(
    nzchar(lb) & nzchar(ub), glue("<{lb},{ub}>"),
    ifelse(
      nzchar(lb) & !nzchar(ub), glue("<{lb}>"),
      ifelse(
        !nzchar(lb) & nzchar(ub), glue("<{ub}>"),
        default
      )
    )
  )
  out
}

# convert parameter bounds in Stan syntax
# TODO: vectorize over a character vector of bounds?
# complicated because of a mix of character and numeric values
# to a named list with elements 'lb' and 'ub'
convert_stan2bounds <- function(bound, default = c(-Inf, Inf)) {
  bound <- as_one_character(bound)
  stopifnot(length(default) == 2L)
  out <- list(lb = default[[1]], ub = default[[2]])
  if (!is.na(bound) && isTRUE(nzchar(bound))) {
    lb <- get_matches("(<|,)lower=[^,>]+", bound)
    if (isTRUE(nzchar(lb))) {
      lb <- substr(lb, 8, nchar(lb))
      lb_num <- SW(as.numeric(lb))
      if (!is.na(lb_num)) {
        lb <- lb_num
      }
      out$lb <- lb
    }
    ub <- get_matches("(<|,)upper=[^,>]+", bound)
    if (isTRUE(nzchar(ub))) {
      ub <- substr(ub, 8, nchar(ub))
      ub_num <- SW(as.numeric(ub))
      if (!is.na(ub_num)) {
        ub <- ub_num
      }
      out$ub <- ub
    }
  }
  out
}

#' Priors of \code{brms} models
#'
#' Extract priors of models fitted with \pkg{brms}.
#'
#' @aliases prior_summary
#'
#' @param object An object of class \code{brmsfit}.
#' @param all Logical; Show all parameters in the model which may have
#'   priors (\code{TRUE}) or only those with proper priors (\code{FALSE})?
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An \code{brmsprior} object.
#'
#' @examples
#' \dontrun{
#' fit <- brm(
#'   count ~ zAge + zBase * Trt + (1|patient) + (1|obs),
#'   data = epilepsy, family = poisson(),
#'   prior = prior(student_t(5,0,10), class = b) +
#'     prior(cauchy(0,2), class = sd)
#' )
#'
#' prior_summary(fit)
#' prior_summary(fit, all = FALSE)
#' print(prior_summary(fit, all = FALSE), show_df = FALSE)
#' }
#'
#' @method prior_summary brmsfit
#' @export
#' @export prior_summary
#' @importFrom rstantools prior_summary
#' @export
prior_summary.brmsfit <- function(object, all = TRUE, ...) {
  object <- restructure(object)
  prior <- object$prior
  if (!all) {
    prior <- prior[nzchar(prior$prior), ]
  }
  prior
}

#' @export
default_prior.brmsfit <- function(object, ...) {
  # just in case people try to apply default_prior to brmsfit objects
  prior_summary.brmsfit(object, ...)
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
print.brmsprior <- function(x, show_df = NULL, ...) {
  if (is.null(show_df)) {
    show_df <- nrow(x) > 1L
  }
  show_df <- as_one_logical(show_df)
  y <- prepare_print_prior(x)
  if (show_df) {
    print.data.frame(y, row.names = FALSE, ...)
  } else {
    cat(collapse(.print_prior(y), "\n"))
  }
  invisible(x)
}

# prepare pretty printing of brmsprior objects
prepare_print_prior <- function(x) {
  stopifnot(is.brmsprior(x))
  if (is.null(x$source)) {
    x$source <- ""
  }
  x$source[!nzchar(x$source)] <- "(unknown)"
  # vectorize priors and bounds for pretty printing
  # TODO: improve efficiency of adding vectorization tags
  for (i in which(!nzchar(x$prior))) {
    base_prior <- stan_base_prior(x, sel_prior = x[i, ])
    if (nzchar(base_prior)) {
      x$prior[i] <- base_prior
      x$source[i] <- "(vectorized)"
    } else {
      x$prior[i] <- "(flat)"
    }
  }
  for (i in which(!nzchar(x$lb) & !nzchar(x$ub))) {
    base_bounds <- stan_base_prior(x, c("lb", "ub"), sel_prior = x[i, ])
    x$lb[i] <- base_bounds[, "lb"]
    x$ub[i] <- base_bounds[, "ub"]
  }
  x
}

# prepare text for print.brmsprior
.print_prior <- function(x) {
  group <-  usc(x$group)
  resp <- usc(x$resp)
  dpar <- usc(x$dpar)
  nlpar <- usc(x$nlpar)
  coef <- usc(x$coef)
  if (any(nzchar(c(resp, dpar, nlpar, coef)))) {
    group <- usc(group, "suffix")
  }
  bound <- convert_bounds2stan(x[c("lb", "ub")])
  bound <- ifelse(nzchar(bound), paste0(bound, " "), "")
  tilde <- ifelse(nzchar(x$class) | nzchar(group) | nzchar(coef), " ~ ", "")
  prior <- ifelse(nzchar(x$prior), x$prior, "(flat)")
  paste0(bound, x$class, group, resp, dpar, nlpar, coef, tilde, prior)
}

# combine multiple brmsprior objects into one brmsprior
#' @export
c.brmsprior <- function(x, ..., replace = FALSE) {
  dots <- list(...)
  if (all(sapply(dots, is.brmsprior))) {
    replace <- as_one_logical(replace)
    # don't use 'c()' here to avoid creating a recursion
    out <- do_call(rbind, list(x, ...))
    if (replace) {
      # update duplicated priors
      out <- unique(out, fromLast = TRUE)
    }
  } else {
    if (length(dots)) {
      stop2("Cannot add '", class(dots[[1]])[1], "' objects to the prior.")
    }
    out <- c(as.data.frame(x))
  }
  out
}

#' @export
"+.brmsprior" <- function(e1, e2) {
  if (is.null(e2)) {
    return(e1)
  }
  if (!is.brmsprior(e2)) {
    stop2("Cannot add '", class(e2)[1], "' objects to the prior.")
  }
  c(e1, e2)
}

#' Transform into a brmsprior object
#'
#' Try to transform an object into a \code{brmsprior} object.
#'
#' @param x An object to be transformed.
#' @return A \code{brmsprior} object if the transformation was possible.
#'
#' @export
as.brmsprior <- function(x) {
  if (is.brmsprior(x)) {
    return(x)
  }
  x <- as.data.frame(x)
  if (!"prior" %in% names(x)) {
    stop2("Column 'prior' is required.")
  }
  x$prior <- as.character(x$prior)

  defaults <- c(
    class = "b", coef = "", group = "", resp = "",
    dpar = "", nlpar = "", lb = NA, ub = NA
  )
  for (v in names(defaults)) {
    if (!v %in% names(x)) {
      x[[v]] <- defaults[v]
    }
    x[[v]] <- as.character(x[[v]])
  }
  x$source <- "user"
  all_vars <- c("prior", names(defaults), "source")
  x <- x[, all_vars, drop = FALSE]
  class(x) <- c("brmsprior", "data.frame")
  x
}

#' @export
duplicated.brmsprior <- function(x, incomparables = FALSE, ...) {
  # compare only specific columns of the brmsprior object
  duplicated.data.frame(x[, rcols_prior()], incomparables, ...)
}

# evaluate the dirichlet prior of simplex parameters
# avoid name clashing with the dirichlet family
# @param prior a character vector of the form 'dirichlet(...)'
# @param len desired length of the prior concentration vector
# @param env environment in which to search for data
# @return a numeric vector of prior concentration values
eval_dirichlet <- function(prior, len = NULL, env = NULL) {
  dirichlet <- function(...) {
    out <- try(as.numeric(c(...)))
    if (is_try_error(out)) {
      stop2("Something went wrong. Did you forget to store ",
            "auxiliary data in the 'data2' argument?")
    }
    if (anyNA(out) || any(out <= 0)) {
      stop2("The dirichlet prior expects positive values.")
    }
    if (!is.null(len)) {
      if (length(out) == 1L) {
        out <- rep(out, len)
      }
      if (length(out) != len) {
        stop2("Invalid Dirichlet prior. Expected input of length ", len, ".")
      }
    }
    return(out)
  }
  prior <- as_one_character(prior)
  if (!nzchar(prior)) {
    prior <- "dirichlet(1)"
  }
  eval2(prior, envir = env, enclos = environment())
}

#' Regularized horseshoe priors in \pkg{brms}
#'
#' Function used to set up regularized horseshoe priors and related
#' hierarchical shrinkage priors for population-level effects in \pkg{brms}. The
#' function does not evaluate its arguments -- it exists purely to help set up
#' the model.
#'
#' @param df Degrees of freedom of student-t prior of the
#'   local shrinkage parameters. Defaults to \code{1}.
#' @param scale_global Scale of the student-t prior of the global shrinkage
#'   parameter. Defaults to \code{1}.
#'   In linear models, \code{scale_global} will internally be
#'   multiplied by the residual standard deviation parameter \code{sigma}.
#' @param df_global Degrees of freedom of student-t prior of the
#'   global shrinkage parameter. Defaults to \code{1}. If \code{df_global}
#'   is greater \code{1}, the shape of the prior will no longer resemble
#'   a horseshoe and it may be more appropriately called an hierarchical
#'   shrinkage prior in this case.
#' @param scale_slab Scale of the Student-t slab. Defaults to \code{2}. The
#'   original unregularized horseshoe prior is obtained by setting
#'   \code{scale_slab} to infinite, which we can approximate in practice by
#'   setting it to a very large real value.
#' @param df_slab Degrees of freedom of the student-t slab.
#'   Defaults to \code{4}.
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
#' @param main Logical (defaults to \code{FALSE}); only relevant if the horseshoe
#'   prior spans multiple parameter classes. In this case, only arguments given
#'   in the single instance where \code{main} is \code{TRUE} will be used.
#'   Arguments given in other instances of the prior will be ignored.
#'   See the Examples section below.
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
#'   Currently, the following classes support the horseshoe prior: \code{b}
#'   (overall regression coefficients), \code{sds} (SDs of smoothing splines),
#'   \code{sdgp} (SDs of Gaussian processes), \code{ar} (autoregressive
#'   coefficients), \code{ma} (moving average coefficients), \code{sderr} (SD of
#'   latent residuals), \code{sdcar} (SD of spatial CAR structures), \code{sd}
#'   (SD of varying coefficients).
#'
#' @references
#' Carvalho, C. M., Polson, N. G., & Scott, J. G. (2009). Handling sparsity via
#' the horseshoe. Artificial Intelligence and Statistics.
#' \url{http://proceedings.mlr.press/v5/carvalho09a}
#'
#' Piironen J. & Vehtari A. (2017). On the Hyperprior Choice for the Global
#' Shrinkage Parameter in the Horseshoe Prior. Artificial Intelligence and
#' Statistics. \url{https://arxiv.org/pdf/1610.05559v1.pdf}
#'
#' Piironen, J., and Vehtari, A. (2017). Sparsity information and regularization
#' in the horseshoe and other shrinkage priors. Electronic Journal of
#' Statistics. \url{https://arxiv.org/abs/1707.01694}
#'
#' @seealso \code{\link{set_prior}}
#'
#' @examples
#' set_prior(horseshoe(df = 3, par_ratio = 0.1))
#'
#' # specify the horseshoe prior across multiple parameter classes
#' set_prior(horseshoe(df = 3, par_ratio = 0.1, main = TRUE), class = "b") +
#'   set_prior(horseshoe(), class = "sd")
#'
#' @export
horseshoe <- function(df = 1, scale_global = 1, df_global = 1,
                      scale_slab = 2, df_slab = 4, par_ratio = NULL,
                      autoscale = TRUE, main = FALSE) {
  out <- deparse0(match.call())
  name <- "horseshoe"
  df <- as.numeric(df)
  df_global <- as.numeric(df_global)
  df_slab <- as.numeric(df_slab)
  scale_global <- as.numeric(scale_global)
  scale_slab <- as.numeric(scale_slab)
  main <- as_one_logical(main)
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
    if (!isTRUE(par_ratio > 0)) {
      stop2("Argument 'par_ratio' must be greater 0.")
    }
  }
  autoscale <- as_one_logical(autoscale)
  att <- nlist(
    name, df, df_global, df_slab, scale_global,
    scale_slab, par_ratio, autoscale, main
  )
  attributes(out)[names(att)] <- att
  out
}

#' R2D2 Priors in \pkg{brms}
#'
#' Function used to set up R2D2 priors for population-level effects in
#' \pkg{brms}. The function does not evaluate its arguments -- it exists purely
#' to help set up the model.
#'
#' @param mean_R2 Mean of the Beta prior on the coefficient of determination R^2.
#' @param prec_R2 Precision of the Beta prior on the coefficient of determination R^2.
#' @param cons_D2 Concentration vector of the Dirichlet prior on the variance
#'   decomposition parameters. Lower values imply more shrinkage.
#' @param autoscale Logical; indicating whether the R2D2
#'   prior should be scaled using the residual standard deviation
#'   \code{sigma} if possible and sensible (defaults to \code{TRUE}).
#'   Autoscaling is not applied for distributional parameters or
#'   when the model does not contain the parameter \code{sigma}.
#' @param main Logical (defaults to \code{FALSE}); only relevant if the R2D2
#'   prior spans multiple parameter classes. In this case, only arguments given
#'   in the single instance where \code{main} is \code{TRUE} will be used.
#'   Arguments given in other instances of the prior will be ignored.
#'   See the Examples section below.
#'
#' @details
#'   Currently, the following classes support the R2D2 prior: \code{b}
#'   (overall regression coefficients), \code{sds} (SDs of smoothing splines),
#'   \code{sdgp} (SDs of Gaussian processes), \code{ar} (autoregressive
#'   coefficients), \code{ma} (moving average coefficients), \code{sderr} (SD of
#'   latent residuals), \code{sdcar} (SD of spatial CAR structures), \code{sd}
#'   (SD of varying coefficients).
#'
#'   Even when the R2D2 prior is applied to multiple parameter classes at once,
#'   the concentration vector (argument \code{cons_D2}) has to be provided
#'   jointly in the the one instance of the prior where \code{main = TRUE}. The
#'   order in which the elements of concentration vector correspond to the
#'   classes' coefficients is the same as the order of the classes provided
#'   above.
#'
#' @references
#' Zhang, Y. D., Naughton, B. P., Bondell, H. D., & Reich, B. J. (2020).
#' Bayesian regression using a prior on the model fit: The R2-D2 shrinkage
#' prior. Journal of the American Statistical Association.
#' \url{https://arxiv.org/pdf/1609.00046.pdf}
#'
#' Aguilar J. E. & Bürkner P. C. (2022). Intuitive Joint Priors for Bayesian
#' Linear Multilevel Models: The R2D2M2 prior. ArXiv preprint.
#' \url{https://arxiv.org/pdf/2208.07132.pdf}
#'
#' @seealso \code{\link{set_prior}}
#'
#' @examples
#' set_prior(R2D2(mean_R2 = 0.8, prec_R2 = 10))
#'
#' # specify the R2D2 prior across multiple parameter classes
#' set_prior(R2D2(mean_R2 = 0.8, prec_R2 = 10, main = TRUE), class = "b") +
#'   set_prior(R2D2(), class = "sd")
#'
#' @export
R2D2 <- function(mean_R2 = 0.5, prec_R2 = 2, cons_D2 = 0.5, autoscale = TRUE,
                 main = FALSE) {
  out <- deparse0(match.call())
  name <- "R2D2"
  mean_R2 <- as_one_numeric(mean_R2)
  prec_R2 <- as_one_numeric(prec_R2)
  cons_D2 <- as.numeric(cons_D2)
  main <- as_one_logical(main)
  if (!(mean_R2 > 0 && mean_R2 < 1)) {
    stop2("Invalid R2D2 prior: Mean of the R2 prior ",
          "must be a single number in (0, 1).")
  }
  if (prec_R2 <= 0) {
    stop2("Invalid R2D2 prior: Precision of the R2 prior ",
          "must be a single positive number.")
  }
  if (any(cons_D2 <= 0)) {
    stop2("Invalid R2D2 prior: Concentration of the D2 prior ",
          "must be a vector of positive numbers.")
  }
  autoscale <- as_one_logical(autoscale)
  att <- nlist(name, mean_R2, prec_R2, cons_D2, autoscale, main)
  attributes(out)[names(att)] <- att
  out
}

#' (Defunct) Set up a lasso prior in \pkg{brms}
#'
#' This functionality is no longer supported as of brms version 2.19.2. Please
#' use the \code{\link{horseshoe}} or \code{\link{R2D2}} shrinkage priors instead.
#'
#' @param df Degrees of freedom of the chi-square prior of the inverse tuning
#'   parameter. Defaults to \code{1}.
#' @param scale Scale of the lasso prior. Defaults to \code{1}.
#'
#' @return An error indicating that the lasso prior is no longer supported.
#'
#' @references
#' Park, T., & Casella, G. (2008). The Bayesian Lasso. Journal of the American
#'    Statistical Association, 103(482), 681-686.
#'
#' @seealso \code{\link{set_prior}}, \code{\link{horseshoe}}, \code{\link{R2D2}}
#'
#' @export
lasso <- function(df = 1, scale = 1) {
  stop2("The lasso prior is no longer supported as of brms version 2.19.2. ",
        "Please use the horseshoe or R2D2 shrinkage priors instead.")
}

# check for the usage of special priors
# @param prior a character vector of priors
# @param target optional special priors to search for
#   if NULL search for all special priors
# @return a logical vector equal to the length of 'prior'
is_special_prior <- function(prior, target = NULL) {
  stopifnot(is.character(prior))
  if (is.null(target)) {
    target <- c("horseshoe", "R2D2", "lasso")
  }
  regex <- paste0("^", regex_or(target), "\\(")
  grepl(regex, prior)
}

# extract special prior information
# @param prior a brmsprior object
# @param class parameter class to be checked. the default ensures that
#.  the presence of any special prior is always detected
# @param px object from which the prefix can be extract
# @param type type of the special prior
get_special_prior <- function(prior, px, class = NULL, main = FALSE) {
  out <- attr(prior, "special")
  prefix <- combine_prefix(px, keep_mu = TRUE)
  out <- out[[prefix]]
  if (!length(out)) {
    return(NULL)
  }
  if (main) {
    # get the main special prior to extract arguments from
    if (length(out) == 1L) {
      # only one class present which must then be main
      out <- out[[1]]
    } else {
      main <- which(ufrom_list(out, "main"))
      if (length(main) != 1L) {
        stop2("If special priors for multiple classes are given, ",
              "exactly one of them must be marked with 'main = TRUE'.")
      }
      out <- out[[main]]
    }
  } else if (!is.null(class)) {
    out <- out[[class]]
  } else {
    # just extract info on any class for example the first
    out <- out[[1]]
  }
  out
}

# is special prior information present?
has_special_prior <- function(prior, px = NULL, class = NULL) {
  if (is.null(px)) {
    # is any special prior present?
    return(length(rmNULL(attr(prior, "special"))) > 0L)
  }
  .has_special_prior <- function(px) {
    !is.null(get_special_prior(prior, px = px, class = class))
  }
  if (is.data.frame(px)) {
    # this case occurs for group-level parameters
    out <- FALSE
    for (i in seq_rows(px)) {
      out <- out || .has_special_prior(px[i, ])
    }
  } else {
    out <- .has_special_prior(px)
  }
  out
}

#' Constant priors in \pkg{brms}
#'
#' Function used to set up constant priors in \pkg{brms}.
#' The function does not evaluate its arguments -- it exists purely
#' to help set up the model.
#'
#' @param const Numeric value, vector, matrix of values to which the parameters
#'   should be fixed to. Can also be a valid Stan variable in the model.
#' @param broadcast Should \code{const} be automatically broadcasted to the
#'   correct size of the parameter? Defaults to \code{TRUE}. If you supply
#'   vectors or matrices in \code{const} or vector/matrix valued Stan variables,
#'   you need to set \code{broadcast} to \code{TRUE} (see Examples).
#'
#' @returns A named list with elements \code{const} and \code{broadcast}.
#'
#' @examples
#' stancode(count ~ Base + Age, data = epilepsy,
#'          prior = prior(constant(1), class = "b"))
#'
#' # will fail parsing because brms will try to broadcast a vector into a vector
#' stancode(count ~ Base + Age, data = epilepsy,
#'          prior = prior(constant(alpha), class = "b"),
#'          stanvars = stanvar(c(1, 0), name = "alpha"))
#'
#' stancode(count ~ Base + Age, data = epilepsy,
#'          prior = prior(constant(alpha, broadcast = FALSE), class = "b"),
#'          stanvars = stanvar(c(1, 0), name = "alpha"))
#'
#' @seealso \code{\link{set_prior}}
#'
#' @export
constant <- function(const, broadcast = TRUE) {
  const <- deparse0(substitute(const))
  const <- rename(const, c("\"", "'"), c("", ""))
  broadcast <- as_one_logical(broadcast)
  nlist(const, broadcast)
}

# check if parameters should be sampled only from the prior
is_prior_only <- function(prior) {
  is_equal(get_sample_prior(prior), "only")
}
