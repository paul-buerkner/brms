#' Set up a model formula for use in \pkg{brms}
#'
#' Set up a model formula for use in the \pkg{brms} package
#' allowing to define (potentially non-linear) additive multilevel
#' models for all parameters of the assumed response distribution.
#'
#' @aliases bf
#'
#' @param formula An object of class \code{formula}
#'   (or one that can be coerced to that class):
#'   a symbolic description of the model to be fitted.
#'   The details of model specification are given in 'Details'.
#' @param ... Additional \code{formula} objects to specify predictors of
#'   non-linear and distributional parameters. Formulas can either be named
#'   directly or contain names on their left-hand side. Alternatively,
#'   it is possible to fix parameters to certain values by passing
#'   numbers or character strings in which case arguments have to be named
#'   to provide the parameter names. See 'Details' for more information.
#' @param flist Optional list of formulas, which are treated in the
#'   same way as formulas passed via the \code{...} argument.
#' @param nl Logical; Indicates whether \code{formula} should be
#'   treated as specifying a non-linear model. By default, \code{formula}
#'   is treated as an ordinary linear model formula.
#' @param loop Logical; Only used in non-linear models.
#'   Indicates if the computation of the non-linear formula should be
#'   done inside (\code{TRUE}) or outside (\code{FALSE}) a loop
#'   over observations. Defaults to \code{TRUE}.
#' @param center Logical; Indicates if the population-level design
#'   matrix should be centered, which usually increases sampling efficiency.
#'   See the 'Details' section for more information.
#'   Defaults to \code{TRUE} for distributional parameters
#'   and to \code{FALSE} for non-linear parameters.
#' @param cmc Logical; Indicates whether automatic cell-mean coding
#'   should be enabled when removing the intercept by adding \code{0}
#'   to the right-hand of model formulas. Defaults to \code{TRUE} to
#'   mirror the behavior of standard \R formula parsing.
#' @param sparse Logical; indicates whether the population-level design matrices
#'   should be treated as sparse (defaults to \code{FALSE}). For design matrices
#'   with many zeros, this can considerably reduce required memory. Sampling
#'   speed is currently not improved or even slightly decreased.
#' @param decomp Optional name of the decomposition used for the
#'   population-level design matrix. Defaults to \code{NULL} that is
#'   no decomposition. Other options currently available are
#'   \code{"QR"} for the QR decomposition that helps in fitting models
#'   with highly correlated predictors.
#' @param family Same argument as in \code{\link{brm}}.
#'   If \code{family} is specified in \code{brmsformula}, it will
#'   overwrite the value specified in other functions.
#' @param autocor An optional \code{formula} which contains
#'   autocorrelation terms as described in \code{\link{autocor-terms}}
#'   or alternatively a \code{\link{cor_brms}} object (deprecated).
#'   If \code{autocor} is specified in \code{brmsformula}, it will
#'   overwrite the value specified in other functions.
#' @param unused An optional \code{formula} which contains variables
#'   that are unused in the model but should still be stored in the
#'   model's data frame. This can be useful, for example,
#'   if those variables are required for post-processing the model.
#'
#' @return An object of class \code{brmsformula}, which
#'   is essentially a \code{list} containing all model
#'   formulas as well as some additional information.
#'
#' @seealso \code{\link{mvbrmsformula}}, \code{\link{brmsformula-helpers}}
#'
#' @details
#'
#'   \bold{General formula structure}
#'
#'   The \code{formula} argument accepts formulas of the following syntax:
#'
#'   \code{response | aterms ~ pterms + (gterms | group)}
#'
#'   The \code{pterms} part contains effects that are assumed to be the same
#'   across observations. We call them 'population-level' or 'overall' effects,
#'   or (adopting frequentist vocabulary) 'fixed' effects. The optional
#'   \code{gterms} part may contain effects that are assumed to vary across
#'   grouping variables specified in \code{group}. We call them 'group-level' or
#'   'varying' effects, or (adopting frequentist vocabulary) 'random' effects,
#'   although the latter name is misleading in a Bayesian context. For more
#'   details type \code{vignette("brms_overview")} and
#'   \code{vignette("brms_multilevel")}.
#'
#'   \bold{Group-level terms}
#'
#'   Multiple grouping factors each with multiple group-level effects are
#'   possible. (Of course we can also run models without any group-level
#'   effects.) Instead of \code{|} you may use \code{||} in grouping terms to
#'   prevent correlations from being modeled. Equivalently, the \code{cor}
#'   argument of the \code{\link{gr}} function can be used for this purpose,
#'   for example, \code{(1 + x || g)} is equivalent to
#'   \code{(1 + x | gr(g, cor = FALSE))}.
#'
#'   It is also possible to model different group-level terms of the same
#'   grouping factor as correlated (even across different formulas, e.g., in
#'   non-linear models) by using \code{|<ID>|} instead of \code{|}. All
#'   group-level terms sharing the same ID will be modeled as correlated. If,
#'   for instance, one specifies the terms \code{(1+x|i|g)} and \code{(1+z|i|g)}
#'   somewhere in the formulas passed to \code{brmsformula}, correlations
#'   between the corresponding group-level effects will be estimated. In the
#'   above example, \code{i} is not a variable in the data but just a symbol to
#'   indicate correlations between multiple group-level terms. Equivalently, the
#'   \code{id} argument of the \code{\link{gr}} function can be used as well,
#'   for example, \code{(1 + x | gr(g, id = "i"))}.
#'
#'   If levels of the grouping factor belong to different sub-populations,
#'   it may be reasonable to assume a different covariance matrix for each
#'   of the sub-populations. For instance, the variation within the
#'   treatment group and within the control group in a randomized control
#'   trial might differ. Suppose that \code{y} is the outcome, and
#'   \code{x} is the factor indicating the treatment and control group.
#'   Then, we could estimate different hyper-parameters of the varying
#'   effects (in this case a varying intercept) for treatment and control
#'   group via \code{y ~ x + (1 | gr(subject, by = x))}.
#'
#'   You can specify multi-membership terms using the \code{\link{mm}}
#'   function. For instance, a multi-membership term with two members
#'   could be \code{(1 | mm(g1, g2))}, where \code{g1} and \code{g2}
#'   specify the first and second member, respectively. Moreover,
#'   if a covariate \code{x} varies across the levels of the grouping-factors
#'   \code{g1} and \code{g2}, we can save the respective covariate values
#'   in the variables \code{x1} and \code{x2} and then model the varying
#'   effect as \code{(1 + mmc(x1, x2) | mm(g1, g2))}.
#'
#'   \bold{Special predictor terms}
#'
#'   Flexible non-linear smooth terms can modeled using the \code{\link{s}}
#'   and \code{\link{t2}} functions in the \code{pterms} part
#'   of the model formula. This allows to fit generalized additive mixed
#'   models (GAMMs) with \pkg{brms}. The implementation is similar to that
#'   used in the \pkg{gamm4} package. For more details on this model class
#'   see \code{\link[mgcv:gam]{gam}} and \code{\link[mgcv:gamm]{gamm}}.
#'
#'   Gaussian process terms can be fitted using the \code{\link{gp}}
#'   function in the \code{pterms} part of the model formula. Similar to
#'   smooth terms, Gaussian processes can be used to model complex non-linear
#'   relationships, for instance temporal or spatial autocorrelation.
#'   However, they are computationally demanding and are thus not recommended
#'   for very large datasets or approximations need to be used.
#'
#'   The \code{pterms} and \code{gterms} parts may contain four non-standard
#'   effect types namely monotonic, measurement error, missing value, and
#'   category specific effects, which can be specified using terms of the
#'   form \code{mo(predictor)}, \code{me(predictor, sd_predictor)},
#'   \code{mi(predictor)}, and \code{cs(<predictors>)}, respectively.
#'   Category specific effects can only be estimated in
#'   ordinal models and are explained in more detail in the package's
#'   main vignette (type \code{vignette("brms_overview")}).
#'   The other three effect types are explained in the following.
#'
#'   A monotonic predictor must either be integer valued or an ordered factor,
#'   which is the first difference to an ordinary continuous predictor.
#'   More importantly, predictor categories (or integers) are not assumed to be
#'   equidistant with respect to their effect on the response variable.
#'   Instead, the distance between adjacent predictor categories (or integers)
#'   is estimated from the data and may vary across categories.
#'   This is realized by parameterizing as follows:
#'   One parameter takes care of the direction and size of the effect similar
#'   to an ordinary regression parameter, while an additional parameter vector
#'   estimates the normalized distances between consecutive predictor categories.
#'   A main application of monotonic effects are ordinal predictors that
#'   can this way be modeled without (falsely) treating them as continuous
#'   or as unordered categorical predictors. For more details and examples
#'   see \code{vignette("brms_monotonic")}.
#'
#'   Quite often, predictors are measured and as such naturally contain
#'   measurement error. Although most researchers are well aware of this problem,
#'   measurement error in predictors is ignored in most
#'   regression analyses, possibly because only few packages allow
#'   for modeling it. Notably, measurement error can be handled in
#'   structural equation models, but many more general regression models
#'   (such as those featured by \pkg{brms}) cannot be transferred
#'   to the SEM framework. In \pkg{brms}, effects of noise-free predictors
#'   can be modeled using the \code{me} (for 'measurement error') function.
#'   If, say, \code{y} is the response variable and
#'   \code{x} is a measured predictor with known measurement error
#'   \code{sdx}, we can simply include it on the right-hand side of the
#'   model formula via \code{y ~ me(x, sdx)}.
#'   This can easily be extended to more general formulas.
#'   If \code{x2} is another measured predictor with corresponding error
#'   \code{sdx2} and \code{z} is a predictor without error
#'   (e.g., an experimental setting), we can model all main effects
#'   and interactions of the three predictors in the well known manner:
#'   \code{y ~ me(x, sdx) * me(x2, sdx2) * z}.
#'   The \code{me} function is soft deprecated in favor of the more flexible
#'   and consistent \code{mi} function (see below).
#'
#'   When a variable contains missing values, the corresponding rows will
#'   be excluded from the data by default (row-wise exclusion). However,
#'   quite often we want to keep these rows and instead estimate the missing values.
#'   There are two approaches for this: (a) Impute missing values before
#'   the model fitting for instance via multiple imputation (see
#'   \code{\link{brm_multiple}} for a way to handle multiple imputed datasets).
#'   (b) Impute missing values on the fly during model fitting. The latter
#'   approach is explained in the following. Using a variable with missing
#'   values as predictors requires two things, First, we need to specify that
#'   the predictor contains missings that should to be imputed.
#'   If, say, \code{y} is the primary response, \code{x} is a
#'   predictor with missings and \code{z} is a predictor without missings,
#'   we go for \code{y ~ mi(x) + z}. Second, we need to model \code{x}
#'   as an additional response with corresponding predictors and the
#'   addition term \code{mi()}. In our example, we could write
#'   \code{x | mi() ~ z}. Measurement error may be included via
#'   the \code{sdy} argument, say, \code{x | mi(sdy = se) ~ z}.
#'   See \code{\link{mi}} for examples with real data.
#'
#'
#'   \bold{Autocorrelation terms}
#'
#'   Autocorrelation terms can be directly specified inside the \code{pterms}
#'   part as well. Details can be found in \code{\link{autocor-terms}}.
#'
#'   \bold{Additional response information}
#'
#'   Another special of the \pkg{brms} formula syntax is the optional
#'   \code{aterms} part, which may contain multiple terms of the form
#'   \code{fun(<variable>)} separated by \code{+} each providing special
#'   information on the response variable. \code{fun} can be replaced with
#'   either \code{se}, \code{weights}, \code{subset}, \code{cens}, \code{trunc},
#'   \code{trials}, \code{cat}, \code{dec}, \code{rate}, \code{vreal}, or
#'   \code{vint}. Their meanings are explained below
#'   (see also \code{\link{addition-terms}}).
#'
#'   For families \code{gaussian}, \code{student} and \code{skew_normal}, it is
#'   possible to specify standard errors of the observations, thus allowing
#'   to perform meta-analysis. Suppose that the variable \code{yi} contains
#'   the effect sizes from the studies and \code{sei} the corresponding
#'   standard errors. Then, fixed and random effects meta-analyses can
#'   be conducted using the formulas \code{yi | se(sei) ~ 1} and
#'   \code{yi | se(sei) ~ 1 + (1|study)}, respectively, where
#'   \code{study} is a variable uniquely identifying every study.
#'   If desired, meta-regression can be performed via
#'   \code{yi | se(sei) ~ 1 + mod1 + mod2 + (1|study)}
#'   or \cr \code{yi | se(sei) ~ 1 + mod1 + mod2 + (1 + mod1 + mod2|study)},
#'   where \code{mod1} and \code{mod2} represent moderator variables.
#'   By default, the standard errors replace the parameter \code{sigma}.
#'   To model \code{sigma} in addition to the known standard errors,
#'   set argument \code{sigma} in function \code{se} to \code{TRUE},
#'   for instance, \code{yi | se(sei, sigma = TRUE) ~ 1}.
#'
#'   For all families, weighted regression may be performed using
#'   \code{weights} in the \code{aterms} part. Internally, this is
#'   implemented by multiplying the log-posterior values of each
#'   observation by their corresponding weights.
#'   Suppose that variable \code{wei} contains the weights
#'   and that \code{yi} is the response variable.
#'   Then, formula \code{yi | weights(wei) ~ predictors}
#'   implements a weighted regression.
#'
#'   For multivariate models, \code{subset} may be used in the \code{aterms}
#'   part, to use different subsets of the data in different univariate
#'   models. For instance, if \code{sub} is a logical variable and
#'   \code{y} is the response of one of the univariate models, we may
#'   write \code{y | subset(sub) ~ predictors} so that \code{y} is
#'   predicted only for those observations for which \code{sub} evaluates
#'   to \code{TRUE}.
#'
#'   For log-linear models such as poisson models, \code{rate} may be used
#'   in the \code{aterms} part to specify the denominator of a response that
#'   is expressed as a rate. The numerator is given by the actual response
#'   variable and has a distribution according to the family as usual. Using
#'   \code{rate(denom)} is equivalent to adding \code{offset(log(denom))} to
#'   the linear predictor of the main parameter but the former is arguably
#'   more convenient and explicit.
#'
#'   With the exception of categorical and ordinal families,
#'   left, right, and interval censoring can be modeled through
#'   \code{y | cens(censored) ~ predictors}. The censoring variable
#'   (named \code{censored} in this example) should contain the values
#'   \code{'left'}, \code{'none'}, \code{'right'}, and \code{'interval'}
#'   (or equivalently \code{-1}, \code{0}, \code{1}, and \code{2}) to indicate that
#'   the corresponding observation is left censored, not censored, right censored,
#'   or interval censored. For interval censored data, a second variable
#'   (let's call it \code{y2}) has to be passed to \code{cens}. In this case,
#'   the formula has the structure \code{y | cens(censored, y2) ~ predictors}.
#'   While the lower bounds are given in \code{y}, the upper bounds are given
#'   in \code{y2} for interval censored data. Intervals are assumed to be open
#'   on the left and closed on the right: \code{(y, y2]}.
#'
#'   With the exception of categorical and ordinal families,
#'   the response distribution can be truncated using the \code{trunc}
#'   function in the addition part. If the response variable is truncated
#'   between, say, 0 and 100, we can specify this via
#'   \code{yi | trunc(lb = 0, ub = 100) ~ predictors}.
#'   Instead of numbers, variables in the data set can also be passed allowing
#'   for varying truncation points across observations. Defining only one of
#'   the two arguments in \code{trunc} leads to one-sided truncation.
#'
#'   For all continuous families, missing values in the responses can be imputed
#'   within Stan by using the addition term \code{mi}. This is mostly
#'   useful in combination with \code{mi} predictor terms as explained
#'   above under 'Special predictor terms'.
#'
#'   For families \code{binomial} and \code{zero_inflated_binomial},
#'   addition should contain a variable indicating the number of trials
#'   underlying each observation. In \code{lme4} syntax, we may write for instance
#'   \code{cbind(success, n - success)}, which is equivalent
#'   to \code{success | trials(n)} in \pkg{brms} syntax. If the number of trials
#'   is constant across all observations, say \code{10},
#'   we may also write \code{success | trials(10)}.
#'   \bold{Please note that the \code{cbind()} syntax will not work
#'   in \pkg{brms} in the expected way because this syntax is reserved
#'   for other purposes.}
#'
#'   For all ordinal families, \code{aterms} may contain a term
#'   \code{thres(number)} to specify the number thresholds (e.g,
#'   \code{thres(6)}), which should be equal to the total number of response
#'   categories - 1. If not given, the number of thresholds is calculated from
#'   the data. If different threshold vectors should be used for different
#'   subsets of the data, the \code{gr} argument can be used to provide the
#'   grouping variable (e.g, \code{thres(6, gr = item)}, if \code{item} is the
#'   grouping variable). In this case, the number of thresholds can also be a
#'   variable in the data with different values per group.
#'
#'   A deprecated quasi alias of \code{thres()} is \code{cat()} with which the
#'   total number of response categories (i.e., number of thresholds + 1) can be
#'   specified.
#'
#'   In Wiener diffusion models (family \code{wiener}) the addition term
#'   \code{dec} is mandatory to specify the (vector of) binary decisions
#'   corresponding to the reaction times. Non-zero values will be treated
#'   as a response on the upper boundary of the diffusion process and zeros
#'   will be treated as a response on the lower boundary. Alternatively,
#'   the variable passed to \code{dec} might also be a character vector
#'   consisting of \code{'lower'} and \code{'upper'}.
#'
#'   All families support the \code{index} addition term to uniquely identify
#'   each observation of the corresponding response variable. Currently,
#'   \code{index} is primarily useful in combination with the \code{subset}
#'   addition and \code{\link{mi}} terms.
#'
#'   For custom families, it is possible to pass an arbitrary number of real and
#'   integer vectors via the addition terms \code{vreal} and \code{vint},
#'   respectively. An example is provided in
#'   \code{vignette('brms_customfamilies')}. To pass multiple vectors of the
#'   same data type, provide them separated by commas inside a single
#'   \code{vreal} or \code{vint} statement.
#'
#'   Multiple addition terms of different types may be specified at the same
#'   time using the \code{+} operator. For example, the formula
#'   \code{formula = yi | se(sei) + cens(censored) ~ 1} implies a censored
#'   meta-analytic model.
#'
#'   The addition argument \code{disp} (short for dispersion)
#'   has been removed in version 2.0. You may instead use the
#'   distributional regression approach by specifying
#'   \code{sigma ~ 1 + offset(log(xdisp))} or
#'   \code{shape ~ 1 + offset(log(xdisp))}, where \code{xdisp} is
#'   the variable being previously passed to \code{disp}.
#'
#'   \bold{Parameterization of the population-level intercept}
#'
#'   By default, the population-level intercept (if incorporated) is estimated
#'   separately and not as part of population-level parameter vector \code{b} As
#'   a result, priors on the intercept also have to be specified separately.
#'   Furthermore, to increase sampling efficiency, the population-level design
#'   matrix \code{X} is centered around its column means \code{X_means} if the
#'   intercept is incorporated. This leads to a temporary bias in the intercept
#'   equal to \code{<X_means, b>}, where \code{<,>} is the scalar product. The
#'   bias is corrected after fitting the model, but be aware that you are
#'   effectively defining a prior on the intercept of the centered design matrix
#'   not on the real intercept. You can turn off this special handling of the
#'   intercept by setting argument \code{center} to \code{FALSE}. For more
#'   details on setting priors on population-level intercepts, see
#'   \code{\link{set_prior}}.
#'
#'   This behavior can be avoided by using the reserved
#'   (and internally generated) variable \code{Intercept}.
#'   Instead of \code{y ~ x}, you may write
#'   \code{y ~ 0 + Intercept + x}. This way, priors can be
#'   defined on the real intercept, directly. In addition,
#'   the intercept is just treated as an ordinary population-level effect
#'   and thus priors defined on \code{b} will also apply to it.
#'   Note that this parameterization may be less efficient
#'   than the default parameterization discussed above.
#'
#'   \bold{Formula syntax for non-linear models}
#'
#'   In \pkg{brms}, it is possible to specify non-linear models
#'   of arbitrary complexity.
#'   The non-linear model can just be specified within the \code{formula}
#'   argument. Suppose, that we want to predict the response \code{y}
#'   through the predictor \code{x}, where \code{x} is linked to \code{y}
#'   through \code{y = alpha - beta * lambda^x}, with parameters
#'   \code{alpha}, \code{beta}, and \code{lambda}. This is certainly a
#'   non-linear model being defined via
#'   \code{formula = y ~ alpha - beta * lambda^x} (addition arguments
#'   can be added in the same way as for ordinary formulas).
#'   To tell \pkg{brms} that this is a non-linear model,
#'   we set argument \code{nl} to \code{TRUE}.
#'   Now we have to specify a model for each of the non-linear parameters.
#'   Let's say we just want to estimate those three parameters
#'   with no further covariates or random effects. Then we can pass
#'   \code{alpha + beta + lambda ~ 1} or equivalently
#'   (and more flexible) \code{alpha ~ 1, beta ~ 1, lambda ~ 1}
#'   to the \code{...} argument.
#'   This can, of course, be extended. If we have another predictor \code{z} and
#'   observations nested within the grouping factor \code{g}, we may write for
#'   instance \code{alpha ~ 1, beta ~ 1 + z + (1|g), lambda ~ 1}.
#'   The formula syntax described above applies here as well.
#'   In this example, we are using \code{z} and \code{g} only for the
#'   prediction of \code{beta}, but we might also use them for the other
#'   non-linear parameters (provided that the resulting model is still
#'   scientifically reasonable).
#'
#'   By default, non-linear covariates are treated as real vectors in Stan.
#'   However, if the data of the covariates is of type `integer` in \R (which
#'   can be enforced by the `as.integer` function), the Stan type will be
#'   changed to an integer array. That way, covariates can also be used
#'   for indexing purposes in Stan.
#'
#'   Non-linear models may not be uniquely identified and / or show bad convergence.
#'   For this reason it is mandatory to specify priors on the non-linear parameters.
#'   For instructions on how to do that, see \code{\link{set_prior}}.
#'   For some examples of non-linear models, see \code{vignette("brms_nonlinear")}.
#'
#'   \bold{Formula syntax for predicting distributional parameters}
#'
#'   It is also possible to predict parameters of the response distribution such
#'   as the residual standard deviation \code{sigma} in gaussian models or the
#'   hurdle probability \code{hu} in hurdle models. The syntax closely resembles
#'   that of a non-linear parameter, for instance \code{sigma ~ x + s(z) +
#'   (1+x|g)}. For some examples of distributional models, see
#'   \code{vignette("brms_distreg")}.
#'
#'   Parameter \code{mu} exists for every family and can be used as an
#'   alternative to specifying terms in \code{formula}. If both \code{mu} and
#'   \code{formula} are given, the right-hand side of \code{formula} is ignored.
#'   Accordingly, specifying terms on the right-hand side of both \code{formula}
#'   and \code{mu} at the same time is deprecated. In future versions,
#'   \code{formula} might be updated by \code{mu}.
#'
#'   The following are
#'   distributional parameters of specific families (all other parameters are
#'   treated as non-linear parameters): \code{sigma} (residual standard
#'   deviation or scale of the \code{gaussian}, \code{student},
#'   \code{skew_normal}, \code{lognormal} \code{exgaussian}, and
#'   \code{asym_laplace} families); \code{shape} (shape parameter of the
#'   \code{Gamma}, \code{weibull}, \code{negbinomial}, and related zero-inflated
#'   / hurdle families); \code{nu} (degrees of freedom parameter of the
#'   \code{student} and \code{frechet} families); \code{phi} (precision
#'   parameter of the \code{beta} and \code{zero_inflated_beta} families);
#'   \code{kappa} (precision parameter of the \code{von_mises} family);
#'   \code{beta} (mean parameter of the exponential component of the
#'   \code{exgaussian} family); \code{quantile} (quantile parameter of the
#'   \code{asym_laplace} family); \code{zi} (zero-inflation probability);
#'   \code{hu} (hurdle probability); \code{zoi} (zero-one-inflation
#'   probability); \code{coi} (conditional one-inflation probability);
#'   \code{disc} (discrimination) for ordinal models; \code{bs}, \code{ndt}, and
#'   \code{bias} (boundary separation, non-decision time, and initial bias of
#'   the \code{wiener} diffusion model). By default, distributional parameters
#'   are modeled on the log scale if they can be positive only or on the logit
#'   scale if the can only be within the unit interval.
#'
#'   Alternatively, one may fix distributional parameters to certain values.
#'   However, this is mainly useful when models become too
#'   complicated and otherwise have convergence issues.
#'   We thus suggest to be generally careful when making use of this option.
#'   The \code{quantile} parameter of the \code{asym_laplace} distribution
#'   is a good example where it is useful. By fixing \code{quantile},
#'   one can perform quantile regression for the specified quantile.
#'   For instance, \code{quantile = 0.25} allows predicting the 25\%-quantile.
#'   Furthermore, the \code{bias} parameter in drift-diffusion models,
#'   is assumed to be \code{0.5} (i.e. no bias) in many applications.
#'   To achieve this, simply write \code{bias = 0.5}.
#'   Other possible applications are the Cauchy distribution as a
#'   special case of the Student-t distribution with
#'   \code{nu = 1}, or the geometric distribution as a special case of
#'   the negative binomial distribution with \code{shape = 1}.
#'   Furthermore, the parameter \code{disc} ('discrimination') in ordinal
#'   models is fixed to \code{1} by default and not estimated,
#'   but may be modeled as any other distributional parameter if desired
#'   (see examples). For reasons of identification, \code{'disc'}
#'   can only be positive, which is achieved by applying the log-link.
#'
#'   In categorical models, distributional parameters do not have
#'   fixed names. Instead, they are named after the response categories
#'   (excluding the first one, which serves as the reference category),
#'   with the prefix \code{'mu'}. If, for instance, categories are named
#'   \code{cat1}, \code{cat2}, and \code{cat3}, the distributional parameters
#'   will be named \code{mucat2} and \code{mucat3}.
#'
#'   Some distributional parameters currently supported by \code{brmsformula}
#'   have to be positive (a negative standard deviation or precision parameter
#'   does not make any sense) or are bounded between 0 and 1 (for zero-inflated /
#'   hurdle probabilities, quantiles, or the initial bias parameter of
#'   drift-diffusion models).
#'   However, linear predictors can be positive or negative, and thus the log link
#'   (for positive parameters) or logit link (for probability parameters) are used
#'   by default to ensure that distributional parameters are within their valid intervals.
#'   This implies that, by default, effects for such distributional parameters are
#'   estimated on the log / logit scale and one has to apply the inverse link
#'   function to get to the effects on the original scale.
#'   Alternatively, it is possible to use the identity link to predict parameters
#'   on their original scale, directly. However, this is much more likely to lead
#'   to problems in the model fitting, if the parameter actually has a restricted range.
#'
#'   See also \code{\link{brmsfamily}} for an overview of valid link functions.
#'
#'   \bold{Formula syntax for mixture models}
#'
#'   The specification of mixture models closely resembles that
#'   of non-mixture models. If not specified otherwise (see below),
#'   all mean parameters of the mixture components are predicted
#'   using the right-hand side of \code{formula}. All types of predictor
#'   terms allowed in non-mixture models are allowed in mixture models
#'   as well.
#'
#'   Distributional parameters of mixture distributions have the same
#'   name as those of the corresponding ordinary distributions, but with
#'   a number at the end to indicate the mixture component. For instance, if
#'   you use family \code{mixture(gaussian, gaussian)}, the distributional
#'   parameters are \code{sigma1} and \code{sigma2}.
#'   Distributional parameters of the same class can be fixed to the same value.
#'   For the above example, we could write \code{sigma2 = "sigma1"} to make
#'   sure that both components have the same residual standard deviation,
#'   which is in turn estimated from the data.
#'
#'   In addition, there are two types of special distributional parameters.
#'   The first are named \code{mu<ID>}, that allow for modeling different
#'   predictors for the mean parameters of different mixture components.
#'   For instance, if you want to predict the mean of the first component
#'   using predictor \code{x} and the mean of the second component using
#'   predictor \code{z}, you can write \code{mu1 ~ x} as well as \code{mu2 ~ z}.
#'   The second are named \code{theta<ID>}, which constitute the mixing
#'   proportions. If the mixing proportions are fixed to certain values,
#'   they are internally normalized to form a probability vector.
#'   If one seeks to predict the mixing proportions, all but
#'   one of the them has to be predicted, while the remaining one is used
#'   as the reference category to identify the model. The so-called 'softmax'
#'   transformation is applied on the linear predictor terms to form a
#'   probability vector.
#'
#'   For more information on mixture models, see
#'   the documentation of \code{\link{mixture}}.
#'
#'   \bold{Formula syntax for multivariate models}
#'
#'   Multivariate models may be specified using \code{mvbind} notation
#'   or with help of the \code{\link{mvbf}} function.
#'   Suppose that \code{y1} and \code{y2} are response variables
#'   and \code{x} is a predictor. Then \code{mvbind(y1, y2) ~ x}
#'   specifies a multivariate model.
#'   The effects of all terms specified at the RHS of the formula
#'   are assumed to vary across response variables.
#'   For instance, two parameters will be estimated for \code{x},
#'   one for the effect on \code{y1} and another for the effect on \code{y2}.
#'   This is also true for group-level effects. When writing, for instance,
#'   \code{mvbind(y1, y2) ~ x + (1+x|g)}, group-level effects will be
#'   estimated separately for each response. To model these effects
#'   as correlated across responses, use the ID syntax (see above).
#'   For the present example, this would look as follows:
#'   \code{mvbind(y1, y2) ~ x + (1+x|2|g)}. Of course, you could also use
#'   any value other than \code{2} as ID.
#'
#'   It is also possible to specify different formulas for different responses.
#'   If, for instance, \code{y1} should be predicted by \code{x} and \code{y2}
#'   should be predicted by \code{z}, we could write \code{mvbf(y1 ~ x, y2 ~ z)}.
#'   Alternatively, multiple \code{brmsformula} objects can be added to
#'   specify a joint multivariate model (see 'Examples').
#'
#' @examples
#' # multilevel model with smoothing terms
#' brmsformula(y ~ x1*x2 + s(z) + (1+x1|1) + (1|g2))
#'
#' # additionally predict 'sigma'
#' brmsformula(y ~ x1*x2 + s(z) + (1+x1|1) + (1|g2),
#'             sigma ~ x1 + (1|g2))
#'
#' # use the shorter alias 'bf'
#' (formula1 <- brmsformula(y ~ x + (x|g)))
#' (formula2 <- bf(y ~ x + (x|g)))
#' # will be TRUE
#' identical(formula1, formula2)
#'
#' # incorporate censoring
#' bf(y | cens(censor_variable) ~ predictors)
#'
#' # define a simple non-linear model
#' bf(y ~ a1 - a2^x, a1 + a2 ~ 1, nl = TRUE)
#'
#' # predict a1 and a2 differently
#' bf(y ~ a1 - a2^x, a1 ~ 1, a2 ~ x + (x|g), nl = TRUE)
#'
#' # correlated group-level effects across parameters
#' bf(y ~ a1 - a2^x, a1 ~ 1 + (1 |2| g), a2 ~ x + (x |2| g), nl = TRUE)
#' # alternative but equivalent way to specify the above model
#' bf(y ~ a1 - a2^x, a1 ~ 1 + (1 | gr(g, id = 2)),
#'    a2 ~ x + (x | gr(g, id = 2)), nl = TRUE)
#'
#' # define a multivariate model
#' bf(mvbind(y1, y2) ~ x * z + (1|g))
#'
#' # define a zero-inflated model
#' # also predicting the zero-inflation part
#' bf(y ~ x * z + (1+x|ID1|g), zi ~ x + (1|ID1|g))
#'
#' # specify a predictor as monotonic
#' bf(y ~ mo(x) + more_predictors)
#'
#' # for ordinal models only
#' # specify a predictor as category specific
#' bf(y ~ cs(x) + more_predictors)
#' # add a category specific group-level intercept
#' bf(y ~ cs(x) + (cs(1)|g))
#' # specify parameter 'disc'
#' bf(y ~ person + item, disc ~ item)
#'
#' # specify variables containing measurement error
#' bf(y ~ me(x, sdx))
#'
#' # specify predictors on all parameters of the wiener diffusion model
#' # the main formula models the drift rate 'delta'
#' bf(rt | dec(decision) ~ x, bs ~ x, ndt ~ x, bias ~ x)
#'
#' # fix the bias parameter to 0.5
#' bf(rt | dec(decision) ~ x, bias = 0.5)
#'
#' # specify different predictors for different mixture components
#' mix <- mixture(gaussian, gaussian)
#' bf(y ~ 1, mu1 ~ x, mu2 ~ z, family = mix)
#'
#' # fix both residual standard deviations to the same value
#' bf(y ~ x, sigma2 = "sigma1", family = mix)
#'
#' # use the '+' operator to specify models
#' bf(y ~ 1) +
#'   nlf(sigma ~ a * exp(b * x), a ~ x) +
#'   lf(b ~ z + (1|g), dpar = "sigma") +
#'   gaussian()
#'
#' # specify a multivariate model using the '+' operator
#' bf(y1 ~ x + (1|g)) +
#'   gaussian() + cor_ar(~1|g) +
#'   bf(y2 ~ z) + poisson()
#'
#' # specify correlated residuals of a gaussian and a poisson model
#' form1 <- bf(y1 ~ 1 + x + (1|c|obs), sigma = 1) + gaussian()
#' form2 <- bf(y2 ~ 1 + x + (1|c|obs)) + poisson()
#'
#' # model missing values in predictors
#' bf(bmi ~ age * mi(chl)) +
#'   bf(chl | mi() ~ age) +
#'   set_rescor(FALSE)
#'
#' # model sigma as a function of the mean
#' bf(y ~ eta, nl = TRUE) +
#'   lf(eta ~ 1 + x) +
#'   nlf(sigma ~ tau * sqrt(eta)) +
#'   lf(tau ~ 1)
#'
#' @export
brmsformula <- function(formula, ..., flist = NULL, family = NULL,
                        autocor = NULL, nl = NULL, loop = NULL,
                        center = NULL, cmc = NULL, sparse = NULL,
                        decomp = NULL, unused = NULL) {
  if (is.brmsformula(formula)) {
    out <- formula
  } else {
    out <- list(formula = as_formula(formula))
    class(out) <- "brmsformula"
  }
  # parse and validate dots arguments
  dots <- c(out$pforms, out$pfix, list(...), flist)
  dots <- lapply(dots, function(x) if (is.list(x)) x else list(x))
  dots <- unlist(dots, recursive = FALSE)
  forms <- list()
  for (i in seq_along(dots)) {
    c(forms) <- validate_par_formula(dots[[i]], par = names(dots)[i])
  }
  is_dupl_pars <- duplicated(names(forms), fromLast = TRUE)
  if (any(is_dupl_pars)) {
    dupl_pars <- collapse_comma(names(forms)[is_dupl_pars])
    message("Replacing initial definitions of parameters ", dupl_pars)
    forms[is_dupl_pars] <- NULL
  }
  not_form <- ulapply(forms, function(x) !is.formula(x))
  fix <- forms[not_form]
  forms[names(fix)] <- NULL
  out$pforms <- forms
  # validate fixed distributional parameters
  fix_theta <- fix[dpar_class(names(fix)) %in% "theta"]
  if (length(fix_theta)) {
    # normalize mixing proportions
    sum_theta <- sum(unlist(fix_theta))
    fix_theta <- lapply(fix_theta, "/", sum_theta)
    fix[names(fix_theta)] <- fix_theta
  }
  out$pfix <- fix
  for (dp in names(out$pfix)) {
    if (is.character(out$pfix[[dp]])) {
      if (identical(dp, out$pfix[[dp]])) {
        stop2("Equating '", dp, "' with itself is not meaningful.")
      }
      ap_class <- dpar_class(dp)
      if (ap_class == "mu") {
        stop2("Equating parameters of class 'mu' is not allowed.")
      }
      if (!identical(ap_class, dpar_class(out$pfix[[dp]]))) {
        stop2("Can only equate parameters of the same class.")
      }
      if (out$pfix[[dp]] %in% names(out$pfix)) {
        stop2("Cannot use fixed parameters on ",
              "the right-hand side of an equation.")
      }
      if (out$pfix[[dp]] %in% names(out$pforms)) {
        stop2("Cannot use predicted parameters on ",
              "the right-hand side of an equation.")
      }
    }
  }
  if (!is.null(nl)) {
    attr(out$formula, "nl") <- as_one_logical(nl)
  } else if (!is.null(out[["nl"]])) {
    # for backwards compatibility with brms <= 1.8.0
    attr(out$formula, "nl") <- out[["nl"]]
    out[["nl"]] <- NULL
  }
  if (is.null(attr(out$formula, "nl"))) {
    attr(out$formula, "nl") <- FALSE
  }
  if (!is.null(loop)) {
    attr(out$formula, "loop") <- as_one_logical(loop)
  }
  if (is.null(attr(out$formula, "loop"))) {
    attr(out$formula, "loop") <- TRUE
  }
  if (!is.null(center)) {
    attr(out$formula, "center") <- as_one_logical(center)
  }
  if (!is.null(cmc)) {
    attr(out$formula, "cmc") <- as_one_logical(cmc)
  }
  if (!is.null(sparse)) {
    attr(out$formula, "sparse") <- as_one_logical(sparse)
  }
  if (!is.null(decomp)) {
    attr(out$formula, "decomp") <- match.arg(decomp, decomp_opts())
  }
  if (!is.null(unused)) {
    attr(out$formula, "unused") <- as.formula(unused)
  }
  if (!is.null(autocor)) {
    attr(out$formula, "autocor") <- validate_autocor(autocor)
  } else if (!is.null(out$autocor)) {
    # for backwards compatibility with brms <= 2.11.0
    attr(out$formula, "autocor") <- validate_autocor(out$autocor)
    out$autocor <- NULL
  }
  if (!is.null(family)) {
    out$family <- validate_family(family)
  }
  if (!is.null(lhs(formula))) {
    out$resp <- terms_resp(formula)
  }
  # add default values for unspecified elements
  defs <- list(pforms = list(), pfix = list(), family = NULL, resp = NULL)
  defs <- defs[setdiff(names(defs), names(rmNULL(out, FALSE)))]
  out[names(defs)] <- defs
  class(out) <- c("brmsformula", "bform")
  split_bf(out)
}

# alias of brmsformula
#' @export
bf <- function(formula, ..., flist = NULL, family = NULL, autocor = NULL,
               nl = NULL, loop = NULL, center = NULL, cmc = NULL,
               sparse = NULL, decomp = NULL) {
  brmsformula(
    formula, ..., flist = flist, family = family, autocor = autocor,
    nl = nl, loop = loop, center = center, cmc = cmc, sparse = sparse,
    decomp = decomp
  )
}

#' Linear and Non-linear formulas in \pkg{brms}
#'
#' Helper functions to specify linear and non-linear
#' formulas for use with \code{\link[brms:brmsformula]{brmsformula}}.
#'
#' @name brmsformula-helpers
#' @aliases bf-helpers nlf lf set_nl set_rescor
#'
#' @param formula Non-linear formula for a distributional parameter.
#'   The name of the distributional parameter can either be specified
#'   on the left-hand side of \code{formula} or via argument \code{dpar}.
#' @param dpar Optional character string specifying the distributional
#'   parameter to which the formulas passed via \code{...} and
#'   \code{flist} belong.
#' @param resp Optional character string specifying the response
#'   variable to which the formulas passed via \code{...} and
#'   \code{flist} belong. Only relevant in multivariate models.
#' @param autocor A one sided formula containing autocorrelation
#'   terms. All none autocorrelation terms in \code{autocor} will
#'   be silently ignored.
#' @param rescor Logical; Indicates if residual correlation between
#'   the response variables should be modeled. Currently this is only
#'   possible in multivariate \code{gaussian} and \code{student} models.
#'   Only relevant in multivariate models.
#' @param mecor Logical; Indicates if correlations between latent variables
#'   defined by \code{\link{me}} terms should be modeled. Defaults to \code{TRUE}.
#' @inheritParams brmsformula
#'
#' @return For \code{lf} and \code{nlf} a \code{list} that can be
#'   passed to \code{\link[brms:brmsformula]{brmsformula}} or added
#'   to an existing \code{brmsformula} or \code{mvbrmsformula} object.
#'   For \code{set_nl} and \code{set_rescor} a logical value that can be
#'   added to an existing \code{brmsformula} or \code{mvbrmsformula} object.
#'
#' @seealso \code{\link{brmsformula}}, \code{\link{mvbrmsformula}}
#'
#' @examples
#' # add more formulas to the model
#' bf(y ~ 1) +
#'   nlf(sigma ~ a * exp(b * x)) +
#'   lf(a ~ x, b ~ z + (1|g)) +
#'   gaussian()
#'
#' # specify 'nl' later on
#' bf(y ~ a * inv_logit(x * b)) +
#'   lf(a + b ~ z) +
#'   set_nl(TRUE)
#'
#' # specify a multivariate model
#' bf(y1 ~ x + (1|g)) +
#'   bf(y2 ~ z) +
#'   set_rescor(TRUE)
#'
#' # add autocorrelation terms
#' bf(y ~ x) + acformula(~ arma(p = 1, q = 1) + car(W))
NULL

#' @rdname brmsformula-helpers
#' @export
nlf <- function(formula, ..., flist = NULL, dpar = NULL,
                resp = NULL, loop = NULL) {
  formula <- as_formula(formula)
  if (is.null(lhs(formula))) {
    stop2("Argument 'formula' must be two-sided.")
  }
  if (length(c(list(...), flist))) {
    warning2(
      "Arguments '...' and 'flist' in nlf() will be reworked ",
      "at some point. Please avoid using them if possible."
    )
  }
  warn_dpar(dpar)
  if (!is.null(resp)) {
    resp <- as_one_character(resp)
  }
  if (!is.null(loop)) {
    attr(formula, "loop") <- as_one_logical(loop)
  }
  if (is.null(attr(formula, "loop"))) {
    attr(formula, "loop") <- TRUE
  }
  out <- c(
    list(structure(formula, nl = TRUE)),
    lf(..., flist = flist)
  )
  structure(out, resp = resp)
}

#' @rdname brmsformula-helpers
#' @export
lf <- function(..., flist = NULL, dpar = NULL, resp = NULL,
               center = NULL, cmc = NULL, sparse = NULL,
               decomp = NULL) {
  out <- c(list(...), flist)
  warn_dpar(dpar)
  if (!is.null(resp)) {
    resp <- as_one_character(resp)
  }
  cmc <- if (!is.null(cmc)) as_one_logical(cmc)
  center <- if (!is.null(center)) as_one_logical(center)
  decomp <- if (!is.null(decomp)) match.arg(decomp, decomp_opts())
  for (i in seq_along(out)) {
    if (!is.null(cmc)) {
      attr(out[[i]], "cmc") <- cmc
    }
    if (!is.null(center)) {
      attr(out[[i]], "center") <- center
    }
    if (!is.null(sparse)) {
      attr(out[[i]], "sparse") <- sparse
    }
    if (!is.null(decomp)) {
      attr(out[[i]], "decomp") <- decomp
    }
  }
  structure(out, resp = resp)
}

#' @rdname brmsformula-helpers
#' @export
acformula <- function(autocor, resp = NULL) {
  autocor <- terms_ac(as.formula(autocor))
  if (!is.formula(autocor)) {
    stop2("'autocor' must contain at least one autocorrelation term.")
  }
  if (!is.null(resp)) {
    resp <- as_one_character(resp)
  }
  structure(autocor, resp = resp, class = c("acformula", "formula"))
}

#' @rdname brmsformula-helpers
#' @export
set_nl <- function(nl = TRUE, dpar = NULL, resp = NULL) {
  nl <- as_one_logical(nl)
  if (!is.null(dpar)) {
    dpar <- as_one_character(dpar)
  }
  if (!is.null(resp)) {
    resp <- as_one_character(resp)
  }
  structure(nl, dpar = dpar, resp = resp, class = "setnl")
}

#' Set up a multivariate model formula for use in \pkg{brms}
#'
#' Set up a multivariate model formula for use in the \pkg{brms} package
#' allowing to define (potentially non-linear) additive multilevel
#' models for all parameters of the assumed response distributions.
#'
#' @aliases mvbf
#'
#' @param ... Objects of class \code{formula} or \code{brmsformula},
#'   each specifying a univariate model. See \code{\link{brmsformula}}
#'   for details on how to specify univariate models.
#' @param flist Optional list of formulas, which are treated in the
#'   same way as formulas passed via the \code{...} argument.
#' @param rescor Logical; Indicates if residual correlation between
#'   the response variables should be modeled. Currently, this is only
#'   possible in multivariate \code{gaussian} and \code{student} models.
#'   If \code{NULL} (the default), \code{rescor} is internally set to
#'   \code{TRUE} when possible.
#'
#' @return An object of class \code{mvbrmsformula}, which
#'   is essentially a \code{list} containing all model formulas
#'   as well as some additional information for multivariate models.
#'
#' @details See \code{vignette("brms_multivariate")} for a case study.
#'
#' @seealso \code{\link{brmsformula}}, \code{\link{brmsformula-helpers}}
#'
#' @examples
#' bf1 <- bf(y1 ~ x + (1|g))
#' bf2 <- bf(y2 ~ s(z))
#' mvbf(bf1, bf2)
#'
#' @export
mvbrmsformula <- function(..., flist = NULL, rescor = NULL) {
  dots <- c(list(...), flist)
  if (!length(dots)) {
    stop2("No objects passed to 'mvbrmsformula'.")
  }
  forms <- list()
  for (i in seq_along(dots)) {
    if (is.mvbrmsformula(dots[[i]])) {
      forms <- c(forms, dots[[i]]$forms)
      if (is.null(rescor)) {
        rescor <- dots[[i]]$rescor
      }
    } else {
      forms <- c(forms, list(bf(dots[[i]])))
    }
  }
  if (!is.null(rescor)) {
    rescor <- as_one_logical(rescor)
  }
  responses <- ulapply(forms, "[[", "resp")
  if (any(duplicated(responses))) {
    stop2("Cannot use the same response variable twice in the same model.")
  }
  names(forms) <- responses
  structure(
    nlist(forms, responses, rescor),
    class = c("mvbrmsformula", "bform")
  )
}

#' @export
mvbf <- function(..., flist = NULL, rescor = NULL) {
  mvbrmsformula(..., flist = flist, rescor = rescor)
}

# build a mvbrmsformula object based on a brmsformula object
# which uses mvbind on the left-hand side to specify MV models
split_bf <- function(x) {
  stopifnot(is.brmsformula(x))
  resp <- terms_resp(x$formula, check_names = FALSE)
  str_adform <- formula2str(x$formula)
  str_adform <- get_matches("\\|[^~]*(?=~)", str_adform, perl = TRUE)
  if (length(resp) > 1L) {
    # mvbind syntax used to specify MV model
    flist <- named_list(resp)
    for (i in seq_along(resp)) {
      flist[[i]] <- x
      str_lhs <- paste0(resp[[i]], str_adform)
      flist[[i]]$formula[[2]] <- parse(text = str_lhs)[[1]]
      flist[[i]]$resp <- resp[[i]]
    }
    x <- mvbf(flist = flist)
  }
  x
}

#' Bind response variables in multivariate models
#'
#' Can be used to specify a multivariate \pkg{brms} model within a single
#' formula. Outside of \code{\link{brmsformula}}, it just behaves like
#' \code{\link{cbind}}.
#'
#' @param ... Same as in \code{\link{cbind}}
#'
#' @seealso \code{\link{brmsformula}}, \code{\link{mvbrmsformula}}
#'
#' @examples
#' bf(mvbind(y1, y2) ~ x)
#'
#' @export
mvbind <- function(...) {
  cbind(...)
}

#' @rdname brmsformula-helpers
#' @export
set_rescor <- function(rescor = TRUE) {
  structure(as_one_logical(rescor), class = "setrescor")
}

allow_rescor <- function(x) {
  # indicate if estimating 'rescor' is allowed for this model
  if (!(is.mvbrmsformula(x) || is.mvbrmsterms(x))) {
    return(FALSE)
  }
  parts <- if (is.mvbrmsformula(x)) x$forms else x$terms
  families <- lapply(parts, "[[", "family")
  has_rescor <- ulapply(families, has_rescor)
  is_mixture <- ulapply(families, is.mixfamily)
  family_names <- ulapply(families, "[[", "family")
  all(has_rescor) && !any(is_mixture) &&
    length(unique(family_names)) == 1L
}

#' @rdname brmsformula-helpers
#' @export
set_mecor <- function(mecor = TRUE) {
  structure(as_one_logical(mecor), class = "setmecor")
}

#' @export
"+.bform" <- function(e1, e2) {
  if (is.brmsformula(e1)) {
    out <- plus_brmsformula(e1, e2)
  } else if (is.mvbrmsformula(e1)) {
    out <- plus_mvbrmsformula(e1, e2)
  } else {
    stop2("Method '+.bform' not implemented for ", class(e1)[1], " objects.")
  }
  out
}

# internal helper function of '+.bform'
plus_brmsformula <- function(e1, e2) {
  if (is.function(e2)) {
    e2 <- try(e2(), silent = TRUE)
    if (!is.family(e2)) {
      stop2("Don't know how to handle non-family functions.")
    }
  }
  if (is.family(e2)) {
    e1 <- bf(e1, family = e2)
  } else if (is.cor_brms(e2) || inherits(e2, "acformula")) {
    e1 <- bf(e1, autocor = e2)
  } else if (inherits(e2, "setnl")) {
    dpar <- attr(e2, "dpar")
    if (is.null(dpar)) {
      e1 <- bf(e1, nl = e2)
    } else {
      if (is.null(e1$pforms[[dpar]])) {
        stop2("Parameter '", dpar, "' has no formula.")
      }
      attr(e1$pforms[[dpar]], "nl") <- e2
      e1 <- bf(e1)
    }
  } else if (inherits(e2, "setmecor")) {
    e1$mecor <- e2[1]
  } else if (is.brmsformula(e2)) {
    e1 <- mvbf(e1, e2)
  } else if (inherits(e2, "setrescor")) {
    stop2("Setting 'rescor' is only possible in multivariate models.")
  } else if (is.ac_term(e2)) {
    stop2("Autocorrelation terms can only be specified on the right-hand ",
          "side of a formula, not added to a 'brmsformula' object.")
  } else if (!is.null(e2)) {
    e1 <- bf(e1, e2)
  }
  e1
}

# internal helper function of '+.bform'
plus_mvbrmsformula <- function(e1, e2) {
  if (is.function(e2)) {
    e2 <- try(e2(), silent = TRUE)
    if (!is.family(e2)) {
      stop2("Don't know how to handle non-family functions.")
    }
  }
  if (is.family(e2) || is.cor_brms(e2)) {
    e1$forms <- lapply(e1$forms, "+", e2)
  } else if (inherits(e2, "setrescor")) {
    e1$rescor <- e2[1]
  } else if (inherits(e2, "setmecor")) {
    e1$mecor <- e2[1]
  } else if (is.brmsformula(e2)) {
    e1 <- mvbf(e1, e2)
  } else if (is.mvbrmsformula(e2)) {
    # TODO: enable this option
    stop2("Cannot add two 'mvbrmsformula' objects together. Instead, ",
          "please add the individual 'brmsformula' objects directly.")
  } else if (is.ac_term(e2)) {
    stop2("Autocorrelation terms can only be specified on the right-hand ",
          "side of a formula, not added to a 'mvbrmsformula' object.")
  } else if (!is.null(e2)) {
    resp <- attr(e2, "resp", TRUE)
    if (is.null(resp)) {
      stop2(
        "Don't know how to add a ", class(e2)[1], " object ",
        "without the response variable name. ",
        "See help('brmsformula-helpers') for more details."
      )
    }
    if (!isTRUE(resp %in% e1$responses)) {
      stop2("'resp' should be one of ", collapse_comma(e1$responses), ".")
    }
    e1$forms[[resp]] <- e1$forms[[resp]] + e2
  }
  e1
}

# extract the 'nl' attribute from a brmsformula object
# @param x object to extract 'nl' from
# @param dpar optional name of a distributional parameter
#   for which 'nl' should be extracted
# @param resp: optional name of a response variable
#   for which 'nl' should be extracted
# @param aol: (as one logical) apply isTRUE to the result?
get_nl <- function(x, dpar = NULL, resp = NULL, aol = TRUE) {
  if (is.mvbrmsformula(x)) {
    resp <- as_one_character(resp)
    x <- x$forms[[resp]]
  }
  if (is.brmsformula(x)) {
    if (is.null(dpar)) {
      x <- x$formula
    } else {
      dpar <- as_one_character(dpar)
      x <- x$pforms[[dpar]]
    }
  }
  nl <- attr(x, "nl", TRUE)
  if (aol) {
    nl <- isTRUE(nl)
  }
  nl
}

# available options for the 'decomp' argument
decomp_opts <- function() {
  c("none", "QR")
}

# validate a formula of an additional parameter
# @param formula an formula object
# @param par optional name of the parameter; if not specified
#   the parameter name will be inferred from the formula
# @param rsv_pars optional character vector of reserved parameter names
# @return a named list of length one containing the formula
validate_par_formula <- function(formula, par = NULL, rsv_pars = NULL) {
  stopifnot(length(par) <= 1L)
  try_formula <- try(as_formula(formula), silent = TRUE)
  if (is(try_formula, "try-error")) {
    if (length(formula) != 1L) {
      stop2("Expecting a single value when fixing parameter '", par, "'.")
    }
    scalar <- SW(as.numeric(formula))
    if (!is.na(scalar)) {
      formula <- scalar
    } else {
      formula <- as.character(formula)
    }
    out <- named_list(par, formula)
  } else {
    formula <- try_formula
    if (!is.null(lhs(formula))) {
      resp_pars <- all.vars(formula[[2]])
      out <- named_list(resp_pars, list(formula))
      for (i in seq_along(out)) {
        out[[i]][[2]] <- eval2(paste("quote(", resp_pars[i], ")"))
      }
    } else {
      if (!isTRUE(nzchar(par))) {
        stop2("Additional formulas must be named.")
      }
      formula <- formula(paste(par, formula2str(formula)))
      out <- named_list(par, list(formula))
    }
  }
  pars <- names(out)
  if (any(grepl("\\.|_", pars))) {
    stop2("Parameter names should not contain dots or underscores.")
  }
  inv_pars <- intersect(pars, rsv_pars)
  if (length(inv_pars)) {
    stop2("The following parameter names are reserved",
          "for this model:\n", collapse_comma(inv_pars))
  }
  out
}

# validate formulas dedicated to response variables
# @param x coerced to a formula object
# @param empty_ok is an empty left-hand-side ok?
# @return a formula of the form <response> ~ 1
validate_resp_formula <- function(x, empty_ok = TRUE) {
  out <- lhs(as_formula(x))
  if (is.null(out)) {
    if (empty_ok) {
      out <- ~ 1
    } else {
      str_x <- formula2str(x, space = "trim")
      stop2("Response variable is missing in formula ", str_x)
    }
  }
  out <- gsub("\\|+[^~]*~", "~", formula2str(out))
  out <- try(formula(out), silent = TRUE)
  if (is(out, "try-error")) {
    str_x <- formula2str(x, space = "trim")
    stop2("Incorrect use of '|' on the left-hand side of ", str_x)
  }
  environment(out) <- environment(x)
  out
}

# incorporate additional arguments into the model formula
validate_formula <- function(formula, ...) {
  UseMethod("validate_formula")
}

#' @export
validate_formula.default <- function(formula, ...) {
  validate_formula(bf(formula), ...)
}

# incorporate additional arguments into the model formula
# @param formula object of class 'formula' of 'brmsformula'
# @param data optional data.frame to validate data related arguments
# @param family optional 'family' object
# @param autocor (deprecated) optional 'cor_brms' object
# @param threshold (deprecated) threshold type for ordinal models
# @param cov_ranef (deprecated) named list of group covariance matrices
# @return a brmsformula object compatible with the current version of brms
#' @export
validate_formula.brmsformula <- function(
  formula, family = gaussian(), autocor = NULL,
  data = NULL, threshold = NULL, sparse = NULL,
  cov_ranef = NULL, ...
) {
  out <- bf(formula)
  if (is.null(out$family) && !is.null(family)) {
    out$family <- validate_family(family)
  }
  # allow the '.' symbol in the formulas
  out$formula <- expand_dot_formula(out$formula, data)
  for (i in seq_along(out$pforms)) {
    out$pforms[[i]] <- expand_dot_formula(out$pforms[[i]], data)
  }
  # allow 'me' terms to be correlated
  out$mecor <- default_mecor(out$mecor)
  if (has_cat(out) && !is.null(data)) {
    # for easy access of response categories
    # allow to update 'cats' with new data
    out$family$cats <- extract_cat_names(out, data)
  }
  if (is_ordinal(out$family)) {
    # thresholds and category names are data dependent
    try_terms <- try(stats::terms(out$formula), silent = TRUE)
    intercept <- attr(try_terms, "intercept", TRUE)
    if (!is(try_terms, "try-error") && isTRUE(intercept == 0)) {
      stop2("Cannot remove the intercept in an ordinal model.")
    }
    if (!is.null(data)) {
      # for easy access of thresholds and response categories (#838)
      # allow to update 'cats' and 'thres' with new data
      out$family$thres <- extract_thres_names(out, data)
      out$family$cats <- extract_cat_names(out, data)
    }
    if (is.mixfamily(out$family)) {
      # every mixture family needs to know about response categories
      for (i in seq_along(out$family$mix)) {
        out$family$mix[[i]]$thres <- out$family$thres
      }
    }
  }
  conv_cats_dpars <- conv_cats_dpars(out$family)
  if (conv_cats_dpars && !is.null(data)) {
    # allow to update 'dpars' with new data
    # define distributional parameters based on response categories
    if (length(out$family$cats) < 2L) {
      stop2("At least 2 response categories are required.")
    }
    if (is.null(out$family$refcat)) {
      # the first level serves as the reference category
      out$family$refcat <- out$family$cats[1]
    }
    if (isNA(out$family$refcat)) {
      # implies predicting all categories
      predcats <- out$family$cats
    } else {
      if (!out$family$refcat %in% out$family$cats) {
        stop2("The reference response category must be one of ",
              collapse_comma(out$family$cats), ".")
      }
      predcats <- setdiff(out$family$cats, out$family$refcat)
    }
    multi_dpars <- valid_dpars(out$family, type = "multi")
    # 'rev' so that mu comes last but appears first in the end
    for (dp in rev(multi_dpars)) {
      dp_dpars <- make_stan_names(paste0(dp, predcats))
      if (any(duplicated(dp_dpars))) {
        stop2("Invalid response category names. Please avoid ",
              "using any special characters in the names.")
      }
      old_dp_dpars <- str_subset(out$family$dpars, paste0("^", dp))
      out$family$dpars <- setdiff(out$family$dpars, old_dp_dpars)
      out$family$dpars <- union(dp_dpars, out$family$dpars)
    }
  }

  # incorporate deprecated arguments
  require_threshold <- is_ordinal(out$family) && is.null(out$family$threshold)
  if (require_threshold && !is.null(threshold)) {
    # slot 'threshold' is deprecated as of brms 1.7.0
    out$family <- validate_family(out$family, threshold = threshold)
  }
  if (!is.null(sparse)) {
    # a global 'sparse' argument is deprecated as of brms 2.8.3
    warning2(
      "Argument 'sparse' should be specified within the ",
      "'formula' argument. See ?brmsformula for help."
    )
    sparse <- as_one_logical(sparse)
    if (is.null(attr(out$formula, "sparse"))) {
      attr(out$formula, "sparse") <- sparse
    }
    for (i in seq_along(out$pforms)) {
      if (is.null(attr(out$pforms[[i]], "sparse"))) {
        attr(out$pforms[[i]], "sparse") <- sparse
      }
    }
  }
  if (is.null(attr(out$formula, "autocor")) && !is.null(autocor)) {
    # 'autocor' interface has been changed in brms 2.11.1
    warning2(
      "Argument 'autocor' should be specified within the ",
      "'formula' argument. See ?brmsformula for help."
    )
    # store 'autocor' as an attribute to carry it around more easily
    attr(out$formula, "autocor") <- validate_autocor(autocor)
  }
  if (!is.null(cov_ranef)) {
    # 'cov_ranef' is deprecated as of brms 2.12.5
    out$cov_ranef <- validate_cov_ranef(cov_ranef)
  }
  bf(out)
}

# incorporate additional arguments into MV model formulas
# allow passing lists of families or autocors
#' @export
validate_formula.mvbrmsformula <- function(
  formula, family = NULL, autocor = NULL, cov_ranef = NULL, ...
) {
  nresp <- length(formula$forms)
  if (!is(family, "list")) {
    family <- replicate(nresp, family, simplify = FALSE)
  } else if (length(family) != nresp) {
    stop2("If 'family' is a list, it has to be of the same ",
          "length as the number of response variables.")
  }
  if (!is(autocor, "list")) {
    autocor <- replicate(nresp, autocor, simplify = FALSE)
  } else if (length(autocor) != nresp) {
    stop2("If 'autocor' is a list, it has to be of the same ",
          "length as the number of response variables.")
  }
  for (i in seq_len(nresp)) {
    formula$forms[[i]] <- validate_formula(
      formula$forms[[i]], family = family[[i]],
      autocor = autocor[[i]], ...
    )
  }
  if (length(formula$forms) < 2L) {
    stop2("Multivariate models require at least two responses.")
  }
  allow_rescor <- allow_rescor(formula)
  if (is.null(formula$rescor)) {
    # with 'mi' terms we usually don't want rescor to be estimated
    miforms <- ulapply(formula$forms, function(f)
      terms_ad(f$formula, f$family, FALSE)[["mi"]]
    )
    formula$rescor <- allow_rescor && !length(miforms)
    message("Setting 'rescor' to ", formula$rescor, " by default for this model")
    if (formula$rescor) {
      warning2(
        "In the future, 'rescor' will be set to FALSE by default for ",
        "all models. It is thus recommended to explicitely set ",
        "'rescor' via 'set_rescor' instead of using the default."
      )
    }
  }
  formula$rescor <- as_one_logical(formula$rescor)
  if (formula$rescor) {
    if (!allow_rescor) {
      stop2("Currently, estimating 'rescor' is only possible ",
            "in multivariate gaussian or student models.")
    }
  }
  # handle default of correlations between 'me' terms
  formula$mecor <- default_mecor(formula$mecor)
  for (i in seq_along(formula$forms)) {
    formula$forms[[i]]$mecor <- formula$mecor
  }
  # incorporate deprecated arguments
  if (!is.null(cov_ranef)) {
    # 'cov_ranef' is deprecated as of brms 2.12.5
    formula$cov_ranef <- validate_cov_ranef(cov_ranef)
  }
  formula
}

# update a brmsformula and / or its attributes
# @param brmsformula object
# @param formula.: formula to update 'object'
# @param mode supports the following options:
#   "update": apply update.formula
#   "replace": replace old formula
#   "keep": keep old formula
#   attributes are always updated
# @param ... currently unused
# @return a brmsformula object
#' @export
update.brmsformula <- function(object, formula.,
                               mode = c("update", "replace", "keep"),
                               ...) {
  mode <- match.arg(mode)
  object <- bf(object)
  up_nl <- get_nl(formula., aol = FALSE)
  if (is.null(up_nl)) {
    up_nl <- get_nl(object)
  }
  # already use up_nl here to avoid ordinary parsing of NL formulas
  formula. <- bf(formula., nl = up_nl)
  up_family <- formula.[["family"]]
  if (is.null(up_family)) {
    up_family <- object[["family"]]
  }
  up_autocor <- attr(formula.$formula, "autocor")
  if (is.null(up_autocor)) {
    up_autocor <- attr(object$formula, "autocor")
  }
  old_form <- object$formula
  up_form <- formula.$formula
  if (mode == "update") {
    new_form <- update(old_form, up_form, ...)
  } else if (mode == "replace") {
    new_form <- up_form
  } else if (mode == "keep") {
    new_form <- old_form
  }
  flist <- c(object$pforms, object$pfix, formula.$pforms, formula.$pfix)
  bf(new_form, flist = flist, family = up_family,
     autocor = up_autocor, nl = up_nl)
}

#' @export
update.mvbrmsformula <- function(object, formula., ...) {
  # temporary until proper updating is implemented
  if (!missing(formula.)) {
    stop2("Updating formulas of multivariate models is not yet possible.")
  }
  object
}

#' Update Formula Addition Terms
#'
#' Update additions terms used in formulas of \pkg{brms}. See
#' \code{\link{addition-terms}} for details.
#'
#' @param formula Two-sided formula to be updated.
#' @param adform One-sided formula containing addition terms to update
#'   \code{formula} with.
#' @param action Indicates what should happen to the existing addition terms in
#'   \code{formula}. If \code{"update"} (the default), old addition terms that
#'   have no corresponding term in \code{adform} will be kept. If
#'   \code{"replace"}, all old addition terms will be removed.
#'
#' @return An object of class \code{formula}.
#'
#' @examples
#' form <- y | trials(size) ~ x
#' update_adterms(form, ~ trials(10))
#' update_adterms(form, ~ weights(w))
#' update_adterms(form, ~ weights(w), action = "replace")
#' update_adterms(y ~ x, ~ trials(10))
#'
#' @export
update_adterms <- function(formula, adform, action = c("update", "replace")) {
  formula <- as_formula(formula)
  adform <- as_formula(adform)
  action <- match.arg(action)
  if (is.null(lhs(formula))) {
    stop2("Can't update a ond-sided formula.")
  }
  str_formula <- formula2str(formula)
  old_ad <- get_matches("(?<=\\|)[^~]*(?=~)", str_formula, perl = TRUE)
  new_ad_terms <- attr(terms(adform), "term.labels")
  if (action == "update" && length(old_ad)) {
    # extract adterms from the original formula
    old_ad <- formula(paste("~", old_ad))
    old_ad_terms <- attr(terms(old_ad), "term.labels")
    old_adnames <- get_matches("^[^\\(]+", old_ad_terms)
    old_adnames <- sub("^resp_", "", old_adnames)
    new_adnames <- get_matches("^[^\\(]+", new_ad_terms)
    new_adnames <- sub("^resp_", "", new_adnames)
    # keep unmatched adterms of the original formula
    keep <- !old_adnames %in% new_adnames
    new_ad_terms <- c(old_ad_terms[keep], new_ad_terms)
  }
  if (length(new_ad_terms)) {
    new_ad_terms <- paste(new_ad_terms, collapse = "+")
    new_ad_terms <- paste("|", new_ad_terms)
  }
  resp <- gsub("\\|.+", "", deparse0(formula[[2]]))
  out <- formula(paste(resp, new_ad_terms, "~1"))
  out[[3]] <- formula[[3]]
  attributes(out) <- attributes(formula)
  out
}

#' @export
print.brmsformula <- function(x, wsp = 0, digits = 2, ...) {
  cat(formula2str(x$formula, space = "trim"), "\n")
  str_wsp <- collapse(rep(" ", wsp))
  autocor <- attr(x$formula, "autocor")
  if (!is.null(autocor)) {
    autocor <- formula2str(autocor, rm = 1, space = "trim")
    cat(paste0(str_wsp, "autocor ~ ", autocor, "\n"))
  }
  pforms <- x$pforms
  if (length(pforms)) {
    pforms <- ulapply(pforms, formula2str, space = "trim")
    cat(collapse(str_wsp, pforms, "\n"))
  }
  pfix <- x$pfix
  if (length(pfix)) {
    pfix <- lapply(pfix, function(x)
      ifelse(is.numeric(x), round(x, digits), x)
    )
    pfix <- paste0(names(pfix), " = ", unlist(pfix))
    cat(collapse(str_wsp, pfix, "\n"))
  }
  invisible(x)
}

#' @export
print.mvbrmsformula <- function(x, wsp = 0, ...) {
  for (i in seq_along(x$forms)) {
    if (i > 1) cat(collapse(rep(" ", wsp)))
    print(x$forms[[i]], wsp = wsp, ...)
  }
  invisible(x)
}

#' Checks if argument is a \code{brmsformula} object
#'
#' @param x An \R object
#'
#' @export
is.brmsformula <- function(x) {
  inherits(x, "brmsformula")
}

#' Checks if argument is a \code{mvbrmsformula} object
#'
#' @param x An \R object
#'
#' @export
is.mvbrmsformula <- function(x) {
  inherits(x, "mvbrmsformula")
}

is_nonlinear <- function(x) {
  stopifnot(is.brmsfit(x))
  get_nl(bf(x$formula))
}

warn_dpar <- function(dpar) {
  # argument 'dpar' in formula helper functions is deprecated as of 2.3.7
  if (!is.null(dpar)) {
    warning2("Argument 'dpar' is no longer necessary and ignored.")
  }
  NULL
}

# return the right-hand side of a formula
rhs <- function(x) {
  attri <- attributes(x)
  x <- as.formula(x)
  x <- if (length(x) == 3) x[-2] else x
  do_call(structure, c(list(x), attri))
}

# return the left-hand side of a formula
lhs <- function(x) {
  x <- as.formula(x)
  if (length(x) == 3L) update(x, . ~ 1) else NULL
}

# convert a string to a formula
# @param x vector of strings to be converted
# @param ... passed to formula()
str2formula <- function(x, env = parent.frame(), collapse = "+") {
  has_chars <- nzchar(x)
  if (length(x) && any(has_chars)) {
    out <- paste(x[has_chars], collapse = collapse)
  } else {
    out <- "1"
  }
  as.formula(paste("~", out), env = env)
}

# convert a formula to a character string
# @param formula a model formula
# @param rm a vector of to elements indicating how many characters
#   should be removed at the beginning and end of the string respectively
# @param space how should whitespaces be treated?
#    option 'rm' is dangerous as it may combine different operators (#1142)
# @return a single character string or NULL
formula2str <- function(formula, rm = c(0, 0), space = c("trim", "rm")) {
  if (is.null(formula)) {
    return(NULL)
  }
  formula <- as.formula(formula)
  space <- match.arg(space)
  if (anyNA(rm[2])) rm[2] <- 0
  x <- Reduce(paste, deparse(formula))
  x <- gsub("[\t\r\n]+", " ", x, perl = TRUE)
  if (space == "trim") {
    x <- trim_wsp(x)
  } else {
    x <- rm_wsp(x)
  }
  substr(x, 1 + rm[1], nchar(x) - rm[2])
}

# right-hand side of a formula as a character string
str_rhs <- function(x) {
  formula2str(rhs(x), rm = c(1, 0))
}

# left-hand side of a formula as a character string
str_lhs <- function(x) {
  formula2str(lhs(x), rm = c(0, 2))
}

is.formula <- function(x) {
  inherits(x, "formula")
}

# wrapper around as.formula with additional checks
as_formula <- function(x) {
  x <- as.formula(x)
  # fixes issue #749
  rhs <- rhs(x)[[2]]
  if (isTRUE(is.call(rhs) && rhs[[1]] == "~")) {
    stop2("Nested formulas are not allowed. Did you use '~~' somewhere?")
  }
  x
}

# expand the '.' variable in formula using stats::terms
expand_dot_formula <- function(formula, data = NULL) {
  if (isTRUE("." %in% all.vars(formula))) {
    att <- attributes(formula)
    try_terms <- try(
      stats::terms(formula, data = data),
      silent = TRUE
    )
    if (!is(try_terms, "try-error")) {
      formula <- formula(try_terms)
    }
    attributes(formula) <- att
  }
  formula
}
