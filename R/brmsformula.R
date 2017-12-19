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
#' @param ... Additional \code{formula} objects to specify 
#'   predictors of non-linear and distributional parameters. 
#'   Formulas can either be named directly or contain
#'   names on their left-hand side. 
#'   The following are distributional parameters of specific families
#'   (all other parameters are treated as non-linear parameters):
#'   \code{sigma} (residual standard deviation or scale of
#'   the \code{gaussian}, \code{student}, \code{skew_normal}, 
#'   \code{lognormal} \code{exgaussian}, and \code{asym_laplace} families);
#'   \code{shape} (shape parameter of the \code{Gamma},
#'   \code{weibull}, \code{negbinomial}, and related
#'   zero-inflated / hurdle families); \code{nu} (degrees of freedom 
#'   parameter of the \code{student} and \code{frechet} families);
#'   \code{phi} (precision parameter of the \code{beta} 
#'   and \code{zero_inflated_beta} families);
#'   \code{kappa} (precision parameter of the \code{von_mises} family);
#'   \code{beta} (mean parameter of the exponential component
#'   of the \code{exgaussian} family);
#'   \code{quantile} (quantile parameter of the \code{asym_laplace} family);
#'   \code{zi} (zero-inflation probability); 
#'   \code{hu} (hurdle probability);
#'   \code{zoi} (zero-one-inflation probability);
#'   \code{coi} (conditional one-inflation probability);
#'   \code{disc} (discrimination) for ordinal models;
#'   \code{bs}, \code{ndt}, and \code{bias} (boundary separation,
#'   non-decision time, and initial bias of the \code{wiener}
#'   diffusion model).
#'   By default, distributional parameters are modeled 
#'   on the log scale if they can be positive only or on the 
#'   logit scale if the can only be within the unit interval.
#'   See 'Details' for more explanation.
#' @param flist Optional list of formulas, which are treated in the 
#'   same way as formulas passed via the \code{...} argument.
#' @param nl Logical; Indicates whether \code{formula} should be
#'   treated as specifying a non-linear model. By default, \code{formula} 
#'   is treated as an ordinary linear model formula.
#' @param family Same argument as in \code{\link{brm}}.
#'   If \code{family} is specified in \code{brmsformula}, it will 
#'   overwrite the value specified in \code{\link{brm}}.
#' @param autocor Same argument as in \code{\link{brm}}.
#'   If \code{autocor} is specified in \code{brmsformula}, it will 
#'   overwrite the value specified in \code{\link{brm}}.
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
#'   The \code{pterms} part contains effects that are assumed to be the 
#'   same across observations. We call them 'population-level' effects
#'   or (adopting frequentist vocabulary) 'fixed' effects. The optional
#'   \code{gterms} part may contain effects that are assumed to vary
#'   across grouping variables specified in \code{group}. We
#'   call them 'group-level' effects or (adopting frequentist 
#'   vocabulary) 'random' effects, although the latter name is misleading
#'   in a Bayesian context. For more details type 
#'   \code{vignette("brms_overview")} and \code{vignette("brms_multilevel")}. 
#'   
#'   \bold{Group-level terms}
#'   
#'   Multiple grouping factors each with multiple group-level effects 
#'   are possible. (Of course we can also run models without any
#'   group-level effects.) 
#'   Instead of \code{|} you may use \code{||} in grouping terms
#'   to prevent correlations from being modeled. 
#'   Alternatively, it is possible to model different group-level terms of 
#'   the same grouping factor as correlated (even across different formulas,
#'   e.g., in non-linear models) by using \code{|<ID>|} instead of \code{|}.
#'   All group-level terms sharing the same ID will be modeled as correlated.
#'   If, for instance, one specifies the terms \code{(1+x|2|g)} and 
#'   \code{(1+z|2|g)} somewhere in the formulas passed to \code{brmsformula},
#'   correlations between the corresponding group-level effects 
#'   will be estimated. 
#'   
#'   You can specify multi-membership terms
#'   using the \code{\link[brms:mm]{mm}} function. For instance, 
#'   a multi-membership term with two members could be
#'   \code{(1|mm(g1, g2))}, where \code{g1} and \code{g2} specify
#'   the first and second member, respectively.
#'   
#'   \bold{Special predictor terms}
#'   
#'   Smoothing terms can modeled using the \code{\link{s}}
#'   and \code{\link{t2}} functions in the \code{pterms} part 
#'   of the model formula. This allows to fit generalized additive mixed
#'   models (GAMMs) with \pkg{brms}. The implementation is similar to that 
#'   used in the \pkg{gamm4} package. For more details on this model class 
#'   see \code{\link[mgcv:gam]{gam}} and \code{\link[mgcv:gamm]{gamm}}.
#'   
#'   Gaussian process terms can be fitted using the \code{\link[brms:gp]{gp}}
#'   function in the \code{pterms} part of the model formula. Similar to
#'   smooth terms, Gaussian processes can be used to model complex non-linear
#'   relationships, for instance temporal or spatial autocorrelation. 
#'   However, they are computationally demanding and are thus not recommended 
#'   for very large datasets.
#'   
#'   The \code{pterms} and \code{gterms} parts may contain three non-standard
#'   effect types namely monotonic, measurement error, and category specific effects,
#'   which can be specified using terms of the form \code{mo(predictor)},
#'   \code{me(predictor, sd_predictor)}, and \code{cs(<predictors>)}, 
#'   respectively. Category specific effects can only be estimated in
#'   ordinal models and are explained in more detail in the package's 
#'   main vignette (type \code{vignette("brms_overview")}). 
#'   The other two effect types are explained in the following.
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
#'   \code{y ~ me(x, sdx) * me(x2, sdx2) * z}. In future version of \pkg{brms},
#'   a vignette will be added to explain more details about these
#'   so called 'error-in-variables' models and provide real world examples.
#'   
#'   \bold{Additional response information}
#'   
#'   Another special of the \pkg{brms} formula syntax is the optional 
#'   \code{aterms} part, which may contain multiple terms of the form 
#'   \code{fun(<variable>)} separated by \code{+} each providing special 
#'   information on the response variable. \code{fun} can be replaced with 
#'   either \code{se}, \code{weights}, \code{cens}, \code{trunc}, 
#'   \code{trials}, \code{cat}, or \code{dec}. Their meanings are explained below.
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
#'   With the exception of categorical, ordinal, and mixture families, 
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
#'   With the exception of categorical, ordinal, and mixture families, 
#'   the response distribution can be truncated using the \code{trunc} 
#'   function in the addition part. If the response variable is truncated 
#'   between, say, 0 and 100, we can specify this via
#'   \code{yi | trunc(lb = 0, ub = 100) ~ predictors}. 
#'   Instead of numbers, variables in the data set can also be passed allowing 
#'   for varying truncation points across observations. Defining only one of 
#'   the two arguments in \code{trunc} leads to one-sided truncation.
#'   
#'   For families \code{binomial} and \code{zero_inflated_binomial}, 
#'   addition should contain a variable indicating the number of trials 
#'   underlying each observation. In \code{lme4} syntax, we may write for instance 
#'   \code{cbind(success, n - success)}, which is equivalent
#'   to \code{success | trials(n)} in \pkg{brms} syntax. If the number of trials
#'   is constant across all observations, say \code{10}, 
#'   we may also write \code{success | trials(10)}. 
#'   
#'   For all ordinal families, \code{aterms} may contain a term 
#'   \code{cat(number)} to specify the number categories (e.g, \code{cat(7)}). 
#'   If not given, the number of categories is calculated from the data.
#'   
#'   In Wiener diffusion models (family \code{wiener}) the addition term
#'   \code{dec} is mandatory to specify the (vector of) binary decisions 
#'   corresponding to the reaction times. Non-zero values will be treated
#'   as a response on the upper boundary of the diffusion process and zeros
#'   will be treated as a response on the lower boundary. Alternatively,
#'   the variable passed to \code{dec} might also be a character vector 
#'   consisting of \code{'lower'} and \code{'upper'}.
#' 
#'   Multiple addition terms may be specified at the same time using 
#'   the \code{+} operator, for instance \cr
#'   \code{formula = yi | se(sei) + cens(censored) ~ 1} 
#'   for a censored meta-analytic model. 
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
#'   The population-level intercept (if incorporated) is estimated separately 
#'   and not as part of population-level parameter vector \code{b}. 
#'   As a result, priors on the intercept also have to be specified separately.
#'   Furthermore, to increase sampling efficiency, the population-level 
#'   design matrix \code{X} is centered around its column means 
#'   \code{X_means} if the intercept is incorporated. 
#'   This leads to a temporary bias in the intercept equal to 
#'   \code{<X_means, b>}, where \code{<,>} is the scalar product. 
#'   The bias is corrected after fitting the model, but be aware 
#'   that you are effectively defining a prior on the intercept 
#'   of the centered design matrix not on the real intercept.
#'   For more details on setting priors on population-level intercepts,
#'   see \code{\link{set_prior}}.
#'   
#'   This behavior can be avoided by using the reserved 
#'   (and internally generated) variable \code{intercept}. 
#'   Instead of \code{y ~ x}, you may write
#'   \code{y ~ 0 + intercept + x}. This way, priors can be
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
#'   To tell \code{brms} that this is a non-linear model, 
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
#'   Non-linear models may not be uniquely identified and / or show bad convergence.
#'   For this reason it is mandatory to specify priors on the non-linear parameters.
#'   For instructions on how to do that, see \code{\link{set_prior}}.
#'   For some examples of non-linear models, see \code{vignette("brms_nonlinear")}.
#'   
#'   \bold{Formula syntax for predicting distributional parameters}
#'   
#'   It is also possible to predict parameters of the response
#'   distribution such as the residual standard deviation \code{sigma} 
#'   in gaussian models or the hurdle probability \code{hu} in hurdle models. 
#'   The syntax closely resembles that of a non-linear 
#'   parameter, for instance \code{sigma ~ x + s(z) + (1+x|g)}. 
#'   For some examples of distributional models, see \code{vignette("brms_distreg")}.
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
#'   distributional parameters of mixture distributions have the same 
#'   name as those of the corresponding ordinary distributions, but with 
#'   a number at the end to indicate the mixture component. For instance, if
#'   you use family \code{mixture(gaussian, gaussian)}, the distributional
#'   parameters are \code{sigma1} and \code{sigma2}.
#'   distributional parameters of the same class can be fixed to the same value. 
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
#'   as the reference category to identify the model. The \code{softmax} 
#'   function is applied on the linear predictor terms to form a 
#'   probability vector.
#'   
#'   For more information on mixture models, see
#'   the documentation of \code{\link{mixture}}.
#'   
#'   \bold{Formula syntax for multivariate models}
#'   
#'   Multivariate models may be specified using \code{cbind} notation
#'   or with help of the \code{\link{mvbf}} function.
#'   Suppose that \code{y1} and \code{y2} are response variables 
#'   and \code{x} is a predictor. Then \code{cbind(y1, y2) ~ x} 
#'   specifies a multivariate model,
#'   The effects of all terms specified at the RHS of the formula 
#'   are assumed to vary across response variables. 
#'   For instance, two parameters will be estimated for \code{x}, 
#'   one for the effect on \code{y1} and another for the effect on \code{y2}.
#'   This is also true for group-level effects. When writing, for instance,
#'   \code{cbind(y1, y2) ~ x + (1+x|g)}, group-level effects will be
#'   estimated separately for each response. To model these effects
#'   as correlated across responses, use the ID syntax (see above).
#'   For the present example, this would look as follows:
#'   \code{cbind(y1, y2) ~ x + (1+x|2|g)}. Of course, you could also use
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
#' bf(y ~ a1 - a2^x, a1 ~ 1 + (1|2|g), a2 ~ x + (x|2|g), nl = TRUE)
#' 
#' # define a multivariate model
#' bf(cbind(y1, y2) ~ x * z + (1|g))
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
#' @export
brmsformula <- function(formula, ..., flist = NULL, family = NULL,
                        autocor = NULL, nl = NULL) {
  if (is.brmsformula(formula)) {
    out <- formula
  } else {
    out <- list(formula = as.formula(formula))
    class(out) <- "brmsformula"
  }
  # parse and validate dots arguments
  dots <- c(out$pforms, out$pfix, list(...), flist)
  dots <- lapply(dots, function(x) if (is.list(x)) x else list(x))
  dots <- unlist(dots, recursive = FALSE)
  forms <- list()
  for (i in seq_along(dots)) {
    forms <- c(forms, prepare_auxformula(dots[[i]], par = names(dots)[i]))
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
  if (!is.null(family)) {
    out$family <- check_family(family)
  }
  if (!is.null(autocor)) {
    out$autocor <- check_autocor(autocor)
  }
  respform <- lhs(formula)
  if (!is.null(respform)) {
    respform <- formula(gsub("\\|+[^~]*~", "~", formula2str(respform)))
    out$resp <- parse_resp(respform)
  }
  if (!is.null(out$family)) {
    # check for the presence of non-linear parameters
    dpars <- names(out$pforms)
    dpars <- dpars[is_dpar_name(dpars, out$family)]
    for (dp in names(out$pforms)) {
      if (!dp %in% dpars) {
        # indicate the correspondence to distributional parameter 
        if (is.null(attr(out$pforms[[dp]], "dpar"))) {
          attr(out$pforms[[dp]], "dpar") <- "mu"
        }
        dpar <- attr(out$pforms[[dp]], "dpar")
        if (!is.null(out$pforms[[dpar]])) {
          nl_allowed <- get_nl(out, dpar = dpar)
        } else {
          if (dpar_class(dpar) == "mu") {
            nl_allowed <- get_nl(out)
          } else {
            nl_allowed <- FALSE
          }
        }
        if (!nl_allowed) {
          stop2(
            "The parameter '", dp, "' is not a valid ", 
            "distributional or non-linear parameter. ",
            "Did you forget to set 'nl = TRUE'?"
          )
        }
      }
    }
  }
  # add default values for unspecified elements
  defs <- list(
    pforms = list(), pfix = list(), family = NULL, 
    autocor = NULL, resp = NULL
  )
  defs <- defs[setdiff(names(defs), names(rmNULL(out, FALSE)))]
  out[names(defs)] <- defs
  class(out) <- c("brmsformula", "bform")
  split_bf(out)
}

#' @export
bf <- function(formula, ..., flist = NULL, family = NULL, 
               autocor = NULL, nl = NULL) {
  # alias of brmsformula
  brmsformula(formula, ..., flist = flist, family = family,
              autocor = autocor, nl = nl)
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
#' @param rescor Logical; Indicates if residual correlation between
#'   the response variables should be modeled. Currently this is only
#'   possible in multivariate \code{gaussian} and \code{student} models.
#'   Only relevant in multivariate models.
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
#'   nlf(sigma ~ a * exp(b * x), a ~ x) + 
#'   lf(b ~ z + (1|g), dpar = "sigma") +
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
NULL

#' @rdname brmsformula-helpers
#' @export
nlf <- function(formula, ..., flist = NULL, dpar = NULL, resp = NULL) {
  formula <- as.formula(formula)
  resp_pars <- all.vars(formula[[2]])
  if (length(resp_pars) == 0L) {
    if (is.null(dpar)) {
      stop2("No parameter name passed via the LHS of ", 
            "'formula' or argument 'dpar'.")
    }
  } else if (length(resp_pars) == 1L) {
    dpar <- resp_pars
  } else {
    stop2("LHS of non-linear formula should contain only one variable.")
  }
  out <- c(
    setNames(list(structure(formula, nl = TRUE)), dpar),
    lf(..., flist = flist, dpar = dpar)
  )
  structure(out, dpar = dpar, resp = resp)
}

#' @rdname brmsformula-helpers
#' @export
lf <- function(..., flist = NULL, dpar = NULL, resp = NULL) {
  out <- c(list(...), flist)
  if (!is.null(dpar)) {
    dpar <- as_one_character(dpar)
    for (i in seq_along(out)) {
      attr(out[[i]], "dpar") <- dpar
    }
  }
  if (!is.null(resp)) {
    resp <- as_one_character(resp)
  }
  structure(out, dpar = dpar, resp = resp)
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
    stop2("No objects passed to 'mvbf'.")
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

split_bf <- function(x) {
  # build a mvbrmsformula object based on a brmsformula object
  # which uses cbind on the left-hand side to specify MV models
  stopifnot(is.brmsformula(x))
  resp <- parse_resp(x$formula, check_names = FALSE)
  str_adform <- get_matches(
    "\\|[^~]*(?=~)", formula2str(x$formula), perl = TRUE
  )
  if (length(resp) > 1L) {
    # cbind syntax used to specify MV model
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
  families <- ulapply(parts, function(f) f$family$family)
  all(families == "gaussian") || all(families == "student")
}

#' @export
"+.bform" <- function(e1, e2) {
  if (is.brmsformula(e1)) {
    out <- plus_brmsformula(e1, e2)
  } else if (is.mvbrmsformula(e1)) {
    out <- plus_mvbrmsformula(e1, e2)
  } else {
    stop2("Method '+.bform' not implemented for ", class(e1), " objects.")
  }
  out
}
  
plus_brmsformula <- function(e1, e2) {
  if (is.function(e2)) {
    e2 <- try(e2(), silent = TRUE)
    if (!is.family(e2)) {
      stop2("Don't know how to handle non-family functions.")
    }
  } 
  if (is.family(e2)) {
    e1 <- bf(e1, family = e2)
  } else if (is.cor_brms(e2)) {
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
  } else if (is.brmsformula(e2)) {
    e1 <- mvbf(e1, e2)
  } else if (inherits(e2, "setrescor")) {
    stop2("Setting 'rescor' is only possible in multivariate models.")
  } else {
    e1 <- bf(e1, e2)
  }
  e1
}

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
    e1$rescor <- e2
  } else if (is.brmsformula(e2)) {
    e1 <- mvbf(e1, e2)
  } else {
    resp <- attr(e2, "resp", TRUE)
    if (is.null(resp)) {
      stop2(
        "Don't know how to add a ", class(e2), " object ",
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

get_nl <- function(x, dpar = NULL, resp = NULL, aol = TRUE) {
  # extract the 'nl' attribute from a (brms)formula object
  # Args:
  #   x: object to extract 'nl' from
  #   dpar: optional name of a distributional parameter
  #     for which 'nl' should be extracted
  #   resp: optional name of a response variable
  #     for which 'nl' should be extracted
  #   aol: (as one logical) apply isTRUE to the result?
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

prepare_auxformula <- function(formula, par = NULL, rsv_pars = NULL) {
  # validate and prepare a formula of an distributional parameter
  # Args:
  #   formula: an object of class formula
  #   par: optional name of the parameter; if not specified
  #        the parameter name will be inferred from the formula
  #   rsv_pars: optional character vector of reserved parameter names
  stopifnot(length(par) <= 1L)
  try_formula <- try(as.formula(formula), silent = TRUE)
  if (is(try_formula, "try-error")) {
    if (length(formula) != 1L) {
      stop2("Expecting a single value when fixing parameters.")
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

dpars <- function() {
  # names of distributional parameters
  c("mu", "sigma", "shape", "nu", "phi", "kappa", "beta", "xi",
    "zi", "hu", "zoi", "coi", "disc", "bs", "ndt", "bias", 
    "quantile", "alpha", "theta")
}

links_dpars <- function(dp) {
  # link functions for distributional parameters
  switch(dp,
    mu = "identity",
    sigma = c("log", "identity"), 
    shape = c("log", "identity"),
    nu = c("logm1", "identity"), 
    phi = c("log", "identity"),
    kappa = c("log", "identity"), 
    beta = c("log", "identity"),
    zi = c("logit", "identity"), 
    hu = c("logit", "identity"),
    zoi = c("logit", "identity"), 
    coi = c("logit", "identity"), 
    disc = c("log", "identity"),
    bs = c("log", "identity"), 
    ndt = c("log", "identity"),
    bias = c("logit", "identity"),
    quantile = c("logit", "identity"),
    xi = c("log1p", "identity"),
    alpha = c("identity", "log"),
    theta = c("identity"), 
    stop2("Parameter '", dp, "' is not supported.")
  )
}

valid_dpars <- function(family, ...) {
  # get valid auxiliary parameters for a family
  UseMethod("valid_dpars")
}

#' @export
valid_dpars.default <- function(family, bterms = NULL, ...) {
  if (is.family(family) && !is.null(family[["dpars"]])) {
    # some families have dpar names that depend on the data
    return(family$dpars)
  }
  x <- c(
    mu = TRUE,
    sigma = has_sigma(family),
    shape = has_shape(family), 
    nu = has_nu(family, bterms = bterms), 
    phi = has_phi(family),
    kappa = has_kappa(family),
    beta = has_beta(family),
    zi = is_zero_inflated(family, zi_beta = TRUE), 
    hu = is_hurdle(family, zi_beta = FALSE),
    zoi = is_zero_one_inflated(family),
    coi = is_zero_one_inflated(family),
    bs = is_wiener(family), 
    ndt = is_wiener(family), 
    bias = is_wiener(family), 
    disc = is_ordinal(family),
    quantile = is_asym_laplace(family),
    alpha = has_alpha(family),
    xi = has_xi(family)
  )
  names(x)[x]
}

#' @export
valid_dpars.mixfamily <- function(family, ...) {
  out <- lapply(family$mix, valid_dpars, ...)
  for (i in seq_along(out)) {
    out[[i]] <- paste0(out[[i]], i)
  }
  c(unlist(out), paste0("theta", seq_along(out)))
}

#' @export
valid_dpars.brmsterms <- function(family, ...) {
  valid_dpars(family$family, bterms = family, ...)
}

is_dpar_name <- function(dpars, family = NULL, ...) {
  # check if provided names of auxiliar parameters are valid
  # Args:
  #   dpars: character vector to be checked
  #   family: the model family
  #   ...: further arguments passed to valid_dpars
  dpars <- as.character(dpars)
  if (!length(dpars)) {
    return(logical(0))
  }
  if (is.null(family)) {
    patterns <- paste0("^", dpars(), "[[:digit:]]*$")
    .is_dpar_name <- function(dpar, ...) {
      any(ulapply(patterns, grepl, x = dpar))
    }
    out <- ulapply(dpars, .is_dpar_name)
  } else {
    out <- dpars %in% valid_dpars(family, ...)
  }
  as.logical(out)
}

dpar_class <- function(dpar) {
  # class of a distributional parameter
  out <- get_matches("^[^[:digit:]]+", dpar, simplify = FALSE)
  ulapply(out, function(x) ifelse(length(x), x, ""))
}

dpar_id <- function(dpar) {
  # id of a distributional parameter
  out <- get_matches("[[:digit:]]+$", dpar, simplify = FALSE)
  ulapply(out, function(x) ifelse(length(x), x, ""))
}

validate_formula <- function(formula, ...) {
  # incorporate additional arguments into the model formula
  UseMethod("validate_formula")
}

#' @export
validate_formula.default <- function(formula, ...) {
  validate_formula(bf(formula), ...)
}

#' @export
validate_formula.brmsformula <- function(
  formula, family = gaussian(), autocor = cor_empty(), 
  data = NULL, threshold = NULL, ...
) {
  # incorporate additional arguments into the model formula
  # Args:
  #   formula: object of class 'formula' of 'brmsformula'
  #   data: optional data.frame
  #   family: optional object of class 'family'
  #   autocor: optional object of class 'cor_brms'
  #   threshold: (deprecated); threshold type for ordinal models
  # Returns:
  #   a brmsformula object compatible with the current version of brms
  out <- bf(formula)
  if (is.null(out$family)) {
    out$family <- check_family(family)
  }
  if (is.null(out$autocor)) {
    out$autocor <- check_autocor(autocor)
  }
  # allow the '.' symbol in the formulas
  out$formula <- expand_dot_formula(out$formula, data)
  for (i in seq_along(out$pforms)) {
    out$pforms[[i]] <- expand_dot_formula(out$pforms[[i]], data)
  }
  if (is_ordinal(out$family)) {
    if (is.null(out$family$threshold) && !is.null(threshold)) {
      # slot 'threshold' is deprecated as of brms > 1.7.0
      out$family <- check_family(out$family, threshold = threshold)
    }
    try_terms <- try(stats::terms(out$formula), silent = TRUE)
    intercept <- attr(try_terms, "intercept", TRUE)
    if (!is(try_terms, "try-error") && isTRUE(intercept == 0)) {
      stop2("Cannot remove the intercept in an ordinal model.")
    }
  }
  needs_cat <- is_categorical(out$family) && is.null(out$family$cats)
  if (needs_cat && !is.null(data)) {
    respform <- formula2str(lhs(out$formula))
    respform <- formula(gsub("\\|+[^~]*~", "~", respform))
    model_response <- model.response(model.frame(respform, data))
    cats <- levels(factor(model_response))
    if (length(cats) <= 2L) {
      stop2("At least 3 response categories are required for family ", 
            "'categorical'.\nPlease use family 'bernoulli' instead.")
    }
    # the first level will serve as the reference category
    out$family$dpars <- make.names(paste0("mu", cats[-1]), unique = TRUE)
    out$family$dpars <- gsub("\\.|_", "", out$family$dpars)
    if (any(duplicated(out$family$dpars))) {
      stop2("Invalid response category names. Please avoid ",
            "using any special characters in the names.")
    }
  }
  bf(out)
}

#' @export
validate_formula.mvbrmsformula <- function(
  formula, family = NULL, autocor = NULL, ...
) {
  # incorporate additional arguments into the MV model formula
  # allow passing lists of families or autocors
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
    formula$rescor <- allow_rescor
    message(
      "Setting 'rescor' to ", formula$rescor, 
      " by default for this combination of families"
    )
  }
  formula$rescor <- as_one_logical(formula$rescor)
  if (formula$rescor) {
    if (!allow_rescor) {
      stop2("Currently, estimating 'rescor' is only possible ", 
            "in multivariate gaussian or student models.")
    }
  }
  formula
}

#' @export
update.brmsformula <- function(object, formula., 
                               mode = c("update", "replace", "keep"), 
                               ...) {
  # update a brmsformula and / or its attributes
  # Args:
  #   object: an object of class 'brmsformula'
  #   formula.: formula to update object
  #   mode: "update": apply update.formula
  #         "replace": replace old formula
  #         "keep": keep old formula
  #         attributes are always updated
  #   ...: currently unused
  # Returns:
  #   a brmsformula object
  mode <- match.arg(mode)
  object <- bf(object)
  up_family <- formula.[["family"]]
  if (is.null(up_family)) {
    up_family <- object[["family"]]
  }
  up_autocor <- formula.[["autocor"]]
  if (is.null(up_autocor)) {
    up_autocor <- object[["autocor"]]
  }
  up_nl <- get_nl(formula, aol = FALSE)
  if (is.null(up_nl)) {
    up_nl <- get_nl(object)
  }
  # already use up_nl here to avoid ordinary parsing of NL formulas
  formula. <- bf(formula., nl = up_nl)
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

#' @export
print.brmsformula <- function(x, wsp = 0, digits = 2, ...) {
  cat(formula2str(x$formula, space = "trim"), "\n")
  str_wsp <- collapse(rep(" ", wsp))
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
