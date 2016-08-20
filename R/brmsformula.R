#' Set up a model formula for use in the \pkg{brms} package
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
#'   The details of model specification are given under 'Details'.
#' @param ... Additional \code{formula} objects to specify 
#'   predictors of special model parts and auxiliary parameters. 
#'   Formulas can either be named directly or contain
#'   names on their left-hand side. Currently, the following
#'   names are accepted: 
#'   \code{sigma} (residual standard deviation of
#'   the \code{gaussian} and \code{student} families);
#'   \code{shape} (shape parameter of the \code{Gamma},
#'   \code{weibull}, \code{negbinomial} and related
#'   zero-inflated / hurdle families); \code{nu}
#'   (degrees of freedom parameter of the \code{student} family);
#'   \code{phi} (precision parameter of the \code{beta} 
#'   and \code{zero_inflated_beta} families);
#'   \code{zi} (zero-inflation probability); 
#'   \code{hu} (hurdle probability).
#'   All auxiliary parameters are modeled 
#'   on the log or logit scale to ensure correct definition
#'   intervals after transformation.
#' @inheritParams brm
#' 
#' @return An object of class \code{brmsformula}, which inherits
#'  from class \code{formula} but contains additional attributes.
#'   
#' @details The \code{formula} argument accepts formulae of the following syntax:
#'   
#'   \code{response | addition ~ Pterms + (Gterms | group)} 
#'   
#'   The \code{Pterms} part contains effects that are assumed to be the 
#'   same across obervations. We call them 'population-level' effects
#'   or (adopting frequentist vocabulary) 'fixed' effects. The optional
#'   \code{Gterms} part may contain effects that are assumed to vary
#'   accross grouping variables specified in \code{group}. We
#'   call them 'group-level' effects or (adopting frequentist 
#'   vocabulary) 'random' effects, although the latter name is misleading
#'   in a Bayesian context (for more details type \code{vignette("brms")}).
#'   Multiple grouping factors each with multiple group-level effects 
#'   are possible. Instead of \code{|} you may use \code{||} in grouping terms
#'   to prevent correlations from being modeled. 
#'   Alternatively, it is possible to model different group-level terms of 
#'   the same grouping factor as correlated (even across different formulae,
#'   e.g. in non-linear models) by using \code{|<ID>|} instead of \code{|}.
#'   All group-level terms sharing the same ID will be modeled as correlated.
#'   If, for instance, one specifies the terms \code{(1+x|2|g)} and 
#'   \code{(1+z|2|g)} somewhere in the formulae passed to \code{brmsformula},
#'   correlations between the corresponding group-level effects 
#'   will be estimated.
#'   
#'   Smoothing terms can modeled using the \code{\link[mgcv:s]{s}}
#'   and \code{\link[mgcv:t2]{t2}} functions of the \pkg{mgcv} package 
#'   in the \code{Pterms} part of the model formula.
#'   This allows to fit generalized additive mixed models (GAMMs) with \pkg{brms}. 
#'   The implementation is similar to that used in the \pkg{gamm4} package.
#'   For more details on this model class see \code{\link[mgcv:gam]{gam}} 
#'   and \code{\link[mgcv:gamm]{gamm}}.
#'   
#'   The \code{Pterms} part may contain two non-standard types
#'   of population-level effects namely monotonic and category specific effects,
#'   which can be specified using terms of the form \code{monotonic(<predictors>)} 
#'   and \code{cse(<predictors>)} respectively. The latter can only be applied in
#'   ordinal models and is explained in more detail in the package's vignette
#'   (type \code{vignette("brms")}). The former effect type is explained here.
#'   A monotonic predictor must either be integer valued or an ordered factor, 
#'   which is the first difference to an ordinary continuous predictor. 
#'   More importantly, predictor categories (or integers) are not assumend to be 
#'   equidistant with respect to their effect on the response variable. 
#'   Instead, the distance between adjacent predictor categories (or integers) 
#'   is estimated from the data and may vary across categories. 
#'   This is realized by parameterizing as follows: 
#'   One parameter takes care of the direction and size of the effect similar 
#'   to an ordinary regression parameter, while an additional parameter vector 
#'   estimates the normalized distances between consecutive predictor categories.     
#'   A main application of monotonic effects are ordinal predictors that
#'   can this way be modeled without (falsely) treating them as continuous
#'   or as unordered categorical predictors.
#'   
#'   The third exception is the optional \code{addition} term, which may contain 
#'   multiple terms of the form \code{fun(variable)} seperated by \code{+} each 
#'   providing special information on the response variable. \code{fun} can be 
#'   replaced with either \code{se}, \code{weights}, \code{disp}, \code{trials},
#'   \code{cat}, \code{cens}, or \code{trunc}. Their meanings are explained below. 
#'   
#'   For families \code{gaussian}, \code{student}, and \code{cauchy} it is 
#'   possible to specify standard errors of the observation, thus allowing 
#'   to perform meta-analysis. Suppose that the variable \code{yi} contains 
#'   the effect sizes from the studies and \code{sei} the corresponding 
#'   standard errors. Then, fixed and random effects meta-analyses can 
#'   be conducted using the formulae \code{yi | se(sei) ~ 1} and 
#'   \code{yi | se(sei) ~ 1 + (1|study)}, respectively, where 
#'   \code{study} is a variable uniquely identifying every study.
#'   If desired, meta-regression can be performed via 
#'   \code{yi | se(sei) ~ 1 + mod1 + mod2 + (1|study)} 
#'   or \cr \code{yi | se(sei) ~ 1 + mod1 + mod2 + (1 + mod1 + mod2|study)}, where
#'   \code{mod1} and \code{mod2} represent moderator variables. 
#'   
#'   For all families, weighted regression may be performed using
#'   \code{weights} in the addition part. Internally, this is 
#'   implemented by multiplying the log-posterior values of each 
#'   observation by their corresponding weights.
#'   Suppose that variable \code{wei} contains the weights 
#'   and that \code{yi} is the response variable. 
#'   Then, formula \code{yi | weights(wei) ~ predictors} 
#'   implements a weighted regression. 
#'   
#'   The addition argument \code{disp} (short for dispersion) serves a
#'   similar purpose than \code{weight}. However, it has a different 
#'   implementation and is less general as it is only usable for the
#'   families \code{gaussian}, \code{student}, \code{cauchy},
#'   \code{lognormal}, \code{Gamma}, \code{weibull}, and \code{negbinomial}.
#'   For the former four families, the residual standard deviation 
#'   \code{sigma} is multiplied by the values given in 
#'   \code{disp}, so that higher values lead to lower weights.
#'   Contrariwise, for the latter three families, the parameter \code{shape}
#'   is multiplied by the values given in \code{disp}. As \code{shape}
#'   can be understood as a precision parameter (inverse of the variance),
#'   higher values will lead to higher weights in this case.
#'   
#'   For families \code{binomial} and \code{zero_inflated_binomial}, 
#'   addition should contain a variable indicating the number of trials 
#'   underlying each observation. In \code{lme4} syntax, we may write for instance 
#'   \code{cbind(success, n - success)}, which is equivalent
#'   to \code{success | trials(n)} in \code{brms} syntax. If the number of trials
#'   is constant across all observation (say \code{10}), 
#'   we may also write \code{success | trials(10)}. 
#'   
#'   For all ordinal families, \code{addition} may contain a term 
#'   \code{cat(number)} to specify the number categories (e.g, \code{cat(7)}). 
#'   If not given, the number of categories is calculated from the data.
#'   
#'   With the expection of \code{categorical} and ordinal families, 
#'   left and right censoring can be modeled through 
#'   \code{yi | cens(censored) ~ predictors}.
#'   The censoring variable (named \code{censored} in this example) should 
#'   contain the values \code{'left'}, \code{'none'}, and \code{'right'}  
#'   (or equivalenty -1, 0, and 1) to indicate that the corresponding observation is 
#'   left censored, not censored, or right censored. 
#'   
#'   With the expection of \code{categorical} and ordinal families, the response 
#'   distribution can be truncated using the \code{trunc} function in the addition part.
#'   If the response variable is truncated between, say, 0 and 100, we can specify this via
#'   \code{yi | trunc(lb = 0, ub = 100) ~ predictors}. 
#'   Instead of numbers, variables in the data set can also be passed allowing 
#'   for varying truncation points across observations. 
#'   Defining only one of the two arguments in \code{trunc} leads to one-sided truncation.
#' 
#'   Mutiple \code{addition} terms may be specified at the same time using 
#'   the \code{+} operator, for instance \code{formula = yi | se(sei) + cens(censored) ~ 1} 
#'   for a censored meta-analytic model. \cr
#'   
#'   For families \code{gaussian}, \code{student}, and \code{cauchy} 
#'   multivariate models may be specified using \code{cbind} notation. 
#'   In \pkg{brms} 1.0.0, the multvariate 'trait' syntax was removed 
#'   from the package as it repeatedly confused users, required much 
#'   special case coding, and was hard to maintain. Below the new 
#'   syntax is described. 
#'   Suppose that \code{y1} and \code{y2} are response variables 
#'   and \code{x} is a predictor. 
#'   Then \code{cbind(y1,y2) ~ x} specifies a multivariate model,
#'   The effects of all terms specified at the RHS of the formula 
#'   are assumed to vary across response variables (this was not the
#'   case by default in \pkg{brms} < 1.0.0). For instance, two parameters will
#'   be estimated for \code{x}, one for the effect
#'   on \code{y1} and another for the effect on \code{y2}.
#'   This is also true for group-level effects. When writing, for instance,
#'   \code{cbind(y1,y2) ~ x + (1+x|g)}, group-level effects will be
#'   estimated separately for each response. To model these effects
#'   as correlated across responses, use the ID syntax (see above).
#'   For the present example, this would look as follows:
#'   \code{cbind(y1,y2) ~ x + (1+x|2|g)}. Of course, you could also use
#'   any value other than \code{2} as ID. It is not yet possible
#'   to model terms as only affecting certain responses (and not others),
#'   but this will comebe implemented in the future.
#'    
#'   Categorical models use the same syntax as multivariate
#'   models. As in most other implementations of categorical models,
#'   values of one category (the first in \pkg{brms}) are fixed 
#'   to identify the model. Thus, all terms on the RHS of 
#'   the formula correspond to \code{K - 1} effects 
#'   (\code{K} = number of categories), one for each non-fixed category.
#'   Group-level effects may be specified as correlated across
#'   categories using the ID syntax.
#'   
#'   As of \pkg{brms} 1.0.0, zero-inflated and hurdle models are specfied 
#'   in the same way as as their non-inflated counterparts. 
#'   However, they have additional auxiliary parameters 
#'   (named \code{zi} and \code{hu} respectively)
#'   modeling the zero-inflation / hurdle probability depending on which 
#'   model you choose. These parameters can also be affected by predictors
#'   in the same way the response variable itself. See the end of the
#'   Details section for information on how to accomplish that.
#'   
#'   \bold{Parameterization of the population-level intercept}
#'   
#'   The population-level intercept (if incorporated) is estimated separately 
#'   and not as part of population-level parameter vector \code{b}. 
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
#'   \code{y ~ 0 + intercept + x}. This way, priors can be
#'   defined on the real intercept, directly. In addition,
#'   the intercept is just treated as an ordinary fixed effect
#'   and thus priors defined on \code{b} will also apply to it. 
#'   Note that this parameterization may be a bit less efficient
#'   than the default parameterization discussed above.  
#'   
#'   \bold{Formula syntax for non-linear models}
#'   
#'   Using the \code{nonlinear} argument, it is possible to specify
#'   non-linear models in \pkg{brms}. Contrary to what the name might suggest,
#'   \code{nonlinear} should not contain the non-linear model itself
#'   but rather information on the non-linear parameters. 
#'   The non-linear model will just be specified within the \code{formula}
#'   argument. Suppose, that we want to predict the response \code{y}
#'   through the predictor \code{x}, where \code{x} is linked to \code{y}
#'   through \code{y = alpha - beta * lambda^x}, with parameters
#'   \code{alpha}, \code{beta}, and \code{lambda}. This is certainly a
#'   non-linear model being defined via
#'   \code{formula = y ~ alpha - beta * lambda^x} (addition arguments 
#'   can be added in the same way as for ordinary formulas).
#'   Now we have to tell \pkg{brms} the names of the non-linear parameters 
#'   and specfiy a (linear mixed) model for each of them using the \code{nonlinear}
#'   argument. Let's say we just want to estimate those three parameters
#'   with no further covariates or random effects. Then we can write
#'   \code{nonlinear = alpha + beta + lambda ~ 1} or equivalently
#'   (and more flexible) \code{nonlinear = list(alpha ~ 1, beta ~ 1, lambda ~ 1)}. 
#'   This can, of course, be extended. If we have another predictor \code{z} and 
#'   observations nested within the grouping factor \code{g}, we may write for 
#'   instance \cr \code{nonlinear = list(alpha ~ 1, beta ~ 1 + z + (1|g), lambda ~ 1)}.
#'   The formula syntax described above applies here as well.
#'   In this example, we are using \code{z} and \code{g} only for the 
#'   prediction of \code{beta}, but we might also use them for the other
#'   non-linear parameters (provided that the resulting model is still 
#'   scientifically reasonable). 
#'   
#'   Non-linear models may not be uniquely identified and / or show bad convergence.
#'   For this reason it is mandatory to specify priors on the non-linear parameters.
#'   For instructions on how to do that, see \code{\link[brms:set_prior]{set_prior}}.
#'   
#'   \bold{Formula syntax for predicting auxiliary parameters}
#'   
#'   It is also possible to predict auxiliary parameters of the response
#'   distribution such as the residual standard deviation \code{sigma} 
#'   in gaussian models or the hurdle probability \code{hu} in hurdle models. 
#'   The syntax closely resembles that of a non-linear 
#'   parameter, for instance \code{sigma ~ x + s(z) + (1+x|g)}.
#'   
#'   All auxiliary parameters currently supported by \code{brmsformula}
#'   have to positive (a negative standard deviation or precision parameter 
#'   doesn't make any sense) or are bounded between 0 and 1 (for zero-inflated / 
#'   hurdle proabilities). 
#'   However, linear predictors can be positive or negative, and thus
#'   the log link (for positive parameters) or logit link (for probability parameters) 
#'   are used to ensure that auxiliary parameters are within their valid intervals.
#'   This implies that effects for auxiliary parameters are estimated on the
#'   log / logit scale and one has to apply the inverse link function to get 
#'   to the effects on the original scale.
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
#' # define a non-linear model
#' bf(y ~ a1 - a2^x, nonlinear = list(a1 ~ 1, a2 ~ x + (x|g)))
#' 
#' # correlated group-level effects across parameters
#' bf(y ~ a1 - a2^x, nonlinear = list(a1 ~ 1 + (1|2|g), a2 ~ x + (x|2|g)))
#' 
#' # define a multivariate model
#' bf(cbind(y1, y2) ~ x * z + (1|g))
#' 
#' # define a zero-inflated model 
#' # also predicting the zero-inflation part
#' bf(y ~ x * z + (1+x|ID1|g), zi ~ x + (1|ID1|g))
#' 
#' # specify a predictor as monotonic
#' bf(y ~ mono(x) + more_predictors)
#' 
#' # specify a predictor as category specific
#' # for ordinal models only
#' bf(y ~ cse(x) + more_predictors)
#' 
#' @export
brmsformula <- function(formula, ..., nonlinear = NULL) {
  # parse and validate dots arguments
  dots <- list(...)
  for (i in seq_along(dots)) {
    dots[[i]] <- prepare_auxformula(dots[[i]], par = names(dots)[i])
    names(dots)[i] <- attr(dots[[i]], "par")
    attr(dots[[i]], "par") <- NULL
  }
  if (any(!nzchar(names(dots)))) {
    stop("'brmsformula' requires named arguments.", call. = FALSE)
  }
  invalid_names <- setdiff(names(dots), auxpars())
  if (length(invalid_names)) {
    stop("The following argument names were invalid: ",
         paste(invalid_names, collapse = ", "), call. = FALSE)
  }
  # add attributes to formula
  if (is.logical(attr(formula, "nonlinear"))) {
    # In brms < 0.10.0 the nonlinear attribute was used differently
    attr(formula, "nonlinear") <- NULL
  }
  auxpars <- auxpars(incl_nl = TRUE)
  old_att <- rmNULL(attributes(formula)[auxpars])
  formula <- as.formula(formula)
  nonlinear <- nonlinear2list(nonlinear, rsv_pars = names(dots))
  new_att <- rmNULL(c(nlist(nonlinear), dots))
  dupl_args <- intersect(names(new_att), names(old_att))
  if (length(dupl_args)) {
    warning("Duplicated definitions of arguments ", 
            paste0("'", dupl_args, "'", collapse = ", "),
            "\nIgnoring definitions outside the formula",
            call. = FALSE)
  }
  null_pars <- setdiff(auxpars, names(old_att))
  new_pars <- intersect(names(new_att), null_pars)
  att <- c(old_att, new_att[new_pars])
  attributes(formula)[names(att)] <- att
  class(formula) <- c("brmsformula", "formula")
  formula
}

#' @export
bf <- function(formula, ..., nonlinear = NULL) {
  # alias of brmsformula
  brmsformula(formula, ..., nonlinear = nonlinear)
}

nonlinear2list <- function(x, rsv_pars = NULL) {
  # convert a single formula into a list of formulas
  # one for each non-linear parameter
  # Args:
  #   x: a formula or a list of formulas
  #   rsv_pars: optional character vector of reserved parameter names
  if (!(is.list(x) || is.null(x))) {
    x <- as.formula(x)
  }
  if (is(x, "formula")) {
    if (length(x) != 3L) {
      stop("Non-linear formulas must be two-sided.", call. = FALSE)
    }
    nlpars <- all.vars(lhs(x))
    x <- lapply(nlpars, function(nlp) update(x, paste(nlp, " ~ .")))
  }
  for (i in seq_along(x)) {
    x[[i]] <- prepare_auxformula(x[[i]], par = names(x)[i],
                                 rsv_pars = rsv_pars)
    names(x)[i] <- attr(x[[i]], "par")
    attr(x[[i]], "par") <- NULL
  }
  x
}

prepare_auxformula <- function(formula, par = NULL, rsv_pars = NULL) {
  # validate and prepare a formula of an auxiliary parameter
  # Args:
  #   formula: an object of class formula
  #   par: optional name of the parameter; if not specified
  #        the parameter name will be inferred from the formula
  #   rsv_pars: optional character vector of reserved parameter names
  stopifnot(length(par) <= 1L)
  formula <- as.formula(formula)
  if (!is.null(lhs(formula))) {
    resp_pars <- all.vars(formula[[2]])
    if (length(resp_pars) != 1L) {
      stop("LHS of additional formulas must contain exactly one variable.",
           call. = FALSE)
    }
    par <- resp_pars
    formula[[2]] <- eval(parse(text = paste("quote(", par, ")")))
  } else {
    if (!isTRUE(nzchar(par))) {
      stop("Additional formulas must be named.", call. = FALSE)
    }
    formula <- formula(paste(par, formula2string(formula)))
  }
  if (any(ulapply(c(".", "_"), grepl, x = par, fixed = TRUE))) {
    stop("Parameter names should not contain dots or underscores.",
         call. = FALSE)
  }
  if (par %in% rsv_pars) {
    stop("Parameter name '", par, "' is reserved for this model.",
         call. = FALSE)
  }
  if (!is.null(attr(terms(formula), "offset"))) {
    stop("Offsets in additional formulas are currently not allowed.", 
         call. = FALSE)
  }
  structure(formula, par = par)
}

auxpars <- function(incl_nl = FALSE) {
  # names of auxiliary parameters
  auxpars <- c("sigma", "shape", "nu", "phi", "zi", "hu")
  if (incl_nl) {
    auxpars <- c(auxpars, "nonlinear")
  }
  auxpars
}

sformula <- function(x, incl_nl = TRUE, flatten = FALSE, ...) {
  # extract special formulas stored in brmsformula objects
  # Args:
  #   x: coerced to a 'brmsformula' object
  #   incl_nl: include the 'nonlinear' argument in the output?
  #   flatten: should formulas stored in lists be put on the first
  #            level of the returned list?
  out <- rmNULL(attributes(bf(x))[auxpars(incl_nl = incl_nl)])
  if (flatten) {
    nonlinear <- out$nonlinear
    out$nonlinear <- NULL
    out <- c(out, nonlinear)
  }
  out
}
