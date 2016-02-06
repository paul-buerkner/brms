#' Prior Definitions for \pkg{brms} Models
#'
#' Define priors for specific parameters or classes of parameters
#'
#' @param prior A character string defining a distribution in \pkg{Stan} language
#' @param class The parameter class. Defaults to \code{"b"} (fixed effects). 
#'   See 'Details' for other valid parameter classes. 
#' @param coef Name of the (fixed, category specific, or random effects) parameter  
#' @param group Grouping factor for random effects parameters.
#' @param nlpar Name of a non-linear parameter. Only used in non-linear models.
#' @param lb Lower bound for parameter restriction. Currently only allowed
#'   if \code{class = "b"}. Defaults to \code{NULL}, that is no restriction.
#' @param ub Upper bound for parameter restriction. Currently only allowed
#'   if \code{class = "b"}. Defaults to \code{NULL}, that is no restriction.
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
#'   \pkg{brms} performs no checks if the priors are written in 
#'   correct Stan language. Instead, Stan will check their correctness 
#'   when the model is parsed to C++ and returns an error if they are not.
#'   Currently, there are five types of parameters in \pkg{brms} models, 
#'   for which the user can specify prior distributions. \cr
#'   
#'   1. Fixed and category specific effects 
#'   
#'   Every fixed (and category specific) effect has its own regression parameter. 
#'   These parameters are internally named as \code{b_<fixed>}, 
#'   where \code{<fixed>} represents 
#'   the name of the corresponding fixed effect. 
#'   Suppose, for instance, that \code{y} is predicted by \code{x1} and \code{x2} 
#'   (i.e. \code{y ~ x1+x2} in formula syntax). 
#'   Then, \code{x1} and \code{x2} have regression parameters 
#'   \code{b_x1} and \code{b_x2} respectively. 
#'   The default prior for fixed and category specific effects is an 
#'   improper flat prior over the reals. Other common options are normal priors 
#'   or student-t priors. If we want to have a normal prior with mean 0 and 
#'   standard deviation 5 for \code{x1}, 
#'   and a unit student-t prior with 10 degrees of freedom for \code{x2}, 
#'   we can specify this via
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
#'   of the other fixed effects.
#'   
#'   A special shrinkage prior to be applied on fixed effects is the horseshoe prior.
#'   It is symmetric around zero with fat tails and an infinitely large spike
#'   at zero. This makes it ideal for sparse models that have 
#'   many regression coefficients,although only a minority of them is non-zero. 
#'   For more details see Carvalho et al. (2009).
#'   The horseshoe prior can be applied on all fixed effects at once 
#'   (excluding the intercept) by using \code{set_prior("horseshoe(<df>)")}.
#'   Replace \code{<df>} with the desired degrees of freedom of the student-t prior
#'   of the local shrinkage parameters. 
#'   In their paper, Carvalho et al. (2009) use one degrees of freedom, but this
#'   my lead to an increased number of divergent transition in \pkg{Stan}
#'   so that slightly higher values may often be a better option. 
#'   Generally, models with horseshoe priors a more likely than other models
#'   to have divergent transitions so that increasing \code{adapt_delta} 
#'   from \code{0.8} to values closer to \code{1} will often be necessary.
#'   See the documentation of \code{\link[brms:brm]{brm}} for instructions
#'   on how to increase \code{adapt_delta}. \cr
#'   
#'   If desired, fixed effects parameters can be restricted to fall only 
#'   within a certain interval using the \code{lb} and \code{ub} arguments
#'   of \code{set_prior}. This is often required when defining priors
#'   that are not defined everywhere on the real line, such as uniform
#'   or gamma priors. When defining a \code{uniform(2,4)} prior, 
#'   you should write \code{set_prior("uniform(2,4)", lb = 2, ub = 4)}. 
#'   When using a prior that is defined on the postive reals only 
#'   (such as a gamma prior) set \code{lb = 0}. 
#'   In most situations, it is not useful to restrict fixed effects
#'   parameters through bounded priors, but if you really want to
#'   this is the way to go.
#'   
#'   3. Autocorrelation parameters
#'   
#'   The autocorrelation parameters currently implemented are named 
#'   \code{ar} (autoregression), \code{ma} (moving average),
#'   and \code{arr} (autoregression of the response).
#'   The default prior for autocorrelation parameters is an 
#'   improper flat prior over the reals. 
#'   Other priors can be defined by \cr
#'   \code{set_prior("<prior>", class = "ar")} 
#'   for \code{ar} effects and similar for \code{ma} and \code{arr} effects.
#'   
#'   4. Standard deviations of random effects
#'   
#'   Each random effect of each grouping factor has a standard deviation named
#'   \code{sd_<group>_<random>}. Consider, for instance, the formula 
#'   \code{y ~ x1+x2+(1+x1|g)}.
#'   We see that the intercept as well as \code{x1} are random effects 
#'   nested in the grouping factor \code{g}. 
#'   The corresponding standard deviation parameters are named as 
#'   \code{sd_g_Intercept} and \code{sd_g_x1} respectively. 
#'   These parameters are restriced to be non-negative and, by default, 
#'   have a half student-t prior with 3 degrees of freedom and a 
#'   scale parameter that depends on the standard deviation of the response 
#'   after applying the link function. Minimally, the scale parameter is 5. 
#'   To define a prior distribution only for standard deviations 
#'   of a specific grouping factor,
#'   use \cr \code{set_prior("<prior>", class = "sd", group = "<group>")}. 
#'   To define a prior distribution only for a specific standard deviation 
#'   of a specific grouping factor, you may write \cr
#'   \code{set_prior("<prior>", class = "sd", group = "<group>", coef = "<coef>")}. 
#'   Recommendations on useful prior distributions for 
#'   standard deviations are given in Gelman (2006). \cr
#'   
#'   5. Correlations of random effects 
#'   
#'   If there is more than one random effect per grouping factor, 
#'   the correlations between those random
#'   effects have to be estimated. 
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
#'   \code{"student_t(3,0,5)"} prior by default. 
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
#' make_stancode(rating ~ period + carry + (1|subject),
#'               data = inhaler, family = sratio(), 
#'               partial = ~ treat, threshold = "equidistant",
#'               prior = prior)
#'               
#' ## use horseshoe priors to model sparsity in fixed effects parameters
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
  if (length(prior) != 1 || length(class) != 1 || length(coef) != 1 || 
      length(group) != 1 || length(nlpar) != 1 || length(lb) > 1 || 
      length(ub) > 1)
    stop("All arguments of set_prior must be of length 1", call. = FALSE)
  valid_classes <- c("Intercept", "b", "sd", "cor", "L", "ar", "ma", "arr",
                     "sigma", "rescor", "Lrescor", "nu", "shape", "delta", "phi")
  if (!class %in% valid_classes)
    stop(paste(class, "is not a valid parameter class"), call. = FALSE)
  if (nchar(group) && !class %in% c("sd", "cor", "L"))
    stop(paste("argument group not meaningful for class", class), 
         call. = FALSE)
  if (nchar(coef) && !class %in% c("b", "sd", "sigma"))
    stop(paste("argument coef not meaningful for class", class))
  if (nchar(nlpar) && !class %in% valid_classes[1:5])
    stop(paste("argument nlpar not meaningful for class", class))
  if (length(lb) || length(ub)) {
    if (!class %in% c("b"))
      stop("currently boundaries are only allowed for fixed effects")
    if (coef != "")
      stop("coef may not be specified when using boundaries")
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
#' ## define a prior on all fixed effects a once
#' prior$prior[1] <- "normal(0,10)"
#' 
#' ## define a specific prior on the fixed effect of Trt_c
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
  formula <- update_formula(formula, data = data)
  family <- check_family(family) 
  link <- family$link
  threshold <- match.arg(threshold)
  autocor <- check_autocor(autocor)
  ee <- extract_effects(formula, partial, family = family,
                        nonlinear = nonlinear)
  data <- update_data(data, family = family, effects = ee)
  
  # ensure that RE and residual SDs only have a weakly informative prior by default
  Y <- unname(model.response(data))
  prior_scale <- 5
  if (link %in% c("identity", "log", "inverse", "sqrt", "1/mu^2")) {
    if (link %in% c("log", "inverse", "1/mu^2")) {
      # avoid Inf in link(Y)
      Y <- ifelse(Y == 0, Y + 0.1, Y)
    }
    suggested_scale <- round(sd(link(Y, link = link)))
    if (!is.nan(suggested_scale)) {
      prior_scale <- max(prior_scale, suggested_scale, na.rm = TRUE)
    } 
  }
  def_scale_prior <- paste0("student_t(3, 0, ", prior_scale, ")")
  
  # initialize output
  prior <- prior_frame(prior = character(0), class = character(0), 
                       coef = character(0), group = character(0),
                       nlpar = character(0), bound = character(0))
  if (length(nonlinear)) {
    prior <- rbind(prior, prior_frame(class = "b"))
    for (i in seq_along(nonlinear)) {
      fixef <- colnames(get_model_matrix(ee$nonlinear[[i]]$fixed, data = data))
      prior <- rbind(prior, prior_frame(class = "b", coef = c("", fixef), 
                                        nlpar = names(ee$nonlinear)[i])) 
    }
  } else {
    # fixed and category specific effects 
    fixef <- colnames(get_model_matrix(ee$fixed, data = data))
    if (length(fixef)) {
      if ("Intercept" %in% fixef) {
        prior <- rbind(prior, prior_frame(class = "Intercept"))
      }
      prior <- rbind(prior, prior_frame(class = "b", coef = c("", fixef))) 
    }
    if (is.formula(partial)) {
      paref <- colnames(get_model_matrix(partial, data = data, 
                                         rm_intercept = TRUE))
      prior <- rbind(prior, prior_frame(class = "b", coef = paref))
    }
  }
  # random effects
  random <- get_random(ee)
  if (nrow(random)) {
    # global sd class
    prior <- rbind(prior, prior_frame(class = "sd", prior = def_scale_prior))  
    gs <- random$group
    nlpars <- if (length(nonlinear)) colnames(random) 
              else rep("", ncol(random))
    for (i in seq_along(gs)) {
      ranef <- colnames(get_model_matrix(random$form[[i]], data = data))
      # include random effects standard deviations
      prior <- rbind(prior, prior_frame(class = "sd", coef = c("", ranef), 
                                        group = gs[i], nlpar = nlpars[i]))
      # detect duplicated random effects
      J <- with(prior, class == "sd" & group == gs[i] & 
                       nlpar == nlpars[i] & nchar(coef))
      dupli <- duplicated(prior[J, ])
      if (any(dupli)) {
        stop(paste("Duplicated random effects detected for group", gs[i]),
             call. = FALSE)
      }
      # include correlation parameters
      if (random$cor[[i]] && length(ranef) > 1) {
        if (internal) {
          prior <- rbind(prior, 
            prior_frame(class = "L", group = c("", gs[i]), nlpar = nlpars[i],
                        prior = c("lkj_corr_cholesky(1)", "")))
        } else {
          prior <- rbind(prior, 
            prior_frame(class = "cor", group = c("", gs[i]),
                        nlpar = nlpars[i], prior = c("lkj(1)", "")))
        }
      }
    }
  }
  # handle additional parameters
  is_ordinal <- is.ordinal(family)
  is_linear <- is.linear(family)
  nresp <- length(ee$response)
  if (get_ar(autocor)) 
    prior <- rbind(prior, prior_frame(class = "ar"))
  if (get_ma(autocor)) 
    prior <- rbind(prior, prior_frame(class = "ma"))
  if (get_arr(autocor)) 
    prior <- rbind(prior, prior_frame(class = "arr"))
  if (has_sigma(family, se = is.formula(ee$se), autocor = autocor)) {
    sigma_prior <- prior_frame(class = "sigma", coef = c("", ee$response),
                               prior = c(def_scale_prior, rep("", nresp)))
    prior <- rbind(prior, sigma_prior)
  }
  if (is_linear && nresp > 1) {
    if (internal) {
      prior <- rbind(prior, prior_frame(class = "Lrescor", 
                                        prior = "lkj_corr_cholesky(1)"))
    } else {
      prior <- rbind(prior, prior_frame(class = "rescor", prior = "lkj(1)"))
    }
  }
  if (family$family == "student") {
    prior <- rbind(prior, prior_frame(class = "nu", 
                                      prior = "gamma(2, 0.1)"))
  }
  if (family$family == "beta") {
    prior <- rbind(prior, prior_frame(class = "phi", 
                                      prior = "gamma(0.01, 0.01)"))
  }
  if (has_shape(family)) {
    prior <- rbind(prior, prior_frame(class = "shape", 
                                      prior = def_scale_prior))
  }
  if (is_ordinal && threshold == "equidistant") {
    prior <- rbind(prior, prior_frame(class = "delta"))
  }
  prior <- unique(prior)
  prior <- prior[with(prior, order(class, group, coef)), ]
  rownames(prior) <- 1:nrow(prior)
  prior
}

check_prior <- function(prior, formula, data = NULL, family = gaussian(), 
                        autocor = NULL,  nonlinear = NULL, partial = NULL, 
                        threshold = "flexible") {
  # check prior input and amend it if needed
  #
  # Args:
  #   same as the respective parameters in brm
  #
  # Returns:
  #   a data.frame of prior specifications to be used in stan_prior (see stan.R)
  if (isTRUE(attr(prior, "checked"))) {
    # prior has already been checked; no need to do it twice
    return(prior)
  }
  ee <- extract_effects(formula, family = family, nonlinear = nonlinear)  
  all_priors <- get_prior(formula = formula, data = data, 
                          family = family, autocor = autocor, 
                          partial = partial, threshold = threshold, 
                          internal = TRUE, nonlinear = nonlinear)
  if (is.null(prior)) {
    prior <- all_priors  
  } else if (is(prior, "brmsprior")) {
    # a single prior may be specified without c(.)
    prior <- c(prior)
  } else if (!is(prior, "prior_frame") && is.list(prior) 
             && !is.null(names(prior))) {
    # deprecated prior specification brms < 0.5.0
    warning(paste("Specifying priors using a named list is deprecated. \n",
                  "We strongly recommend to use the set_prior function instead. \n",
                  "See help(set_prior) for further information."))
    prior <- update_prior(prior)
  } else if (!is(prior, "prior_frame")) {
    stop(paste("Invalid input for argument prior. See help(set_prior)", 
               "for further information."), call. = FALSE)
  }
  
  # exclude prior using increment_log_prob to readd the at the end
  has_incr_lp <- grepl("^increment_log_prob\\(", prior$prior)
  prior_incr_lp <- prior[has_incr_lp, ]
  prior <- prior[!has_incr_lp, ]
  
  prior$class <- rename(prior$class, symbols = c("^cor$", "^rescor$"), 
                        subs = c("L", "Lrescor"), fixed = FALSE)
  duplicated_input <- duplicated(prior[, 2:5])
  if (any(duplicated_input)) {
    stop("Duplicated prior specifications are not allowed.", call. = FALSE)
  }
  
  # handle special priors that are not explictly coded as functions in Stan
  temp <- handle_special_priors(prior)  
  prior <- temp$prior
  attrib <- temp$attrib 
  
  # check if parameters in prior are valid
  if (nrow(prior)) {
    valid <- which(duplicated(rbind(all_priors[, 2:5], prior[, 2:5])))
    invalid <- which(!1:nrow(prior) %in% (valid - nrow(all_priors)))
    if (length(invalid)) {
      message(paste("Prior element", paste(invalid, collapse = ", "),
                    "is invalid and will be removed."))
      prior <- prior[-invalid, ]
    }
  }
  # merge prior with all_priors
  prior <- rbind(prior, all_priors)
  rm <- which(duplicated(prior[, 2:5]))
  if (length(rm)) { 
    # else it may happen that all rows a removed...
    prior <- prior[-rm, ]
  }
  
  rows2remove <- NULL
  # special treatment of fixed effects Intercept(s)
  Int_index <- which(prior$class == "Intercept")
  if (length(Int_index)) {
    # if an intercept is present
    rows2remove <- c(rows2remove, Int_index)
    Int_prior <- prior[Int_index, ]
    old_index <- which(prior$class == "b" & prior$coef == "Intercept")
    rows2remove <- c(rows2remove, old_index)
    if (length(old_index) && nchar(prior$prior[old_index])) {
      # for backwards compatibility
      Int_prior$prior <- prior$prior[old_index]
      warning(paste("Using class = 'b' with coef = 'Intercept' is deprecated.", 
                    "See help(set_prior) for further details."), call. = FALSE)
    }
    # (temporary) Intercepts have their own internal parameter class
    res_thres <- is.ordinal(family) && threshold == "equidistant"
    Int_prior$class <- ifelse(res_thres, "temp_Intercept1", "temp_Intercept")
    Int_prior$coef <- ""
    prior <- rbind(prior, Int_prior)
  }
  # get category specific priors out of fixef priors
  if (is.categorical(family) || is.formula(partial)) {
    paref <- colnames(get_model_matrix(partial, data = data, 
                                       rm_intercept = TRUE))
    b_index <- which(prior$class == "b" & !nchar(prior$coef))
    partial_index <- which(prior$class == "b" & prior$coef %in% paref)
    rows2remove <- c(rows2remove, partial_index)
    partial_prior <- prior[c(b_index, partial_index), ]
    partial_prior$class <- "bp"  # the category specific effects class
    prior <- rbind(prior, partial_prior)
  }
  # check if priors for non-linear parameters are defined
  if (length(nonlinear)) {
    nlpars <- names(ee$nonlinear)
    for (nlp in nlpars) {
      nlp_prior <- prior$prior[with(prior, nlpar == nlp & class == "b")]
      if (!any(as.logical(nchar(nlp_prior)))) {
        stop(paste0("Priors for non-linear parameters are required, ",
                    "but no prior found for parameter '", nlp, "'. \n",
                    "See help(set_prior) for more details."), 
             call. = FALSE)
      }
    }
  }
  # rename random effects priors to match names in Stan code
  random <- get_random(ee)
  group_indices <- which(nchar(prior$group) > 0)
  for (i in group_indices) {
    if (!prior$group[i] %in% random$group) { 
      stop(paste("grouping factor", prior$group[i], "not found in the model"),
           call. = FALSE)
    } else if (sum(prior$group[i] == random$group) == 1) {
      # matches only one grouping factor in the model
      prior$group[i] <- match(prior$group[i], random$group)
    } else {
      # matches multiple grouping factors in the model
      rows2remove <- c(rows2remove, i)
      which_match <- which(prior$group[i] == random$group)
      new_rows <- lapply(which_match, function(j) {
        new_row <- prior[i, ]
        new_row$group <- j
        new_row
      })
      prior <- rbind(prior, do.call(rbind, new_rows))  # add new rows
    }
  }
  # remove unnecessary rows
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
  attr(prior, "checked") <- TRUE
  prior
}

handle_special_priors <- function(prior) {
  # look for special priors such as horseshoe and process them appropriately
  #
  # Args:
  #   prior: an object of class prior_frame
  #
  # Returns:
  #   an named list of two objects: 
  #   prior: an updated version of prior
  #   attrib: a named list containing future attributes of prior
  attrib <- list()
  b_index <- which(prior$class == "b" & !nchar(prior$coef))
  if (length(b_index) && grepl("^horseshoe\\(.+\\)$", prior$prior[b_index])) {
    # horseshoe prior for fixed effects parameters
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
  prior$prior <- sub("^lkj\\(", "lkj_corr_cholesky(", prior$prior)
  list(prior = prior, attrib = attrib)
}

update_prior <- function(prior) {
  # update deprecated prior specifications from brms < 0.5.0
  #
  # Args:
  #   prior: A named list
  #
  # Returns:
  #   a data.frame compatible with check_prior of brms >= 0.5.0
  if (!is.list(prior) || is.null(names(prior))) {
    stop("Only named lists can be updated")
  }
  prior_names <- names(prior)
  class <- regmatches(prior_names, regexpr("^[^_]+", prior_names))
  
  # try to separate group from coef
  group_coef <- substr(prior_names, nchar(class) + 2, nchar(prior_names))
  group <- rep("", length(prior))
  for (i in 1:length(prior)) {
    if (class[i] == "sd" && prior_names[i] != "sd") {
      s <- regmatches(group_coef[i], regexpr("^[^_]+", group_coef[i]))
      group[i] <- ifelse(length(s), s, "")
    } else if (class[i] == "cor") {
      group[i] <- group_coef[i]
    }
  }
  coef <- substr(group_coef, nchar(group) + ifelse(nchar(group), 2, 1), 
                 nchar(group_coef))
  
  prior_frame(prior = unlist(prior, use.names = FALSE),
              class = class, coef = coef, group = group)
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

#' @export
print.brmsprior <- function(x, ...) {
  group <- ifelse(nchar(x$group), paste0("_", x$group), "")
  coef <- ifelse(nchar(x$coef), paste0("_", x$coef), "")
  nlpar <- ifelse(nchar(x$nlpar), paste0("_", x$nlpar), "")
  tilde <- ifelse(nchar(x$class) + nchar(group) + nchar(coef), " ~ ", "")
  bound <- ifelse(nchar(x$bound), paste0(x$bound, " "), "")
  cat(paste0("Prior: ", bound, x$class, nlpar, group, coef, tilde, x$prior))
  invisible(x)
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