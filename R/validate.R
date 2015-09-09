extract_effects <- function(formula, ..., family = "none") {
  # Extract fixed and random effects from a formula
  # 
  # Args:
  #   formula: An object of class "formula" using mostly the syntax of the \code{lme4} package
  #   ...: Additional objects of class "formula"
  #   family: the model family
  # 
  # Returns: 
  #   A named list of the following elements: 
  #   fixed: An object of class "formula" that contains the fixed effects including the dependent variable. 
  #   random: A list of formulas containing the random effects per grouping variable. 
  #   group: A vector of names of the grouping variables. 
  #   weights, se, cens, trials, cat: information on possible addition arguments
  #   all: A formula that contains every variable mentioned in formula and ...
  formula <- formula2string(formula)  
  fixed <- gsub(paste0("\\([^(\\||~)]*\\|[^\\)]*\\)\\+|\\+\\([^(\\||~)]*\\|[^\\)]*\\)",
                       "|\\([^(\\||~)]*\\|[^\\)]*\\)"), "", formula)
  fixed <- gsub("\\|+[^~]*~", "~", fixed)
  if (substr(fixed, nchar(fixed), nchar(fixed)) == "~") 
    fixed <- paste0(fixed, "1")
  fixed <- formula(fixed)
  if (family %in% c("cumulative", "sratio", "cratio", "acat"))
    fixed <- update.formula(fixed, . ~ . +1)
  if (length(fixed) < 3) 
    stop("invalid formula: response variable is missing")
  
  # extract random effects part
  rg <- unlist(regmatches(formula, gregexpr("\\([^\\|\\)]*\\|[^\\)]*\\)", formula)))
  random <- lapply(regmatches(rg, gregexpr("\\([^\\|]*", rg)), 
                   function(r) formula(paste0("~ ",substr(r, 2, nchar(r)))))
  cor <- unlist(lapply(regmatches(rg, gregexpr("\\|[^\\)]*", rg)), 
                       function(g) substr(g, 1, 2) != "||"))
  group <- regmatches(rg, gregexpr("\\|[^\\)]*", rg))
  group_formula <- lapply(group, get_group_formula)
  group <- unlist(lapply(group_formula, function(g) 
                         paste0(all.vars(g), collapse = ":")))
  
  # ordering is to ensure that all REs of the same grouping factor are next to each other
  x <- list(fixed = fixed, 
            random = if (length(group)) random[order(group)] else random, 
            cor = if (length(group)) cor[order(group)] else cor,
            group = if (length(group)) group[order(group)] else group)
  
  # handle addition arguments
  fun <- c("se", "weights", "trials", "cat", "cens")
  add_vars <- list()
  if (family != "none") {
    add <- unlist(regmatches(formula, gregexpr("\\|[^~]*~", formula)))[1]
    add <- substr(add, 2, nchar(add)-1)
    families <- list(se = c("gaussian","student","cauchy"),
                     weights = "all",
                     trials = "binomial", 
                     cat = c("categorical", "cumulative", "cratio", "sratio", "acat"), 
                     cens = c("gaussian","student","cauchy","binomial","poisson",
                              "geometric","negbinomial","exponential", "weibull","gamma"))
    for (f in fun) {
      x[[f]] <- unlist(regmatches(add, gregexpr(paste0(f,"\\([^\\|]*\\)"), add)))[1]
      add <- gsub(paste0(f,"\\([^~|\\|]*\\)\\|*"), "", add)
      if (is.na(x[[f]])) {
        x[[f]] <- NULL
      } else if (family %in% families[[f]] || families[[f]][1] == "all") {
        args <- substr(x[[f]], nchar(f) + 2, nchar(x[[f]]) -1)
        if (is.na(suppressWarnings(as.numeric(args)))) {
          x[[f]] <- as.formula(paste0("~ .", x[[f]]))
          add_vars[[f]] <- as.formula(paste("~", paste(all.vars(x[[f]]), collapse = "+")))
        } else {
          x[[f]] <- as.numeric(args)
        }
      } else {
        stop(paste("Argument", f, "in formula is not supported by family", family))
      }
    }
    if (nchar(gsub("\\|", "", add)) > 0 && !is.na(add))
      stop(paste0("Invalid addition part of formula. Please see the 'Details' section of help(brm)"))
  }
  
  # make a formula containing all required variables (element 'all')
  new_formula <- unlist(lapply(c(random, group_formula, add_vars, ...), 
                               function(x) paste0("+", Reduce(paste, deparse(x[[2]])))))
  new_formula <- paste0("update(",Reduce(paste, deparse(fixed)),
                        ", . ~ .",paste0(new_formula, collapse=""),")")
  x$all <- eval(parse(text = new_formula))
  environment(x$all) <- globalenv()
  
  # extract response variables
  x$response = all.vars(x$all[[2]])
  if (length(x$response) > 1) {
    if (!is.null(x$cens) || !is.null(x$se))
      stop("multivariate models currently allow only weights as addition arguments")
    x$fixed <- eval(parse(text = paste0("update(x$fixed, ", x$response[1], " ~ .)"))) 
    x$all <- eval(parse(text = paste0("update(x$all, ", x$response[1], " ~ .)"))) 
  }  
  x
} 

extract_time <- function(formula) {
  # extract time and grouping variabels for correlation structure
  # 
  # Args:
  #   formula: a one sided formula of the form ~ time|group typically taken from a cor_brms object
  # 
  # Returns: 
  #   a list with elements time, group, and all, where all contains a formula with all variables in formula
  if (is.null(formula)) 
    return(NULL)
  formula <- gsub(" ","",Reduce(paste, deparse(formula))) 
  time <- all.vars(as.formula(paste("~", gsub("~|\\|[[:print:]]*", "", formula))))
  if (length(time) > 1) 
    stop("Autocorrelation structures may only contain 1 time variable")
  x <- list(time = ifelse(length(time), time, ""))
  group <- get_group_formula(sub("~[^\\|]*", "", formula))
  x$group <- paste0(all.vars(group), collapse = ":")
  x$all <- formula(paste("~",paste(c("1", time, all.vars(group)), collapse = "+")))
  x
}

update_formula <- function(formula, addition = NULL, partial = NULL) {
  # incorporate addition and partial arguments into formula 
  # 
  # Args:
  #   formula: a model formula 
  #   addition: a list with one sided formulas taken from the addition arguments in brm
  #   partial: a one sided formula containing partial effects
  #
  # Returns:
  #   an updated formula containing the addition and partial effects
  var_names <- names(addition)
  addition <- lapply(addition, formula2string, rm = 1)
  fnew <- "."
  if (length(addition)) 
    for (i in 1:length(addition))
      fnew <- paste0(fnew, " | ", var_names[i], "(", addition[[i]], ")")
  fnew <- paste(fnew, "~ .")
  if (is.formula(partial)) {
    partial <- formula2string(partial, rm=1)
    fnew <- paste(fnew, "+ partial(", partial, ")")
  }
  update.formula(formula, formula(fnew))
}

get_group_formula <- function(g) {
  # transform grouping term in formula
  # 
  # Args: 
  #   g: a grouping term 
  #
  # Returns:
  #   the formula ~ g if g is valid and else an error
  g <- sub("^\\|*", "", g)
  if (nchar(gsub(":|[^([:digit:]|[:punct:])][[:alnum:]_\\.]*", "", g)))
    stop(paste("Illegal grouping term:", g, "\n",
               "may contain only variable names combined by the symbol ':'"))
  if (nchar(g)) {
    return(formula(paste("~", g)))
  } else {
    return(~1)
  }
}

gather_ranef <- function(effects, data = NULL) {
  # gathers helpful information on the random effects
  #
  # Args:
  #   effects: output of extract_effects
  #   data: data passed to brm after updating
  #
  # Returns: 
  #   A named list with one element per grouping factor
  Z <- lapply(effects$random, get_model_matrix, data = data)
  ranef <- setNames(lapply(Z, colnames), effects$group)
  if (length(ranef)) {
    for (i in 1:length(ranef)) {
      attr(ranef[[i]], "levels") <- 
        levels(as.factor(get(effects$group[[i]], data)))  
    }
  }
  ranef
}

check_family <- function(family) {
  # check validity of model family
  if (family == "normal")
    family <- "gaussian"
  if (family == "multigaussian") 
    stop("family 'multigaussian' is deprecated. Use family 'gaussian' instead")
  if (!family %in% c("gaussian", "student", "cauchy", "binomial", 
                     "bernoulli", "categorical", "poisson", 
                     "negbinomial", "geometric", "gamma", 
                     "weibull", "exponential", "cumulative", 
                     "cratio", "sratio", "acat"))
    stop(paste(family, "is not a valid family"))
  family
}

link4family <- function(family) {
  # check validity of the link function and return default links if no link is specified
  link <- family[2]
  family <- family[1]
  is_linear <- family %in% c("gaussian", "student", "cauchy")
  is_skew <- family %in% c("gamma", "weibull", "exponential")
  is_cat <- family %in% c("cumulative", "cratio", "sratio", 
                          "acat", "binomial", "bernoulli")                    
  is_count <- family %in% c("poisson", "negbinomial", "geometric")
  
  if (is.na(link)) {
    if (is_linear) link <- "identity"
    else if (is_skew || is_count) link <- "log"
    else if (is_cat || family == "categorical") link <- "logit"
  }
  else if (is_linear && !link %in% c("identity", "log", "inverse") ||
           is_count && !link %in% c("log", "identity", "sqrt") ||
           is_cat && !link %in% c("logit", "probit", "probit_approx", "cloglog") ||
           family == "categorical" && link != "logit" ||
           is_skew && !link %in% c("log", "identity", "inverse"))
    stop(paste(link, "is not a valid link for family", family))
  else if (is_count && link == "sqrt") 
    warning(paste(family, "model with sqrt link may not be uniquely identified"))
  link
}

exclude_pars <- function(formula, ranef = TRUE) {
  # list irrelevant parameters not to be saved by Stan
  # 
  # Args:
  #   formula: a model formula
  #   ranef: logical; should random effects parameters of each levels be saved?
  #
  # Returns:
  #   a vector of parameters to be excluded
  ee <- extract_effects(formula)
  out <- c("eta", "etam", "etap", "b_Intercept1", "Lrescor", "Rescor",
           "p", "q", "e", "Ema", "lp_pre")
  if (length(ee$group)) {
    for (i in 1:length(ee$group)) {
      out <- c(out, paste0("pre_",i), paste0("L_",i), paste0("Cor_",i))
      if (!ranef) out <- c(out, paste0("r_",i))
    }
  }
  out
}

remove_chains <- function(i, sflist) {
  # remove chains that produce errors leaving the other chains untouched
  #
  # Args:
  #   i: an index between 1 and length(sflist) 
  #   sflist: list of stanfit objects as returned by parLapply
  if (!is(sflist[[i]], "stanfit") || length(sflist[[i]]@sim$samples) == 0) {
    warning(paste("chain", i, "did not contain samples and was removed from the fitted model"))
    return(NULL)  
  } else {
    return(sflist[[i]])
  }
}

#' Prior Deinitions for \pkg{brms} Models
#'
#' Define priors for specific parameters or classes of parameters
#'
#' @param prior A character string defining a distribution in \pkg{Stan} language
#' @param class The parameter class. Defaults to \code{"b"} (fixed effects). 
#'   See 'Details' for other valid parameter classes. 
#' @param coef Name of the (fixed, partial, or random effects) parameter  
#' @param group Grouping factor for random effects parameters.
#' 
#' @return An object of class \code{brmsprior} to be used in the \code{prior}
#'   argument of \code{\link[brms:brm]{brm}}.
#' 
#' @details 
#'   \code{set_prior} is used to define prior distributions for parameters in \pkg{brms} models.
#'   Below, we explain its usage and list some common prior distributions for parameters. 
#'   A complete overview on possible prior distributions is given in the Stan Reference Manual available at 
#'   \url{http://mc-stan.org/}.
#'   
#'   \pkg{brms} performs no checks if the priors are written in correct Stan language.
#'   Instead, Stan will check their correctness when the model is parsed to C++ and returns an error if they are not.
#'   Currently, there are five types of parameters in \pkg{brms} models, 
#'   for which the user can specify prior distributions. \cr
#'   
#'   1. Fixed and partial effects 
#'   
#'   Every fixed (and partial) effect has its corresponding regression parameter. These parameters are internally named as
#'   \code{b_<fixed>}, where \code{<fixed>} represents the name of the corresponding fixed effect. 
#'   Suppose, for instance, that \code{y} is predicted by \code{x1} and \code{x2} 
#'   (i.e. \code{y ~ x1+x2} in formula syntax). 
#'   Then, \code{x1} and \code{x2} have regression parameters \code{b_x1} and \code{b_x2} respectively. 
#'   The default prior for fixed and partial effects is an improper flat prior over the reals. 
#'   Other common options are normal priors or uniform priors over a finite interval.
#'   If we want to have a normal prior with mean 0 and standard deviation 5 for \code{x1}, 
#'   and a uniform prior between -10 and 10 for \code{x2}, we can specify this via \cr
#'   \code{set_prior("normal(0,5)", class = "b", coef = "x1")} and \cr
#'   \code{set_prior("uniform(-10,10)", class = "b", coef = "x2")}.
#'   To put the same prior on all fixed effects at once, 
#'   we may write as a shortcut \code{set_prior("<prior>", class = "b")}. This also
#'   leads to faster sampling, because priors can be vectorized in this case. \cr
#'   
#'   2. Autocorrelation parameters
#'   
#'   The autocorrelation parameters currently implemented are named \code{ar} (autoregression) and \code{ma} (moving average).
#'   The default prior for autocorrelation parameters is an improper flat prior over the reals. 
#'   Other priors can be defined with \code{set_prior("<prior>", class = "ar")} 
#'   or \cr \code{set_prior("<prior>", class = "ma")}. It should be noted that \code{ar} will
#'   only take one values between -1 and 1 if the response variable is wide-sence stationay, 
#'   i.e. if there is no drift in the responses. \cr
#'   
#'   3. Standard deviations of random effects
#'   
#'   Each random effect of each grouping factor has a standard deviation named
#'   \code{sd_<group>_<random>}. Consider, for instance, the formula \code{y ~ x1+x2+(1+x1|g)}.
#'   We see that the intercept as well as \code{x1} are random effects nested in the grouping factor \code{g}. 
#'   The corresponding standard deviation parameters are named as \code{sd_g_Intercept} and \code{sd_g_x1} respectively. 
#'   These parameters are restriced to be non-negative and, by default, 
#'   have a half cauchy prior with scale parameter 5. 
#'   We could make this explicit by writing \code{set_prior("cauchy(0,5)", class = "sd")}. 
#'   To define a prior distribution only for standard deviations of a specific grouping factor,
#'   use \cr \code{set_prior("<prior>", class = "sd", group = "<group>")}. 
#'   To define a prior distribution only for a specific standard deviation 
#'   of a specific grouping factor, you may write \cr
#'   \code{set_prior("<prior>", class = "sd", group = "<group>", coef = "<coef>")}. 
#'   Recommendations on useful prior distributions for standard deviations are given in Gelman (2006). \cr
#'   
#'   4. Correlations of random effects 
#'   
#'   If there is more than one random effect per grouping factor, the correlations between those random
#'   effects have to be estimated. 
#'   The prior \code{"lkj_corr_cholesky(eta)"} or in short \code{"lkj(eta)"} with \code{eta > 0} is essentially the only prior 
#'   for (choelsky factors) of correlation matrices. If \code{eta = 1} (the default) all correlations matrices 
#'   are equally likely a priori. If \code{eta > 1}, extreme correlations become less likely, 
#'   whereas \code{0 < eta < 1} results in higher probabilities for extreme correlations. 
#'   Correlation matrix parameters in \code{brms} models are named as 
#'   \code{cor_(group)}, (e.g., \code{cor_g} if \code{g} is the grouping factor).
#'   To set the same prior on every correlation matrix, use for instance \code{set_prior("lkj(2)", class = "cor")}.
#'   
#'   5. Parameters for specific families 
#'   
#'   Some families need additional parameters to be estimated. 
#'   Families \code{gaussian}, \code{student}, and \code{cauchy} need the parameter \code{sigma} 
#'   to account for the standard deviation of the response variable around the regression line
#'   (not to be confused with the standard deviations of random effects). 
#'   By default, \code{sigma} has a half cauchy prior with 'mean' 0 and 'standard deviation' 5. 
#'   Furthermore, family \code{student} needs the parameter \code{nu} representing 
#'   the degrees of freedom of students t distribution. 
#'   By default, \code{nu} has prior \code{"uniform(1,100)"}. 
#'   Families \code{gamma} and \code{weibull} need the parameter \code{shape} 
#'   that has a \code{"gamma(0.01,0.01)"} prior by default. For families \code{cumulative}, \code{cratio}, \code{sratio}, 
#'   and \code{acat}, and only if \code{threshold = "equidistant"}, the parameter \code{delta} is used to model the distance
#'   between to adjacent thresholds. By default, \code{delta} has an improper flat prior over the reals. \cr
#'   Every family specific parameter has its own prior class, so that \cr
#'   \code{set_prior("<prior>", class = "<parameter>")} it the right way to go. \cr
#' 
#'   Often, it may not be immediately clear, which parameters are present in the model.
#'   To get a full list of parameters and parameter classes for which priors can be specified (depending on the model) 
#'   use function \code{\link[brms:get_prior]{get_prior}}
#'
#' @seealso \code{\link[brms:get_prior]{get_prior}}
#' 
#' @references 
#' Gelman A (2006). Prior distributions for variance parameters in hierarchical models."
#'    Bayesian analysis, 1(3), 515 -- 534.
#' 
#' @examples
#' \dontrun{
#' ## check which parameters can have priors
#' get_prior(rating ~ treat + period + carry + (1|subject),
#'           data = inhaler, family = "cumulative", 
#'           threshold = "equidistant")
#'          
#' ## define some priors          
#' prior <- c(set_prior("normal(0,10)", class = "b"),
#'            set_prior("normal(1,2)", class = "b", coef = "treat"),
#'            set_prior("cauchy(0,2)", class = "sd", 
#'                      group = "subject", coef = "Intercept"),
#'            set_prior("uniform(-5,5)", class = "delta"))
#'               
#' ## use the defined priors in the model
#' fit <- brm(rating ~ period + carry + (1|subject),
#'            data = inhaler, family = "sratio", 
#'            partial = ~ treat, threshold = "equidistant",
#'            prior = prior, n.iter = 1000, n.cluster = 2)
#'            
#' ## check that the priors found their way into Stan's model code
#' fit$model             
#' }
#'
#' @export
set_prior <- function(prior, class = "b", coef = "", group = "") {
  prior <- as.character(prior)
  class <- as.character(class)
  group <- as.character(group)
  coef <- as.character(coef)
  if (length(prior) != 1 || length(class) != 1 
      || length(coef) != 1 || length(group) != 1)
    stop("All arguments of set_prior must be of length 1")
  valid_classes <- c("b", "sd", "cor", "L", "ar", "ma", "sigma", 
                     "rescor", "Lrescor", "nu", "shape", "delta")
  if (!class %in% valid_classes)
    stop(paste(class, "is not a valid paramter class"))
  if (nchar(group) && !class %in% c("sd", "cor", "L"))
    stop(paste("argument group not meaningful for class", class))
  if (nchar(coef) && !class %in% c("b", "sd", "sigma"))
    stop(paste("argument coef not meaningful for class", class))
  out <- list(prior = prior, class = class, coef = coef, group = group)
  class(out) <- c("brmsprior", "list")
  out
}

#' Overview on Priors for \pkg{brms} Models
#' 
#' Get information on all parameters (and parameter classes) for which priors may be specified including default priors.
#' 
#' @inheritParams brm
#' @param internal A flag indicating if the names of additional internal parameters should be displayed. 
#'   Setting priors on these parameters is not recommended
#' 
#' @return A data.frame with columns \code{prior}, \code{class}, \code{coef}, and \code{group}
#'   and several rows, each providing information on a paramter (or parameter class) on which
#'   priors can be specified. The prior column is empty except for internal default priors.
#'   
#' @seealso \code{\link[brms:set_prior]{set_prior}}
#' 
#' @examples 
#' \dontrun{
#' ## get all parameters and parameters classes to define priors on
#' (prior <- get_prior(count ~ log_Age_c + log_Base4_c * Trt_c
#'                     + (1|patient) + (1|visit),
#'                     data = epilepsy, family = "poisson"))   
#'          
#' ## define a prior on all fixed effects a once
#' prior$prior[1] <- "normal(0,10)"
#' 
#' ## define a specific prior on the fixed effect of Trt_c
#' prior$prior[5] <- "uniform(-5,5)"       
#' 
#' ## fit a model using the priors above
#' fit <- brm(count ~ log_Age_c + log_Base4_c * Trt_c 
#'            + (1|patient) + (1|visit),
#'            data = epilepsy, family = "poisson", prior = prior)
#'            
#' ## check that priors indeed found their way into Stan's model code
#' fit$model
#' }
#' 
#' @export
get_prior <- function(formula, data = NULL, family = "gaussian", addition = NULL, 
                      autocor = NULL, partial = NULL, threshold = c("flexible", "equidistant"), 
                      internal = FALSE) {
  # note that default priors are stored in this function
  if (is.null(autocor)) 
    autocor <- cor_arma()
  if (!is(autocor, "cor_brms")) 
    stop("cor must be of class cor_brms")
  threshold <- match.arg(threshold)
  family <- check_family(family[1])
  formula <- update_formula(formula, addition = addition)
  ee <- extract_effects(formula, partial, family = family)
  data <- update_data(data, family = family, effects = ee)
  
  # initialize output
  prior <- prior_frame(prior = character(0), class = character(0), 
                     coef = character(0), group = character(0))
  # fixed and partial effects
  fixef <- colnames(get_model_matrix(ee$fixed, data = data))
  if (length(fixef)) {
    prior <- rbind(prior, prior_frame(class = "b", coef = c("", fixef)))
  }
  if (is.formula(partial)) {
    paref <- colnames(get_model_matrix(partial, data = data, rm_intercept = TRUE))
    prior <- rbind(prior, prior_frame(class = "b", coef = paref))
  }
  # random effects
  if (length(ee$group)) {
    prior <- rbind(prior, prior_frame(class = "sd", prior = "cauchy(0,5)"))  # global sd class
    gs <- unlist(ee$group)
    for (i in 1:length(gs)) {
      ranef <- colnames(get_model_matrix(ee$random[[i]], data = data))
      # include random effects standard deviations
      prior <- rbind(prior, prior_frame(class = "sd", coef = c("", ranef), group = gs[i]))
      # detect duplicated random effects
      J <- with(prior, class == "sd" & group == gs[i] & nchar(coef))
      dupli <- duplicated(prior[J, ])
      if (any(dupli)) {
        stop(paste("Duplicated random effects detected for group", gs[i]))
      }
      # include correlation parameters
      if (ee$cor[[i]] && length(ranef) > 1) {
        if (internal) {
          prior <- rbind(prior, prior_frame(class = "L", group = c("", gs[i]),
                                            prior = c("lkj_corr_cholesky(1)", "")))
        } else {
          prior <- rbind(prior, prior_frame(class = "cor", group = c("", gs[i]),
                                            prior = c("lkj(1)", "")))
        }
      }
    }
  }
  # handle additional parameters
  is_ordinal <- family %in% c("cumulative", "sratio", "cratio", "acat") 
  if (is(autocor, "cor_arma") && autocor$p) 
    prior <- rbind(prior, prior_frame(class = "ar"))
  if (is(autocor, "cor_arma") && autocor$q) 
    prior <- rbind(prior, prior_frame(class = "ma"))
  if (family %in% c("gaussian", "student", "cauchy") && !is.formula(ee$se))
    prior <- rbind(prior, prior_frame(class = "sigma", coef = c("", ee$response),
                                      prior = c("cauchy(0,5)", rep("", length(ee$response)))))
  if (family == "gaussian" && length(ee$response) > 1) {
    if (internal) {
      prior <- rbind(prior, prior_frame(class = "Lrescor", prior = "lkj_corr_cholesky(1)"))
    } else {
      prior <- rbind(prior, prior_frame(class = "rescor", prior = "lkj(1)"))
    }
  }
  if (family == "student") 
    prior <- rbind(prior, prior_frame(class = "nu", prior = "uniform(1,100)"))
  if (family %in% c("gamma", "weibull", "negbinomial")) 
    prior <- rbind(prior, prior_frame(class = "shape", prior = "gamma(0.01,0.01)"))
  if (is_ordinal && threshold == "equidistant")
    prior <- rbind(prior, prior_frame(class = "delta"))
  prior <- unique(prior)
  prior <- prior[with(prior, order(class, group, coef)), ]
  rownames(prior) <- 1:nrow(prior)
  prior
}

check_prior <- function(prior, formula, data = NULL, family = "gaussian", 
                        autocor = NULL, partial = NULL, threshold = "flexible") {
  # check prior input and amend it if needed
  #
  # Args:
  #   same as the respective parameters in brm
  #
  # Returns:
  #   a data.frame of prior specifications to be used in stan_prior (see stan.R)
  ee <- extract_effects(formula, family = family)  
  all_prior <- get_prior(formula = formula, data = data, family = family, 
                         autocor = autocor, partial = partial, 
                         threshold = threshold, internal = TRUE)
  if (is.null(prior)) {
    return(all_prior)  # nothing to check
  } else if (is(prior, "brmsprior")) {
    # a single prior may be specified without c(.)
    prior <- c(prior)
  } else if (!is(prior, "prior_frame") && is.list(prior) && !is.null(names(prior))) {
    # deprecated prior specification brms < 0.5.0
    warning(paste("Specifying priors using a named list is deprecated. \n",
                  "We strongly recommend to use the set_prior function instead. \n",
                  "See help(set_prior) for further information."))
    prior <- update_prior(prior)
  } else if (!is(prior, "prior_frame")) {
    stop("Invalid input for argument prior. See help(set_prior) for further information.")
  }
  
  prior$class <- rename(prior$class, symbols = c("^cor$", "^rescor$"), 
                        subs = c("L", "Lrescor"), fixed = FALSE)
  duplicated_input <- duplicated(prior[, 2:4])
  if (any(duplicated_input)) {
    stop("Duplicated prior specifications are not allowed. \n")
  }
  # expand lkj correlation prior to full name
  prior$prior <- sub("^lkj\\(", "lkj_corr_cholesky(", prior$prior)
  
  # check if parameters in prior are valid
  valid <- which(duplicated(rbind(all_prior[, 2:4], prior[, 2:4])))
  invalid <- which(!1:nrow(prior) %in% (valid - nrow(all_prior)))
  if (length(invalid)) {
    message(paste("Prior element", paste(invalid, collapse = ", "),
                  "is invalid and will be removed."))
    prior <- prior[-invalid, ]
  }
  
  # merge prior with all_prior
  prior <- rbind(prior, all_prior)
  rm <- which(duplicated(prior[, 2:4]))
  if (length(rm))  # else it may happen that all rows a removed...
    prior <- prior[-rm, ]
  
  # rename parameter groups
  rows2remove <- NULL
  group_indices <- which(nchar(prior$group) > 0)
  for (i in group_indices) {
    if (!prior$group[i] %in% ee$group) { 
      stop(paste("grouping factor", prior$group[i], "not found in the model"))
    } else if (sum(prior$group[i] == ee$group) == 1) {
      # matches only one grouping factor in the model
      prior$group[i] <- match(prior$group[i], ee$group)
    } else {
      # matches multiple grouping factors in the model
      rows2remove <- c(rows2remove, i)
      which_match <- which(prior$group[i] == ee$group)
      new_rows <- lapply(which_match, function(j) {
        new_row <- prior[i, ]
        new_row$group <- j
        new_row
      })
      prior <- rbind(prior, do.call(rbind, new_rows))  # add new rows
    }
  }
  # get partial priors out of fixef priors
  if (family == "categorical" || is.formula(partial)) {
    paref <- colnames(get_model_matrix(partial, data = data, rm_intercept = TRUE))
    b_index <- which(prior$class == "b" & !nchar(prior$coef))
    partial_index <- which(prior$class == "b" & prior$coef %in% paref)
    rows2remove <- c(rows2remove, partial_index)
    partial_prior <- prior[c(b_index, partial_index), ]
    partial_prior$class <- "bp"  # the partial effects class
    prior <- rbind(prior, partial_prior)
  }
  # special treatment of thresholds in ordinal models
  if (family %in% c("cumulative", "sratio", "cratio", "acat")) {
    # take specific fixed effects Intercept prior
    Int_index <- which(prior$class == "b" & prior$coef == "Intercept")
    rows2remove <- c(rows2remove, Int_index)
    if (!length(Int_index)) {  # take global fixed effects prior
      Int_index <- which(prior$class == "b" & !nchar(prior$coef))
    }
    if (length(Int_index)) {
      Int_prior <- prior[Int_index, ] 
      # thresholds have their own internal parameter class
      Int_prior$class <- ifelse(threshold == "equidistant", 
                                "b_Intercept1", "b_Intercept")
      Int_prior$coef <- ""
      prior <- rbind(prior, Int_prior)
    }
  }
  # remove unnecessary rows
  if (length(rows2remove)) {   
    prior <- prior[-rows2remove, ]
  }
  prior <- prior[with(prior, order(class, group, coef)), ]
  rownames(prior) <- 1:nrow(prior)
  prior
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

prior_frame <- function(prior = "", class = "", coef = "", group = "") {
  # helper function to create data.frames containing prior information 
  out <- data.frame(prior = prior, class = class, coef = coef, group = group,
                    stringsAsFactors = FALSE)
  class(out) <- c("prior_frame", "data.frame")
  out
}

#' @export
print.brmsprior <- function(x, ...) {
  group <- ifelse(nchar(x$group), paste0("_", x$group), "")
  coef <- ifelse(nchar(x$coef), paste0("_", x$coef), "")
  cat(paste0("Prior: ", x$class, group, coef, " ~ ", x$prior))    
}

#' @export
c.brmsprior <- function(x, ...) {
  # combines multiple brmsprior objects into one prior_frame
  if(any(!sapply(list(...), is, class2 = "brmsprior")))
    stop("All arguments must be of class brmsprior")
  prior <- data.frame(matrix(unlist(list(x, ...)), ncol = 4, byrow = TRUE),
                      stringsAsFactors = FALSE)
  names(prior) <- c("prior", "class", "coef", "group") 
  class(prior) <- c("prior_frame", "data.frame")
  prior
}