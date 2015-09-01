# Extract fixed and random effects from a formula
# 
# @param formula An object of class "formula" using mostly the syntax of the \code{lme4} package
# @param ... Additional objects of class "formula"
# 
# @return A named list of the following six objects: \cr
#   \code{fixed}:    An object of class "formula" that contains the fixed effects including the dependent variable. \cr
#   \code{random}:   A list of formulas containing the random effects per grouping variable. \cr
#   \code{group}:    A list of names of the grouping variables. \cr
#   \code{add}:      A one sided formula containing the \code{add} part of \code{formula = y | add ~ predictors} if present. \cr
#   \code{add2}:  A one sided formula containing the \code{add2} part of \code{formula = y || add2 ~ predictors} if present. \cr
#   \code{all}:      A formula that contains every variable mentioned in \code{formula} and \code{...} 
# 
# @examples
# \dontrun{ 
# # fixed effects model
# extract_effects(response ~ I(1/a) + b)
# 
# # mixed effects model
# extract_effects(response ~ I(1/a) + b + (1 + c | d))
# 
# # mixed effects model with additional information on the response variable 
# # in this case standard errors in a gaussian linear model
# extract_effects(response | se(sei) ~ I(1/a) + b + (1 + c | d), family = "gaussian")
# }
extract_effects <- function(formula, ..., family = "none") {
  formula <- formula2string(formula)  
  fixed <- gsub(paste0("\\([^(\\||~)]*\\|[^\\)]*\\)\\+|\\+\\([^(\\||~)]*\\|[^\\)]*\\)",
                       "|\\([^(\\||~)]*\\|[^\\)]*\\)"),"",formula)
  fixed <- gsub("\\|+[^~]*~", "~", fixed)
  if (substr(fixed, nchar(fixed), nchar(fixed)) == "~") fixed <- paste0(fixed, "1")
  fixed <- formula(fixed)
  if (family %in% c("cumulative", "sratio", "cratio", "acat"))
    fixed <- update.formula(fixed, . ~ . +1)
  if (length(fixed) < 3) stop("invalid formula: response variable is missing")
  
  rg <- unlist(regmatches(formula, gregexpr("\\([^\\|\\)]*\\|[^\\)]*\\)", formula)))
  random <- lapply(regmatches(rg, gregexpr("\\([^\\|]*", rg)), function(r) 
    formula(paste0("~ ",substr(r, 2, nchar(r)))))
  cor <- unlist(lapply(regmatches(rg, gregexpr("\\|[^\\)]*", rg)), function(g) substr(g, 1, 2) != "||"))
  group_formulas <- lapply(regmatches(rg, gregexpr("\\|[^\\)]*", rg)), function(g) {
    g <- ifelse(substr(g, 1, 2) == "||", substr(g, 3, nchar(g)), substr(g, 2, nchar(g)))
    if (nchar(gsub(":", "", gsub("[^([:digit:]|[:punct:])][[:alnum:]_\\.]*", "", g))))
      stop(paste("Illegal grouping term:",g,"\nGrouping terms may contain only variable names",
                 "combined by the interaction symbol ':'"))
    return(formula(paste("~",g)))})
  group <- unlist(lapply(group_formulas, function(g) paste0(all.vars(g), collapse = ":")))
  # ordering is to ensure that all REs of the same grouping factor are next to each other
  x <- list(fixed = fixed, 
            random = if (length(group)) random[order(group)] else random, 
            cor = if (length(group)) cor[order(group)] else cor,
            group = if (length(group)) group[order(group)] else group)
  
  fun <- c("se", "weights", "trials", "cat", "cens")
  add_vars <- list()
  if (family != "none") {
    add <- unlist(regmatches(formula, gregexpr("\\|[^~]*~", formula)))[1]
    add <- substr(add, 2, nchar(add)-1)
    families <- list(se = c("gaussian","student","cauchy"), weights = c("all"),
                     trials = c("binomial"), cat = c("categorical", "cumulative", "cratio", "sratio", "acat"), 
                     cens = c("gaussian","student","cauchy","binomial","poisson","geometric","negbinomial","exponential",
                              "weibull","gamma"))
    for (f in fun) {
      x[[f]] <- unlist(regmatches(add, gregexpr(paste0(f,"\\([^\\|]*\\)"), add)))[1]
      add <- gsub(paste0(f,"\\([^~|\\|]*\\)\\|*"), "", add)
      if (is.na(x[[f]])) x[[f]] <- NULL
      else if (family %in% families[[f]] || families[[f]][1] == "all") {
        args <- substr(x[[f]], nchar(f) + 2, nchar(x[[f]]) -1)
        if (is.na(suppressWarnings(as.numeric(args)))) {
          x[[f]] <- as.formula(paste0("~ .", x[[f]]))
          add_vars[[f]] <- as.formula(paste("~", paste(all.vars(x[[f]]), collapse = "+")))
        }  
        else x[[f]] <- as.numeric(args)
      }  
      else stop(paste("Argument",f,"in formula is not supported by family",family))
    }
    if (nchar(gsub("\\|", "", add)) > 0 && !is.na(add))
      stop(paste0("Invalid addition part of formula. Please see the 'Details' section of help(brm) ",
                  "for further information. \nNote that the syntax of addition has changed in brms 0.2.1 as ",
                  "the old one was not flexible enough."))
  }
  
  new_formula <- unlist(lapply(c(random, group_formulas, add_vars, ...), 
                               function(x) paste0("+", Reduce(paste, deparse(x[[2]])))))
  new_formula <- paste0("update(",Reduce(paste, deparse(fixed)),", . ~ .",paste0(new_formula, collapse=""),")")
  x$all <- eval(parse(text = new_formula))
  environment(x$all) <- globalenv()
  x$response = all.vars(x$all[[2]])
  if (length(x$response) > 1) {
    if (!is.null(x$cens) || !is.null(x$se))
      stop("multivariate models currently allow only weights as addition arguments")
    x$fixed <- eval(parse(text = paste0("update(x$fixed, ", x$response[1], " ~ .)"))) 
    x$all <- eval(parse(text = paste0("update(x$all, ", x$response[1], " ~ .)"))) 
  }  
  x
} 

# extract time and grouping variabels for correlation structure
# 
# @formula a one sided formula of the form ~ time|group typically taken from a cor_brms object
# 
# @return a list with elements time, group, and all, where all contains a formula with all variables in formula
extract_time <- function(formula) {
  if (is.null(formula)) return(NULL)
  formula <- gsub(" ","",Reduce(paste, deparse(formula))) 
  time <- all.vars(as.formula(paste("~", gsub("~|\\|[[:print:]]*", "", formula))))
  if (length(time) > 1) stop("Autocorrelation structures may only contain 1 time variable")
  x <- list(time = ifelse(length(time), time, ""), group = "")
  
  group <- gsub("~[^\\|]*|\\|", "", formula)
  if (nchar(group)) {
    if (nchar(gsub(":", "", gsub("[^([:digit:]|[:punct:])][[:alnum:]_\\.]*", "", group))))
      stop(paste("Illegal grouping term:",group,"\nGrouping terms may contain only variable names",
                 "combined by the interaction symbol ':'"))
    group <- formula(paste("~", group))
    x$group <- paste0(all.vars(group), collapse = ":")
  }
  x$all <- formula(paste("~",paste(c("1", time, all.vars(group)), collapse = "+")))
  x
}

# incorporate addition and partial arguments into formula 
# 
# @param formula a model formula 
# @param addition a list with one sided formulas taken from the addition arguments in brm
# @param partial a one sided formula containing partial effects
#
# @return an updated formula containing the addition and partial effects
update_formula <- function(formula, addition = NULL, partial = NULL) {
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

# check validity of family
check_family <- function(family) {
  if (family == "normal") family <- "gaussian"
  if (family == "multigaussian") 
    stop("family 'multigaussian' is deprecated. Use family 'gaussian' instead")
  if (!family %in% c("gaussian", "student", "cauchy", "binomial", "bernoulli", "categorical",
                     "poisson", "negbinomial", "geometric", "gamma", "weibull", "exponential",
                     "cumulative", "cratio", "sratio", "acat"))
    stop(paste(family, "is not a valid family"))
  family
}

# Links for \code{brms} families
# 
# @param family A vector of one or two character strings. The first string indicates the distribution of the dependent variable (the 'family'). Currently, the following distributions are supported:
#  \code{"gaussian"}, \code{"student"}, \code{"cauchy"}, \code{"poisson"}, \code{"binomial"}, \code{"categorical"}, 
#  \code{"gamma"}, \code{"exponential"}, \code{"weibull"}, \code{"cumulative"}, \cr
#  \code{"cratio"}, \code{"sratio"}, and \code{"acat"}.
#  The second string indicates the link function, which must supported by the distribution of the dependent variable. If not specified, default link functions are used (see 'Details').
# @return The second element of \code{family} (if present) or else the default link of the specified family
#   
# @details The families \code{gaussian}, \code{student}, and \code{cauchy} accept the links (as names) \code{identity}, \code{log}, and \code{inverse};
# the \code{poisson} family the links \code{log}, \code{identity}, and \code{sqrt}; 
# families \code{binomial}, \code{cumulative}, \code{cratio}, \code{sratio}, and \code{acat} the links \code{logit}, \code{probit}, \code{probit_approx}, and \code{cloglog};
# family  \code{categorical} the link \code{logit}; families \code{gamma}, \code{weibull}, and \code{exponential} the links \code{log}, \code{identity}, and \code{inverse}. 
# The first link mentioned for each family is the default.     
# 
# @examples brm.link("gaussian")
# brm.link(c("gaussian","log"))
link4family <- function(family) {
  link <- family[2]
  family <- family[1]
  is_linear <- family %in% c("gaussian", "student", "cauchy")
  is_skew <- family %in% c("gamma", "weibull", "exponential")
  is_cat <- family %in% c("cumulative", "cratio", "sratio", "acat", "binomial", "bernoulli")                    
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

# list irrelevant parameters not to be saved by Stan
# 
# @param formula a model formula
# @param ranef logical; should random effects parameters of each levels be saved?
#
# @return a vector of parameters to be excluded
exclude_pars <- function(formula, ranef = TRUE) {
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

# check prior input by and amend it if needed
#
# @param same as the respective parameters in brm
#
# @return a list of prior specifications adjusted to be used in stan_prior (see stan.R)
check_prior <- function(prior, formula, data = NULL, family = "gaussian", autocor = NULL, 
                        partial = NULL, threshold = "flexible") {
  
  #expand lkj correlation prior to full name
  prior <- lapply(prior, function(p) sub("^lkj\\(", "lkj_corr_cholesky(", p))
  
  #check if parameter names in prior are correct
  ee <- extract_effects(formula, family = family)  
  possible_priors <- unlist(parnames(formula, data = data, family = family, autocor = autocor,
                                     partial = partial, threshold = threshold, internal = TRUE), use.names = FALSE)
  meta_priors <- unlist(regmatches(possible_priors, gregexpr("^[^_]+", possible_priors)))
  if ("sd" %in% meta_priors)
    meta_priors <- c(meta_priors, paste0("sd_",ee$group))
  possible_priors <- unique(c(possible_priors, meta_priors))
  wrong_priors <- names(prior)[!names(prior) %in% possible_priors]
  if (length(wrong_priors))
    warning(paste("Some parameter names in prior cannot be found in the model:", 
                  paste0(wrong_priors, collapse = ", ")))
  
  #rename certain parameters
  names(prior) <- rename(names(prior), symbols = c("^cor_", "^cor$", "^rescor$"), 
                         subs = c("L_", "L", "Lrescor"), fixed = FALSE)
  if (any(grepl("^sd_.+", names(prior))))
    for (i in 1:length(ee$group)) 
      names(prior) <-  rename(names(prior), symbols = paste0("^sd_",ee$group[[i]]),
                              subs = paste0("sd_",i), fixed = FALSE)
  if (family %in% c("cumulative", "sratio", "cratio", "acat") && threshold == "equidistant")
    names(prior) <- rename(names(prior), symbols = "^b_Intercept$", subs = "b_Intercept1",
                           fixed = FALSE)
  prior
}

#' Extract parameter names
#' 
#' Extract all parameter names for which priors may be specified
#' 
#' @param x An object of class \code{formula}
#' @inheritParams brm
#' @param internal A flag indicating if the names of additional internal parameters should be displayed. 
#'   Setting priors on these parameters is not recommended
#' @param ... Currently ignored
#' 
#' @return A list of character vectors containing the parameter names for which priors may be specified
#' 
#' @examples 
#' parnames(rating ~ treat + period + carry + (1+carry|subject), 
#'           data = inhaler, family = "student")
#'           
#' parnames(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|visit),
#'           data = epilepsy, family = "poisson")          
#' 
#' @export
parnames.formula <- function(x, data = NULL, family = "gaussian", addition = NULL, autocor = NULL, 
                             partial = NULL, threshold = "flexible", internal = FALSE, ...) {
  if (is.null(autocor)) autocor <- cor_arma()
  if (!is(autocor, "cor_brms")) stop("cor must be of class cor_brms")
  if (!threshold %in% c("flexible","equidistant")) 
    stop("threshold must be either flexible or equidistant")
  family <- check_family(family[1])
  x <- update_formula(x, addition = addition)
  ee <- extract_effects(x, partial, family = family)
  data <- update_data(data, family = family, effects = ee)
  out <- list(fixef = paste0("b_",colnames(get_model_matrix(ee$fixed, data = data))),
              ranef = list(), other = NULL)
  if (is.formula(partial)) 
    out$fixef <- c(out$fixef, paste0("b_",colnames(get_model_matrix(partial, data = data, rm_intercept = TRUE))))
  if (length(ee$group)) {
    gs <- unlist(ee$group)
    for (i in 1:length(gs)) {
      ranef <- colnames(get_model_matrix(ee$random[[i]], data = data))
      out$ranef[[gs[i]]] <- sort(unique(c(out$ranef[[gs[i]]], paste0("sd_",gs[i],"_",ranef),
                                          if (ee$cor[[i]] && length(ranef) > 1) 
                                            c(paste0("cor_",gs[i]), if(internal) paste0("L_",gs[i]))))) 
    }
  }
  if (is(autocor, "cor_arma") && autocor$p) out$other <- c(out$other, "ar")
  if (is(autocor, "cor_arma") && autocor$q) out$other <- c(out$other, "ma")
  if (family %in% c("gaussian", "student", "cauchy") && !is.formula(ee$se))
    out$other <- c(out$other, paste0("sigma_",ee$response))
  if (family == "gaussian" && length(ee$response) > 1)
    out$other <- c(out$other, "rescor", if (internal) "Lrescor")
  if (family == "student") out$other <- c(out$other, "nu")
  if (family %in% c("gamma", "weibull", "negbinomial")) 
    out$other <- c(out$other, "shape")
  if (family %in% c("cumulative", "sratio", "cratio", "acat") && threshold == "equidistant")
    out$other <- c(out$other, "delta")
  out
}