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

prior_frame <- function(prior = "", class = "", coef = "", group = "") {
  # helper function to easily create data.frames containing prior information 
  data.frame(prior = prior, class = class, coef = coef, group = group,
             stringsAsFactors = FALSE)
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
#'   arguments of \code{\link[brms:brm]{brm}}.
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
#'   \code{set_prior("normal(0,5)", class = "b", coef = "x1")} and 
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
#'   or \code{set_prior("<prior>", class = "ma")}. It should be noted that \code{ar} will
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
#'   use \code{set_prior("<prior>", class = "sd", group = "<group>")}. 
#'   To define a prior distribution only for a specific standard deviation 
#'   of a specific grouping factor, you may write
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
#'   Every family specific parameter has its own prior class, so that 
#'   \code{set_prior("<prior>", class = "<parameter>")} it the right way to go. \cr
#' 
#'   Often, it may not be immediately clear, which parameters are present in the model.
#'   To get a full list of parameters for which priors can be specified (depending on the model) 
#'   use method \code{\link[brms:parnames.formula]{parnames.formula}}.
#'
#' @references 
#' Gelman A (2006). Prior distributions for variance parameters in hierarchical models."
#'    Bayesian analysis, 1(3), 515 -- 534.
#' 
#' @examples
#' \dontrun{
#' ## check which parameters can have priors
#' parnames(rating ~ period + carry + (1|subject),
#'          data = inhaler, family = "sratio", 
#'          partial = ~ treat, threshold = "equidistant")
#'          
#' ## define some priors          
#' prior <- list(set_prior("normal(0,10)", class = "b"),
#'               set_prior("normal(1,2)", class = "b", coef = "treat"),
#'               set_prior("cauchy(0,2)", class = "sd", 
#'                         group = "subject", coef = "Intercept"),
#'               set_prior("uniform(-5,5)", class = "delta"))
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
  if (nchar(group) && !class %in% c("sd", "cor", "L"))
    stop(paste("Group attribute not meaningful for class", class))
  valid_classes <- c("b", "bp", "sd", "cor", "L", "ar", "ma", "sigma", 
                     "rescor", "Lrescor", "nu", "shape", "delta")
  if (!class %in% valid_classes)
    stop(paste(class, "is not a valid paramter class"))
  out <- list(prior = prior, class = class, coef = coef, group = group)
  class(out) <- c("brmsprior", "list")
  out
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
  if (is.null(prior)) {
    return(prior)  # nothing to check
  } else if (is(prior, "brmsprior")) {
    # a single prior may be specified without list(.)
    prior <- as.data.frame(prior, stringsAsFactors = FALSE)
  } else if (!is.null(names(prior))) {
    # deprecated prior specification brms < 0.5.0
    warning(paste("Specifying priors using a named list is deprecated. \n",
                  "See help(set_prior) for further information."))
    prior <- update_deprecated_prior(prior)
  } else {
    # the usual case since brms 0.5.0
    lapply(prior, function(x) if (!is(x, "brmsprior")) 
      stop(paste("Elements of prior must be of class brmsprior. \n",
                 "See help(set_prior) for further information.")))
    prior <- data.frame(matrix(unlist(prior), ncol = 4, byrow = TRUE),
                        stringsAsFactors = FALSE)
    names(prior) <- c("prior", "class", "coef", "group")
  }
  dupli <- duplicated(prior[, 2:4])
  if (any(dupli)) {
    stop("Duplicated prior specifications are not allowed. \n")
  }
  
  # check if parameters in prior are valid
  ee <- extract_effects(formula, family = family)  
  possible_prior <- parnames(formula, data = data, family = family, 
                             autocor = autocor, partial = partial, 
                             threshold = threshold, internal = TRUE)
  # index of valid rows in prior
  valid <- which(duplicated(rbind(possible_prior, prior[, 2:4]))) - 
           nrow(possible_prior)
  invalid <- which(!1:nrow(prior) %in% valid)
  if (length(invalid)) {
    message(paste("Prior element", paste(invalid, collapse = ", "),
                  "is invalid and will be ignored."))
  }
  
  # expand lkj correlation prior to full name
  prior$prior <- unlist(lapply(prior$prior, function(p) 
    sub("^lkj\\(", "lkj_corr_cholesky(", p)))
  # rename parameter clases
  prior$class <- rename(prior$class, symbols = c("^cor$", "^rescor$"), 
                        subs = c("L", "Lrescor"), fixed = FALSE)
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
      prior <- do.call(rbind, c(prior, new_rows))  # add new rows
    }
  }
  # get partial priors out of fixef priors
  if (is.formula(partial)) {
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
  rownames(prior) <- 1:nrow(prior)
  prior
}

update_deprecated_prior <- function(prior) {
  # update prior specifications from before 0.5.0
  #
  # Args:
  #   prior: A named list
  #
  # Returns:
  #   a data.frame compatible with check_prior of brms >= 0.5.0
  if (is.null(names(prior)) | length(prior) == 0) {
    stop("Only named lists can be updated")
  }
  prior_names <- names(prior)
  class <- regmatches(prior_names, regexpr("^[^_]+", prior_names))
  group_coef <- substr(prior_names, nchar(class) + 2, nchar(prior_names))
  group <- ifelse(prior_names == "sd", "", ifelse(class == "sd",
                                                  regmatches(group_coef, regexpr("^[^_]+", group_coef)), 
                                                  ifelse(class == "cor", group_coef, "")))
  coef <- substr(group_coef, nchar(group) + ifelse(nchar(group), 2, 1), 
                 nchar(group_coef))
  prior_frame(prior = unlist(prior, use.names = FALSE),
              class = class, coef = coef, group = group)
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
#' @return A data.frame with columns \code{class}, \code{coef}, and \code{group}
#'   and several rows, each providing information on a paramter (or parameter class) on which
#'   priors can be specified. See also \code{\link[brms:set_prior]{set_prior}}.
#' 
#' @examples 
#' parnames(rating ~ treat + period + carry + (1+carry|subject), 
#'          data = inhaler, family = "student")
#'           
#' parnames(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|visit),
#'          data = epilepsy, family = "poisson")          
#' 
#' @export
parnames.formula <- function(x, data = NULL, family = "gaussian", addition = NULL, 
                             autocor = NULL, partial = NULL, 
                             threshold = c("flexible", "equidistant"), 
                             internal = FALSE, ...) {
  
  if (is.null(autocor)) 
    autocor <- cor_arma()
  if (!is(autocor, "cor_brms")) 
    stop("cor must be of class cor_brms")
  threshold <- match.arg(threshold)
  family <- check_family(family[1])
  x <- update_formula(x, addition = addition)
  ee <- extract_effects(x, partial, family = family)
  data <- update_data(data, family = family, effects = ee)
  
  # initialize output
  out <- prior_frame(prior = character(0), class = character(0), 
                     coef = character(0), group = character(0))
  # fixed and partial effects
  fixef <- colnames(get_model_matrix(ee$fixed, data = data))
  if (length(fixef)) {
    out <- rbind(out, prior_frame(class = "b", coef = c("", fixef)))
  }
  if (is.formula(partial)) {
    paref <- colnames(get_model_matrix(partial, data = data, rm_intercept = TRUE))
    out <- rbind(out, prior_frame(class = "b", coef = paref))
  }
  # random effects
  if (length(ee$group)) {
    out <- rbind(out, prior_frame(class = "sd"))  # global sd class
    gs <- unlist(ee$group)
    for (i in 1:length(gs)) {
      ranef <- colnames(get_model_matrix(ee$random[[i]], data = data))
      # include random effects standard deviations
      out <- rbind(out, prior_frame(class = "sd", coef = c("", ranef), group = gs[i]))
      # detect duplicated random effects
      J <- with(out, class == "sd" & group == gs[i] & nchar(coef))
      dupli <- duplicated(out[J, ])
      if (any(dupli)) {
        stop(paste("Duplicated random effects detected for group", gs[i]))
      }
      # include correlation parameters
      if (ee$cor[[i]] && length(ranef) > 1) {
        out <- rbind(out, prior_frame(class = "cor"))  # global cor class  
        out <- rbind(out, prior_frame(class = "cor", group = gs[i]))
        if (internal) 
          out <- rbind(out, prior_frame(class = "L", group = gs[i]))
      }
    }
  }
  # handle additional parameters
  is_ordinal <- family %in% c("cumulative", "sratio", "cratio", "acat") 
  if (is(autocor, "cor_arma") && autocor$p) 
    out <- rbind(out, prior_frame(class = "ar"))
  if (is(autocor, "cor_arma") && autocor$q) 
    out <- rbind(out, prior_frame(class = "ma"))
  if (family %in% c("gaussian", "student", "cauchy") && !is.formula(ee$se))
    out <- rbind(out, prior_frame(class = "sigma", coef = c("", ee$response)))
  if (family == "gaussian" && length(ee$response) > 1)
    out <- rbind(out, prior_frame(class = c("rescor", if (internal) "Lrescor")))
  if (family == "student") 
    out <- rbind(out, prior_frame(class = "nu"))
  if (family %in% c("gamma", "weibull", "negbinomial")) 
    out <- rbind(out, prior_frame(class = "shape"))
  if (is_ordinal && threshold == "equidistant")
    out <- rbind(out, prior_frame(class = "delta"))
  unique(out[,2:4])  # do not return prior column as it is empty anyway
}