extract_effects <- function(formula, ..., family = NA, check_response = TRUE) {
  # Extract fixed and random effects from a formula
  # 
  # Args:
  #   formula: An object of class "formula" using mostly the syntax of the \code{lme4} package
  #   ...: Additional objects of class "formula"
  #   family: the model family
  #   check_response: check if the response part is non-empty?
  # 
  # Returns: 
  #   A named list of the following elements: 
  #   fixed: An object of class "formula" that contains the fixed effects 
  #          including the dependent variable. 
  #   random: A list of formulas containing the random effects per grouping variable. 
  #   group: A vector of names of the grouping variables. 
  #   weights, se, cens, trials, cat: information on possible addition arguments
  #   all: A formula that contains every variable mentioned in formula and ...
  if (class(family) == "family") {
    family <- family$family
  }
  formula <- formula2string(formula)  
  fixed <- gsub(paste0("\\([^(\\||~)]*\\|[^\\)]*\\)\\+",
                       "|\\+\\([^(\\||~)]*\\|[^\\)]*\\)",
                       "|\\([^(\\||~)]*\\|[^\\)]*\\)"), 
                "", formula)
  fixed <- gsub("\\|+[^~]*~", "~", fixed)
  if (substr(fixed, nchar(fixed), nchar(fixed)) == "~") 
    fixed <- paste0(fixed, "1")
  fixed <- formula(fixed)
  if (family %in% c("cumulative", "sratio", "cratio", "acat"))
    fixed <- update.formula(fixed, . ~ . + 1)
  if (check_response && length(fixed) < 3) 
    stop("invalid formula: response variable is missing")
  
  # extract random effects part
  rg <- gregexpr("\\([^\\|\\)]*\\|[^\\)]*\\)", formula)
  rg <- unlist(regmatches(formula, rg))
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
  fun <- c("se", "weights", "trials", "cat", "cens", "trunc")
  add_vars <- list()
  if (!is.na(family)) {
    add <- unlist(regmatches(formula, gregexpr("\\|[^~]*~", formula)))[1]
    add <- substr(add, 2, nchar(add)-1)
    families <- list(se = c("gaussian", "student", "cauchy"),
                     weights = "all",
                     trials = c("binomial", "binomial_2pl"),
                     cat = c("categorical", "cumulative", "cratio", "sratio", "acat"), 
                     cens = c("gaussian", "student", "cauchy", "inverse.gaussian", 
                              "binomial", "binomial_2pl",
                              "poisson", "geometric", "negbinomial", 
                              "exponential", "weibull", "gamma"),
                     trunc = c("gaussian", "student", "cauchy",
                               "binomial", "binomial_2pl",
                               "poisson", "geometric", "negbinomial", 
                               "exponential", "weibull", "gamma"))
    for (f in fun) {
      x[[f]] <- unlist(regmatches(add, gregexpr(paste0(f,"\\([^\\|]*\\)"), add)))[1]
      add <- gsub(paste0(f,"\\([^~|\\|]*\\)\\|*"), "", add)
      if (is.na(x[[f]])) {
        x[[f]] <- NULL
      } else if (family %in% families[[f]] || families[[f]][1] == "all") {
        args <- substr(x[[f]], nchar(f) + 2, nchar(x[[f]]) - 1)
        try_numeric <- suppressWarnings(as.numeric(args))
        if (f %in% c("trials", "cat") && !is.na(try_numeric)) {
          x[[f]] <- try_numeric
        } else {
          x[[f]] <- as.formula(paste0("~ .", x[[f]]))
          if (length(all.vars(x[[f]]))) {
            form <- paste("~", paste(all.vars(x[[f]]), collapse = "+"))
            add_vars[[f]] <- as.formula(form)
          }
        }
      } else {
        stop(paste("Argument", f, "in formula", 
                   "is not supported by family", family))
      }
    }
    if (nchar(gsub("\\|", "", add)) > 0 && !is.na(add))
      stop(paste("Invalid addition part of formula.", 
                 "Please see the 'Details' section of help(brm)"))
  }
  
  # make a formula containing all required variables (element 'all')
  plus_rh <- function(x) {
    # take the right hand side of a formula and add a +
    paste0("+", Reduce(paste, deparse(x[[2]])))
  }
  formula_list <- c(random, group_formula, add_vars, ...)
  new_formula <- ulapply(formula_list, plus_rh)
  new_formula <- paste0("update(",Reduce(paste, deparse(fixed)),
                        ", ~ .", collapse(new_formula), ")")
  x$all <- eval(parse(text = new_formula))
  environment(x$all) <- globalenv()
  
  # extract response variables
  if (check_response) {
    x$response <- gather_response(x$all[[2]])
    if (is.hurdle(family)) {
      x$response <- c(x$response, paste0("hu_", x$response))
    } else if (is.zero_inflated(family)) {
      x$response <- c(x$response, paste0("zi_", x$response))
    } else if (is.2pl(family)) {
      x$response <- c(x$response, paste0("logDisc_", x$response))
    }
    if (length(x$response) > 1) {
      if (!(is.null(x$cens) && is.null(x$se) && is.null(x$trunc))
          && is.linear(family)) {
        stop("multivariate models currently allow only weights as addition arguments")
      }
      first_var <- all.vars(x$all[[2]])[1]
      x$fixed <- eval(parse(text = paste0("update(x$fixed, ", first_var, " ~ .)"))) 
      x$all <- eval(parse(text = paste0("update(x$all, ", first_var, " ~ .)"))) 
    }  
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
  #   a list with elements time, group, and all, where all contains a 
  #   formula with all variables in formula
  if (is.null(formula)) 
    return(NULL)
  formula <- gsub(" ","",Reduce(paste, deparse(formula))) 
  time <- all.vars(as.formula(paste("~", gsub("~|\\|[[:print:]]*", "", formula))))
  if (length(time) > 1) {
    stop("Autocorrelation structures may only contain 1 time variable")
  }
  x <- list(time = ifelse(length(time), time, ""))
  group <- get_group_formula(sub("~[^\\|]*", "", formula))
  x$group <- paste0(all.vars(group), collapse = ":")
  x$all <- formula(paste("~",paste(c("1", time, all.vars(group)), collapse = "+")))
  x
}

update_formula <- function(formula, addition = NULL, partial = NULL) {
  # incorporate addition arguments and category specific effects into formula 
  # 
  # Args:
  #   formula: a model formula 
  #   addition: a list with one sided formulas taken from the addition arguments in brm
  #   partial: a one sided formula containing category specific effects
  #
  # Returns:
  #   an updated formula containing the addition and category specific effects
  var_names <- names(addition)
  addition <- lapply(addition, formula2string, rm = 1)
  fnew <- "."
  if (length(addition)) {
    warning("Argument addition is deprecated. See help(brm) for further details.")
    for (i in 1:length(addition)) {
      fnew <- paste0(fnew, " | ", var_names[i], "(", addition[[i]], ")")
    }
  }
  fnew <- paste(fnew, "~ .")
  if (is.formula(partial)) {
    partial <- formula2string(partial, rm=1)
    fnew <- paste(fnew, "+ partial(", partial, ")")
  }
  if (fnew == ". ~ .") {
    formula
  } else {
    update.formula(formula, formula(fnew))
  }
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

check_re_formula <- function(re_formula, old_ranef, data) {
  # validate the re_formula argument as passed to predict and fitted
  #
  # Args:
  #   re_formula: see predict.brmsfit for documentation
  #   old_ranef: named list containing the RE names 
  #              of each grouping factor in the original model
  #   data: data supplied by the user
  #
  # Returns:
  #   named list containing the RE names of each grouping factor
  #   as defined in re_formula; or NULL if re_formula is NA or ~ 1
  if (is.null(re_formula)) {
    new_ranef <- old_ranef
  } else if (is.formula(re_formula)) {
    if (!is.data.frame(data)) {
      stop("argument re_formula requires models fitted with brms > 0.5.0")
    }
    if (length(re_formula) == 3) {
      stop("re_formula must be one-sided")
    }
    ee <- extract_effects(re_formula, check_response = FALSE)
    if (length(all.vars(ee$fixed))) {
      stop("fixed effects are not allowed in re_formula")
    }
    if (!length(ee$group)) {
      # if no RE terms are present in re_formula
      return(NULL)
    }
    # the true family doesn't matter here
    data <- update_data(data, family = "none", effects = ee)
    new_ranef <- combine_duplicates(gather_ranef(effects = ee, data = data))
    invalid_gf <- setdiff(names(new_ranef), names(old_ranef))
    if (length(invalid_gf)) {
      stop(paste("Invalid grouping factors detected:", 
                 paste(invalid_gf, collapse = ", ")))
    }
    for (gf in names(new_ranef)) {
      invalid_re <- setdiff(new_ranef[[gf]], old_ranef[[gf]])
      if (length(invalid_re)) {
        stop(paste0("Invalid random effects detected for grouping factor ", 
                    gf, ": ", paste(invalid_re, collapse = ", ")))
      } 
    }
  } else if (is.na(re_formula)) {
    new_ranef <- NULL
  } else {
    stop("invalid re_formula argument")
  }
  new_ranef
}

update_re_terms <- function(formula, re_formula = NULL) {
  # remove RE terms in formula and add RE terms of re_formula
  #
  # Args:
  #   formula: model formula to be updated
  #   re_formula: formula containing new RE terms
  #
  # Returns:
  #  a formula with updated RE terms
  if (suppressWarnings(anyNA(re_formula))) {
    re_formula <- ~ 1
  }
  if (is.formula(re_formula)) {
    formula <- formula2string(formula)  
    re_formula <- formula2string(re_formula)
    fixef_formula <- gsub(paste0("\\([^(\\||~)]*\\|[^\\)]*\\)\\+",
                                 "|\\+\\([^(\\||~)]*\\|[^\\)]*\\)",
                                 "|\\([^(\\||~)]*\\|[^\\)]*\\)"), 
                          "", formula)
    new_re_terms <- gregexpr("\\([^\\|\\)]*\\|[^\\)]*\\)", re_formula)
    new_re_terms <- unlist(regmatches(re_formula, new_re_terms))
    new_formula <- paste(c(fixef_formula, new_re_terms), collapse = "+")
    new_formula <- formula(new_formula)   
  } else if (is.null(re_formula)) {
    new_formula <- formula
  } else {
    stop("invalid re_formula argument")
  } 
  new_formula
}

gather_response <- function(expr) {
  # gather response variable names
  # Args:
  #   expr: an expression containing the reponse part
  #         of the model formula
  # Returns:
  #   a vector of names of the response variables (columns)
  all_vars <- all.vars(expr)
  if (length(all_vars) == 0) {
    stop("formula must contain at least one response variable")
  }
  pseudo_data <- as.data.frame(setNames(as.list(rep(1, length(all_vars))), 
                                        all_vars))
  resp_formula <- formula(paste(deparse(expr), " ~ 1"))
  pseudo_resp <- model.response(model.frame(resp_formula, data = pseudo_data))
  if (is.null(dim(pseudo_resp))) {
    # response is a vector
    response <- all_vars[1]
  } else if (length(dim(pseudo_resp)) == 2) {
    # response is a matrix
    response <- colnames(pseudo_resp)
    empty_names <- which(!nchar(response))
    if (length(empty_names)) {
      response[empty_names] <- paste0("response", empty_names)
    }
  }
  response
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

family.character <- function(object, link = NA, ...) {
  # build a family object
  # 
  # Args:
  #   object: A character string defining the family
  #   link: A character string defining the link
  family <- tolower(object)
  # check validity of family
  if (family == "normal")
    family <- "gaussian"
  if (family == "multigaussian") 
    stop("family 'multigaussian' is deprecated. Use family 'gaussian' instead")
  okFamilies <- c("gaussian", "student", "cauchy", 
                  "binomial", "bernoulli", "categorical", 
                  "poisson", "negbinomial", "geometric", 
                  "gamma", "weibull", "exponential", "inverse.gaussian", 
                  "cumulative", "cratio", "sratio", "acat",
                  "hurdle_poisson", "hurdle_negbinomial", "hurdle_gamma",
                  "zero_inflated_poisson", "zero_inflated_negbinomial",
                  "beta", "bernoulli_2pl", "binomial_2pl")
  if (!family %in% okFamilies)
    stop(paste(family, "is not a supported family. Supported families are: \n",
               paste(okFamilies, collapse = ", ")))
  
  # check validity of link
  if (is.linear(family)) {
    okLinks <- c("identity", "log", "inverse")
  } else if (family == "inverse.gaussian") {
    okLinks <- c("1/mu^2", "inverse", "identity", "log")
  } else if (is.count(family)) {
    okLinks <- c("log", "identity", "sqrt")
  } else if (is.binary(family) || is.ordinal(family) || family == "beta") {
    okLinks <- c("logit", "probit", "probit_approx", "cloglog", "cauchit")
  } else if (family == "categorical") {
    okLinks <- c("logit")
  } else if (is.skewed(family)) {
    okLinks <- c("log", "identity", "inverse")
  } else if (is.hurdle(family) || is.zero_inflated(family)) {
    okLinks <- c("log")
  } else if (is.2pl(family)) {
    okLinks <- c("logit")
  } 
  if (is.na(link)) {
    link <- okLinks[1]
  }
  if (!link %in% okLinks)
    stop(paste0(link, " is not a supported link for family ", family, ". ", 
                "Supported links are: \n", paste(okLinks, collapse = ", ")))
  structure(list(family = family, link = link), class = "family")
}

check_family <- function(family) {
  # checks and corrects validity of the model family
  #
  # Args:
  #   family: Either a function, an object of class 'family' of a character string
  if (is.function(family)) {
    family <- family()   
  }
  if (is(family, "family")) {
    family <- family(family$family, link = family$link)
  } else if (is.character(family)) {
    family <- family(family[1], link = family[2])
  } else {
    stop("family argument is invalid")
  }
  family
}

exclude_pars <- function(formula, ranef = TRUE) {
  # list irrelevant parameters NOT to be saved by Stan
  # 
  # Args:
  #   formula: a model formula
  #   ranef: logical; should random effects of each level be saved?
  #
  # Returns:
  #   a vector of parameters to be excluded
  ee <- extract_effects(formula)
  out <- c("eta", "etam", "etap", "eta_2pl", "b_Intercept1", 
           "Lrescor", "Rescor", "Sigma", "LSigma",
           "p", "q", "e", "E", "res_cov_matrix", 
           "lp_pre", "hs_local", "hs_global")
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