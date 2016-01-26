extract_effects <- function(formula, ..., family = NA, 
                            check_response = TRUE) {
  # Extract fixed and random effects from a formula
  # 
  # Args:
  #   formula: An object of class "formula" using mostly the syntax 
  #            of the \code{lme4} package
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
  term_labels <- rename(attr(terms(formula), "term.labels"), " ", "")
  tformula <- formula2string(formula) 
  tfixed <- gsub("\\|+[^~]*~", "~", tformula)
  re_terms <- term_labels[grepl("\\|", term_labels)]
  if (length(re_terms)) {
    re_terms <- paste0("(", re_terms, ")")
    # make sure that + before random terms are also removed
    extended_re_terms <- c(paste0("+", re_terms), re_terms)
    tfixed <- rename(tfixed, extended_re_terms, "")
  } 
  if (substr(tfixed, nchar(tfixed), nchar(tfixed)) == "~") {
    tfixed <- paste0(tfixed, "1")
  }
  if (grepl("|", x = tfixed, fixed = TRUE)) {
    stop("Random effects terms should be enclosed in brackets", call. = FALSE)
  }
  fixed <- formula(tfixed)
  rsv_intercept <- has_rsv_intercept(fixed)
  if (!is.na(family[[1]])) 
    family <- check_family(family)
  if (is.ordinal(family) || rsv_intercept)
    fixed <- update.formula(fixed, . ~ . + 1)
  if (rsv_intercept)
    attr(fixed, "rsv_intercept") <- TRUE
  if (check_response && length(fixed) < 3) 
    stop("Invalid formula: response variable is missing", call. = FALSE)
  
  # extract random effects parts
  form <- lapply(get_matches("\\([^\\|]*", re_terms), function(r) 
                 formula(paste0("~ ", substr(r, 2, nchar(r)))))
  group <- get_matches("\\|[^\\)]*", re_terms)
  group_formula <- lapply(group, get_group_formula)
  group <- ulapply(group_formula, function(g) 
                   paste0(all.vars(g), collapse = ":"))
  cor <- ulapply(get_matches("\\|[^\\)]*", re_terms), 
                 function(g) substr(g, 1, 2) != "||")
  random <- data.frame(group = group, cor = cor, 
                       stringsAsFactors = FALSE)
  # ensure that all REs of the same gf are next to each other
  if (nrow(random)) {
    random$form <- form
    random <- random[order(random$group), ]
  }
  x <- nlist(fixed, random)
  
  # handle addition arguments
  fun <- c("se", "weights", "trials", "cat", "cens", "trunc")
  add_vars <- list()
  if (!is.na(family[[1]])) {
    add <- get_matches("\\|[^~]*~", tformula)[1]
    add <- substr(add, 2, nchar(add)-1)
    families <- list(se = c("gaussian", "student", "cauchy"),
                     weights = "all",
                     trials = c("binomial", "zero_inflated_binomial"),
                     cat = c("categorical", "cumulative", 
                             "cratio", "sratio", "acat"), 
                     cens = c("gaussian", "student", "cauchy", 
                              "inverse.gaussian", "binomial",
                              "poisson", "geometric", "negbinomial", 
                              "exponential", "weibull", "gamma"),
                     trunc = c("gaussian", "student", "cauchy", "binomial",
                               "poisson", "geometric", "negbinomial", 
                               "exponential", "weibull", "gamma"))
    for (f in fun) {
      x[[f]] <- get_matches(paste0(f, "\\([^\\|]*\\)"), add)[1]
      add <- gsub(paste0(f,"\\([^~|\\|]*\\)\\|*"), "", add)
      add_present <- 
      if (is.na(x[[f]])) {
        x[[f]] <- NULL
      } else if (family$family %in% families[[f]] || 
                 families[[f]][1] == "all") {
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
        stop(paste("Argument", f, "in formula is not supported", 
                   "by family", family$family), call. = FALSE)
      }
    }
    if (nchar(gsub("\\|", "", add)) > 0 && !is.na(add))
      stop(paste("Invalid addition part of formula.", 
                 "Please see the 'Details' section of help(brm)"),
           call. = FALSE)
  }
  
  # make a formula containing all required variables (element 'all')
  plus_rh <- function(x) {
    # take the right hand side of a formula and add a +
    if (is.formula(x)) {
      paste0("+", Reduce(paste, deparse(x[[2]])))
    } else ""
  }
  formula_list <- c(random$form, group_formula, add_vars, ...)
  new_formula <- collapse(ulapply(formula_list, plus_rh))
  x$all <- paste0("update(", tfixed, ", ~ .", new_formula, ")")
  x$all <- eval(parse(text = x$all))
  environment(x$all) <- globalenv()
  
  # extract response variables
  if (check_response) {
    x$respform <- update(x$all, . ~ 1)
    x$response <- gather_response(x$respform)
    if (is.hurdle(family)) {
      x$response <- c(x$response, paste0("hu_", x$response))
    } else if (is.zero_inflated(family)) {
      x$response <- c(x$response, paste0("zi_", x$response))
    } else if (is.2PL(family)) {
      x$response <- c(x$response, paste0("logDisc_", x$response))
    }
    if (length(x$response) > 1) {
      if (!(is.null(x$cens) && is.null(x$se) && is.null(x$trunc))
          && is.linear(family)) {
        stop(paste("Multivariate models currently allow", 
                   "only weights as addition arguments"), 
             call. = FALSE)
      }
      x$fixed <- update(x$fixed, response ~ .)
      x$all <- update(x$all, response ~ .)
    }  
  }
  x
} 

extract_time <- function(formula) {
  # extract time and grouping variabels for correlation structure
  # 
  # Args:
  #   formula: a one sided formula of the form ~ time|group 
  #            typically taken from a cor_brms object
  # 
  # Returns: 
  #   a list with elements time, group, and all, where all contains a 
  #   formula with all variables in formula
  if (is.null(formula)) 
    return(NULL)
  formula <- gsub(" ","",Reduce(paste, deparse(formula))) 
  time <- all.vars(as.formula(paste("~", gsub("~|\\|[[:print:]]*", "", formula))))
  if (length(time) > 1) {
    stop("Autocorrelation structures may only contain 1 time variable", 
         call. = FALSE)
  }
  x <- list(time = ifelse(length(time), time, ""))
  group <- get_group_formula(sub("~[^\\|]*", "", formula))
  x$group <- paste0(all.vars(group), collapse = ":")
  x$all <- formula(paste("~",paste(c("1", time, all.vars(group)), collapse = "+")))
  x
}

update_formula <- function(formula, data = NULL, addition = NULL, 
                           partial = NULL) {
  # incorporate addition arguments and category specific effects into formula 
  # 
  # Args:
  #   formula: a model formula 
  #   data: a data.frame or NULL 
  #   addition: a list with one sided formulas taken from the addition arguments in brm
  #   partial: a one sided formula containing category specific effects
  #
  # Returns:
  #   an updated formula containing the addition and category specific effects
  var_names <- names(addition)
  addition <- lapply(addition, formula2string, rm = 1)
  fnew <- "."
  if (length(addition)) {
    warning("Argument addition is deprecated. See help(brm) for further details.",
            call. = FALSE)
    for (i in 1:length(addition)) {
      fnew <- paste0(fnew, " | ", var_names[i], "(", addition[[i]], ")")
    }
  }
  fnew <- paste(fnew, "~ .")
  if (is.formula(partial)) {
    partial <- formula2string(partial, rm = 1)
    fnew <- paste(fnew, "+ partial(", partial, ")")
  }
  # to allow the '.' symbol in formula
  formula <- formula(terms(formula, data = data))
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
               "may contain only variable names combined by the symbol ':'"),
         call. = FALSE)
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
      stop("argument re_formula requires models fitted with brms > 0.5.0",
           call. = FALSE)
    }
    if (length(re_formula) == 3) {
      stop("re_formula must be one-sided", call. = FALSE)
    }
    ee <- extract_effects(re_formula, check_response = FALSE)
    if (length(all.vars(ee$fixed))) {
      stop("fixed effects are not allowed in re_formula", call. = FALSE)
    }
    if (!nrow(ee$random)) {
      # if no RE terms are present in re_formula
      return(NULL)
    }
    # the true family doesn't matter here
    data <- update_data(data, family = NA, effects = ee)
    new_ranef <- gather_ranef(random = ee$random, data = data)
    new_ranef <- combine_duplicates(new_ranef)
    invalid_gf <- setdiff(names(new_ranef), names(old_ranef))
    if (length(invalid_gf)) {
      stop(paste("Invalid grouping factors detected:", 
                 paste(invalid_gf, collapse = ", ")), call. = FALSE)
    }
    for (gf in names(new_ranef)) {
      invalid_re <- setdiff(new_ranef[[gf]], old_ranef[[gf]])
      if (length(invalid_re)) {
        stop(paste0("Invalid random effects detected for grouping factor ", 
                    gf, ": ", paste(invalid_re, collapse = ", ")),
             call. = FALSE)
      } 
    }
  } else if (is.na(re_formula)) {
    new_ranef <- NULL
  } else {
    stop("invalid re_formula argument", call. = FALSE)
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
    new_re_terms <- get_matches("\\([^\\|\\)]*\\|[^\\)]*\\)", re_formula)
    new_formula <- paste(c(fixef_formula, new_re_terms), collapse = "+")
    new_formula <- formula(new_formula)   
  } else if (is.null(re_formula)) {
    new_formula <- formula
  } else {
    stop("invalid re_formula argument", call. = FALSE)
  } 
  new_formula
}

amend_terms <- function(x, rm_intercept = FALSE, is_forked = FALSE) {
  # amend a terms object (or one that can be coerced to it)
  # to be used in get_model_matrix
  # Args:
  #   x: any R object; if not a formula or terms, NULL is returned
  #   rm_intercept: a flag indicating if the intercept column 
  #                 should be removed from the model.matrix. 
  #                 Primarily useful for ordinal models
  #   is_forked: a flag indicating if the model is forked into
  #              two parts (e.g., a hurdle model)
  # Returns:
  #   a (possibly amended) terms object or NULL
  if (is.formula(x) || is(x, "terms")) {
    x <- terms(x)
  } else {
    return(NULL)
  }
  attr(x, "rm_intercept") <- as.logical(rm_intercept)
  if (is_forked) {
    # ensure that interactions with main and spec won't
    # cause automatic cell mean coding of factors
    term_labels <- attr(x, "term.labels")
    if (any(grepl("(^|:)(main|spec)($|:)", term_labels))) {
      if (any(grepl("(^|:)trait($|:)", term_labels))) {
        stop(paste("formula may not contain variable 'trait'",
                   "when using variables 'main' or 'spec'"),
             call. = FALSE)
      }
      if (attr(x, "intercept")) {
        stop(paste("formula may not contain an intercept",
                   "when using variables 'main' or 'spec'"),
             call. = FALSE)
      }
      attr(x, "intercept") <- 1
      attr(x, "rm_intercept") <- TRUE
    }
  }
  x
}

check_intercept <- function(names) {
  # check if model contains fixed effects intercept
  #
  # Args:
  #   names: The names of the design matrix
  #          to be checked for an intercept
  # Returns:
  #   a list containing the updated effect names
  #   as well as an indicator if the model has an intercept
  if (!is.null(names)) {
    has_intercept <- "Intercept" == names[1]
    if (has_intercept) names <- names[-1]
  } else {
    has_intercept <- FALSE
  }
  nlist(names, has_intercept)
}

has_rsv_intercept <- function(formula) {
  # check if model makes use of the reserved variable 'intercept'
  # Args:
  #   formula: a model formula
  formula <- as.formula(formula)
  has_intercept <- attr(terms(formula), "intercept")
  rhs <- if (length(formula) == 2L) formula[[2]]
         else formula[[3]]
  !has_intercept && "intercept" %in% all.vars(rhs)
}

gather_response <- function(formula) {
  # gather response variable names
  # Args:
  #   formula: a formula containing only the model reponse
  # Returns:
  #   a vector of names of the response variables (columns)
  stopifnot(is.formula(formula))
  all_vars <- all.vars(formula)
  if (length(all_vars) == 0) {
    stop("formula must contain at least one response variable", call. = FALSE)
  }
  mf <- as.data.frame(setNames(as.list(rep(1, length(all_vars))), 
                               all_vars))
  mf <- model.frame(formula, data = mf, na.action = NULL)
  pseudo_resp <- model.response(mf)
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

gather_ranef <- function(random, data = NULL, ...) {
  # gathers helpful information on the random effects
  #
  # Args:
  #   effects: output of extract_effects
  #   data: data passed to brm after updating
  #   ...: Further arguments passed to get_model_matrix
  #
  # Returns: 
  #   A named list with one element per grouping factor
  Z <- lapply(random$form, get_model_matrix, data = data, ...)
  ranef <- setNames(lapply(Z, colnames), random$group)
  for (i in seq_along(ranef)) {
    attr(ranef[[i]], "levels") <- 
      levels(as.factor(get(random$group[[i]], data)))
    attr(ranef[[i]], "group") <- names(ranef)[i]
    attr(ranef[[i]], "cor") <- random$cor[[i]]
  }
  ranef
}

rsv_vars <- function(family, nresp = 1) {
  # returns reservered variables for the family
  if (is.forked(family)) {
    rsv <- c("trait", "main", "spec")
  } else if (is.linear(family) && nresp > 1) {
    rsv <- "trait"
  } else {
    rsv <- NULL
  }
  rsv
}

check_brm_input <- function(x) {
  # misc checks on brm arguments 
  # Args:
  #   x: A named list
  if (x$chains %% x$cluster != 0) {
    stop("chains must be a multiple of cluster", call. = FALSE)
  }
  family <- check_family(x$family) 
  if (family$family %in% c("exponential", "weibull") && 
      x$inits == "random") {
    warning(paste("Families exponential and weibull may not work well",
                  "with default initial values. \n",
                  " It is thus recommended to set inits = '0'"), 
            call. = FALSE)
  }
  if (family$family == "cauchy") {
    warning(paste("Family cauchy is deprecated and will be removed soon",
                  "as it often has convergence issues and not much",
                  "practical application anyway."), call. = FALSE)
  }
  if (family$family == "inverse.gaussian") {
    warning(paste("Inverse gaussian models require carefully chosen", 
                  "prior distributions to ensure convergence of the chains."),
            call. = FALSE)
  }
  if (family$link == "sqrt") {
    warning(paste(family$family, "model with sqrt link may not be", 
                  "uniquely identified"), call. = FALSE)
  }
  invisible(NULL)
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
  out <- c("eta", "etap", "eta_2PL", "Eta", 
           "temp_Intercept1", "temp_Intercept", 
           "Lrescor", "Rescor", "Sigma", "LSigma",
           "p", "q", "e", "E", "res_cov_matrix", 
           "lp_pre", "hs_local", "hs_global")
  for (i in seq_along(ee$random$group)) {
    out <- c(out, paste0("pre_",i), paste0("L_",i), paste0("Cor_",i))
    if (!ranef) out <- c(out, paste0("r_",i))
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
    warning(paste("chain", i, "did not contain samples", 
                  "and was removed from the fitted model"))
    return(NULL)  
  } else {
    return(sflist[[i]])
  }
}