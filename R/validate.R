extract_effects <- function(formula, ..., family = NA, nonlinear = NULL, 
                            check_response = TRUE, resp_rhs_all = TRUE) {
  # Parse the model formula and related arguments
  # 
  # Args:
  #   formula: An object of class "formula" using mostly the syntax 
  #            of the \code{lme4} package
  #   ...: Additional objects of class "formula"
  #   family: the model family
  #   nonlinear: a list of formulas specifying non-linear effects
  #   check_response: check if the response part is non-empty?
  #   resp_rhs_all: include response variables on the RHS of $all? 
  # 
  # Returns: 
  #   A named list whose elements depend on the formula input 
  tformula <- formula2string(formula) 
  tfixed <- gsub("\\|+[^~]*~", "~", tformula)
  if (length(nonlinear)) {
    if (grepl("|", tfixed, fixed = TRUE)) {
      stop(paste("Random effects in non-linear models should be specified", 
                 "in the 'nonlinear' argument."), call. = FALSE)
    }
    if (is.ordinal(family) || is.categorical(family)) {
      stop("Non-linear effects are not yet allowed for this family.", 
           call. = FALSE)
    }
    re_terms <- cse_terms <- NULL
  } else {
    # terms() doesn't like non-linear formulas
    term_labels <- rename(attr(terms(formula), "term.labels"), " ", "")
    re_terms <- term_labels[grepl("\\|", term_labels)]
    if (length(re_terms)) {
      re_terms <- paste0("(", re_terms, ")")
      tfixed <- rename(tfixed, c(paste0("+", re_terms), re_terms), "")
    } 
    # category specific effects in ordinal models
    cse_terms <- term_labels[grepl("^cse\\(", term_labels)]
    if (length(cse_terms)) {
      if (!is.na(family[[1]]) && !allows_cse(family)) {
        stop(paste("Category specific effects are only meaningful for", 
                   "families 'sratio', 'cratio', and 'acat'."), 
             call. = FALSE)
      }
      tfixed <- rename(tfixed, c(paste0("+", cse_terms), cse_terms), "")
      cse_terms <- substr(cse_terms, 5, nchar(cse_terms) - 1)
      cse_terms <- formula(paste("~", paste(cse_terms, collapse = "+")))
      if (!length(all.vars(cse_terms))) {
        stop("invalid input to function 'cse'", call. = FALSE)
      }
    }
    if (substr(tfixed, nchar(tfixed), nchar(tfixed)) == "~") {
      tfixed <- paste0(tfixed, "1")
    }
    if (grepl("|", x = tfixed, fixed = TRUE)) {
      stop("Random effects terms should be enclosed in brackets", call. = FALSE)
    }
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
  if (length(cse_terms)) x$cse <- cse_terms
  x$nonlinear <- nonlinear_effects(nonlinear, model = x$fixed)
  
  # handle addition arguments
  fun <- c("se", "weights", "trials", "cat", "cens", "trunc", "disp")
  add_vars <- list()
  if (!is.na(family[[1]])) {
    add <- get_matches("\\|[^~]*~", tformula)
    if (length(add)) {
      # replace deprecated '|' by '+'
      add <- paste("~", rename(substr(add, 2, nchar(add) - 1), "|", "+"))
      add_terms <- attr(terms(formula(add)), "term.labels")
      for (f in fun) {
        matches <- grep(paste0("^", f, "\\(.+\\)$"), add_terms)
        if (length(matches) == 1) {
          x[[f]] <- add_terms[matches]
          add_terms <- add_terms[-matches]
          if (!is.na(x[[f]]) && (add_families(f)[1] == "all" ||
              family$family %in% add_families(f))) {
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
        } else if (length(matches) > 1) {
          stop("Addition arguments may be only defined once.", call. = FALSE)
        } 
      }
      if (length(add_terms))
        stop(paste("Invalid addition part of formula.", 
                   "Please see the 'Details' section of help(brm)"),
             call. = FALSE)
      if (is.formula(x$se) && is.formula(x$disp)) {
        stop(paste("Addition arguments 'se' and 'disp' cannot be used", 
                   "at the same time."), call. = FALSE)
      }
    }
  }
  
  covars <- setdiff(all.vars(rhs(x$fixed)), names(x$nonlinear))
  x$covars <- formula(paste("~", paste(c("1", covars), collapse = "+")))
  # make a formula containing all required variables (element 'all')
  formula_list <- c(if (resp_rhs_all) all.vars(lhs(x$fixed)), add_vars, 
                    x$covars, x$cse, if (!length(x$nonlinear)) rhs(x$fixed), 
                    x$random$form, group_formula, get_offset(x$fixed), 
                    lapply(x$nonlinear, function(nl) nl$all), ...)
  new_formula <- collapse(ulapply(formula_list, plus_rhs))
  x$all <- paste0("update(", tfixed, ", ~ ", new_formula, ")")
  x$all <- eval(parse(text = x$all))
  environment(x$all) <- globalenv()
  
  # extract response variables
  if (check_response) {
    x$respform <- lhs(x$all)
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
      # don't use update on a formula that is possibly non-linear
      x$fixed[[2]] <- quote(response)
      x$all <- update(x$all, response ~ .)
    }  
  }
  x
} 

extract_time <- function(formula) {
  # extract time and grouping variables for autocorrelation structures
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
  formula <- gsub(" ", "", Reduce(paste, deparse(formula))) 
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

nonlinear_effects <- function(x, model = ~ 1) {
  # prepare nonlinear formulas
  # Args:
  #   x: a list for formulas specifying linear predictors for 
  #      non-linear parameters
  #   model: formula of the non-linear model
  # Returns:
  #   A list of objects each returned by extract_effects
  if (length(x)) {
    nleffects <- vector("list", length = length(x))
    for (i in seq_along(x)) {
      if (!is(x[[i]], "formula")) {
        stop("Argument 'nonlinear' must be a list of formulas.")
      }
      if (length(x[[i]]) != 3) {
        stop("Non-linear formulas must be two-sided.")
      }
      nlresp <- all.vars(x[[i]][[2]])
      if (length(nlresp) != 1) {
        stop("LHS of non-linear formula must contain exactly one variable.")
      }
      if (any(ulapply(c(".", "_"), grepl, x = nlresp, fixed = TRUE))) {
        stop("Non-linear parameters should not contain dots or underscores.")
      }
      x[[i]] <- rhs(x[[i]])
      nleffects[[i]] <- extract_effects(x[[i]], check_response = FALSE)
      nleffects[[i]]$nonlinear <- NULL
      names(nleffects)[[i]] <- nlresp
    }
    model_vars <- all.vars(rhs(model))
    missing_pars <- setdiff(names(nleffects), model_vars)
    if (length(missing_pars)) {
      stop(paste("Some non-linear parameters are missing in formula:", 
                 paste(missing_pars, collapse = ", ")), call. = FALSE)
    }
  } else {
    nleffects <- NULL 
  }
  nleffects
}

nonlinear2list <- function(x) {
  # convert a single formula into a list of formulas
  # one for each non-linear parameter
  if (is(x, "formula")) {
    if (length(x) != 3) {
      stop("Non-linear formulas must be two-sided.")
    }
    nlpars <- all.vars(lhs(x))
    x <- lapply(nlpars, function(nlp) update(x, paste(nlp, " ~ .")))
  } else if (!(is.list(x) || is.null(x))) {
    stop("invalid 'nonlinear' argument", call. = FALSE)
  }
  x
}

update_formula <- function(formula, data = NULL, partial = NULL, 
                           nonlinear = NULL) {
  # incorporate addition arguments and category specific effects into formula 
  # 
  # Args:
  #   formula: a model formula 
  #   data: a data.frame or NULL 
  #   partial: a one sided formula containing category specific effects
  #
  # Returns:
  #   an updated formula containing the addition and category specific effects
  fnew <- ". ~ ."
  if (is(partial, "formula")) {
    warning(paste("Argument 'partial' is deprecated. Please use the 'cse'", 
                  "function inside the model formula instead."), 
            call. = FALSE)
    partial <- formula2string(partial, rm = 1)
    fnew <- paste(fnew, "+ cse(", partial, ")")
  } else if (!is.null(partial)) {
    stop("argument 'partial' must be of class formula")
  }
  # to allow the '.' symbol in formula
  try_terms <- try(terms(formula, data = data), silent = TRUE)
  if (!is(try_terms, "try-error")) {
    formula <- formula(try_terms)
  }
  if (fnew != ". ~ .") {
    formula <- update.formula(formula, formula(fnew))
  }
  if (!isTRUE(attr(formula, "nonlinear"))) {
    attr(formula, "nonlinear") <- length(nonlinear) > 0
  }
  formula
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
    new_ranef <- gather_ranef(ee, data = data)
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
  if (isTRUE(attr(formula, "nonlinear"))) {
    # non-linear formulae may cause errors when passed to terms
    # and do not contain random effects anyway
    return(formula)
  }
  if (suppressWarnings(anyNA(re_formula))) {
    re_formula <- ~ 1
  }
  if (is.formula(re_formula)) {
    new_formula <- formula2string(formula)
    old_term_labels <- rename(attr(terms(formula), "term.labels"), " ", "")
    old_re_terms <- old_term_labels[grepl("\\|", old_term_labels)]
    if (length(old_re_terms)) {
      old_re_terms <- paste0("(", old_re_terms, ")")
      # make sure that + before random terms are also removed
      old_re_terms <- c(paste0("+", old_re_terms), old_re_terms)
      new_formula <- rename(new_formula, old_re_terms, "")
      if (grepl("~$", new_formula)) {
        # lhs only formulas are not allowed
        new_formula <- paste(new_formula, "1")
      }
    }
    new_term_labels <- rename(attr(terms(re_formula), "term.labels"),  " ", "")
    new_re_terms <- new_term_labels[grepl("\\|", new_term_labels)]
    if (length(new_re_terms)) {
      new_re_terms <- paste0("(", new_re_terms, ")")
      new_formula <- paste(c(new_formula, new_re_terms), collapse = "+")
    }
    new_formula <- formula(new_formula)
  } else if (is.null(re_formula)) {
    new_formula <- formula
  } else {
    stop("invalid re_formula argument", call. = FALSE)
  } 
  new_formula
}

plus_rhs <- function(x) {
  # take the right hand side of a formula and add a +
  if (is.formula(x)) x <- Reduce(paste, deparse(x[[2]]))
  if (!is.null(x) && nchar(x)) paste("+", paste(x, collapse = "+"))
  else "+ 1"
}

get_fixed <- function(effects) {
  # get fixed effects formulas in a list
  # Args:
  #   effects: object returned by extract_effects
  out <- list(effects$fixed)
  if (!is.null(effects$nonlinear)) {
    out <- c(out, lapply(effects$nonlinear, function(par) par$fixed))
    attr(out, "nonlinear") <- TRUE
  }
  out
}

get_random <- function(effects) {
  # get random effects information in a data.frame
  # Args:
  #   effects: object returned by extract_effects
  if (!is.null(effects$nonlinear)) {
    out <- do.call(rbind, lapply(effects$nonlinear, function(par) par$random))
    # R does not allow duplicated rownames and adds "." causing Stan to fail
    out$nlpar <- get_matches("^[^\\.]+", rownames(out))
    attr(out, "nonlinear") <- TRUE
  } else {
    out <- effects$random
  }
  out
}

get_re_index <- function(i, random) {
  # get random effects index for the ith row of random
  # Args:
  #   i: an index
  #   random: data.frame returned by get_random
  rn <- random$nlpar
  if (isTRUE(attr(random, "nonlinear"))) {
    # each non-linear parameter may has its own random effects
    i <- which(which(rn == rn[i]) == i)
  }
  i
}

get_offset <- function(x) {
  # extract offset terms from a formula
  x <- try(terms(as.formula(x)), silent = TRUE)
  if (is(x, "try-error")) {
    # terms doesn't like non-linear formulas
    offset <- NULL
  } else {
    offset_pos <- attr(x, "offset")
    if (!is.null(offset_pos)) {
      vars <- attr(x, "variables")
      offset <- ulapply(offset_pos, function(i) deparse(vars[[i+1]]))
      offset <- paste(offset, collapse = "+")
    } else {
      offset <- NULL
    }
  }
  offset
}

get_var_combs <- function(x) {
  # get all variable combinations occuring in elements of x
  # Args:
  #   x: a character vector
  stopifnot(is.character(x))
  unique(lapply(x, function(y) all.vars(parse(text = y))))
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
  try_terms <- try(terms(formula), silent = TRUE)
  if (is(try_terms, "try-error")) {
    out <- FALSE
  } else {
    has_intercept <- attr(try_terms, "intercept")
    out <- !has_intercept && "intercept" %in% all.vars(rhs(formula))
  }
  out
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

gather_ranef <- function(effects, data = NULL, ...) {
  # gathers helpful information on the random effects
  # 
  # Args:
  #   effects: output of extract_effects
  #   data: data passed to brm after updating
  #   ...: Further arguments passed to get_model_matrix
  #
  # Returns: 
  #   A named list with one element per grouping factor
  random <- get_random(effects)
  Z <- lapply(random$form, get_model_matrix, data = data, ...)
  ranef <- setNames(lapply(Z, colnames), random$group)
  for (i in seq_along(ranef)) {
    attr(ranef[[i]], "levels") <- 
      levels(as.factor(get(random$group[[i]], data)))
    attr(ranef[[i]], "group") <- names(ranef)[i]
    attr(ranef[[i]], "cor") <- random$cor[[i]]
    if (length(effects$nonlinear))
      attr(ranef[[i]], "nlpar") <- random$nlpar[i]
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

cse <- function(...) {
  stop("inappropriate use of function 'cse'", call. = FALSE)
}

add_families <- function(x) {
  # return names of valid families for addition argument x
  switch(x, weights = "all",
         se = c("gaussian", "student", "cauchy"),
         trials = c("binomial", "zero_inflated_binomial"),
         cat = c("categorical", "cumulative", 
               "cratio", "sratio", "acat"), 
         cens = c("gaussian", "student", "cauchy", 
                "inverse.gaussian", "binomial",
                "poisson", "geometric", "negbinomial", 
                "exponential", "weibull", "gamma"),
         trunc = c("gaussian", "student", "cauchy", "binomial",
                 "poisson", "geometric", "negbinomial", 
                 "exponential", "weibull", "gamma"),
         disp = c("gaussian", "student", "cauchy", "gamma",
                  "weibull", "negbinomial"),
         stop(paste("addition argument", x, "is not supported")))
}

is.formula <- function(x, or = TRUE) {
  # checks if x is formula (or list of formulas)
  #
  # Returns:
  #   x: a formula or a list of formulas
  #   or: logical; indicates if any element must be a formula (or = TRUE) 
  #       or if all elements must be formulas
  if (!is.list(x)) x <- list(x)
  out <- sapply(x, function(y) is(y, "formula"))
  if (or) {
    out <- any(out)
  } else out <- all(out)
  out
}

formula2string <- function(formula, rm = c(0, 0)) {
  # converts formula to string
  #
  # Args:
  #   formula: a model formula
  #   rm: a vector of to elements indicating how many characters 
  #       should be removed at the beginning
  #       and end of the string respectively
  #
  # Returns:
  #    the formula as string 
  if (!is.formula(formula))
    stop(paste(deparse(substitute(formula)), "must be of class formula"))
  if (is.na(rm[2])) rm[2] <- 0
  x <- gsub(" ","", Reduce(paste, deparse(formula)))
  x <- substr(x, 1 + rm[1], nchar(x) - rm[2])
  x
} 

get_boundaries <- function(trunc) {
  # extract truncation boundaries out of a formula
  # that is known to contain the .trunc function
  # Returns:
  #   a list containing two numbers named lb and ub
  if (is.formula(trunc)) {
    .addition(trunc)
  } else {
    .trunc()
  }
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

exclude_pars <- function(ranef = list(), save_ranef = TRUE) {
  # list irrelevant parameters NOT to be saved by Stan
  # 
  # Args:
  #   ranef: output of gather_ranef
  #   save_ranef: should random effects of each level be saved?
  #
  # Returns:
  #   a vector of parameters to be excluded
  out <- c("eta", "etap", "eta_2PL", "Eta", "temp_Intercept1", 
           "temp_Intercept",  "Lrescor", "Rescor", "Sigma", "LSigma",
           "disp_sigma", "e", "E", "res_cov_matrix", 
           "lp_pre", "hs_local", "hs_global")
  if (length(ranef)) {
    rm_re_pars <- c("z", "L", "Cor", if (!save_ranef) "r")
    if (!is.null(attr(ranef[[1]], "nlpar"))) {
      nlpars <- ulapply(ranef, function(r) attr(r, "nlpar"))
      out <- c(out, unique(paste0("eta_", nlpars)))
      for (k in seq_along(ranef)) {
        i <- which(which(nlpars == nlpars[k]) == k)
        out <- c(out, paste0(rm_re_pars, "_", nlpars[k], "_", i))
        neff <- length(ranef[[k]])
        if (neff > 1L) {
          out <- c(out, paste0("r_", nlpars[k], "_", i, "_", 1:neff))
        }
      }
    } else {
      for (k in seq_along(ranef)) {
        out <- c(out, paste0(rm_re_pars, "_", k))
        neff <- length(ranef[[k]])
        if (neff > 1L) {
          out <- c(out, paste0("r_", k, "_", 1:neff))
        }
      }
    }
  }
  out
}
