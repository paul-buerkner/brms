extract_effects <- function(formula, ..., family = NA, nonlinear = NULL, 
                            check_response = TRUE, resp_rhs_all = TRUE) {
  # Parse the model formula and related arguments
  # Args:
  #   formula: An object of class 'formula' or 'brmsformula'
  #   ...: Additional objects of class "formula"
  #   family: the model family
  #   nonlinear: a list of formulas specifying non-linear effects
  #   check_response: check if the response part is non-empty
  #   resp_rhs_all: include response variables on the RHS of $all? 
  # Returns: 
  #   A named list whose elements depend on the formula input
  old_mv <- isTRUE(attr(formula, "old_mv"))
  formula <- bf(formula, nonlinear = nonlinear)
  nonlinear <- attr(formula, "nonlinear")
  if (!is.na(family[[1]])) {
    family <- check_family(family)
  }
  tformula <- formula2string(formula) 
  tfixed <- gsub("\\|+[^~]*~", "~", tformula)
  x <- nlist(formula)
  if (length(nonlinear)) {
    if (grepl("|", tfixed, fixed = TRUE)) {
      stop(paste("Random effects in non-linear models should be specified", 
                 "in the 'nonlinear' argument."), call. = FALSE)
    }
    if (is.ordinal(family) || is.categorical(family) || is.forked(family)) {
      stop("Non-linear effects are not yet allowed for this family.", 
           call. = FALSE)
    }
    x$fixed <- formula(tfixed)
    x$nonlinear <- nonlinear_effects(nonlinear, x$fixed, family = family)
    re_terms <- NULL
  } else {
    # terms() doesn't like non-linear formulas
    terms <- terms(formula)
    all_terms <- gsub("[ \t\r\n]+", "", perl = TRUE,
                      x = attr(terms, "term.labels"))
    pos_re_terms <- grep("\\|", all_terms)
    re_terms <- all_terms[pos_re_terms]
    # monotonic effects
    pos_mono_terms <- grep("^mono(|tonic|tonous)\\(", all_terms)
    mono_terms <- all_terms[pos_mono_terms]
    if (length(mono_terms)) {
      mono_terms <- sub("^mono(|tonous)\\(", "monotonic(", mono_terms)
      mono_terms <- substr(mono_terms, 11, nchar(mono_terms) - 1)
      mono_terms <- formula(paste("~", paste(mono_terms, collapse = "+")))
      attr(mono_terms, "rsv_intercept") <- TRUE
      if (!length(all.vars(mono_terms))) {
        stop("invalid input to function 'monotonic'", call. = FALSE)
      }
      if (any(grepl(":", attr(terms(mono_terms), "term.labels")))) {
        stop("Interactions cannot be modeled as monotonic effects.",
             call. = FALSE)
      }
      x$mono <- mono_terms
    }
    # category specific effects in ordinal models
    pos_cse_terms <- grep("^cse\\(", all_terms)
    cse_terms <- all_terms[pos_cse_terms]
    if (length(cse_terms)) {
      if (!is.na(family[[1]]) && !allows_cse(family)) {
        stop(paste("Category specific effects are only meaningful for", 
                   "families 'sratio', 'cratio', and 'acat'."), 
             call. = FALSE)
      }
      cse_terms <- substr(cse_terms, 5, nchar(cse_terms) - 1)
      cse_terms <- formula(paste("~", paste(cse_terms, collapse = "+")))
      attr(cse_terms, "rsv_intercept") <- TRUE
      if (!length(all.vars(cse_terms))) {
        stop("invalid input to function 'cse'", call. = FALSE)
      }
      x$cse <- cse_terms
    }
    # parse spline expression for GAMMs
    pos_spline_terms <- grep("^(s|t2|te|ti)\\(", all_terms)
    spline_terms <- all_terms[pos_spline_terms]
    if (length(spline_terms)) {
      if (is.mv(family) || is.categorical(family)) {
        stop("Splines are not yet implemented for this family.", 
             call. = FALSE)
      }
      if (any(grepl("^(te|ti)\\(", spline_terms))) {
        stop(paste("Tensor product splines 'te' and 'ti' are not yet", 
                   "implemented in brms. Consider using 't2' instead."),
             call. = FALSE)
      }
      x$gam <- formula(paste("~", paste(spline_terms, collapse = "+")))
    }
    rm_terms <- c(pos_re_terms, pos_mono_terms, 
                  pos_cse_terms, pos_spline_terms)
    fe_terms <- all_terms
    if (length(rm_terms)) {
      fe_terms <- fe_terms[-rm_terms]
    }
    int_term <- ifelse(attr(terms, "intercept") == 1, "1", "0")
    fe_terms <- paste(c(int_term, fe_terms, get_offset(formula)), 
                      collapse = "+")
    tfixed <- paste(sub("~.*", "", tfixed), "~", fe_terms)
    x$fixed <- formula(tfixed)
    if (is.ordinal(family)) {
      x$fixed <- update.formula(x$fixed, . ~ . + 1)
    }
    if (has_rsv_intercept(x$fixed)) {
      attr(x$fixed, "rsv_intercept") <- TRUE
    }
    if (is.forked(family)) {
      attr(x$fixed, "forked") <- TRUE
    }
  }
  if (check_response && length(x$fixed) < 3L) { 
    stop("Invalid formula: response variable is missing", call. = FALSE)
  }
  # make sure to store the plain names of all predictors
  covars <- setdiff(all.vars(rhs(x$fixed)), names(x$nonlinear))
  x$covars <- formula(paste("~", paste(c("1", covars), collapse = "+")))
  attr(x$covars, "rsv_intercept") <- TRUE
  # extract random effects parts
  x$random <- extract_random(re_terms)
  
  # evaluate formulas for auxiliary parameters
  auxpars <- sformula(formula, incl_nl = FALSE)
  for (ap in names(auxpars)) {
    x[[ap]] <- extract_effects(rhs(auxpars[[ap]]), family = family,
                               check_response = FALSE)
  }
  
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
        } else if (length(matches) > 1L) {
          stop("Addition arguments may be only defined once.", call. = FALSE)
        } 
      }
      if (length(add_terms)) {
        stop(paste("Invalid addition part of formula.", 
                   "Please see the 'Details' section of help(brm)"),
             call. = FALSE)
      }
      if (is.formula(x$se) && is.formula(x$disp)) {
        stop(paste("Addition arguments 'se' and 'disp' cannot be used", 
                   "at the same time."), call. = FALSE)
      }
    }
  }
  
  # make a formula containing all required variables (element 'all')
  formula_list <- c(
    if (resp_rhs_all) all.vars(lhs(x$fixed)), 
    add_vars, x[c("covars", "cse", "mono")], all.vars(rhs(x$gam)),  
    if (!length(x$nonlinear)) c(rhs(x$fixed), all.vars(rhs(x$fixed))), 
    x$random$form, lapply(x$random$form, all.vars), x$random$group, 
    get_offset(x$fixed), lapply(x$nonlinear, function(e) e$all),
    lapply(x[auxpars()], function(e) e$all), ...)
  new_formula <- collapse(ulapply(formula_list, plus_rhs))
  x$all <- paste0("update(", tfixed, ", ~ ", new_formula, ")")
  x$all <- eval(parse(text = x$all))
  environment(x$all) <- environment(formula)
  
  # extract response variables
  if (check_response) {
    x$respform <- lhs(x$all)
    if (!is.null(attr(formula, "response"))) {
      x$response <- attr(formula, "response")
    } else { 
      x$response <- gather_response(x$respform)
    }
    if (old_mv) {
      # multivariate ('trait') syntax is deprecated as of brms 1.0.0
      if (is.hurdle(family)) {
        x$response <- c(x$response, paste0("hu_", x$response))
      } else if (is.zero_inflated(family)) {
        x$response <- c(x$response, paste0("zi_", x$response))
      } else if (is.2PL(family)) {
        x$response <- c(x$response, paste0("logDisc_", x$response))
      } 
      if (length(x$response) > 1L) {
        if (is.linear(family) && length(rmNULL(x[c("se", "cens", "trunc")]))) {
          stop(paste("Multivariate models currently allow", 
                     "only weights as addition arguments"), 
               call. = FALSE)
        }
        # don't use update on a formula that is possibly non-linear
        x$fixed[[2]] <- quote(response)
        x$all <- update(x$all, response ~ .)
        attr(x$formula, "old_mv") <- TRUE
      }   
    }
  }
  # check validity of auxiliary parameters
  if (!is.na(family[[1]])) {
    # don't do this earlier to allow exlusion of MV models
    inv_auxpars <- setdiff(auxpars(), valid_auxpars(family, effects = x))
    inv_auxpars <- intersect(inv_auxpars, names(x))
    if (length(inv_auxpars)) {
      stop("Prediction of the parameter(s) ", 
           paste0("'", inv_auxpars, "'", collapse = ", "),
           " is not allowed for this model.", call. = FALSE)
    }
  }
  x
} 

extract_time <- function(formula) {
  # extract time and grouping variables for autocorrelation structures
  # Args:
  #   formula: a one sided formula of the form ~ time|group 
  #            typically taken from a cor_brms object
  # Returns: 
  #   a list with elements time, group, and all, where all contains a 
  #   formula with all variables in formula
  if (is.null(formula)) return(NULL)
  if (!is.null(lhs(formula))) {
    stop("autocorrelation formula must be one-sided", call. = FALSE)
  }
  formula <- formula2string(formula)
  time <- as.formula(paste("~", gsub("~|\\|[[:print:]]*", "", formula)))
  time <- all.vars(time)
  if (length(time) > 1L) {
    stop("Autocorrelation structures may only contain 1 time variable", 
         call. = FALSE)
  }
  x <- list(time = ifelse(length(time), time, ""))
  group <- sub("^\\|*", "", sub("~[^\\|]*", "", formula))
  if (illegal_group_expr(group, bs_valid = FALSE)) {
    stop(paste("Illegal grouping term:", group, "\n may contain only", 
               "variable names combined by the symbols ':'"),
         call. = FALSE)
  }
  group <- formula(paste("~", ifelse(nchar(group), group, "1")))
  x$group <- paste0(all.vars(group), collapse = ":")
  x$all <- formula(paste("~", paste(c("1", time, all.vars(group)), 
                                    collapse = "+")))
  x
}

nonlinear_effects <- function(x, model = ~1, family = NA) {
  # prepare nonlinear formulas
  # Args:
  #   x: a list for formulas specifying linear predictors for 
  #      non-linear parameters
  #   model: formula of the non-linear model
  # Returns:
  #   A list of objects each returned by extract_effects
  stopifnot(is.list(x), is(model, "formula"))
  if (length(x)) {
    nleffects <- named_list(names(x))
    for (i in seq_along(x)) {
      nleffects[[i]] <- extract_effects(rhs(x[[i]]), family = family, 
                                        check_response = FALSE)
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

valid_auxpars <- function(family, effects = list(), autocor = cor_arma()) {
  # convenience function to find relevant auxiliary parameters
  x <- c(sigma = has_sigma(family, effects = effects, autocor = autocor),
         shape = has_shape(family), nu = has_nu(family), phi = has_phi(family),
         zi = is.zero_inflated(family, zi_beta = TRUE), 
         hu = is.hurdle(family, zi_beta = FALSE))
  names(x)[x]
}

avoid_auxpars <- function(names, effects) {
  # avoid ambiguous parameter names
  # Args:
  #   names: names to check for ambiguity
  #   effects: output of extract_effects
  auxpars <- intersect(auxpars(), names(effects))
  if (length(auxpars)) {
    auxpars_prefix <- paste0("^", auxpars, "_")
    invalid <- any(ulapply(auxpars_prefix, grepl, names))
    if (invalid) {
      auxpars <- paste0(auxpars, "_", collapse = ", ")
      stop("Variable names starting with ", auxpars,
           " are not allowed for this model.", 
           call. = FALSE)
    }
  }
  invisible(NULL)
}

update_formula <- function(formula, data = NULL, family = gaussian(),
                           nonlinear = NULL, partial = NULL) {
  # incorporate addition arguments and category specific effects into formula 
  # Args:
  #   formula: a model formula 
  #   data: a data.frame or NULL 
  #   partial: a one sided formula containing category specific effects
  # Returns:
  #   an updated formula containing the addition and category specific effects
  formula <- bf(formula, nonlinear = nonlinear)
  old_attributes <- attributes(formula)
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
  attributes(formula) <- old_attributes
  if (is.categorical(family) && is.null(attr(formula, "response"))) {
    respform <- extract_effects(formula)$respform
    model_response <- model.response(model.frame(respform, data = data))
    response <- levels(factor(model_response))
    if (length(response) <= 2L) {
      stop(paste("At least 3 response categories are required",
                 "for categorical models"), call. = FALSE)
    }
    # the first level will serve as the reference category
    attr(formula, "response") <- response[-1]
  }
  formula
}

extract_random <- function(re_terms) {
  # generate a data.frame with all information about the group-level terms
  # Args:
  #   re_terms: A vector of random effects terms in lme4 syntax
  stopifnot(!length(re_terms) || is.character(re_terms))
  lhs_terms <- get_matches("^[^\\|]*", re_terms)
  mid_terms <- get_matches("\\|([^\\|]*\\||)", re_terms)
  rhs_terms <- sub("^\\|", "", get_matches("\\|[^\\|]*$", re_terms))
  random <- vector("list", length(re_terms))
  for (i in seq_along(re_terms)) {
    form <- formula(paste("~", lhs_terms[i]))
    cor <- substr(mid_terms[i], 1, 2) != "||"
    id <- gsub("\\|", "", mid_terms[i])
    if (!nzchar(id)) id <- NA
    groups <- unlist(strsplit(rhs_terms[i], "/", fixed = TRUE))
    new_groups <- c(groups[1], rep("", length(groups) - 1L))
    for (j in seq_along(groups)) {
      if (illegal_group_expr(groups[j])) {
        stop(paste("Illegal grouping term:", rhs_terms[i], "\n may contain",
                   "only variable names combined by the symbols ':' or '/'"),
             call. = FALSE)
      }
      if (j > 1L) {
        new_groups[j] <- paste0(new_groups[j - 1], ":", groups[j])
      }
    }
    random[[i]] <- data.frame(group = new_groups, gn = i, id = id,
                              cor = rep(cor, length(groups)),
                              stringsAsFactors = FALSE)
    random[[i]]$form <- replicate(length(new_groups), form)
  }
  if (length(random)) {
    random <- do.call(rbind, random)
    # ensure that all REs of the same gf are next to each other
    random <- random[order(random$group), ]
  } else {
    random <- data.frame(group = character(0), gn = numeric(0),
                         id = numeric(0), cor = logical(0), 
                         form = character(0))
  }
  random
}

illegal_group_expr <- function(group, bs_valid = TRUE) {
  # check if the group part of a group-level term is invalid
  # Args:
  #  g: the group part of a group-level term
  #  bs_valid: is '/' valid? 
  valid_expr <- ":|[^([:digit:]|[:punct:])][[:alnum:]_\\.]*"
  rsv_signs <- c(".", "+", "*", "|", "?", "::", "\\", if (!bs_valid) "/")
  nchar(gsub(valid_expr, "", group)) ||
    any(ulapply(rsv_signs, grepl, x = group, fixed = TRUE))
}

get_re_terms <- function(x, formula = FALSE, brackets = TRUE) {
  # extract RE terms from a formula of character vector
  # Args:
  #   x: formula or character vector
  #   formula: return a formula containing only ranefs?
  if (is(x, "formula")) {
    x <- gsub("[ \t\r\n]+", "", attr(terms(x), "term.labels"), perl = TRUE)
  }
  re_terms <- x[grepl("\\|", x)]
  if (brackets && length(re_terms)) {
    re_terms <- paste0("(", re_terms, ")")
  } 
  if (formula) {
    if (length(re_terms)) {
      re_terms <- formula(paste("~ 1", collapse("+", re_terms)))
    } else {
      re_terms <- ~ 1
    }
  }
  re_terms
}

check_re_formula <- function(re_formula, formula) {
  # validate the re_formula argument as passed to predict and fitted
  # Args:
  #   re_formula: see predict.brmsfit for documentation
  #   formula: formula to match re_formula with
  # Returns:
  #   updated re_formula containign only terms existent in formula
  old_re_formula <- get_re_terms(formula, formula = TRUE)
  if (is.null(re_formula)) {
    re_formula <- old_re_formula
  } else if (SW(anyNA(re_formula))) {
    re_formula <- ~ 1
  } else {
    re_formula <- get_re_terms(as.formula(re_formula), formula = TRUE)
    new <- extract_effects(re_formula, check_response = FALSE)$random
    old <- extract_effects(old_re_formula, check_response = FALSE)$random
    if (nrow(new)) {
      new_terms <- lapply(new$form, terms)
      found <- rep(FALSE, nrow(new))
      for (i in 1:nrow(new)) {
        group <- new$group[[i]]
        old_terms <- lapply(old$form[old$group == group], terms)
        j <- 1
        while (!found[i] && j <= length(old_terms)) {
          found[i] <- all(attr(new_terms[[i]], "term.labels") %in% 
                       attr(old_terms[[j]], "term.labels")) &&
                      attr(new_terms[[i]], "intercept") <=
                       attr(old_terms[[j]], "intercept")
          j <- j + 1
        }
      }  
      new <- new[found, ]
      if (nrow(new)) {
        forms <- ulapply(new$form, formula2string, rm = 1)
        re_terms <- paste("(", forms, "|", new$group, ")")
        re_formula <- formula(paste("~", paste(re_terms, collapse = "+")))
      } else {
        re_formula <- ~ 1
      }
    } else {
      re_formula <- ~ 1
    }
  }
  re_formula
}

update_re_terms <- function(formula, re_formula = NULL) {
  # remove RE terms in formula and add valid RE terms of re_formula
  # can handle brmsformula objects
  # Args:
  #   formula: model formula to be updated
  #   re_formula: formula containing new RE terms
  # Returns:
  #  a formula with updated RE terms
  old_attr <- attributes(formula)[c("response", "old_mv")]
  spars <- sformula(formula)
  for (ap in setdiff(names(spars), "nonlinear")) {
    spars[[ap]] <- update_re_terms(spars[[ap]], re_formula = re_formula)
  }
  if (!is.null(spars$nonlinear)) {
    # non-linear formulae may cause errors when passed to terms
    # and do not contain random effects anyway
    new_formula <- formula
    spars$nonlinear <- lapply(spars$nonlinear, update_re_terms, 
                              re_formula = re_formula)
  } else {
    re_formula <- check_re_formula(re_formula, formula)
    new_formula <- formula2string(formula)
    old_re_terms <- get_re_terms(formula)
    if (length(old_re_terms)) {
      # make sure that + before random terms are also removed
      rm_terms <- c(paste0("+", old_re_terms), old_re_terms)
      new_formula <- rename(new_formula, rm_terms, "")
      if (grepl("~$", new_formula)) {
        # lhs only formulas are not allowed
        new_formula <- paste(new_formula, "1")
      }
    }
    new_re_terms <- get_re_terms(re_formula)
    new_formula <- paste(c(new_formula, new_re_terms), collapse = "+")
    new_formula <- formula(new_formula)
  }
  attributes(new_formula)[names(spars)] <- spars
  attributes(new_formula)[names(old_attr)] <- old_attr
  environment(new_formula) <- environment(formula)
  new_formula
}

plus_rhs <- function(x) {
  # take the right hand side of a formula and add a +
  if (is.formula(x)) x <- Reduce(paste, deparse(x[[2]]))
  if (length(x) && nchar(x)) paste("+", paste(x, collapse = "+"))
  else "+ 1"
}

get_effect <- function(effects, target = c("fixed", "mono", "cse", "gam"),
                       all = TRUE) {
  # get formulas of certain effects in a list
  # Args:
  #   effects: object returned by extract_effects
  #   target: type of effects to return
  #   all: logical; include effects of nl and auxpars?
  target <- match.arg(target)
  out <- list(effects[[target]])
  if (all) {
    el <- rmNULL(c(effects[auxpars()], effects$nonlinear), recursive = FALSE)
    if (length(el)) {
      out <- c(out, lapply(el, function(par) par[[target]]))
    }
  }
  rmNULL(out)
}

get_random <- function(effects, all = TRUE) {
  # get random effects information in a data.frame
  # Args:
  #   effects: object returned by extract_effects
  #   all: logical; include ranefs of nl and aux parameters?
  if (!is.null(effects$random)) {
    stopifnot(is.data.frame(effects$random))
    out <- effects$random
    nresp <- length(effects$response)
    if (nresp > 1L && nrow(out)) {
      # mv models are also using the 'nlpar' argument
      out <- replicate(nresp, out, simplify = FALSE)
      for (i in seq_len(nresp)) {
        out[[i]]$nlpar <- rep(effects$response[i], nrow(out[[i]]))
      }
      out <- do.call(rbind, out)
    } else {
      out$nlpar <- rep("", nrow(out)) 
    }
  } else {
    out <- NULL
  }
  if (all) {
    el <- rmNULL(c(effects[auxpars()], effects$nonlinear), recursive = FALSE)
    if (length(el)) {
      rand <- do.call(rbind, lapply(el, function(par) par$random))
      # R does not allow duplicated rownames and adds "." causing Stan to fail
      rand$nlpar <- get_matches("^[^\\.]+", rownames(rand))
      out <- rbind(out, rand)
    } 
  }
  out
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
    } else {
      offset <- NULL
    }
  }
  offset
}

get_var_combs <- function(..., alist = list()) {
  # get all variable combinations occuring in elements of ...
  # Args:
  #   ...: character vectors or formulas
  #   alist: a list of character vectors or formulas
  dots <- c(list(...), alist)
  for (i in seq_along(dots)) {
    if (is(dots[[i]], "formula")) {
      dots[[i]] <- attr(terms(dots[[i]]), "term.labels")
    }
    dots[[i]] <- lapply(dots[[i]], function(y) 
      all.vars(parse(text = y)))
  }
  unique(unlist(dots, recursive = FALSE))
}

get_spline_labels <- function(x) {
  # extract labels of splines for GAMMs
  # Args:
  #   x: either a formula or a list containing an element "gam"
  if (is.list(x)) x <- x$gam
  if (is.null(x)) return(NULL)
  term_labels <- rename(attr(terms(x), "term.labels"), " ", "")
  term_labels[grepl("^(s|t2|te|ti)\\(", term_labels)]
}

eval_spline <- function(spline) {
  eval(parse(text = paste0("mgcv::", spline)))
}

amend_terms <- function(x) {
  # amend a terms object (or one that can be coerced to it)
  # to be used in get_model_matrix
  # Args:
  #   x: any R object; if not a formula or terms, NULL is returned
  # Returns:
  #   a (possibly amended) terms object or NULL
  if (is.formula(x) || is(x, "terms")) {
    y <- terms(x)
  } else {
    return(NULL)
  }
  if (isTRUE(attr(x, "forked")) && isTRUE(attr(x, "old_mv"))) {
    # ensure that interactions with main and spec won't
    # cause automatic cell mean coding of factors
    term_labels <- attr(y, "term.labels")
    if (any(grepl("(^|:)(main|spec)($|:)", term_labels))) {
      if (any(grepl("(^|:)trait($|:)", term_labels))) {
        stop(paste("formula may not contain variable 'trait'",
                   "when using variables 'main' or 'spec'"),
             call. = FALSE)
      }
      if (attr(y, "intercept")) {
        stop(paste("formula may not contain an intercept",
                   "when using variables 'main' or 'spec'"),
             call. = FALSE)
      }
      attr(x, "rsv_intercept") <- TRUE
    }
  }
  if (isTRUE(attr(x, "rsv_intercept"))) {
    attr(y, "intercept") <- 1
    attr(y, "rm_intercept") <- TRUE
  }
  y
}

has_intercept <- function(formula) {
  # checks if the formula contains an intercept
  # can handle non-linear formulae
  # Args:
  #   formula: a formula object
  formula <- as.formula(formula)
  try_terms <- try(terms(formula), silent = TRUE)
  if (is(try_terms, "try-error")) {
    out <- FALSE
  } else {
    out <- as.logical(attr(try_terms, "intercept"))
  }
  out
}

has_rsv_intercept <- function(formula) {
  # check if model makes use of the reserved variable 'intercept'
  # can handle non-linear formulae
  # Args:
  #   formula: a formula object
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
  #   formula: a formula containing only the model response
  # Returns:
  #   a vector of names of the response variables (columns)
  formula <- lhs(as.formula(formula))
  all_vars <- all.vars(formula)
  if (!length(all_vars)) {
    stop("formula contains no response variables", call. = FALSE)
  }
  mf <- as.data.frame(named_list(all_vars, values = 1))
  mf <- model.frame(formula, data = mf, na.action = NULL)
  pseudo_resp <- model.response(mf)
  if (is.null(dim(pseudo_resp))) {
    # response is a vector
    response <- all_vars[1]
  } else if (length(dim(pseudo_resp)) == 2L) {
    # response is a matrix
    response <- colnames(pseudo_resp)
    empty_names <- which(!nchar(response))
    if (length(empty_names)) {
      response[empty_names] <- paste0("response", empty_names)
    }
  } else {
    stop("invalid response part of 'formula'", call. = FALSE)
  }
  gsub("\\.|_", "", make.names(response, unique = TRUE))
}

has_splines <- function(effects) {
  # check if splines are present in the model
  # Args:
  #   effects: output of extract_effects
  if (length(effects$nonlinear)) {
    out <- any(ulapply(effects$nonlinear, has_splines))
  } else {
    out <- !is.null(effects[["gam"]])
  }
  ee_auxpars <- rmNULL(effects[auxpars()], recursive = FALSE)
  out || any(ulapply(ee_auxpars, has_splines))
}

gather_ranef <- function(effects, data = NULL, all = TRUE) {
  # gathers helpful information on the random effects
  # Args:
  #   effects: output of extract_effects
  #   data: data passed to brm after updating
  #   all: include REs of non-linear and auxiliary parameters?
  #   combined: should 
  # Returns: 
  #   A named list with one element per grouping factor
  random <- get_random(effects, all = all)
  ranef <- vector("list", nrow(random))
  used_ids <- new_ids <- id_groups <- NULL
  j <- 1
  for (i in seq_len(nrow(random))) {
    Z <- get_model_matrix(random$form[[i]], data = data)
    rdat <- data.frame(id = random$id[[i]], 
                       group = random$group[[i]], 
                       gn = random$gn[[i]],
                       coef = colnames(Z), cn = NA,
                       nlpar = random$nlpar[[i]],
                       cor = random$cor[[i]], 
                       stringsAsFactors = FALSE)
    rdat$form <- replicate(nrow(rdat), random$form[[i]])
    id <- random$id[[i]]
    if (is.na(id)) {
      rdat$id <- j
      j <- j + 1
    } else {
      if (id %in% used_ids) {
        k <- match(id, used_ids)
        rdat$id <- new_ids[k]
        if (!identical(random$group[[i]], id_groups[k])) {
          stop("Can only combine group-level terms of the ",
               "same grouping factor.", call. = FALSE)
        }
      } else {
        used_ids <- c(used_ids, id)
        k <- length(used_ids)
        rdat$id <- new_ids[k] <- j
        id_groups[k] <- random$group[[i]]
        j <- j + 1
      }
    }
    ranef[[i]] <- rdat 
  }
  ranef <- do.call(rbind, c(list(empty_ranef()), ranef))
  if (nrow(ranef)) {
    for (id in unique(ranef$id)) {
      ranef$cn[ranef$id == id] <- seq_len(sum(ranef$id == id))
    }
    levels <- lapply(unique(random$group), function(g) 
      levels(factor(get(g, data))))
    names(levels) <- unique(random$group)
    attr(ranef, "levels") <- levels
  }
  ranef
}

empty_ranef <- function() {
  data.frame(id = numeric(0), group = character(0), gn = numeric(0),
             coef = character(0), cn = numeric(0), nlpar = character(0),
             cor = logical(0), form = character(0), stringsAsFactors = FALSE)
}

rsv_vars <- function(family, nresp = 1, rsv_intercept = FALSE,
                     old_mv = FALSE) {
  # returns names of reserved variables
  # Args:
  #   family: the model family
  #   nresp: number of response variables
  #   rsv_intercept: is the reserved variable "intercept" used?
  if (isTRUE(old_mv)) {
    if (is.linear(family) && nresp > 1L || is.categorical(family)) {
      rsv <- "trait"
    } else if (is.forked(family)) {
      rsv <- c("trait", "main", "spec")
    } else {
      rsv <- NULL
    }
  } else {
    rsv <- NULL
  }
  if (isTRUE(rsv_intercept)) {
    rsv <- c(rsv, "intercept")
  }
  rsv
}

add_families <- function(x) {
  # return names of valid families for addition argument x
  switch(x, weights = "all",
         se = c("gaussian", "student", "cauchy"),
         trials = c("binomial", "zero_inflated_binomial"),
         cat = c("cumulative", "cratio", "sratio", "acat"), 
         cens = c("gaussian", "student", "cauchy", "lognormal",
                  "inverse.gaussian", "binomial", "poisson", 
                  "geometric", "negbinomial", "exponential", 
                  "weibull", "gamma"),
         trunc = c("gaussian", "student", "cauchy", "lognormal", 
                   "binomial", "poisson", "geometric", "negbinomial",
                   "exponential", "weibull", "gamma"),
         disp = c("gaussian", "student", "cauchy", "lognormal", 
                  "gamma", "weibull", "negbinomial"),
         stop(paste("addition argument", x, "is not supported")))
}

formula2string <- function(formula, rm = c(0, 0)) {
  # converts formula to string
  # Args:
  #   formula: a model formula
  #   rm: a vector of to elements indicating how many characters 
  #       should be removed at the beginning
  #       and end of the string respectively
  # Returns:
  #    the formula as string 
  if (!is(formula, "formula")) {
    formula <- as.formula(formula)
  }
  if (is.na(rm[2])) rm[2] <- 0
  x <- gsub("[ \t\r\n]+", "", Reduce(paste, deparse(formula)), perl = TRUE)
  x <- substr(x, 1 + rm[1], nchar(x) - rm[2])
  x
}

get_bounds <- function(trunc, data = NULL) {
   # extract truncation boundaries out of a formula
   # that is known to contain the .trunc function
   # Returns:
   #   a list containing two numbers named lb and ub
   if (is(trunc, "formula")) {
     .addition(trunc, data = data)
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

exclude_pars <- function(effects, ranef = empty_ranef(),
                         save_ranef = TRUE) {
  # list irrelevant parameters NOT to be saved by Stan
  # Args:
  #   effects: output of extract_effects
  #   ranef: output of gather_ranef
  #   save_ranef: should random effects of each level be saved?
  # Returns:
  #   a vector of parameters to be excluded
  out <- c("temp_Intercept1", "temp_Intercept", "Lrescor", 
           "Rescor", "Sigma", "LSigma", "res_cov_matrix", 
           "hs_local", "hs_global",
           intersect(auxpars(), names(effects)))
  # TODO: correctly remove spline helper parameters
  # TODO: correctly remove temporary intercepts
  nlpars <- names(effects$nonlinear)
  # exclude spline helper parameters
  for (i in seq_along(nlpars)) {
    splines <- get_spline_labels(effects$nonlinear[[i]])
    if (length(splines)) {
      out <- c(out, paste0("zs_", nlpars[i], "_", seq_along(splines)))
    }
  }
  splines <- get_spline_labels(effects)
  if (length(splines)) {
    out <- c(out, paste0("zs_", seq_along(splines)))
  }
  # exclude group-level helper parameters
  if (nrow(ranef)) {
    rm_re_pars <- c("z", "L", "Cor", "r")
    for (id in unique(ranef$id)) {
      out <- c(out, paste0(rm_re_pars, "_", id))
    }
    if (!save_ranef) {
      usc_nlpar <- usc(ranef$nlpar, "prefix")
      out <- c(out, paste0("r_", ranef$gn, usc_nlpar, "_", ranef$cn))
    }
  }
  out
}
