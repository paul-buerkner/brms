extract_effects <- function(formula, family = NA, nonlinear = NULL, 
                            autocor = NULL, check_response = TRUE, 
                            resp_rhs_all = TRUE, ...) {
  # Parse the model formula and related arguments
  # Args:
  #   formula: An object of class 'formula' or 'brmsformula'
  #   family: the model family
  #   nonlinear: a list of formulas specifying non-linear effects
  #   autocor: object of class 'cor_brms' or NULL
  #   check_response: check if the response part is non-empty?
  #   resp_rhs_all: include response variables on the RHS of x$allvars? 
  #   ...: Additional objects of class 'formula'
  # Returns: 
  #   A named list whose elements depend on the formula input
  old_mv <- isTRUE(attr(formula, "old_mv"))
  formula <- bf(formula, nonlinear = nonlinear)
  nonlinear <- attr(formula, "nonlinear")
  if (!is.na(family[[1]])) {
    family <- check_family(family)
  }
  tformula <- formula2str(formula) 
  tfixed <- gsub("\\|+[^~]*~", "~", tformula)
  x <- nlist(formula)
  if (length(nonlinear)) {
    if (grepl("|", tfixed, fixed = TRUE)) {
      stop2("Group-level effects in non-linear models have ", 
            "to be specified in the 'nonlinear' argument.")
    }
    if (is.ordinal(family) || is.categorical(family)) {
      stop2("Non-linear effects are not yet allowed for this family.")
    }
    x$fixed <- formula(tfixed)
    x$nonlinear <- extract_nonlinear(nonlinear, x$fixed, family = family)
    re_terms <- NULL
  } else {
    # terms() doesn't like non-linear formulas
    terms <- terms(formula)
    all_terms <- all_terms(formula)
    pos_re_terms <- grepl("\\|", all_terms)
    re_terms <- all_terms[pos_re_terms]
    mo_form <- extract_mo(formula)
    if (is.formula(mo_form)) {
      x[["mo"]] <- mo_form
    }
    cs_form <- extract_cs(formula, family = family)
    if (is.formula(cs_form)) {
      x[["cs"]] <- cs_form
    }
    me_form <- extract_me(formula)
    if (is.formula(me_form)) {
      x[["me"]] <- me_form
    }
    gam_form <- extract_gam(formula)
    if (is.formula(gam_form)) {
      x[["gam"]] <- gam_form
    }
    rm_pos <- lapply(list(mo_form, cs_form, me_form, gam_form), attr, "pos")
    rm_pos <- c(rm_pos, list(pos_re_terms))
    fe_terms <- all_terms[!Reduce("|", rm_pos)]
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
    stop2("Invalid formula: response variable is missing.")
  }
  
  # make sure to store the plain names of all predictors
  covars <- setdiff(all.vars(rhs(x$fixed)), names(x$nonlinear))
  x$covars <- formula(paste("~", paste(c("1", covars), collapse = "+")))
  attr(x$covars, "rsv_intercept") <- TRUE
  # extract group-level effects parts
  x$random <- extract_random(re_terms)
  # evaluate autocor formula
  x$time <- extract_time(autocor$formula)
  # evaluate formulas for auxiliary parameters
  auxpars <- sformula(formula, incl_nl = FALSE)
  for (ap in names(auxpars)) {
    x[[ap]] <- extract_effects(rhs(auxpars[[ap]]), family = family,
                               check_response = FALSE)
  }
  
  # handle addition arguments
  add_funs <- c("se", "weights", "trials", "cat", 
                "cens", "trunc", "disp", "dec")
  add_vars <- list()
  if (!is.na(family[[1]])) {
    add <- get_matches("\\|[^~]*~", tformula)
    if (length(add)) {
      # replace deprecated '|' by '+'
      add <- paste("~", rename(substr(add, 2, nchar(add) - 1), "|", "+"))
      add_terms <- attr(terms(formula(add)), "term.labels")
      for (af in add_funs) {
        matches <- grep(paste0("^", af, "\\(.+\\)$"), add_terms)
        if (length(matches) == 1L) {
          x[[af]] <- add_terms[matches]
          add_terms <- add_terms[-matches]
          if (!is.na(x[[af]]) && (add_families(af)[1] == "all" ||
              family$family %in% add_families(af))) {
            args <- substr(x[[af]], nchar(af) + 2, nchar(x[[af]]) - 1)
            try_numeric <- SW(as.numeric(args))
            if (af %in% c("trials", "cat") && !is.na(try_numeric)) {
              x[[af]] <- try_numeric
            } else {
              x[[af]] <- as.formula(paste0("~ .", x[[af]]))
              if (length(all.vars(x[[af]]))) {
                form <- paste("~", paste(all.vars(x[[af]]), collapse = "+"))
                add_vars[[af]] <- as.formula(form)
              }
            }
          } else {
            stop2("Argument '", af, "' is not supported for ", 
                  "family '", family$family, "'.")
          } 
        } else if (length(matches) > 1L) {
          stop2("Addition arguments may be only defined once.")
        } 
      }
      if (length(add_terms)) {
        stop2("Invalid addition part of formula. \nPlease see the ", 
              "'Details' section of help(brmsformula).")
      }
      if (is.formula(x$se) && is.formula(x$disp)) {
        stop2("Addition arguments 'se' and 'disp' cannot ", 
              "be used at the same time.")
      }
      cens_or_weights <- is.formula(x$cens) || is.formula(x$weights)
      if (is.formula(x$trunc) && cens_or_weights) {
        stop2("Truncation is not yet possible in censored or weighted models.")
      }
    }
    if (is.wiener(family) && check_response && !is.formula(x$dec)) {
      stop2("Addition argument 'dec' is required for family 'wiener'.")
    }
  }
  
  # make a formula containing all required variables
  lhs_vars <- if (resp_rhs_all) all.vars(lhs(x$fixed))
  fixed_vars <- if (!length(x$nonlinear)) rhs(x$fixed)
  formula_list <- c(lhs_vars, add_vars, fixed_vars,
                    x[c("covars", "cs", "mo", "me")], 
                    attr(x$gam, "allvars"), x$random$form,
                    lapply(x$random$gcall, "[[", "allvars"), 
                    lapply(x$nonlinear, "[[", "allvars"),
                    lapply(x[auxpars()], "[[", "allvars"), 
                    get_offset(x$fixed), x$time$allvars)
  new_formula <- collapse(ulapply(rmNULL(formula_list), plus_rhs))
  all_vars <- c("1", all.vars(parse(text = new_formula)))
  new_formula <- paste(new_formula, "+", paste0(all_vars, collapse = "+"))
  x$allvars <- eval2(paste0("update(", tfixed, ", ~ ", new_formula, ")"))
  environment(x$allvars) <- environment(formula)
  
  # extract response variables
  if (check_response) {
    x$respform <- lhs(x$allvars)
    if (!is.null(attr(formula, "response"))) {
      x$response <- attr(formula, "response")
    } else { 
      x$response <- extract_response(x$respform, keep_dot_usc = old_mv)
    }
    if (is.linear(family) && length(x$response) > 1L &&
        length(rmNULL(x[c("se", "cens", "trunc")]))) {
      stop2("Multivariate models currently allow only ",
            "addition argument 'weights'.")
    }
    if (old_mv) {
      # multivariate ('trait') syntax is deprecated as of brms 1.0.0
      if (is.hurdle(family)) {
        x$response <- c(x$response, paste0("hu_", x$response))
      } else if (is.zero_inflated(family)) {
        x$response <- c(x$response, paste0("zi_", x$response))
      }
      if (length(x$response) > 1L) {
        # don't use update on a formula that is possibly non-linear
        x$fixed[[2]] <- quote(response)
        x$allvars[[2]] <- quote(response)
        attr(x$formula, "old_mv") <- TRUE
      }
    }
  }
  # check validity of auxiliary parameters
  if (!is.na(family[[1]])) {
    # don't do this earlier to allow exlusion of MV models
    inv_auxpars <- setdiff(auxpars(), valid_auxpars(family, x))
    inv_auxpars <- intersect(inv_auxpars, names(x))
    if (length(inv_auxpars)) {
      inv_auxpars <- paste0("'", inv_auxpars, "'", collapse = ", ")
      stop2("Prediction of parameter(s) ", inv_auxpars,
            " is not allowed for this model.")
    }
  }
  x
} 

extract_mo <- function(formula) {
  # extract monotonic effects
  # Args:
  #   formula: a formula object
  all_terms <- all_terms(formula)
  pos_mo_terms <- grepl("^mo((no)?|(notonic)?)\\([^\\|]+$", all_terms)
  mo_terms <- all_terms[pos_mo_terms]
  if (length(mo_terms)) {
    mo_terms <- formula(paste("~", paste(eval2(mo_terms), collapse = "+")))
    if (!length(all.vars(mo_terms))) {
      stop2("No variable supplied to function 'mo'.")
    }
    attr(mo_terms, "rsv_intercept") <- TRUE
  }
  structure(mo_terms, pos = pos_mo_terms)
}

extract_cs <- function(formula, family = NA) {
  # category specific effects in ordinal models
  all_terms <- all_terms(formula)
  pos_cs_terms <- grepl("^cse?\\([^\\|]+$", all_terms)
  cs_terms <- all_terms[pos_cs_terms]
  if (length(cs_terms)) {
    if (!is.na(family[[1]]) && !allows_cs(family)) {
      stop2("Category specific effects are only meaningful for ", 
            "families 'sratio', 'cratio', and 'acat'.")
    }
    cs_terms <- formula(paste("~", paste(eval2(cs_terms), collapse = "+")))
    # do not test whether variables were supplied to 'cs'
    # to allow category specific group-level intercepts
    attr(cs_terms, "rsv_intercept") <- TRUE
  }
  structure(cs_terms, pos = pos_cs_terms)
}

extract_gam <- function(formula) {
  # parse spline expression for GAMMs
  all_terms <- all_terms(formula)
  pos_gam_terms <- grepl("^(s|t2|te|ti)\\(", all_terms)
  gam_terms <- all_terms[pos_gam_terms]
  if (length(gam_terms)) {
    if (any(grepl("^(te|ti)\\(", gam_terms))) {
      stop2("Tensor product splines 'te' and 'ti' are not yet ", 
            "implemented in brms. Consider using 't2' instead.")
    }
    covars <- byvars <- named_list(gam_terms)
    for (i in seq_along(gam_terms)) {
      es <- eval_spline(gam_terms[i])
      covars[[i]] <- es$term
      if (es$by != "NA") {
        byvars[[i]] <- es$by 
      }
    }
    gam_terms <- formula(paste("~", paste(gam_terms, collapse = "+")))
    allvars <- mgcv::interpret.gam(gam_terms)$fake.formula
  } else {
    covars <- byvars <- NULL
    allvars <- ~ 1
  }
  structure(gam_terms, pos = pos_gam_terms, 
            covars = covars, byvars = byvars,  
            allvars = allvars)
}

extract_me <- function(formula) {
  # extract variables modeled with measurement error
  # Args:
  #   formula: a formula object
  all_terms <- all_terms(formula)
  # stop group-level terms from being included
  all_terms[grepl("\\|", all_terms)] <- ""
  pos_me_terms <- grepl_expr("^me\\([^:]*\\)$", all_terms)
  me_terms <- all_terms[pos_me_terms]
  if (length(me_terms)) {
    me_terms <- formula(paste("~", paste(me_terms, collapse = "+")))
    if (!length(all.vars(me_terms))) {
      stop2("No variable supplied to function 'me'.")
    }
    attr(me_terms, "rsv_intercept") <- TRUE
  }
  structure(me_terms, pos = pos_me_terms)
}

extract_random <- function(re_terms) {
  # generate a data.frame with all information about the group-level terms
  # Args:
  #   re_terms: A vector of group-level effects terms in extended lme4 syntax
  re_terms <- split_re_terms(re_terms)
  lhs_terms <- lhs_terms(re_terms)
  mid_terms <- mid_terms(re_terms)
  rhs_terms <- rhs_terms(re_terms)
  random <- vector("list", length(re_terms))
  type <- attr(re_terms, "type")
  for (i in seq_along(re_terms)) {
    id <- gsub("\\|", "", mid_terms[i])
    if (!nzchar(id)) id <- NA
    gcall <- eval2(rhs_terms[i])
    group <- paste0(gcall$type, collapse(gcall$groups))
    random[[i]] <- data.frame(group = group, gtype = gcall$type,
                              gn = i, id = id, type = type[i],
                              cor = substr(mid_terms[i], 1, 2) != "||",
                              stringsAsFactors = FALSE)
    random[[i]]$gcall <- list(gcall)
    random[[i]]$form <- list(formula(paste("~", lhs_terms[i])))
  }
  if (length(random)) {
    random <- do.call(rbind, random)
    # ensure that all REs of the same gf are next to each other
    random <- random[order(random$group), ]
  } else {
    random <- data.frame(group = character(0), gn = numeric(0),
                         id = numeric(0), cor = logical(0), 
                         type = character(0), form = character(0))
  }
  random
}

extract_response <- function(formula, keep_dot_usc = FALSE) {
  # extract response variable names
  # Args:
  #   formula: a formula containing only the model response
  #   keep_dot_usc: keep dots and underscores in the names?
  # Returns:
  #   a vector of names of the response variables (columns)
  formula <- lhs(as.formula(formula))
  all_vars <- all.vars(formula)
  if (!length(all_vars)) {
    stop2("The formula contains no response variables.")
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
    stop2("Response part of 'formula' is invalid.")
  }
  response <- make.names(response, unique = TRUE)
  if (!keep_dot_usc) {
    response <- gsub("\\.|_", "", response)
  }
  response
}

extract_time <- function(formula) {
  # extract time and grouping variables for autocorrelation structures
  # Args:
  #   formula: a one sided formula of the form ~ time|group 
  #            typically taken from a cor_brms object
  # Returns: 
  #   a list with elements time, group, and all, where all contains a 
  #   formula with all variables in formula
  if (is.null(formula)) {
    formula <- ~ 1
  }
  formula <- as.formula(formula)
  if (!is.null(lhs(formula))) {
    stop2("Autocorrelation formula must be one-sided.")
  }
  formula <- formula2str(formula)
  time <- as.formula(paste("~", gsub("~|\\|[[:print:]]*", "", formula)))
  time <- all.vars(time)
  if (length(time) > 1L) {
    stop2("Autocorrelation structures may only contain 1 time variable.")
  }
  x <- list(time = ifelse(length(time), time, ""))
  group <- sub("^\\|*", "", sub("~[^\\|]*", "", formula))
  if (illegal_group_expr(group)) {
    stop2("Illegal grouping term: ", group, "\nIt may contain only ", 
          "variable names combined by the symbol ':'")
  }
  group <- formula(paste("~", ifelse(nchar(group), group, "1")))
  x$group <- paste0(all.vars(group), collapse = ":")
  x$allvars <- paste(c("1", time, all.vars(group)), collapse = "+")
  x$allvars <- formula(paste("~", x$allvars))
  x
}

extract_nonlinear <- function(x, model = ~1, family = NA) {
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
      stop2("Some non-linear parameters are missing in formula: ", 
            paste0("'", missing_pars, "'", collapse = ", "))
    }
  } else {
    nleffects <- list() 
  }
  nleffects
}

avoid_auxpars <- function(names, effects) {
  # avoid ambiguous parameter names
  # Args:
  #   names: names to check for ambiguity
  #   effects: output of extract_effects
  auxpars <- c(auxpars(), "mo", "cs")
  auxpars <- intersect(auxpars, names(effects))
  if (length(auxpars)) {
    auxpars_prefix <- paste0("^", auxpars, "_")
    invalid <- any(ulapply(auxpars_prefix, grepl, names))
    if (invalid) {
      auxpars <- paste0("'", auxpars, "_'", collapse = ", ")
      stop2("Variable names starting with ", auxpars,
            " are not allowed for this model.")
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
  if (!is.null(partial)) {
    warning2("Argument 'partial' is deprecated. Please use the 'cs' ", 
             "function inside the model formula instead.")
    partial <- formula2str(partial, rm = 1)
    fnew <- paste(fnew, "+ cs(", partial, ")")
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
      stop2("At least 3 response categories are required ",
            "for family 'categorical'.")
    }
    # the first level will serve as the reference category
    attr(formula, "response") <- response[-1]
  }
  formula
}

illegal_group_expr <- function(group) {
  # check if the group part of a group-level term is invalid
  # Args:
  #  g: the group part of a group-level term
  valid_expr <- ":|[^([:digit:]|[:punct:])][[:alnum:]_\\.]*"
  rsv_signs <- c("+", "-", "*", "/", "|", "::")
  nzchar(gsub(valid_expr, "", group)) ||
    any(ulapply(rsv_signs, grepl, x = group, fixed = TRUE))
}

get_re_terms <- function(x, formula = FALSE, brackets = TRUE) {
  # extract RE terms from a formula of character vector
  # Args:
  #   x: formula or character vector
  #   formula: return a formula containing only ranefs?
  if (is(x, "formula")) {
    x <- all_terms(x)
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

split_re_terms <- function(re_terms) {
  # split nested group-level terms by the '/' sign
  # and check for special effects terms
  # Args:
  #   re_terms: group-level terms in extended lme4 syntax
  stopifnot(!length(re_terms) || is.character(re_terms))
  lhs_terms <- lhs_terms(re_terms)
  mid_terms <- mid_terms(re_terms)
  rhs_terms <- rhs_terms(re_terms)
  new_re_terms <- vector("list", length(re_terms))
  type <- as.list(rep("", length(re_terms)))
  for (i in seq_along(re_terms)) {
    # check for special terms
    lhs_form <- formula(paste("~", lhs_terms[i]))
    lhs_all_terms <- all_terms(lhs_form)
    lhs_form_mo <- extract_mo(lhs_form)
    if (is.formula(lhs_form_mo)) {
      pos_mo <- attr(lhs_form_mo, "pos")
      if (!all(pos_mo)) {
        stop2("Please specify monotonic effects ", 
              "in separate group-level terms.")
      }
      lhs_terms[i] <- formula2str(lhs_form_mo, rm = 1)
      type[[i]] <- "mo"
    }
    lhs_form_cs <- extract_cs(lhs_form)
    if (is.formula(lhs_form_cs)) {
      pos_cs <- attr(lhs_form_cs, "pos")
      if (!all(pos_cs)) {
        stop2("Please specify category specific effects ", 
              "in separate group-level terms.")
      }
      lhs_terms[i] <- formula2str(lhs_form_cs, rm = 1)
      type[[i]] <- "cs"
    }
    lhs_form_me <- extract_me(lhs_form)
    if (is.formula(lhs_form_me)) {
      pos_me <- attr(lhs_form_me, "pos")
      if (!all(pos_me)) {
        stop2("Please specify terms of noisy variables ", 
              "in separate group-level terms.")
      }
      lhs_terms[i] <- formula2str(lhs_form_me, rm = 1)
      type[[i]] <- "me"
    }
    # expand grouping factor terms
    groups <- terms(formula(paste0("~", rhs_terms[i])))
    groups <- attr(groups, "term.labels")
    groups <- ifelse(!grepl("^(gr|mm)\\(", groups), 
                     paste0("gr(", groups, ")"), groups)
    new_re_terms[[i]] <- paste0(lhs_terms[i], mid_terms[i], groups)
    type[[i]] <- rep(type[[i]], length(new_re_terms[[i]]))
  }
  structure(unlist(new_re_terms), type = unlist(type))
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
        forms <- ulapply(new$form, formula2str, rm = 1)
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
    # and do not contain group-level effects anyway
    new_formula <- formula
    spars$nonlinear <- lapply(spars$nonlinear, update_re_terms, 
                              re_formula = re_formula)
  } else {
    re_formula <- check_re_formula(re_formula, formula)
    new_formula <- formula2str(formula)
    old_re_terms <- get_re_terms(formula)
    if (length(old_re_terms)) {
      # make sure that + before group-level terms are also removed
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
  if (is.formula(x)) {
    x <- Reduce(paste, deparse(rhs(x)[[2]]))
  }
  if (length(x) && all(nzchar(x))) {
    out <- paste0(" + ", paste(x, collapse = "+"))
  } else {
    out <- " + 1"
  }
  out
}

get_effect <- function(effects, 
                       target = c("fixed", "mo", "me", "cs", "gam"),
                       all = TRUE) {
  # get formulas of certain effects in a list
  # Args:
  #   effects: object returned by extract_effects
  #   target: type of effects to return
  #   all: logical; include effects of nl and auxpars?
  target <- match.arg(target)
  eff <- effects[[target]]
  out <- list(eff)
  if (all) {
    el <- c(effects[auxpars()], effects$nonlinear)
    el <- rmNULL(el, recursive = FALSE)
    if (length(el)) {
      out <- c(out, lapply(el, function(e) e[[target]]))
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
    old_mv <- isTRUE(attr(effects$formula, "old_mv"))
    if (!old_mv && nresp > 1L && nrow(out)) {
      # new MV models are also using the 'nlpar' argument
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
    # non-linear and auxiliary parameters
    ep <- rmNULL(c(effects[auxpars()], effects$nonlinear), recursive = FALSE)
    if (length(ep)) {
      rand <- do.call(rbind, lapply(ep, function(e) get_random(e)))
      # R does not allow duplicated rownames and adds "."
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
    dots[[i]] <- lapply(dots[[i]], function(y) all.vars(parse(text = y)))
  }
  unique(unlist(dots, recursive = FALSE))
}

get_all_effects <- function(effects, rsv_vars = NULL) {
  # get all effects for use in marginal_effects
  # Args:
  #   effects: output of extract_effects
  #   rsv_vars: character vector of reserved variables
  # Returns:
  #   a list with one element per valid effect / effects combination
  #   excludes all 3-way or higher interactions
  stopifnot(is.list(effects))
  stopifnot(is.atomic(rsv_vars))
  .get_all_effects <- function(ee) {
    alist <- lapply(attr(ee$gam, "covars"), function(x) 
      formula(paste("~", paste(x, collapse = "*"))))
    get_var_combs(ee$fixed, ee$mo, ee$cs, ee$me, alist = alist)
  }
  if (length(effects$nonlinear)) {
    # allow covariates as well as fixed effects of non-linear parameters
    covars <- setdiff(all.vars(rhs(effects$fixed)), names(effects$nonlinear))
    covars_comb <- as.list(covars)
    if (length(covars) > 1L) {
      covars_comb <- c(covars_comb, utils::combn(covars, 2, simplify = FALSE))
    } 
    nl_effects <- lapply(effects$nonlinear, .get_all_effects)
    nl_effects <- unlist(nl_effects, recursive = FALSE)
    all_effects <- unique(c(covars_comb, nl_effects))
  } else {
    all_effects <- .get_all_effects(effects)
  }
  # make sure to also include effects only present in auxpars
  effects_auxpars <- rmNULL(effects[auxpars()])
  if (length(effects_auxpars)) {
    ap_effects <- lapply(effects_auxpars, .get_all_effects)
    ap_effects <- unlist(ap_effects, recursive = FALSE)
    all_effects <- unique(c(all_effects, ap_effects))
  }
  all_effects <- rmNULL(lapply(all_effects, setdiff, y = rsv_vars))
  all_effects[lengths(all_effects) <= 2L]
}

get_spline_labels <- function(x, data = NULL, covars = FALSE,
                              combine = TRUE) {
  # extract labels of splines for GAMMs
  # Args:
  #   x: either a formula or a list containing an element "gam"
  #   data: optional data frame containing the covariates
  #   covars: should the names of the covariates be returned
  #           instead of the full term names?
  #   combine: combine names of the covariates (TRUE) 
  #            or just return the covariate names (FALSE)?
  if (is.formula(x)) {
    x <- extract_effects(x, check_response = FALSE)
  }
  if (!is.formula(x$gam)) {
    return(NULL)
  }
  term_labels <- rename(attr(terms(x$gam), "term.labels"), " ", "")
  splines <- term_labels[grepl("^(s|t2|te|ti)\\(", term_labels)] 
  if (covars) {
    sfuns <- get_matches("^[^\\(]+", splines)
    covars <- attr(x$gam, "covars")
    byvars <- attr(x$gam, "byvars")
    var_labels <- named_list(attr(x$gam, "covars"))
    for (i in seq_along(var_labels)) {
      var_labels[[i]] <- c(covars[[i]], byvars[[i]])
    }
    if (combine) {
      splines <- paste0(sfuns, ulapply(var_labels, collapse))
    } else {
      splines <- var_labels
    }
  }
  if (length(splines) && !is.null(data)) {
    # one spline term may contain multiple spline matrices
    sdata_fixef <- data_fixef(x, data, knots = attr(data, "knots"))
    knots <- sdata_fixef[grepl("^knots_", names(sdata_fixef))]
    attr(splines, "nbases") <- setNames(ulapply(knots, length), splines)
  }
  splines
}

eval_spline <- function(spline) {
  eval2(paste0("mgcv::", spline))
}

get_me_labels <- function(x, data) {
  if (is.formula(x)) {
    x <- extract_effects(x, check_response = FALSE)
  }
  if (!is.formula(x$me)) {
    return(NULL)
  }
  mm <- get_model_matrix(x$me, data, rename = FALSE)
  not_one <- apply(mm, 2, function(x) any(x != 1))
  uni_me <- get_matches_expr("^me\\([^:]*\\)$", colnames(mm))
  uni_me <- unique(gsub("[[:space:]]", "", uni_me))
  structure(colnames(mm), not_one = not_one, uni_me = uni_me)
}

all_terms <- function(formula) {
  if (is.null(formula)) return(NULL)
  terms <- terms(as.formula(formula))
  gsub("[ \t\r\n]+", "", x = attr(terms, "term.labels"), perl = TRUE)
}

lhs_terms <- function(re_terms) {
  out <- get_matches("^[^\\|]*", re_terms) 
  if (length(out) != length(re_terms)) {
    stop2("One or more group-levels terms were invalid.")
  }
  out
}

mid_terms <- function(re_terms) {
  out <- get_matches("\\|([^\\|]*\\||)", re_terms)
  if (length(out) != length(re_terms)) {
    stop2("One or more group-levels terms were invalid.")
  }
  out
}

rhs_terms <- function(re_terms) {
  out <- sub("^\\|", "", get_matches("\\|[^\\|]*$", re_terms))
  if (length(out) != length(re_terms)) {
    stop2("One or more group-levels terms were invalid.")
  }
  out
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
        stop2("formula may not contain variable 'trait' when ",
              "using variables 'main' or 'spec'")
      }
      if (attr(y, "intercept")) {
        stop2("formula may not contain an intercept when ",
              "using variables 'main' or 'spec'")
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

has_cs <- function(effects) {
  length(all_terms(effects$cs)) ||
    any(get_random(effects)$type %in% "cs")
}

tidy_ranef <- function(effects, data = NULL, all = TRUE, 
                       ncat = NULL, old_levels = NULL) {
  # combines helpful information on the group-level effects
  # Args:
  #   effects: output of extract_effects
  #   data: data passed to brm after updating
  #   all: include REs of non-linear and auxiliary parameters?
  #   ncat: optional number of response categories
  #         only used for category specific group-level effects
  # Returns: 
  #   A tidy data.frame with the following columns:
  #     id: ID of the group-level effect 
  #     group: name of the grouping factor
  #     gn: number of the grouping term within the respective formula
  #     coef: name of the group-level effect
  #     cn: number of the effect within the ID
  #     nlpar: name of the corresponding non-linear parameter
  #     cor: are correlations modeled for this effect?
  #     type: special effects type; can be "mo", "cs", or "me"
  #     gcall: output of functions 'gr' or 'mm'
  #     form: formula used to compute the effects
  random <- get_random(effects, all = all)
  ranef <- vector("list", nrow(random))
  used_ids <- new_ids <- NULL
  id_groups <- list()
  j <- 1
  for (i in seq_len(nrow(random))) {
    if (random$type[[i]] == "mo") {
      coef <- prepare_mo_vars(random$form[[i]], data, check = FALSE)
      coef <- colnames(coef)
    } else if (random$type[[i]] == "cs") {
      coef <- colnames(get_model_matrix(random$form[[i]], data = data))
      if (is.null(ncat)) {
        # try to infer ncat from the data
        Y <- as.numeric(model.response(data))
        ncat <- max(Y) - min(Y) + 1
      }
      indices <- paste0("[", seq_len(ncat - 1), "]")
      coef <- as.vector(t(outer(coef, indices, paste0)))
    } else if (random$type[[i]] == "me") {
      coef <- rename(get_me_labels(random$form[[i]], data))
    } else {
      coef <- colnames(get_model_matrix(random$form[[i]], data = data)) 
    }
    avoid_auxpars(coef, effects = effects)
    rdat <- data.frame(id = random$id[[i]],
                       group = random$group[[i]],
                       gn = random$gn[[i]],
                       gtype = random$gtype[[i]],
                       coef = coef, cn = NA,
                       nlpar = random$nlpar[[i]],
                       cor = random$cor[[i]],
                       type = random$type[[i]],
                       stringsAsFactors = FALSE)
    rdat$gcall <- replicate(nrow(rdat), random$gcall[i]) 
    rdat$form <- replicate(nrow(rdat), random$form[[i]])
    id <- random$id[[i]]
    if (is.na(id)) {
      rdat$id <- j
      j <- j + 1
    } else {
      if (id %in% used_ids) {
        k <- match(id, used_ids)
        rdat$id <- new_ids[k]
        new_id_groups <- c(random$group[[i]], random$gcall[[i]]$groups)
        if (!identical(new_id_groups, id_groups[[k]])) {
          stop2("Can only combine group-level terms of the ",
                "same grouping factors.")
        }
      } else {
        used_ids <- c(used_ids, id)
        k <- length(used_ids)
        rdat$id <- new_ids[k] <- j
        id_groups[[k]] <- c(random$group[[i]], random$gcall[[i]]$groups)
        j <- j + 1
      }
    }
    ranef[[i]] <- rdat 
  }
  ranef <- do.call(rbind, c(list(empty_ranef()), ranef))
  # check for overlap between different group types
  rsv_groups <- ranef[nzchar(ranef$gtype), "group"]
  other_groups <- ranef[!nzchar(ranef$gtype), "group"]
  inv_groups <- intersect(rsv_groups, other_groups)
  if (length(inv_groups)) {
    inv_groups <- paste0("'", inv_groups, "'", collapse = ", ")
    stop2("Grouping factor names ", inv_groups, " are resevered.")
  }
  # check for duplicated and thus not identified effects
  dup <- duplicated(ranef[, c("group", "coef", "nlpar")])
  if (any(dup)) {
    stop2("Duplicated group-level effects are not allowed.")
  }
  if (nrow(ranef)) {
    for (id in unique(ranef$id)) {
      ranef$cn[ranef$id == id] <- seq_len(sum(ranef$id == id))
    }
    if (is.null(old_levels)) {
      un_random <- random[!duplicated(random$group), ]
      levels <- named_list(un_random$group)
      for (i in seq_along(levels)) {
        # combine levels of all grouping factors within one grouping term
        levels[[i]] <- ulapply(un_random$gcall[[i]]$groups, 
                               function(g) levels(factor(get(g, data))))
        levels[[i]] <- unique(levels[[i]])
      }
      attr(ranef, "levels") <- levels 
    } else {
      # for newdata numeration has to depend on the original levels
      attr(ranef, "levels") <- old_levels
    }
  }
  ranef
}

empty_ranef <- function() {
  data.frame(id = numeric(0), group = character(0), gn = numeric(0),
             coef = character(0), cn = numeric(0), nlpar = character(0),
             cor = logical(0), type = character(0), form = character(0), 
             stringsAsFactors = FALSE)
}

rsv_vars <- function(family, nresp = 1L, rsv_intercept = FALSE,
                     old_mv = FALSE) {
  # returns names of reserved variables
  # Args:
  #   family: the model family
  #   nresp: number of response variables
  #   rsv_intercept: is the reserved variable "intercept" used?
  if (isTRUE(old_mv)) {
    if (is.linear(family) && nresp > 1L || is.categorical(family)) {
      rsv <- c("trait", "response")
    } else if (is.forked(family)) {
      rsv <- c("trait", "response", "main", "spec")
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
                  "weibull", "gamma", "exgaussian"),
         trunc = c("gaussian", "student", "cauchy", "lognormal", 
                   "binomial", "poisson", "geometric", "negbinomial",
                   "exponential", "weibull", "gamma", "inverse.gaussian",
                   "exgaussian"),
         disp = c("gaussian", "student", "cauchy", "lognormal", 
                  "gamma", "weibull", "negbinomial", "exgaussian"),
         dec = c("wiener"),
         stop2(paste("Addition argument '", x, "' is not supported.")))
}

get_bounds <- function(formula, data = NULL) {
  # extract truncation boundaries out of a formula
  # that is known to contain the .trunc function
  # Returns:
  #   a list containing two numbers named lb and ub
  if (is.formula(formula)) {
    term <- attr(terms(formula), "term.labels")
    stopifnot(length(term) == 1L && grepl("\\.trunc\\(", term))
    trunc <- eval_rhs(formula, data = data)
  } else {
    trunc <- .trunc()
  }
  trunc
}

has_cens <- function(formula, data = NULL) {
  # indicate if the model is (possibly interval) censored
  if (is.formula(formula)) {
    term <- attr(terms(formula), "term.labels")
    stopifnot(length(term) == 1L && grepl("\\.cens\\(", term))
    cens <- eval_rhs(formula, data = data)
    cens <- structure(TRUE, interval = !is.null(attr(cens, "y2")))
  } else {
    cens <- FALSE
  }
  cens
}

check_brm_input <- function(x) {
  # misc checks on brm arguments 
  # Args:
  #   x: A named list
  if (x$chains %% x$cluster != 0L) {
    stop2("'chains' must be a multiple of 'cluster'")
  }
  family <- check_family(x$family) 
  if (family$family == "inverse.gaussian") {
    warning2("Inverse gaussian models require carefully chosen ", 
             "prior distributions to ensure convergence of the chains.")
  }
  if (family$link == "sqrt") {
    warning2(family$family, " model with sqrt link may not be ", 
             "uniquely identified")
  }
  invisible(NULL)
}

exclude_pars <- function(effects, data = NULL, ranef = empty_ranef(),
                         save_ranef = TRUE, save_mevars = FALSE) {
  # list irrelevant parameters NOT to be saved by Stan
  # Args:
  #   effects: output of extract_effects
  #   data: data passed by the user
  #   ranef: output of tidy_ranef
  #   save_ranef: should group-level effects be saved?
  #   save_mevars: should samples of noise-free variables be saved?
  # Returns:
  #   a vector of parameters to be excluded
  out <- c("temp_Intercept1", "temp_Intercept", "Lrescor", "Rescor",
           "Sigma", "LSigma", "res_cov_matrix", "hs_local",
           intersect(auxpars(), names(effects)))
  # exclude spline helper parameters and temporary Intercepts
  par_effects <- rmNULL(c(effects[auxpars()], effects$nonlinear))
  for (par in names(par_effects)) {
    out <- c(out, paste0("temp_", par, "_Intercept"))
    splines <- get_spline_labels(par_effects[[par]], data)
    if (length(splines) && !is.null(data)) {
      for (i in seq_along(splines)) {
        nb <- seq_len(attr(splines, "nbases")[[i]])
        out <- c(out, paste0("zs_", par, "_", i, "_", nb))
      } 
    }
    meef <- get_me_labels(par_effects[[par]], data)
    if (!save_mevars && length(meef)) {
      out <- c(out, paste0("Xme_", par, "_", seq_along(meef)))
    }
  }
  splines <- get_spline_labels(effects, data)
  if (length(splines) && !is.null(data)) {
    for (i in seq_along(splines)) {
      nb <- seq_len(attr(splines, "nbases")[[i]])
      out <- c(out, paste0("zs_", i, "_", nb))
    }
  }
  meef <- get_me_labels(effects, data)
  if (!save_mevars && length(meef)) {
    out <- c(out, paste0("Xme_", seq_along(meef)))
  }
  # exclude group-level helper parameters
  if (nrow(ranef)) {
    rm_re_pars <- c("z", "L", "Cor", "r")
    for (id in unique(ranef$id)) {
      out <- c(out, paste0(rm_re_pars, "_", id))
    }
    if (!save_ranef) {
      usc_nlpar <- usc(ranef$nlpar, "prefix")
      out <- c(out, paste0("r_", ranef$id, usc_nlpar, "_", ranef$cn))
    }
  }
  out
}
