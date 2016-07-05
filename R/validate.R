extract_effects <- function(formula, ..., family = NA, nonlinear = NULL, 
                            check_response = TRUE, resp_rhs_all = TRUE,
                            lhs_char = "") {
  # Parse the model formula and related arguments
  # Args:
  #   formula: An object of class 'formula' or 'brmsformula'
  #   ...: Additional objects of class "formula"
  #   family: the model family
  #   nonlinear: a list of formulas specifying non-linear effects
  #   check_response: check if the response part is non-empty
  #   resp_rhs_all: include response variables on the RHS of $all? 
  #   lhs_char: the response part of the model as a character string;
  #             currently only used for splines in non-linear models
  # Returns: 
  #   A named list whose elements depend on the formula input 
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
    rsv_pars <- names(sformula(formula))
    x$nonlinear <- nonlinear_effects(nonlinear, model = x$fixed,
                                     family = family, rsv_pars = rsv_pars)
    re_terms <- NULL
  } else {
    # terms() doesn't like non-linear formulas
    re_terms <- get_re_terms(formula)
    if (length(re_terms)) {
      tfixed <- rename(tfixed, c(paste0("+", re_terms), re_terms), "")
    } 
    # monotonic effects
    term_labels <- gsub(" ", "", attr(terms(formula), "term.labels"))
    mono_terms <- term_labels[grepl("^mono(|tonic|tonous)\\(", term_labels)]
    if (length(mono_terms)) {
      tfixed <- rename(tfixed, c(paste0("+", mono_terms), mono_terms), "")
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
      attr(cse_terms, "rsv_intercept") <- TRUE
      if (!length(all.vars(cse_terms))) {
        stop("invalid input to function 'cse'", call. = FALSE)
      }
      x$cse <- cse_terms
    }
    # parse spline expression for GAMMs
    splines <- term_labels[grepl("^(s|t2|te|ti)\\(", term_labels)]
    if (length(splines)) {
      if (is.mv(family) || is.forked(family) || is.categorical(family)) {
        stop("Splines are not yet implemented for this family.", 
             call. = FALSE)
      }
      if (any(grepl("^(te|ti)\\(", splines))) {
        stop(paste("Tensor product splines 'te' and 'ti' are not yet", 
                   "implemented in brms. Consider using 't2' instead."),
             call. = FALSE)
      }
      tfixed <- rename(tfixed, c(paste0("+", splines), splines), "")
      if (!nchar(lhs_char)) {
        lhs_char <- get_matches("^[^~]*", tfixed)
      }
      stopifnot(nchar(lhs_char) > 0L)
      x$gam <- formula(paste(lhs_char, "~", paste(splines, collapse = "+")))
    }
    if (substr(tfixed, nchar(tfixed), nchar(tfixed)) == "~") {
      tfixed <- paste0(tfixed, "1")
    }
    if (grepl("|", x = tfixed, fixed = TRUE)) {
      stop("Random effects terms should be enclosed in brackets", call. = FALSE)
    }
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
  # TODO include validity checks
  auxpars <- sformula(formula, incl_nl = FALSE)
  for (ap in names(auxpars)) {
    x[[ap]] <- extract_effects(auxpars[[ap]], family = family,
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
    for (i in seq_along(x$nonlinear)) {
      # currently only required for GAMMs
      x$nonlinear[[i]]$respform <- x$respform
    }
    if (!is.null(attr(formula, "response"))) {
      x$response <- attr(formula, "response")
    } else { 
      x$response <- gather_response(x$respform)
    }
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
  time <- all.vars(as.formula(paste("~", gsub("~|\\|[[:print:]]*", "", formula))))
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

nonlinear_effects <- function(x, model = ~1, family = NA, rsv_pars = NULL) {
  # prepare nonlinear formulas
  # Args:
  #   x: a list for formulas specifying linear predictors for 
  #      non-linear parameters
  #   model: formula of the non-linear model
  # Returns:
  #   A list of objects each returned by extract_effects
  stopifnot(is.list(x), is(model, "formula"))
  if (length(x)) {
    lhs_char <- as.character(model[[2]])
    nleffects <- vector("list", length = length(x))
    for (i in seq_along(x)) {
      x[[i]] <- as.formula(x[[i]])
      if (length(x[[i]]) != 3L) {
        stop("Non-linear formulas must be two-sided.", call. = FALSE)
      }
      nlresp <- all.vars(x[[i]][[2]])
      if (length(nlresp) != 1L) {
        stop("LHS of non-linear formula must contain exactly one variable.",
             call. = FALSE)
      }
      if (any(ulapply(c(".", "_"), grepl, x = nlresp, fixed = TRUE))) {
        stop("Non-linear parameters should not contain dots or underscores.",
             call. = FALSE)
      }
      if (nlresp %in% rsv_pars) {
        stop("Parameter name '", nlresp, "' is reserved for this model.",
             call. = FALSE)
      }
      if (!is.null(attr(terms(x[[i]]), "offset"))) {
        stop("Offsets are currently not allowed in non-linear models.",
             call. = FALSE)
      }
      x[[i]] <- rhs(x[[i]])
      nleffects[[i]] <- extract_effects(x[[i]], family = family, 
                                        check_response = FALSE,
                                        lhs_char = lhs_char)
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
  x
}

update_formula <- function(formula, data = NULL, family = gaussian(),
                           nonlinear = NULL, partial = NULL) {
  # incorporate addition arguments and category specific effects into formula 
  # 
  # Args:
  #   formula: a model formula 
  #   data: a data.frame or NULL 
  #   partial: a one sided formula containing category specific effects
  #
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
    response <- levels(as.factor(model_response))
    if (length(response) <= 2L) {
      stop(paste("At least 3 response categories are required",
                 "for categorical models"), call. = FALSE)
    }
    attr(formula, "response") <- response
  }
  formula
}

#' Set up a model formula for use in the \pkg{brms} package
#' 
#' Set up a model formula for use in the \pkg{brms} package
#' allowing to define additive multilevel models for all parameters
#' of the assumed distribution of the response.
#' 
#' @inheritParams brm
#' @param ... Additional \code{formula} objects to specify 
#'   predictors of special model parts and auxiliary parameters. 
#'   Formulas can either be named directly or contain
#'   names on their left-hand side. Currently, the following
#'   names corresponding are accepted: 
#'   \code{sigma} (residual standard deviation of
#'   the \code{gaussian} and \code{student} families);
#'   \code{shape} (shape parameter of the \code{Gamma},
#'   \code{weibull}, \code{negbinomial} and related
#'   zero-inflated / hurdle families); \code{nu}
#'   (degrees of freedom parameter of the \code{student} family);
#'   \code{phi} (precision parameter of the \code{beta} 
#'   and \code{zero_inflated_beta} families).
#' 
#' @export
bf <- function(formula, ..., nonlinear = NULL) {
  # parse and validate dots arguments
  dots <- list(...)
  parnames <- names(dots)
  if (is.null(parnames)) {
    parnames <- rep("", length(dots))
  }
  for (i in seq_along(dots)) {
    dots[[i]] <- as.formula(dots[[i]])
    if (length(dots[[i]]) == 3L && !nchar(parnames[i])) {
      resp_par <- all.vars(dots[[i]][[2]])
      if (length(resp_par) != 1L) {
        stop("LHS of additional formulas must contain exactly one variable.",
             call. = FALSE)
      }
      parnames[i] <- resp_par
    }
    if (!is.null(attr(terms(dots[[i]]), "offset"))) {
      stop("Offsets are currently not allowed.", call. = FALSE)
    }
    dots[[i]] <- rhs(dots[[i]])
  }
  names(dots) <- parnames
  if (any(!nchar(names(dots)))) {
    stop("Function 'bf' requires named arguments.", call. = FALSE)
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
  nonlinear <- nonlinear2list(nonlinear)
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

sformula <- function(x, incl_nl = TRUE, ...) {
  # special formulas stored in brmsformula objects
  # Args:
  #   x: coerced to a 'brmsformula' object
  #   incl_nl: include the 'nonlinear' argument in the output?
  rmNULL(attributes(bf(x))[auxpars(incl_nl = incl_nl)])
}

auxpars <- function(incl_nl = FALSE) {
  auxpars <- c("sigma", "shape", "nu", "phi")
  if (incl_nl) {
    auxpars <- c(auxpars, "nonlinear")
  }
  auxpars
}

valid_auxpars <- function(family, effects = list(), autocor = cor_arma()) {
  # convenience function to find relevant auxiliary parameters
  x <- c(sigma = has_sigma(family, effects = effects, autocor = autocor),
         shape = has_shape(family), nu = has_nu(family), phi = has_phi(family))
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

#' @export
update.brmsformula <- function(object, formula., 
                               mode = c("update", "replace", "keep"), 
                               ...) {
  # update a brmsformula and / or its attributes
  # Args:
  #   object: an object of class 'brmsformula'
  #   formula.: formula to update object
  #   mode: "update": apply update.formula
  #         "replace": replace old formula
  #         "keep": keep old formula
  #   ...: currently unused
  mode <- match.arg(mode)
  new_att <- sformula(formula.)
  old_att <- sformula(object)
  if (mode == "update") {
    new_formula <- update.formula(object, formula., ...)
  } else if (mode == "replace") {
    new_formula <- formula.
  } else {
    new_formula <- object
  }
  attributes(new_formula)[union(names(new_att), names(old_att))] <- NULL
  new_formula <- do.call(bf, c(new_formula, new_att))
  new_formula <- SW(do.call(bf, c(new_formula, old_att)))
  new_formula
}

extract_random <- function(re_terms) {
  # generate a data.frame with all information about the group-level terms
  # Args:
  #   re_terms: A vector of random effects terms in lme4 syntax
  stopifnot(!length(re_terms) || is.character(re_terms))
  lhs_terms <- get_matches("\\([^\\|]*", re_terms)
  rhs_terms <- get_matches("\\|[^\\)]*", re_terms)
  random <- vector("list", length(re_terms))
  for (i in seq_along(re_terms)) {
    form <- formula(paste("~", substring(lhs_terms[i], 2)))
    cor <- substr(rhs_terms[i], 1, 2) != "||"
    rhs_terms[i] <- sub("^\\|*", "", rhs_terms[i])
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
    random[[i]] <- data.frame(group = new_groups, 
                              cor = rep(cor, length(groups)),
                              stringsAsFactors = FALSE)
    random[[i]]$form <- replicate(length(new_groups), form)
  }
  if (length(random)) {
    random <- do.call(rbind, random)
    # ensure that all REs of the same gf are next to each other
    # to allow combining them in rename_pars later on
    random <- random[order(random$group), ]
  } else {
    random <- data.frame(group = character(0), 
                         cor = logical(0), 
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

get_re_terms <- function(x, formula = FALSE) {
  # extract RE terms from a formula of character vector
  # Args:
  #   x: formula or character vector
  #   formula: return a formula containing only ranefs?
  if (is(x, "formula")) {
    x <- gsub(" ", "", attr(terms(x), "term.labels"))
  }
  re_terms <- x[grepl("\\|", x)]
  if (length(re_terms)) {
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
  environment(new_formula) <- environment(formula)
  new_formula
}

plus_rhs <- function(x) {
  # take the right hand side of a formula and add a +
  if (is.formula(x)) x <- Reduce(paste, deparse(x[[2]]))
  if (length(x) && nchar(x)) paste("+", paste(x, collapse = "+"))
  else "+ 1"
}

get_effect <- function(effects, target = c("fixed", "mono", "cse", "gam")) {
  # get fixed effects formulas in a list
  # Args:
  #   effects: object returned by extract_effects
  target <- match.arg(target)
  out <- list(effects[[target]])
  if (!is.null(effects$nonlinear)) {
    out <- c(out, lapply(effects$nonlinear, function(par) par[[target]]))
    attr(out, "nonlinear") <- TRUE
  }
  out
}

get_random <- function(effects, all = TRUE) {
  # get random effects information in a data.frame
  # Args:
  #   effects: object returned by extract_effects
  #   all: logical; include ranefs of nl and aux parameters?
  if (!is.null(effects$random)) {
    stopifnot(is.data.frame(effects$random))
    out <- effects$random
    out$nlpar <- rep("", nrow(out))
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
      attr(out, "nonlinear") <- TRUE
    } 
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
    # each non-linear parameter may have its own random effects
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

amend_terms <- function(x, forked = FALSE) {
  # amend a terms object (or one that can be coerced to it)
  # to be used in get_model_matrix
  # Args:
  #   x: any R object; if not a formula or terms, NULL is returned
  #   forked: a flag indicating if the model is forked into
  #           two parts (e.g., a hurdle model).
  # Returns:
  #   a (possibly amended) terms object or NULL
  if (is.formula(x) || is(x, "terms")) {
    y <- terms(x)
  } else {
    return(NULL)
  }
  if (forked) {
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

check_intercept <- function(names) {
  # check if model contains fixed effects intercept
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

gather_ranef <- function(effects, data = NULL, all = TRUE,
                         combine = FALSE) {
  # gathers helpful information on the random effects
  # Args:
  #   effects: output of extract_effects
  #   data: data passed to brm after updating
  #   all: include REs of non-linear and auxiliary parameters?
  # Returns: 
  #   A named list with one element per grouping factor
  random <- get_random(effects, all = all)
  Z <- lapply(random$form, get_model_matrix, data = data, 
              forked = isTRUE(attr(effects$fixed, "forked")))
  ranef <- setNames(lapply(Z, colnames), random$group)
  for (i in seq_along(ranef)) {
    attr(ranef[[i]], "levels") <- 
      levels(factor(get(random$group[[i]], data)))
    attr(ranef[[i]], "group") <- names(ranef)[i]
    attr(ranef[[i]], "cor") <- random$cor[[i]]
    attr(ranef[[i]], "nlpar") <- random$nlpar[i]
  }
  if (combine) {
    ranef <- combine_duplicates(ranef, sep = c("nlpar", "levels"))
  }
  ranef
}

rsv_vars <- function(family, nresp = 1, rsv_intercept = NULL) {
  # returns names of reserved variables
  # Args:
  #   family: the model family
  #   nresp: number of response variables
  #   rsv_intercept: is the reserved variable "intercept" used?
  if (is.forked(family)) {
    rsv <- c("trait", "main", "spec")
  } else if (is.linear(family) && nresp > 1L || is.categorical(family)) {
    rsv <- "trait"
  } else {
    rsv <- NULL
  }
  if (isTRUE(rsv_intercept)) {
    rsv <- c(rsv, "intercept")
  }
  rsv
}

check_mv_formula <- function(family, effects) {
  # check if reserved mv variables were (correctly) used
  # TODO: add more / better checks
  # Args:
  #   family: the model family
  #   effects: output of extract effects
  rsv_vars <- rsv_vars(family, nresp = length(effects$response))
  if (!is.null(rsv_vars)) {
    if (!any(rsv_vars %in% all.vars(rhs(effects$formula)))) {
      warning(paste(
        "Apparently, you did not use any of the variables specific",
        "to multivariate models. \nThis will lead to parameter",
        "estimates which are pooled over all model parts, \nwhich is",
        "probably not what you wanted. \nPlease see the Details section",
        "of help(brm) for more information."), call. = FALSE)
    }
  } 
  invisible(NULL)
}

cse <- function(...) {
  stop("inappropriate use of function 'cse'", call. = FALSE)
}

monotonic <- function(...) {
  stop("inappropriate use of function 'monotonic'", call. = FALSE)
}

mono <- function(...) {
  # abbreviation of monotonic
  stop("inappropriate use of function 'monotonic'", call. = FALSE)
}

monotonous <- function(...) {
  # abbreviation of monotonic
  stop("please use function 'monotonic' instead", call. = FALSE)
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

is.formula <- function(x, or = TRUE) {
  # checks if x is formula (or list of formulas)
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
  # Args:
  #   formula: a model formula
  #   rm: a vector of to elements indicating how many characters 
  #       should be removed at the beginning
  #       and end of the string respectively
  # Returns:
  #    the formula as string 
  if (!is(formula, "formula")) {
    stop(paste(deparse(substitute(formula)), "must be of class formula"))
  }
  if (is.na(rm[2])) rm[2] <- 0
  x <- gsub(" ", "", Reduce(paste, deparse(formula)))
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

exclude_pars <- function(effects, ranef = list(), 
                         save_ranef = TRUE) {
  # list irrelevant parameters NOT to be saved by Stan
  # Args:
  #   effects: output of extract_effects
  #   ranef: output of gather_ranef
  #   save_ranef: should random effects of each level be saved?
  # Returns:
  #   a vector of parameters to be excluded
  out <- c("eta", "etap", "eta_2PL", "Eta", "temp_Intercept1", 
           "temp_Intercept",  "Lrescor", "Rescor", "Sigma", 
           "LSigma", "disp_sigma", "e", "E", "res_cov_matrix", 
           "lp_pre", "hs_local", "hs_global",
           intersect(auxpars(), names(effects)))
  nlpars <- names(effects$nonlinear)
  if (length(nlpars)) {
    out <- c(out, unique(paste0("eta_", nlpars)))
    for (i in seq_along(nlpars)) {
      splines <- get_spline_labels(effects$nonlinear[[i]])
      if (length(splines)) {
        out <- c(out, paste0("zs_", nlpars[i], "_", seq_along(splines)))
      }
    }
  } else {
    splines <- get_spline_labels(effects)
    if (length(splines)) {
      out <- c(out, paste0("zs_", seq_along(splines)))
    }
  }
  if (length(ranef)) {
    rm_re_pars <- c("z", "L", "Cor", if (!save_ranef) "r")
    # names of NL-parameters must be computed based on ranef here
    nlp <- ulapply(ranef, attr, "nlpar")
    if (length(nlp)) {
      stopifnot(length(nlp) == length(ranef))
      nlp <- ifelse(nchar(nlp), paste0(nlp, "_"), nlp)
      for (k in seq_along(ranef)) {
        i <- which(which(nlp == nlp[k]) == k)
        out <- c(out, paste0(rm_re_pars, "_", nlp[k], i))
        neff <- length(ranef[[k]])
        if (neff > 1L) {
          out <- c(out, paste0("r_", nlp[k], i, "_", 1:neff))
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
