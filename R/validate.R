#' Parse Formulas of \pkg{brms} Models
#' 
#' Parse \code{formula} or \code{brmsformula} objects for use 
#' in \pkg{brms}.
#' 
#' @inheritParams brm
#' @param check_response Logical; Indicates whether the left-hand side 
#'   of \code{formula} (i.e. response variables and addition arguments) 
#'   should be parsed. If \code{FALSE}, \code{formula} may also be one-sided.
#' @param resp_rhs_all Logical; Indicates whether to also include response 
#'   variables on the right-hand side of formula \code{.$allvars}, 
#'   where \code{.} represents the output of \code{parse_bf}. 
#'  
#' @return An object of class \code{brmsterms}, which is a \code{list}
#'   containing all required information initially stored in \code{formula} 
#'   in an easier to use format, basically a list of formulas 
#'   (not an abstract syntax tree).
#' 
#' @details This is the main formula parsing function of \pkg{brms}.
#'   It should usually not be called directly, but is exported to allow
#'   package developers making use of the formula syntax implemented 
#'   in \pkg{brms}. As long as no other packages depend on this functions,
#'   it may be changed without deprecation warnings, when new features make
#'   this necessary.
#'   
#' @seealso 
#'   \code{\link[brms:brm]{brm}}, 
#'   \code{\link[brms:brmsformula]{brmsformula}}
#' 
#' @export
parse_bf <- function(formula, family = NULL, autocor = NULL, 
                     check_response = TRUE, resp_rhs_all = TRUE) {
  x <- bf(formula, family = family)
  old_mv <- isTRUE(x[["old_mv"]])
  if (!(is.null(autocor) || is.cor_brms(autocor))) {
    stop2("Argument 'autocor' has to be of class 'cor_brms'")
  }
  
  formula <- x$formula
  family <- x$family
  y <- nlist(formula)
  add_forms <- parse_add(formula, family, check_response)
  add_vars <- str2formula(ulapply(add_forms, all.vars))
  y[names(add_forms)] <- add_forms
  
  pforms <- pforms(x)
  auxpars <- intersect(names(pforms), valid_auxpars(family, y))
  nlpars <- setdiff(names(pforms), auxpars)
  tfixed <- gsub("\\|+[^~]*~", "~", formula2str(formula))
  if (isTRUE(x[["nl"]])) {
    if (!length(nlpars)) {
      stop2("No non-linear parameters specified.")
    }
    if (grepl("|", tfixed, fixed = TRUE)) {
      stop2("Group-level terms cannot be specified ", 
            "in the non-linear formula itself.")
    }
    if (is_ordinal(family) || is_categorical(family)) {
      stop2("Non-linear formulas are not yet allowed for this family.")
    }
    y$fixed <- formula(tfixed)
    y$nlpars <- parse_nonlinear(pforms[nlpars], y$fixed)
    re_terms <- NULL
  } else {
    if (length(nlpars)) {
      nlpars <- collapse_comma(nlpars)
      stop2("Prediction of parameter(s) ", nlpars,
            " is not allowed for this model.")
    }
    terms <- terms(formula)
    all_terms <- all_terms(formula)
    pos_re_terms <- grepl("\\|", all_terms)
    re_terms <- all_terms[pos_re_terms]
    mo_form <- parse_mo(formula)
    if (is.formula(mo_form)) {
      y[["mo"]] <- mo_form
    }
    cs_form <- parse_cs(formula, family = family)
    if (is.formula(cs_form)) {
      y[["cs"]] <- cs_form
    }
    me_form <- parse_me(formula)
    if (is.formula(me_form)) {
      y[["me"]] <- me_form
    }
    gam_form <- parse_gam(formula)
    if (is.formula(gam_form)) {
      y[["gam"]] <- gam_form
    }
    rm_pos <- lapply(list(mo_form, cs_form, me_form, gam_form), attr, "pos")
    rm_pos <- c(rm_pos, list(pos_re_terms))
    fe_terms <- all_terms[!Reduce("|", rm_pos)]
    int_term <- ifelse(attr(terms, "intercept") == 1, "1", "0")
    fe_terms <- paste(c(int_term, fe_terms, get_offset(formula)), 
                      collapse = "+")
    tfixed <- paste(sub("~.*", "", tfixed), "~", fe_terms)
    y$fixed <- formula(tfixed)
    if (is_ordinal(family)) {
      y$fixed <- update.formula(y$fixed, . ~ . + 1)
    }
    if (has_rsv_intercept(y$fixed)) {
      attr(y$fixed, "rsv_intercept") <- TRUE
    }
    if (is_forked(family)) {
      attr(y$fixed, "forked") <- TRUE
    }
  }
  if (check_response && length(y$fixed) < 3L) { 
    stop2("Invalid formula: response variable is missing.")
  }
  
  # make sure to store the plain names of all predictors
  covars <- setdiff(all.vars(rhs(y$fixed)), names(y$nlpars))
  y$covars <- str2formula(covars)
  attr(y$covars, "rsv_intercept") <- TRUE
  # extract group-level effects parts
  y$random <- parse_random(re_terms)
  # evaluate autocor formula
  y$time <- parse_time(autocor$formula)
  # evaluate formulas for auxiliary parameters
  for (ap in auxpars) {
    y$auxpars[[ap]] <- parse_bf(rhs(pforms[[ap]]), check_response = FALSE)
  }
  # fixed auxiliary parameters
  pfix <- pfix(x)
  inv_fauxpars <- setdiff(names(pfix), valid_auxpars(family, y))
  if (length(inv_fauxpars)) {
    stop2("Invalid auxiliary parameters: ", collapse_comma(inv_fauxpars))
  }
  y$fauxpars <- pfix
  
  # make a formula containing all required variables
  lhs_vars <- if (resp_rhs_all) all.vars(lhs(y$fixed))
  fixed_vars <- if (!length(y$nlpars)) rhs(y$fixed)
  formula_list <- c(lhs_vars, add_vars, fixed_vars,
                    y[c("covars", "cs", "mo", "me")], 
                    attr(y$gam, "allvars"), y$random$form,
                    lapply(y$random$gcall, "[[", "allvars"), 
                    lapply(y$nlpars, "[[", "allvars"),
                    lapply(y$auxpars, "[[", "allvars"), 
                    get_offset(y$fixed), y$time$allvars)
  new_formula <- collapse(ulapply(rmNULL(formula_list), plus_rhs))
  all_vars <- c("1", all.vars(parse(text = new_formula)))
  new_formula <- paste(new_formula, "+", paste0(all_vars, collapse = "+"))
  y$allvars <- eval2(paste0("update(", tfixed, ", ~ ", new_formula, ")"))
  environment(y$allvars) <- environment(formula)
  
  # extract response variables
  if (check_response) {
    y$respform <- lhs(y$allvars)
    if (!is.null(x$response)) {
      y$response <- x$response
    } else { 
      y$response <- parse_response(y$respform, keep_dot_usc = old_mv)
    }
    if (is_linear(family) && length(y$response) > 1L) {
      if (length(rmNULL(y[c("se", "cens", "trunc")]))) {
        stop2("Multivariate models currently allow only ",
              "addition argument 'weights'.")
      }
      if (length(y$auxpars)) {
        stop2("Auxiliary parameter cannot yet be ", 
              "predicted in multivariate models.")
      }
    }
    if (old_mv) {
      # multivariate ('trait') syntax is deprecated as of brms 1.0.0
      if (is_hurdle(family)) {
        y$response <- c(y$response, paste0("hu_", y$response))
      } else if (is_zero_inflated(family)) {
        y$response <- c(y$response, paste0("zi_", y$response))
      }
      if (length(y$response) > 1L) {
        # don't use update on a formula that is possibly non-linear
        y$fixed[[2]] <- quote(response)
        y$allvars[[2]] <- quote(response)
      }
      attr(y$formula, "old_mv") <- TRUE
    }
  }
  class(y) <- "brmsterms"
  y
}

parse_add <- function(formula, family = NULL, check_response = TRUE) {
  # extract addition arguments out formula
  # Args:
  #   see parse_bf
  # Returns:
  #   A list of formulas each containg a single addition term
  x <- list()
  add_funs <- lsp("brms", what = "exports", pattern = "^resp_")
  add_funs <- sub("^resp_", "", add_funs)
  if (!is.null(family) && !isNA(family$family)) {
    add <- get_matches("\\|[^~]*~", formula2str(formula))
    if (length(add)) {
      # replace deprecated '|' by '+'
      add <- paste("~", rename(substr(add, 2, nchar(add) - 1), "|", "+"))
      add_terms <- attr(terms(formula(add)), "term.labels")
      for (af in add_funs) {
        matches <- grep(paste0("^(resp_)?", af, "\\(.+\\)$"), add_terms)
        if (length(matches) == 1L) {
          x[[af]] <- add_terms[matches]
          if (!grepl("^resp_", x[[af]])) {
            x[[af]] <- paste0("resp_", x[[af]])
          }
          add_terms <- add_terms[-matches]
          add_fams <- add_families(af)
          valid <- add_fams[1] == "all" || family$family %in% add_fams
          if (!is.na(x[[af]]) && valid) {
            x[[af]] <- str2formula(x[[af]])
          } else {
            stop2("Argument '", af, "' is not supported for ", 
                  "family '", family$family, "'.")
          } 
        } else if (length(matches) > 1L) {
          stop2("Each addition argument may only be defined once.")
        } 
      }
      if (length(add_terms)) {
        stop2("The following addition terms are invalid:\n",
              collapse_comma(add_terms))
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
    if (is_wiener(family) && check_response && !is.formula(x$dec)) {
      stop2("Addition argument 'dec' is required for family 'wiener'.")
    }
  }
  x
}

parse_mo <- function(formula) {
  # extract monotonic terms
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

parse_cs <- function(formula, family = NULL) {
  # category specific terms for ordinal models
  all_terms <- all_terms(formula)
  pos_cs_terms <- grepl("^cse?\\([^\\|]+$", all_terms)
  cs_terms <- all_terms[pos_cs_terms]
  if (length(cs_terms)) {
    if (!is.null(family) && !allows_cs(family)) {
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

parse_gam <- function(formula) {
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

parse_me <- function(formula) {
  # extract variables modeled with measurement error
  # Args:
  #   formula: a formula object
  all_terms <- all_terms(formula)
  # do not include group-level terms
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

parse_random <- function(re_terms) {
  # generate a data.frame with all information about the group-level terms
  # Args:
  #   re_terms: A vector of group-level terms in extended lme4 syntax
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
    random <- random[order(random$group), ]
  } else {
    random <- data.frame(group = character(0), gn = numeric(0),
                         id = numeric(0), cor = logical(0), 
                         type = character(0), form = character(0))
  }
  random
}

parse_response <- function(formula, keep_dot_usc = FALSE) {
  # extract response variable names
  # Args:
  #   formula: a two-sided formula
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

parse_time <- function(formula) {
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

parse_nonlinear <- function(x, model = ~1) {
  # prepare non-linear formulas
  # Args:
  #   x: a list for formulas specifying linear predictors for 
  #      non-linear parameters
  #   model: formula of the non-linear model
  # Returns:
  #   A list of brmsterms objects
  stopifnot(is.list(x), is.formula(model))
  if (length(x)) {
    out <- named_list(names(x))
    for (i in seq_along(x)) {
      out[[i]] <- parse_bf(rhs(x[[i]]), check_response = FALSE)
    }
    model_vars <- all.vars(rhs(model))
    missing_pars <- setdiff(names(out), model_vars)
    if (length(missing_pars)) {
      stop2("Some non-linear parameters are missing in formula: ", 
            collapse_comma(missing_pars))
    }
  } else {
    out <- list() 
  }
  out
}

#' Checks if argument is a \code{brmsterms} object
#' 
#' @param x An \R object
#' 
#' @seealso \code{\link[brms:parse_bf]{parse_bf}}
#' 
#' @export
is.brmsterms <- function(x) {
  inherits(x, "brmsterms")
}

avoid_auxpars <- function(names, bterms) {
  # avoid ambiguous parameter names
  # Args:
  #   names: names to check for ambiguity
  #   bterms: object of class brmsterms
  stopifnot(is.brmsterms(bterms))
  auxpars <- c(names(bterms$auxpars), "mo", "cs", "me")
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
    lhs_form_mo <- parse_mo(lhs_form)
    if (is.formula(lhs_form_mo)) {
      pos_mo <- attr(lhs_form_mo, "pos")
      if (!all(pos_mo)) {
        stop2("Please specify monotonic effects ", 
              "in separate group-level terms.")
      }
      lhs_terms[i] <- formula2str(lhs_form_mo, rm = 1)
      type[[i]] <- "mo"
    }
    lhs_form_cs <- parse_cs(lhs_form)
    if (is.formula(lhs_form_cs)) {
      pos_cs <- attr(lhs_form_cs, "pos")
      if (!all(pos_cs)) {
        stop2("Please specify category specific effects ", 
              "in separate group-level terms.")
      }
      lhs_terms[i] <- formula2str(lhs_form_cs, rm = 1)
      type[[i]] <- "cs"
    }
    lhs_form_me <- parse_me(lhs_form)
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
  structure_not_null(unlist(new_re_terms), type = unlist(type))
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
    new <- parse_bf(re_formula, check_response = FALSE)$random
    old <- parse_bf(old_re_formula, check_response = FALSE)$random
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

update_re_terms <- function(x, re_formula = NULL) {
  # update group-level terms
  # Args:
  #   x: Either 'formula' or 'brmsformula' object
  #   re_formula: formula containing new RE terms
  if (is.formula(x)) {
    x <- .update_re_terms(x, re_formula) 
  } else if (is.brmsformula(x)) {
    if (!x[["nl"]]) {
      x$formula <- .update_re_terms(x$formula, re_formula)
    }
    x$pforms <- lapply(pforms(x), .update_re_terms, re_formula)
  } else {
    stop("Don't know how to handle objects of class ",
         collapse_comma(class(x)))
  }
  x
}

.update_re_terms <- function(formula, re_formula = NULL) {
  # remove existing group-level terms in formula and
  # add valid group-level terms of re_formula
  # Args:
  #   formula: object of class 'formula'
  formula <- as.formula(formula)
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

get_effect <- function(bterms, 
                       target = c("fixed", "mo", "me", "cs", "gam"),
                       all = TRUE) {
  # get formulas of certain effects in a list
  # Args:
  #   bterms: object of class brmsterms
  #   target: type of effects to return
  #   all: logical; include effects of nl and auxpars?
  stopifnot(is.brmsterms(bterms))
  target <- match.arg(target)
  eff <- bterms[[target]]
  out <- list(eff)
  if (all) {
    el <- c(bterms$auxpars, bterms$nlpars)
    if (length(el)) {
      out <- c(out, lapply(el, function(e) e[[target]]))
    }
  }
  rmNULL(out)
}

get_random <- function(bterms, all = TRUE) {
  # get group-level information in a data.frame
  # Args:
  #   bterms: object of class brmsterms
  #   all: logical; include ranefs of nl and aux parameters?
  stopifnot(is.brmsterms(bterms))
  if (!is.null(bterms$random)) {
    stopifnot(is.data.frame(bterms$random))
    out <- bterms$random
    nresp <- length(bterms$response)
    old_mv <- isTRUE(attr(bterms$formula, "old_mv"))
    if (!old_mv && nresp > 1L && nrow(out)) {
      # new MV models are also using the 'nlpar' argument
      out <- replicate(nresp, out, simplify = FALSE)
      for (i in seq_len(nresp)) {
        out[[i]]$nlpar <- rep(bterms$response[i], nrow(out[[i]]))
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
    ep <- c(bterms$auxpars, bterms$nlpars)
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

get_all_effects <- function(bterms, rsv_vars = NULL, 
                            comb_all = FALSE) {
  # get all effects for use in marginal_effects
  # Args:
  #   bterms: object of class brmsterms
  #   rsv_vars: character vector of reserved variables
  #   comb_all: include all main effects and two-way interactions?
  # Returns:
  #   a list with one element per valid effect / effects combination
  #   excludes all 3-way or higher interactions
  stopifnot(is.brmsterms(bterms))
  stopifnot(is.atomic(rsv_vars))
  int_formula <- function(x) {
    formula(paste("~", paste(x, collapse = "*")))
  }
  .get_all_effects <- function(bt) {
    covars <- attr(bt$gam, "covars")
    byvars <- attr(bt$gam, "byvars")
    svars <- mapply(c, covars, byvars, SIMPLIFY = FALSE)
    alist <- lapply(svars, int_formula)
    get_var_combs(bt$fixed, bt$mo, bt$cs, bt$me, alist = alist)
  }
  if (length(bterms$nlpars)) {
    # allow covariates as well as pop-level effects of non-linear parameters
    covars <- setdiff(all.vars(rhs(bterms$fixed)), names(bterms$nlpars))
    covars_comb <- as.list(covars)
    if (length(covars) > 1L) {
      covars_comb <- c(covars_comb, utils::combn(covars, 2, simplify = FALSE))
    } 
    nl_effects <- lapply(bterms$nlpars, .get_all_effects)
    nl_effects <- unlist(nl_effects, recursive = FALSE)
    out <- unique(c(covars_comb, nl_effects))
  } else {
    out <- .get_all_effects(bterms)
  }
  # make sure to also include effects only present in auxpars
  if (length(bterms$auxpars)) {
    ap_effects <- lapply(bterms$auxpars, .get_all_effects)
    ap_effects <- unlist(ap_effects, recursive = FALSE)
    out <- unique(c(out, ap_effects))
  }
  out <- rmNULL(lapply(out, setdiff, y = rsv_vars))
  if (comb_all) {
    out <- unique(unlist(out))
    int <- expand.grid(out, out, stringsAsFactors = FALSE)
    int <- int[int[, 1] != int[, 2], ]
    int <- as.list(as.data.frame(t(int), stringsAsFactors = FALSE))
    int <- unique(unname(lapply(int, sort)))
    out <- c(as.list(out), int)
  }
  out[lengths(out) <= 2L] 
}

get_spline_labels <- function(x, data = NULL, covars = FALSE,
                              combine = TRUE) {
  # extract labels of splines for GAMMs
  # Args:
  #   x: either a formula or a list containing an element "gam"
  #   data: optional data frame containing the covariates
  #   covars: should the names of the covariates be returned
  #           instead of the full term names?
  #   byvars_as_covars: should byvars be returned as covars?
  #   combine: combine names of the covariates (TRUE) 
  #            or just return the covariate names (FALSE)?
  if (is.formula(x)) {
    x <- parse_bf(x, check_response = FALSE)
  }
  if (!is.formula(x$gam)) {
    return(NULL)
  }
  return_covars <- covars
  term_labels <- rename(attr(terms(x$gam), "term.labels"), " ", "")
  splines <- term_labels[grepl("^(s|t2|te|ti)\\(", term_labels)] 
  byvars <- attr(x$gam, "byvars")
  if (return_covars) {
    sfuns <- get_matches("^[^\\(]+", splines)
    covars <- attr(x$gam, "covars")
    for (i in seq_along(covars)) {
      covars[[i]] <- c(covars[[i]], byvars[[i]])
    }
    if (combine) {
      splines <- paste0(sfuns, ulapply(covars, collapse))
    } else {
      splines <- covars
    }
  }
  if (length(splines) && !is.null(data)) {
    # one spline term may contain multiple spline matrices
    if (return_covars && !combine) {
      stop("Invalid combination of arguments. Please report a bug.")
    }
    sdata_fixef <- data_fixef(x, data, knots = attr(data, "knots"))
    by_levels <- attr(sdata_fixef[["X"]], "by_levels")
    nby <- lengths(by_levels)
    splines <- as.list(splines)
    for (i in seq_along(splines)) {
      if (nby[i] > 0L) {
        splines[[i]] <- paste0(splines[[i]], rename(by_levels[[i]]))
      }
    }
    nby[nby == 0L] <- 1L
    termnum <- rep(seq_along(splines), nby)
    splines <- unlist(splines)
    knots <- sdata_fixef[grepl("^knots_", names(sdata_fixef))]
    nbases <- setNames(ulapply(knots, length), splines)
    alist <- nlist(nbases, by_levels, termnum)
    attributes(splines)[names(alist)] <- alist
  }
  structure_not_null(splines, byvars = byvars)
}

eval_spline <- function(spline) {
  eval2(paste0("mgcv::", spline))
}

get_me_labels <- function(x, data) {
  # get labels of measurement error terms
  # Args:
  #   x: either a formula or a list containing an element "gam"
  #   data: data frame containing the noisy variables
  if (is.formula(x)) {
    x <- parse_bf(x, check_response = FALSE)
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

has_splines <- function(bterms) {
  # check if splines are present in the model
  # Args:
  #   bterms: object of class brmsterms
  stopifnot(is.brmsterms(bterms))
  if (length(bterms$nlpars)) {
    out <- any(ulapply(bterms$nlpars, has_splines))
  } else {
    out <- !is.null(bterms[["gam"]])
  }
  out || any(ulapply(bterms$auxpars, has_splines))
}

has_cs <- function(bterms) {
  stopifnot(is.brmsterms(bterms))
  length(all_terms(bterms$cs)) ||
    any(get_random(bterms)$type %in% "cs")
}

tidy_ranef <- function(bterms, data = NULL, all = TRUE, 
                       ncat = NULL, old_levels = NULL) {
  # combines helpful information on the group-level effects
  # Args:
  #   bterms: object of class brmsterms
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
  stopifnot(is.brmsterms(bterms))
  random <- get_random(bterms, all = all)
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
    avoid_auxpars(coef, bterms = bterms)
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
    if (is_linear(family) && nresp > 1L || is_categorical(family)) {
      rsv <- c("trait", "response")
    } else if (is_forked(family)) {
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
  switch(x, 
    weights = "all",
    se = c("gaussian", "student", "cauchy"),
    trials = c("binomial", "zero_inflated_binomial"),
    cat = c("cumulative", "cratio", "sratio", "acat"), 
    cens = c("gaussian", "student", "cauchy", "lognormal",
             "inverse.gaussian", "binomial", "poisson", 
             "geometric", "negbinomial", "exponential", 
             "weibull", "gamma", "exgaussian", "frechet",
             "asym_laplace"),
    trunc = c("gaussian", "student", "cauchy", "lognormal", 
              "binomial", "poisson", "geometric", "negbinomial",
              "exponential", "weibull", "gamma", "inverse.gaussian",
              "exgaussian", "frechet", "asym_laplace"),
    disp = c("gaussian", "student", "cauchy", "lognormal", 
             "gamma", "weibull", "negbinomial", "exgaussian",
             "asym_laplace"),
    dec = c("wiener"),
    stop2("Addition argument '", x, "' is not supported.")
  )
}

get_bounds <- function(formula, data = NULL) {
  # extract truncation boundaries
  # Returns:
  #   a list containing two numbers named lb and ub
  if (is.formula(formula)) {
    term <- attr(terms(formula), "term.labels")
    stopifnot(length(term) == 1L && grepl("resp_trunc\\(", term))
    trunc <- eval_rhs(formula, data = data)
  } else {
    trunc <- resp_trunc()
  }
  trunc
}

has_cens <- function(formula, data = NULL) {
  # indicate if the model is (possibly interval) censored
  if (is.formula(formula)) {
    term <- attr(terms(formula), "term.labels")
    stopifnot(length(term) == 1L && grepl("resp_cens\\(", term))
    cens <- eval_rhs(formula, data = data)
    cens <- structure(TRUE, interval = !is.null(attr(cens, "y2")))
  } else {
    cens <- FALSE
  }
  cens
}

exclude_pars <- function(bterms, data = NULL, ranef = empty_ranef(),
                         save_ranef = TRUE, save_mevars = FALSE) {
  # list irrelevant parameters NOT to be saved by Stan
  # Args:
  #   bterms: object of class brmsterms
  #   data: data passed by the user
  #   ranef: output of tidy_ranef
  #   save_ranef: should group-level effects be saved?
  #   save_mevars: should samples of noise-free variables be saved?
  # Returns:
  #   a vector of parameters to be excluded
  stopifnot(is.brmsterms(bterms))
  out <- c("temp_Intercept1", "temp_Intercept", "Lrescor", "Rescor",
           "Sigma", "LSigma", "res_cov_matrix", "hs_local",
           intersect(auxpars(), names(bterms)))
  # exclude spline helper parameters and temporary Intercepts
  par_effects <- c(bterms$auxpars, bterms$nlpars)
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
  splines <- get_spline_labels(bterms, data)
  if (length(splines) && !is.null(data)) {
    for (i in seq_along(splines)) {
      nb <- seq_len(attr(splines, "nbases")[[i]])
      out <- c(out, paste0("zs_", i, "_", nb))
    }
  }
  meef <- get_me_labels(bterms, data)
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
