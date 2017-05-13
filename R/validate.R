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
  y <- nlist(formula, family, autocor) 
  class(y) <- "brmsterms"
  
  if (check_response) {
    # extract response variables
    y$respform <- lhs(formula)
    if (is.null(y$respform)) {
      stop2("Invalid formula: response variable is missing.")
    }
    y$respform <- formula(gsub("\\|+[^~]*~", "~", formula2str(y$respform)))
    if (!is.null(x$response)) {
      y$response <- x$response
    } else { 
      y$response <- parse_resp(y$respform, keep_dot_usc = old_mv)
    }
  }
  
  # extract addition arguments
  adforms <- parse_ad(formula, family, check_response)
  advars <- str2formula(ulapply(adforms, all.vars))
  y$adforms[names(adforms)] <- adforms
  
  str_rhs_form <- formula2str(rhs(formula))
  terms <- try(terms(rhs(formula)), silent = TRUE)
  has_terms <- is(terms, "try-error") || 
    length(attr(terms, "term.labels")) ||
    length(attr(terms, "offset"))
  rhs_needed <- FALSE
  if (is.mixfamily(family)) {
    for (i in seq_along(family$mix)) {
      mui <- paste0("mu", i)
      if (!is.formula(x$pforms[[mui]])) {
        x$pforms[[mui]] <- eval2(paste0(mui, str_rhs_form))
        rhs_needed <- TRUE
      }
    }
  } else {
    if (!is.formula(x$pforms[["mu"]])) { 
      x$pforms[["mu"]] <- eval2(paste0("mu", str_rhs_form))
      rhs_needed <- TRUE
    }
    x$pforms <- x$pforms[c("mu", setdiff(names(x$pforms), "mu"))]
  }
  if (!rhs_needed && has_terms) {
    stop2("All 'mu' parameters are specified so that ",
          "the right-hand side of 'formula' is unused.")
  }
  auxpars <- is_auxpar_name(names(x$pforms), family, bterms = y)
  auxpars <- names(x$pforms)[auxpars]
  # amend when generalizing non-linear models to auxiliary parameters
  nlpars <- setdiff(names(x$pforms), auxpars)
  if (isTRUE(x[["nl"]])) {
    if (is.mixfamily(family) || is_ordinal(family) || is_categorical(family)) {
      stop2("Non-linear formulas are not yet allowed for this family.")
    }
    y$auxpars[["mu"]] <- parse_nlf(x$pforms[["mu"]], x$pforms[nlpars])
    y$auxpars[["mu"]]$family <- auxpar_family(family, "mu")
    auxpars <- setdiff(auxpars, "mu")
  } else {
    if (length(nlpars)) {
      stop2("Parameter '", nlpars[1], "' is not part of the model.")
    }
  }
  
  # predicted auxiliary parameters
  for (ap in auxpars) {
    y$auxpars[[ap]] <- parse_lf(x$pforms[[ap]], family = family)
    y$auxpars[[ap]]$family <- auxpar_family(family, ap)
  }
  if (!is.null(y$auxpars[["mu"]])) {
    y$auxpars[["mu"]][["autocor"]] <- autocor
  }
  # fixed auxiliary parameters
  inv_fauxpars <- setdiff(names(x$pfix), valid_auxpars(family, y))
  if (length(inv_fauxpars)) {
    stop2("Invalid auxiliary parameters: ", collapse_comma(inv_fauxpars))
  }
  y$fauxpars <- x$pfix
  check_fauxpars(y$fauxpars)
  # check for illegal use of cs terms
  if (has_cs(y) && !(is.null(family) || allows_cs(family))) {
    stop2("Category specific effects are only meaningful for ", 
          "families 'sratio', 'cratio', and 'acat'.")
  }
  # parse autocor formula
  y$time <- parse_time(autocor$formula)
  
  # make a formula containing all required variables
  lhsvars <- if (resp_rhs_all) all.vars(y$respform)
  lformula <- c(
    lhsvars, advars, 
    lapply(y$auxpars, "[[", "allvars"), 
    y$time$allvars
  )
  y$allvars <- allvars_formula(lformula)
  if (check_response) {
    y$allvars <- update(y$respform, y$allvars) 
  }
  environment(y$allvars) <- environment(formula)
  
  if (is_linear(family) && length(y$response) > 1L) {
    if (any(!names(y$adforms) %in% "weights")) {
      stop2("Multivariate models currently allow only ",
            "addition argument 'weights'.")
    }
    if (length(y$auxpars) > 1L) {
      stop2("Auxiliary parameters cannot yet be ", 
            "predicted in multivariate models.")
    }
  }
  if (check_response && old_mv) {
    # multivariate ('trait') syntax is deprecated as of brms 1.0.0
    if (is_hurdle(family)) {
      y$response <- c(y$response, paste0("hu_", y$response))
    } else if (is_zero_inflated(family)) {
      y$response <- c(y$response, paste0("zi_", y$response))
    }
    if (length(y$response) > 1L) {
      y$allvars[[2]] <- quote(response)
    }
    attr(y$formula, "old_mv") <- TRUE
  }
  y
}

parse_lf <- function(formula, family = NULL) {
  # parse linear formulas
  # Args:
  #   formula: an ordinary R formula
  #   family: the model family
  # Returns:
  #   object of class 'btl'
  formula <- rhs(as.formula(formula))
  y <- nlist(formula)
  
  terms <- terms(formula)
  all_terms <- all_terms(formula)
  pos_re_terms <- grepl("\\|", all_terms)
  re_terms <- all_terms[pos_re_terms]
  mo_form <- parse_mo(formula)
  if (is.formula(mo_form)) {
    y[["mo"]] <- mo_form
  }
  cs_form <- parse_cs(formula)
  if (is.formula(cs_form)) {
    y[["cs"]] <- cs_form
  }
  me_form <- parse_me(formula)
  if (is.formula(me_form)) {
    y[["me"]] <- me_form
  }
  sm_form <- parse_sm(formula)
  if (is.formula(sm_form)) {
    y[["sm"]] <- sm_form
  }
  gp_form <- parse_gp(formula)
  if (is.formula(gp_form)) {
    y[["gp"]] <- gp_form
  }
  offset_form <- parse_offset(formula)
  if (is.formula(offset_form)) {
    y[["offset"]] <- offset_form
  }
  rm_pos <- list(mo_form, cs_form, me_form, sm_form, gp_form)
  rm_pos <- c(lapply(rm_pos, attr, "pos"), list(pos_re_terms))
  fe_terms <- all_terms[!Reduce("|", rm_pos)]
  int_term <- ifelse(attr(terms, "intercept") == 1, "1", "0")
  fe_terms <- paste(c(int_term, fe_terms), collapse = "+")
  y[["fe"]] <- str2formula(fe_terms)
  if (is_ordinal(family)) {
    y[["fe"]] <- update.formula(y[["fe"]], . ~ . + 1)
  }
  if (has_rsv_intercept(y[["fe"]])) {
    attr(y[["fe"]], "rsv_intercept") <- TRUE
  }
  if (is_forked(family)) {
    attr(y[["fe"]], "forked") <- TRUE
  }
  # parse group-level terms
  y$re <- parse_re(re_terms)
  lformula <- c(
    y[c("fe", "cs", "mo", "me")], 
    attr(y$sm, "allvars"), attr(y$gp, "allvars"),
    y$re$form, lapply(y$re$gcall, "[[", "allvars"),
    str2formula(all.vars(y$offset))
  )
  y$allvars <- allvars_formula(lformula)
  environment(y$allvars) <- environment(formula)
  structure(y, class = "btl")
}

parse_nlf <- function(formula, nlpar_forms) {
  # parse non-linear formulas
  # Args:
  #   formula: formula of the non-linear model
  #   nlpar_forms: a list for formulas specifying linear predictors 
  #                for non-linear parameters
  # Returns:
  #   object of class 'btnl'
  stopifnot(is.list(nlpar_forms))
  formula <- rhs(as.formula(formula))
  y <- list()
  if (length(nlpar_forms)) {
    y$formula <- formula
    nlpars <- names(nlpar_forms)
    y$nlpars <- named_list(nlpars)
    for (nlp in nlpars) {
      y$nlpars[[nlp]] <- parse_lf(nlpar_forms[[nlp]])
    }
    model_vars <- all.vars(formula)
    missing_pars <- setdiff(names(y$nlpars), model_vars)
    if (length(missing_pars)) {
      stop2("Some non-linear parameters are missing in formula: ", 
            collapse_comma(missing_pars))
    }
    covars <- str2formula(setdiff(all.vars(formula), nlpars))
    y$covars <- structure(covars, rsv_intercept = TRUE)
    lformula <- c(covars, lapply(y$nlpars, "[[", "allvars"))
    y$allvars <- allvars_formula(lformula)
    environment(y$allvars) <- environment(formula)
  } else {
    stop2("No non-linear parameters specified.")
  }
  structure(y, class = "btnl")
}

parse_ad <- function(formula, family = NULL, check_response = TRUE) {
  # extract addition arguments out of formula
  # Args:
  #   see parse_bf
  # Returns:
  #   A list of formulas each containg a single addition term
  x <- list()
  ad_funs <- lsp("brms", what = "exports", pattern = "^resp_")
  ad_funs <- sub("^resp_", "", ad_funs)
  families <- family_names(family)
  if (is.family(family) && any(nzchar(families))) {
    ad <- get_matches("\\|[^~]*~", formula2str(formula))
    if (length(ad)) {
      # replace deprecated '|' by '+'
      ad <- paste("~", rename(substr(ad, 2, nchar(ad) - 1), "|", "+"))
      ad_terms <- attr(terms(formula(ad)), "term.labels")
      for (a in ad_funs) {
        matches <- grep(paste0("^(resp_)?", a, "\\(.+\\)$"), ad_terms)
        if (length(matches) == 1L) {
          x[[a]] <- ad_terms[matches]
          if (!grepl("^resp_", x[[a]])) {
            x[[a]] <- paste0("resp_", x[[a]])
          }
          ad_terms <- ad_terms[-matches]
          ad_fams <- ad_families(a)
          valid <- ad_fams[1] == "all" || all(families %in% ad_fams)
          if (!is.na(x[[a]]) && valid) {
            x[[a]] <- str2formula(x[[a]])
          } else {
            stop2("Argument '", a, "' is not supported for ", 
                  "family '", summary(family), "'.")
          } 
        } else if (length(matches) > 1L) {
          stop2("Each addition argument may only be defined once.")
        } 
      }
      if (length(ad_terms)) {
        stop2("The following addition terms are invalid:\n",
              collapse_comma(ad_terms))
      }
      if (is.formula(x$se) && is.formula(x$disp)) {
        stop2("Addition arguments 'se' and 'disp' cannot ", 
              "be used at the same time.")
      }
    }
    if (is_wiener(family) && check_response && !is.formula(x$dec)) {
      stop2("Addition argument 'dec' is required for family 'wiener'.")
    }
    if (is.mixfamily(family) && (is.formula(x$cens) || is.formula(x$trunc))) {
      stop2("Censoring or truncation is not yet allowed in mixture models.")
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
    mo_terms <- ulapply(mo_terms, eval2)
    mo_terms <- str2formula(mo_terms)
    if (!length(all.vars(mo_terms))) {
      stop2("No variable supplied to function 'mo'.")
    }
    attr(mo_terms, "rsv_intercept") <- TRUE
  }
  structure(mo_terms, pos = pos_mo_terms)
}

parse_cs <- function(formula) {
  # category specific terms for ordinal models
  all_terms <- all_terms(formula)
  pos_cs_terms <- grepl("^cse?\\([^\\|]+$", all_terms)
  cs_terms <- all_terms[pos_cs_terms]
  if (length(cs_terms)) {
    cs_terms <- ulapply(cs_terms, eval2)
    cs_terms <- str2formula(cs_terms)
    # do not test whether variables were supplied to 'cs'
    # to allow category specific group-level intercepts
    attr(cs_terms, "rsv_intercept") <- TRUE
  }
  structure(cs_terms, pos = pos_cs_terms)
}

parse_sm <- function(formula) {
  # parse smooth functions
  all_terms <- all_terms(formula)
  pos_sm_terms <- grepl("^(s|t2|te|ti)\\(", all_terms)
  sm_terms <- all_terms[pos_sm_terms]
  if (length(sm_terms)) {
    if (any(grepl("^(te|ti)\\(", sm_terms))) {
      stop2("Tensor product smooths 'te' and 'ti' are not yet ", 
            "implemented in brms. Consider using 't2' instead.")
    }
    covars <- byvars <- named_list(sm_terms)
    for (i in seq_along(sm_terms)) {
      es <- eval_smooth(sm_terms[i])
      covars[[i]] <- es$term
      if (es$by != "NA") {
        byvars[[i]] <- es$by 
      }
    }
    sm_terms <- str2formula(sm_terms)
    allvars <- mgcv::interpret.gam(sm_terms)$fake.formula
  } else {
    covars <- byvars <- NULL
    allvars <- ~ 1
  }
  structure_not_null(
    sm_terms, pos = pos_sm_terms, covars = covars, 
    byvars = byvars, allvars = allvars
  )
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
    me_terms <- str2formula(me_terms)
    if (!length(all.vars(me_terms))) {
      stop2("No variable supplied to function 'me'.")
    }
    attr(me_terms, "rsv_intercept") <- TRUE
  }
  structure(me_terms, pos = pos_me_terms)
}

parse_gp <- function(formula) {
  # extract terms for gaussian processes
  # Args:
  #   formula: a formula object
  all_terms <- all_terms(formula)
  pos_gp_terms <- grepl_expr("^gp\\([^:]*\\)$", all_terms)
  gp_terms <- all_terms[pos_gp_terms]
  if (length(gp_terms)) {
    eterms <- lapply(gp_terms, eval2)
    covars <- lapply(eterms, "[[", "term")
    byvars <- lapply(eterms, "[[", "by")
    allvars <- str2formula(unlist(c(covars, byvars)))
    allvars <- str2formula(all.vars(allvars))
    if (!length(all.vars(allvars))) {
      stop2("No variable supplied to function 'gp'.")
    }
    gp_terms <- str2formula(gp_terms)
  } else {
    byvars <- NULL
    allvars <- ~ 1
  }
  structure_not_null(
    gp_terms, pos = pos_gp_terms, 
    byvars = byvars, allvars = allvars
  )
}

parse_offset <- function(formula) {
  # extract offset terms from a formula
  terms <- terms(as.formula(formula))
  pos_offset_terms <- attr(terms, "offset")
  if (!is.null(pos_offset_terms)) {
    vars <- attr(terms, "variables")
    offset_terms <- ulapply(
      pos_offset_terms, function(i) deparse(vars[[i+1]])
    )
    offset_terms <- str2formula(offset_terms)
  } else {
    offset_terms <- character(0)
  }
  offset_terms
}

parse_re <- function(re_terms) {
  # generate a data.frame with all information about the group-level terms
  # Args:
  #   re_terms: A vector of group-level terms in extended lme4 syntax
  re_terms <- split_re_terms(re_terms)
  lhs_terms <- lhs_terms(re_terms)
  mid_terms <- mid_terms(re_terms)
  rhs_terms <- rhs_terms(re_terms)
  out <- vector("list", length(re_terms))
  type <- attr(re_terms, "type")
  for (i in seq_along(re_terms)) {
    id <- gsub("\\|", "", mid_terms[i])
    if (!nzchar(id)) id <- NA
    gcall <- eval2(rhs_terms[i])
    group <- paste0(gcall$type, collapse(gcall$groups))
    out[[i]] <- data.frame(
      group = group, gtype = gcall$type,
      gn = i, id = id, type = type[i],
      cor = substr(mid_terms[i], 1, 2) != "||",
      stringsAsFactors = FALSE
    )
    out[[i]]$gcall <- list(gcall)
    out[[i]]$form <- list(formula(paste("~", lhs_terms[i])))
  }
  if (length(out)) {
    out <- do.call(rbind, out)
    out <- out[order(out$group), ]
  } else {
    out <- data.frame(
      group = character(0), gn = numeric(0),
      id = numeric(0), cor = logical(0), 
      type = character(0), form = character(0)
    )
  }
  out
}

parse_resp <- function(formula, keep_dot_usc = FALSE) {
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
    out <- all_vars[1]
  } else if (length(dim(pseudo_resp)) == 2L) {
    # response is a matrix
    out <- colnames(pseudo_resp)
    empty_names <- which(!nchar(out))
    if (length(empty_names)) {
      out[empty_names] <- paste0("response", empty_names)
    }
  } else {
    stop2("Response part of 'formula' is invalid.")
  }
  out <- make.names(out, unique = TRUE)
  if (!keep_dot_usc) {
    out <- gsub("\\.|_", "", out)
  }
  out
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
  x$allvars <- str2formula(c(time, all.vars(group)))
  x
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

is.btl <- function(x) {
  inherits(x, "btl")
}

is.btnl <- function(x) {
  inherits(x, "btnl")
}

avoid_auxpars <- function(names, bterms) {
  # avoid ambiguous parameter names
  # Args:
  #   names: names to check for ambiguity
  #   bterms: object of class brmsterms
  #stopifnot(is.brmsterms(bterms))
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

check_nlpar <- function(nlpar) {
  stopifnot(length(nlpar) == 1L)
  nlpar <- as.character(nlpar)
  ifelse(nlpar == "mu", "", nlpar)
}

check_fauxpars <- function(x) {
  # check validity of fixed auxiliary parameters
  stopifnot(is.null(x) || is.list(x))
  pos_pars <- c(
    "sigma", "shape", "nu", "phi", "kappa", 
    "beta", "disc", "bs", "ndt", "theta"
  )
  prob_pars <- c("zi", "hu", "bias", "quantile")
  for (ap in names(x)) {
    apc <- auxpar_class(ap)
    if (apc %in% pos_pars && x[[ap]] < 0) {
      stop2("Parameter '", ap, "' must be positive.")
    }
    if (apc %in% prob_pars && (x[[ap]] < 0 || x[[ap]] > 1)) {
      stop2("Parameter '", ap, "' must be between 0 and 1.")
    }
  }
  invisible(TRUE)
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
    new <- parse_bf(re_formula, check_response = FALSE)$auxpars$mu$re
    old <- parse_bf(old_re_formula, check_response = FALSE)$auxpars$mu$re
    if (nrow(new)) {
      new_terms <- lapply(new$form, terms)
      found <- rep(FALSE, nrow(new))
      for (i in 1:nrow(new)) {
        group <- new$group[[i]]
        old_terms <- lapply(old$form[old$group == group], terms)
        j <- 1
        while (!found[i] && j <= length(old_terms)) {
          found[i] <- isTRUE(
            all(attr(new_terms[[i]], "term.labels") %in% 
                attr(old_terms[[j]], "term.labels")) &&
            attr(new_terms[[i]], "intercept") <=
              attr(old_terms[[j]], "intercept")
          )
          j <- j + 1
        }
      }  
      new <- new[found, ]
      if (nrow(new)) {
        forms <- ulapply(new$form, formula2str, rm = 1)
        groups <- ulapply(new$gcall, "[[", "label")
        re_terms <- paste("(", forms, "|", groups, ")")
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
    return(new_formula)
  }
  
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

ad_families <- function(x) {
  # names of valid families for addition arguments
  switch(x, 
    weights = "all",
    se = c("gaussian", "student", "cauchy"),
    trials = c("binomial", "zero_inflated_binomial"),
    cat = c("cumulative", "cratio", "sratio", "acat"), 
    cens = c(
      "gaussian", "student", "cauchy", "lognormal",
      "inverse.gaussian", "binomial", "poisson", 
      "geometric", "negbinomial", "exponential", 
      "weibull", "gamma", "exgaussian", "frechet",
      "asym_laplace", "gen_extreme_value"
    ),
    trunc = c(
      "gaussian", "student", "cauchy", "lognormal", 
      "binomial", "poisson", "geometric", "negbinomial",
      "exponential", "weibull", "gamma", "inverse.gaussian",
      "exgaussian", "frechet", "asym_laplace",
      "gen_extreme_value"
    ),
    disp = c(
      "gaussian", "student", "cauchy", "lognormal", 
      "gamma", "weibull", "negbinomial", "exgaussian",
      "asym_laplace"
    ),
    dec = c("wiener"),
    stop2("Addition argument '", x, "' is not supported.")
  )
}

allvars_formula <- function(lformula) {
  # combine all variables in one formuula
  # Args:
  #   lformula: list of formulas or character strings
  out <- collapse(ulapply(rmNULL(lformula), plus_rhs))
  out <- str2formula(c(out, all.vars(parse(text = out))))
  update(out, ~ .)
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

#' @export
get_re.brmsterms <- function(x, all = TRUE, ...) {
  # get group-level information in a data.frame
  # Args:
  #   bterms: object of class brmsterms
  #   all: logical; include ranefs of nl and aux parameters?
  old_mv <- isTRUE(attr(x$formula, "old_mv"))
  if (all) {
    re <- named_list(names(x$auxpars))
    for (ap in names(re)) {
      re[[ap]] <- get_re(
        x$auxpars[[ap]], nlpar = ap,
        response = x$response, old_mv = old_mv
      )
    }
    re <- do.call(rbind, re)
  } else {
    x$auxpars[["mu"]]$nlpars <- NULL
    re <- get_re(x$auxpars[["mu"]])
  }
  re
}

#' @export
get_re.btl <- function(x, nlpar = "", response = "", old_mv = FALSE, ...) {
  stopifnot(is.data.frame(x$re))
  nlpar <- check_nlpar(nlpar)
  re <- x$re
  nresp <- length(response)
  if (!old_mv && nresp > 1L && nrow(re)) {
    # new MV models are also using the 'nlpar' argument
    re <- replicate(nresp, re, simplify = FALSE)
    for (i in seq_len(nresp)) {
      re[[i]]$nlpar <- rep(response[i], nrow(re[[i]]))
    }
    re <- do.call(rbind, re)
  } else {
    re$nlpar <- rep(nlpar, nrow(re)) 
  }
  re
}

#' @export
get_re.btnl <- function(x, ...) {
  re <- named_list(names(x$nlpars))
  for (nlp in names(re)) {
    re[[nlp]] <- get_re(x$nlpars[[nlp]], nlpar = nlp)
  }
  do.call(rbind, re)
}

#' @export
get_effect.brmsterms <- function(x, target = "fe", all = TRUE, ...) {
  # get formulas of certain effects in a list
  # Args:
  #   target: type of effects to return
  #   all: logical; include effects of nlpars and auxpars?
  if (all) {
    out <- named_list(names(x$auxpars))
    for (ap in names(out)) {
      out[[ap]] <- get_effect(x$auxpars[[ap]], target = target)
    }
  } else {
    x$auxpars[["mu"]]$nlpars <- NULL
    out <- get_effect(x$auxpars[["mu"]], target = target)
  }
  unlist(out, recursive = FALSE)
}

#' @export
get_effect.btl <- function(x, target = "fe", ...) {
  x[[target]]
}

#' @export
get_effect.btnl <- function(x, target = "fe", ...) {
  out <- named_list(names(x$nlpars))
  for (nlp in names(out)) {
    out[[nlp]] <- get_effect(x$nlpars[[nlp]], target = target)
  }
  rmNULL(out)
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

#' @export
get_all_effects.brmsterms <- function(x, rsv_vars = NULL, 
                                      comb_all = FALSE) {
  # get all effects for use in marginal_effects
  # Args:
  #   bterms: object of class brmsterms
  #   rsv_vars: character vector of reserved variables
  #   comb_all: include all main effects and two-way interactions?
  # Returns:
  #   a list with one element per valid effect / effects combination
  #   excludes all 3-way or higher interactions
  #stopifnot(is.brmsterms(bterms))
  stopifnot(is.atomic(rsv_vars))
  out <- list()
  for (ap in names(x$auxpars)) {
    out <- c(out, get_all_effects(x$auxpars[[ap]]))
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

#' @export
get_all_effects.btl <- function(x, ...) {
  int_formula <- function(x) {
    formula(paste("~", paste(x, collapse = "*")))
  }
  covars <- attr(x$sm, "covars")
  byvars <- attr(x$sm, "byvars")
  svars <- mapply(c, covars, byvars, SIMPLIFY = FALSE)
  alist <- lapply(svars, int_formula)
  get_var_combs(x$fe, x$mo, x$cs, x$me, x$gp, alist = alist)
}

#' @export
get_all_effects.btnl <- function(x, ...) {
  covars <- all.vars(rhs(x$covars))
  covars_comb <- as.list(covars)
  if (length(covars) > 1L) {
    covars_comb <- c(covars_comb, 
      utils::combn(covars, 2, simplify = FALSE)
    )
  }
  nl_effects <- lapply(x$nlpars, get_all_effects)
  nl_effects <- unlist(nl_effects, recursive = FALSE)
  unique(c(covars_comb, nl_effects))
}

get_sm_labels <- function(x, data = NULL, covars = FALSE,
                          combine = TRUE) {
  # extract labels of smooth terms
  # Args:
  #   x: either a formula or a list containing an element "sm"
  #   data: optional data frame containing the covariates
  #   covars: should the names of the covariates be returned
  #           instead of the full term names?
  #   combine: combine names of the covariates (TRUE) 
  #            or just return the covariate names (FALSE)?
  if (is.formula(x)) {
    x <- parse_bf(x, check_response = FALSE)
    sm_form <- x$auxpars$mu[["sm"]]
  } else {
    sm_form <- x[["sm"]] 
  }
  if (!is.formula(sm_form)) {
    return(character(0))
  }
  return_covars <- covars
  term_labels <- rename(attr(terms(sm_form), "term.labels"), " ", "")
  sms <- term_labels[grepl("^(s|t2|te|ti)\\(", term_labels)] 
  byvars <- attr(sm_form, "byvars")
  if (return_covars) {
    sfuns <- get_matches("^[^\\(]+", sms)
    covars <- attr(sm_form, "covars")
    for (i in seq_along(covars)) {
      covars[[i]] <- c(covars[[i]], byvars[[i]])
    }
    if (combine) {
      sms <- paste0(sfuns, rename(ulapply(covars, collapse)))
    } else {
      sms <- covars
    }
  }
  if (length(sms) && !is.null(data)) {
    # one smooth term may contain multiple design matrices
    if (return_covars && !combine) {
      stop("Invalid combination of arguments. Please report a bug.")
    }
    sdata_fe <- data_fe(x, data, knots = attr(data, "knots"))
    by_levels <- attr(sdata_fe[["X"]], "by_levels")
    nby <- lengths(by_levels)
    sms <- as.list(sms)
    for (i in seq_along(sms)) {
      if (nby[i] > 0L) {
        sms[[i]] <- paste0(sms[[i]], rename(by_levels[[i]]))
      }
    }
    nby[nby == 0L] <- 1L
    termnum <- rep(seq_along(sms), nby)
    sms <- unlist(sms)
    knots <- sdata_fe[grepl("^knots_", names(sdata_fe))]
    nbases <- setNames(ulapply(knots, length), sms)
    alist <- nlist(nbases, by_levels, termnum)
    attributes(sms)[names(alist)] <- alist
  }
  structure_not_null(sms, byvars = byvars)
}

get_me_labels <- function(x, data) {
  # get labels of measurement error terms
  # Args:
  #   x: either a formula or a list containing an element "me"
  #   data: data frame containing the noisy variables
  if (is.formula(x)) {
    x <- parse_bf(x, check_response = FALSE)
    me_form <- x$auxpars$mu[["me"]]
  } else {
    me_form <- x[["me"]]
  }
  if (!is.formula(me_form)) {
    return(character(0))
  }
  mm <- get_model_matrix(me_form, data, rename = FALSE)
  not_one <- apply(mm, 2, function(x) any(x != 1))
  uni_me <- get_matches_expr("^me\\([^:]*\\)$", colnames(mm))
  uni_me <- unique(gsub("[[:space:]]", "", uni_me))
  structure(colnames(mm), not_one = not_one, uni_me = uni_me)
}

get_gp_labels <- function(x, data = NULL, covars = FALSE) {
  # get labels of gaussian process terms
  # Args:
  #   x: either a formula or a list containing an element "gp"
  #   covars: should the conbined names of the covariates be 
  #           returned instead of the full term names?
  if (is.formula(x)) {
    x <- parse_bf(x, check_response = FALSE)
    gp_form <- x$auxpars$mu[["gp"]]
  } else {
    gp_form <- x[["gp"]]
  }
  if (!is.formula(gp_form)) {
    return(character(0))
  }
  byvars = attr(gp_form, "byvars")
  gp_terms <- all_terms(gp_form)
  by_levels <- named_list(gp_terms)
  out <- rep(NA, length(gp_terms))
  for (i in seq_along(gp_terms)) {
    gp <- eval2(gp_terms[i])
    if (covars) {
      out[i] <- paste0("gp", collapse(gp$term))
      if (gp$by != "NA") {
        out[i] <- paste0(out[i], gp$by)
      }
    } else {
      out[i] <- gp$label
    }
    if (!is.null(data)) {
      if (gp$by != "NA") {
        Cgp <- get(gp$by, data)
        if (!is.numeric(Cgp)) {
          by_levels[[i]] <- levels(factor(Cgp))
        }
      }
    }
  }
  if (covars) {
    out <- rename(out)
  }
  structure_not_null(out, byvars = byvars, by_levels = by_levels)
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
  formula <- try(as.formula(formula), silent = TRUE)
  if (is(formula, "try-error")) {
    out <- FALSE
  } else {
    try_terms <- try(terms(formula), silent = TRUE)
    if (is(try_terms, "try-error")) {
      out <- FALSE
    } else {
      has_intercept <- attr(try_terms, "intercept")
      out <- !has_intercept && "intercept" %in% all.vars(rhs(formula))
    }
  }
  out
}

has_smooths <- function(bterms) {
  # check if smooths are present in the model
  length(get_effect(bterms, target = "sm")) > 0L
}

has_cs <- function(bterms) {
  # check if category specific effects are present in the model
  length(get_effect(bterms, target = "cs")) > 0L ||
    any(get_re(bterms)$type %in% "cs")
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
  re <- get_re(bterms, all = all)
  ranef <- vector("list", nrow(re))
  used_ids <- new_ids <- NULL
  id_groups <- list()
  j <- 1
  for (i in seq_len(nrow(re))) {
    if (re$type[[i]] == "mo") {
      coef <- prepare_mo_vars(re$form[[i]], data, check = FALSE)
      coef <- colnames(coef)
    } else if (re$type[[i]] == "cs") {
      coef <- colnames(get_model_matrix(re$form[[i]], data = data))
      if (is.null(ncat)) {
        # try to infer ncat from the data
        Y <- as.numeric(model.response(data))
        ncat <- max(Y) - min(Y) + 1
      }
      indices <- paste0("[", seq_len(ncat - 1), "]")
      coef <- as.vector(t(outer(coef, indices, paste0)))
    } else if (re$type[[i]] == "me") {
      coef <- rename(get_me_labels(re$form[[i]], data))
    } else {
      coef <- colnames(get_model_matrix(re$form[[i]], data = data)) 
    }
    avoid_auxpars(coef, bterms = bterms)
    rdat <- data.frame(
      id = re$id[[i]],
      group = re$group[[i]],
      gn = re$gn[[i]],
      gtype = re$gtype[[i]],
      coef = coef, cn = NA,
      nlpar = re$nlpar[[i]],
      cor = re$cor[[i]],
      type = re$type[[i]],
      stringsAsFactors = FALSE
    )
    rdat$gcall <- replicate(nrow(rdat), re$gcall[i]) 
    rdat$form <- replicate(nrow(rdat), re$form[[i]])
    id <- re$id[[i]]
    if (is.na(id)) {
      rdat$id <- j
      j <- j + 1
    } else {
      if (id %in% used_ids) {
        k <- match(id, used_ids)
        rdat$id <- new_ids[k]
        new_id_groups <- c(re$group[[i]], re$gcall[[i]]$groups)
        if (!identical(new_id_groups, id_groups[[k]])) {
          stop2("Can only combine group-level terms of the ",
                "same grouping factors.")
        }
      } else {
        used_ids <- c(used_ids, id)
        k <- length(used_ids)
        rdat$id <- new_ids[k] <- j
        id_groups[[k]] <- c(re$group[[i]], re$gcall[[i]]$groups)
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
      un_re <- re[!duplicated(re$group), ]
      levels <- named_list(un_re$group)
      for (i in seq_along(levels)) {
        # combine levels of all grouping factors within one grouping term
        levels[[i]] <- ulapply(un_re$gcall[[i]]$groups, 
          function(g) levels(factor(get(g, data)))
        )
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
  data.frame(
    id = numeric(0), group = character(0), gn = numeric(0),
    coef = character(0), cn = numeric(0), nlpar = character(0),
    cor = logical(0), type = character(0), form = character(0), 
    stringsAsFactors = FALSE
  )
}

rsv_vars <- function(bterms, incl_intercept = TRUE) {
  # returns names of reserved variables
  # Args:
  #   bterms: object of class brmsterms
  #   incl_intercept: treat variable 'intercept' as reserved?
  stopifnot(is.brmsterms(bterms))
  nresp <- length(bterms$response)
  family <- bterms$family
  old_mv <- isTRUE(attr(bterms$formula, "old_mv"))
  rsv_intercept <- any(ulapply(bterms$auxpars, has_rsv_intercept))
  if (old_mv) {
    if (is_linear(family) && nresp > 1L || is_categorical(family)) {
      rsv <- c("trait", "response")
    } else if (is_forked(family)) {
      rsv <- c("trait", "response", "main", "spec")
    } else {
      rsv <- character(0)
    }
  } else {
    rsv <- character(0)
  }
  if (incl_intercept && rsv_intercept) {
    rsv <- c(rsv, "intercept")
  }
  rsv
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
  .exclude_pars <- function(bt, nlpar = "") {
    stopifnot(is.btl(bt))
    nlpar <- usc(check_nlpar(nlpar))
    out <- c(
      paste0("temp", nlpar, "_Intercept"),
      paste0(c("hs_local", "hs_global", "zb"), nlpar)
    )
    sms <- get_sm_labels(bt, data)
    if (length(sms) && !is.null(data)) {
      for (i in seq_along(sms)) {
        nb <- seq_len(attr(sms, "nbases")[[i]])
        out <- c(out, paste0("zs", nlpar, "_", i, "_", nb))
      } 
    }
    meef <- get_me_labels(bt, data)
    if (!save_mevars && length(meef)) {
      out <- c(out, paste0("Xme", nlpar, "_", seq_along(meef)))
    }
    return(out)
  }
  out <- c(
    "temp_Intercept1", "ordered_Intercept", "Rescor", "Lrescor", 
    "Sigma", "LSigma", "res_cov_matrix", "theta",
    intersect(auxpars(), names(bterms$auxpars))
  )
  if (length(bterms$response) > 1L) {
    for (r in bterms$response) {
      out <- c(out, .exclude_pars(bterms$auxpars$mu, nlpar = r))
    }
    bterms$auxpars$mu <- NULL
  }
  for (ap in names(bterms$auxpars)) {
    bt <- bterms$auxpars[[ap]]
    if (length(bt$nlpars)) {
      for (nlp in names(bt$nlpars)) {
        out <- c(out, .exclude_pars(bt$nlpars[[nlp]], nlpar = nlp))
      }
    } else {
      out <- c(out, .exclude_pars(bt, nlpar = ap))
    }
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
