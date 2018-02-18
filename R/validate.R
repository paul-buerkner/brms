#' Parse Formulas of \pkg{brms} Models
#' 
#' Parse formulas objects for use in \pkg{brms}.
#' 
#' @inheritParams brm
#' @param check_response Logical; Indicates whether the left-hand side 
#'   of \code{formula} (i.e. response variables and addition arguments) 
#'   should be parsed. If \code{FALSE}, \code{formula} may also be one-sided.
#' @param resp_rhs_all Logical; Indicates whether to also include response 
#'   variables on the right-hand side of formula \code{.$allvars}, 
#'   where \code{.} represents the output of \code{parse_bf}.
#' @param mv Indicates if the univariate model is part of a multivariate model.
#' @param rescor Indicates if residual correlations should be estimated.
#'   Only relevant in multivariate models.
#' @param ... Further arguments passed to or from other methods.
#'  
#' @return An object of class \code{brmsterms} or \code{mvbrmsterms} 
#'   (for multivariate models), which is a \code{list} containing all 
#'   required information initially stored in \code{formula} 
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
#'   \code{\link{brm}}, 
#'   \code{\link{brmsformula}},
#'   \code{\link{mvbrmsformula}}
#' 
#' @export
parse_bf <- function(formula, ...) {
  UseMethod("parse_bf")
}

#' @rdname parse_bf
#' @export
parse_bf.default <- function(formula, family = NULL, autocor = NULL, ...) {
  x <- validate_formula(formula, family = family, autocor = autocor)
  parse_bf(x, ...)
}

#' @rdname parse_bf
#' @export
parse_bf.brmsformula <- function(formula, family = NULL, autocor = NULL, 
                                 check_response = TRUE, resp_rhs_all = TRUE,
                                 mv = FALSE, rescor = FALSE, ...) {
  x <- validate_formula(formula, family = family, autocor = autocor)
  mv <- as_one_logical(mv)
  rescor <- as_one_logical(rescor)
  formula <- x$formula
  family <- x$family
  autocor <- x$autocor
  y <- nlist(formula, family, autocor, mv, rescor) 
  class(y) <- "brmsterms"
  
  if (check_response) {
    # extract response variables
    y$respform <- lhs(formula)
    if (is.null(y$respform)) {
      stop2("Invalid formula: response variable is missing.")
    }
    y$respform <- formula(gsub("\\|+[^~]*~", "~", formula2str(y$respform)))
    if (mv) {
      y$resp <- parse_resp(y$respform) 
    } else {
      y$resp <- ""
    }
  }
  
  # extract addition arguments
  adforms <- parse_ad(formula, family, check_response)
  advars <- str2formula(ulapply(adforms, all.vars))
  y$adforms[names(adforms)] <- adforms
  
  # copy stuff from the formula to parameter 'mu'
  str_rhs_form <- formula2str(rhs(formula))
  rhs_needed <- FALSE
  if (is.mixfamily(family)) {
    for (i in seq_along(family$mix)) {
      mui <- paste0("mu", i)
      if (!is.formula(x$pforms[[mui]])) {
        x$pforms[[mui]] <- eval2(paste0(mui, str_rhs_form))
        attr(x$pforms[[mui]], "nl") <- attr(formula, "nl")
        rhs_needed <- TRUE
      }
    }
  } else if (is_categorical(x$family)) {
    for (dp in x$family$dpars) {
      if (!is.formula(x$pforms[[dp]])) {
        x$pforms[[dp]] <- eval2(paste0(dp, str_rhs_form))
        attr(x$pforms[[dp]], "nl") <- attr(formula, "nl")
        rhs_needed <- TRUE
      }
    }
  } else {
    if (!is.formula(x$pforms[["mu"]])) { 
      x$pforms$mu <- eval2(paste0("mu", str_rhs_form))
      attr(x$pforms$mu, "nl") <- attr(formula, "nl")
      rhs_needed <- TRUE
    }
    x$pforms <- x$pforms[c("mu", setdiff(names(x$pforms), "mu"))]
  }
  terms <- try(terms(rhs(formula)), silent = TRUE)
  has_terms <- is(terms, "try-error") || 
    length(attr(terms, "term.labels")) ||
    length(attr(terms, "offset"))
  if (!rhs_needed && has_terms) {
    stop2("All 'mu' parameters are specified so that ",
          "the right-hand side of 'formula' is unused.")
  }
  
  # predicted distributional parameters
  resp <- ifelse(mv && !is.null(y$resp), y$resp, "")
  dpars <- is_dpar_name(names(x$pforms), family, bterms = y)
  dpars <- names(x$pforms)[dpars]
  dpar_forms <- x$pforms[dpars]
  nlpars <- setdiff(names(x$pforms), dpars)
  nlpar_forms <- x$pforms[nlpars]
  dpar_of_nlpars <- rep(NA, length(nlpars))
  for (i in seq_along(nlpar_forms)) {
    dpar <- attr(nlpar_forms[[i]], "dpar", TRUE)
    if (length(dpar) != 1L) {
      stop2("Parameter '", nlpars[i], "' is not part of the model")
    }
    dpar_of_nlpars[i] <- dpar
  }
  for (dp in dpars) {
    if (get_nl(dpar_forms[[dp]])) {
      if (is.mixfamily(family) || is_ordinal(family)) {
        stop2("Non-linear formulas are not yet allowed for this family.")
      }
      dpar_nlpar_forms <- nlpar_forms[dpar_of_nlpars %in% dp]
      y$dpars[[dp]] <- parse_nlf(dpar_forms[[dp]], dpar_nlpar_forms, dp, resp)
    } else {
      y$dpars[[dp]] <- parse_lf(dpar_forms[[dp]], family = family)
    }
    y$dpars[[dp]]$family <- dpar_family(family, dp)
    y$dpars[[dp]]$dpar <- dp
    y$dpars[[dp]]$resp <- resp
  }
  y <- store_uni_me(y)
  # fixed distributional parameters
  valid_dpars <- valid_dpars(y)
  inv_fixed_dpars <- setdiff(names(x$pfix), valid_dpars)
  if (length(inv_fixed_dpars)) {
    stop2("Invalid fixed parameters: ", collapse_comma(inv_fixed_dpars))
  }
  if ("sigma" %in% valid_dpars && no_sigma(y)) {
    # some models require setting sigma to 0
    if ("sigma" %in% c(names(x$pforms), names(x$pfix))) {
      stop2("Cannot predict or fix 'sigma' in this model.")
    }
    x$pfix$sigma <- 0
  }
  if ("disc" %in% valid_dpars) {
    # 'disc' is set to 1 and not estimated by default
    if (!"disc" %in% c(names(x$pforms), names(x$pfix))) {
      x$pfix$disc <- 1
    }
  }
  for (dp in names(x$pfix)) {
    y$fdpars[[dp]] <- list(value = x$pfix[[dp]], dpar = dp)
  }
  check_fdpars(y$fdpars)
  # check for illegal use of cs terms
  if (has_cs(y) && !(is.null(family) || allows_cs(family))) {
    stop2("Category specific effects are only meaningful for ", 
          "families 'sratio', 'cratio', and 'acat'.")
  }
  # parse autocor formula
  if (!is.null(y$dpars[["mu"]])) {
    y$dpars$mu$autocor <- autocor
  }
  y$time <- parse_time(autocor)
  
  # make a formula containing all required variables
  lhsvars <- if (resp_rhs_all) all.vars(y$respform)
  lformula <- c(
    lhsvars, advars, 
    lapply(y$dpars, "[[", "allvars"), 
    y$time$allvars
  )
  y$allvars <- allvars_formula(lformula)
  if (check_response) {
    y$allvars <- update(y$respform, y$allvars) 
  }
  environment(y$allvars) <- environment(formula)
  y
}

#' @rdname parse_bf
#' @export
parse_bf.mvbrmsformula <- function(formula, family = NULL, autocor = NULL, ...) {
  x <- validate_formula(formula, family = family, autocor = autocor)
  x$rescor <- isTRUE(x$rescor)
  out <- structure(list(), class = "mvbrmsterms")
  out$terms <- lapply(x$forms, parse_bf, mv = TRUE, rescor = x$rescor, ...)
  out$allvars <- allvars_formula(lapply(out$terms, "[[", "allvars"))
  # required to find variables used solely in the response part
  lhs_resp <- function(x) deparse_combine(lhs(x$respform)[[2]])
  out$respform <- paste0(ulapply(out$terms, lhs_resp), collapse = ",")
  out$respform <- formula(paste0("cbind(", out$respform, ") ~ 1"))
  out$responses <- ulapply(out$terms, "[[", "resp")
  out$rescor <- x$rescor
  out
}

parse_lf <- function(formula, family = NULL) {
  # parse linear formulas
  # Args:
  #   formula: an ordinary R formula
  #   family: the model family
  # Returns:
  #   object of class 'btl'
  formula <- rhs(as.formula(formula))
  check_multiple_special_terms(formula)
  y <- nlist(formula)
  types <- c("re", "sp", "cs", "sm", "gp", "offset")
  for (t in types) {
    tmp <- do.call(paste0("parse_", t), list(formula))
    if (is.data.frame(tmp) || is.formula(tmp)) {
      y[[t]] <- tmp 
    }
  }
  pos_special <- Reduce("|", rmNULL(lapply(y[types], attr, "pos")))
  y$fe <- parse_fe(formula, pos_special)
  lformula <- c(
    y[c("fe", "sp", "cs")],
    attr(y$sm, "allvars"), attr(y$gp, "allvars"),
    y$re$form, lapply(y$re$gcall, "[[", "allvars"),
    str2formula(all.vars(y$offset))
  )
  y$allvars <- allvars_formula(lformula)
  environment(y$allvars) <- environment(formula)
  structure(y, class = "btl")
}

parse_nlf <- function(formula, nlpar_forms, dpar, resp = "") {
  # parse non-linear formulas
  # Args:
  #   formula: formula of the non-linear model
  #   nlpar_forms: a list for formulas specifying linear predictors 
  #                for non-linear parameters
  # Returns:
  #   object of class 'btnl'
  stopifnot(is.list(nlpar_forms))
  stopifnot(length(dpar) == 1L)
  formula <- rhs(as.formula(formula))
  y <- list()
  if (length(nlpar_forms)) {
    y$formula <- formula
    nlpars <- names(nlpar_forms)
    y$nlpars <- named_list(nlpars)
    for (nlp in nlpars) {
      y$nlpars[[nlp]] <- parse_lf(nlpar_forms[[nlp]])
      y$nlpars[[nlp]]$nlpar <- nlp
      y$nlpars[[nlp]]$dpar <- dpar
      y$nlpars[[nlp]]$resp <- resp
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
    str_formula <- formula2str(formula)
    ad <- get_matches("(?<=\\|)[^~]*(?=~)", str_formula, perl = TRUE)
    if (length(ad)) {
      ad_terms <- attr(terms(formula(paste("~", ad))), "term.labels")
      for (a in ad_funs) {
        matches <- grep(paste0("^(resp_)?", a, "\\(.*\\)$"), ad_terms)
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

parse_fe <- function(formula, pos_special = NULL) {
  # extract fixed effects terms
  # Args:
  #   pos_special: logical position vector of non-FE terms
  terms <- terms(formula)
  all_terms <- all_terms(terms)
  if (length(pos_special)) {
    stopifnot(length(pos_special) == length(all_terms))
    fe_terms <- all_terms[!pos_special]
  } else {
    fe_terms <- all_terms
  }
  int_term <- attr(terms, "intercept")
  fe_terms <- paste(c(int_term, fe_terms), collapse = "+")
  fe_form <- str2formula(fe_terms)
  if (has_rsv_intercept(fe_form)) {
    attr(fe_form, "rsv_intercept") <- TRUE
  }
  fe_form
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

parse_sp <- function(formula) {
  # extract special effects terms 
  # Args:
  #   formula: a formula object
  all_terms <- all_terms(formula)
  # do not include group-level terms
  all_terms[grepl("\\|", all_terms)] <- ""
  pos_terms <- grepl_expr(regex_sp(), all_terms)
  out <- all_terms[pos_terms]
  if (length(out)) {
    uni_mo <- rm_wsp(get_matches_expr(regex_sp("mo"), out))
    uni_me <- rm_wsp(get_matches_expr(regex_sp("me"), out))
    uni_mi <- rm_wsp(get_matches_expr(regex_sp("mi"), out))
    out <- str2formula(out)
    attr(out, "rsv_intercept") <- TRUE
    attr(out, "uni_mo") <- uni_mo
    attr(out, "uni_me") <- uni_me
    attr(out, "uni_mi") <- uni_mi
  }
  structure(out, pos = pos_terms)
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
      es <- eval2(sm_terms[i])
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

parse_re <- function(formula) {
  # generate a data.frame with all information about the group-level terms
  re_terms <- get_re_terms(formula, brackets = FALSE)
  re_pos <- attr(re_terms, "pos")
  re_terms <- split_re_terms(re_terms)
  re_parts <- re_parts(re_terms)
  out <- vector("list", length(re_terms))
  type <- attr(re_terms, "type")
  for (i in seq_along(re_terms)) {
    id <- gsub("\\|", "", re_parts$mid[i])
    if (!nzchar(id)) id <- NA
    gcall <- eval2(re_parts$rhs[i])
    group <- paste0(gcall$type, collapse(gcall$groups))
    out[[i]] <- data.frame(
      group = group, gtype = gcall$type, 
      gn = i, id = id, type = type[i],
      cor = substr(re_parts$mid[i], 1, 2) != "||",
      stringsAsFactors = FALSE
    )
    out[[i]]$gcall <- list(gcall)
    out[[i]]$form <- list(formula(paste("~", re_parts$lhs[i])))
  }
  if (length(out)) {
    out <- do.call(rbind, out)
    out <- out[order(out$group), ]
  } else {
    out <- data.frame(
      group = character(0), gtype = character(0),
      gn = numeric(0), id = numeric(0), type = character(0), 
      cor = logical(0), form = character(0)
    )
  }
  structure(out, pos = re_pos)
}

parse_resp <- function(formula, check_names = TRUE) {
  # extract response variable names
  # assumes multiple response variables to be combined via cbind
  formula <- lhs(as.formula(formula))
  if (is.null(formula)) {
    return(NULL)
  }
  str_formula <- gsub("\\|+[^~]*~", "~", formula2str(formula))
  expr <- formula(str_formula)[[2]]
  if (length(expr) <= 1L) {
    out <- deparse_no_string(expr)
  } else {
    str_fun <- deparse_no_string(expr[[1]]) 
    use_cbind <- identical(str_fun, "cbind")
    if (use_cbind) {
      out <- ulapply(expr[-1], deparse_no_string)
    } else {
      out <- deparse_no_string(expr) 
    }
  }
  if (check_names) {
    out <- gsub("\\.|_", "", make.names(out, unique = TRUE))
  }
  out
}

parse_time <- function(autocor) {
  # extract time and grouping variables for autocorrelation structures
  # Args:
  #   autocor: object of class 'cor_brms'
  # Returns: 
  #   a list with elements time, group, and all, where all contains a 
  #   formula with all variables in formula
  out <- list()
  formula <- autocor$formula
  if (is.null(formula)) {
    formula <- ~ 1 
  }
  if (!is.null(lhs(formula))) {
    stop2("Autocorrelation formulas must be one-sided.")
  }
  formula <- formula2str(formula)
  time <- as.formula(paste("~", gsub("~|\\|[[:print:]]*", "", formula)))
  time_vars <- all.vars(time)
  if (is.cor_car(autocor) && length(time_vars) > 0L) {
    stop2("The CAR structure should not contain a 'time' variable.")
  }
  if (length(time_vars) > 1L) {
    stop2("Autocorrelation structures may only contain 1 time variable.")
  }
  if (length(time_vars)) {
    out$time <- time_vars
  }
  group <- sub("^\\|*", "", sub("~[^\\|]*", "", formula))
  if (illegal_group_expr(group)) {
    stop2("Illegal grouping term: ", group, "\nIt may contain only ", 
          "variable names combined by the symbol ':'")
  }
  group <- formula(paste("~", ifelse(nchar(group), group, "1")))
  group_vars <- all.vars(group)
  if (length(group_vars)) {
    out$group <- paste0(group_vars, collapse = ":")
  }
  out$allvars <- str2formula(c(time_vars, group_vars))
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

#' Checks if argument is a \code{mvbrmsterms} object
#' 
#' @param x An \R object
#' 
#' @seealso \code{\link[brms:parse_bf]{parse_bf}}
#' 
#' @export
is.mvbrmsterms <- function(x) {
  inherits(x, "mvbrmsterms")
}

is.btl <- function(x) {
  inherits(x, "btl")
}

is.btnl <- function(x) {
  inherits(x, "btnl")
}

as.brmsterms <- function(x) {
  # transform mvbrmsterms objects for use in stan_llh.brmsterms
  stopifnot(is.mvbrmsterms(x), x$rescor)
  families <- ulapply(x$terms, function(y) y$family$family)
  stopifnot(all(families == families[1]))
  out <- structure(list(), class = "brmsterms")
  out$family <- structure(
    list(family = paste0(families[1], "_mv"), link = "identity"),
    class = c("brmsfamily", "family")
  )
  out$sigma_pred <- any(ulapply(x$terms, 
    function(x) "sigma" %in% names(x$dpar) || is.formula(x$adforms$se)
  ))
  weight_forms <- rmNULL(lapply(x$terms, function(x) x$adforms$weights))
  if (length(weight_forms)) {
    str_wf <- unique(ulapply(weight_forms, formula2str))
    if (length(str_wf) > 1L) {
      stop2("All responses should use the same", 
            "weights if 'rescor' is estimated.")
    }
    out$adforms$weights <- weight_forms[[1]]
  }
  miforms <- rmNULL(lapply(x$terms, function(x) x$adforms$mi))
  if (length(miforms)) {
    out$adforms$mi <- miforms[[1]]
  }
  out
}

avoid_dpars <- function(names, bterms) {
  # avoid ambiguous parameter names
  # Args:
  #   names: names to check for ambiguity
  #   bterms: object of class brmsterms
  dpars <- c(names(bterms$dpars), "sp", "cs")
  if (length(dpars)) {
    dpars_prefix <- paste0("^", dpars, "_")
    invalid <- any(ulapply(dpars_prefix, grepl, names))
    if (invalid) {
      dpars <- paste0("'", dpars, "_'", collapse = ", ")
      stop2("Variable names starting with ", dpars,
            " are not allowed for this model.")
    }
  }
  invisible(NULL)
}

vars_prefix <- function() {
  c("dpar", "resp", "nlpar") 
}

check_prefix <- function(x, keep_mu = FALSE) {
  vpx <- vars_prefix()
  if (is.data.frame(x) && nrow(x) == 0) {
    # avoids a bug in data.frames with zero rows
    x <- list()
  }
  x[setdiff(vpx, names(x))] <- ""
  x <- x[vpx]
  for (i in seq_along(x)) {
    x[[i]] <- as.character(x[[i]])
    if (!length(x[[i]])) {
      x[[i]] <- ""
    }
    x[[i]] <- ifelse(
      !keep_mu & names(x)[i] == "dpar" & x[[i]] %in% "mu", 
      yes = "", no = x[[i]]
    )
    x[[i]] <- ifelse(
      keep_mu & names(x)[i] == "dpar" & x[[i]] %in% "", 
      yes = "mu", no = x[[i]]
    )
  }
  x
}

combine_prefix <- function(prefix, keep_mu = FALSE) {
  prefix <- check_prefix(prefix, keep_mu = keep_mu)
  prefix <- lapply(prefix, function(x)
    ifelse(nzchar(x), paste0("_", x), "")
  )
  sub("^_", "", do.call(paste0, prefix))
}

check_fdpars <- function(x) {
  # check validity of fixed distributional parameters
  stopifnot(is.null(x) || is.list(x))
  pos_pars <- c(
    "sigma", "shape", "nu", "phi", "kappa", 
    "beta", "disc", "bs", "ndt", "theta"
  )
  prob_pars <- c("zi", "hu", "bias", "quantile")
  for (dp in names(x)) {
    apc <- dpar_class(dp)
    value <- x[[dp]]$value
    if (apc %in% pos_pars && value < 0) {
      stop2("Parameter '", dp, "' must be positive.")
    }
    if (apc %in% prob_pars && (value < 0 || value > 1)) {
      stop2("Parameter '", dp, "' must be between 0 and 1.")
    }
  }
  invisible(TRUE)
}

check_multiple_special_terms <- function(x) {
  # check if one formula term contains multiple special terms
  # Args:
  #   x: (coerced to) a character vector of formula terms
  # Returns:
  #   This function is called for its side effects (errors)
  if (is.formula(x)) {
    x <- all_terms(x)
  }
  x <- as.character(x)
  sterms <- c(
    "(mo((no)?|(notonic)?))|(me)|(mi)", 
    "cse?", "(s|(t2)|(te)|(ti))", "gp"
  )
  sterms <- paste0("^(", sterms, ")\\([^:]*\\)$") 
  smatches <- matrix(NA, nrow = length(x), ncol = length(sterms))
  for (i in seq_along(sterms)) {
    smatches[, i] <- grepl_expr(sterms[i], x)
  }
  invalid <- x[rowSums(smatches) > 1]
  if (length(invalid)) {
    stop2("Cannot use multiple special terms within one term.\n",
          "Occured for: ", collapse_comma(invalid))
  }
  invisible(TRUE)
}

ad_families <- function(x) {
  # names of valid families for addition arguments
  switch(x, 
    weights = "all",
    se = c("gaussian", "student", "skew_normal"),
    trials = c("binomial", "zero_inflated_binomial"),
    cat = c("cumulative", "cratio", "sratio", "acat"), 
    cens = c(
      "gaussian", "student", "lognormal", "skew_normal",
      "inverse.gaussian", "binomial", "poisson", 
      "geometric", "negbinomial", "exponential", "beta",
      "weibull", "gamma", "exgaussian", "frechet",
      "asym_laplace", "gen_extreme_value", "shifted_lognormal"
    ),
    trunc = c(
      "gaussian", "student", "lognormal", "skew_normal",
      "binomial", "poisson", "geometric", "negbinomial",
      "exponential", "weibull", "gamma", "inverse.gaussian",
      "exgaussian", "frechet", "asym_laplace", "beta",
      "gen_extreme_value", "shifted_lognormal"
    ),
    mi = c(
      "gaussian", "student", "lognormal", "skew_normal",
      "inverse.gaussian", "exponential", "weibull", 
      "gamma", "exgaussian", "frechet", "beta",
      "asym_laplace", "gen_extreme_value"
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

get_effect <- function(x, ...) {
  # extract various kind of effects
  UseMethod("get_effect")
}

#' @export
get_effect.mvbrmsterms <- function(x, ...) {
  unlist(lapply(x$terms, get_effect, ...), recursive = FALSE)
}

#' @export
get_effect.brmsterms <- function(x, target = "fe", all = TRUE, ...) {
  # get formulas of certain effects in a list
  # Args:
  #   target: type of effects to return
  #   all: logical; include effects of nlpars and dpars?
  if (all) {
    out <- named_list(names(x$dpars))
    for (dp in names(out)) {
      out[[dp]] <- get_effect(x$dpars[[dp]], target = target)
    }
  } else {
    x$dpars$mu$nlpars <- NULL
    out <- get_effect(x$dpars[["mu"]], target = target)
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

get_advars <- function(x, ...) {
  # extract variable names used in addition terms
  UseMethod("get_advars")
}

#' @export
get_advars.brmsterms <- function(x, ad, ...) {
  ad <- as_one_character(ad)
  all.vars(x$adforms[[ad]])
}

#' @export
get_advars.mvbrmsterms <- function(x, ad, ...) {
  unique(ulapply(x$terms, get_advars, ad = ad, ...))
}

get_sm_labels <- function(x, data = NULL, covars = FALSE, combine = TRUE) {
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
    sm_form <- x$dpars$mu[["sm"]]
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
    by_levels <- attr(sdata_fe$X, "by_levels")
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

get_gp_labels <- function(x, data = NULL, covars = FALSE) {
  # get labels of gaussian process terms
  # Args:
  #   x: either a formula or a list containing an element "gp"
  #   covars: should the conbined names of the covariates be 
  #           returned instead of the full term names?
  if (is.formula(x)) {
    x <- parse_bf(x, check_response = FALSE)
    gp_form <- x$dpars$mu[["gp"]]
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

all_terms <- function(x) {
  if (!length(x)) {
    return(character(0))
  }
  if (!inherits(x, "terms")) {
    x <- terms(as.formula(x))
  }
  rm_wsp(attr(x, "term.labels"))
}

validate_terms <- function(x) {
  # validate a terms object (or one that can be coerced to it)
  # to be used in get_model_matrix
  # Args:
  #   x: any R object; if not a formula or terms, NULL is returned
  # Returns:
  #   a (possibly amended) terms object or NULL
  rsv_intercept <- isTRUE(attr(x, "rsv_intercept"))
  if (is.formula(x) && !inherits(x, "terms")) {
    x <- terms(x)
  }
  if (!inherits(x, "terms")) {
    return(NULL)
  }
  if (rsv_intercept) {
    # allows to remove the intercept without causing cell mean coding
    attr(x, "intercept") <- 1
    attr(x, "rm_intercept") <- TRUE
  }
  x
}

has_intercept <- function(formula) {
  # checks if the formula contains an intercept
  # can handle non-linear formulas
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

rsv_vars <- function(bterms) {
  # returns names of reserved variables
  # Args:
  #   bterms: object of class brmsterms
  stopifnot(is.brmsterms(bterms) || is.mvbrmsterms(bterms))
  .rsv_vars <- function(x) {
    rsv_int <- any(ulapply(x$dpars, has_rsv_intercept))
    if (rsv_int) "intercept" else NULL
  }
  if (is.mvbrmsterms(bterms)) {
    out <- unique(ulapply(bterms$terms, .rsv_vars))
  } else {
    out <- .rsv_vars(bterms)
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

get_autocor_vars <- function(x, ...) {
  # extract variable names used in autocor structures
  UseMethod("get_autocor_vars")
}

#' @export
get_autocor_vars.cor_brms <- function(x, var = "time", incl_car = TRUE, ...) {
  if (incl_car || !is.cor_car(x)) parse_time(x)[[var]]
}

#' @export
get_autocor_vars.brmsterms <- function(x, var = "time", incl_car = TRUE, ...) {
  if (incl_car || !is.cor_car(x$autocor)) x$time[[var]]
}

#' @export
get_autocor_vars.brmsformula <- function(x, ...) {
  get_autocor_vars(x$autocor, ...)
}

#' @export
get_autocor_vars.mvbrmsformula <- function(x, ...) {
  unique(ulapply(x$forms, get_autocor_vars, ...))
}

#' @export
get_autocor_vars.mvbrmsterms <- function(x, ...) {
  unique(ulapply(x$terms, get_autocor_vars, ...))
}

#' @export
get_autocor_vars.brmsfit <- function(x, ...) {
  get_autocor_vars(x$formula, ...)
}

get_bounds <- function(bterms, data = NULL, incl_family = FALSE, 
                       stan = FALSE) {
  # extract truncation boundaries
  # Returns:
  #   incl_family: include the family in the derivation of bounds?
  #   stan: return bounds in form of Stan syntax?
  stopifnot(is.brmsterms(bterms))
  formula <- bterms$adforms$trunc
  if (is.formula(formula)) {
    term <- attr(terms(formula), "term.labels")
    stopifnot(length(term) == 1L && grepl("resp_trunc\\(", term))
    out <- eval_rhs(formula, data = data)
  } else {
    out <- resp_trunc()
  }
  if (incl_family) {
    family_bounds <- get_family_bounds(bterms)
    out$lb <- max(out$lb, family_bounds$lb)
    out$ub <- min(out$ub, family_bounds$ub)
  }
  if (stan) {
    if (any(out$lb > -Inf | out$ub < Inf)) {
      tmp <- c(
        if (out$lb > -Inf) paste0("lower=", out$lb),
        if (out$ub < Inf) paste0("upper=", out$ub)
      )
      out <- paste0("<", paste0(tmp, collapse = ","), ">")
    } else {
      out <- ""
    }
  }
  out
}

get_family_bounds <- function(bterms) {
  # get boundaries of response distribution
  stopifnot(is.brmsterms(bterms))
  family <- bterms$family$family
  if (is.null(family)) {
    return(list(lb = -Inf, ub = Inf))
  }
  resp <- usc(bterms$resp)
  pos_families <- c(
    "poisson", "negbinomial", "geometric", "gamma", "weibull", 
    "exponential", "lognormal", "frechet", "inverse.gaussian", 
    "hurdle_poisson", "hurdle_negbinomial", "hurdle_gamma",
    "hurdle_lognormal", "zero_inflated_poisson", 
    "zero_inflated_negbinomial"
  )
  beta_families <- c("beta", "zero_inflated_beta", "zero_one_inflated_beta")
  ordinal_families <- c("cumulative", "cratio", "sratio", "acat")
  if (family %in% pos_families) {
    out <- list(lb = 0, ub = Inf)
  } else if (family %in% c("bernoulli", beta_families)) {
    out <- list(lb = 0, ub = 1)
  } else if (family %in% c("categorical", ordinal_families)) {
    out <- list(lb = 1, ub = paste0("ncat", resp))
  } else if (family %in% c("binomial", "zero_inflated_binomial")) {
    out <- list(lb = 0, ub = paste0("trials", resp))
  } else if (family %in% "von_mises") {
    out <- list(lb = -pi, ub = pi)
  } else if (family %in% c("wiener", "shifted_lognormal")) {
    out <- list(lb = paste("min_Y", resp), ub = Inf)
  } else {
    out <- list(lb = -Inf, ub = Inf)
  }
  out
}

has_cens <- function(bterms, data = NULL) {
  # indicate if the model is (possibly interval) censored
  stopifnot(is.brmsterms(bterms))
  formula <- bterms$adforms$cens
  if (is.formula(formula)) {
    term <- attr(terms(formula), "term.labels")
    stopifnot(length(term) == 1L && grepl("resp_cens\\(", term))
    out <- eval_rhs(formula, data = data)
    out <- structure(TRUE, interval = !is.null(attr(out, "y2")))
  } else {
    out <- FALSE
  }
  out
}

get_sdy <- function(x, data = NULL) {
  stopifnot(is.brmsterms(x))
  miform <- x$adforms[["mi"]]
  if (is.formula(miform)) {
    sdy <- eval_rhs(miform, data = data)
  } else {
    sdy <- NULL
  }
  sdy
}
