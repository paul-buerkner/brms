#' Parse Formulas of \pkg{brms} Models
#'
#' Parse formulas objects for use in \pkg{brms}.
#'
#' @aliases parse_bf
#'
#' @inheritParams brm
#' @param check_response Logical; Indicates whether the left-hand side
#'   of \code{formula} (i.e. response variables and addition arguments)
#'   should be parsed. If \code{FALSE}, \code{formula} may also be one-sided.
#' @param resp_rhs_all Logical; Indicates whether to also include response
#'   variables on the right-hand side of formula \code{.$allvars},
#'   where \code{.} represents the output of \code{brmsterms}.
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
brmsterms <- function(formula, ...) {
  UseMethod("brmsterms")
}

# the name 'parse_bf' is deprecated as of brms 2.12.4
# remove it eventually in brms 3.0
#' @export
parse_bf <- function(x, ...) {
  warning2("Method 'parse_bf' is deprecated. Please use 'brmsterms' instead.")
  UseMethod("brmsterms")
}

#' @rdname brmsterms
#' @export
brmsterms.default <- function(formula, ...) {
  brmsterms(validate_formula(formula), ...)
}

#' @rdname brmsterms
#' @export
brmsterms.brmsformula <- function(formula, check_response = TRUE,
                                  resp_rhs_all = TRUE, ...) {
  x <- validate_formula(formula)
  mv <- isTRUE(x$mv)
  rescor <- mv && isTRUE(x$rescor)
  mecor <- isTRUE(x$mecor)
  formula <- x$formula
  family <- x$family
  y <- nlist(formula, family, mv, rescor, mecor)
  y$cov_ranef <- x$cov_ranef
  class(y) <- "brmsterms"

  y$resp <- ""
  if (check_response) {
    # extract response variables
    y$respform <- validate_resp_formula(formula, empty_ok = FALSE)
    if (mv) {
      y$resp <- terms_resp(y$respform)
    }
  }

  # extract addition arguments
  adforms <- terms_ad(formula, family, check_response)
  advars <- str2formula(ulapply(adforms, all_vars))
  y$adforms[names(adforms)] <- adforms

  # centering would lead to incorrect results for grouped threshold vectors
  # as each threshold vector only affects a subset of observations
  if (!is.null(get_ad_expr(y, "thres", "gr"))) {
    attr(formula, "center") <- FALSE
    dp_classes <- dpar_class(names(x$pforms))
    mu_names <- names(x$pforms)[dp_classes == "mu"]
    for (dp in mu_names) {
      attr(x$pforms[[dp]], "center") <- FALSE
    }
  }

  # combine the main formula with formulas for the 'mu' parameters
  if (is.mixfamily(family)) {
    mu_dpars <- paste0("mu", seq_along(family$mix))
    for (dp in mu_dpars) {
      x$pforms[[dp]] <- combine_formulas(formula, x$pforms[[dp]], dp)
    }
    x$pforms <- move2start(x$pforms, mu_dpars)
    for (i in seq_along(family$mix)) {
      # store the respective mixture index in each mixture component
      # this enables them to be easily passed along, e.g. in stan_log_lik
      y$family$mix[[i]]$mix <- i
    }
  } else if (conv_cats_dpars(x$family)) {
    mu_dpars <- str_subset(x$family$dpars, "^mu")
    for (dp in mu_dpars) {
      x$pforms[[dp]] <- combine_formulas(formula, x$pforms[[dp]], dp)
    }
    x$pforms <- move2start(x$pforms, mu_dpars)
  } else {
    x$pforms[["mu"]] <- combine_formulas(formula, x$pforms[["mu"]], "mu")
    x$pforms <- move2start(x$pforms, "mu")
  }

  # predicted distributional parameters
  dpars <- intersect(names(x$pforms), valid_dpars(family))
  dpar_forms <- x$pforms[dpars]
  nlpars <- setdiff(names(x$pforms), dpars)

  y$dpars <- named_list(dpars)
  for (dp in dpars) {
    if (get_nl(dpar_forms[[dp]])) {
      y$dpars[[dp]] <- terms_nlf(dpar_forms[[dp]], nlpars, y$resp)
    } else {
      y$dpars[[dp]] <- terms_lf(dpar_forms[[dp]])
    }
    y$dpars[[dp]]$family <- dpar_family(family, dp)
    y$dpars[[dp]]$dpar <- dp
    y$dpars[[dp]]$resp <- y$resp
    if (dpar_class(dp) == "mu") {
      y$dpars[[dp]]$respform <- y$respform
      y$dpars[[dp]]$adforms <- y$adforms
    }
    y$dpars[[dp]]$transform <- stan_eta_transform(y, y$dpars[[dp]]$family)
    check_cs(y$dpars[[dp]])
  }

  y$nlpars <- named_list(nlpars)
  if (length(nlpars)) {
    nlpar_forms <- x$pforms[nlpars]
    for (nlp in nlpars) {
      if (is.null(attr(nlpar_forms[[nlp]], "center"))) {
        # design matrices of non-linear parameters will not be
        # centered by default to make prior specification easier
        attr(nlpar_forms[[nlp]], "center") <- FALSE
      }
      if (get_nl(nlpar_forms[[nlp]])) {
        y$nlpars[[nlp]] <- terms_nlf(nlpar_forms[[nlp]], nlpars, y$resp)
      } else {
        y$nlpars[[nlp]] <- terms_lf(nlpar_forms[[nlp]])
      }
      y$nlpars[[nlp]]$nlpar <- nlp
      y$nlpars[[nlp]]$resp <- y$resp
      check_cs(y$nlpars[[nlp]])
    }
    used_nlpars <- ufrom_list(c(y$dpars, y$nlpars), "used_nlpars")
    unused_nlpars <- setdiff(nlpars, used_nlpars)
    if (length(unused_nlpars)) {
      stop2(
        "The parameter '", unused_nlpars[1], "' is not a ",
        "valid distributional or non-linear parameter. ",
        "Did you forget to set 'nl = TRUE'?"
      )
    }
    # sort non-linear parameters after dependency
    used_nlpars <- from_list(y$nlpars, "used_nlpars")
    sorted_nlpars <- sort_dependencies(used_nlpars)
    y$nlpars <- y$nlpars[sorted_nlpars]
  }

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
  if ("nu" %in% valid_dpars && no_nu(y)) {
    if ("nu" %in% c(names(x$pforms), names(x$pfix))) {
      stop2("Cannot predict or fix 'nu' in this model.")
    }
    x$pfix$nu <- 1
  }
  disc_pars <- valid_dpars[dpar_class(valid_dpars) %in% "disc"]
  for (dp in disc_pars) {
    # 'disc' is set to 1 and not estimated by default
    if (!dp %in% c(names(x$pforms), names(x$pfix))) {
      x$pfix[[dp]] <- 1
    }
  }
  for (dp in names(x$pfix)) {
    y$fdpars[[dp]] <- list(value = x$pfix[[dp]], dpar = dp)
  }
  check_fdpars(y$fdpars)

  # make a formula containing all required variables
  y$unused <- attr(x$formula, "unused")
  lhsvars <- if (resp_rhs_all) all_vars(y$respform)
  y$allvars <- allvars_formula(
    lhsvars, advars, lapply(y$dpars, get_allvars),
    lapply(y$nlpars, get_allvars), y$time$allvars,
    get_unused_arg_vars(y),
    .env = environment(formula)
  )
  if (check_response) {
    # add y$respform to the left-hand side of y$allvars
    # avoid using update.formula as it is inefficient for longer formulas
    formula_allvars <- y$respform
    formula_allvars[[3]] <- y$allvars[[2]]
    environment(formula_allvars) <- environment(y$allvars)
    y$allvars <- formula_allvars
  }
  y
}

#' @rdname brmsterms
#' @export
brmsterms.mvbrmsformula <- function(formula, ...) {
  x <- validate_formula(formula)
  x$rescor <- isTRUE(x$rescor)
  x$mecor <- isTRUE(x$mecor)
  out <- structure(list(), class = "mvbrmsterms")
  out$terms <- named_list(names(x$forms))
  for (i in seq_along(out$terms)) {
    x$forms[[i]]$rescor <- x$rescor
    x$forms[[i]]$mecor <- x$mecor
    x$forms[[i]]$mv <- TRUE
    out$terms[[i]] <- brmsterms(x$forms[[i]], ...)
  }
  list_allvars <- lapply(out$terms, get_allvars)
  out$allvars <- allvars_formula(
    list_allvars, .env = environment(list_allvars[[1]])
  )
  # required to find variables used solely in the response part
  lhs_resp <- function(x) deparse0(lhs(x$respform)[[2]])
  out$respform <- paste0(ulapply(out$terms, lhs_resp), collapse = ",")
  out$respform <- formula(paste0("mvbind(", out$respform, ") ~ 1"))
  out$responses <- ufrom_list(out$terms, "resp")
  out$rescor <- x$rescor
  out$mecor <- x$mecor
  out$cov_ranef <- x$cov_ranef
  out
}

# parse linear/additive formulas
# @param formula an ordinary model formula
# @return a 'btl' object
terms_lf <- function(formula) {
  formula <- rhs(as.formula(formula))
  y <- nlist(formula)
  formula <- terms(formula)
  check_accidental_helper_functions(formula)
  types <- setdiff(all_term_types(), excluded_term_types(formula))
  for (t in types) {
    tmp <- do_call(paste0("terms_", t), list(formula))
    if (is.data.frame(tmp) || is.formula(tmp)) {
      y[[t]] <- tmp
    }
  }
  y$allvars <- allvars_formula(
    get_allvars(y$fe), get_allvars(y$re),
    get_allvars(y$cs), get_allvars(y$sp),
    get_allvars(y$sm), get_allvars(y$gp),
    get_allvars(y$ac), get_allvars(y$offset)
  )
  structure(y, class = "btl")
}

# parse non-linear formulas
# @param formula non-linear model formula
# @param nlpars names of all non-linear parameters
# @param resp optional name of a response variable
# @return a 'btnl' object
terms_nlf <- function(formula, nlpars, resp = "") {
  if (!length(nlpars)) {
    stop2("No non-linear parameters specified.")
  }
  loop <- !isFALSE(attr(formula, "loop"))
  formula <- rhs(as.formula(formula))
  y <- nlist(formula)
  all_vars <- all_vars(formula)
  y$used_nlpars <- intersect(all_vars, nlpars)
  covars <- setdiff(all_vars, nlpars)
  y$covars <- structure(str2formula(covars), int = FALSE)
  if (!"ac" %in% excluded_term_types(formula)) {
    y$ac <- terms_ac(attr(formula, "autocor"))
  }
  y$allvars <- allvars_formula(covars, get_allvars(y$ac))
  y$loop <- loop
  structure(y, class = "btnl")
}

# extract addition arguments out of formula
# @return a list of formulas each containg a single addition term
terms_ad <- function(formula, family = NULL, check_response = TRUE) {
  x <- list()
  ad_funs <- lsp("brms", what = "exports", pattern = "^resp_")
  ad_funs <- sub("^resp_", "", ad_funs)
  families <- family_names(family)
  if (is.family(family) && any(nzchar(families))) {
    str_formula <- formula2str(formula)
    ad <- get_matches("(?<=\\|)[^~]*(?=~)", str_formula, perl = TRUE)
    valid_ads <- family_info(family, "ad")
    if (length(ad)) {
      ad_terms <- terms(str2formula(ad))
      if (length(attr(ad_terms, "offset"))) {
        stop2("Offsets are not allowed in addition terms.")
      }
      ad_terms <- attr(ad_terms, "term.labels")
      for (a in ad_funs) {
        matches <- grep(paste0("^(resp_)?", a, "\\(.*\\)$"), ad_terms)
        if (length(matches) == 1L) {
          x[[a]] <- ad_terms[matches]
          if (!grepl("^resp_", x[[a]])) {
            x[[a]] <- paste0("resp_", x[[a]])
          }
          ad_terms <- ad_terms[-matches]
          if (!is.na(x[[a]]) && a %in% valid_ads) {
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
    if (check_response && "wiener" %in% families && !is.formula(x$dec)) {
      stop2("Addition argument 'dec' is required for family 'wiener'.")
    }
    if (is.formula(x$cat)) {
      # 'cat' was replaced by 'thres' in brms 2.10.5
      x$thres <- x$cat
    }
  }
  x
}

# extract fixed effects terms
terms_fe <- function(formula) {
  if (!is.terms(formula)) {
    formula <- terms(formula)
  }
  all_terms <- all_terms(formula)
  sp_terms <- find_terms(all_terms, "all", complete = FALSE)
  re_terms <- all_terms[grepl("\\|", all_terms)]
  int_term <- attr(formula, "intercept")
  fe_terms <- setdiff(all_terms, c(sp_terms, re_terms))
  out <- paste(c(int_term, fe_terms), collapse = "+")
  out <- str2formula(out)
  attr(out, "allvars") <- allvars_formula(out)
  attr(out, "decomp") <- get_decomp(formula)
  if (has_rsv_intercept(out, has_intercept(formula))) {
    attr(out, "int") <- FALSE
  }
  if (no_cmc(formula)) {
    attr(out, "cmc") <- FALSE
  }
  if (no_center(formula)) {
    attr(out, "center") <- FALSE
  }
  if (is_sparse(formula)) {
    attr(out, "sparse") <- TRUE
  }
  out
}

# gather information of group-level terms
# @return a data.frame with one row per group-level term
terms_re <- function(formula) {
  re_terms <- get_re_terms(formula, brackets = FALSE)
  if (!length(re_terms)) {
    return(NULL)
  }
  re_terms <- split_re_terms(re_terms)
  re_parts <- re_parts(re_terms)
  out <- allvars <- vector("list", length(re_terms))
  type <- attr(re_terms, "type")
  for (i in seq_along(re_terms)) {
    gcall <- eval2(re_parts$rhs[i])
    form <- str2formula(re_parts$lhs[i])
    group <- paste0(gcall$type, collapse(gcall$groups))
    out[[i]] <- data.frame(
      group = group, gtype = gcall$type, gn = i,
      id = gcall$id, type = type[i], cor = gcall$cor,
      stringsAsFactors = FALSE
    )
    out[[i]]$gcall <- list(gcall)
    out[[i]]$form <- list(form)
    # gather all variables used in the group-level term
    # at this point 'cs' terms are no longer recognized as such
    ftype <- str_if(type[i] %in% "cs", "", type[i])
    re_allvars <- get_allvars(form, type = ftype)
    allvars[[i]] <- allvars_formula(re_allvars, gcall$allvars)
  }
  out <- do_call(rbind, out)
  out <- out[order(out$group), ]
  attr(out, "allvars") <- allvars_formula(allvars)
  if (no_cmc(formula)) {
    # disabling cell-mean coding in all group-level terms
    # has to come last to avoid removal of attributes
    for (i in seq_rows(out)) {
      attr(out$form[[i]], "cmc") <- FALSE
    }
  }
  out
}

# extract category specific terms for ordinal models
terms_cs <- function(formula) {
  out <- find_terms(formula, "cs")
  if (!length(out)) {
    return(NULL)
  }
  out <- ulapply(out, eval2, envir = environment())
  out <- str2formula(out)
  attr(out, "allvars") <- allvars_formula(out)
  # do not test whether variables were supplied to 'cs'
  # to allow category specific group-level intercepts
  attr(out, "int") <- FALSE
  out
}

# extract special effects terms
terms_sp <- function(formula) {
  types <- c("mo", "me", "mi")
  out <- find_terms(formula, types, complete = FALSE)
  if (!length(out)) {
    return(NULL)
  }
  uni_mo <- get_matches_expr(regex_sp("mo"), out)
  uni_me <- get_matches_expr(regex_sp("me"), out)
  uni_mi <- get_matches_expr(regex_sp("mi"), out)
  # remove the intercept as it is handled separately
  out <- str2formula(c("0", out))
  attr(out, "int") <- FALSE
  attr(out, "uni_mo") <- uni_mo
  attr(out, "uni_me") <- uni_me
  attr(out, "uni_mi") <- uni_mi
  attr(out, "allvars") <- str2formula(all_vars(out))
  # TODO: do we need sp_fake_formula at all?
  # attr(out, "allvars") <- sp_fake_formula(uni_mo, uni_me, uni_mi)
  out
}

# extract spline terms
terms_sm <- function(formula) {
  out <- find_terms(formula, "sm")
  if (!length(out)) {
    return(NULL)
  }
  if (any(grepl("^(te|ti)\\(", out))) {
    stop2("Tensor product smooths 'te' and 'ti' are not yet ",
          "implemented in brms. Consider using 't2' instead.")
  }
  out <- str2formula(out)
  attr(out, "allvars") <- mgcv::interpret.gam(out)$fake.formula
  out
}

# extract gaussian process terms
terms_gp <- function(formula) {
  out <- find_terms(formula, "gp")
  if (!length(out)) {
    return(NULL)
  }
  eterms <- lapply(out, eval2, envir = environment())
  covars <- from_list(eterms, "term")
  byvars <- from_list(eterms, "by")
  allvars <- str2formula(unlist(c(covars, byvars)))
  allvars <- str2formula(all_vars(allvars))
  if (!length(all_vars(allvars))) {
    stop2("No variable supplied to function 'gp'.")
  }
  out <- str2formula(out)
  attr(out, "allvars") <- allvars
  out
}

# extract autocorrelation terms
terms_ac <- function(formula) {
  autocor <- attr(formula, "autocor")
  out <- c(find_terms(formula, "ac"), find_terms(autocor, "ac"))
  if (!length(out)) {
    return(NULL)
  }
  eterms <- lapply(out, eval2, envir = environment())
  allvars <- unlist(c(
    from_list(eterms, "time"),
    from_list(eterms, "gr")
  ))
  allvars <- str2formula(all_vars(allvars))
  out <- str2formula(out)
  attr(out, "allvars") <- allvars
  out
}

# extract offset terms
terms_offset <- function(formula) {
  if (!is.terms(formula)) {
    formula <- terms(as.formula(formula))
  }
  pos <- attr(formula, "offset")
  if (is.null(pos)) {
    return(NULL)
  }
  vars <- attr(formula, "variables")
  out <- ulapply(pos, function(i) deparse0(vars[[i + 1]]))
  out <- str2formula(out)
  attr(out, "allvars") <- str2formula(all_vars(out))
  out
}

# extract multiple covariates in multi-membership terms
terms_mmc <- function(formula) {
  out <- find_terms(formula, "mmc")
  if (!length(out)) {
    return(NULL)
  }
  # remove the intercept as it is handled separately
  out <- str2formula(c("0", out))
  attr(out, "allvars") <- allvars_formula(out)
  attr(out, "int") <- FALSE
  out
}

# extract response variable names
# assumes multiple response variables to be combined via mvbind
terms_resp <- function(formula, check_names = TRUE) {
  formula <- lhs(as.formula(formula))
  if (is.null(formula)) {
    return(NULL)
  }
  used_mvbind <- FALSE
  expr <- validate_resp_formula(formula)[[2]]
  if (length(expr) <= 1L) {
    out <- deparse_no_string(expr)
  } else {
    str_fun <- deparse_no_string(expr[[1]])
    used_mvbind <- grepl("^(brms:::?)?mvbind$", str_fun)
    if (used_mvbind) {
      out <- ulapply(expr[-1], deparse_no_string)
    } else {
      out <- deparse_no_string(expr)
    }
  }
  if (check_names) {
    out <- make_stan_names(out)
  }
  attr(out, "mvbind") <- used_mvbind
  out
}

#' Checks if argument is a \code{brmsterms} object
#'
#' @param x An \R object
#'
#' @seealso \code{\link[brms:brmsterms]{brmsterms}}
#'
#' @export
is.brmsterms <- function(x) {
  inherits(x, "brmsterms")
}

#' Checks if argument is a \code{mvbrmsterms} object
#'
#' @param x An \R object
#'
#' @seealso \code{\link[brms:brmsterms]{brmsterms}}
#'
#' @export
is.mvbrmsterms <- function(x) {
  inherits(x, "mvbrmsterms")
}

# useful for functions that require either of the two objects
is.anybrmsterms <- function(x) {
  is.brmsterms(x) || is.mvbrmsterms(x)
}

is.btl <- function(x) {
  inherits(x, "btl")
}

is.btnl <- function(x) {
  inherits(x, "btnl")
}

# figure out if a certain distributional parameter is predicted
is_pred_dpar <- function(bterms, dpar) {
  stopifnot(is.brmsterms(bterms))
  if (!length(dpar)) {
    return(FALSE)
  }
  mix <- get_mix_id(bterms)
  any(paste0(dpar, mix) %in% names(bterms$dpars))
}

# transform mvbrmsterms objects for use in stan_llh.brmsterms
as.brmsterms <- function(x) {
  stopifnot(is.mvbrmsterms(x), x$rescor)
  families <- ulapply(x$terms, function(y) y$family$family)
  stopifnot(all(families == families[1]))
  out <- structure(list(), class = "brmsterms")
  out$family <- structure(
    list(family = families[1], link = "identity"),
    class = c("brmsfamily", "family")
  )
  out$family$fun <- paste0(out$family$family, "_mv")
  info <- get(paste0(".family_", families[1]))()
  out$family[names(info)] <- info
  out$sigma_pred <- any(ulapply(x$terms,
    function(x) is_pred_dpar(x, "sigma") || has_ad_terms(x, "se")
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

# names of supported term types
all_term_types <- function() {
  c("fe", "re", "sp", "cs", "sm", "gp", "ac", "offset")
}

# avoid ambiguous parameter names
# @param names names to check for ambiguity
# @param bterms a brmsterms object
avoid_dpars <- function(names, bterms) {
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

# check and tidy parameter prefixes
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

# combined parameter prefixes
# @param prefix object from which to extract prefixes
# @param keep_mu keep the 'mu' prefix if available or remove it?
# @param nlp include the 'nlp' prefix for non-linear parameters?
combine_prefix <- function(prefix, keep_mu = FALSE, nlp = FALSE) {
  prefix <- check_prefix(prefix, keep_mu = keep_mu)
  if (is_nlpar(prefix) && nlp) {
    prefix$dpar <- "nlp"
  }
  prefix <- lapply(prefix, usc)
  sub("^_", "", do_call(paste0, prefix))
}

# check validity of fixed distributional parameters
check_fdpars <- function(x) {
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

# combine all variables in one formuula
# @param x (list of) formulas or character strings
# @return a formula with all variables on the right-hand side
allvars_formula <- function(..., .env = parent.frame()) {
  out <- rmNULL(c(...))
  out <- collapse(ulapply(out, plus_rhs))
  all_vars <- all_vars(out)
  invalid_vars <- setdiff(all_vars, make.names(all_vars))
  if (length(invalid_vars)) {
    stop2("The following variable names are invalid: ",
          collapse_comma(invalid_vars))
  }
  out <- str2formula(c(out, all_vars), env = .env)
  # remove redundant terms to simplify the formula (#1786)
  formula(terms.formula(out, simplify = TRUE))
}

# conveniently extract a formula of all relevant variables
# @param x any object from which to extract 'allvars'
# @param type predictor type; requires a 'parse_<type>' function
# @return a formula with all variables on the right-hand side
#   or NULL if 'allvars' cannot be found
get_allvars <- function(x, type = "") {
  out <- attr(x, "allvars", TRUE)
  if (is.null(out) && "allvars" %in% names(x)) {
    out <- x[["allvars"]]
  }
  if (is.null(out) && is.formula(x)) {
    type <- as_one_character(type)
    type <- str_if(nzchar(type), type, "fe")
    terms_fun <- get(paste0("terms_", type), mode = "function")
    out <- attr(terms_fun(x), "allvars")
  }
  stopifnot(is.null(out) || is.formula(out))
  out
}

# add 'x' to the right-hand side of a formula
plus_rhs <- function(x) {
  if (is.formula(x)) {
    x <- sub("^[^~]*~", "", formula2str(x))
  }
  if (length(x) && all(nzchar(x))) {
    out <- paste0(" + ", paste(x, collapse = "+"))
  } else {
    out <- " + 1"
  }
  out
}

# like stats::terms but keeps attributes if possible
terms <- function(formula, ...) {
  old_attributes <- attributes(formula)
  formula <- stats::terms(formula, ...)
  new_attributes <- attributes(formula)
  sel_names <- setdiff(names(old_attributes), names(new_attributes))
  attributes(formula)[sel_names] <- old_attributes[sel_names]
  formula
}

is.terms <- function(x) {
  inherits(x, "terms")
}

# combine formulas for distributional parameters
# @param formula1 primary formula from which to take the RHS
# @param formula2 secondary formula used to update the RHS of formula1
# @param lhs character string to define the left-hand side of the output
# @param update a flag to indicate whether updating should be allowed.
#   Defaults to FALSE to maintain backwards compatibility
# @return a formula object
combine_formulas <- function(formula1, formula2, lhs = "", update = FALSE) {
  stopifnot(is.formula(formula1))
  stopifnot(is.null(formula2) || is.formula(formula2))
  lhs <- as_one_character(lhs)
  update <- as_one_logical(update)
  if (is.null(formula2)) {
    rhs <- str_rhs(formula1)
    att <- attributes(formula1)
  } else if (update && has_terms(formula1)) {
    # TODO: decide about intuitive updating behavior
    if (get_nl(formula1) || get_nl(formula2)) {
      stop2("Cannot combine non-linear formulas.")
    }
    old_formula <- eval2(paste0("~ ", str_rhs(formula1)))
    new_formula <- eval2(paste0("~ . + ", str_rhs(formula2)))
    rhs <- str_rhs(update(old_formula, new_formula))
    att <- attributes(formula1)
    att[names(attributes(formula2))] <- attributes(formula2)
  } else {
    rhs <- str_rhs(formula2)
    att <- attributes(formula2)
  }
  out <- eval2(paste0(lhs, " ~ ", rhs))
  attributes(out)[names(att)] <- att
  out
}

# does the formula contain any terms?
# @return TRUE or FALSE
has_terms <- function(formula) {
  stopifnot(is.formula(formula))
  terms <- try(terms(rhs(formula)), silent = TRUE)
  is_try_error(terms) ||
    length(attr(terms, "term.labels")) ||
    length(attr(terms, "offset"))
}

# has a linear formula any terms except overall effects?
has_special_terms <- function(x) {
  if (!is.btl(x)) {
    return(FALSE)
  }
  special_terms <- c("sp", "sm", "gp", "ac", "cs", "offset")
  NROW(x[["re"]]) > 0 || any(lengths(x[special_terms]))
}

# indicate if the predictor term belongs to a non-linear parameter
is_nlpar <- function(x) {
  isTRUE(nzchar(x[["nlpar"]]))
}

# indicate if the intercept should be removed
no_int <- function(x) {
  isFALSE(attr(x, "int", exact = TRUE))
}

# indicate if cell mean coding should be disabled
no_cmc <- function(x) {
  isFALSE(attr(x, "cmc", exact = TRUE))
}

# indicate if centering of the design matrix should be disabled
no_center <- function(x) {
  isFALSE(attr(x, "center", exact = TRUE))
}

# indicate if the design matrix should be handled as sparse
is_sparse <- function(x) {
  isTRUE(attr(x, "sparse", exact = TRUE))
}

# get the decomposition type of the design matrix
get_decomp <- function(x) {
  out <- attr(x, "decomp", exact = TRUE)
  if (is.null(out)) {
    out <- "none"
  }
  as_one_character(out)
}

# extract different types of effects
get_effect <- function(x, ...) {
  UseMethod("get_effect")
}

#' @export
get_effect.default <- function(x, ...) {
  NULL
}

#' @export
get_effect.brmsfit <- function(x, ...) {
  get_effect(x$formula, ...)
}

#' @export
get_effect.brmsformula <- function(x, ...) {
  get_effect(brmsterms(x), ...)
}

#' @export
get_effect.mvbrmsformula <- function(x, ...) {
  get_effect(brmsterms(x), ...)
}

#' @export
get_effect.mvbrmsterms <- function(x, ...) {
  ulapply(x$terms, get_effect, recursive = FALSE, ...)
}

# extract formulas of a certain effect type
# @param target effect type to return
# @param all logical; include effects of nlpars and dpars?
# @return a list of formulas
#' @export
get_effect.brmsterms <- function(x, target = "fe", ...) {
  out <- named_list(c(names(x$dpars), names(x$nlpars)))
  for (dp in names(x$dpars)) {
    out[[dp]] <- get_effect(x$dpars[[dp]], target = target)
  }
  for (nlp in names(x$nlpars)) {
    out[[nlp]] <- get_effect(x$nlpars[[nlp]], target = target)
  }
  unlist(out, recursive = FALSE)
}

#' @export
get_effect.btl <- function(x, target = "fe", ...) {
  x[[target]]
}

#' @export
get_effect.btnl <- function(x, target = "fe", ...) {
  x[[target]]
}

all_terms <- function(x) {
  if (!length(x)) {
    return(character(0))
  }
  if (!is.terms(x)) {
    x <- terms(as.formula(x))
  }
  trim_wsp(attr(x, "term.labels"))
}

# generate a regular expression to extract special terms
# @param type one or more special term types to be extracted
# TODO: rule out expressions such as mi(y) + mi(x)
regex_sp <- function(type = "all") {
  choices <- c("all", "sp", "sm", "gp", "cs", "mmc", "ac", all_sp_types())
  type <- unique(match.arg(type, choices, several.ok = TRUE))
  funs <- c(
    sm = "(s|(t2)|(te)|(ti))",
    gp = "gp", cs = "cse?", mmc = "mmc",
    ac = "((arma)|(ar)|(ma)|(cosy)|(unstr)|(sar)|(car)|(fcor))"
  )
  funs[all_sp_types()] <- all_sp_types()
  if ("sp" %in% type) {
    # allows extracting all 'sp' terms at once
    type <- setdiff(type, "sp")
    type <- union(type, all_sp_types())
  }
  if ("all" %in% type) {
    # allows extracting all special terms at once
    type <- names(funs)
  }
  funs <- funs[type]
  allow_colon <- c("cs", "mmc", "ac")
  inner <- ifelse(names(funs) %in% allow_colon, ".*", "[^:]*")
  out <- paste0("^(", funs, ")\\(", inner, "\\)$")
  paste0("(", out, ")", collapse = "|")
}

# find special terms of a certain type
# @param x formula object of character vector from which to extract terms
# @param type special terms type to be extracted. see regex_sp()
# @param complete check if terms consist completely of single special terms?
# @param ranef include group-level terms?
# @return a character vector of matching terms
find_terms <- function(x, type, complete = TRUE, ranef = FALSE) {
  if (is.formula(x)) {
    x <- all_terms(x)
  } else {
    x <- trim_wsp(as.character(x))
  }
  complete <- as_one_logical(complete)
  ranef <- as_one_logical(ranef)
  regex <- regex_sp(type)
  is_match <- grepl_expr(regex, x)
  if (!ranef) {
    is_match <- is_match & !grepl("\\|", x)
  }
  out <- x[is_match]
  if (complete) {
    matches <- lapply(out, get_matches_expr, pattern = regex)
    # each term may contain only one special function call
    invalid <- out[lengths(matches) > 1L]
    if (!length(invalid)) {
      # each term must be exactly equal to the special function call
      invalid <- out[unlist(matches) != out]
    }
    # TODO: some terms can be part of I() calls (#1520); reflect this here?
    if (length(invalid)) {
      stop2("The term '", invalid[1], "' is invalid in brms syntax.")
    }
  }
  out
}

# validate a terms object (or one that can be coerced to it)
# for use primarily in 'get_model_matrix'
# @param x any R object
# @return a (possibly amended) terms object or NULL
#   if 'x' could not be coerced to a terms object
validate_terms <- function(x) {
  no_int <- no_int(x)
  no_cmc <- no_cmc(x)
  if (is.formula(x) && !is.terms(x)) {
    x <- terms(x)
  }
  if (!is.terms(x)) {
    return(NULL)
  }
  if (no_int || !has_intercept(x) && no_cmc) {
    # allows to remove the intercept without causing cell mean coding
    attr(x, "intercept") <- 1
    attr(x, "int") <- FALSE
  }
  x
}

# checks if the formula contains an intercept
has_intercept <- function(formula) {
  if (is.terms(formula)) {
    out <- as.logical(attr(formula, "intercept"))
  } else {
    formula <- as.formula(formula)
    try_terms <- try(terms(formula), silent = TRUE)
    if (is_try_error(try_terms)) {
      out <- FALSE
    } else {
      out <- as.logical(attr(try_terms, "intercept"))
    }
  }
  out
}

# check if model makes use of the reserved intercept variables
# @param has_intercept does the model have an intercept?
#   if NULL this will be inferred from formula itself
has_rsv_intercept <- function(formula, has_intercept = NULL) {
  .has_rsv_intercept <- function(terms, has_intercept) {
    has_intercept <- as_one_logical(has_intercept)
    intercepts <- c("intercept", "Intercept")
    out <- !has_intercept && any(intercepts %in% all_vars(rhs(terms)))
    return(out)
  }
  if (is.terms(formula)) {
    if (is.null(has_intercept)) {
      has_intercept <- has_intercept(formula)
    }
    return(.has_rsv_intercept(formula, has_intercept))
  }
  formula <- try(as.formula(formula), silent = TRUE)
  if (is_try_error(formula)) {
    return(FALSE)
  }
  if (is.null(has_intercept)) {
    try_terms <- try(terms(formula), silent = TRUE)
    if (is_try_error(try_terms)) {
      return(FALSE)
    }
    has_intercept <- has_intercept(try_terms)
  }
  .has_rsv_intercept(formula, has_intercept)
}

# names of reserved variables
rsv_vars <- function(bterms) {
  stopifnot(is.brmsterms(bterms) || is.mvbrmsterms(bterms))
  .rsv_vars <- function(x) {
    rsv_int <- any(ulapply(x$dpars, has_rsv_intercept))
    if (rsv_int) c("intercept", "Intercept") else NULL
  }
  if (is.mvbrmsterms(bterms)) {
    out <- unique(ulapply(bterms$terms, .rsv_vars))
  } else {
    out <- .rsv_vars(bterms)
  }
  out
}

# are category specific effects present?
has_cs <- function(bterms) {
  length(get_effect(bterms, target = "cs")) > 0L ||
    any(get_re(bterms)$type %in% "cs")
}

# check if category specific effects are allowed
check_cs <- function(bterms) {
  stopifnot(is.btl(bterms) || is.btnl(bterms))
  if (has_cs(bterms)) {
    if (!is_equal(dpar_class(bterms$dpar), "mu")) {
      stop2("Category specific effects are only supported ",
            "for the main parameter 'mu'.")
    }
    if (!(is.null(bterms$family) || allow_cs(bterms$family))) {
      stop2("Category specific effects are not supported for this family.")
    }
    if (needs_ordered_cs(bterms$family)) {
      warning2("Category specific effects for this family should be ",
               "considered experimental and may have convergence issues.")
    }
  }
  invisible(NULL)
}

# check for the presence of helper functions accidentally used
# within a formula instead of added to bf(). See #1103
check_accidental_helper_functions <- function(formula) {
  terms <- all_terms(formula)
  # see help("brmsformula-helpers") for the list of functions
  funs <- c("nlf", "lf", "acformula", "set_nl", "set_rescor", "set_mecor")
  regex <- paste0("(", funs, ")", collapse = "|")
  regex <- paste0("^(", regex, ")\\(")
  matches <- get_matches(regex, terms, first = TRUE)
  matches <- sub("\\($", "", matches)
  matches <- unique(matches)
  matches <- matches[nzchar(matches)]
  for (m in matches) {
    loc <- utils::find(m, mode = "function")
    if (is_equal(loc[1], "package:brms")) {
      stop2("Function '", m, "' should not be part of the right-hand side ",
            "of a formula. See help('brmsformula-helpers') for the correct syntax.")
    }
  }
  invisible(TRUE)
}

# extract names of variables added via the 'unused' argument
get_unused_arg_vars <- function(x, ...) {
  UseMethod("get_unused_arg_vars")
}

#' @export
get_unused_arg_vars.brmsformula <- function(x, ...) {
  all_vars(attr(x$formula, "unused"))
}

#' @export
get_unused_arg_vars.mvbrmsformula <- function(x, ...) {
  unique(ulapply(x$forms, get_unused_arg_vars, ...))
}

#' @export
get_unused_arg_vars.brmsterms <- function(x, ...) {
  all_vars(x$unused)
}

#' @export
get_unused_arg_vars.mvbrmsterms <- function(x, ...) {
  unique(ulapply(x$terms, get_unused_arg_vars, ...))
}

# extract elements from objects
# @param x object from which to extract elements
# @param name name of the element to be extracted
get_element <- function(x, name, ...) {
  UseMethod("get_element")
}

#' @export
get_element.default <- function(x, name, ...) {
  x[[name]]
}

#' @export
get_element.mvbrmsformula <- function(x, name, ...) {
  lapply(x$forms, get_element, name = name, ...)
}

#' @export
get_element.mvbrmsterms <- function(x, name, ...) {
  lapply(x$terms, get_element, name = name, ...)
}
