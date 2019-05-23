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
                                 mv = FALSE, ...) {
  x <- validate_formula(formula, family = family, autocor = autocor)
  mv <- as_one_logical(mv)
  rescor <- mv && isTRUE(x$rescor)
  mecor <- isTRUE(x$mecor)
  formula <- x$formula
  family <- x$family
  autocor <- x$autocor
  y <- nlist(formula, family, autocor, mv, rescor, mecor) 
  class(y) <- "brmsterms"
  
  if (check_response) {
    # extract response variables
    y$respform <- validate_resp_formula(formula, empty_ok = FALSE)
    if (mv) {
      y$resp <- parse_resp(y$respform) 
    } else {
      y$resp <- ""
    }
  }
  
  # extract addition arguments
  adforms <- parse_ad(formula, family, check_response)
  advars <- str2formula(ulapply(adforms, all_vars))
  y$adforms[names(adforms)] <- adforms
  
  # copy stuff from the formula to parameter 'mu'
  str_rhs_form <- formula2str(rhs(formula))
  rhs_needed <- FALSE
  if (is.mixfamily(family)) {
    for (i in seq_along(family$mix)) {
      mui <- paste0("mu", i)
      if (!is.formula(x$pforms[[mui]])) {
        x$pforms[[mui]] <- eval2(paste0(mui, str_rhs_form))
        attributes(x$pforms[[mui]]) <- attributes(formula)
        rhs_needed <- TRUE
      }
    }
  } else if (conv_cats_dpars(x$family)) {
    mu_dpars <- str_subset(x$family$dpars, "^mu")
    for (dp in mu_dpars) {
      if (!is.formula(x$pforms[[dp]])) {
        x$pforms[[dp]] <- eval2(paste0(dp, str_rhs_form))
        attributes(x$pforms[[dp]]) <- attributes(formula)
        rhs_needed <- TRUE
      }
    }
  } else {
    if (!is.formula(x$pforms[["mu"]])) { 
      x$pforms$mu <- eval2(paste0("mu", str_rhs_form))
      attributes(x$pforms$mu) <- attributes(formula)
      rhs_needed <- TRUE
    }
    x$pforms <- move2start(x$pforms, "mu")
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
  dpars <- intersect(names(x$pforms), valid_dpars(family))
  dpar_forms <- x$pforms[dpars]
  nlpars <- setdiff(names(x$pforms), dpars)
  
  y$dpars <- named_list(dpars)
  for (dp in dpars) {
    if (get_nl(dpar_forms[[dp]])) {
      y$dpars[[dp]] <- parse_nlf(dpar_forms[[dp]], nlpars, resp)
    } else {
      y$dpars[[dp]] <- parse_lf(dpar_forms[[dp]])
    }
    y$dpars[[dp]]$family <- dpar_family(family, dp)
    y$dpars[[dp]]$dpar <- dp
    y$dpars[[dp]]$resp <- resp
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
        y$nlpars[[nlp]] <- parse_nlf(nlpar_forms[[nlp]], nlpars, resp)
      } else {
        y$nlpars[[nlp]] <- parse_lf(nlpar_forms[[nlp]])
      }
      y$nlpars[[nlp]]$nlpar <- nlp
      y$nlpars[[nlp]]$resp <- resp
    }
    used_nlpars <- ulapply(c(y$dpars, y$nlpars), "[[", "used_nlpars")
    unused_nlpars <- setdiff(nlpars, used_nlpars)
    if (length(unused_nlpars)) {
      stop2(
        "The parameter '", unused_nlpars[1], "' is not a ", 
        "valid distributional or non-linear parameter. ",
        "Did you forget to set 'nl = TRUE'?"
      )
    }
    # sort non-linear parameters after dependency
    used_nlpars <- lapply(y$nlpars, "[[", "used_nlpars")
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
  # check for illegal use of cs terms
  if (has_cs(y) && !(is.null(family) || allow_cs(family))) {
    stop2("Category specific effects require families ", 
          "'sratio', 'cratio', or 'acat'.")
  }
  # parse autocor formula
  y$time <- parse_time(autocor)
  if (!is.null(y$dpars[["mu"]])) {
    y$dpars$mu$autocor <- autocor
    y$dpars$mu$time <- y$time
  }
  
  # make a formula containing all required variables
  lhsvars <- if (resp_rhs_all) all_vars(y$respform)
  y$allvars <- allvars_formula(
    lhsvars, advars, lapply(y$dpars, get_allvars), 
    lapply(y$nlpars, get_allvars), y$time$allvars
  )
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
  x$mecor <- isTRUE(x$mecor)
  out <- structure(list(), class = "mvbrmsterms")
  out$terms <- named_list(names(x$forms))
  for (i in seq_along(out$terms)) {
    x$forms[[i]]$rescor <- x$rescor
    x$forms[[i]]$mecor <- x$mecor
    out$terms[[i]] <- parse_bf(x$forms[[i]], mv = TRUE, ...)
  }
  out$allvars <- allvars_formula(lapply(out$terms, get_allvars))
  # required to find variables used solely in the response part
  lhs_resp <- function(x) deparse_combine(lhs(x$respform)[[2]])
  out$respform <- paste0(ulapply(out$terms, lhs_resp), collapse = ",")
  out$respform <- formula(paste0("cbind(", out$respform, ") ~ 1"))
  out$responses <- ulapply(out$terms, "[[", "resp")
  out$rescor <- x$rescor
  out$mecor <- x$mecor
  out
}

# parse linear/additive formulas
# @param formula an ordinary model formula
# @return a 'btl' object
parse_lf <- function(formula) {
  formula <- rhs(as.formula(formula))
  y <- nlist(formula)
  types <- c("fe", "re", "sp", "cs", "sm", "gp", "offset")
  for (t in types) {
    tmp <- do_call(paste0("parse_", t), list(formula))
    if (is.data.frame(tmp) || is.formula(tmp)) {
      y[[t]] <- tmp 
    }
  }
  y$allvars <- allvars_formula(
    get_allvars(y$fe), get_allvars(y$re),
    get_allvars(y$cs), get_allvars(y$sp),
    get_allvars(y$sm), get_allvars(y$gp),
    get_allvars(y$offset)
  )
  environment(y$allvars) <- environment(formula)
  structure(y, class = "btl")
}

# parse non-linear formulas
# @param formula non-linear model formula
# @param nlpars names of all non-linear parameters
# @param resp optional name of a response variable
# @return a 'btnl' object
parse_nlf <- function(formula, nlpars, resp = "") {
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
  y$allvars <- allvars_formula(covars)
  environment(y$allvars) <- environment(formula)
  y$loop <- loop
  structure(y, class = "btnl")
}

# extract addition arguments out of formula
# @return a list of formulas each containg a single addition term
parse_ad <- function(formula, family = NULL, check_response = TRUE) {
  x <- list()
  ad_funs <- lsp("brms", what = "exports", pattern = "^resp_")
  ad_funs <- sub("^resp_", "", ad_funs)
  families <- family_names(family)
  if (is.family(family) && any(nzchar(families))) {
    str_formula <- formula2str(formula)
    ad <- get_matches("(?<=\\|)[^~]*(?=~)", str_formula, perl = TRUE)
    valid_ads <- family_info(family, "ad")
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
  }
  x
}

# extract fixed effects terms
parse_fe <- function(formula) {
  all_terms <- all_terms(formula)
  sp_terms <- find_terms(all_terms, "all", complete = FALSE)
  re_terms <- all_terms[grepl("\\|", all_terms)]
  int_term <- attr(terms(formula), "intercept")
  out <- setdiff(all_terms, c(sp_terms, re_terms))
  out <- paste(c(int_term, out), collapse = "+")
  out <- str2formula(out)
  attr(out, "allvars") <- allvars_formula(out)
  attr(out, "decomp") <- get_decomp(formula)
  if (has_rsv_intercept(out)) {
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
parse_re <- function(formula) {
  re_terms <- get_re_terms(formula, brackets = FALSE)
  re_terms <- split_re_terms(re_terms)
  re_parts <- re_parts(re_terms)
  out <- allvars <- vector("list", length(re_terms))
  type <- attr(re_terms, "type")
  for (i in seq_along(re_terms)) {
    id <- gsub("\\|", "", re_parts$mid[i])
    if (!nzchar(id)) id <- NA
    gcall <- eval2(re_parts$rhs[i])
    form <- str2formula(re_parts$lhs[i])
    group <- paste0(gcall$type, collapse(gcall$groups))
    out[[i]] <- data.frame(
      group = group, gtype = gcall$type, 
      gn = i, id = id, type = type[i],
      cor = substr(re_parts$mid[i], 1, 2) != "||",
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
  if (length(out)) {
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
  } else {
    out <- data.frame(
      group = character(0), gtype = character(0),
      gn = numeric(0), id = numeric(0), type = character(0), 
      cor = logical(0), form = character(0)
    )
  }
  out
}

# extract category specific terms for ordinal models
parse_cs <- function(formula) {
  out <- find_terms(formula, "cs")
  if (length(out)) {
    out <- ulapply(out, eval2)
    out <- str2formula(out)
    attr(out, "allvars") <- allvars_formula(out)
    # do not test whether variables were supplied to 'cs'
    # to allow category specific group-level intercepts
    attr(out, "int") <- FALSE
  }
  out
}

# extract special effects terms 
parse_sp <- function(formula) {
  types <- c("mo", "me", "mi")
  out <- find_terms(formula, types, complete = FALSE)
  if (length(out)) {
    uni_mo <- rm_wsp(get_matches_expr(regex_sp("mo"), out))
    uni_me <- rm_wsp(get_matches_expr(regex_sp("me"), out))
    uni_mi <- rm_wsp(get_matches_expr(regex_sp("mi"), out))
    out <- str2formula(out)
    attr(out, "int") <- FALSE
    attr(out, "uni_mo") <- uni_mo
    attr(out, "uni_me") <- uni_me
    attr(out, "uni_mi") <- uni_mi
    uni_terms <- c(uni_mo, uni_me, uni_mi)
    attr(out, "allvars") <- sp_fake_formula(uni_terms)
  }
  out
}

# extract spline terms
parse_sm <- function(formula) {
  out <- find_terms(formula, "sm")
  if (length(out)) {
    if (any(grepl("^(te|ti)\\(", out))) {
      stop2("Tensor product smooths 'te' and 'ti' are not yet ", 
            "implemented in brms. Consider using 't2' instead.")
    }
    out <- str2formula(out)
    attr(out, "allvars") <- mgcv::interpret.gam(out)$fake.formula
  }
  out
}

# extract gaussian process terms
parse_gp <- function(formula) {
  out <- find_terms(formula, "gp")
  if (length(out)) {
    eterms <- lapply(out, eval2)
    covars <- lapply(eterms, "[[", "term")
    byvars <- lapply(eterms, "[[", "by")
    allvars <- str2formula(unlist(c(covars, byvars)))
    allvars <- str2formula(all_vars(allvars))
    if (!length(all_vars(allvars))) {
      stop2("No variable supplied to function 'gp'.")
    }
    out <- str2formula(out)
    attr(out, "allvars") <- allvars
  }
  out
}

# extract offset terms
parse_offset <- function(formula) {
  out <- character(0)
  terms <- terms(as.formula(formula))
  pos <- attr(terms, "offset")
  if (!is.null(pos)) {
    vars <- attr(terms, "variables")
    out <- ulapply(pos, function(i) deparse(vars[[i + 1]]))
    out <- str2formula(out)
    attr(out, "allvars") <- str2formula(all_vars(out))
  }
  out
}

# extract multiple covariates in multi-membership terms
parse_mmc <- function(formula) {
  out <- find_terms(formula, "mmc")
  if (length(out)) {
    out <- str2formula(out)
    attr(out, "allvars") <- allvars_formula(out)
    attr(out, "int") <- FALSE
  }
  out
}

# extract response variable names
# assumes multiple response variables to be combined via mvbind
parse_resp <- function(formula, check_names = TRUE) {
  formula <- lhs(as.formula(formula))
  if (is.null(formula)) {
    return(NULL)
  }
  expr <- validate_resp_formula(formula)[[2]]
  if (length(expr) <= 1L) {
    out <- deparse_no_string(expr)
  } else {
    str_fun <- deparse_no_string(expr[[1]]) 
    use_mvbind <- identical(str_fun, "mvbind")
    use_cbind <- identical(str_fun, "cbind")
    if (use_mvbind) {
      out <- ulapply(expr[-1], deparse_no_string)
    } else if (use_cbind) {
      # deprecated as of brms 2.7.2
      warning2("Using 'cbind' for multivariate models is ", 
               "deprecated. Please use 'mvbind' instead.")
      out <- ulapply(expr[-1], deparse_no_string)
    } else {
      out <- deparse_no_string(expr) 
    }
  }
  if (check_names) {
    out <- make_stan_names(out)
  }
  out
}

# extract time and grouping variables for autocorrelation structures
# @param autocor object of class 'cor_brms'
# @return a list with elements 'time', 'group', and 'allvars'
parse_time <- function(autocor) {
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
  time_vars <- all_vars(time)
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
  stopif_illegal_group(group)
  group_vars <- all_vars(group)
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

# transform mvbrmsterms objects for use in stan_llh.brmsterms
as.brmsterms <- function(x) {
  stopifnot(is.mvbrmsterms(x), x$rescor)
  families <- ulapply(x$terms, function(y) y$family$family)
  stopifnot(all(families == families[1]))
  out <- structure(list(), class = "brmsterms")
  out$family <- structure(
    list(family = paste0(families[1], "_mv"), link = "identity"),
    class = c("brmsfamily", "family")
  )
  info <- get(paste0(".family_", families[1]))()
  out$family[names(info)] <- info
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
allvars_formula <- function(...) {
  out <- rmNULL(c(...))
  out <- collapse(ulapply(out, plus_rhs))
  out <- str2formula(c(out, all_vars(out)))
  update(out, ~ .)
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
    parse_fun <- get(paste0("parse_", type), mode = "function")
    out <- attr(parse_fun(x), "allvars")
  }
  stopifnot(is.null(out) || is.formula(out))
  out
}

# add 'x' to the right-hand side of a formula
plus_rhs <- function(x) {
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
  get_effect(parse_bf(x), ...)
}

#' @export
get_effect.mvbrmsformula <- function(x, ...) {
  get_effect(parse_bf(x), ...)
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
get_effect.brmsterms <- function(x, target = "fe", all = TRUE, ...) {
  if (all) {
    out <- named_list(c(names(x$dpars), names(x$nlpars)))
    for (dp in names(x$dpars)) {
      out[[dp]] <- get_effect(x$dpars[[dp]], target = target)
    }
    for (nlp in names(x$nlpars)) {
      out[[nlp]] <- get_effect(x$nlpars[[nlp]], target = target)
    }
  } else {
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
  NULL
}

# extract variable names used in addition terms
get_advars <- function(x, ...) {
  UseMethod("get_advars")
}

#' @export
get_advars.brmsterms <- function(x, ad, ...) {
  ad <- as_one_character(ad)
  all_vars(x$adforms[[ad]])
}

#' @export
get_advars.mvbrmsterms <- function(x, ad, ...) {
  unique(ulapply(x$terms, get_advars, ad = ad, ...))
}

# extract information about smooth terms
# @param x either a formula or a list containing an element "sm"
# @param data data.frame containing the covariates
tidy_smef <- function(x, data) {
  if (is.formula(x)) {
    x <- parse_bf(x, check_response = FALSE)$dpars$mu
  }
  form <- x[["sm"]] 
  if (!is.formula(form)) {
    return(empty_data_frame())
  }
  out <- data.frame(term = all_terms(form), stringsAsFactors = FALSE)
  nterms <- nrow(out)
  out$sfun <- get_matches("^[^\\(]+", out$term)
  out$vars <- out$byvars <- out$covars <- vector("list", nterms)
  for (i in seq_len(nterms)) {
    sm <- eval2(out$term[i])
    out$covars[[i]] <- sm$term
    if (sm$by != "NA") {
      out$byvars[[i]] <- sm$by
    }
    out$vars[[i]] <- c(out$covars[[i]], out$byvars[[i]])
  }
  out$label <- paste0(out$sfun, rename(ulapply(out$vars, collapse)))
  # prepare information inferred from the data
  sdata <- data_sm(x, data, knots = attr(data, "knots"))
  bylevels <- attr(sdata$Xs, "bylevels")
  nby <- lengths(bylevels)
  tmp <- vector("list", nterms)
  for (i in seq_len(nterms)) {
    tmp[[i]] <- out[i, , drop = FALSE]
    tmp[[i]]$termnum <- i
    if (nby[i] > 0L) {
      tmp[[i]] <- do_call(rbind, repl(tmp[[i]], nby[i]))
      tmp[[i]]$bylevel <- rm_wsp(bylevels[[i]])
      tmp[[i]]$byterm <- paste0(tmp[[i]]$term, tmp[[i]]$bylevel)
      str_add(tmp[[i]]$label) <- rename(tmp[[i]]$bylevel)
    } else {
      tmp[[i]]$bylevel <- NA
      tmp[[i]]$byterm <- tmp[[i]]$term
    }
  }
  out <- do_call(rbind, tmp)
  out$knots <- sdata[grepl("^knots_", names(sdata))]
  out$nbases <- lengths(out$knots)
  attr(out, "Xs_names") <- colnames(sdata$Xs)
  rownames(out) <- NULL
  out
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

# generate a regular expression to extract special terms
# @param type one or more special term types to be extracted 
regex_sp <- function(type = "all") {
  choices <- c("all", "sp", "sm", "gp", "cs", "mmc", all_sp_types())
  type <- unique(match.arg(type, choices, several.ok = TRUE))
  funs <- c(sm = "(s|(t2)|(te)|(ti))", gp = "gp", cs = "cse?", mmc = "mmc")
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
  allow_colon <- c("cs", "mmc")
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
    x <- as.character(x)
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
    inv <- out[lengths(matches) > 1L]
    if (!length(inv)) {
      # each term must be exactly equal to the special function call
      inv <- out[rm_wsp(unlist(matches)) != out]
    }
    if (length(inv)) {
      stop2("The term '", inv[1], "' is invalid in brms syntax.")
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
  if (is.formula(x) && !inherits(x, "terms")) {
    x <- terms(x)
  }
  if (!inherits(x, "terms")) {
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
  formula <- as.formula(formula)
  try_terms <- try(terms(formula), silent = TRUE)
  if (is(try_terms, "try-error")) {
    out <- FALSE
  } else {
    out <- as.logical(attr(try_terms, "intercept"))
  }
  out
}

# check if model makes use of the reserved intercept variables
has_rsv_intercept <- function(formula) {
  formula <- try(as.formula(formula), silent = TRUE)
  if (is(formula, "try-error")) {
    out <- FALSE
  } else {
    try_terms <- try(terms(formula), silent = TRUE)
    if (is(try_terms, "try-error")) {
      out <- FALSE
    } else {
      has_intercept <- attr(try_terms, "intercept")
      intercepts <- c("intercept", "Intercept")
      out <- !has_intercept && any(intercepts %in% all_vars(rhs(formula)))
    }
  }
  out
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

# check if smooths are present in the model
has_smooths <- function(bterms) {
  length(get_effect(bterms, target = "sm")) > 0L
}

# check if category specific effects are present in the model
has_cs <- function(bterms) {
  length(get_effect(bterms, target = "cs")) > 0L ||
    any(get_re(bterms)$type %in% "cs")
}

# extract names of response categories
extract_cat_names <- function(x, data) {
  stopifnot(is.brmsformula(x) || is.brmsterms(x))
  respform <- validate_resp_formula(x$formula)
  mr <- model.response(model.frame(respform, data))
  if (is_ordinal(x) && is.numeric(mr)) {
    out <- as.character(seq_len(max(mr)))
  } else if (has_multicol(x)) {
    mr <- as.matrix(mr)
    out <- as.character(colnames(mr))
    if (!length(out)) {
      out <- as.character(seq_cols(mr))
    }
  } else {
    out <- levels(factor(mr))
  }
  out
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

# extract variable names used in autocor structures
get_autocor_vars <- function(x, ...) {
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

# convenient wrapper around 'get_autocor_vars'
get_ac_groups <- function(x, ...) {
  get_autocor_vars(x, var = "group", ...)
}

# extract truncation boundaries
# @param incl_family include the family in the derivation of the bounds?
# @param stan return bounds in form of Stan syntax?
get_bounds <- function(bterms, data = NULL, incl_family = FALSE, 
                       stan = FALSE) {
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

# get boundaries of response distributions
get_family_bounds <- function(bterms) {
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

# indicate if the model is censored
has_cens <- function(bterms, data = NULL) {
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

# check if addition argument 'subset' ist used in the model
has_subset <- function(bterms) {
  .has_subset <- function(x) {
    is.formula(x$adforms$subset)
  }
  if (is.brmsterms(bterms)) {
    out <- .has_subset(bterms)
  } else if (is.mvbrmsterms(bterms)) {
    out <- any(ulapply(bterms$terms, .has_subset))
  } else {
    out <- FALSE
  }
  out 
}
