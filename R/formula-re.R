# This file contains functions dealing with the extended
# lme4-like formula syntax to specify group-level terms

#' Set up basic grouping terms in \pkg{brms}
#'
#' Function used to set up a basic grouping term in \pkg{brms}.
#' The function does not evaluate its arguments --
#' it exists purely to help set up a model with grouping terms.
#' \code{gr} is called implicitly inside the package
#' and there is usually no need to call it directly.
#'
#' @param ... One or more terms containing grouping factors.
#' @param by An optional factor variable, specifying sub-populations of the
#'   groups. For each level of the \code{by} variable, a separate
#'   variance-covariance matrix will be fitted. Levels of the grouping factor
#'   must be nested in levels of the \code{by} variable.
#' @param cor Logical. If \code{TRUE} (the default), group-level terms will be
#'   modelled as correlated.
#' @param id Optional character string. All group-level terms across the model
#'   with the same \code{id} will be modeled as correlated (if \code{cor} is
#'   \code{TRUE}). See \code{\link{brmsformula}} for more details.
#' @param cov An optional matrix which is proportional to the withon-group
#'   covariance matrix of the group-level effects. All levels of the grouping
#'   factor should appear as rownames of the corresponding matrix. This argument
#'   can be used, among others, to model pedigrees and phylogenetic effects. See
#'   \code{vignette("brms_phylogenetics")} for more details. By default, levels
#'   of the same grouping factor are modeled as independent of each other.
#' @param dist Name of the distribution of the group-level effects.
#' Currently \code{"gaussian"} is the only option.
#'
#' @seealso \code{\link{brmsformula}}
#'
#' @examples
#' \dontrun{
#' # model using basic lme4-style formula
#' fit1 <- brm(count ~ Trt + (1|patient), data = epilepsy)
#' summary(fit1)
#'
#' # equivalent model using 'gr' which is called anyway internally
#' fit2 <- brm(count ~ Trt + (1|gr(patient)), data = epilepsy)
#' summary(fit2)
#'
#' # include Trt as a by variable
#' fit3 <- brm(count ~ Trt + (1|gr(patient, by = Trt)), data = epilepsy)
#' summary(fit3)
#' }
#'
#' @export
gr <- function(..., by = NULL, cor = TRUE, id = NA,
               cov = NULL, dist = "gaussian") {
  label <- deparse0(match.call())
  groups <- as.character(as.list(substitute(list(...)))[-1])
  if (length(groups) > 1L) {
    stop2("Grouping structure 'gr' expects only a single grouping term")
  }
  stopif_illegal_group(groups[1])
  cor <- as_one_logical(cor)
  id <- as_one_character(id, allow_na = TRUE)
  by <- substitute(by)
  if (!is.null(by)) {
    by <- deparse0(by)
  } else {
    by <- ""
  }
  cov <- substitute(cov)
  if (!is.null(cov)) {
    cov <- all.vars(cov)
    if (length(cov) != 1L) {
      stop2("Argument 'cov' must contain exactly one variable.")
    }
  } else {
    cov <- ""
  }
  dist <- match.arg(dist, c("gaussian", "student"))
  byvars <- all_vars(by)
  allvars <- str2formula(c(groups, byvars))
  nlist(groups, allvars, label, by, cor, id, cov, dist, type = "")
}

#' Set up multi-membership grouping terms in \pkg{brms}
#'
#' Function to set up a multi-membership grouping term in \pkg{brms}.
#' The function does not evaluate its arguments --
#' it exists purely to help set up a model with grouping terms.
#'
#' @inheritParams gr
#' @param weights A matrix specifying the weights of each member.
#'  It should have as many columns as grouping terms specified in \code{...}.
#'  If \code{NULL} (the default), equally weights are used.
#' @param by An optional factor matrix, specifying sub-populations of the
#'   groups. It should have as many columns as grouping terms specified in
#'   \code{...}. For each level of the \code{by} variable, a separate
#'   variance-covariance matrix will be fitted. Levels of the grouping factor
#'   must be nested in levels of the \code{by} variable matrix.
#' @param scale Logical; if \code{TRUE} (the default),
#'  weights are standardized in order to sum to one per row.
#'  If negative weights are specified, \code{scale} needs
#'  to be set to \code{FALSE}.
#'
#' @seealso \code{\link{brmsformula}}, \code{\link{mmc}}
#'
#' @examples
#' \dontrun{
#' # simulate some data
#' dat <- data.frame(
#'  y = rnorm(100), x1 = rnorm(100), x2 = rnorm(100),
#'  g1 = sample(1:10, 100, TRUE), g2 = sample(1:10, 100, TRUE)
#' )
#'
#' # multi-membership model with two members per group and equal weights
#' fit1 <- brm(y ~ x1 + (1|mm(g1, g2)), data = dat)
#' summary(fit1)
#'
#' # weight the first member two times for than the second member
#' dat$w1 <- rep(2, 100)
#' dat$w2 <- rep(1, 100)
#' fit2 <- brm(y ~ x1 + (1|mm(g1, g2, weights = cbind(w1, w2))), data = dat)
#' summary(fit2)
#'
#' # multi-membership model with level specific covariate values
#' dat$xc <- (dat$x1 + dat$x2) / 2
#' fit3 <- brm(y ~ xc + (1 + mmc(x1, x2) | mm(g1, g2)), data = dat)
#' summary(fit3)
#' }
#'
#' @export
mm <- function(..., weights = NULL, scale = TRUE, by = NULL, cor = TRUE,
               id = NA, cov = NULL, dist = "gaussian") {
  label <- deparse0(match.call())
  groups <- as.character(as.list(substitute(list(...)))[-1])
  if (length(groups) < 2) {
    stop2("Multi-membership terms require at least two grouping variables.")
  }
  for (i in seq_along(groups)) {
    stopif_illegal_group(groups[i])
  }
  cor <- as_one_logical(cor)
  id <- as_one_character(id, allow_na = TRUE)
  by <- substitute(by)
  if (!is.null(by)) {
    by <- deparse0(by)
  } else {
    by <- ""
  }
  cov <- substitute(cov)
  if (!is.null(cov)) {
    cov <- all.vars(cov)
    if (length(cov) != 1L) {
      stop2("Argument 'cov' must contain exactly one variable.")
    }
  } else {
    cov <- ""
  }
  dist <- match.arg(dist, c("gaussian", "student"))
  scale <- as_one_logical(scale)
  weights <- substitute(weights)
  weightvars <- all_vars(weights)
  byvars <- all_vars(by)
  allvars <- str2formula(c(groups, weightvars, byvars))
  if (!is.null(weights)) {
    weights <- str2formula(deparse_no_string(weights))
    attr(weights, "scale") <- scale
    weightvars <- str2formula(weightvars)
  }
  nlist(
    groups, weights, weightvars, allvars, label,
    by, cor, id, cov, dist, type = "mm"
  )
}

#' Multi-Membership Covariates
#'
#' Specify covariates that vary over different levels
#' of multi-membership grouping factors thus requiring
#' special treatment. This function is almost solely useful,
#' when called in combination with \code{\link{mm}}.
#' Outside of multi-membership terms it will behave
#' very much like \code{\link{cbind}}.
#'
#' @param ... One or more terms containing covariates
#'   corresponding to the grouping levels specified in \code{\link{mm}}.
#'
#' @return A matrix with covariates as columns.
#'
#' @seealso \code{\link{mm}}
#'
#' @examples
#' \dontrun{
#' # simulate some data
#' dat <- data.frame(
#'   y = rnorm(100), x1 = rnorm(100), x2 = rnorm(100),
#'   g1 = sample(1:10, 100, TRUE), g2 = sample(1:10, 100, TRUE)
#' )
#'
#' # multi-membership model with level specific covariate values
#' dat$xc <- (dat$x1 + dat$x2) / 2
#' fit <- brm(y ~ xc + (1 + mmc(x1, x2) | mm(g1, g2)), data = dat)
#' summary(fit)
#' }
#'
#' @export
mmc <- function(...) {
  dots <- list(...)
  if (any(ulapply(dots, is_like_factor))) {
    stop2("'mmc' requires numeric variables.")
  }
  out <- cbind(...)
  colnames(out) <- paste0("?", colnames(out))
  out
}

# check if the group part of a group-level term is invalid
# @param group the group part of a group-level term
illegal_group_expr <- function(group) {
  group <- as_one_character(group)
  valid_expr <- ":|([^([:digit:]|[:punct:])]|\\.)[[:alnum:]_\\.]*"
  rsv_signs <- c("+", "-", "*", "/", "|", "::")
  nzchar(gsub(valid_expr, "", group)) ||
    any(ulapply(rsv_signs, grepl, x = group, fixed = TRUE))
}

stopif_illegal_group <- function(group) {
  if (illegal_group_expr(group)) {
    stop2(
      "Illegal grouping term '", group, "'. It may contain ",
      "only variable names combined by the symbol ':'"
    )
  }
  invisible(NULL)
}

re_lhs <- function(re_terms) {
  get_matches("^[^\\|]*", re_terms)
}

re_mid <- function(re_terms) {
  get_matches("\\|([^\\|]*\\||)", re_terms)
}

re_rhs <- function(re_terms) {
  sub("^\\|", "", get_matches("\\|[^\\|]*$", re_terms))
}

# extract the three parts of group-level terms
# @param re_terms character vector of RE terms in lme4 syntax
# @return a data.frame with one row per group-level term
re_parts <- function(re_terms) {
  lhs <- re_lhs(re_terms)
  mid <- re_mid(re_terms)
  rhs <- re_rhs(re_terms)
  out <- nlist(lhs, mid, rhs)
  if (any(lengths(out) != length(re_terms))) {
    stop2("Invalid syntax used in group-level terms.")
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

# split nested group-level terms and check for special effects terms
# @param re_terms character vector of RE terms in extended lme4 syntax
split_re_terms <- function(re_terms) {
  if (!length(re_terms)) {
    return(re_terms)
  }
  stopifnot(is.character(re_terms))

  # split after grouping factor terms
  re_parts <- re_parts(re_terms)
  new_re_terms <- vector("list", length(re_terms))
  for (i in seq_along(re_terms)) {
    new_re_rhs <- terms(formula(paste0("~", re_parts$rhs[i])))
    new_re_rhs <- attr(new_re_rhs, "term.labels")
    new_re_rhs <- ifelse(
      !grepl("^(gr|mm)\\(", new_re_rhs),
      paste0("gr(", new_re_rhs, ")"), new_re_rhs
    )
    new_re_terms[[i]] <- paste0(
      re_parts$lhs[i], re_parts$mid[i], new_re_rhs
    )
  }
  re_terms <- unlist(new_re_terms)

  # split after coefficient types
  re_parts <- re_parts(re_terms)
  new_re_terms <- type <- vector("list", length(re_terms))
  for (i in seq_along(re_terms)) {
    lhs_form <- formula(paste("~", re_parts$lhs[i]))
    lhs_all_terms <- all_terms(lhs_form)
    # otherwise varying intercepts cannot be handled reliably
    is_cs_term <- grepl_expr(regex_sp("cs"), lhs_all_terms)
    if (any(is_cs_term) && !all(is_cs_term)) {
      stop2("Please specify category specific effects ",
            "in separate group-level terms.")
    }
    new_lhs <- NULL
    # prepare effects of special terms
    valid_types <- c("sp", "cs", "mmc")
    invalid_types <- c("sm", "gp")
    for (t in c(valid_types, invalid_types)) {
      lhs_tform <- do_call(paste0("terms_", t), list(lhs_form))
      if (is.formula(lhs_tform)) {
        if (t %in% invalid_types) {
          stop2("Cannot handle splines or GPs in group-level terms.")
        }
        new_lhs <- c(new_lhs, formula2str(lhs_tform, rm = 1))
        type[[i]] <- c(type[[i]], t)
      }
    }
    # prepare effects of basic terms
    lhs_terms <- terms(lhs_form)
    fe_form <- terms_fe(lhs_terms)
    fe_terms <- all_terms(fe_form)
    # the intercept lives within not outside of 'cs' terms
    has_intercept <- has_intercept(lhs_terms) && !"cs" %in% type[[i]]
    if (length(fe_terms) || has_intercept) {
      new_lhs <- c(new_lhs, formula2str(fe_form, rm = 1))
      type[[i]] <- c(type[[i]], "")
    }
    # extract information from the mid section of the terms
    rhs_call <- str2lang(re_parts$rhs[i])
    if (re_parts$mid[i] == "||") {
      # ||-syntax overwrites the 'cor' argument
      rhs_call$cor <- FALSE
    }
    gcall <- eval(rhs_call)
    if (gcall$cor) {
      id <- gsub("\\|", "", re_parts$mid[i])
      if (nzchar(id)) {
        # ID-syntax overwrites the 'id' argument
        rhs_call$id <- id
      } else {
        id <- gcall$id
      }
      if (length(new_lhs) > 1 && isNA(id)) {
        # ID is required to model coefficients as correlated
        # if multiple types are provided within the same term
        rhs_call$id <- collapse(sample(0:9, 10, TRUE))
      }
    }
    re_parts$mid[i] <- "|"
    re_parts$rhs[i] <- deparse0(rhs_call)
    new_re_terms[[i]] <- paste0(new_lhs, re_parts$mid[i], re_parts$rhs[i])
    new_re_terms[[i]] <- new_re_terms[[i]][order(type[[i]])]
    type[[i]] <- sort(type[[i]])
  }
  re_terms <- unlist(new_re_terms)
  structure(re_terms, type = unlist(type))
}

# extract group-level terms from a formula of character vector
# @param x formula or character vector
# @param formula return a formula rather than a character string?
# @param brackets include group-level terms in brackets?
get_re_terms <- function(x, formula = FALSE, brackets = TRUE) {
  if (is.formula(x)) {
    x <- all_terms(x)
  }
  re_pos <- grepl("\\|", x)
  out <- x[re_pos]
  if (brackets && length(out)) {
    out <- paste0("(", out, ")")
  }
  if (formula) {
    out <- str2formula(out)
  }
  out
}

# validate the re_formula argument
# @inheritParams extract_draws.brmsfit
# @param formula: formula to match re_formula with
# @return updated re_formula containing only terms existent in formula
check_re_formula <- function(re_formula, formula) {
  old_re_formula <- get_re_terms(formula, formula = TRUE)
  if (is.null(re_formula)) {
    re_formula <- old_re_formula
  } else if (SW(anyNA(re_formula))) {
    re_formula <- ~1
  } else {
    re_formula <- get_re_terms(as.formula(re_formula), formula = TRUE)
    new <- brmsterms(re_formula, check_response = FALSE)$dpars$mu[["re"]]
    old <- brmsterms(old_re_formula, check_response = FALSE)$dpars$mu[["re"]]
    if (NROW(new) && NROW(old)) {
      # compare old and new ranefs
      new_terms <- lapply(new$form, terms)
      found <- rep(FALSE, NROW(new))
      for (i in seq_rows(new)) {
        group <- new$group[[i]]
        old_terms <- lapply(old$form[old$group == group], terms)
        j <- 1
        while (!found[i] && j <= length(old_terms)) {
          new_term_labels <- attr(new_terms[[i]], "term.labels")
          old_term_labels <- attr(old_terms[[j]], "term.labels")
          new_intercept <- attr(new_terms[[i]], "intercept")
          old_intercept <- attr(old_terms[[j]], "intercept")
          found[i] <- isTRUE(
            all(new_term_labels %in% old_term_labels) &&
              new_intercept <= old_intercept
          )
          if (found[i]) {
            # terms have to maintain the original order so that Z_* data
            # and r_* parameters match in 'extract_draws' (fixes issue #844)
            term_matches <- match(new_term_labels, old_term_labels)
            if (is.unsorted(term_matches)) {
              stop2("Order of terms in 're_formula' should match the original order.")
            }
          }
          j <- j + 1
        }
      }
      new <- new[found, ]
      if (NROW(new)) {
        forms <- ulapply(new$form, formula2str, rm = 1)
        groups <- ufrom_list(new$gcall, "label")
        re_terms <- paste("(", forms, "|", groups, ")")
        re_formula <- formula(paste("~", paste(re_terms, collapse = "+")))
      } else {
        re_formula <- ~1
      }
    } else {
      re_formula <- ~1
    }
  }
  re_formula
}

# remove existing group-level terms in formula and
# add valid group-level terms of re_formula
update_re_terms <- function(formula, re_formula) {
  UseMethod("update_re_terms")
}

#' @export
update_re_terms.mvbrmsformula <- function(formula, re_formula) {
  formula$forms <- lapply(formula$forms, update_re_terms, re_formula)
  formula
}

#' @export
update_re_terms.brmsformula <- function(formula, re_formula) {
  formula$formula <- update_re_terms(formula$formula, re_formula)
  formula$pforms <- lapply(formula$pforms, update_re_terms, re_formula)
  formula
}

#' @export
update_re_terms.formula <- function(formula, re_formula = NULL) {
  if (is.null(re_formula) || get_nl(formula)) {
    return(formula)
  }
  re_formula <- check_re_formula(re_formula, formula)
  new_formula <- formula2str(formula)
  old_re_terms <- get_re_terms(formula, brackets = FALSE)
  if (length(old_re_terms)) {
    # remove old group-level terms
    rm_terms <- c(
      paste0("+ (", old_re_terms, ")"),
      paste0("(", old_re_terms, ")"),
      old_re_terms
    )
    new_formula <- rename(new_formula, rm_terms, "")
    if (grepl("~( *\\+*)*$", new_formula)) {
      # lhs only formulas are syntactically invalid
      # also check for trailing '+' signs (#769)
      new_formula <- paste(new_formula, "1")
    }
  }
  # add new group-level terms
  new_re_terms <- get_re_terms(re_formula)
  new_formula <- paste(c(new_formula, new_re_terms), collapse = "+")
  new_formula <- formula(new_formula)
  attributes(new_formula) <- attributes(formula)
  new_formula
}

# extract group-level terms
get_re <- function(x, ...) {
  UseMethod("get_re")
}

#' @export
get_re.default <- function(x, ...) {
  NULL
}

# get group-level information in a data.frame
# @param bterms object of class 'brmsterms'
# @param all logical; include ranefs of additional parameters?
#' @export
get_re.brmsterms <- function(x, ...) {
  re <- named_list(c(names(x$dpars), names(x$nlpars)))
  for (dp in names(x$dpars)) {
    re[[dp]] <- get_re(x$dpars[[dp]])
  }
  for (nlp in names(x$nlpars)) {
    re[[nlp]] <- get_re(x$nlpars[[nlp]])
  }
  do_call(rbind, re)
}

#' @export
get_re.mvbrmsterms <- function(x, ...) {
  do_call(rbind, lapply(x$terms, get_re, ...))
}

#' @export
get_re.btl <- function(x, ...) {
  px <- check_prefix(x)
  re <- x[["re"]]
  if (is.null(re)) {
    re <- empty_re()
  }
  re$resp <- rep(px$resp, nrow(re))
  re$dpar <- rep(px$dpar, nrow(re))
  re$nlpar <- rep(px$nlpar, nrow(re))
  re
}

# gather information on group-level effects
# @param bterms object of class brmsterms
# @param data data.frame containing all model variables
# @param old_levels optional original levels of the grouping factors
# @return a tidy data.frame with the following columns:
#   id: ID of the group-level effect
#   group: name of the grouping factor
#   gn: number of the grouping term within the respective formula
#   coef: name of the group-level effect
#   cn: number of the effect within the ID
#   resp: name of the response variable
#   dpar: name of the distributional parameter
#   nlpar: name of the non-linear parameter
#   cor: are correlations modeled for this effect?
#   ggn: global number of the grouping factor
#   type: special effects type; can be 'sp' or 'cs'
#   gcall: output of functions 'gr' or 'mm'
#   form: formula used to compute the effects
frame_re <- function(bterms, data, old_levels = NULL) {
  data <- combine_groups(data, get_group_vars(bterms))
  re <- get_re(bterms)
  out <- vector("list", nrow(re))
  used_ids <- new_ids <- NULL
  id_groups <- list()
  j <- 1
  for (i in seq_rows(re)) {
    if (!nzchar(re$type[i])) {
      coef <- colnames(get_model_matrix(re$form[[i]], data))
    } else if (re$type[i] == "sp") {
      # TODO: try to avoid having to call frame_sp here
      coef <- frame_sp(re$form[[i]], data)$coef
    } else if (re$type[i] == "mmc") {
      coef <- rename(all_terms(re$form[[i]]))
    } else if (re$type[i] == "cs") {
      resp <- re$resp[i]
      if (nzchar(resp)) {
        stopifnot(is.mvbrmsterms(bterms))
        nthres <- max(get_thres(bterms$terms[[resp]]))
      } else {
        stopifnot(is.brmsterms(bterms))
        nthres <- max(get_thres(bterms))
      }
      indices <- paste0("[", seq_len(nthres), "]")
      coef <- colnames(get_model_matrix(re$form[[i]], data = data))
      coef <- as.vector(t(outer(coef, indices, paste0)))
    }
    avoid_dpars(coef, bterms)
    rdat <- data.frame(
      id = re$id[[i]],
      group = re$group[[i]],
      gn = re$gn[[i]],
      gtype = re$gtype[[i]],
      coef = coef,
      cn = NA,
      resp = re$resp[[i]],
      dpar = re$dpar[[i]],
      nlpar = re$nlpar[[i]],
      ggn = NA,
      cor = re$cor[[i]],
      type = re$type[[i]],
      by = re$gcall[[i]]$by,
      cov = re$gcall[[i]]$cov,
      dist = re$gcall[[i]]$dist,
      stringsAsFactors = FALSE
    )
    bylevels <- NULL
    if (nzchar(rdat$by[1])) {
      bylevels <- eval2(rdat$by[1], data)
      bylevels <- rm_wsp(extract_levels(bylevels))
    }
    rdat$bylevels <- repl(bylevels, nrow(rdat))
    rdat$form <- repl(re$form[[i]], nrow(rdat))
    rdat$gcall <- repl(re$gcall[[i]], nrow(rdat))
    # prepare group-level IDs
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
    out[[i]] <- rdat
  }
  out <- do_call(rbind, c(list(empty_reframe()), out))
  # check for overlap between different group types
  rsv_groups <- out[nzchar(out$gtype), "group"]
  other_groups <- out[!nzchar(out$gtype), "group"]
  inv_groups <- intersect(rsv_groups, other_groups)
  if (length(inv_groups)) {
    inv_groups <- paste0("'", inv_groups, "'", collapse = ", ")
    stop2("Grouping factor names ", inv_groups, " are resevered.")
  }
  # check for duplicated and thus not identified effects
  dup <- duplicated(out[, c("group", "coef", vars_prefix())])
  if (any(dup)) {
    dr <- out[which(dup)[1], ]
    stop2(
      "Duplicated group-level effects are not allowed.\n",
      "Occured for effect '", dr$coef, "' of group '", dr$group, "'."
    )
  }
  if (has_rows(out)) {
    for (id in unique(out$id)) {
      out$cn[out$id == id] <- seq_len(sum(out$id == id))
    }
    out$ggn <- match(out$group, unique(out$group))
    # compute random effects levels
    rsub <- out[!duplicated(out$group), ]
    levels <- named_list(rsub$group)
    for (i in seq_along(levels)) {
      # combine levels of all grouping factors within one grouping term
      levels[[i]] <- unique(ulapply(
        rsub$gcall[[i]]$groups,
        function(g) extract_levels(get(g, data))
      ))
      # fixes issue #1353
      bysel <- out$group == names(levels)[i] &
        nzchar(out$by) & !duplicated(out$by)
      bysel <- which(bysel)
      if (length(bysel) > 1L) {
        stop2("Each grouping factor can only be associated with one 'by' variable.")
      }
      # ensure that a non-NULL by-variable is found if present
      if (length(bysel) == 1L) {
        rsub[i, ] <- out[bysel, ]
      }
      # store information of corresponding by-levels
      if (nzchar(rsub$by[i])) {
        stopifnot(rsub$type[i] %in% c("", "mmc"))
        by <- rsub$by[i]
        bylevels <- rsub$bylevels[[i]]
        byvar <- rm_wsp(eval2(by, data))
        groups <- rsub$gcall[[i]]$groups
        if (rsub$gtype[i] == "mm") {
          byvar <- as.matrix(byvar)
          if (!identical(dim(byvar), c(nrow(data), length(groups)))) {
            stop2(
              "Grouping structure 'mm' expects 'by' to be ",
              "a matrix with as many columns as grouping factors."
            )
          }
          df <- J <- named_list(groups)
          for (k in seq_along(groups)) {
            J[[k]] <- match(get(groups[k], data), levels[[i]])
            df[[k]] <- data.frame(J = J[[k]], by = byvar[, k])
          }
          J <- unlist(J)
          df <- do_call(rbind, df)
        } else {
          J <- match(get(groups, data), levels[[i]])
          df <- data.frame(J = J, by = byvar)
        }
        df <- unique(df)
        if (nrow(df) > length(unique(J))) {
          stop2("Some levels of ", collapse_comma(groups),
                " correspond to multiple levels of '", by, "'.")
        }
        df <- df[order(df$J), ]
        by_per_level <- bylevels[match(df$by, bylevels)]
        attr(levels[[i]], "by") <- by_per_level
      }
    }
    if (!is.null(old_levels)) {
      # for newdata numeration has to depend on the original levels
      set_levels(out) <- old_levels
      set_levels(out, "used") <- levels
    } else {
      set_levels(out) <- levels
    }
    # incorporate deprecated 'cov_ranef' argument
    out <- update_ranef_cov(out, bterms)
  }
  # ordering after IDs matches the order of the posterior draws
  # if multiple IDs are used for the same grouping factor (#835)
  out <- out[order(out$id), , drop = FALSE]
  class(out) <- reframe_class()
  out
}

empty_reframe <- function() {
  out <- data.frame(
    id = numeric(0), group = character(0), gn = numeric(0),
    coef = character(0), cn = numeric(0), resp = character(0),
    dpar = character(0), nlpar = character(0), ggn = numeric(0),
    cor = logical(0), type = character(0), form = character(0),
    stringsAsFactors = FALSE
  )
  class(out) <- reframe_class()
  out
}

empty_re <- function() {
  data.frame(
    group = character(0), gtype = character(0),
    gn = numeric(0), id = numeric(0), type = character(0),
    cor = logical(0), form = character(0)
  )
}

reframe_class <- function() {
  c("reframe", "data.frame")
}

is.reframe <- function(x) {
  inherits(x, "reframe")
}

# extract names of all grouping variables
get_group_vars <- function(x, ...) {
  UseMethod("get_group_vars")
}

#' @export
get_group_vars.brmsfit <- function(x, ...) {
  get_group_vars(x$formula, ...)
}

#' @export
get_group_vars.default <- function(x, ...) {
  get_group_vars(brmsterms(x), ...)
}

#' @export
get_group_vars.brmsterms <- function(x, ...) {
  .get_group_vars(x, ...)
}

#' @export
get_group_vars.mvbrmsterms <- function(x, ...) {
  .get_group_vars(x, ...)
}

.get_group_vars <- function(x, ...) {
  out <- c(get_re_groups(x), get_me_groups(x), get_ac_groups(x))
  out <- out[nzchar(out)]
  if (length(out)) {
    c(out) <- unlist(strsplit(out, ":"))
    out <- sort(unique(out))
  }
  out
}

# get names of grouping variables of re terms
get_re_groups <- function(x, ...) {
  ufrom_list(get_re(x)$gcall, "groups")
}

# extract information about groups with a certain distribution
get_dist_groups <- function(reframe, dist) {
  out <- subset2(reframe, dist = dist)
  out[!duplicated(out$group), c("group", "ggn", "id")]
}

# extract names of group-level effects
# @param reframe output of frame_re()
# @param group optional name of a grouping factor for
#   which to extract effect names
# @param bylevels optional names of 'by' levels for
#    which to extract effect names
# @return a vector of character strings
get_rnames <- function(reframe, group = NULL, bylevels = NULL) {
  stopifnot(is.data.frame(reframe))
  if (!is.null(group)) {
    group <- as_one_character(group)
    reframe <- subset2(reframe, group = group)
  }
  stopifnot(length(unique(reframe$group)) == 1L)
  out <- paste0(usc(combine_prefix(reframe), "suffix"), reframe$coef)
  if (isTRUE(nzchar(reframe$by[1]))) {
    if (!is.null(bylevels)) {
      stopifnot(all(bylevels %in% reframe$bylevels[[1]]))
    } else {
      bylevels <- reframe$bylevels[[1]]
    }
    bylabels <- paste0(reframe$by[1], bylevels)
    out <- outer(out, bylabels, paste, sep = ":")
  }
  out
}

# validate within-group covariance matrices
# @param M a matrix to be validated
validate_recov_matrix <- function(M) {
  M <- as.matrix(M)
  if (!isSymmetric(unname(M))) {
    stop2("Within-group covariance matrices must be symmetric.")
  }
  found_levels <- rownames(M)
  if (is.null(found_levels)) {
    found_levels <- colnames(M)
  }
  if (is.null(found_levels)) {
    stop2("Row or column names are required for within-group covariance matrices.")
  }
  rownames(M) <- colnames(M) <- found_levels
  evs <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  if (min(evs) <= 0) {
    stop2("Within-group covariance matrices must be positive definite.")
  }
  M
}

# check validity of the 'cov_ranef' argument
# argument 'cov_ranef' is deprecated as of version 2.12.5
validate_cov_ranef <- function(cov_ranef) {
  if (is.null(cov_ranef)) {
    return(cov_ranef)
  }
  warning2(
    "Argument 'cov_ranef' is deprecated and will be removed in the future. ",
    "Please use argument 'cov' in function 'gr' instead."
  )
  cr_names <- names(cov_ranef)
  cr_is_named <- length(cr_names) && all(nzchar(cr_names))
  if (!is.list(cov_ranef) || !cr_is_named) {
    stop2("'cov_ranef' must be a named list.")
  }
  if (any(duplicated(cr_names))) {
    stop2("Names of 'cov_ranef' must be unique.")
  }
  cov_ranef
}

# update 'reframe' according to information in 'cov_ranef'
# argument 'cov_ranef' is deprecated as of version 2.12.5
update_ranef_cov <- function(reframe, bterms) {
  cr_names <- names(bterms$cov_ranef)
  if (!length(cr_names)) {
    return(reframe)
  }
  unused_names <- setdiff(cr_names, reframe$group)
  if (length(unused_names)) {
    stop2("The following elements of 'cov_ranef' are unused: ",
          collapse_comma(unused_names))
  }
  has_cov <- reframe$group %in% cr_names
  reframe$cov[has_cov] <- reframe$group[has_cov]
  reframe
}

# extract 'cov_ranef' for storage in 'data2'
# @param x a list-like object
get_data2_cov_ranef <- function(x) {
  x[["cov_ranef"]]
}
