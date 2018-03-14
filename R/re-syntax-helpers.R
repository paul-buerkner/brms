# This file contains functions dealing with the extended 
# lme4-like formula syntax to specify group-level terms

illegal_group_expr <- function(group) {
  # check if the group part of a group-level term is invalid
  # Args:
  #  group: the group part of a group-level term
  valid_expr <- ":|[^([:digit:]|[:punct:])][[:alnum:]_\\.]*"
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

get_groups <- function(x) {
  # TODO: merge with get_group_vars
  if (!(is.brmsterms(x) || is.mvbrmsterms(x))) {
    x <- parse_bf(x)
  }
  out <- ulapply(get_re(x)$gcall, "[[", "groups")
  out <- c(out, get_autocor_vars(x, "group"))
  unique(out[nzchar(out)])
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

re_parts <- function(re_terms) {
  # extract the three parts of group-level terms
  # Returns:
  #   A data.frame with one row per group-level term
  lhs <- re_lhs(re_terms)
  mid <- re_mid(re_terms)
  rhs <- re_rhs(re_terms)
  out <- nlist(lhs, mid, rhs)
  if (any(lengths(out) != length(re_terms))) {
    stop2("Invalid syntax used in group-level terms.")
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

split_re_terms <- function(re_terms) {
  # split nested group-level terms 
  # and check for special effects terms
  # Args:
  #   re_terms: group-level terms in extended lme4 syntax
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
      lhs_tform <- do.call(paste0("parse_", t), list(lhs_form))
      if (is.formula(lhs_tform)) {
        if (t %in% invalid_types) {
          stop2("Cannot handle splines or GPs in group-level terms.")
        }
        new_lhs <- c(new_lhs, formula2str(lhs_tform, rm = 1))
        type[[i]] <- c(type[[i]], t)
      }
    }
    # prepare effects of basic terms
    fe_form <- parse_fe(lhs_form)
    fe_terms <- all_terms(fe_form)
    has_intercept <- attr(terms(fe_form), "intercept")
    if (length(fe_terms) || has_intercept && !"cs" %in% type[[i]]) {
      new_lhs <- c(new_lhs, formula2str(fe_form, rm = 1))
      type[[i]] <- c(type[[i]], "")
    }
    if (length(new_lhs) > 1 && re_parts$mid[i] != "||") {
      id <- gsub("\\|", "", re_parts$mid[i])
      if (!nzchar(id)) {
        # ID is required to model coefficients as correlated 
        # if multiple types are provided within the same term
        id <- collapse(sample(0:9, 10, TRUE))
        re_parts$mid[i] <- paste0("|", id, "|")
      }
    }
    new_re_terms[[i]] <- paste0(new_lhs, re_parts$mid[i], re_parts$rhs[i])
    new_re_terms[[i]] <- new_re_terms[[i]][order(type[[i]])]
    type[[i]] <- sort(type[[i]])
  }
  re_terms <- unlist(new_re_terms)
  structure(re_terms, type = unlist(type))
}

get_re_terms <- function(x, formula = FALSE, brackets = TRUE) {
  # extract group-level terms from a formula of character vector
  # Args:
  #   x: formula or character vector
  #   formula: return a formula rather than a character string?
  #   brackets: include group-level terms in brackets?
  if (is.formula(x)) {
    x <- all_terms(x)
  }
  re_pos <- grepl("\\|", x)
  out <- x[re_pos]
  if (brackets && length(out)) {
    out <- paste0("(", out, ")")
  } 
  if (formula) {
    if (length(out)) {
      out <- formula(paste("~ 1", collapse("+", out)))
    } else {
      out <- ~ 1
    }
  }
  out
}

check_re_formula <- function(re_formula, formula) {
  # validate the re_formula argument
  # Args:
  #   re_formula: see predict.brmsfit for documentation
  #   formula: formula to match re_formula with
  # Returns:
  #   updated re_formula containing only terms existent in formula
  old_re_formula <- get_re_terms(formula, formula = TRUE)
  if (is.null(re_formula)) {
    re_formula <- old_re_formula
  } else if (SW(anyNA(re_formula))) {
    re_formula <- ~ 1
  } else {
    re_formula <- get_re_terms(as.formula(re_formula), formula = TRUE)
    new <- parse_bf(re_formula, check_response = FALSE)$dpars$mu$re
    old <- parse_bf(old_re_formula, check_response = FALSE)$dpars$mu$re
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

update_re_terms <- function(formula, re_formula) {
  # remove existing group-level terms in formula and
  # add valid group-level terms of re_formula
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
      paste0("+(", old_re_terms, ")"),
      paste0("(", old_re_terms, ")"),
      old_re_terms
    )
    new_formula <- rename(new_formula, rm_terms, "")
    if (grepl("~$", new_formula)) {
      # lhs only formulas are not allowed
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

get_re <- function(x, ...) {
  # extract group-level terms
  UseMethod("get_re")
}

#' @export
get_re.brmsterms <- function(x, all = TRUE, ...) {
  # get group-level information in a data.frame
  # Args:
  #   bterms: object of class brmsterms
  #   all: logical; include ranefs of nl and aux parameters?
  if (all) {
    re <- named_list(names(x$dpars))
    for (dp in names(re)) {
      re[[dp]] <- get_re(x$dpars[[dp]])
    }
    re <- do.call(rbind, re)
  } else {
    x$dpars[["mu"]]$nlpars <- NULL
    re <- get_re(x$dpars[["mu"]])
  }
  re
}

#' @export
get_re.mvbrmsterms <- function(x, ...) {
  do.call(rbind, lapply(x$terms, get_re, ...))
}

#' @export
get_re.btl <- function(x, ...) {
  stopifnot(is.data.frame(x$re))
  px <- check_prefix(x)
  re <- x$re
  re$resp <- rep(px$resp, nrow(re)) 
  re$dpar <- rep(px$dpar, nrow(re))
  re$nlpar <- rep(px$nlpar, nrow(re)) 
  re
}

#' @export
get_re.btnl <- function(x, ...) {
  re <- named_list(names(x$nlpars))
  for (nlp in names(re)) {
    re[[nlp]] <- get_re(x$nlpars[[nlp]])
  }
  do.call(rbind, re)
}

tidy_ranef <- function(bterms, data, all = TRUE, 
                       old_levels = NULL, old_sdata = NULL) {
  # combines helpful information on the group-level effects
  # Args:
  #   bterms: object of class brmsterms
  #   data: data.frame containing all model variables
  #   all: include REs of all parameters?
  #   old_levels: optional original levels of the grouping factors
  #   old_sdata: optional; see 'extract_old_standata'
  #     only used for category specific group-level effects
  # Returns: 
  #   A tidy data.frame with the following columns:
  #     id: ID of the group-level effect 
  #     group: name of the grouping factor
  #     gn: number of the grouping term within the respective formula
  #     coef: name of the group-level effect
  #     cn: number of the effect within the ID
  #     resp: name of the response variable
  #     dpar: name of the distributional parameter
  #     nlpar: name of the non-linear parameter
  #     cor: are correlations modeled for this effect?
  #     type: special effects type; can be 'sp' or 'cs'
  #     gcall: output of functions 'gr' or 'mm'
  #     form: formula used to compute the effects
  data <- combine_groups(data, get_groups(bterms))
  re <- get_re(bterms, all = all)
  ranef <- vector("list", nrow(re))
  used_ids <- new_ids <- NULL
  id_groups <- list()
  j <- 1
  for (i in seq_len(nrow(re))) {
    if (!nzchar(re$type[i])) {
      coef <- colnames(get_model_matrix(re$form[[i]], data)) 
    } else if (re$type[i] == "sp") {
      coef <- tidy_spef(re$form[[i]], data)$coef
    } else if (re$type[i] == "mmc") {
      coef <- rename(all_terms(re$form[[i]]))
    } else if (re$type[i] == "cs") {
      resp <- re$resp[i]
      if (!is.null(old_sdata)) {
        # extract ncat from the original data
        if (nzchar(resp)) {
          ncat <- old_sdata[[resp]][["ncat"]]
        } else {
          ncat <- old_sdata[["ncat"]]
        }
      } else {
        # infer ncat from the current data
        if (nzchar(resp)) {
          respform <- bterms$terms[[resp]]$respform
        } else {
          respform <- bterms$respform
        }
        Y <- as.numeric(model.response(model.frame(respform, data)))
        ncat <- max(Y) - min(Y) + 1
      }
      indices <- paste0("[", seq_len(ncat - 1), "]")
      coef <- colnames(get_model_matrix(re$form[[i]], data = data))
      coef <- as.vector(t(outer(coef, indices, paste0)))
    }
    avoid_dpars(coef, bterms = bterms)
    rdat <- data.frame(
      id = re$id[[i]],
      group = re$group[[i]],
      gn = re$gn[[i]],
      gtype = re$gtype[[i]],
      coef = coef, cn = NA,
      resp = re$resp[[i]],
      dpar = re$dpar[[i]],
      nlpar = re$nlpar[[i]],
      cor = re$cor[[i]],
      type = re$type[[i]],
      by = re$gcall[[i]]$by,
      stringsAsFactors = FALSE
    )
    bylevels <- NULL
    if (nzchar(rdat$by[1])) {
      bylevels <- levels(factor(get(rdat$by[1], data)))
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
  dup <- duplicated(ranef[, c("group", "coef", vars_prefix())])
  if (any(dup)) {
    dr <- ranef[which(dup)[1], ]
    stop2(
      "Duplicated group-level effects are not allowed.\n",
      "Occured for effect '", dr$coef, "' of group '", dr$group, "'."
    )
  }
  if (nrow(ranef)) {
    for (id in unique(ranef$id)) {
      ranef$cn[ranef$id == id] <- seq_len(sum(ranef$id == id))
    }
    if (is.null(old_levels)) {
      rsub <- ranef[!duplicated(ranef$group), ]
      levels <- named_list(rsub$group)
      for (i in seq_along(levels)) {
        # combine levels of all grouping factors within one grouping term
        levels[[i]] <- unique(ulapply(
          rsub$gcall[[i]]$groups, 
          function(g) levels(factor(get(g, data)))
        ))
        by <- rsub$by[i]
        if (nzchar(by)) {
          # store information of corresponding by levels
          stopifnot(!nzchar(rsub$type[i]))
          bylevels <- rsub$bylevels[[i]]
          byvar <- get(by, data)
          g <- rsub$gcall[[i]]$groups
          gvar <- get(g, data)
          not_dupl <- !duplicated(data.frame(gvar, byvar))
          if (sum(not_dupl) > length(levels[[i]])) {
            stop2(
              "Some levels of '", g, "' correspond ", 
              "to multiple levels of '", by, "'. "
            )
          }
          by_per_level <- bylevels[match(byvar[not_dupl], bylevels)]
          attr(levels[[i]], "by") <- by_per_level
        }
      }
      attr(ranef, "levels") <- levels 
    } else {
      # for newdata numeration has to depend on the original levels
      attr(ranef, "levels") <- old_levels
    }
  }
  structure(ranef, class = c("ranef_frame", "data.frame"))
}

empty_ranef <- function() {
  structure(
    data.frame(
      id = numeric(0), group = character(0), gn = numeric(0),
      coef = character(0), cn = numeric(0), resp = character(0),
      dpar = character(0), nlpar = character(0), cor = logical(0), 
      type = character(0), form = character(0), 
      stringsAsFactors = FALSE
    ),
    class = c("ranef_frame", "data.frame")
  )
}

is.ranef_frame <- function(x) {
  inherits(x, "ranef_frame")
}

get_group_vars <- function(x, ...) {
  # extract names of grouping variables
  UseMethod("get_group_vars") 
}

#' @export
get_group_vars.brmsfit <- function(x, ...) {
  bterms <- parse_bf(x$formula)
  unique(c(
    get_group_vars(x$ranef),
    get_group_vars(tidy_meef(bterms, x$data)),
    get_autocor_vars(x, var = "group")
  ))
}

#' @export
get_group_vars.ranef_frame <- function(x, ...) {
  unique(ulapply(x$gcall, "[[", "groups"))
}

#' @export
get_group_vars.meef_frame <- function(x, ...) {
  if (nrow(x)) {
    out <- unique(x$grname[!is.na(x$grname)])
  } else {
    out <- NULL
  }
  out
}
