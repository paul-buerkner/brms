# This file contains functions dealing with the extended 
# formula syntax to specify special effects terms

vars_keep_na <- function(x, ...) {
  # find variable names for which to keep NAs
  UseMethod("vars_keep_na")
}

#' @export
vars_keep_na.mvbrmsterms <- function(x, ...) {
  resps <- get_element(x, "respform")
  resps <- ulapply(resps, parse_resp, check_names = FALSE)
  out <- lapply(x$terms, vars_keep_na, responses = resps, ...)
  vars_mi <- unique(ulapply(out, attr, "vars_mi"))
  out <- unique(unlist(out))
  miss_mi <- setdiff(vars_mi, out)
  if (length(miss_mi)) {
    stop2(
      "Response models of variables in 'mi' terms require " ,
      "specification of the addition argument 'mi'. See ?mi. ", 
      "Error occured for ", collapse_comma(miss_mi), "."
    )
  }
  out
}

#' @export
vars_keep_na.brmsterms <- function(x, responses = NULL, ...) {
  if (is.formula(x$adforms$mi)) {
    mi_respvars <- parse_resp(x$respform, check_names = FALSE)
    mi_advars <- all.vars(x$adforms$mi)
    out <- unique(c(mi_respvars, mi_advars))
  } else {
    out <- character(0)
  }
  uni_mi <- ulapply(get_effect(x, "sp"), attr, "uni_mi")
  if (length(uni_mi)) {
    vars_mi <- gsub("(^mi\\()|(\\)$)", "", uni_mi)
    miss_mi <- setdiff(vars_mi, responses)
    if (length(miss_mi)) {
      stop2(
        "Variables in 'mi' terms should also be specified " ,
        "as response variables in the model. See ?mi. ", 
        "Error occured for ", collapse_comma(miss_mi), "."
      )
    }
    attr(out, "vars_mi") <- vars_mi
  }
  out
}

get_uni_me <- function(x) {
  # extract unique names of noise-free terms 
  uni_me <- ulapply(get_effect(x, "sp"), attr, "uni_me")
  if (!length(uni_me)) {
    return(NULL)
  }
  all_vars <- all.vars(parse(text = uni_me))
  elist <- named_list(all_vars, values = NA_real_)
  xname <- ulapply(uni_me, function(x) attr(eval2(x, elist), "xname"))
  df <- data.frame(xname, uni_me)
  df <- df[!duplicated(df), ]
  xdupl <- df$xname[duplicated(df$xname)]
  if (length(xdupl)) {
    calls <- df$uni_me[df$xname == xdupl[1]]
    stop2(
      "Variable '", xdupl[1], "' is used in different calls to 'me'.\n",
      "Associated calls are: ", collapse_comma(calls)
    )
  }
  unique(uni_me)
}

tidy_meef <- function(bterms, data, old_levels = NULL) {
  # save all me-terms within a tidy data.frame
  uni_me <- get_uni_me(bterms)
  if (length(uni_me)) {
    out <- data.frame(
      term = uni_me, xname = "", grname = "", 
      stringsAsFactors = FALSE
    )
    levels <- list("vector", nrow(out))
    for (i in seq_rows(out)) {
      att <- attributes(eval2(out$term[i], data))
      out$xname[i] <- att$xname
      if (isTRUE(nzchar(att$grname))) {
        out$grname[i] <- att$grname
      }
      if (is.null(old_levels)) {
        levels[[i]] <- levels(factor(att$gr))
      } else {
        levels[[i]] <- old_levels[[att$grname]]
      }
    }
    out$coef <- rename(paste0("me", out$xname))
    out$cor <- isTRUE(bterms$mecor)
    names(levels) <- out$grname
    levels <- levels[!is.na(names(levels))]
    if (length(levels)) {
      levels <- levels[!duplicated(names(levels))]
      attr(out, "levels") <- levels
    }
  } else {
    out <- data.frame(
      terms = character(0), xname = character(0),
      grname = character(0), cor = logical(0),
      stringsAsFactors = FALSE
    )
  }
  structure(out, class = c("meef_frame", "data.frame"))
}

is.meef_frame <- function(x) {
  inherits(x, "meef_frame")
}

default_mecor <- function(mecor = NULL) {
  # handle default of correlations between 'me' terms
  if (is.null(mecor)) TRUE else as_one_logical(mecor)
}

get_sp_vars <- function(x, type) {
  # find names of all variables used in a special effects type
  sp_terms <- ulapply(get_effect(x, "sp"), all_terms)
  all.vars(str2formula(get_matches_expr(regex_sp(type), sp_terms)))
}

tidy_spef <- function(x, data) {
  # get labels of special effects terms
  # Args:
  #   x: either a formula or a list containing an element "sp"
  #   data: data frame containing the monotonic variables
  if (is.formula(x)) {
    x <- parse_bf(x, check_response = FALSE)$dpars$mu
  }
  form <- x[["sp"]]
  if (!is.formula(form)) {
    return(empty_data_frame())
  }
  mm <- get_model_matrix(form, data, rename = FALSE)
  out <- data.frame(term = rm_wsp(colnames(mm)), stringsAsFactors = FALSE)
  out$coef <- rename(out$term)
  call_cols <- paste0("call_", c("mo", "me", "mi"))
  for (col in c(call_cols, "call_prod", "vars_mi", "Imo")) {
    out[[col]] <- vector("list", nrow(out))
  }
  kmo <- 0
  terms_split <- strsplit(out$term, ":")
  for (i in seq_rows(out)) {
    # prepare mo terms
    take_mo <- grepl_expr(regex_sp("mo"), terms_split[[i]])
    if (sum(take_mo)) {
      out$call_mo[[i]] <- terms_split[[i]][take_mo]
      nmo <- length(out$call_mo[[i]])
      out$Imo[[i]] <- (kmo + 1):(kmo + nmo)
      kmo <- kmo + nmo
      for (j in seq_along(out$call_mo[[i]])) {
        mo_term <- out$call_mo[[i]][[j]]
        mo_match <- get_matches_expr(regex_sp("mo"), mo_term)
        if (length(mo_match) > 1L || nchar(mo_match) < nchar(mo_term)) {
          stop2("The monotonic term '",  mo_term, "' is invalid.")
        }
      }
    }
    # prepare me terms
    take_me <- grepl_expr(regex_sp("me"), terms_split[[i]])
    if (sum(take_me)) {
      out$call_me[[i]] <- terms_split[[i]][take_me]
      # remove 'I' (identity) function calls that 
      # were used solely to separate formula terms
      out$call_me[[i]] <- gsub("^I\\(", "(", out$call_me[[i]])
    }
    # prepare mi terms 
    take_mi <- grepl_expr(regex_sp("mi"), terms_split[[i]])
    if (sum(take_mi)) {
      out$call_mi[[i]] <- terms_split[[i]][take_mi]
      out$call_mi[[i]] <- gsub("^I\\(", "(", out$call_mi[[i]])
      out$vars_mi[[i]] <- get_matches_expr(regex_sp("mi"), out$call_mi[[i]])
      out$vars_mi[[i]] <- gsub("(^mi\\()|(\\)$)", "", out$vars_mi[[i]])
      # do it like parse_resp to ensure correct matching
      out$vars_mi[[i]] <- gsub("\\.|_", "", make.names(out$vars_mi[[i]]))
    }
    out$call_prod[[i]] <- paste0(unlist(out[i, call_cols]), collapse = " * ")
  }
  not_one <- apply(mm, 2, function(x) any(x != 1))
  out$Ic <- ulapply(seq_along(not_one), function(i) sum(not_one[1:i]))
  out
}

get_simo_labels <- function(spef) {
  # extract names of monotonic simplex parameters 
  # Args:
  #   spef: output of tidy_spef
  fun <- function(i) paste0(spef$coef[i], seq_along(spef$Imo[[i]]))
  ulapply(which(lengths(spef$Imo) > 0), fun)
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

get_me_groups <- function(x) {
  # get names of me grouping variables
  uni_me <- get_uni_me(x)
  out <- lapply(uni_me, eval_NA) 
  out <- ulapply(out, attr, "grname")
  out[nzchar(out)]
}
