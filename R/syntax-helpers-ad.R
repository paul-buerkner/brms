# has a certain addition term a certain variable?
# @param x list with potentail $adforms elements
# @param data data passed by the user
# @return TRUE or FALSE
has_advar <- function(x, ad, var) {
  as <- as_one_character(ad)
  var <- as_one_character(var)
  ad <- x$adforms[[ad]]
  if (is.null(ad)) {
    return(FALSE)
  }
  ad <- eval_rhs(ad)
  isTRUE(!is.null(ad$vars[[var]]) && ad$vars[[var]] != "NA")
}

# get values of a variable used in an addition term
# @param x list with potentail $adforms elements
# @param data data passed by the user
# @return a vector of threshold groups or NULL
get_advalues <- function(x, ad, var, data) {
  as <- as_one_character(ad)
  var <- as_one_character(var)
  if (!has_advar(x, ad, var)) {
    return(NULL)
  }
  ad <- eval_rhs(x$adforms[[ad]])
  eval2(ad$vars[[var]], data)
}

# coerce censored values into the right format
# @param x vector of censoring indicators
# @return transformed vector of censoring indicators
prepare_cens <- function(x) {
  .prepare_cens <- function(x) {  
    stopifnot(length(x) == 1L)
    regx <- paste0("^", x)
    if (grepl(regx, "left")) {
      x <- -1
    } else if (grepl(regx, "none") || isFALSE(x)) {
      x <- 0
    } else if (grepl(regx, "right") || isTRUE(x)) {
      x <- 1
    } else if (grepl(regx, "interval")) {
      x <- 2
    }
    return(x)
  }
  x <- unname(x)
  if (is.factor(x)) {
    x <- as.character(x)
  }
  ulapply(x, .prepare_cens)
}

# extract information on censoring of the response variable
# @param x a brmsfit object
# @param resp optional names of response variables for which to extract values
# @return vector of censoring indicators or NULL in case of no censoring
get_cens <- function(x, resp = NULL, newdata = NULL) {
  stopifnot(is.brmsfit(x))
  resp <- validate_resp(resp, x, multiple = FALSE)
  bterms <- parse_bf(x$formula)
  if (!is.null(resp)) {
    bterms <- bterms$terms[[resp]]
  }
  if (is.null(newdata)) {
    newdata <- model.frame(x)
  }
  out <- NULL
  if (is.formula(bterms$adforms$cens)) {
    cens <- eval_rhs(bterms$adforms$cens)
    out <- eval2(cens$vars$cens, newdata)
    out <- prepare_cens(out)
  }
  out
}

# extract truncation boundaries
trunc_bounds <- function(x, ...) {
  UseMethod("trunc_bounds")
}

# @return a named list with one element per response variable
#' @export
trunc_bounds.mvbrmsterms <- function(x, ...) {
  lapply(x$terms, trunc_bounds, ...)
}

# @param data data.frame containing the truncation variables
# @param incl_family include the family in the derivation of the bounds?
# @param stan return bounds in form of Stan syntax?
# @return a list with elements 'lb' and 'ub'
#' @export
trunc_bounds.brmsterms <- function(x, data = NULL, incl_family = FALSE, 
                                   stan = FALSE, ...) {
  if (is.formula(x$adforms$trunc)) {
    trunc <- eval_rhs(x$adforms$trunc)
  } else {
    trunc <- resp_trunc()
  }
  out <- list(
    lb = eval2(trunc$vars$lb, data),
    ub = eval2(trunc$vars$ub, data)
  )
  if (incl_family) {
    family_bounds <- family_bounds(x)
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
