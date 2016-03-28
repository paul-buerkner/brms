isNULL <- function(x) {
  # check if an object is NULL
  is.null(x) || ifelse(is.vector(x), all(sapply(x, is.null)), FALSE)
}

rmNULL <- function(x) {
  # recursively removes NULL entries from an object
  x <- Filter(Negate(isNULL), x)
  lapply(x, function(x) if (is.list(x)) rmNULL(x) else x)
}

isFALSE <- function(x) {
  identical(FALSE, x)
}

is_equal <- function(x, y, ...) {
  isTRUE(all.equal(x, y, ...))
}

rmNum <- function(x) {
  # remove all numeric elements from an object
  x[sapply(x, Negate(is.numeric))]
}

rmMatch <- function(x, ...) {
  # remove all elements in x that also appear in ... while keeping all attributes
  att <- attributes(x)
  keep <- which(!(x %in% c(...)))
  x <- x[keep]
  attributes(x) <- att
  attr(x, "match.length") <- att$match.length[keep] 
  x
}

subset_attr <- function(x, y) {
  # take a subset of vector, list, etc. 
  # while keeping all attributes except for names
  att <- attributes(x)
  x <- x[y]
  att[["names"]] <- names(x)
  attributes(x) <- att
  x
} 

is.wholenumber <- function(x, tol = .Machine$double.eps) {  
  # check if x is a whole number (integer)
  if (!is.numeric(x)) {
    return(FALSE)
  } else {
    return(abs(x - round(x)) < tol)
  }
} 

ulapply <- function(X, FUN, ...) {
  # short for unlist(lapply(.))
  unlist(lapply(X = X, FUN = FUN, ...))
}

as_matrix <- function(x, ...) {
  # wrapper around as.matrix that can handle NULL
  if (is.null(x)) {
    NULL
  } else {
    as.matrix(x, ...) 
  }
}

lc <- function(l, x) {
  # append x to l
  l[[length(l) + 1]] <- x
  l
}

collapse <- function(..., sep = "") {
  # wrapper for paste with collapse
  paste(..., sep = sep, collapse = "")
}

collapse_lists <- function(ls) {
  # collapse strings having the same name in different lists
  #
  # Args:
  #  ls: a list of named lists
  # 
  # Returns:
  #  a named list containg the collapsed strings
  elements <- unique(unlist(lapply(ls, names)))
  out <- do.call(mapply, c(FUN = collapse, lapply(ls, "[", elements), 
                           SIMPLIFY = FALSE))
  names(out) <- elements
  out
}

nlist <- function(...) {
  # create a named list using object names
  m <- match.call()
  dots <- list(...)
  no_names <- is.null(names(dots))
  has_name <- if (no_names) FALSE 
              else nzchar(names(dots))
  if (all(has_name)) return(dots)
  nms <- as.character(m)[-1]
  if (no_names) {
    names(dots) <- nms
  } else {
    names(dots)[!has_name] <- nms[!has_name]
  }
  dots
}

get_arg <- function(x, ...) {
  # find first occurrence of x in ... objects
  # Args:
  #  x: The name of the required element
  #  ...: R objects that may contain x
  dots <- list(...)
  i <- 1
  out <- NULL
  while(i <= length(dots) && is.null(out)) {
    if (!is.null(dots[[i]][[x]])) {
      out <- dots[[i]][[x]]
    } else i <- i + 1
  }
  out
}

rhs <- function(x) {
  # return the righthand side of a formula
  attri <- attributes(x)
  x <- as.formula(x)
  x <- if (length(x) == 3) x[-2] else x
  do.call(structure, c(list(x), attri))
}

lhs <- function(x) {
  # return the lefthand side of a formula
  x <- as.formula(x)
  if (length(x) == 3) update(x, . ~ 1) else NULL
}

SW <- function(expr) {
  # just a short form for suppressWarnings
  base::suppressWarnings(expr)
}

get_matches <- function(pattern, text, simplify = TRUE, ...) {
  # get pattern matches in text as vector
  x <- regmatches(text, gregexpr(pattern, text, ...))
  if (simplify) x <- unlist(x)
  x
}

logit <- function(p) {
  # compute the logit
  log(p / (1 - p))
}

inv_logit <- function(x) { 
  # compute the inverse of logit
  1 / (1 + exp(-x))
}

incgamma <- function(x, a) {
  # incomplete gamma funcion
  pgamma(x, shape = a) * gamma(a)
}

wsp <- function(x, nsp = 1) {
  # add leading and trailing whitespaces
  # Args:
  #   x: object accepted by paste
  #   nsp: number of whitespaces to add
  sp <- collapse(rep(" ", nsp))
  if (length(x)) paste0(sp, x, sp)
  else NULL
}

limit_chars <- function(x, chars = NULL, lsuffix = 4) {
  # limit the number of characters of a vector
  # Args:
  #   x: a character vector
  #   chars: maximum number of characters to show
  #   lsuffix: number of characters to keep 
  #            at the end of the strings
  stopifnot(is.character(x))
  if (!is.null(chars)) {
    chars_x <- nchar(x) - lsuffix
    suffix <- substr(x, chars_x + 1, chars_x + lsuffix)
    x <- substr(x, 1, chars_x)
    x <- ifelse(chars_x <= chars, x, paste0(substr(x, 1, chars - 3), "..."))
    x <- paste0(x, suffix)
  }
  x
}

use_alias <- function(arg, alias = NULL, warn = TRUE) {
  # ensure that deprecated arguments still work
  # Args:
  #   arg: input to the new argument
  #   alias: input to the deprecated argument
  arg_name <- Reduce(paste, deparse(substitute(arg)))
  alias_name <- Reduce(paste, deparse(substitute(alias)))
  if (!is.null(alias)) {
    arg <- alias
    if (substr(alias_name, 1, 5) == "dots$") {
      alias_name <- substr(alias_name, 6, nchar(alias_name))
    }
    if (warn) {
      warning(paste0("Argument '", alias_name, "' is deprecated. ", 
                     "Please use argument '", arg_name, "' instead."), 
              call. = FALSE)
    }
  }
  arg
}

.addition <- function(formula, data = NULL) {
  # computes data for addition arguments
  if (!is.formula(formula))
    formula <- as.formula(formula)
  eval(formula[[2]], data, environment(formula))
}

.se <- function(x) {
  # standard errors for meta-analysis
  if (!is.numeric(x)) 
    stop("SEs must be numeric")
  if (min(x) < 0) 
    stop("standard errors must be non-negative", call. = FALSE)
  x  
}

.weights <- function(x) {
  # weights to be applied on any model
  if (!is.numeric(x)) 
    stop("weights must be numeric")
  if (min(x) < 0) 
    stop("weights must be non-negative", call. = FALSE)
  x
}

.disp <- function(x) {
  # dispersion factors
  if (!is.numeric(x)) 
    stop("dispersion factors must be numeric")
  if (min(x) < 0) 
    stop("dispersion factors must be non-negative", call. = FALSE)
  x  
}

.trials <- function(x) {
  # trials for binomial models
  if (any(!is.wholenumber(x) || x < 1))
    stop("number of trials must be positive integers", call. = FALSE)
  x
}

.cat <- function(x) {
  # number of categories for categorical and ordinal models
  if (any(!is.wholenumber(x) || x < 1))
    stop("number of categories must be positive integers", call. = FALSE)
  x
}

.cens <- function(x) {
  # indicator for censoring
  if (is.factor(x)) x <- as.character(x)
  cens <- unname(sapply(x, function(x) {
    if (grepl(paste0("^", x), "right") || isTRUE(x)) x <- 1
    else if (grepl(paste0("^", x), "none") || isFALSE(x)) x <- 0
    else if (grepl(paste0("^", x), "left")) x <- -1
    else x
  }))
  if (!all(unique(cens) %in% c(-1:1)))
    stop (paste0("Invalid censoring data. Accepted values are ", 
                 "'left', 'none', and 'right' \n(abbreviations are allowed) ", 
                 "or -1, 0, and 1. TRUE and FALSE are also accepted \n",
                 "and refer to 'right' and 'none' respectively."),
          call. = FALSE)
  cens
}

.trunc <- function(lb = -Inf, ub = Inf) {
  lb <- as.numeric(lb)
  ub <- as.numeric(ub)
  if (length(lb) != 1 || length(ub) != 1) {
    stop("invalid truncation values", call. = FALSE)
  }
  nlist(lb, ub)
}

# startup messages for brms
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0(
    "Loading 'brms' package (version ", utils::packageVersion("brms"), "). ",
    "Useful instructions \n", 
    "can be found by typing help('brms'). A more detailed introduction \n", 
    "to the package is available through vignette('brms')."))
}