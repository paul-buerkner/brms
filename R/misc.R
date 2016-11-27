p <- function(x, i = NULL, row = TRUE) {
  # flexible indexing of vector and matrix type objects
  # Args:
  #   x: an R object typically a vector or matrix
  #   i: optional index; if NULL, x is returned unchanged
  #   row: indicating if rows or cols should be indexed
  #        only relevant if x has two dimensions
  if (!length(i)) {
    x
  } else if (length(dim(x)) == 2L) {
    if (row) {
      x[i, , drop = FALSE]
    } else {
      x[, i, drop = FALSE]
    }
  } else {
    x[i]
  }
}

isNULL <- function(x) {
  # check if an object is NULL
  is.null(x) || ifelse(is.vector(x), all(sapply(x, is.null)), FALSE)
}

rmNULL <- function(x, recursive = TRUE) {
  # recursively removes NULL entries from an object
  x <- Filter(Negate(isNULL), x)
  if (recursive) {
    x <- lapply(x, function(x) if (is.list(x)) rmNULL(x) else x)
  }
  x
}

first_not_null <- function(...) {
  # find the first argument that is not NULL
  dots <- list(...)
  out <- NULL
  i <- 1L
  while(isNULL(out) && i <= length(dots)) {
    if (!isNULL(dots[[i]])) {
      out <- dots[[i]]
    }
    i <- i + 1L
  }
  out
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
  # remove all elements in x that also appear in ... 
  # while keeping all attributes
  att <- attributes(x)
  keep <- which(!(x %in% c(...)))
  x <- x[keep]
  attributes(x) <- att
  attr(x, "match.length") <- att$match.length[keep] 
  x
}

rm_attr <- function(x, attr) {
  # remove certain attributes
  attributes(x)[attr] <- NULL
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
  # Args:
  #  ls: a list of named lists
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

named_list <- function(names, values = NULL) {
  # initialize a named list
  # Args:
  #   names: names of the elements
  #   values: values of the elements
  if (!is.null(values)) {
    if (length(values) <= 1L) {
      values <- replicate(length(names), values)
    }
    values <- as.list(values)
    stopifnot(length(values) == length(names))
  } else {
    values <- vector("list", length(names))
  }
  setNames(values, names)
} 

deparse_no_string <- function(x) {
  # deparse x if it is no string
  if (!is.character(x)) {
    x <- deparse(x)
  } 
  x
}

eval2 <- function(text, ...) {
  # evaluate a string
  eval(parse(text = text), ...)
}

stop2 <- function(...) {
  stop(..., call. = FALSE)
}

warning2 <- function(...) {
  warning(..., call. = FALSE)
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
    } else {
      i <- i + 1
    }
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
  if (length(x) == 3L) update(x, . ~ 1) else NULL
}

is.formula <- function(x) {
  is(x, "formula")
}

SW <- function(expr) {
  # just a short form for suppressWarnings
  base::suppressWarnings(expr)
}

get_matches <- function(pattern, text, simplify = TRUE, ...) {
  # get pattern matches in text as vector
  x <- regmatches(text, gregexpr(pattern, text, ...))
  if (simplify) {
    x <- unlist(x)
  }
  x
}

usc <- function(x, pos = c("prefix", "suffix")) {
  # add an underscore to non-empty character strings
  # Args:
  #   x: a character vector
  #   pos: position of the underscore
  pos <- match.arg(pos)
  x <- as.character(x)
  if (pos == "prefix") {
    x <- ifelse(nzchar(x), paste0("_", x), "")
  } else {
    x <- ifelse(nzchar(x), paste0(x, "_"), "")
  }
  x
}

logit <- function(p) {
  # logit link
  log(p / (1 - p))
}

inv_logit <- function(x) { 
  # inverse of logit link
  1 / (1 + exp(-x))
}

cloglog <- function(x) {
  # cloglog link
  log(-log(1-x))
}

inv_cloglog <- function(x) {
  # inverse of the cloglog link
  1 - exp(-exp(x))
}

Phi <- function(x) {
  pnorm(x)
}

incgamma <- function(x, a) {
  # incomplete gamma funcion
  pgamma(x, shape = a) * gamma(a)
}

square <- function(x) {
  x^2
}

cbrt <- function(x) {
  x^(1/3)
}

exp2 <- function(x) {
  2^x
}

pow <- function(x, y) {
  x^y
}

inv <- function(x) {
  1/x
}

inv_sqrt <- function(x) {
  1/sqrt(x)
}

inv_square <- function(x) {
  1/x^2
}

hypot <- function(x, y) {
  stopifnot(all(x >= 0))
  stopifnot(all(y >= 0))
  sqrt(x^2 + y^2)
}

log1p <- function(x) {
  log(1 + x)
}

log1m <- function(x) {
  log(1 - x)
}

expm1 <- function(x) {
  exp(x) - 1
}

multiply_log <- function(x, y) {
  ifelse(x == y & x == 0, 0, x * log(y))
}

log1p_exp <- function(x) {
  log(1 + exp(x))
}

log1m_exp <- function(x) {
  ifelse(x < 0, log(1 - exp(x)), NaN)
}

log_diff_exp <- function(x, y) {
  stopifnot(length(x) == length(y))
  ifelse(x > y, log(exp(x) - exp(y)), NaN)
}

log_sum_exp <- function(x, y) {
  max <- max(x, y)
  max + log(exp(x - max) + exp(y - max))
}

log_inv_logit <- function(x) {
  log(inv_logit(x))
}

log1m_inv_logit <- function(x) {
  log(1 - inv_logit(x))
}

fabs <- function(x) {
  abs(x)
}

wsp <- function(x, nsp = 1) {
  # add leading and trailing whitespaces
  # Args:
  #   x: object accepted by paste
  #   nsp: number of whitespaces to add
  sp <- collapse(rep(" ", nsp))
  if (length(x)) {
    paste0(sp, x, sp)
  } else {
    NULL
  } 
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

use_alias <- function(arg, alias = NULL, default = NULL,
                      warn = TRUE) {
  # ensure that deprecated arguments still work
  # Args:
  #   arg: input to the new argument
  #   alias: input to the deprecated argument
  #   default: the default value of alias
  #   warn: should a warning be printed if alias is specified?
  arg_name <- Reduce(paste, deparse(substitute(arg)))
  alias_name <- Reduce(paste, deparse(substitute(alias)))
  if (!is_equal(alias, default)) {
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
  formula <- as.formula(formula)
  eval(formula[[2]], data, environment(formula))
}

.se <- function(x) {
  # standard errors for meta-analysis
  if (!is.numeric(x)) 
    stop2("Standard errors must be numeric.")
  if (min(x) < 0) 
    stop2("Standard errors must be non-negative.")
  x  
}

.weights <- function(x) {
  # weights to be applied on any model
  if (!is.numeric(x)) 
    stop2("Weights must be numeric.")
  if (min(x) < 0) 
    stop2("Weights must be non-negative.")
  x
}

.disp <- function(x) {
  # dispersion factors
  if (!is.numeric(x)) 
    stop2("Dispersion factors must be numeric.")
  if (min(x) < 0) 
    stop2("Dispersion factors must be non-negative.")
  x  
}

.dec <- function(x) {
  # decisions for the wiener diffusion model
  if (is.character(x) || is.factor(x)) {
    x <- ifelse(x == "lower", 0, ifelse(x == "upper", 1, x))
    if (!is.numeric(x)) {
      stop2("Decisions should be 'lower' or 'upper' ", 
            "when supplied as characters or factors.")
    }
  } else {
    x <- as.numeric(as.logical(x))
  }
  x
}

.trials <- function(x) {
  # trials for binomial models
  if (any(!is.wholenumber(x) || x < 1))
    stop2("Number of trials must be positive integers.")
  x
}

.cat <- function(x) {
  # number of categories for categorical and ordinal models
  if (any(!is.wholenumber(x) || x < 1))
    stop2("Number of categories must be positive integers.")
  x
}

.cens <- function(x, y2 = NULL) {
  # indicator for censoring
  if (is.factor(x)) {
    x <- as.character(x)
  }
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
    x
  }
  cens <- unname(ulapply(x, .prepare_cens))
  if (!all(is.wholenumber(cens) & cens %in% -1:2)) {
    stop2("Invalid censoring data. Accepted values are ", 
          "'left', 'none', 'right', and 'interval'\n",
          "(abbreviations are allowed) or -1, 0, 1, and 2.\n",
          "TRUE and FALSE are also accepted ",
          "and refer to 'right' and 'none' respectively.")
  }
  if (any(cens %in% 2)) {
    if (length(y2) != length(cens)) {
      stop2("Argument 'y2' is required for interval censored data.")
    }
    attr(cens, "y2") <- unname(y2)
  }
  cens
}

.trunc <- function(lb = -Inf, ub = Inf) {
  lb <- as.numeric(lb)
  ub <- as.numeric(ub)
  if (any(lb >= ub)) {
    stop2("Truncation bounds are invalid: lb >= ub")
  }
  nlist(lb, ub)
}

cse <- function(...) {
  stop2("Inappropriate use of function 'cse'")
}

monotonic <- function(...) {
  stop2("Inappropriate use of function 'monotonic'")
}

mono <- function(...) {
  # abbreviation of monotonic
  stop2("Inappropriate use of function 'monotonic'")
}

monotonous <- function(...) {
  # abbreviation of monotonic
  stop2("Please use function 'monotonic' instead.")
}

# startup messages for brms
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0(
    "Loading 'brms' package (version ", utils::packageVersion("brms"), "). ",
    "Useful instructions \n", 
    "can be found by typing help('brms'). A more detailed introduction \n", 
    "to the package is available through vignette('brms')."))
}