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

keep_attr <- function(x, y) {
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

get_matches <- function(pattern, text, ...) {
  # get pattern matches in text as vector
  unlist(regmatches(text, gregexpr(pattern, text, ...)))
}

logit <- function(p) {
  # compute the logit
  log(p / (1 - p))
}

ilogit <- function(x) { 
  # compute the inverse of logit
  exp(x) / (1 + exp(x))
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
  sp <- paste(rep(" ", nsp), collapse = "")
  paste0(sp, x, sp)
}

is.formula <- function(x, or = TRUE) {
  # checks if x is formula (or list of formulas)
  #
  # Returns:
  #   x: a formula or a list of formulas
  #   or: logical; indicates if any element must be a formula (or = TRUE) 
  #       or if all elements must be formulas
  if (!is.list(x)) x <- list(x)
  out <- sapply(x, function(y) is(y, "formula"))
  if (or) {
    out <- any(out)
  } else out <- all(out)
  out
}

formula2string <- function(formula, rm = c(0, 0)) {
  # converts formula to string
  #
  # Args:
  #   formula: a model formula
  #   rm: a vector of to elements indicating how many characters 
  #       should be removed at the beginning
  #       and end of the string respectively
  #
  # Returns:
  #    the formula as string 
  if (!is.formula(formula))
    stop(paste(deparse(substitute(formula)), "must be of class formula"))
  if (is.na(rm[2])) rm[2] <- 0
  x <- gsub(" ","", Reduce(paste, deparse(formula)))
  x <- substr(x, 1 + rm[1], nchar(x) - rm[2])
  x
} 

is.linear <- function(family) {
  # indicate if family is for a linear model
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("gaussian", "student", "cauchy")
}

is.lognormal <- function(family, link = "identity", nresp = 1) {
  # indicate transformation to lognormal model
  # Args:
  #   link: A character string
  #   nresp: number of response variables
  if (is(family, "family")) {
    link <- family$link
    family <- family$family
  }
  family == "gaussian" && link == "log" && nresp == 1
}

is.binary <- function(family) {
  # indicate if family is bernoulli or binomial
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("binomial", "bernoulli")
}

is.ordinal <- function(family) {
  # indicate if family is for an ordinal model
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("cumulative", "cratio", "sratio", "acat") 
}

is.categorical <- function(family) {
  if (is(family, "family")) {
    family <- family$family
  }
  family == "categorical" 
}

is.skewed <- function(family) {
  # indicate if family is for model with postive skewed response
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("gamma", "weibull", "exponential")
}

is.count <- function(family) {
  # indicate if family is for a count model
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("poisson", "negbinomial", "geometric")
}

is.hurdle <- function(family) {
  # indicate if family is for a hurdle model
  if (is(family, "family")) {
    family <- family$family
  }
  # zi_beta is technically a hurdle model
  family %in% c("hurdle_poisson", "hurdle_negbinomial", "hurdle_gamma",
                "zero_inflated_beta")
}

is.zero_inflated <- function(family) {
  # indicate if family is for a zero inflated model
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("zero_inflated_poisson", "zero_inflated_negbinomial",
                "zero_inflated_binomial")
}

is.2PL <- function(family) {
  if (!is(family, "brmsfamily")) {
    out <- FALSE
  } else {
    out <- family$family == "bernoulli" && identical(family$type, "2PL")
  }
  out
}

is.forked <- function(family) {
  # indicate if family has two separate model parts
  is.hurdle(family) || is.zero_inflated(family) || is.2PL(family)
}

use_real <- function(family) {
  # indicate if family uses real responses
  if (is(family, "family")) {
    family <- family$family
  }
  is.linear(family) || is.skewed(family) || 
    family %in% c("inverse.gaussian", "beta", "zero_inflated_beta", 
                  "hurdle_gamma")
}

use_int <- function(family) {
  # indicate if family uses integer responses
  if (is(family, "family")) {
    family <- family$family
  }
  is.binary(family) || has_cat(family) || 
    is.count(family) || is.zero_inflated(family) || 
    family %in% c("hurdle_poisson", "hurdle_negbinomial")
}

has_trials <- function(family) {
  # indicate if family makes use of argument trials
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("binomial", "zero_inflated_binomial")
}

has_cat <- function(family) {
  # indicate if family makes use of argument cat
  if (is(family, "family")) {
    family <- family$family
  }
  family == "categorical" || is.ordinal(family)
}

has_shape <- function(family) {
  # indicate if family needs a shape parameter
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("gamma", "weibull", "inverse.gaussian", 
                "negbinomial", "hurdle_negbinomial", 
                "hurdle_gamma", "zero_inflated_negbinomial")
}

has_sigma <- function(family, autocor = cor_arma(), se = FALSE,
                      is_multi = FALSE) {
  # indicate if the model needs a sigma parameter
  # Args:
  #  family: model family
  #  se: does the model contain user defined SEs?
  #  autocor: object of class cor_arma
  #  is_multi: is the model multivariate?
  if (is.null(se)) se <- FALSE
  if (is.formula(se)) se <- TRUE
  is.linear(family) && !is_multi && 
    (!se || get_ar(autocor) || get_ma(autocor)) 
}

needs_kronecker <- function(ranef, names_cov_ranef) {
  # checks if a model needs the kronecker product
  # Args: 
  #   ranef: named list returned by gather_ranef
  #   names_cov_ranef: names of the grouping factors that
  #                    have a cov.ranef matrix 
  .fun <- function(x, names) {
    length(x) > 1 && attr(x, "group") %in% names && attr(x, "cor")
  }
  any(sapply(ranef, .fun, names = names_cov_ranef))
}

get_boundaries <- function(trunc) {
  # extract truncation boundaries out of a formula
  # that is known to contain the .trunc function
  # Returns:
  #   a list containing two numbers named lb and ub
  if (is.formula(trunc)) {
    .addition(trunc)
  } else {
    .trunc()
  }
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

# startup messages for brms
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0(
    "Loading 'brms' package (version ", utils::packageVersion("brms"), "). ",
    "Useful instructions \n", 
    "can be found by typing help('brms'). A more detailed introduction \n", 
    "to the package is available through vignette('brms')."))
}