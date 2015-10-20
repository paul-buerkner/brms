isNULL <- function(x) {
  # check if an object is NULL
  is.null(x) || ifelse(is.vector(x), all(sapply(x, is.null)), FALSE)
}

rmNULL <- function(x) {
  # recursively removes NULL entries from an object
  x <- Filter(Negate(isNULL), x)
  lapply(x, function(x) if (is.list(x)) rmNULL(x) else x)
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

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {  
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
  #   rm: a vector of to elements indicating how many characters should be removed at the beginning
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

indicate_linear <- function(family) {
  # indicate if family is for a linear model
  if (class(family) == "family") {
    family <- family$family
  }
  family %in% c("gaussian", "student", "cauchy")
}

indicate_binary <- function(family) {
  # indicate if family is bernoulli or binomial
  if (class(family) == "family") {
    family <- family$family
  }
  family %in% c("binomial", "bernoulli")
}

indicate_ordinal <- function(family) {
  # indicate if family is for an ordinal model
  if (class(family) == "family") {
    family <- family$family
  }
  family %in% c("cumulative", "cratio", "sratio", "acat") 
}

indicate_skewed <- function(family) {
  # indicate if family is for model with postive skewed response
  if (class(family) == "family") {
    family <- family$family
  }
  family %in% c("gamma", "weibull", "exponential")
}

indicate_count <- function(family) {
  # indicate if family is for a count model
  if (class(family) == "family") {
    family <- family$family
  }
  family %in% c("poisson", "negbinomial", "geometric")
}

indicate_hurdle <- function(family) {
  # indicate if family is for a hurdle model
  if (class(family) == "family") {
    family <- family$family
  }
  family %in% c("hurdle_poisson", "hurdle_negbinomial", "hurdle_gamma")
}

indicate_zero_inflated <- function(family) {
  # indicate if family is for a zero inflated model
  if (class(family) == "family") {
    family <- family$family
  }
  family %in% c("zero_inflated_poisson", "zero_inflated_negbinomial")
}

indicate_shape <- function(family) {
  # indicate if family needs a shape parameter
  if (class(family) == "family") {
    family <- family$family
  }
  family %in% c("gamma", "weibull", "inverse.gaussian", 
                "negbinomial", "hurdle_negbinomial", 
                "hurdle_gamma", "zero_inflated_negbinomial")
}

indicate_sigma <- function(family, se, autocor) {
  # indicate if the model needs a sigma parameter
  # Args:
  #  family: a character string
  #  se: does the model contain user defined SEs?
  #      may be a formula in which case se is treated as TRUE
  #  autocor: object of class cor_arma
  is_linear <- indicate_linear(family)
  if (is.null(se)) se <- FALSE
  if (is.formula(se)) se <- TRUE
  is_linear && (!se || get_ar(autocor) || get_ma(autocor))
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

# startup messages for brms
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0(
    "Loading 'brms' package (version ", utils::packageVersion("brms"), "). ",
    "Useful instructions \n", 
    "can be found by typing help('brms'). A more detailed introduction \n", 
    "to the package is available through vignette('brms')."))
}