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

eval_rhs <- function(formula, data = NULL) {
  # computes data for addition arguments
  formula <- as.formula(formula)
  eval(rhs(formula)[[2]], data, environment(formula))
}

.se <- function(x, sigma = FALSE) {
  # standard errors for meta-analysis
  if (!is.numeric(x)) {
    stop2("Standard errors must be numeric.")
  }
  if (min(x) < 0) {
    stop2("Standard errors must be non-negative.")
  }
  sigma <- as.logical(sigma)
  if (length(sigma) != 1L) {
    stop2("Argument 'sigma' must be either TRUE or FALSE.")
  }
  structure(x, sigma = sigma)  
}

.se_no_data <- function(x, sigma = FALSE) {
  # only evaluate the sigma argument
  .se(1, sigma = sigma)
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

#' Predictors with Measurement Error in \pkg{brms} Models
#' 
#' @param x The variable measured with error.
#' @param sdx Known measurement error of \code{x}
#'   treated as standard deviation.
#' 
#' @details For detailed documentation see \code{help(brmsformula)}. 
#' 
#' This function is almost solely useful when
#' called in formulas passed to the \pkg{brms} package.
#' 
#' @seealso \code{\link[brms:brmsformula]{brmsformula}}
#'   
#' @examples 
#' \dontrun{
#' # sample some data
#' N <- 100
#' dat <- data.frame(y = rnorm(N), x = rnorm(N), sdx = abs(rnorm(N, 1)))
#' # fit a simple error-in-variables model 
#' fit <- brm(y ~ me(x, sdx), data = dat, save_mevars = TRUE)
#' summary(fit)
#' } 
#' 
#' @export
me <- function(x, sdx = NULL) {
  x <- as.vector(x)
  sdx <- as.vector(sdx)
  if (length(sdx) == 0L) {
    stop2("Argument 'sdx' is missing in function 'me'.")
  } else if (length(sdx) == 1L) {
    sdx <- rep(sdx, length(x))
  }
  if (!is.numeric(x)) {
    stop2("Noisy variables should be numeric.")
  }
  if (!is.numeric(sdx)) {
    stop2("Measurement error should be numeric.")
  }
  if (any(sdx <= 0)) {
    stop2("Measurement error should be positive.")
  }
  out <- rep(1, length(x))
  structure(out, var = x, noise = sdx) 
}

#' Category Specific Predictors in \pkg{brms} Models
#' 
#' @aliases cse
#' 
#' @param expr Expression containing predictors,
#'  for which category specific effects should
#'  be estimated. For evaluation, \R formula syntax is applied.
#'  
#' @details For detailed documentation see \code{help(brmsformula)}
#'   as well as \code{vignette("brms_overview")}.
#' 
#' This function is almost solely useful when
#' called in formulas passed to the \pkg{brms} package.
#' 
#' @seealso \code{\link[brms:brmsformula]{brmsformula}}
#'   
#' @examples   
#' \dontrun{
#' fit <- brm(rating ~ period + carry + cs(treat), 
#'            data = inhaler, family = sratio("cloglog"), 
#'            prior = set_prior("normal(0,5)"), chains = 2)
#' summary(fit)
#' plot(fit, ask = FALSE)
#' } 
#'  
#' @export
cs <- function(expr) {
  expr <- substitute(expr)
  if (!is.character(expr)) {
    expr <- deparse(expr)
  }
  expr
}

#' @export
cse <- function(expr) {
  # alias of function 'cs'
  expr <- substitute(expr)
  if (!is.character(expr)) {
    expr <- deparse(expr)
  }
  expr
}

#' Monotonic Predictors in \pkg{brms} Models
#' 
#' @aliases mono monotonic
#' 
#' @param expr Expression containing predictors,
#'  for which monotonic effects should
#'  be estimated. For evaluation, \R formula syntax is applied.
#'  
#' @details For detailed documentation see \code{help(brmsformula)}
#'   as well as \code{vignette("brms_monotonic")}.
#' 
#' This function is almost solely useful when
#' called in formulas passed to the \pkg{brms} package.
#' 
#' @seealso \code{\link[brms:brmsformula]{brmsformula}}
#'   
#' @examples   
#' \dontrun{
#' # generate some data
#' income_options <- c("below_20", "20_to_40", "40_to_100", "greater_100")
#' income <- factor(sample(income_options, 100, TRUE), 
#'                  levels = income_options, ordered = TRUE)
#' mean_ls <- c(30, 60, 70, 75)
#' ls <- mean_ls[income] + rnorm(100, sd = 7)
#' dat <- data.frame(income, ls)
#' 
#' # fit a simple monotonic model
#' fit <- brm(ls ~ mo(income), data = dat)
#' 
#' # summarise the model
#' summary(fit)
#' plot(fit, N = 6)
#' plot(marginal_effects(fit), points = TRUE)
#' } 
#'  
#' @export
mo <- function(expr) {
  expr <- substitute(expr)
  if (!is.character(expr)) {
    expr <- deparse(expr)
  }
  expr
}

#' @export
mono <- function(expr) {
  # alias of function 'mo'
  expr <- substitute(expr)
  if (!is.character(expr)) {
    expr <- deparse(expr)
  }
  expr
}

#' @export
monotonic <- function(expr) {
  # alias of function 'mo'
  expr <- substitute(expr)
  if (!is.character(expr)) {
    expr <- deparse(expr)
  }
  expr
}

#' Set up basic grouping terms in \pkg{brms}
#' 
#' Function used to set up a basic grouping term in \pkg{brms}.
#' The function does not evaluate its arguments --
#' it exists purely to help set up a model with grouping terms.
#' \code{gr} is called implicitely inside the package
#' and there is usually no need to call it directly.
#' 
#' @param ... One or more terms containing grouping factors.
#' 
#' @seealso \code{\link[brms:brmsformula]{brmsformula}}
#' 
#' @examples 
#' \dontrun{
#' # model using basic lme4-style formula
#' fit1 <- brm(count ~ Trt_c + (1|patient), data = epilepsy)
#' summary(fit1)
#' # equivalent model using 'gr' which is called anyway internally
#' fit2 <- brm(count ~ Trt_c + (1|gr(patient)), data = epilepsy)
#' summary(fit2)
#' WAIC(fit1, fit2)
#' }
#' 
#' @export
gr <- function(...) {
  groups <- as.character(as.list(substitute(list(...)))[-1])
  if (length(groups) > 1L) {
    stop2("Grouping structure 'gr' expects only a single grouping term")
  }
  if (illegal_group_expr(groups[1])) {
    stop2("Illegal grouping term: ", groups[1], "\nIt may contain ",
          "only variable names combined by the symbol ':'")
  }
  allvars <- str2formula(groups)
  nlist(groups, allvars, type = "")
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
#'  Weights are standardized in order to sum to one per row.
#'  If \code{NULL} (the default), equally weights are used. 
#'  
#' @seealso \code{\link[brms:brmsformula]{brmsformula}}
#'  
#' @examples 
#' \dontrun{
#' # simulate some data
#' dat <- data.frame(y = rnorm(100), x = rnorm(100), 
#'                   g1 = sample(1:10, 100, TRUE),
#'                   g2 = sample(1:10, 100, TRUE))
#' # multi-membership model with two members per group and equal weights
#' fit1 <- brm(y ~ x + (1|mm(g1, g2)), data = dat)
#' summary(fit1)
#' 
#' # weight the first member two times for than the second member
#' dat$w1 <- rep(2, 100)
#' dat$w2 <- rep(1, 100)
#' fit2 <- brm(y ~ x + (1|mm(g1, g2, weights = cbind(w1, w2))), data = dat)
#' summary(fit2)
#' }
#'   
#' @export
mm <- function(..., weights = NULL) {
  groups <- as.character(as.list(substitute(list(...)))[-1])
  for (i in seq_along(groups)) {
    if (illegal_group_expr(groups[i])) {
      stop2("Illegal grouping term: ", groups[i], "\nIt may contain ",
            "only variable names combined by the symbol ':'")
    }
  }
  weights <- substitute(weights)
  weightvars <- all.vars(weights)
  allvars <- str2formula(c(groups, weightvars))
  if (!is.null(weights)) {
    weights <- str2formula(deparse_no_string(weights))
    weightvars <- str2formula(weightvars)
  }
  nlist(groups, weights, weightvars, allvars, type = "mm")
}

amend_formula <- function(formula, data = NULL, family = gaussian(),
                          nonlinear = NULL, partial = NULL) {
  # incorporate additional arguments into formula
  # Args:
  #   formula: object of class 'formula'
  #   data: optional data.frame
  #   family: optional object of class 'family'
  #   nonlinear, partial: deprecated arguments of brm
  formula <- bf(formula, nonlinear = nonlinear)
  attr(formula, "family") <- family
  old_attr <- attributes(formula)
  fnew <- ". ~ ."
  if (!is.null(partial)) {
    warning2("Argument 'partial' is deprecated. Please use the 'cs' ", 
             "function inside the model formula instead.")
    partial <- formula2str(partial, rm = 1)
    fnew <- paste(fnew, "+ cs(", partial, ")")
  }
  # to allow the '.' symbol in formula
  try_terms <- try(terms(formula, data = data), silent = TRUE)
  if (!is(try_terms, "try-error")) {
    formula <- formula(try_terms)
  }
  if (fnew != ". ~ .") {
    formula <- update.formula(formula, formula(fnew))
  }
  attributes(formula) <- old_attr
  if (is.categorical(family) && is.null(attr(formula, "response"))) {
    respform <- extract_effects(formula)$respform
    model_response <- model.response(model.frame(respform, data = data))
    response <- levels(factor(model_response))
    if (length(response) <= 2L) {
      stop2("At least 3 response categories are required ",
            "for family 'categorical'.")
    }
    # the first level will serve as the reference category
    attr(formula, "response") <- response[-1]
  }
  bf(formula)
}

str2formula <- function(x, ...) {
  # converts a string to a formula
  # Args:
  #   x: vector of strings to be converted
  #   ...: passed to formula(.)
  if (length(x)) {
    x <- paste(x, collapse = "+") 
  } else {
    x <- "1"
  }
  formula(paste("~", x), ...)
}

formula2str <- function(formula, rm = c(0, 0)) {
  # converts a formula to a string
  # Args:
  #   formula: a model formula
  #   rm: a vector of to elements indicating how many characters 
  #       should be removed at the beginning
  #       and end of the string respectively
  if (!is.formula(formula)) {
    formula <- as.formula(formula)
  }
  if (is.na(rm[2])) rm[2] <- 0
  x <- gsub("[ \t\r\n]+", "", Reduce(paste, deparse(formula)), perl = TRUE)
  x <- substr(x, 1 + rm[1], nchar(x) - rm[2])
  x
}

is.formula <- function(x) {
  is(x, "formula")
}
