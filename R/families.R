#' Special Family Functions for \pkg{brms} Models
#' 
#' Family objects provide a convenient way to specify the details of the models 
#' used by many model fitting functions. The families present here are currently 
#' for use with \pkg{brms} only and will NOT work with other model fitting functions such as 
#' \code{glm} or \code{glmer}. However, the standard family functions as decribed in
#' \code{\link[stats:family]{family}} will work with \pkg{brms}. For a full list 
#' of families and link functions supported by \pkg{brms}, see the documentation
#' (in particular the 'Details' section) of \code{\link[brms:brm]{brm}}.
#' 
#' @param link A specification for the model link function. 
#'   This can be a name/expression or character string. 
#'   The following list only refers to \pkg{brms} specific family functions.
#'   Families \code{student}, and \code{cauchy} (deprecated) accept the links 
#'   (as names) \code{identity}, \code{log}, and \code{inverse};
#'   families \code{negbinomial}, and \code{geometric} the links 
#'   \code{log}, \code{identity}, and \code{sqrt}; 
#'   families \code{bernoulli}, \code{Beta}, \code{cumulative}, 
#'   \code{cratio}, \code{sratio}, and \code{acat} 
#'   the links \code{logit}, \code{probit}, \code{probit_approx}, 
#'   \code{cloglog}, and \code{cauchit};
#'   family \code{categorical}, the link \code{logit}; 
#'   families \code{weibull}, and \code{exponential} 
#'   the links \code{log}, \code{identity}, and \code{inverse};
#'   families \code{hurdle_poisson}, \code{hurdle_gamma}, 
#'   \code{hurdle_negbinomial}, \code{zero_inflated_poisson}, 
#'   and \code{zero_inflated_negbinomial} the link \code{log};
#'   families \code{zero_inflated_binomial} and 
#'   \code{zero_inflated_beta} the link \code{logit}. 
#'   The first link mentioned for each family is the default.
#'   A full list of families and link functions supported by \pkg{brms}, 
#'   is provided in the 'Details' section of \code{\link[brms:brm]{brm}}.
#' @param type An optional character string allowing to specify advanced 
#'   models implemented through certain families. Currently, 
#'   only the \code{bernoulli} family uses this argument to define 
#'   2PL models (applied in IRT) by setting \code{type = "2PL"}.
#'   Further options will follow in the future.
#' 
#' @name brmsfamily
NULL

#' @rdname brmsfamily
#' @export
student <- function(link = "identity") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("identity", "log", "inverse")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link for family student.", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "student", link = linktemp), 
            class = c("brmsfamily", "family"))
}

#' @rdname brmsfamily
#' @export
cauchy <- function(link = "identity") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("identity", "log", "inverse")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link for family cauchy.", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "cauchy", link = linktemp), 
            class = c("brmsfamily", "family"))
}

#' @rdname brmsfamily
#' @export
bernoulli <- function(link = "logit", type = NULL) {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  if (!is.null(type)) {
    type <- match.arg(type, choices = "2PL")
    okLinks <- "logit"
  } else {
    okLinks <- c("logit", "probit", "probit_approx", 
                 "cloglog", "cauchit", "identity")
  }
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link for family bernoulli.", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  if (!is.null(type)) {
    type <- match.arg(type, c("2PL"))
  }
  structure(list(family = "bernoulli", link = linktemp, type = type), 
            class = c("brmsfamily", "family"))
}

#' @rdname brmsfamily
#' @export
negbinomial <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("log", "identity", "sqrt")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link for family negbinimial.", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "negbinomial", link = linktemp), 
            class = c("brmsfamily", "family"))
}

#' @rdname brmsfamily
#' @export
geometric <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("log", "identity", "sqrt")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link for family geometric.", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "geometric", link = linktemp), 
            class = c("brmsfamily", "family"))
}

#' @rdname brmsfamily
#' @export
exponential <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("log", "identity", "inverse")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link for family exponential. ", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "exponential", link = linktemp), 
            class = c("brmsfamily", "family"))
}

#' @rdname brmsfamily
#' @export
weibull <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("log", "identity", "inverse")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link for family weibull. ", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "weibull", link = linktemp), 
            class = c("brmsfamily", "family"))
}

#' @rdname brmsfamily
#' @export
Beta <- function(link = "logit") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("logit", "probit", "probit_approx", "cloglog", 
               "cauchit", "identity")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link for family beta.", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "beta", link = linktemp), 
            class = c("brmsfamily", "family"))
}

#' @rdname brmsfamily
#' @export
hurdle_poisson <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("log")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link for family hurdle_poisson. ", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "hurdle_poisson", link = linktemp), 
            class = c("brmsfamily", "family"))
}

#' @rdname brmsfamily
#' @export
hurdle_negbinomial <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("log")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link for family hurdle_negbinomial. ", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "hurdle_negbinomial", link = linktemp), 
            class = c("brmsfamily", "family"))
}

#' @rdname brmsfamily
#' @export
hurdle_gamma <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("log")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link for family hurdle_gamma. ", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "hurdle_gamma", link = linktemp), 
            class = c("brmsfamily", "family"))
}

#' @rdname brmsfamily
#' @export
zero_inflated_beta <- function(link = "logit") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("logit")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link", 
               "for family zero_inflated_beta. ", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "zero_inflated_beta", link = linktemp), 
            class = c("brmsfamily", "family"))
}

#' @rdname brmsfamily
#' @export
zero_inflated_poisson <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("log")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link", 
               "for family zero_inflated_poisson. ", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "zero_inflated_poisson", link = linktemp), 
            class = c("brmsfamily", "family"))
}

#' @rdname brmsfamily
#' @export
zero_inflated_negbinomial <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("log")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link", 
               "for family zero_inflated_negbinomial. ", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "zero_inflated_negbinomial", link = linktemp), 
            class = c("brmsfamily", "family"))
}

#' @rdname brmsfamily
#' @export
zero_inflated_binomial <- function(link = "logit") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("logit")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link", 
               "for family zero_inflated_binomial. ", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "zero_inflated_binomial", link = linktemp), 
            class = c("brmsfamily", "family"))
}

#' @rdname brmsfamily
#' @export
categorical <- function(link = "logit") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("logit")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link for family categorical.", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "categorical", link = linktemp), 
            class = c("brmsfamily", "family"))
}

#' @rdname brmsfamily
#' @export
cumulative <- function(link = "logit") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("logit", "probit", "probit_approx", "cloglog", "cauchit")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link for family cumulative.", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "cumulative", link = linktemp), 
            class = c("brmsfamily", "family"))
}

#' @rdname brmsfamily
#' @export
sratio <- function(link = "logit") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("logit", "probit", "probit_approx", "cloglog", "cauchit")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link for family sratio.",
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "sratio", link = linktemp), 
            class = c("brmsfamily", "family"))
}

#' @rdname brmsfamily
#' @export
cratio <- function(link = "logit") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("logit", "probit", "probit_approx", "cloglog", "cauchit")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link for family cratio.", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "cratio", link = linktemp), 
            class = c("brmsfamily", "family"))
}

#' @rdname brmsfamily
#' @export
acat <- function(link = "logit") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("logit", "probit", "probit_approx", "cloglog", "cauchit")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link for family acat.", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "acat", link = linktemp), 
            class = c("brmsfamily", "family"))
}

family.character <- function(object, link = NA, type = NULL, ...) {
  # build a family object
  # 
  # Args:
  #   object: A character string defining the family
  #   link: A character string defining the link
  family <- tolower(object)
  # check validity of family
  if (family == "normal")
    family <- "gaussian"
  if (family == "multigaussian") 
    stop("family 'multigaussian' is deprecated. Use family 'gaussian' instead",
         call. = FALSE)
  okFamilies <- c("gaussian", "student", "cauchy", 
                  "binomial", "bernoulli", "categorical", "beta",
                  "poisson", "negbinomial", "geometric", 
                  "gamma", "weibull", "exponential", "inverse.gaussian", 
                  "cumulative", "cratio", "sratio", "acat",
                  "hurdle_poisson", "hurdle_negbinomial", "hurdle_gamma",
                  "zero_inflated_poisson", "zero_inflated_negbinomial",
                  "zero_inflated_binomial", "zero_inflated_beta")
  if (!family %in% okFamilies) {
    stop(paste(family, "is not a supported family. Supported families are: \n",
               paste(okFamilies, collapse = ", ")), call. = FALSE)
  }
  
  # check validity of link
  if (is.linear(family)) {
    okLinks <- c("identity", "log", "inverse")
  } else if (family == "inverse.gaussian") {
    okLinks <- c("1/mu^2", "inverse", "identity", "log")
  } else if (is.count(family)) {
    okLinks <- c("log", "identity", "sqrt")
  } else if (is.ordinal(family) || family %in% "zero_inflated_beta") {
    okLinks <- c("logit", "probit", "probit_approx", "cloglog", "cauchit")
  } else if (is.binary(family) || family %in% "beta") {
    okLinks <- c("logit", "probit", "probit_approx", 
                 "cloglog", "cauchit", "identity")
  } else if (family %in% c("categorical", "zero_inflated_binomial")) {
    okLinks <- c("logit")
  } else if (is.skewed(family)) {
    okLinks <- c("log", "identity", "inverse")
  } else if (is.hurdle(family) || is.zero_inflated(family)) {
    # does not include zi_binomial or zi_beta
    okLinks <- c("log")
  } 
  if (is.na(link)) {
    link <- okLinks[1]
  }
  if (!link %in% okLinks)
    stop(paste0(link, " is not a supported link for family ", family, ". ", 
                "Supported links are: \n", paste(okLinks, collapse = ", ")))
  structure(nlist(family, link, type), class = c("brmsfamily", "family"))
}

check_family <- function(family, link = NULL) {
  # checks and corrects validity of the model family
  #
  # Args:
  #   family: Either a function, an object of class 'family' 
  #   or a character string
  #   link: an optional character string naming the link function
  #         ignored if family is a function or a family object
  if (is.function(family)) {
    family <- family()   
  }
  if (is(family, "family")) {
    family <- family(family$family, link = family$link, type = family$type)
  } else if (is.character(family)) {
    if (is.null(link)) {
      link <- family[2]
    }
    family <- family(family[1], link = link)
  } else {
    stop("family argument is invalid")
  }
  family
}

#' @export
print.brmsfamily <- function(x, ...) {
  cat("\nFamily:", x$family, "\n")
  cat("Link function:", x$link, "\n")
  if (!is.null(x$type)) {
    cat("Type:", x$type, "\n") 
  }
  cat("\n")
  invisible(x)
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
  family %in% "gaussian" && link == "log" && nresp == 1
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
  family %in% "categorical" 
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
    out <- family$family %in% "bernoulli" && identical(family$type, "2PL")
  }
  out
}

is.forked <- function(family) {
  # indicate if family has two separate model parts
  is.hurdle(family) || is.zero_inflated(family) || is.2PL(family)
}

is.mv <- function(family, response = NULL) {
  # indicate if the model uses multivariate formula syntax
  nresp <- length(response)
  is_mv <- nresp > 1L && is.linear(family) || is.categorical(family) || 
           nresp == 2L && is.forked(family)
  if (nresp > 1L && !is_mv) {
    stop("invalid multivariate model", call. = FALSE)
  }
  is_mv
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
  is.categorical(family) || is.ordinal(family)
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

allows_cse <- function(family) {
  # checks if category specific effects are allowed
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("sratio", "cratio", "acat")
}

is.old_categorical <- function(x) {
  # indicate if the model is categorical fitted with brms <= 0.8.0
  stopifnot(is(x, "brmsfit"))
  is(x$fit, "stanfit") && is.categorical(x$family) && "bp" %in% x$fit@model_pars
}
