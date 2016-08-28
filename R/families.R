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
#'   Familiy \code{student} accept the links 
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
#'   family \code{lognormal} the links \code{identity} and \code{inverse};
#'   families \code{hurdle_poisson}, \code{hurdle_gamma}, 
#'   \code{hurdle_negbinomial}, \code{zero_inflated_poisson}, 
#'   and \code{zero_inflated_negbinomial} the link \code{log};
#'   families \code{zero_inflated_binomial} and 
#'   \code{zero_inflated_beta} the link \code{logit}. 
#'   The first link mentioned for each family is the default.
#'   A full list of families and link functions supported by \pkg{brms}, 
#'   is provided in the 'Details' section of \code{\link[brms:brm]{brm}}.
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
bernoulli <- function(link = "logit") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("logit", "probit", "probit_approx", 
               "cloglog", "cauchit", "identity")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link for family bernoulli.", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "bernoulli", link = linktemp), 
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
lognormal <- function(link = "identity") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("identity", "inverse")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link for family lognormal.", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "lognormal", link = linktemp), 
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
von_mises <- function(link = "tan_half") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  }
  okLinks <- c("tan_half")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link for family von_mises.", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "von_mises", link = linktemp), 
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

family.character <- function(object, link = NA, ...) {
  # build a family object
  # Args:
  #   object: A character string defining the family
  #   link: A character string defining the link
  family <- tolower(object)
  # check validity of family
  if (family == "normal") {
    family <- "gaussian"
  }
  okFamilies <- c("gaussian", "student", "lognormal", 
                  "binomial", "bernoulli", "categorical", 
                  "poisson", "negbinomial", "geometric", 
                  "gamma", "weibull", "exponential", 
                  "inverse.gaussian", "beta", "von_mises",
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
  } else if (family %in% "lognormal") {
    okLinks <- c("identity", "inverse")
  } else if (family %in% "von_mises") {
    okLinks <- c("tan_half")
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
  structure(nlist(family, link), class = c("brmsfamily", "family"))
}

check_family <- function(family, link = NULL) {
  # checks and corrects validity of the model family
  # Args:
  #   family: Either a function, an object of class 'family' 
  #   or a character string
  #   link: an optional character string naming the link function
  #         ignored if family is a function or a family object
  if (is.function(family)) {
    family <- family()   
  }
  if (is(family, "family")) {
    family <- family(family$family, link = family$link)
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

is.lognormal <- function(family) {
  # indicate if family is lognormal
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% "lognormal"
}

is.count <- function(family) {
  # indicate if family is for a count model
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("poisson", "negbinomial", "geometric")
}

is.hurdle <- function(family, zi_beta = TRUE) {
  # indicate if family is for a hurdle model
  if (is(family, "family")) {
    family <- family$family
  }
  # zi_beta is technically a hurdle model
  family %in% c("hurdle_poisson", "hurdle_negbinomial", "hurdle_gamma",
                if (zi_beta) "zero_inflated_beta")
}

is.zero_inflated <- function(family, zi_beta = FALSE) {
  # indicate if family is for a zero inflated model
  if (is(family, "family")) {
    family <- family$family
  }
  # zi_beta is technically a hurdle model
  family %in% c("zero_inflated_poisson", "zero_inflated_negbinomial",
                "zero_inflated_binomial", if (zi_beta) "zero_inflated_beta")
}

is.2PL <- function(family) {
  if (!is(family, "brmsfamily")) {
    out <- FALSE
  } else {
    out <- family$family %in% "bernoulli" && identical(family$type, "2PL")
  }
  if (out) {
    stop("The special implementation of 2PL models has been removed.\n",
         "You can now use argument 'nonlinear' to fit such models.",
         call. = FALSE)
  }
  out
}

is.forked <- function(family) {
  # indicate if family has two separate model parts
  is.hurdle(family) || is.zero_inflated(family) || is.2PL(family)
}

is.mv <- function(family, response = NULL) {
  # indicate if the model uses multiple responses
  nresp <- length(response)
  is_mv <- nresp > 1L && is.linear(family) || is.categorical(family) || 
           nresp == 2L && is.forked(family)
  if (nresp > 1L && !is_mv) {
    stop("Invalid multivariate model", call. = FALSE)
  }
  is_mv
}

use_real <- function(family) {
  # indicate if family uses real responses
  if (is(family, "family")) {
    family <- family$family
  }
  is.linear(family) || is.skewed(family) || 
    family %in% c("lognormal", "inverse.gaussian", "beta", "von_mises",
                  "zero_inflated_beta", "hurdle_gamma")
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

has_nu <- function(family) {
  # indicate if family needs a nu parameter
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("student")
}

has_phi <- function(family) {
  # indicate if family needs a phi parameter
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("beta", "zero_inflated_beta")
}

has_kappa <- function(family) {
  # indicate if family needs a kappa parameter
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("von_mises")
}

has_sigma <- function(family, effects = NULL, autocor = cor_arma(),
                      incmv = FALSE) {
  # indicate if the model needs a sigma parameter
  # Args:
  #  family: model family
  #  effects: list returned by extract_effects
  #  autocor: object of class cor_arma
  #  incmv: should MV (linear) models be treated as having sigma? 
  has_se <- !is.null(effects$se)
  out <- (is.linear(family) || is.lognormal(family)) && 
         (!has_se || use_cov(autocor)) && !is(autocor, "cov_fixed")
  if (!incmv) {
    is_multi <- is.linear(family) && length(effects$response) > 1L
    out <- out && !is_multi
  }
  out
}

allows_cse <- function(family) {
  # checks if category specific effects are allowed
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("sratio", "cratio", "acat")
}

is.old_lognormal <- function(family, link = "identity", nresp = 1,
                             version = utils::packageVersion("brms")) {
  # indicate transformation to lognormal models
  # Args:
  #   link: A character string; ignored if family is of class family
  #   nresp: number of response variables
  #   version: brms version with which the model was fitted
  if (is(family, "family")) {
    link <- family$link
    family <- family$family
  }
  family %in% "gaussian" && link == "log" && nresp == 1 &&
    (is.null(version) || version <= "0.9.1")
}

is.old_categorical <- function(x) {
  # indicate if the model is and old categorical model
  stopifnot(is(x, "brmsfit"))
  if (is(x$fit, "stanfit") && is.categorical(x$family)) {
    if ("bp" %in% x$fit@model_pars) {
      # fitted with brms <= 0.8.0
      out <- 1L
    } else if (is.old_mv(x)) {
      # fitted with brms <= 1.0.0
      out <- 2L
    } else {
      out <- 0L
    }
  } else {
    out <- 0L
  }
  out
}

is.old_mv <- function(x) {
  # indicate if the model uses the old multivariate syntax 
  # from brms < 1.0.0
  stopifnot(is(x, "brmsfit"))
  ee <- extract_effects(formula(x), family = family(x))
  (is.null(x$version) || x$version <= "0.10.0.9000") &&
    (is.mv(family(x), ee$response) || is.forked(family(x)))
}
