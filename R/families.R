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
#'   Families \code{student}, and \code{cauchy} accept the links (as names) 
#'   \code{identity}, \code{log}, and \code{inverse};
#'   families \code{negbinomial}, and \code{geometric} the links 
#'   \code{log}, \code{identity}, and \code{sqrt}; 
#'   families \code{bernoulli}, \code{cumulative}, \code{cratio}, \code{sratio}, 
#'   and \code{acat} the links \code{logit}, \code{probit}, \code{probit_approx}, 
#'   \code{cloglog}, and \code{cauchit};
#'   family \code{categorical} the link \code{logit}; families \code{weibull}, 
#'   and \code{exponential} the links \code{log}, \code{identity}, and \code{inverse};
#'   families \code{hurdle_poisson}, \code{hurdle_gamma}, 
#'   \code{hurdle_negbinomial}, \code{zero_inflated_poisson}, 
#'   and \code{zero_inflated_negbinomial} the link \code{log}.  
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
  structure(list(family = "student", link = linktemp), class = "family")
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
  structure(list(family = "cauchy", link = linktemp), class = "family")
}

#' @rdname brmsfamily
#' @export
bernoulli <- function(link = "logit") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  } 
  okLinks <- c("logit", "probit", "probit_approx", "cloglog", "cauchit")
  if (!linktemp %in% okLinks && is.character(link)) {
    linktemp <- link
  }
  if (!linktemp %in% okLinks) {
    stop(paste(linktemp, "is not a supported link for family bernoulli.", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "bernoulli", link = linktemp), class = "family")
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
  structure(list(family = "negbinomial", link = linktemp), class = "family")
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
  structure(list(family = "geometric", link = linktemp), class = "family")
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
  structure(list(family = "exponential", link = linktemp), class = "family")
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
  structure(list(family = "weibull", link = linktemp), class = "family")
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
    stop(paste(linktemp, "is not a supported link for family weibull. ", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "hurdle_poisson", link = linktemp), class = "family")
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
    stop(paste(linktemp, "is not a supported link for family weibull. ", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "hurdle_negbinomial", link = linktemp), class = "family")
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
    stop(paste(linktemp, "is not a supported link for family weibull. ", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "hurdle_gamma", link = linktemp), class = "family")
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
    stop(paste(linktemp, "is not a supported link for family weibull. ", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "zero_inflated_poisson", link = linktemp), 
            class = "family")
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
    stop(paste(linktemp, "is not a supported link for family weibull. ", 
               "Supported links are: \n", paste(okLinks, collapse = ", ")))
  }
  structure(list(family = "zero_inflated_negbinomial", link = linktemp), 
            class = "family")
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
  structure(list(family = "categorical", link = linktemp), class = "family")
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
  structure(list(family = "cumulative", link = linktemp), class = "family")
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
  structure(list(family = "sratio", link = linktemp), class = "family")
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
  structure(list(family = "cratio", link = linktemp), class = "family")
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
  structure(list(family = "acat", link = linktemp), class = "family")
}
 
 