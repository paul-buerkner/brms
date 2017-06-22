#' Special Family Functions for \pkg{brms} Models
#' 
#' Family objects provide a convenient way to specify the details of the models 
#' used by many model fitting functions. The familiy functions present here are 
#' currently for use with \pkg{brms} only and will NOT work with other model 
#' fitting functions such as \code{glm} or \code{glmer}. 
#' However, the standard family functions as decribed in
#' \code{\link[stats:family]{family}} will work with \pkg{brms}.
#' 
#' @param family A character string naming the distribution
#'   of the response variable be used in the model.
#'   Currently, the following families are supported:
#'   \code{gaussian}, \code{student}, \code{binomial}, 
#'   \code{bernoulli}, \code{poisson}, \code{negbinomial}, 
#'   \code{geometric}, \code{Gamma}, \code{lognormal}, 
#'   \code{exgaussian}, \code{skew_normal}, \code{wiener}, 
#'   \code{inverse.gaussian}, \code{exponential}, \code{weibull}, 
#'   \code{frechet}, \code{Beta}, \code{von_mises}, \code{asym_laplace},
#'   \code{gen_extreme_value}, \code{categorical}, \code{cumulative}, 
#'   \code{cratio}, \code{sratio}, \code{acat}, \code{hurdle_poisson}, 
#'   \code{hurdle_negbinomial}, \code{hurdle_gamma}, \code{hurdle_lognormal},
#'   \code{zero_inflated_binomial}, \code{zero_inflated_beta},
#'   \code{zero_inflated_negbinomial}, \code{zero_inflated_poisson},
#'   and \code{zero_one_inflated_beta}.
#' @param link A specification for the model link function. 
#'   This can be a name/expression or character string. 
#'   See the 'Details' section for more information on link
#'   functions supported by each family.
#' @param link_sigma Link of auxiliary parameter \code{sigma} if being predicted.
#' @param link_shape Link of auxiliary parameter \code{shape} if being predicted.
#' @param link_nu Link of auxiliary parameter \code{nu} if being predicted.
#' @param link_phi Link of auxiliary parameter \code{phi} if being predicted.
#' @param link_kappa Link of auxiliary parameter \code{kappa} if being predicted.
#' @param link_beta Link of auxiliary parameter \code{beta} if being predicted.
#' @param link_zi Link of auxiliary parameter \code{zi} if being predicted.
#' @param link_hu Link of auxiliary parameter \code{hu} if being predicted.
#' @param link_zoi Link of auxiliary parameter \code{zoi} if being predicted.
#' @param link_coi Link of auxiliary parameter \code{coi} if being predicted.
#' @param link_disc Link of auxiliary parameter \code{disc} if being predicted.
#' @param link_bs Link of auxiliary parameter \code{bs} if being predicted.
#' @param link_ndt Link of auxiliary parameter \code{ndt} if being predicted.
#' @param link_bias Link of auxiliary parameter \code{bias} if being predicted.
#' @param link_alpha Link of auxiliary parameter \code{alpha} if being predicted.
#' @param link_quantile Link of auxiliary parameter \code{quantile} if being predicted.
#' @param link_xi Link of auxiliary parameter \code{xi} if being predicted.
#' 
#' @details 
#'   Family \code{gaussian} with \code{identity} link leads to linear regression. 
#'   Family \code{student} with \code{identity} link leads to 
#'   robust linear regression that is less influenced by outliers. 
#'   Family \code{skew_normal} can handle skewed responses in linear regression.
#'   Families \code{poisson}, \code{negbinomial}, and \code{geometric} 
#'   with \code{log} link lead to regression models for count data. 
#'   Families \code{binomial} and \code{bernoulli} with \code{logit} link leads to 
#'   logistic regression and family \code{categorical} to multi-logistic regression 
#'   when there are more than two possible outcomes.
#'   Families \code{cumulative}, \code{cratio} ('continuation ratio'), 
#'   \code{sratio} ('stopping ratio'), and \code{acat} ('adjacent category') 
#'   leads to ordinal regression. Families \code{Gamma}, \code{weibull}, 
#'   \code{exponential}, \code{lognormal}, \code{frechet}, and 
#'   \code{inverse.gaussian} can be used (among others) for survival regression.
#'   Families \code{weibull}, \code{frechet}, and \code{gen_extreme_value}
#'   ('generalized extreme value') allow for modeling extremes.
#'   Family \code{asym_laplace} allows for quantile regression when fixing
#'   the auxiliary \code{quantile} parameter to the quantile of interest.
#'   Family \code{exgaussian} ('exponentially modified Gaussian') is especially
#'   suited to model reaction times and the \code{wiener} family provides
#'   an implementation of the Wiener diffusion model. For this family,
#'   the main formula predicts the drift parameter 'delta' and
#'   all other parameters are modeled as auxiliary parameters 
#'   (see \code{\link[brms:brmsformula]{brmsformula}} for details).
#'   Families \code{hurdle_poisson}, \code{hurdle_negbinomial}, 
#'   \code{hurdle_gamma}, \code{hurdle_lognormal}, \code{zero_inflated_poisson},
#'   \code{zero_inflated_negbinomial}, \code{zero_inflated_binomial},
#'   \code{zero_inflated_beta}, and \code{zero_one_inflated_beta} 
#'   allow to estimate zero-inflated and hurdle models. 
#'   These models can be very helpful when there are many zeros in the data 
#'   (or ones in case of one-inflated models)
#'   that cannot be explained by the primary distribution of the response. 
#'   Families \code{hurdle_lognormal} and \code{hurdle_gamma} are 
#'   especially useful, as traditional \code{lognormal} or \code{Gamma}
#'   models cannot be reasonably fitted for data containing zeros in the response.
#'   
#'   In the following, we list all possible links for each family.
#'   The families \code{gaussian}, \code{student}, \code{skew_normal},
#'   \code{exgaussian}, \code{asym_laplace}, and \code{gen_extreme_value} 
#'   accept the links (as names) \code{identity}, \code{log}, and \code{inverse};
#'   families \code{poisson}, \code{negbinomial}, \code{geometric},
#'   \code{zero_inflated_poisson}, \code{zero_inflated_negbinomial},
#'   \code{hurdle_poisson}, and \code{hurdle_negbinomial} the links 
#'   \code{log}, \code{identity}, and \code{sqrt}; 
#'   families \code{binomial}, \code{bernoulli}, \code{Beta}, 
#'   \code{zero_inflated_binomial}, \code{zero_inflated_beta}, 
#'   and \code{zero_one_inflated_beta} the links \code{logit}, 
#'   \code{probit}, \code{probit_approx}, \code{cloglog}, 
#'   \code{cauchit}, and \code{identity}; 
#'   families \code{cumulative}, \code{cratio}, \code{sratio}, 
#'   and \code{acat} the links \code{logit}, \code{probit}, 
#'   \code{probit_approx}, \code{cloglog}, and \code{cauchit}; 
#'   family \code{categorical} the link \code{logit};
#'   families \code{Gamma}, \code{weibull}, \code{exponential}, 
#'   \code{frechet}, and \code{hurdle_gamma} the links 
#'   \code{log}, \code{identity}, and \code{inverse};
#'   families \code{lognormal} and \code{hurdle_lognormal} 
#'   the links \code{identity} and \code{inverse};
#'   family \code{inverse.gaussian} the links \code{1/mu^2}, 
#'   \code{inverse}, \code{identity} and \code{log};
#'   family \code{von_mises} the link \code{tan_half};
#'   family \code{wiener} the link \code{identity}.
#'   The first link mentioned for each family is the default.     
#'   
#'   Please note that when calling the \code{\link[stats:family]{Gamma}} 
#'   family function, the default link will be \code{inverse} not \code{log}. 
#'   Also, the \code{probit_approx} link cannot be used when calling the
#'   \code{\link[stats:family]{binomial}} family function. 
#'   
#'   The current implementation of \code{inverse.gaussian} models has some 
#'   convergence problems and requires carefully chosen prior distributions 
#'   to work efficiently. For this reason, we currently do not recommend
#'   to use the \code{inverse.gaussian} family, unless you really feel
#'   that your data requires exactly this type of model.
#'
#' @seealso \code{\link[brms:brm]{brm}}, 
#'   \code{\link[stats:family]{family}}
#'   
#' @examples 
#'  # create a family object
#'  (fam1 <- student("log"))
#'  # alternatively use the brmsfamily function
#'  (fam2 <- brmsfamily("student", "log"))
#'  # both leads to the same object
#'  identical(fam1, fam2) 
#' 
#' @export
brmsfamily <- function(family, link = NULL, link_sigma = "log", 
                       link_shape = "log", link_nu = "logm1",
                       link_phi = "log", link_kappa = "log",
                       link_beta = "log", link_zi = "logit", 
                       link_hu = "logit", link_zoi = "logit",
                       link_coi = "logit", link_disc = "log",
                       link_bs = "log", link_ndt = "log",
                       link_bias = "logit", link_xi = "log1p",
                       link_alpha = "identity", 
                       link_quantile = "logit") {
  slink <- substitute(link)
  .brmsfamily(
    family, link = link, slink = slink,
    link_sigma = link_sigma, link_shape = link_shape, 
    link_nu = link_nu, link_phi = link_phi, 
    link_kappa = link_kappa, link_beta = link_beta, 
    link_zi = link_zi, link_hu = link_hu, 
    link_zoi = link_zoi, link_coi = link_coi,
    link_disc = link_disc, link_bs = link_bs, 
    link_ndt = link_ndt, link_bias = link_bias,
    link_alpha = link_alpha, link_xi = link_xi,
    link_quantile = link_quantile
  )
}

.brmsfamily <- function(family, link = NULL, slink = link, ...) {
  # helper function to prepare brmsfamily objects
  # Args:
  #   family: character string naming the model family
  #   link: character string naming the link function
  #   slink: can be used with substitute(link) for 
  #          non-standard evaluation of the link function
  #   ...: link functions (as character strings) of auxiliary parameters
  # returns:
  #  An object of class = c(brmsfamily, family) to be used
  #  only insided the brms package
  family <- tolower(as.character(family))
  aux_links <- list(...)
  if (length(family) != 1L) {
    stop2("Argument 'family' must be of length 1.")
  }
  family <- rename(family, symbols = c("^normal$", "^zi_", "^hu_"),
                   subs = c("gaussian", "zero_inflated_", "hurdle_"),
                   fixed = FALSE)
  ok_families <- c(
    "gaussian", "student", "lognormal", "skew_normal",
    "binomial", "bernoulli", "categorical", 
    "poisson", "negbinomial", "geometric", 
    "gamma", "weibull", "exponential", "exgaussian", 
    "frechet", "gen_extreme_value", "inverse.gaussian", 
    "wiener", "beta", "von_mises", "asym_laplace",
    "cumulative", "cratio", "sratio", "acat",
    "hurdle_poisson", "hurdle_negbinomial", "hurdle_gamma",
    "hurdle_lognormal", "zero_inflated_poisson", 
    "zero_inflated_negbinomial", "zero_inflated_binomial", 
    "zero_inflated_beta", "zero_one_inflated_beta"
  )
  if (!family %in% ok_families) {
    stop(family, " is not a supported family. Supported families are: \n",
         paste(ok_families, collapse = ", "), call. = FALSE)
  }
  
  # check validity of link
  is_linear <- family %in% c(
    "gaussian", "student", "skew_normal", "exgaussian", 
    "asym_laplace", "gen_extreme_value"
  )  
  is_count <- family %in% c(
    "poisson", "negbinomial", "geometric", "hurdle_poisson",
    "hurdle_negbinomial", "zero_inflated_poisson",
    "zero_inflated_negbinomial"
  )
  is_binomial <- family %in% c(
    "binomial", "bernoulli", "zero_inflated_binomial"
  )
  is_beta <- family %in% c(
    "beta", "zero_inflated_beta", "zero_one_inflated_beta"
  )
  is_skewed <- family %in% c(
    "gamma", "weibull", "exponential", "frechet", "hurdle_gamma"
  )
  is_lognormal <- family %in% c("lognormal", "hurdle_lognormal")
  if (is_linear) {
    ok_links <- c("identity", "log", "inverse")
  } else if (is_count) {
    ok_links <- c("log", "identity", "sqrt")
  } else if (is_binomial || is_beta) {
    ok_links <- c(
      "logit", "probit", "probit_approx", 
      "cloglog", "cauchit", "identity"
    )
  } else if (is_skewed) {
    ok_links <- c("log", "identity", "inverse")
  } else if (is_lognormal) {
    ok_links <- c("identity", "inverse")
  } else if (is_categorical(family)) {
    ok_links <- c("logit")
  } else if (is_ordinal(family)) {
    ok_links <- c(
      "logit", "probit", "probit_approx", "cloglog", "cauchit"
    )
  } else if (family %in% "inverse.gaussian") {
    ok_links <- c("1/mu^2", "inverse", "identity", "log")
  } else if (is_wiener(family)) {
    ok_links <- c("identity")
  } else if (family %in% "von_mises") {
    ok_links <- c("tan_half")
  }
  # non-standard evaluation of link
  if (!is.character(slink)) {
    slink <- deparse(slink)
  } 
  if (!slink %in% ok_links) {
    if (is.character(link)) {
      slink <- link
    } else if (!length(link) || identical(link, NA)) {
      slink <- NA
    }
  }
  if (length(slink) != 1L) {
    stop2("Argument 'link' must be of length 1.")
  }
  if (is.na(slink)) {
    slink <- ok_links[1]
  } 
  if (!slink %in% ok_links) {
    stop2("'", slink, "' is not a supported link ", 
          "for family '", family, "'.\nSupported links are: ",
          collapse_comma(ok_links))
  }
  out <- structure(
    list(family = family, link = slink), 
    class = c("brmsfamily", "family")
  )
  for (ap in valid_auxpars(out$family)) {
    alink <- as.character(aux_links[[paste0("link_", ap)]])
    if (length(alink)) {
      if (length(alink) > 1L) {
        stop2("Link functions must be of length 1.")
      }
      valid_links <- links_auxpars(ap)
      if (!alink %in% valid_links) {
        stop2("'", alink, "' is not a supported link ", 
              "for parameter '", ap, "'.\nSupported links are: ", 
              collapse_comma(valid_links))
      }
      out[[paste0("link_", ap)]] <- alink
    }
  }
  out
}

#' @rdname brmsfamily
#' @export
student <- function(link = "identity", link_sigma = "log", link_nu = "logm1") {
  slink <- substitute(link)
  .brmsfamily("student", link = link, slink = slink, 
              link_sigma = link_sigma, link_nu = link_nu)
}

#' @rdname brmsfamily
#' @export
bernoulli <- function(link = "logit") {
  slink <- substitute(link)
  .brmsfamily("bernoulli", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
negbinomial <- function(link = "log", link_shape = "log") {
  slink <- substitute(link)
  .brmsfamily("negbinomial", link = link, slink = slink,
              link_shape = link_shape)
}

#' @rdname brmsfamily
#' @export
geometric <- function(link = "log") {
  slink <- substitute(link)
  .brmsfamily("geometric", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
lognormal <- function(link = "identity", link_sigma = "log") {
  slink <- substitute(link)
  .brmsfamily("lognormal", link = link, slink = slink,
              link_sigma = link_sigma)
}

#' @rdname brmsfamily
#' @export
skew_normal <- function(link = "identity", link_sigma = "log", 
                        link_alpha = "identity") {
  slink <- substitute(link)
  .brmsfamily("skew_normal", link = link, slink = slink,
              link_sigma = link_sigma, link_alpha = link_alpha)
}

#' @rdname brmsfamily
#' @export
exponential <- function(link = "log") {
  slink <- substitute(link)
  .brmsfamily("exponential", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
weibull <- function(link = "log", link_shape = "log") {
  slink <- substitute(link)
  .brmsfamily("weibull", link = link, slink = slink,
              link_shape = link_shape)
}

#' @rdname brmsfamily
#' @export
frechet <- function(link = "log", link_nu = "logm1") {
  slink <- substitute(link)
  .brmsfamily("frechet", link = link, slink = slink,
              link_nu = link_nu)
}

#' @rdname brmsfamily
#' @export
gen_extreme_value <- function(link = "identity", link_sigma = "log",
                              link_xi = "log1p") {
  slink <- substitute(link)
  .brmsfamily("gen_extreme_value", link = link, slink = slink,
              link_sigma = link_sigma, link_xi = link_xi)
}

#' @rdname brmsfamily
#' @export
exgaussian <- function(link = "identity", link_sigma = "log",
                       link_beta = "log") {
  slink <- substitute(link)
  .brmsfamily("exgaussian", link = link, slink = slink,
              link_sigma = link_sigma, link_beta = link_beta)
}

#' @rdname brmsfamily
#' @export
wiener <- function(link = "identity", link_bs = "log", 
                   link_ndt = "log", link_bias = "logit") {
  slink <- substitute(link)
  .brmsfamily("wiener", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
Beta <- function(link = "logit", link_phi = "log") {
  slink <- substitute(link)
  .brmsfamily("beta", link = link, slink = slink,
              link_phi = link_phi)
}

#' @rdname brmsfamily
#' @export
von_mises <- function(link = "tan_half", link_kappa = "log") {
  slink <- substitute(link)
  .brmsfamily("von_mises", link = link, slink = slink,
              link_kappa = link_kappa)
}

#' @rdname brmsfamily
#' @export
asym_laplace <- function(link = "identity", link_sigma = "log",
                         link_quantile = "logit") {
  slink <- substitute(link)
  .brmsfamily("asym_laplace", link = link, slink = slink,
              link_sigma = link_sigma, link_quantile = link_quantile)
}

#' @rdname brmsfamily
#' @export
hurdle_poisson <- function(link = "log") {
  slink <- substitute(link)
  .brmsfamily("hurdle_poisson", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
hurdle_negbinomial <- function(link = "log", link_shape = "log",
                               link_hu = "logit") {
  slink <- substitute(link)
  .brmsfamily("hurdle_negbinomial", link = link, slink = slink,
              link_shape = link_shape, link_hu = link_hu)
}

#' @rdname brmsfamily
#' @export
hurdle_gamma <- function(link = "log", link_shape = "log",
                         link_hu = "logit") {
  slink <- substitute(link)
  .brmsfamily("hurdle_gamma", link = link, slink = slink,
              link_shape = link_shape, link_hu = link_hu)
}

#' @rdname brmsfamily
#' @export
hurdle_lognormal <- function(link = "identity", link_sigma = "log",
                             link_hu = "logit") {
  slink <- substitute(link)
  .brmsfamily("hurdle_lognormal", link = link, slink = slink,
              link_sigma = link_sigma, link_hu = link_hu)
}

#' @rdname brmsfamily
#' @export
zero_inflated_beta <- function(link = "logit", link_phi = "log",
                               link_zi = "logit") {
  slink <- substitute(link)
  .brmsfamily("zero_inflated_beta", link = link, slink = slink,
              link_phi = link_phi, link_zi = link_zi)
}

#' @rdname brmsfamily
#' @export
zero_one_inflated_beta <- function(link = "logit", link_phi = "log",
                                   link_zoi = "logit", link_coi = "logit") {
  slink <- substitute(link)
  .brmsfamily("zero_one_inflated_beta", link = link, slink = slink,
              link_phi = link_phi, link_zoi = link_zoi,
              link_coi = link_coi)
}

#' @rdname brmsfamily
#' @export
zero_inflated_poisson <- function(link = "log", link_zi = "logit") {
  slink <- substitute(link)
  .brmsfamily("zero_inflated_poisson", link = link, slink = slink,
              link_zi = link_zi)
}

#' @rdname brmsfamily
#' @export
zero_inflated_negbinomial <- function(link = "log", link_shape = "log",
                                      link_zi = "logit") {
  slink <- substitute(link)
  .brmsfamily("zero_inflated_negbinomial", link = link, slink = slink,
              link_shape = link_shape, link_zi = link_zi)
}

#' @rdname brmsfamily
#' @export
zero_inflated_binomial <- function(link = "logit", link_zi = "logit") {
  slink <- substitute(link)
  .brmsfamily("zero_inflated_binomial", link = link, slink = slink,
              link_zi = link_zi)
}

#' @rdname brmsfamily
#' @export
categorical <- function(link = "logit") {
  slink <- substitute(link)
  .brmsfamily("categorical", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
cumulative <- function(link = "logit", link_disc = "log") {
  slink <- substitute(link)
  .brmsfamily("cumulative", link = link, slink = slink,
              link_disc = link_disc)
}

#' @rdname brmsfamily
#' @export
sratio <- function(link = "logit", link_disc = "log") {
  slink <- substitute(link)
  .brmsfamily("sratio", link = link, slink = slink,
              link_disc = link_disc)
}

#' @rdname brmsfamily
#' @export
cratio <- function(link = "logit", link_disc = "log") {
  slink <- substitute(link)
  .brmsfamily("cratio", link = link, slink = slink,
              link_disc = link_disc)
}

#' @rdname brmsfamily
#' @export
acat <- function(link = "logit", link_disc = "log") {
  slink <- substitute(link)
  .brmsfamily("acat", link = link, slink = slink,
              link_disc = link_disc)
}

par_family <- function(par = NULL, link = NULL) {
  # set up special family objects for parameters
  # Args:
  #   par: name of the parameter
  #   link: optional link function of the parameter
  ap_class <- auxpar_class(par)
  if (!is.null(par) && !isTRUE(ap_class %in% auxpars())) {
    stop2("Parameter '", par, "' is invalid.")
  }
  if (is.null(par)) {
    link <- "identity"
  } else {
    links <- links_auxpars(ap_class)
    if (is.null(link)) {
      link <- links[1]
    } else {
      if (!isTRUE(link %in% links)) {
        stop2("Link '", link, "' is invalid for parameter '", par, "'.")
      }
    }
  }
  structure(
    list(family = "", link = link, par = par),
    class = c("brmsfamily", "family")
  )
}

check_family <- function(family, link = NULL) {
  # checks and corrects validity of the model family
  # Args:
  #   family: Either a function, an object of class 'family' 
  #           or a character string of length one or two
  #   link: an optional character string naming the link function
  #         ignored if family is a function or a family object
  if (is.function(family)) {
    family <- family()   
  }
  if (!is(family, "brmsfamily")) {
    if (is.family(family)) {
      link <- family$link
      family <- family$family
    } 
    if (is.character(family)) {
      if (is.null(link)) {
        link <- family[2]
      }
      family <- .brmsfamily(family[1], link = link)
    } else {
      stop2("Argument 'family' is invalid.")
    } 
  }
  family
}

#' Finite Mixture Families in \pkg{brms}
#' 
#' Set up a finite mixture family for use in \pkg{brms}.
#' 
#' @param ... One or more objects providing a description of the 
#'   response distributions to be combined in the mixture model. 
#'   These can be family functions, calls to family functions or 
#'   character strings naming the families.
#'   For details of supported families see 
#'   \code{\link[brms:brmsfamily]{brmsfamily}}.
#' @param flist Optional list of objects, which are treated in the 
#'   same way as objects passed via the \code{...} argument.
#' @param nmix Optional numeric vector specifying the number of times
#'   each family is repeated. If specified, it must have the same length 
#'   as the number of families passed via \code{...} or \code{flist}.
#' @param order Ordering constraint to identify mixture components.
#'   If \code{'mu'} or \code{TRUE}, population-level intercepts
#'   of the mean parameters are ordered. 
#'   If \code{'none'} or \code{FALSE}, no ordering constraint is applied.
#'   If \code{NULL} (the default), \code{order} is set to \code{'mu'}
#'   if all families are the same and \code{'none'} otherwise.
#'   Other ordering constraints may be implemented in the future.
#'
#' @return An object of class \code{mixfamily}.
#' 
#' @details
#' 
#' Most families supported by \pkg{brms} can be used to form 
#' mixtures. The response variable has to be valid for all components
#' of the mixture family. Currently, the number of mixture components 
#' has to be specified by the user. It is not yet possible to estimate 
#' the number of mixture components from the data.
#' 
#' For most mixture models, you may want to specify priors on the population-level
#' intercepts via \code{\link[brms:set_prior]{set_prior}} to improve convergence. 
#' In addition, it is sometimes necessary to set \code{inits = 0} in the call to 
#' \code{\link[brms:brms]{brm}} to allow chains to initialize properly.
#' 
#' For more details on the specification of mixture
#' models, see \code{\link[brms:brmsformula]{brmsformula}}.
#' 
#' @examples
#' \dontrun{
#' ## simulate some data
#' set.seed(1234)
#' dat <- data.frame(
#'   y = c(rnorm(200), rnorm(100, 6)), 
#'   x = rnorm(300),
#'   z = sample(0:1, 300, TRUE)
#' )
#' 
#' ## fit a simple normal mixture model
#' mix <- mixture(gaussian, gaussian)
#' prior <- c(
#'   prior(normal(0, 7), Intercept, nlpar = mu1),
#'   prior(normal(5, 7), Intercept, nlpar = mu2)
#' )
#' fit1 <- brm(bf(y ~ x + z), dat, family = mix,
#'             prior = prior, chains = 2) 
#' summary(fit1)
#' pp_check(fit1)
#' 
#' ## use different predictors for the components
#' fit2 <- brm(bf(y ~ 1, mu1 ~ x, mu2 ~ z), dat, family = mix,
#'             prior = prior, chains = 2) 
#' summary(fit2)
#' 
#' ## fix the mixing proportions
#' fit3 <- brm(bf(y ~ x + z, theta1 = 1, theta2 = 2), 
#'             dat, family = mix, prior = prior, 
#'             inits = 0, chains = 2)
#' summary(fit3)
#' pp_check(fit3)    
#' 
#' ## predict the mixing proportions
#' fit4 <- brm(bf(y ~ x + z, theta2 ~ x), 
#'             dat, family = mix, prior = prior, 
#'             inits = 0, chains = 2)
#' summary(fit4)
#' pp_check(fit4)           
#'
#' ## compare model fit
#' LOO(fit1, fit2, fit3, fit4)  
#' }
#' 
#' @export
mixture <- function(..., flist = NULL, nmix = 1, order = NULL) {
  dots <- c(list(...), flist)
  if (length(nmix) == 1L) {
    nmix <- rep(nmix, length(dots))
  }
  if (length(dots) != length(nmix)) {
    stop2("The length of 'nmix' should be the same ", 
          "as the number of mixture components.")
  }
  dots <- dots[rep(seq_along(dots), nmix)]
  family <- list(
    family = "mixture", 
    link = "identity",
    mix = lapply(dots, check_family)
  )
  class(family) <- c("mixfamily", "brmsfamily", "family")
  # validity checks
  families <- family_names(family)
  if (length(families) < 2L) {
    stop2("Expecting at least 2 mixture components.")
  }
  non_mix_families <- c("categorical")
  non_mix_families <- intersect(families, non_mix_families)
  if (length(non_mix_families)) {
    stop2("Families ", collapse_comma(non_mix_families), 
          " are currently not allowed in mixture models.")
  }
  ordinal_families <- c("cumulative", "sratio", "cratio", "acat")
  if (is_ordinal(family) && any(!families %in% ordinal_families)) {
    stop2("Cannot mix ordinal and non-ordinal families.")
  }
  if (use_real(family) && use_int(family)) {
    stop2("Cannot mix families with real and integer support.")
  }
  if (is.null(order)) {
    if (length(unique(families)) == 1L) {
      family$order <- "mu"
      message("Setting order = 'mu' for mixtures of the same family.")
    } else {
      family$order <- "none"
      message("Setting order = 'none' for mixtures of different families.")
    }
  } else {
    if (length(order) != 1L) {
      stop2("Argument 'order' must be of length 1.")
    }
    if (is.character(order)) {
      valid_order <- c("none", "mu")
      if (!order %in% valid_order) {
        stop2("Argument 'order' is invalid. Valid options are: ",
              collapse_comma(valid_order))
      }
      family$order <- order
    } else {
      family$order <- ifelse(as.logical(order), "mu", "none")
    }
  }
  family
}

#' @export
family_names.default <- function(family, ...) {
  family
}

#' @export
family_names.family <- function(family, ...) {
  family$family
}

#' @export
family_names.mixfamily <- function(family, ...) {
  ulapply(family$mix, "[[", "family")
}

#' @export
auxpar_family.default <- function(family, auxpar, ...) {
  ap_class <- auxpar_class(auxpar)
  if (!identical(ap_class, "mu")) {
    link <- family[[paste0("link_", ap_class)]]
    family <- par_family(auxpar, link)
  }
  family
}

#' @export
auxpar_family.mixfamily <- function(family, auxpar, ...) {
  ap_id <- as.numeric(auxpar_id(auxpar))
  if (!(length(ap_id) == 1L && is.numeric(ap_id))) {
    stop2("Parameter '", auxpar, "' is not a valid mixture parameter.")
  }
  auxpar_family(family$mix[[ap_id]], auxpar, ...)
}

#' @export
print.brmsfamily <- function(x, links = FALSE, newline = TRUE, ...) {
  cat("\nFamily:", x$family, "\n")
  cat("Link function:", x$link, "\n")
  if (!is.null(x$type)) {
    cat("Type:", x$type, "\n") 
  }
  if (isTRUE(links) || is.character(links)) {
    ap_links <- x[grepl("^link_", names(x))]
    names(ap_links) <- sub("^link_", "", names(ap_links))
    if (is.character(links)) {
      ap_links <- rmNULL(ap_links[links])
    }
    for (ap in names(ap_links)) {
      cat(paste0("Link function of ", ap, ": ", ap_links[[ap]], " \n")) 
    }
  }
  if (newline) {
    cat("\n")
  }
  invisible(x)
}

#' @export
print.mixfamily <- function(x, newline = TRUE, ...) {
  cat("\nMixture\n")
  for (i in seq_along(x$mix)) {
    print(x$mix[[i]], newline = FALSE, ...)
  }
  if (newline) {
    cat("\n")
  }
  invisible(x)
}

#' @method summary family
#' @export
summary.family <- function(object, link = TRUE, ...) {
  out <- object$family
  if (link) {
    out <- paste0(out, "(", object$link, ")")
  }
  out
}

#' @method summary mixfamily
#' @export
summary.mixfamily <- function(object, link = FALSE, ...) {
  families <- ulapply(object$mix, summary, link = link, ...)
  paste0("mixture(", paste0(families, collapse = ", "), ")")
}

is.family <- function(x) {
  inherits(x, "family")
}

is.mixfamily <- function(x) {
  inherits(x, "mixfamily")
}

is_linear <- function(family) {
  # indicate if family is for a linear model
  any(family_names(family) %in% c("gaussian", "student", "cauchy"))
}

is_binary <- function(family) {
  # indicate if family is bernoulli or binomial
  any(family_names(family) %in% c("binomial", "bernoulli"))
}

is_ordinal <- function(family) {
  # indicate if family is for an ordinal model
  any(family_names(family) %in% c("cumulative", "cratio", "sratio", "acat"))
}

is_categorical <- function(family) {
  any(family_names(family) %in% "categorical")
}

is_skewed <- function(family) {
  # indicate if family is for model with postive skewed response
  any(family_names(family) %in% 
        c("gamma", "weibull", "exponential", "frechet"))
}

is_lognormal <- function(family) {
  # indicate if family is lognormal
  any(family_names(family) %in% c("lognormal"))
}

is_exgaussian <- function(family) {
  # indicate if family is exgaussian
  any(family_names(family) %in% c("exgaussian"))
}

is_wiener <- function(family) {
  # indicate if family is the wiener diffusion model
  any(family_names(family) %in% c("wiener"))
}

is_asym_laplace <- function(family) {
  # indicates if family is asymmetric laplace
  any(family_names(family) %in% c("asym_laplace"))
}

is_gev <- function(family) {
  # indicates if family is generalized extreme value
  any(family_names(family) %in% c("gen_extreme_value"))
}

is_count <- function(family) {
  # indicate if family is for a count model
  any(family_names(family) %in% c("poisson", "negbinomial", "geometric"))
}

is_hurdle <- function(family, zi_beta = TRUE) {
  # indicate if family is for a hurdle model
  # zi_beta is technically a hurdle model
  any(family_names(family) %in% 
        c("hurdle_poisson", "hurdle_negbinomial", "hurdle_gamma",
          "hurdle_lognormal", if (zi_beta) "zero_inflated_beta"))
}

is_zero_inflated <- function(family, zi_beta = FALSE) {
  # indicate if family is for a zero inflated model
  # zi_beta is technically a hurdle model
  any(family_names(family) %in% 
        c("zero_inflated_poisson", "zero_inflated_negbinomial",
          "zero_inflated_binomial", if (zi_beta) "zero_inflated_beta"))
}

is_zero_one_inflated <- function(family, zi_beta = FALSE) {
  # indicate if family is for a zero one inflated model
  any(family_names(family) %in% "zero_one_inflated_beta")
}

is_2PL <- function(family) {
  # do not remove to provide an informative error message
  # why the special 2PL implementation is not supported anymore
  if (!is(family, "brmsfamily")) {
    out <- FALSE
  } else {
    out <- any(family_names(family) %in% "bernoulli") && 
      identical(family$type, "2PL")
  }
  if (out) {
    stop2("The special implementation of 2PL models has been removed.\n",
          "You can now use argument 'nonlinear' to fit such models.")
  }
  out
}

is_forked <- function(family) {
  # indicate if family has two separate model parts
  is_hurdle(family) || is_zero_inflated(family) || is_2PL(family)
}

is_mv <- function(family, response = NULL) {
  # indicate if the model uses multiple responses
  nresp <- length(response)
  is_mv <- nresp > 1L && is_linear(family) || is_categorical(family) || 
           nresp == 2L && is_forked(family)
  if (nresp > 1L && !is_mv) {
    stop2("Invalid multivariate model")
  }
  is_mv
}

use_real <- function(family) {
  # indicate if family uses real responses
  families <- family_names(family)
  is_linear(families) || is_skewed(families) || 
    any(families %in% 
      c("lognormal", "exgaussian", "inverse.gaussian", "beta", 
        "von_mises", "zero_inflated_beta", "hurdle_gamma", 
        "hurdle_lognormal", "wiener", "asym_laplace", "skew_normal",
        "gen_extreme_value", "zero_one_inflated_beta")
    )
}

use_int <- function(family) {
  # indicate if family uses integer responses
  families <- family_names(family)
  is_binary(families) || has_cat(families) || 
    is_count(families) || is_zero_inflated(families) || 
    any(families %in% c("hurdle_poisson", "hurdle_negbinomial"))
}

has_trials <- function(family) {
  # indicate if family makes use of argument trials
  any(family_names(family) %in% c("binomial", "zero_inflated_binomial"))
}

has_cat <- function(family) {
  # indicate if family makes use of argument cat
  families <- family_names(family)
  is_categorical(families) || is_ordinal(families)
}

has_shape <- function(family) {
  # indicate if family needs a shape parameter
  any(family_names(family) %in% 
        c("gamma", "weibull", "inverse.gaussian", 
          "negbinomial", "hurdle_negbinomial", 
          "hurdle_gamma", "zero_inflated_negbinomial"))
}

has_nu <- function(family) {
  # indicate if family needs a nu parameter
  any(family_names(family) %in% c("student", "frechet"))
}

has_phi <- function(family) {
  # indicate if family needs a phi parameter
  any(family_names(family) %in% 
        c("beta", "zero_inflated_beta", "zero_one_inflated_beta"))
}

has_kappa <- function(family) {
  # indicate if family needs a kappa parameter
  any(family_names(family) %in% c("von_mises"))
}

has_beta <- function(family) {
  # indicate if family needs a beta parameter
  any(family_names(family) %in% c("exgaussian"))
}

has_alpha <- function(family) {
  # indicate if family needs an alpha parameter
  any(family_names(family) %in% c("skew_normal"))
}

has_xi <- function(family) {
  # indicate if family needs a xi parameter
  any(family_names(family) %in% c("gen_extreme_value"))
}

has_sigma <- function(family, bterms = NULL, incmv = FALSE) {
  # indicate if the model needs a sigma parameter
  # Args:
  #  family: model family
  #  bterms: object of class brmsterms
  #  incmv: should MV (linear) models be treated as having sigma? 
  families <- family_names(family)
  is_ln_eg <- any(families %in% 
    c("lognormal", "hurdle_lognormal", "exgaussian",
      "asym_laplace", "gen_extreme_value", "skew_normal")
  )
  if (is.formula(bterms$adforms$se)) {
    # call .se without evaluating the x argument 
    cl <- rhs(bterms$adforms$se)[[2]]
    cl[[1]] <- quote(resp_se_no_data)
    se_only <- isFALSE(attr(eval(cl), "sigma")) 
    if (se_only && use_cov(bterms$autocor)) {
      stop2("Please set argument 'sigma' of function 'se' ",  
            "to TRUE when modeling ARMA covariance matrices.")
    }
  } else {
    se_only <- FALSE
  }
  out <- (is_linear(family) || is_ln_eg) && 
          !se_only && !is(bterms$autocor, "cov_fixed")
  if (!incmv) {
    is_multi <- is_linear(family) && length(bterms$response) > 1L
    out <- out && !is_multi
  }
  out
}

allows_cs <- function(family) {
  # checks if category specific effects are allowed
  all(family_names(family) %in% c("sratio", "cratio", "acat"))
}

is_old_lognormal <- function(family, link = "identity", nresp = 1L,
                             version = utils::packageVersion("brms")) {
  # indicate transformation to lognormal models
  # Args:
  #   link: A character string; ignored if family is of class family
  #   nresp: number of response variables
  #   version: brms version with which the model was fitted
  if (is.family(family)) {
    link <- family$link
    family <- family$family
  }
  isTRUE(family %in% "gaussian") && link == "log" && nresp == 1L &&
    (is.null(version) || version <= "0.9.1")
}

is_old_categorical <- function(x) {
  # indicate if the model is and old categorical model
  stopifnot(is.brmsfit(x))
  if (is(x$fit, "stanfit") && is_categorical(x$family)) {
    if ("bp" %in% x$fit@model_pars) {
      # fitted with brms <= 0.8.0
      out <- 1L
    } else if (is_old_mv(x)) {
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

is_old_mv <- function(x) {
  # indicate if the model uses the old multivariate syntax 
  # from brms < 1.0.0
  stopifnot(is.brmsfit(x))
  bterms <- parse_bf(formula(x), family = family(x))
  if (is.package_version(x$version)) {
    x$version <- list(brms = x$version)
  }
  (is.null(x$version) || x$version$brms <= "0.10.0.9000") &&
    (is_mv(family(x), bterms$response) || is_forked(family(x)))
}
