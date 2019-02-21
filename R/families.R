#' Special Family Functions for \pkg{brms} Models
#' 
#' Family objects provide a convenient way to specify the details of the models 
#' used by many model fitting functions. The family functions presented here are 
#' for use with \pkg{brms} only and will **not** work with other model 
#' fitting functions such as \code{glm} or \code{glmer}. 
#' However, the standard family functions as described in
#' \code{\link[stats:family]{family}} will work with \pkg{brms}.
#' You can also specify custom families for use in \pkg{brms} with
#' the \code{\link{custom_family}} function.
#' 
#' @param family A character string naming the distribution of the response
#'   variable be used in the model. Currently, the following families are
#'   supported: \code{gaussian}, \code{student}, \code{binomial},
#'   \code{bernoulli}, \code{poisson}, \code{negbinomial}, \code{geometric},
#'   \code{discrete_weibull}, \code{Gamma}, \code{skew_normal}, \code{lognormal},
#'   \code{shifted_lognormal}, \code{exgaussian}, \code{wiener},
#'   \code{inverse.gaussian}, \code{exponential}, \code{weibull},
#'   \code{frechet}, \code{Beta}, \code{dirichlet}, \code{von_mises},
#'   \code{asym_laplace}, \code{gen_extreme_value}, \code{categorical},
#'   \code{multinomial}, \code{cumulative}, \code{cratio}, \code{sratio},
#'   \code{acat}, \code{hurdle_poisson}, \code{hurdle_negbinomial},
#'   \code{hurdle_gamma}, \code{hurdle_lognormal},
#'   \code{zero_inflated_binomial}, \code{zero_inflated_beta},
#'   \code{zero_inflated_negbinomial}, \code{zero_inflated_poisson}, and
#'   \code{zero_one_inflated_beta}.
#' @param link A specification for the model link function. This can be a
#'   name/expression or character string. See the 'Details' section for more
#'   information on link functions supported by each family.
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
#' @param threshold A character string indicating the type 
#'   of thresholds (i.e. intercepts) used in an ordinal model. 
#'   \code{"flexible"} provides the standard unstructured thresholds and 
#'   \code{"equidistant"} restricts the distance between 
#'   consecutive thresholds to the same value.
#' @param refcat Optional name of the reference response category used in
#'   categorical, multinomial, and dirichlet models. If \code{NULL} (the
#'   default), the first category is used as the reference. If \code{NA}, all
#'   categories will be predicted, which requires strong priors or carefully
#'   specified predictor terms in order to lead to an identified model.
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
#'   Family \code{exgaussian} ('exponentially modified Gaussian') and 
#'   \code{shifted_lognormal} are especially suited to model reaction times.
#'   The \code{wiener} family provides an implementation of the Wiener 
#'   diffusion model. For this family, the main formula predicts the drift 
#'   parameter 'delta' and all other parameters are modeled as auxiliary parameters 
#'   (see \code{\link{brmsformula}} for details).
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
#'   \code{\link[stats:family]{family}},
#'   \code{\link{customfamily}}
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
                       link_quantile = "logit",
                       threshold = c("flexible", "equidistant"),
                       refcat = NULL) {
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
    link_quantile = link_quantile,
    threshold = threshold, refcat = refcat
  )
}

.brmsfamily <- function(family, link = NULL, slink = link,
                        threshold = c("flexible", "equidistant"),
                        refcat = NULL, ...) {
  # helper function to prepare brmsfamily objects
  # Args:
  #   family: character string naming the model family
  #   link: character string naming the link function
  #   slink: can be used with substitute(link) for 
  #          non-standard evaluation of the link function
  #   threshold: threshold type for ordinal models
  #   ...: link functions (as character strings) of parameters
  # Returns:
  #   An object of class = c('brmsfamily', 'family') to be used
  #   only inside the brms package
  family <- tolower(as_one_character(family))
  aux_links <- list(...)
  pattern <- c("^normal$", "^zi_", "^hu_")
  replacement <- c("gaussian", "zero_inflated_", "hurdle_")
  family <- rename(family, pattern, replacement, fixed = FALSE)
  ok_families <- lsp("brms", pattern = "^\\.family_")
  ok_families <- sub("^\\.family_", "", ok_families)
  if (!family %in% ok_families) {
    stop2(family, " is not a supported family. Supported ", 
          "families are:\n", collapse_comma(ok_families))
  }
  family_info <- get(paste0(".family_", family))()
  ok_links <- family_info$links
  family_info$links <- NULL
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
  out <- list(
    family = family, link = slink,
    linkfun = function(mu) link(mu, link = slink),
    linkinv = function(eta) ilink(eta, link = slink)
  )
  out[names(family_info)] <- family_info
  class(out) <- c("brmsfamily", "family")
  for (dp in valid_dpars(out)) {
    alink <- as.character(aux_links[[paste0("link_", dp)]])
    if (length(alink)) {
      alink <- as_one_character(alink)
      valid_links <- links_dpars(dp)
      if (!alink %in% valid_links) {
        stop2(
          "'", alink, "' is not a supported link ", 
          "for parameter '", dp, "'.\nSupported links are: ", 
          collapse_comma(valid_links)
        )
      }
      out[[paste0("link_", dp)]] <- alink
    }
  }
  if (is_ordinal(out$family)) {
    out$threshold <- match.arg(threshold)
  }
  if (conv_cats_dpars(out$family)) {
    if (!is.null(refcat)) {
      out$refcat <- as_one_character(refcat, allow_na = TRUE) 
    }
  }
  out
}

check_family <- function(family, link = NULL, threshold = NULL) {
  # checks and corrects validity of the model family
  # Args:
  #   family: Either a function, an object of class 'family' 
  #           or a character string of length one or two
  #   link: an optional character string naming the link function
  #         ignored if family is a function or a family object
  #   threshold: optional character string specifying the threshold
  #              type in ordinal models
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
  if (is_ordinal(family) && !is.null(threshold)) {
    # slot 'threshold' deprecated as of brms > 1.7.0
    threshold <- match.arg(threshold, c("flexible", "equidistant"))
    family$threshold <- threshold
  }
  family
}

family_info <- function(x, y, ...) {
  # extract special information of families
  # Args: 
  #   x: object from which to extract
  #   y: name of the component to extract
  UseMethod("family_info")
}

#' @export
family_info.default <- function(x, y, ...) {
  x <- as.character(x)
  ulapply(x, .family_info, y = y, ...)
}

.family_info <- function(x, y, ...) {
  x <- as_one_character(x)
  y <- as_one_character(y)
  if (y == "family") {
    return(x)
  }
  if (!nzchar(x)) {
    return(NULL)
  }
  info <- get(paste0(".family_", x))()
  if (y == "link") {
    out <- info$links[1]  # default link
  } else {
    info$links <- NULL
    out <- info[[y]] 
  }
  out
}

family_info.NULL <- function(x, y, ...) {
  NULL
}

#' @export
family_info.family <- function(x, y, ...) {
  family_info(x$family, y = y, ...)
}

#' @export
family_info.brmsfamily <- function(x, y, ...) {
  y <- as_one_character(y)
  out <- x[[y]]
  if (is.null(out)) {
    # required for models fitted with brms 2.2 or earlier
    out <- family_info(x$family, y = y, ...)
  }
  out
}

#' @export
family_info.mixfamily <- function(x, y, ...) {
  out <- lapply(x$mix, family_info, y = y, ...)
  combine_family_info(out, y = y)
}

#' @export
family_info.brmsformula <- function(x, y, ...) {
  family_info(x$family, y = y, ...)
}

#' @export
family_info.mvbrmsformula <- function(x, y, ...) {
  out <- lapply(x$forms, family_info, y = y, ...)
  combine_family_info(out, y = y)
}

#' @export
family_info.brmsterms <- function(x, y, ...) {
  family_info(x$family, y = y, ...)
}

#' @export
family_info.mvbrmsterms <- function(x, y, ...) {
  out <- lapply(x$terms, family_info, y = y, ...)
  combine_family_info(out, y = y)
}

#' @export
family_info.btl <- function(x, y, ...) {
  family_info(x$family, y = y, ...)
}

#' @export
family_info.btnl <- function(x, y, ...) {
  family_info(x$family, y = y, ...)
}

#' @export
family_info.brmsfit <- function(x, y, ...) {
  family_info(x$formula, y = y, ...)
}

combine_family_info <- function(x, y, ...) {
  # combine information from multiple families
  # provides special handling for certain elements
  y <- as_one_character(y)
  unite <- c("dpars", "type", "specials", "include", "const")
  if (y %in% c("family", "link")) {
    x <- unlist(x)
  } else if (y %in% unite) {
    x <- Reduce("union", x)
  } else if (y == "ad") {
    x <- Reduce("intersect", x)
  } else if (y == "ybounds") {
    x <- do_call(rbind, x)
    x <- c(max(x[, 1]), min(x[, 2]))
  } else if (y == "closed") {
    # closed only if no bounds are open
    x <- do_call(rbind, x)
    clb <- !any(ulapply(x[, 1], isFALSE))
    cub <- !any(ulapply(x[, 2], isFALSE))
    x <- c(clb, cub)
  }
  x
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
discrete_weibull <- function(link = "logit", link_shape = "log") {
  slink <- substitute(link)
  .brmsfamily("discrete_weibull", link = link, slink = slink,
              link_shape = link_shape)
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
shifted_lognormal <- function(link = "identity", link_sigma = "log",
                              link_ndt = "log") {
  slink <- substitute(link)
  .brmsfamily("shifted_lognormal", link = link, slink = slink,
              link_sigma = link_sigma, link_ndt = link_ndt)
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
  .brmsfamily("wiener", link = link, slink = slink,
              link_bs = link_bs, link_ndt = link_ndt,
              link_bias = link_bias)
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
dirichlet <- function(link = "logit", link_phi = "log", refcat = NULL) {
  slink <- substitute(link)
  .brmsfamily("dirichlet", link = link, slink = slink,
              link_phi = link_phi, refcat = refcat)
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
categorical <- function(link = "logit", refcat = NULL) {
  slink <- substitute(link)
  .brmsfamily("categorical", link = link, slink = slink, refcat = refcat)
}

#' @rdname brmsfamily
#' @export
multinomial <- function(link = "logit", refcat = NULL) {
  slink <- substitute(link)
  .brmsfamily("multinomial", link = link, slink = slink, refcat = refcat)
}

#' @rdname brmsfamily
#' @export
cumulative <- function(link = "logit", link_disc = "log",
                       threshold = c("flexible", "equidistant")) {
  slink <- substitute(link)
  .brmsfamily("cumulative", link = link, slink = slink,
              link_disc = link_disc, threshold = threshold)
}

#' @rdname brmsfamily
#' @export
sratio <- function(link = "logit", link_disc = "log",
                   threshold = c("flexible", "equidistant")) {
  slink <- substitute(link)
  .brmsfamily("sratio", link = link, slink = slink,
              link_disc = link_disc, threshold = threshold)
}

#' @rdname brmsfamily
#' @export
cratio <- function(link = "logit", link_disc = "log",
                   threshold = c("flexible", "equidistant")) {
  slink <- substitute(link)
  .brmsfamily("cratio", link = link, slink = slink,
              link_disc = link_disc, threshold = threshold)
}

#' @rdname brmsfamily
#' @export
acat <- function(link = "logit", link_disc = "log",
                 threshold = c("flexible", "equidistant")) {
  slink <- substitute(link)
  .brmsfamily("acat", link = link, slink = slink,
              link_disc = link_disc, threshold = threshold)
}

#' Finite Mixture Families in \pkg{brms}
#' 
#' Set up a finite mixture family for use in \pkg{brms}.
#' 
#' @param ... One or more objects providing a description of the 
#'   response distributions to be combined in the mixture model. 
#'   These can be family functions, calls to family functions or 
#'   character strings naming the families. For details of supported 
#'   families see \code{\link{brmsfamily}}.
#' @param flist Optional list of objects, which are treated in the 
#'   same way as objects passed via the \code{...} argument.
#' @param nmix Optional numeric vector specifying the number of times
#'   each family is repeated. If specified, it must have the same length 
#'   as the number of families passed via \code{...} and \code{flist}.
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
#' intercepts via \code{\link{set_prior}} to improve convergence. 
#' In addition, it is sometimes necessary to set \code{inits = 0} in the call to 
#' \code{\link{brm}} to allow chains to initialize properly.
#' 
#' For more details on the specification of mixture
#' models, see \code{\link{brmsformula}}.
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
#'   prior(normal(0, 7), Intercept, dpar = mu1),
#'   prior(normal(5, 7), Intercept, dpar = mu2)
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
  if (length(family$mix) < 2L) {
    stop2("Expecting at least 2 mixture components.")
  }
  if (use_real(family) && use_int(family)) {
    stop2("Cannot mix families with real and integer support.")
  }
  is_ordinal <- ulapply(family$mix, is_ordinal)
  if (any(is_ordinal) && any(!is_ordinal)) {
    stop2("Cannot mix ordinal and non-ordinal families.")
  }
  no_mixture <- ulapply(family$mix, no_mixture)
  if (any(no_mixture)) {
    stop2("Some of the families are not allowed in mixture models.")
  }
  for (fam in family$mix) {
    if (is.customfamily(fam) && "theta" %in% fam$dpars) {
      stop2("Parameter name 'theta' is reserved in mixture models.")
    }
  }
  if (is.null(order)) {
    if (any(is_ordinal)) {
      family$order <- "none"
      message("Setting order = 'none' for mixtures of ordinal families.")
    } else if (length(unique(family_names(family))) == 1L) {
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
    if (any(is_ordinal) && family$order != "none") {
      stop2("Ordinal mixture models only support order = 'none'.")
    }
  }
  family
}

#' Custom Families in \pkg{brms} Models
#' 
#' Define custom families (i.e. response distribution) for use in 
#' \pkg{brms} models. It allows users to benefit from the modeling 
#' flexibility of \pkg{brms}, while applying their self-defined likelihood
#' functions. All of the post-processing methods for \code{brmsfit} 
#' objects can be made compatible with custom families. 
#' See \code{vignette("brms_customfamilies")} for more details.
#' For a list of built-in families see \code{\link{brmsfamily}}.
#' 
#' @aliases customfamily
#' 
#' @param name Name of the custom family.
#' @param dpars Names of the distributional parameters of
#'   the family. One parameter must be named \code{"mu"} and
#'   the main formula of the model will correspond to that
#'   parameter.
#' @param links Names of the link functions of the 
#'   distributional parameters.
#' @param type Indicates if the response distribution is
#'   continuous (\code{"real"}) or discrete (\code{"int"}).
#' @param lb Vector of lower bounds of the distributional 
#'   parameters. Defaults to \code{NA} that is no lower bound.
#' @param ub Vector of upper bounds of the distributional 
#'   parameters. Defaults to \code{NA} that is no upper bound.
#' @param vars Names of variables, which are part of the likelihood
#'   function without being distributional parameters. That is,
#'   \code{vars} can be used to pass data to the likelihood. 
#'   See \code{\link{stanvar}} for details about adding self-defined
#'   data to the generated \pkg{Stan} model.
#' @param specials A character vector of special options to enable
#'   for this custom family. Currently for internal use only.
#' @param log_lik Optional function to compute log-likelihood values of
#'   the model in \R. This is only relevant if one wants to ensure 
#'   compatibility with method \code{\link[brms:log_lik.brmsfit]{log_lik}}.
#' @param predict Optional function to compute predicted values of
#'   the model in \R. This is only relevant if one wants to ensure 
#'   compatibility with method \code{\link[brms:predict.brmsfit]{predict}}.  
#' @param fitted Optional function to compute fitted values of
#'   the model in \R. This is only relevant if one wants to ensure 
#'   compatibility with method \code{\link[brms:fitted.brmsfit]{fitted}}.     
#' @param env An \code{\link{environment}} in which certain post-processing 
#'   functions related to the custom family can be found, if there were not 
#'   directly passed to \code{custom_family}. This is only
#'   relevant if one wants to ensure compatibility with the methods
#'   \code{\link[brms:predict.brmsfit]{predict}}, 
#'   \code{\link[brms:fitted.brmsfit]{fitted}}, or
#'   \code{\link[brms:log_lik.brmsfit]{log_lik}}.
#'   By default, \code{env} is the enviroment from which 
#'   \code{custom_family} is called.
#'   
#' @details The corresponding probability density or mass \code{Stan} 
#'   functions need to have the same name as the custom family.
#'   That is if a family is called \code{myfamily}, then the 
#'   \pkg{Stan} functions should be called \code{myfamily_lpdf} or
#'   \code{myfamily_lpmf} depending on whether it defines a 
#'   continuous or discrete distribution.
#'   
#' @return An object of class \code{customfamily} inheriting
#'   from class \code{\link{brmsfamily}}.
#'   
#' @seealso \code{\link{brmsfamily}}, \code{\link{stanvar}}
#' 
#' @examples
#' \dontrun{
#' ## demonstrate how to fit a beta-binomial model
#' ## generate some fake data
#' phi <- 0.7
#' n <- 300
#' z <- rnorm(n, sd = 0.2)
#' ntrials <- sample(1:10, n, replace = TRUE)
#' eta <- 1 + z
#' mu <- exp(eta) / (1 + exp(eta))
#' a <- mu * phi
#' b <- (1 - mu) * phi
#' p <- rbeta(n, a, b)
#' y <- rbinom(n, ntrials, p)
#' dat <- data.frame(y, z, ntrials)
#' 
#' # define a custom family
#' beta_binomial2 <- custom_family(
#'   "beta_binomial2", dpars = c("mu", "phi"),
#'   links = c("logit", "log"), lb = c(NA, 0),
#'   type = "int", vars = "trials[n]"
#' )
#' 
#' # define the corresponding Stan density function
#' stan_funs <- "
#'   real beta_binomial2_lpmf(int y, real mu, real phi, int N) {
#'     return beta_binomial_lpmf(y | N, mu * phi, (1 - mu) * phi);
#'   }
#' "
#' 
#' # fit the model
#' fit <- brm(y | trials(ntrials) ~ z, data = dat, 
#'            family = beta_binomial2, stan_funs = stan_funs)
#' summary(fit)
#' }
#' 
#' @export
custom_family <- function(name, dpars = "mu", links = "identity",
                          type = c("real", "int"), lb = NA, ub = NA,
                          vars = NULL, specials = NULL, 
                          log_lik = NULL, predict = NULL, 
                          fitted = NULL, env = parent.frame()) {
  name <- as_one_character(name)
  dpars <- as.character(dpars)
  links <- as.character(links)
  type <- match.arg(type)
  lb <- as.character(lb)
  ub <- as.character(ub)
  vars <- as.character(vars)
  specials <- as.character(specials)
  env <- as.environment(env)
  if (any(duplicated(dpars))) {
    stop2("Duplicated 'dpars' are not allowed.")
  }
  if (!"mu" %in% dpars) {
    stop2("All families must have a 'mu' parameter.")
  }
  if (any(grepl("_|\\.", dpars))) {
    stop2("Dots or underscores are not allowed in 'dpars'.")
  }
  if (any(grepl("[[:digit:]]+$", dpars))) {
    stop2("'dpars' should not end with a number.")
  }
  for (arg in c("links", "lb", "ub")) {
    obj <- get(arg)
    if (length(obj) == 1L) {
      obj <- rep(obj, length(dpars))
      assign(arg, obj)
    }
    if (length(dpars) != length(obj)) {
      stop2("'", arg, "' must be of the same length as 'dpars'.")
    }
  }
  if (!is.null(log_lik)) {
    log_lik <- as.function(log_lik)
    args <- names(formals(log_lik))
    if (!is_equal(args[1:2], c("i", "draws"))) {
      stop2("The first two arguments of 'log_lik' ", 
            "should be 'i' and 'draws'.")
    }
  }
  if (!is.null(predict)) {
    predict <- as.function(predict)
    args <- names(formals(predict))
    if (!is_equal(args[1:3], c("i", "draws", "..."))) {
      stop2("The first three arguments of 'predict' ", 
            "should be 'i', 'draws', and '...'.")
    }
  }
  if (!is.null(fitted)) {
    fitted <- as.function(fitted)
    args <- names(formals(fitted))
    if (!is_equal(args[1], "draws")) {
      stop2("The first argument of 'fitted' should be 'draws'.")
    }
  }
  lb <- named_list(dpars, lb)
  ub <- named_list(dpars, ub)
  is_mu <- "mu" == dpars
  link <- links[is_mu]
  out <- nlist(
    family = "custom", link, name, 
    dpars, lb, ub, type, vars, specials,
    log_lik, predict, fitted, env
  )
  if (length(dpars) > 1L) {
    out[paste0("link_", dpars[!is_mu])] <- links[!is_mu]
  }
  structure(out, class = c("customfamily", "brmsfamily", "family"))
}

valid_dpars <- function(family, ...) {
  # get valid distributional parameters for a family
  UseMethod("valid_dpars")
}

#' @export
valid_dpars.default <- function(family, ...) {
  if (!length(family)) {
    return("mu")
  }
  family <- check_family(family) 
  family_info(family, "dpars", ...)
}

#' @export
valid_dpars.mixfamily <- function(family, ...) {
  out <- lapply(family$mix, valid_dpars, ...)
  for (i in seq_along(out)) {
    out[[i]] <- paste0(out[[i]], i)
  }
  c(unlist(out), paste0("theta", seq_along(out)))
}

#' @export
valid_dpars.brmsformula <- function(family, ...) {
  valid_dpars(family$family, ...)
}

#' @export
valid_dpars.mvbrmsformula <- function(family, ...) {
  ulapply(family$forms, valid_dpars, ...)
}

#' @export
valid_dpars.brmsterms <- function(family, ...) {
  valid_dpars(family$family, ...)
}

#' @export
valid_dpars.mvbrmsterms <- function(family, ...) {
  ulapply(family$terms, valid_dpars, ...)
}

#' @export
valid_dpars.brmsfit <- function(family, ...) {
  valid_dpars(family$formula, ...)
}

dpar_class <- function(dpar) {
  # class of a distributional parameter
  sub("[[:digit:]]*$", "", dpar)
}

dpar_id <- function(dpar) {
  # id of a distributional parameter
  out <- get_matches("[[:digit:]]+$", dpar, simplify = FALSE)
  ulapply(out, function(x) ifelse(length(x), x, ""))
}

links_dpars <- function(dpar) {
  # link functions for distributional parameters
  if (!length(dpar)) dpar <- ""
  switch(dpar,
    character(0),
    mu = "identity",  # not actually used
    sigma = c("log", "identity"), 
    shape = c("log", "identity"),
    nu = c("logm1", "identity"), 
    phi = c("log", "identity"),
    kappa = c("log", "identity"), 
    beta = c("log", "identity"),
    zi = c("logit", "identity"), 
    hu = c("logit", "identity"),
    zoi = c("logit", "identity"), 
    coi = c("logit", "identity"), 
    disc = c("log", "identity"),
    bs = c("log", "identity"), 
    ndt = c("log", "identity"),
    bias = c("logit", "identity"),
    quantile = c("logit", "identity"),
    xi = c("log1p", "identity"),
    alpha = c("identity", "log"),
    theta = c("identity")
  )
}

dpar_family <- function(family, dpar, ...) {
  # generate a family object of a distributional parameter
  UseMethod("dpar_family")
}

#' @export
dpar_family.default <- function(family, dpar, ...) {
  dp_class <- dpar_class(dpar)
  if (dp_class != "mu" || conv_cats_dpars(family)) {
    link <- family[[paste0("link_", dp_class)]]
    family <- .dpar_family(dpar, link)
  }
  family
}

#' @export
dpar_family.mixfamily <- function(family, dpar, ...) {
  dp_id <- as.numeric(dpar_id(dpar))
  if (!(length(dp_id) == 1L && is.numeric(dp_id))) {
    stop2("Parameter '", dpar, "' is not a valid mixture parameter.")
  }
  dpar_family(family$mix[[dp_id]], dpar, ...)
}

.dpar_family <- function(dpar = NULL, link = NULL) {
  # set up special family objects for distributional parameters
  # Args:
  #   dpar: name of the distributional parameter
  #   link: optional link function of the parameter
  links <- links_dpars(dpar_class(dpar))
  if (!length(link)) {
    if (!length(links)) {
      link <- "identity" 
    } else {
      link <- links[1]
    }
  }
  link <- as_one_character(link)
  structure(
    nlist(family = "", link, dpar),
    class = c("brmsfamily", "family")
  )
}

#' @export
print.brmsfamily <- function(x, links = FALSE, newline = TRUE, ...) {
  cat("\nFamily:", x$family, "\n")
  cat("Link function:", x$link, "\n")
  if (!is.null(x$threshold)) {
    cat("Threshold:", x$threshold, "\n")
  }
  if (isTRUE(links) || is.character(links)) {
    dp_links <- x[grepl("^link_", names(x))]
    names(dp_links) <- sub("^link_", "", names(dp_links))
    if (is.character(links)) {
      dp_links <- rmNULL(dp_links[links])
    }
    for (dp in names(dp_links)) {
      cat(paste0(
        "Link function of '", dp, "' (if predicted): ", 
        dp_links[[dp]], "\n"
      )) 
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

#' @export
print.customfamily <- function(x, links = FALSE, newline = TRUE, ...) {
  cat("\nCustom family:", x$name, "\n")
  cat("Link function:", x$link, "\n")
  cat("Parameters:", paste0(x$dpars, collapse = ", "), "\n")
  if (isTRUE(links) || is.character(links)) {
    dp_links <- x[grepl("^link_", names(x))]
    names(dp_links) <- sub("^link_", "", names(dp_links))
    if (is.character(links)) {
      dp_links <- rmNULL(dp_links[links])
    }
    for (dp in names(dp_links)) {
      cat(paste0(
        "Link function of '", dp, "' (if predicted): ", 
        dp_links[[dp]], "\n"
      )) 
    }
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

#' @method summary customfamily
#' @export
summary.customfamily <- function(object, link = TRUE, ...) {
  object$family <- object$name
  summary.family(object, link = link, ...)
}

summarise_families <- function(x) {
  # summary of families used in summary.brmsfit
  UseMethod("summarise_families")
}

#' @export
summarise_families.mvbrmsformula <- function(x, ...) {
  out <- ulapply(x$forms, summarise_families, ...)
  paste0("MV(", paste0(out, collapse = ", "), ")")
}

#' @export
summarise_families.brmsformula <- function(x, ...) {
  summary(x$family, link = FALSE, ...)
}

summarise_links <- function(x, ...) {
  # summary of link functions used in summary.brmsfit
  UseMethod("summarise_links")
}

#' @export
summarise_links.mvbrmsformula <- function(x, wsp = 0, ...) {
  str_wsp <- collapse(rep(" ", wsp))
  links <- ulapply(x$forms, summarise_links, mv = TRUE, ...)
  paste0(links, collapse = paste0("\n", str_wsp))
}

#' @export
summarise_links.brmsformula <- function(x, mv = FALSE, ...) {
  x <- parse_bf(x)
  dpars <- valid_dpars(x)
  links <- setNames(rep("identity", length(dpars)), dpars)
  links_pred <- ulapply(x$dpars, function(x) x$family$link)
  links[names(links_pred)] <- links_pred
  if (conv_cats_dpars(x)) {
    links[grepl("^mu", names(links))] <- x$family$link
  }
  resp <- if (mv) usc(combine_prefix(x))
  names(links) <- paste0(names(links), resp)
  paste0(names(links), " = ", links, collapse = "; ")
}

is.family <- function(x) {
  inherits(x, "family")
}

is.brmsfamily <- function(x) {
  inherits(x, "brmsfamily")
}

is.mixfamily <- function(x) {
  inherits(x, "mixfamily")
}

is.customfamily <- function(x) {
  inherits(x, "customfamily")
}

family_names <- function(x) {
  family_info(x, "family")
}

use_real <- function(family) {
  # indicate if family uses real responses
  "real" %in% family_info(family, "type")
}

use_int <- function(family) {
  # indicate if family uses integer responses
  "int" %in% family_info(family, "type")
}

is_binary <- function(family) {
  "binary" %in% family_info(family, "specials")
}

is_categorical <- function(family) {
  "categorical" %in% family_info(family, "specials")
}

is_ordinal <- function(family) {
  "ordinal" %in% family_info(family, "specials")
}

is_multinomial <- function(family) {
  "multinomial" %in% family_info(family, "specials")
}

is_dirichlet <- function(family) {
  "dirichlet" %in% family_info(family, "specials")
}

allow_factors <- function(family) {
  specials <- c("binary", "categorical", "ordinal")
  any(specials %in% family_info(family, "specials"))
}

allow_autocor <- function(family) {
  # checks if autocorrelation structures are allowed
  "autocor" %in% family_info(family, "specials")
}

allow_cs <- function(family) {
  # checks if category specific effects are allowed
  "cs" %in% family_info(family, "specials")
}

conv_cats_dpars <- function(family) {
  # choose dpar names based on categories?
  is_categorical(family) || is_multinomial(family) || is_dirichlet(family)
}

no_mixture <- function(family) {
  # families not allowed in mixture models
  is_categorical(family) || is_multinomial(family) || is_dirichlet(family)
}

has_multicol <- function(family) {
  # indicate if the response should consist of multiple columns
  is_multinomial(family) || is_dirichlet(family)
}

has_logscale <- function(family) {
  # indicate if the response is modeled on the log-scale
  # even if formally the link function is not 'log'
  "logscale" %in% family_info(family, "specials")
}

has_trials <- function(family) {
  # indicate if family makes use of argument trials
  "trials" %in% family_info(family, "ad") &&
    !"custom" %in% family_names(family)
}

has_cat <- function(family) {
  # indicate if family has more than two response categories
  is_categorical(family) || is_ordinal(family) || 
    is_multinomial(family) || is_dirichlet(family)
}

has_ndt <- function(family) {
  "ndt" %in% dpar_class(family_info(family, "dpars"))
}

has_sigma <- function(family) {
  "sigma" %in% dpar_class(family_info(family, "dpars"))
}

no_sigma <- function(bterms) {
  # check if sigma should be explicitely set to 0
  stopifnot(is.brmsterms(bterms))
  if (is.formula(bterms$adforms$se)) {
    # call resp_se without evaluating the x argument
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
  se_only || is.cor_fixed(bterms$autocor)
}

simple_sigma <- function(bterms) {
  # has the model a non-predicted but estimated sigma parameter?
  stopifnot(is.brmsterms(bterms))
  has_sigma(bterms) && is.null(bterms$dpars$sigma)
}

pred_sigma <- function(bterms) {
  # has the model a predicted sigma parameter?
  stopifnot(is.brmsterms(bterms))
  "sigma" %in% dpar_class(names(bterms$dpars))
}

no_nu <- function(bterms) {
  # the multi_student_t family only has a single 'nu' parameter
  isTRUE(bterms$rescor) && "student" %in% family_names(bterms)
}
