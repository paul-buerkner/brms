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
#'   \code{bernoulli}, \code{beta-binomial}, \code{poisson}, \code{negbinomial},
#'   \code{geometric}, \code{Gamma}, \code{skew_normal}, \code{lognormal},
#'   \code{shifted_lognormal}, \code{exgaussian}, \code{wiener},
#'   \code{inverse.gaussian}, \code{exponential}, \code{weibull},
#'   \code{frechet}, \code{Beta}, \code{dirichlet}, \code{von_mises},
#'   \code{asym_laplace}, \code{gen_extreme_value}, \code{categorical},
#'   \code{multinomial}, \code{cumulative}, \code{cratio}, \code{sratio},
#'   \code{acat}, \code{hurdle_poisson}, \code{hurdle_negbinomial},
#'   \code{hurdle_gamma}, \code{hurdle_lognormal},
#'   \code{zero_inflated_binomial}, \code{zero_inflated_beta_binomial},
#'   \code{zero_inflated_beta}, \code{zero_inflated_negbinomial},
#'   \code{zero_inflated_poisson}, and \code{zero_one_inflated_beta}.
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
#'   \code{"flexible"} provides the standard unstructured thresholds,
#'   \code{"equidistant"} restricts the distance between
#'   consecutive thresholds to the same value, and
#'   \code{"sum_to_zero"} ensures the thresholds sum to zero.
#' @param refcat Optional name of the reference response category used in
#'   \code{categorical}, \code{multinomial}, \code{dirichlet} and
#'   \code{logistic_normal} models. If \code{NULL} (the default), the first
#'   category is used as the reference. If \code{NA}, all categories will be
#'   predicted, which requires strong priors or carefully specified predictor
#'   terms in order to lead to an identified model.
#' @param bhaz Currently for experimental purposes only.
#'
#' @details
#'   Below, we list common use cases for the different families.
#'   This list is not ment to be exhaustive.
#'   \itemize{
#'   \item{Family \code{gaussian} can be used for linear regression.}
#'
#'   \item{Family \code{student} can be used for robust linear regression
#'   that is less influenced by outliers.}
#'
#'   \item{Family \code{skew_normal} can handle skewed responses in linear
#'   regression.}
#'
#'   \item{Families \code{poisson}, \code{negbinomial}, and \code{geometric}
#'   can be used for regression of unbounded count data.}
#'
#'   \item{Families \code{bernoulli}, \code{binomial}, and \code{beta_binomial}
#'   can be used for binary regression (i.e., most commonly logistic
#'   regression).}
#'
#'   \item{Families \code{categorical} and \code{multinomial} can be used for
#'   multi-logistic regression when there are more than two possible outcomes.}
#'
#'   \item{Families \code{cumulative}, \code{cratio} ('continuation ratio'),
#'   \code{sratio} ('stopping ratio'), and \code{acat} ('adjacent category')
#'   leads to ordinal regression.}
#'
#'   \item{Families \code{Gamma}, \code{weibull}, \code{exponential},
#'   \code{lognormal}, \code{frechet}, \code{inverse.gaussian}, and \code{cox}
#'   (Cox proportional hazards model) can be used (among others) for
#'   time-to-event regression also known as survival regression.}
#'
#'   \item{Families \code{weibull}, \code{frechet}, and \code{gen_extreme_value}
#'   ('generalized extreme value') allow for modeling extremes.}
#'
#'   \item{Families \code{beta}, \code{dirichlet}, and \code{logistic_normal}
#'   can be used to model responses representing rates or probabilities.}
#'
#'   \item{Family \code{asym_laplace} allows for quantile regression when fixing
#'   the auxiliary \code{quantile} parameter to the quantile of interest.}
#'
#'   \item{Family \code{exgaussian} ('exponentially modified Gaussian') and
#'   \code{shifted_lognormal} are especially suited to model reaction times.}
#'
#'   \item{Family \code{wiener} provides an implementation of the Wiener
#'   diffusion model. For this family, the main formula predicts the drift
#'   parameter 'delta' and all other parameters are modeled as auxiliary parameters
#'   (see \code{\link{brmsformula}} for details).}
#'
#'   \item{Families \code{hurdle_poisson}, \code{hurdle_negbinomial},
#'   \code{hurdle_gamma}, \code{hurdle_lognormal}, \code{zero_inflated_poisson},
#'   \code{zero_inflated_negbinomial}, \code{zero_inflated_binomial},
#'   \code{zero_inflated_beta_binomial}, \code{zero_inflated_beta}, and
#'   \code{zero_one_inflated_beta} allow to estimate zero-inflated and hurdle
#'   models. These models can be very helpful when there are many zeros in the
#'   data (or ones in case of one-inflated models)
#'   that cannot be explained by the primary distribution of the response.}
#'   }
#'
#'   Below, we list all possible links for each family.
#'   The first link mentioned for each family is the default.
#'   \itemize{
#'   \item{Families \code{gaussian}, \code{student}, \code{skew_normal},
#'   \code{exgaussian}, \code{asym_laplace}, and \code{gen_extreme_value}
#'   support the links (as names) \code{identity}, \code{log}, \code{inverse},
#'   and \code{softplus}.}
#'
#'   \item{Families \code{poisson}, \code{negbinomial}, \code{geometric},
#'   \code{zero_inflated_poisson}, \code{zero_inflated_negbinomial},
#'   \code{hurdle_poisson}, and \code{hurdle_negbinomial} support
#'   \code{log}, \code{identity}, \code{sqrt}, and \code{softplus}.}
#'
#'   \item{Families \code{binomial}, \code{bernoulli}, \code{beta_binomial},
#'   \code{zero_inflated_binomial}, \code{zero_inflated_beta_binomial}, 
#'   \code{Beta}, \code{zero_inflated_beta}, and \code{zero_one_inflated_beta}
#'   support \code{logit}, \code{probit}, \code{probit_approx}, \code{cloglog}, 
#'   \code{cauchit}, \code{identity}, and \code{log}.}
#'
#'   \item{Families \code{cumulative}, \code{cratio}, \code{sratio},
#'   and \code{acat} support \code{logit}, \code{probit},
#'   \code{probit_approx}, \code{cloglog}, and \code{cauchit}.}
#'
#'   \item{Families \code{categorical}, \code{multinomial}, and \code{dirichlet}
#'   support \code{logit}.}
#'
#'   \item{Families \code{Gamma}, \code{weibull}, \code{exponential},
#'   \code{frechet}, and \code{hurdle_gamma} support
#'   \code{log}, \code{identity}, \code{inverse}, and \code{softplus}.}
#'
#'   \item{Families \code{lognormal} and \code{hurdle_lognormal}
#'   support \code{identity} and \code{inverse}.}
#'
#'   \item{Family \code{logistic_normal} supports \code{identity}.}
#'
#'   \item{Family \code{inverse.gaussian} supports \code{1/mu^2},
#'   \code{inverse}, \code{identity}, \code{log}, and \code{softplus}.}
#'
#'   \item{Family \code{von_mises} supports \code{tan_half} and
#'   \code{identity}.}
#'
#'   \item{Family \code{cox} supports \code{log}, \code{identity},
#'   and \code{softplus} for the proportional hazards parameter.}
#'
#'   \item{Family \code{wiener} supports \code{identity}, \code{log},
#'   and \code{softplus} for the main parameter which represents the
#'   drift rate.}
#'   }
#'
#'   Please note that when calling the \code{\link[stats:family]{Gamma}} family
#'   function of the \pkg{stats} package, the default link will be
#'   \code{inverse} instead of \code{log} although the latter is the default in
#'   \pkg{brms}. Also, when using the family functions \code{gaussian},
#'   \code{binomial}, \code{poisson}, and \code{Gamma} of the \pkg{stats}
#'   package (see \code{\link[stats:family]{family}}), special link functions
#'   such as \code{softplus} or \code{cauchit} won't work. In this case, you
#'   have to use \code{brmsfamily} to specify the family with corresponding link
#'   function.
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
                       threshold = "flexible",
                       refcat = NULL, bhaz = NULL) {
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
    threshold = threshold, refcat = refcat,
    bhaz = bhaz
  )
}

# helper function to prepare brmsfamily objects
# @param family character string naming the model family
# @param link character string naming the link function
# @param slink can be used with substitute(link) for
#   non-standard evaluation of the link function
# @param threshold threshold type for ordinal models
# @param ... link functions (as character strings) of parameters
# @return an object of 'brmsfamily' which inherits from 'family'
.brmsfamily <- function(family, link = NULL, slink = link,
                        threshold = "flexible",
                        refcat = NULL, bhaz = NULL, ...) {
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
    linkinv = function(eta) inv_link(eta, link = slink)
  )
  out[names(family_info)] <- family_info
  class(out) <- c("brmsfamily", "family")
  all_valid_dpars <- c(valid_dpars(out), valid_dpars(out, type = "multi"))
  for (dp in all_valid_dpars) {
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
    # TODO: move specification of 'threshold' to the 'resp_thres' function?
    thres_options <- c("flexible", "equidistant", "sum_to_zero")
    out$threshold <- match.arg(threshold, thres_options)
  }
  if (conv_cats_dpars(out$family)) {
    if (!has_joint_link(out$family)) {
      out$refcat <- NA
    } else if (!is.null(refcat)) {
      allow_na_ref <- !is_logistic_normal(out$family)
      out$refcat <- as_one_character(refcat, allow_na = allow_na_ref)
    }
  }
  if (is_cox(out$family)) {
    if (!is.null(bhaz)) {
      if (!is.list(bhaz)) {
        stop2("'bhaz' should be a list.")
      }
      out$bhaz <- bhaz
    } else {
      out$bhaz <- list()
    }
    # set default arguments
    if (is.null(out$bhaz$df)) {
      out$bhaz$df <- 5L
    }
    if (is.null(out$bhaz$intercept)) {
      out$bhaz$intercept <- TRUE
    }
  }
  out
}

# checks and corrects validity of the model family
# @param family Either a function, an object of class 'family'
#   or a character string of length one or two
# @param link an optional character string naming the link function
#   ignored if family is a function or a family object
# @param threshold optional character string specifying the threshold
#   type in ordinal models
validate_family <- function(family, link = NULL, threshold = NULL) {
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

# extract special information of families
# @param x object from which to extract
# @param y name of the component to extract
family_info <- function(x, y, ...) {
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
family_info.list <- function(x, y, ...) {
  ulapply(x, family_info, y = y, ...)
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

# combine information from multiple families
# provides special handling for certain elements
combine_family_info <- function(x, y, ...) {
  y <- as_one_character(y)
  unite <- c(
    "dpars", "type", "specials", "include",
    "const", "cats", "ad", "normalized"
  )
  if (y %in% c("family", "link")) {
    x <- unlist(x)
  } else if (y %in% unite) {
    x <- Reduce("union", x)
  } else if (y == "ybounds") {
    x <- do_call(rbind, x)
    x <- c(max(x[, 1]), min(x[, 2]))
  } else if (y == "closed") {
    # closed only if no bounds are open
    x <- do_call(rbind, x)
    clb <- !any(ulapply(x[, 1], isFALSE))
    cub <- !any(ulapply(x[, 2], isFALSE))
    x <- c(clb, cub)
  } else if (y == "thres") {
    # thresholds are the same across mixture components
    x <- x[[1]]
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
beta_binomial <- function(link = "logit", link_phi = "log") {
  slink <- substitute(link)
  .brmsfamily("beta_binomial", link = link, slink = slink, link_phi = link_phi)
}

#' @rdname brmsfamily
#' @export
negbinomial <- function(link = "log", link_shape = "log") {
  slink <- substitute(link)
  .brmsfamily("negbinomial", link = link, slink = slink,
              link_shape = link_shape)
}

# not yet officially supported
# @rdname brmsfamily
# @export
negbinomial2 <- function(link = "log", link_sigma = "log") {
  slink <- substitute(link)
  .brmsfamily("negbinomial2", link = link, slink = slink,
              link_sigma = link_sigma)
}

#' @rdname brmsfamily
#' @export
geometric <- function(link = "log") {
  slink <- substitute(link)
  .brmsfamily("geometric", link = link, slink = slink)
}

# do not export yet!
# @rdname brmsfamily
# @export
discrete_weibull <- function(link = "logit", link_shape = "log") {
  slink <- substitute(link)
  .brmsfamily("discrete_weibull", link = link, slink = slink,
              link_shape = link_shape)
}

# do not export yet!
# @rdname brmsfamily
# @export
com_poisson <- function(link = "log", link_shape = "log") {
  slink <- substitute(link)
  .brmsfamily("com_poisson", link = link, slink = slink,
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

# not yet exported
# @rdname brmsfamily
# @export
dirichlet2 <- function(link = "log") {
  slink <- substitute(link)
  .brmsfamily("dirichlet2", link = link, slink = slink, refcat = NA)
}

#' @rdname brmsfamily
#' @export
logistic_normal <- function(link = "identity", link_sigma = "log",
                            refcat = NULL) {
  slink <- substitute(link)
  .brmsfamily("logistic_normal", link = link, slink = slink,
              link_sigma = link_sigma, refcat = refcat)
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

# do not export yet!
# @rdname brmsfamily
# @export
zero_inflated_asym_laplace <- function(link = "identity", link_sigma = "log",
                                       link_quantile = "logit",
                                       link_zi = "logit") {
  slink <- substitute(link)
  .brmsfamily("zero_inflated_asym_laplace", link = link, slink = slink,
              link_sigma = link_sigma, link_quantile = link_quantile,
              link_zi = link_zi)
}

#' @rdname brmsfamily
#' @export
cox <- function(link = "log", bhaz = NULL) {
  slink <- substitute(link)
  .brmsfamily("cox", link = link, bhaz = bhaz)
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
zero_inflated_beta_binomial <- function(link = "logit", link_phi = "log",
                                        link_zi = "logit") {
  slink <- substitute(link)
  .brmsfamily("zero_inflated_beta_binomial", link = link, slink = slink,
              link_phi = link_phi, link_zi = link_zi)
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
                       threshold = "flexible") {
  slink <- substitute(link)
  .brmsfamily("cumulative", link = link, slink = slink,
              link_disc = link_disc, threshold = threshold)
}

#' @rdname brmsfamily
#' @export
sratio <- function(link = "logit", link_disc = "log",
                   threshold = "flexible") {
  slink <- substitute(link)
  .brmsfamily("sratio", link = link, slink = slink,
              link_disc = link_disc, threshold = threshold)
}

#' @rdname brmsfamily
#' @export
cratio <- function(link = "logit", link_disc = "log",
                   threshold = "flexible") {
  slink <- substitute(link)
  .brmsfamily("cratio", link = link, slink = slink,
              link_disc = link_disc, threshold = threshold)
}

#' @rdname brmsfamily
#' @export
acat <- function(link = "logit", link_disc = "log",
                 threshold = "flexible") {
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
#'   of the mean parameters are ordered in non-ordinal models
#'   and fixed to the same value in ordinal models (see details).
#'   If \code{'none'} or \code{FALSE}, no ordering constraint is applied.
#'   If \code{NULL} (the default), \code{order} is set to \code{'mu'}
#'   if all families are the same and \code{'none'} otherwise.
#'   Other ordering constraints may be implemented in the future.
#'
#' @return An object of class \code{mixfamily}.
#'
#' @details
#'
#' Most families supported by \pkg{brms} can be used to form mixtures. The
#' response variable has to be valid for all components of the mixture family.
#' Currently, the number of mixture components has to be specified by the user.
#' It is not yet possible to estimate the number of mixture components from the
#' data.
#'
#' Ordering intercepts in mixtures of ordinal families is not possible as each
#' family has itself a set of vector of intercepts (i.e. ordinal thresholds).
#' Instead, \pkg{brms} will fix the vector of intercepts across components in
#' ordinal mixtures, if desired, so that users can try to identify the mixture
#' model via selective inclusion of predictors.
#'
#' For most mixture models, you may want to specify priors on the
#' population-level intercepts via \code{\link{set_prior}} to improve
#' convergence. In addition, it is sometimes necessary to set \code{init = 0}
#' in the call to \code{\link{brm}} to allow chains to initialize properly.
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
#'             init = 0, chains = 2)
#' summary(fit3)
#' pp_check(fit3)
#'
#' ## predict the mixing proportions
#' fit4 <- brm(bf(y ~ x + z, theta2 ~ x),
#'             dat, family = mix, prior = prior,
#'             init = 0, chains = 2)
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
    mix = lapply(dots, validate_family)
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
#'   continuous (\code{"real"}) or discrete (\code{"int"}). This controls
#'   if the corresponding density function will be named with
#'   \code{<name>_lpdf} or \code{<name>_lpmf}.
#' @param lb Vector of lower bounds of the distributional
#'   parameters. Defaults to \code{NA} that is no lower bound.
#' @param ub Vector of upper bounds of the distributional
#'   parameters. Defaults to \code{NA} that is no upper bound.
#' @param vars Names of variables that are part of the likelihood function
#'   without being distributional parameters. That is, \code{vars} can be used
#'   to pass data to the likelihood. Such arguments will be added to the list of
#'   function arguments at the end, after the distributional parameters. See
#'   \code{\link{stanvar}} for details about adding self-defined data to the
#'   generated \pkg{Stan} model. Addition arguments \code{vreal} and \code{vint}
#'   may be used for this purpose as well (see Examples below). See also
#'   \code{\link{brmsformula}} and \code{\link{addition-terms}} for more
#'   details.
#' @param loop Logical; Should the likelihood be evaluated via a loop
#'   (\code{TRUE}; the default) over observations in Stan?
#'   If \code{FALSE}, the Stan code will be written in a vectorized
#'   manner over observations if possible.
#' @param specials A character vector of special options to enable
#'   for this custom family. Currently for internal use only.
#' @param threshold Optional threshold type for custom ordinal families.
#'   Ignored for non-ordinal families.
#' @param log_lik Optional function to compute log-likelihood values of
#'   the model in \R. This is only relevant if one wants to ensure
#'   compatibility with method \code{\link[brms:log_lik.brmsfit]{log_lik}}.
#' @param posterior_predict Optional function to compute posterior prediction of
#'   the model in \R. This is only relevant if one wants to ensure compatibility
#'   with method \code{\link[brms:posterior_predict.brmsfit]{posterior_predict}}.
#' @param posterior_epred Optional function to compute expected values of the
#'   posterior predictive distribution of the model in \R. This is only relevant
#'   if one wants to ensure compatibility with method
#'   \code{\link[brms:posterior_epred.brmsfit]{posterior_epred}}.
#' @param predict Deprecated alias of `posterior_predict`.
#' @param fitted Deprecated alias of `posterior_epred`.
#' @param env An \code{\link{environment}} in which certain post-processing
#'   functions related to the custom family can be found, if there were not
#'   directly passed to \code{custom_family}. This is only
#'   relevant if one wants to ensure compatibility with the methods
#'   \code{\link[brms:log_lik.brmsfit]{log_lik}},
#'   \code{\link[brms:posterior_predict.brmsfit]{posterior_predict}}, or
#'   \code{\link[brms:posterior_epred.brmsfit]{posterior_epred}}.
#'   By default, \code{env} is the environment from which
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
#' @seealso \code{\link{brmsfamily}}, \code{\link{brmsformula}},
#'    \code{\link{stanvar}}
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
#'   type = "int", vars = "vint1[n]"
#' )
#'
#' # define the corresponding Stan density function
#' stan_density <- "
#'   real beta_binomial2_lpmf(int y, real mu, real phi, int N) {
#'     return beta_binomial_lpmf(y | N, mu * phi, (1 - mu) * phi);
#'   }
#' "
#' stanvars <- stanvar(scode = stan_density, block = "functions")
#'
#' # fit the model
#' fit <- brm(y | vint(ntrials) ~ z, data = dat,
#'            family = beta_binomial2, stanvars = stanvars)
#' summary(fit)
#'
#'
#' # define a *vectorized* custom family (no loop over observations)
#' # notice also that 'vint' no longer has an observation index
#' beta_binomial2_vec <- custom_family(
#'   "beta_binomial2", dpars = c("mu", "phi"),
#'   links = c("logit", "log"), lb = c(NA, 0),
#'   type = "int", vars = "vint1", loop = FALSE
#' )
#'
#' # define the corresponding Stan density function
#' stan_density_vec <- "
#'   real beta_binomial2_lpmf(int[] y, vector mu, real phi, int[] N) {
#'     return beta_binomial_lpmf(y | N, mu * phi, (1 - mu) * phi);
#'   }
#' "
#' stanvars_vec <- stanvar(scode = stan_density_vec, block = "functions")
#'
#' # fit the model
#' fit_vec <- brm(y | vint(ntrials) ~ z, data = dat,
#'            family = beta_binomial2_vec,
#'            stanvars = stanvars_vec)
#' summary(fit_vec)
#' }
#'
#' @export
custom_family <- function(name, dpars = "mu", links = "identity",
                          type = c("real", "int"), lb = NA, ub = NA,
                          vars = NULL, loop = TRUE, specials = NULL,
                          threshold = "flexible",
                          log_lik = NULL, posterior_predict = NULL,
                          posterior_epred = NULL, predict = NULL,
                          fitted = NULL, env = parent.frame()) {
  name <- as_one_character(name)
  dpars <- as.character(dpars)
  links <- as.character(links)
  type <- match.arg(type)
  lb <- as.character(lb)
  ub <- as.character(ub)
  vars <- as.character(vars)
  loop <- as_one_logical(loop)
  specials <- as.character(specials)
  env <- as.environment(env)
  posterior_predict <- use_alias(posterior_predict, predict)
  posterior_epred <- use_alias(posterior_epred, fitted)
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
    if (!is_equal(args[1:2], c("i", "prep"))) {
      stop2("The first two arguments of 'log_lik' ",
            "should be 'i' and 'prep'.")
    }
  }
  if (!is.null(posterior_predict)) {
    posterior_predict <- as.function(posterior_predict)
    args <- names(formals(posterior_predict))
    if (!is_equal(args[1:3], c("i", "prep", "..."))) {
      stop2("The first three arguments of 'posterior_predict' ",
            "should be 'i', 'prep', and '...'.")
    }
  }
  if (!is.null(posterior_epred)) {
    posterior_epred <- as.function(posterior_epred)
    args <- names(formals(posterior_epred))
    if (!is_equal(args[1], "prep")) {
      stop2("The first argument of 'posterior_epred' should be 'prep'.")
    }
  }
  lb <- named_list(dpars, lb)
  ub <- named_list(dpars, ub)
  is_mu <- "mu" == dpars
  link <- links[is_mu]
  normalized <- ""
  out <- nlist(
    family = "custom", link, name,
    dpars, lb, ub, type, vars, loop, specials,
    log_lik, posterior_predict, posterior_epred, env,
    normalized
  )
  if (length(dpars) > 1L) {
    out[paste0("link_", dpars[!is_mu])] <- links[!is_mu]
  }
  class(out) <- c("customfamily", "brmsfamily", "family")
  if (is_ordinal(out)) {
    threshold <- match.arg(threshold)
    out$threshold <- threshold
  }
  out
}

# get post-processing methods for custom families
custom_family_method <- function(family, name) {
  if (!is.customfamily(family)) {
    return(NULL)
  }
  out <- family[[name]]
  if (!is.function(out)) {
    out <- paste0(name, "_", family$name)
    out <- get(out, family$env)
  }
  out
}

# get valid distributional parameters for a family
valid_dpars <- function(family, ...) {
  UseMethod("valid_dpars")
}

#' @export
valid_dpars.default <- function(family, type = NULL, ...) {
  if (!length(family)) {
    if (is.null(type)) {
      return("mu")
    } else {
      return(NULL)
    }
  }
  family <- validate_family(family)
  info <- paste0(usc(type, "suffix"), "dpars")
  family_info(family, info, ...)
}

#' @export
valid_dpars.mixfamily <- function(family, type = NULL, ...) {
  out <- lapply(family$mix, valid_dpars, type = type, ...)
  for (i in seq_along(out)) {
    if (length(out[[i]])) {
      out[[i]] <- paste0(out[[i]], i) 
    }
  }
  out <- unlist(out)
  if (is.null(type)) {
    c(out) <- paste0("theta", seq_along(family$mix))
  }
  out
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

# class of a distributional parameter
dpar_class <- function(dpar, family = NULL) {
  out <- sub("[[:digit:]]*$", "", dpar)
  if (!is.null(family)) {
    # TODO: avoid these special cases by changing naming conventions
    # perhaps add a protected "C" before category names
    # and a protected "M" for mixture components
    if (conv_cats_dpars(family)) {
      # categorical-like models have non-integer suffixes
      # that will not be caught by the standard procedure
      multi_dpars <- valid_dpars(family, type = "multi")
      for (dp in multi_dpars) {
        sel <- grepl(paste0("^", dp), out)
        out[sel] <- dp
      }
    }
  }
  out
}

# id of a distributional parameter
dpar_id <- function(dpar) {
  out <- get_matches("[[:digit:]]+$", dpar, simplify = FALSE)
  ulapply(out, function(x) ifelse(length(x), x, ""))
}

# link functions for distributional parameters
links_dpars <- function(dpar) {
  if (!length(dpar)) dpar <- ""
  switch(dpar,
    character(0),
    mu = "identity",  # not actually used
    sigma = c("log", "identity", "softplus", "squareplus"),
    shape = c("log", "identity", "softplus", "squareplus"),
    nu = c("logm1", "identity"),
    phi = c("log", "identity", "softplus", "squareplus"),
    kappa = c("log", "identity", "softplus", "squareplus"),
    beta = c("log", "identity", "softplus", "squareplus"),
    zi = c("logit", "identity"),
    hu = c("logit", "identity"),
    zoi = c("logit", "identity"),
    coi = c("logit", "identity"),
    disc = c("log", "identity", "softplus", "squareplus"),
    bs = c("log", "identity", "softplus", "squareplus"),
    ndt = c("log", "identity", "softplus", "squareplus"),
    bias = c("logit", "identity"),
    quantile = c("logit", "identity"),
    xi = c("log1p", "identity"),
    alpha = c("identity", "log", "softplus", "squareplus"),
    theta = c("identity")
  )
}

# is a distributional parameter a mixture proportion?
is_mix_proportion <- function(dpar, family) {
  dpar_class <- dpar_class(dpar, family)
  dpar_class %in% "theta" & is.mixfamily(family)
}

# generate a family object of a distributional parameter
dpar_family <- function(family, dpar, ...) {
  UseMethod("dpar_family")
}

#' @export
dpar_family.default <- function(family, dpar, ...) {
  dp_class <- dpar_class(dpar, family)
  if (dp_class == "mu") {
    if (conv_cats_dpars(family)) {
      link <- NULL
      if (!has_joint_link(family)) {
        link <- family$link
      }
      # joint links are applied directly in the likelihood function
      # so link is treated as 'identity'
      out <- .dpar_family(dpar, link)
    } else {
      # standard single mu parameters just store the original family
      out <- family
    }
  } else {
    # link_<dp_class> is always defined for non-mu parameters
    link <- family[[paste0("link_", dp_class)]]
    out <- .dpar_family(dpar, link)
  }
  out
}

#' @export
dpar_family.mixfamily <- function(family, dpar, ...) {
  dp_id <- as.numeric(dpar_id(dpar))
  if (!(length(dp_id) == 1L && is.numeric(dp_id))) {
    stop2("Parameter '", dpar, "' is not a valid mixture parameter.")
  }
  out <- dpar_family(family$mix[[dp_id]], dpar, ...)
  out$order <- family$order
  out
}

# set up special family objects for distributional parameters
# @param dpar name of the distributional parameter
# @param link optional link function of the parameter
.dpar_family <- function(dpar = NULL, link = NULL) {
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
  x <- brmsterms(x)
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

# indicate if family uses real responses
use_real <- function(family) {
  "real" %in% family_info(family, "type")
}

# indicate if family uses integer responses
use_int <- function(family) {
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

is_logistic_normal <- function(family) {
  "logistic_normal" %in% family_info(family, "specials")
}

is_simplex <- function(family) {
  "simplex" %in% family_info(family, "specials")
}

is_polytomous <- function(family) {
  is_categorical(family) || is_ordinal(family) ||
    is_multinomial(family) || is_simplex(family)
}

is_cox <- function(family) {
  "cox" %in% family_info(family, "specials")
}

# has joint link function over multiple inputs
has_joint_link <- function(family) {
  "joint_link" %in% family_info(family, "specials")
}

allow_factors <- function(family) {
  specials <- c("binary", "categorical", "ordinal")
  any(specials %in% family_info(family, "specials"))
}

# check if the family has natural residuals
has_natural_residuals <- function(family) {
  "residuals" %in% family_info(family, "specials")
}

# check if the family allows for residual correlations
has_rescor <- function(family) {
  "rescor" %in% family_info(family, "specials")
}

# check if category specific effects are allowed
allow_cs <- function(family) {
  any(c("cs", "ocs") %in% family_info(family, "specials"))
}

# check if category specific effects should be ordered
needs_ordered_cs <- function(family) {
  "ocs" %in% family_info(family, "specials")
}

# choose dpar names based on categories?
conv_cats_dpars <- function(family) {
  is_categorical(family) || is_multinomial(family) || is_simplex(family)
}

# check if mixtures of the given families are allowed
no_mixture <- function(family) {
  is_categorical(family) || is_multinomial(family) || is_simplex(family)
}

# indicate if the response should consist of multiple columns
has_multicol <- function(family) {
  is_multinomial(family) || is_simplex(family)
}

# indicate if the response is modeled on the log-scale
# even if formally the link function is not 'log'
has_logscale <- function(family) {
  "logscale" %in% family_info(family, "specials")
}

# indicate if family makes use of argument trials
has_trials <- function(family) {
  "trials" %in% family_info(family, "ad") &&
    !"custom" %in% family_names(family)
}

# indicate if family has more than two response categories
has_cat <- function(family) {
  is_categorical(family) || is_multinomial(family) || is_simplex(family)
}

# indicate if family has thresholds
has_thres <- function(family) {
  is_ordinal(family)
}

# indicate if family has equidistant thresholds
has_equidistant_thres <- function(family) {
  "equidistant" %in% family_info(family, "threshold")
}

# indicate if family has sum-to-zero thresholds
has_sum_to_zero_thres <- function(family) {
  "sum_to_zero" %in% family_info(family, "threshold")
}

# indicate if family has ordered thresholds
has_ordered_thres <- function(family) {
  "ordered_thres" %in% family_info(family, "specials")
}

# compute threshold - eta in the likelihood
has_thres_minus_eta <- function(family) {
  "thres_minus_eta" %in% family_info(family, "specials")
}

# compute eta - threshold in the likelihood
has_eta_minus_thres <- function(family) {
  "eta_minus_thres" %in% family_info(family, "specials")
}

# get names of response categories
get_cats <- function(family) {
  family_info(family, "cats")
}

# get reference category categorical-like models
get_refcat <- function(family, int = FALSE) {
  refcat <- family_info(family, "refcat")
  if (int) {
    cats <- family_info(family, "cats")
    refcat <- match(refcat, cats)
  }
  refcat
}

# get names of predicted categories categorical-like models
get_predcats <- function(family) {
  refcat <- family_info(family, "refcat")
  cats <- family_info(family, "cats")
  setdiff(cats, refcat)
}

# get names of ordinal thresholds for prior specification
# @param group name of a group for which to extract categories
get_thres <- function(family, group = "") {
  group <- as_one_character(group)
  thres <- family_info(family, "thres")
  subset2(thres, group = group)$thres
}

# get group names of ordinal thresholds
get_thres_groups <- function(family) {
  thres <- family_info(family, "thres")
  unique(thres$group)
}

# has the model group specific thresholds?
has_thres_groups <- function(family) {
  groups <- get_thres_groups(family)
  any(nzchar(groups))
}

has_ndt <- function(family) {
  "ndt" %in% dpar_class(family_info(family, "dpars"))
}

has_sigma <- function(family) {
  "sigma" %in% dpar_class(family_info(family, "dpars"))
}

# check if sigma should be explicitely set to 0
no_sigma <- function(bterms) {
  stopifnot(is.brmsterms(bterms))
  if (is.formula(bterms$adforms$se)) {
    se <- eval_rhs(bterms$adforms$se)
    se_only <- isFALSE(se$flags$sigma)
    if (se_only && use_ac_cov_time(bterms)) {
      stop2("Please set argument 'sigma' of function 'se' ",
            "to TRUE when modeling time-series covariance matrices.")
    }
  } else {
    se_only <- FALSE
  }
  se_only
}

# has the model a non-predicted but estimated sigma parameter?
simple_sigma <- function(bterms) {
  stopifnot(is.brmsterms(bterms))
  has_sigma(bterms) && !no_sigma(bterms) && !pred_sigma(bterms)
}

# has the model a predicted sigma parameter?
pred_sigma <- function(bterms) {
  stopifnot(is.brmsterms(bterms))
  "sigma" %in% dpar_class(names(bterms$dpars))
}

# do not include a 'nu' parameter in a univariate model?
no_nu <- function(bterms) {
  # the multi_student_t family only has a single 'nu' parameter
  isTRUE(bterms$rescor) && "student" %in% family_names(bterms)
}

# does the family-link combination have a built-in Stan function?
has_built_in_fun <- function(family, link = NULL, dpar = NULL, cdf = FALSE) {
  link <- link %||% family$link
  glm_special <- paste0("sbi", usc(dpar), "_", link, str_if(cdf, "_cdf"))
  all(glm_special %in% family_info(family, "specials"))
}

# suffixes of Stan lpdfs or lpmfs for which only a normalized version exists
always_normalized <- function(family) {
  family_info(family, "normalized")
}

# prepare for calling family specific post-processing functions
prepare_family <- function(x) {
  stopifnot(is.brmsformula(x) || is.brmsterms(x))
  family <- x$family
  acef <- tidy_acef(x)
  if (use_ac_cov_time(acef) && has_natural_residuals(x) && !parameterize_ac_effects(x)) {
    family$fun <- paste0(family$family, "_time")
  } else if (has_ac_class(acef, "sar")) {
    acef_sar <- subset2(acef, class = "sar")
    if (has_ac_subset(acef_sar, type = "lag")) {
      family$fun <- paste0(family$family, "_lagsar")
    } else if (has_ac_subset(acef_sar, type = "error")) {
      family$fun <- paste0(family$family, "_errorsar")
    }
  } else if (has_ac_class(acef, "fcor")) {
    family$fun <- paste0(family$family, "_fcor")
  } else {
    family$fun <- family$family
  }
  family
}

# order intercepts to help identifying mixture components?
# does not work in ordinal models as they have vectors of intercepts
order_intercepts <- function(bterms) {
  dpar <- dpar_class(bterms[["dpar"]])
  if (!length(dpar)) dpar <- "mu"
  isTRUE(!is_ordinal(bterms) && dpar %in% bterms$family[["order"]])
}

# fix intercepts to help identifying mixture components?
# currently enabled only in ordinal models
fix_intercepts <- function(bterms) {
  dpar <- dpar_class(bterms[["dpar"]])
  if (!length(dpar)) dpar <- "mu"
  isTRUE(is_ordinal(bterms) && dpar %in% bterms$family[["order"]])
}

# does the mixture have a joint parameter vector 'theta'
has_joint_theta <- function(bterms) {
  stopifnot(is.brmsterms(bterms))
  is.mixfamily(bterms$family) &&
    !"theta" %in% dpar_class(names(c(bterms$dpars, bterms$fdpars)))
}

# extract family boundaries
family_bounds <- function(x, ...) {
  UseMethod("family_bounds")
}

# @return a named list with one element per response variable
#' @export
family_bounds.mvbrmsterms <- function(x, ...) {
  lapply(x$terms, family_bounds, ...)
}

# bounds of likelihood families
# @return a list with elements 'lb' and 'ub'
#' @export
family_bounds.brmsterms <- function(x, ...) {
  family <- x$family$family
  if (is.null(family)) {
    return(list(lb = -Inf, ub = Inf))
  }
  resp <- usc(x$resp)
  # TODO: define in family-lists.R
  pos_families <- c(
    "poisson", "negbinomial", "negbinomial2", "geometric",
    "gamma", "weibull", "exponential", "lognormal",
    "frechet", "inverse.gaussian",
    "hurdle_poisson", "hurdle_negbinomial", "hurdle_gamma",
    "hurdle_lognormal", "zero_inflated_poisson",
    "zero_inflated_negbinomial"
  )
  beta_families <- c("beta", "zero_inflated_beta", "zero_one_inflated_beta")
  ordinal_families <- c("cumulative", "cratio", "sratio", "acat")
  if (family %in% pos_families) {
    out <- list(lb = 0, ub = Inf)
  } else if (family %in% c("bernoulli", beta_families)) {
    out <- list(lb = 0, ub = 1)
  } else if (family %in% c("categorical", ordinal_families)) {
    out <- list(lb = 1, ub = paste0("ncat", resp))
  } else if (family %in% c("binomial", "zero_inflated_binomial",
                           "beta_binomial", "zero_inflated_beta_binomial")) {
    out <- list(lb = 0, ub = paste0("trials", resp))
  } else if (family %in% "von_mises") {
    out <- list(lb = -pi, ub = pi)
  } else if (family %in% c("wiener", "shifted_lognormal")) {
    out <- list(lb = paste("min_Y", resp), ub = Inf)
  } else {
    out <- list(lb = -Inf, ub = Inf)
  }
  out
}
