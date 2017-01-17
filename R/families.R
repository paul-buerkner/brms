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
#'   \code{exgaussian}, \code{wiener}, \code{inverse.gaussian}, 
#'   \code{exponential}, \code{weibull}, \code{frechet},
#'   \code{Beta}, \code{von_mises},
#'   \code{categorical}, \code{cumulative}, \code{cratio}, \code{sratio},
#'   \code{acat}, \code{hurdle_poisson}, \code{hurdle_negbinomial}, 
#'   \code{hurdle_gamma}, \code{hurdle_lognormal}, 
#'   \code{zero_inflated_binomial}, \code{zero_inflated_beta},
#'   \code{zero_inflated_negbinomial}, 
#'   and \code{zero_inflated_poisson}.
#' @param link A specification for the model link function. 
#'   This can be a name/expression or character string. 
#'   See the 'Details' section for more information on link
#'   functions supported by each family.
#'   
#' @details 
#'   Family \code{gaussian} with \code{identity} link leads to linear regression. 
#'   Family \code{student} with \code{identity} link leads to 
#'   robust linear regression that is less influenced by outliers. 
#'   Families \code{poisson}, \code{negbinomial}, and \code{geometric} 
#'   with \code{log} link lead to regression models for count data. 
#'   Families \code{binomial} and \code{bernoulli} with \code{logit} link leads to 
#'   logistic regression and family \code{categorical} to multi-logistic regression 
#'   when there are more than two possible outcomes.
#'   Families \code{cumulative}, \code{cratio} ('contiuation ratio'), 
#'   \code{sratio} ('stopping ratio'), and \code{acat} ('adjacent category') 
#'   leads to ordinal regression. Families \code{Gamma}, \code{weibull}, 
#'   \code{exponential}, \code{lognormal}, \code{frechet}, and 
#'   \code{inverse.gaussian} can be used (among others) for survival regression.
#'   Family \code{exgaussian} ('exponentially modified Gaussian') is especially
#'   suited to model reaction times and the \code{wiener} family provides
#'   an implementation of the Wiener diffusion model. For this family,
#'   the main formula predicts the drift parameter 'delta' and
#'   all other parameters are modeled as auxiliary parameters 
#'   (see \code{\link[brms:brmsformula]{brmsformula}} for details).
#'   Families \code{hurdle_poisson}, \code{hurdle_negbinomial}, \code{hurdle_gamma}, 
#'   \code{hurdle_lognormal}, \code{zero_inflated_poisson},
#'   \code{zero_inflated_negbinomial}, \code{zero_inflated_binomial}, and
#'   \code{zero_inflated_beta} allow to estimate zero-inflated and hurdle models. 
#'   These models can be very helpful when there are many zeros in the data 
#'   that cannot be explained by the primary distribution of the response. 
#'   Families \code{hurdle_lognormal} and \code{hurdle_gamma} are 
#'   especially useful, as traditional \code{lognormal} or \code{Gamma}
#'   models cannot be reasonably fitted for data containing zeros in the response.
#'   
#'   In the following, we list all possible links for each family.
#'   The families \code{gaussian}, \code{student}, and \code{exgaussian}
#'   accept the links (as names) \code{identity}, \code{log}, and \code{inverse};
#'   families \code{poisson}, \code{negbinomial}, and \code{geometric} the links 
#'   \code{log}, \code{identity}, and \code{sqrt}; 
#'   families \code{binomial}, \code{bernoulli}, \code{Beta},
#'   \code{cumulative}, \code{cratio}, \code{sratio}, and \code{acat} 
#'   the links \code{logit}, \code{probit}, \code{probit_approx}, 
#'   \code{cloglog}, and \code{cauchit}; 
#'   family \code{categorical} the link \code{logit};
#'   families \code{Gamma}, \code{weibull}, \code{exponential}, and 
#'   \code{frechet} the links \code{log}, \code{identity}, and \code{inverse};
#'   family \code{lognormal} the links \code{identity} and \code{inverse};
#'   family \code{inverse.gaussian} the links \code{1/mu^2}, 
#'   \code{inverse}, \code{identity} and \code{log}; 
#'   families \code{hurdle_poisson}, \code{hurdle_negbinomial},
#'   \code{hurdle_gamma}, \code{zero_inflated_poisson}, and
#'   \code{zero_inflated_negbinomial} the link \code{log};
#'   families \code{wiener} and \code{hurdle_lognormal} the link \code{identity}.
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
#'   that your data requires exactly this type of model. \cr
#'   
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
brmsfamily <- function(family, link = NULL) {
  slink <- substitute(link)
  .brmsfamily(family, link = link, slink = slink)
}

.brmsfamily <- function(family, link = NULL, slink = link) {
  # helper function to prepare brmsfamily objects
  # Args:
  #   family: character string naming the model family
  #   link: character string naming the link function
  #   slink: can be used with substitute(link) for 
  #          non-standard evaluation of the link function
  # returns:
  #  An object of class = c(brmsfamily, family) to be used
  #  only insided the brms package
  family <- tolower(as.character(family))
  if (length(family) != 1L) {
    stop("Argument 'family' must be of length 1.", call. = FALSE)
  }
  family <- rename(family, symbols = c("^normal$", "^zi_", "^hu_"),
                   subs = c("gaussian", "zero_inflated_", "hurdle_"),
                   fixed = FALSE)
  ok_families <- c(
    "gaussian", "student", "lognormal", 
    "binomial", "bernoulli", "categorical", 
    "poisson", "negbinomial", "geometric", 
    "gamma", "weibull", "exponential", 
    "exgaussian", "frechet", "inverse.gaussian", 
    "wiener", "beta", "von_mises",
    "cumulative", "cratio", "sratio", "acat",
    "hurdle_poisson", "hurdle_negbinomial", "hurdle_gamma",
    "hurdle_lognormal", "zero_inflated_poisson", 
    "zero_inflated_negbinomial", "zero_inflated_binomial", 
    "zero_inflated_beta")
  if (!family %in% ok_families) {
    stop(family, " is not a supported family. Supported families are: \n",
         paste(ok_families, collapse = ", "), call. = FALSE)
  }
  
  # check validity of link
  if (is_linear(family) || family %in% "exgaussian") {
    ok_links <- c("identity", "log", "inverse")
  } else if (family == "inverse.gaussian") {
    ok_links <- c("1/mu^2", "inverse", "identity", "log")
  } else if (is_count(family)) {
    ok_links <- c("log", "identity", "sqrt")
  } else if (is_ordinal(family) || family %in% "zero_inflated_beta") {
    ok_links <- c("logit", "probit", "probit_approx", "cloglog", "cauchit")
  } else if (is_binary(family) || family %in% "beta") {
    ok_links <- c("logit", "probit", "probit_approx", 
                 "cloglog", "cauchit", "identity")
  } else if (family %in% c("categorical", "zero_inflated_binomial")) {
    ok_links <- c("logit")
  } else if (is_skewed(family)) {
    ok_links <- c("log", "identity", "inverse")
  } else if (is_lognormal(family)) {
    ok_links <- c("identity", "inverse")
  } else if (family %in% c("hurdle_lognormal", "wiener")) {
    ok_links <- c("identity")
  } else if (family %in% c("von_mises")) {
    ok_links <- c("tan_half")
  } else if (is_hurdle(family) || is_zero_inflated(family)) {
    # does not include zi_binomial, zi_beta, or hu_lognormal
    ok_links <- c("log")
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
    stop("Argument 'link' must be of length 1.", call. = FALSE)
  }
  if (is.na(slink)) {
    slink <- ok_links[1]
  } 
  if (!slink %in% ok_links) {
    stop("Link '", slink, "' is not a supported link for family '", 
         family, "'. \nSupported links are: ", 
         paste(ok_links, collapse = ", "), call. = FALSE) 
  }
  structure(list(family = family, link = slink), 
            class = c("brmsfamily", "family"))
}

#' @rdname brmsfamily
#' @export
student <- function(link = "identity") {
  slink <- substitute(link)
  .brmsfamily("student", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
bernoulli <- function(link = "logit") {
  slink <- substitute(link)
  .brmsfamily("bernoulli", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
negbinomial <- function(link = "log") {
  slink <- substitute(link)
  .brmsfamily("negbinomial", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
geometric <- function(link = "log") {
  slink <- substitute(link)
  .brmsfamily("geometric", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
lognormal <- function(link = "identity") {
  slink <- substitute(link)
  .brmsfamily("lognormal", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
exponential <- function(link = "log") {
  slink <- substitute(link)
  .brmsfamily("exponential", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
weibull <- function(link = "log") {
  slink <- substitute(link)
  .brmsfamily("weibull", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
frechet <- function(link = "log") {
  slink <- substitute(link)
  .brmsfamily("frechet", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
exgaussian <- function(link = "identity") {
  slink <- substitute(link)
  .brmsfamily("exgaussian", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
wiener <- function(link = "identity") {
  slink <- substitute(link)
  .brmsfamily("wiener", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
Beta <- function(link = "logit") {
  slink <- substitute(link)
  .brmsfamily("beta", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
von_mises <- function(link = "tan_half") {
  slink <- substitute(link)
  .brmsfamily("von_mises", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
hurdle_poisson <- function(link = "log") {
  slink <- substitute(link)
  .brmsfamily("hurdle_poisson", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
hurdle_negbinomial <- function(link = "log") {
  slink <- substitute(link)
  .brmsfamily("hurdle_negbinomial", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
hurdle_gamma <- function(link = "log") {
  slink <- substitute(link)
  .brmsfamily("hurdle_gamma", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
hurdle_lognormal <- function(link = "identity") {
  slink <- substitute(link)
  .brmsfamily("hurdle_lognormal", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
zero_inflated_beta <- function(link = "logit") {
  slink <- substitute(link)
  .brmsfamily("zero_inflated_beta", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
zero_inflated_poisson <- function(link = "log") {
  slink <- substitute(link)
  .brmsfamily("zero_inflated_poisson", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
zero_inflated_negbinomial <- function(link = "log") {
  slink <- substitute(link)
  .brmsfamily("zero_inflated_negbinomial", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
zero_inflated_binomial <- function(link = "logit") {
  slink <- substitute(link)
  .brmsfamily("zero_inflated_binomial", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
categorical <- function(link = "logit") {
  slink <- substitute(link)
  .brmsfamily("categorical", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
cumulative <- function(link = "logit") {
  slink <- substitute(link)
  .brmsfamily("cumulative", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
sratio <- function(link = "logit") {
  slink <- substitute(link)
  .brmsfamily("sratio", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
cratio <- function(link = "logit") {
  slink <- substitute(link)
  .brmsfamily("cratio", link = link, slink = slink)
}

#' @rdname brmsfamily
#' @export
acat <- function(link = "logit") {
  slink <- substitute(link)
  .brmsfamily("acat", link = link, slink = slink)
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
    if (is(family, "family")) {
      link <- family$link
      family <- family$family
    } 
    if (is.character(family)) {
      if (is.null(link)) {
        link <- family[2]
      }
      family <- .brmsfamily(family[1], link = link)
    } else {
      stop("Argument 'family' is invalid.", call. = FALSE)
    } 
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

is.family <- function(x) {
  inherits(x, "family")
}

is_linear <- function(family) {
  # indicate if family is for a linear model
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("gaussian", "student", "cauchy")
}

is_binary <- function(family) {
  # indicate if family is bernoulli or binomial
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("binomial", "bernoulli")
}

is_ordinal <- function(family) {
  # indicate if family is for an ordinal model
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("cumulative", "cratio", "sratio", "acat") 
}

is_categorical <- function(family) {
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% "categorical" 
}

is_skewed <- function(family) {
  # indicate if family is for model with postive skewed response
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("gamma", "weibull", "exponential", "frechet")
}

is_lognormal <- function(family) {
  # indicate if family is lognormal
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("lognormal")
}

is_exgaussian <- function(family) {
  # indicate if family is exgaussian
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("exgaussian")
}

is_wiener <- function(family) {
  # indicate if family is the wiener diffusion model
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("wiener")
}

is_count <- function(family) {
  # indicate if family is for a count model
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("poisson", "negbinomial", "geometric")
}

is_hurdle <- function(family, zi_beta = TRUE) {
  # indicate if family is for a hurdle model
  if (is(family, "family")) {
    family <- family$family
  }
  # zi_beta is technically a hurdle model
  family %in% c("hurdle_poisson", "hurdle_negbinomial", "hurdle_gamma",
                "hurdle_lognormal", if (zi_beta) "zero_inflated_beta")
}

is_zero_inflated <- function(family, zi_beta = FALSE) {
  # indicate if family is for a zero inflated model
  if (is(family, "family")) {
    family <- family$family
  }
  # zi_beta is technically a hurdle model
  family %in% c("zero_inflated_poisson", "zero_inflated_negbinomial",
                "zero_inflated_binomial", if (zi_beta) "zero_inflated_beta")
}

is_2PL <- function(family) {
  # do not remove to provide an informative error message
  # why the special 2PL implementation is not supported anymore
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
  if (is(family, "family")) {
    family <- family$family
  }
  is_linear(family) || is_skewed(family) || 
    family %in% c("lognormal", "exgaussian", "inverse.gaussian", "beta", 
                  "von_mises", "zero_inflated_beta", "hurdle_gamma", 
                  "hurdle_lognormal", "wiener")
}

use_int <- function(family) {
  # indicate if family uses integer responses
  if (is(family, "family")) {
    family <- family$family
  }
  is_binary(family) || has_cat(family) || 
    is_count(family) || is_zero_inflated(family) || 
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
  is_categorical(family) || is_ordinal(family)
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
  family %in% c("student", "frechet")
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

has_beta <- function(family) {
  # indicate if family needs a kappa parameter
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("exgaussian")
}

has_sigma <- function(family, effects = NULL, 
                      autocor = cor_arma(), incmv = FALSE) {
  # indicate if the model needs a sigma parameter
  # Args:
  #  family: model family
  #  effects: list returned by extract_effects
  #  autocor: object of class cor_arma
  #  incmv: should MV (linear) models be treated as having sigma? 
  if (is(family, "family")) {
    family <- family$family
  }
  is_ln_eg <- family %in% c("lognormal", "hurdle_lognormal", "exgaussian")
  if (is.formula(effects$se)) {
    # call .se without evaluating the x argument 
    cl <- rhs(effects$se)[[2]]
    cl[[1]] <- quote(resp_se_no_data)
    se_only <- isFALSE(attr(eval(cl), "sigma")) 
    if (se_only && use_cov(autocor)) {
      stop2("Please set argument 'sigma' of function 'se' ",  
            "to TRUE when modeling ARMA covariance matrices.")
    }
  } else {
    se_only <- FALSE
  }
  out <- (is_linear(family) || is_ln_eg) && 
           !se_only && !is(autocor, "cov_fixed")
  if (!incmv) {
    is_multi <- is_linear(family) && length(effects$response) > 1L
    out <- out && !is_multi
  }
  out
}

allows_cs <- function(family) {
  # checks if category specific effects are allowed
  if (is(family, "family")) {
    family <- family$family
  }
  family %in% c("sratio", "cratio", "acat")
}

is_old_lognormal <- function(family, link = "identity", nresp = 1,
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

is_old_categorical <- function(x) {
  # indicate if the model is and old categorical model
  stopifnot(is(x, "brmsfit"))
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
  ee <- extract_effects(formula(x), family = family(x))
  (is.null(x$version) || x$version <= "0.10.0.9000") &&
    (is_mv(family(x), ee$response) || is_forked(family(x)))
}
