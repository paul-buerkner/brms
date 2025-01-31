#' Additional Response Information
#'
#' Provide additional information on the response variable
#' in \pkg{brms} models, such as censoring, truncation, or
#' known measurement error. Detailed documentation on the use
#' of each of these functions can be found in the Details section
#' of \code{\link{brmsformula}} (under "Additional response information").
#'
#' @name addition-terms
#' @aliases se weights trials thres cat dec cens trunc
#' @aliases index rate subset vreal vint
#'
#' @param x A vector; Ideally a single variable defined in the data (see
#'   Details). Allowed values depend on the function: \code{resp_se} and
#'   \code{resp_weights} require positive numeric values. \code{resp_trials},
#'   \code{resp_thres}, and \code{resp_cat} require positive integers.
#'   \code{resp_dec} requires \code{0} and \code{1}, or alternatively
#'   \code{'lower'} and \code{'upper'}. \code{resp_cens}
#'   requires \code{'left'}, \code{'none'}, \code{'right'}, and
#'   \code{'interval'} (or equivalently \code{-1}, \code{0}, \code{1}, and
#'   \code{2}) to indicate left, no, right, or interval censoring.
#' @param sigma Logical; Indicates whether the residual standard deviation
#'  parameter \code{sigma} should be included in addition to the known
#'  measurement error. Defaults to \code{FALSE} for backwards compatibility,
#'  but setting it to \code{TRUE} is usually the better choice.
#' @param scale Logical; Indicates whether weights should be scaled
#'  so that the average weight equals one. Defaults to \code{FALSE}.
#' @param y2 A vector specifying the upper bounds in interval censoring.
#'  Will be ignored for non-interval censored observations. However, it
#'  should NOT be \code{NA} even for non-interval censored observations to
#'  avoid accidental exclusion of these observations.
#' @param lb A numeric vector or single numeric value specifying
#'   the lower truncation bound.
#' @param ub A numeric vector or single numeric value specifying
#'   the upper truncation bound.
#' @param denom A vector of positive numeric values specifying
#'   the denominator values from which the response rates are computed.
#' @param gr A vector of grouping indicators.
#' @param df Degrees of freedom of baseline hazard splines for Cox models.
#' @param ... For \code{resp_vreal}, vectors of real values.
#'   For \code{resp_vint}, vectors of integer values. In Stan,
#'   these variables will be named \code{vreal1}, \code{vreal2}, ...,
#'   and \code{vint1}, \code{vint2}, ..., respectively.
#'
#' @return A list of additional response information to be processed further
#'   by \pkg{brms}.
#'
#' @details
#'   These functions are almost solely useful when
#'   called in formulas passed to the \pkg{brms} package.
#'   Within formulas, the \code{resp_} prefix may be omitted.
#'   More information is given in the 'Details' section
#'   of \code{\link{brmsformula}} (under "Additional response information").
#'
#'   It is highly recommended to use a single data variable as input
#'   for \code{x} (instead of a more complicated expression) to make sure all
#'   post-processing functions work as expected.
#'
#' @seealso
#'   \code{\link{brm}},
#'   \code{\link{brmsformula}}
#'
#' @examples
#' \dontrun{
#' ## Random effects meta-analysis
#' nstudies <- 20
#' true_effects <- rnorm(nstudies, 0.5, 0.2)
#' sei <- runif(nstudies, 0.05, 0.3)
#' outcomes <- rnorm(nstudies, true_effects, sei)
#' data1 <- data.frame(outcomes, sei)
#' fit1 <- brm(outcomes | se(sei, sigma = TRUE) ~ 1,
#'             data = data1)
#' summary(fit1)
#'
#' ## Probit regression using the binomial family
#' n <- sample(1:10, 100, TRUE)  # number of trials
#' success <- rbinom(100, size = n, prob = 0.4)
#' x <- rnorm(100)
#' data2 <- data.frame(n, success, x)
#' fit2 <- brm(success | trials(n) ~ x, data = data2,
#'             family = binomial("probit"))
#' summary(fit2)
#'
#' ## Survival regression modeling the time between the first
#' ## and second recurrence of an infection in kidney patients.
#' fit3 <- brm(time | cens(censored) ~ age * sex + disease + (1|patient),
#'             data = kidney, family = lognormal())
#' summary(fit3)
#'
#' ## Poisson model with truncated counts
#' fit4 <- brm(count | trunc(ub = 104) ~ zBase * Trt,
#'             data = epilepsy, family = poisson())
#' summary(fit4)
#' }
#'
NULL

# TODO: split into separate docs for each function

#' @rdname addition-terms
#' @export
resp_se <- function(x, sigma = FALSE) {
  se <- deparse0(substitute(x))
  sigma <- as_one_logical(sigma)
  class_resp_special(
    "se", call = match.call(),
    vars = nlist(se), flags = nlist(sigma)
  )
}

#' @rdname addition-terms
#' @export
resp_weights <- function(x, scale = FALSE) {
  weights <- deparse0(substitute(x))
  scale <- as_one_logical(scale)
  class_resp_special(
    "weights", call = match.call(),
    vars = nlist(weights), flags = nlist(scale)
  )
}

#' @rdname addition-terms
#' @export
resp_trials <- function(x) {
  trials <- deparse0(substitute(x))
  class_resp_special("trials", call = match.call(), vars = nlist(trials))
}

#' @rdname addition-terms
#' @export
resp_thres <- function(x, gr = NA) {
  thres <- deparse0(substitute(x))
  gr <- deparse0(substitute(gr))
  class_resp_special("thres", call = match.call(), vars = nlist(thres, gr))
}

#' @rdname addition-terms
#' @export
resp_cat <- function(x) {
  # deprecated as of brms 2.10.5
  # number of thresholds = number of response categories - 1
  thres <- deparse0(substitute(x))
  str_add(thres) <- " - 1"
  class_resp_special(
    "thres", call = match.call(),
    vars = nlist(thres, gr = "NA")
  )
}

#' @rdname addition-terms
#' @export
resp_dec <- function(x) {
  dec <- deparse0(substitute(x))
  class_resp_special("dec", call = match.call(), vars = nlist(dec))
}

#' @rdname addition-terms
#' @export
resp_bhaz <- function(gr = NA, df = 5, ...) {
  gr <- deparse0(substitute(gr))
  df <- as_one_integer(df)
  args <- nlist(df, ...)
  # non-power users shouldn't know they can change 'intercept'
  args$intercept <- args$intercept %||% TRUE
  class_resp_special("bhaz", call = match.call(), vars = nlist(gr), flags = args)
}

#' @rdname addition-terms
#' @export
resp_cens <- function(x, y2 = NA) {
  cens <- deparse0(substitute(x))
  y2 <- deparse0(substitute(y2))
  class_resp_special("cens", call = match.call(), vars = nlist(cens, y2))
}

#' @rdname addition-terms
#' @export
resp_trunc <- function(lb = -Inf, ub = Inf) {
  lb <- deparse0(substitute(lb))
  ub <- deparse0(substitute(ub))
  class_resp_special("trunc", call = match.call(), vars = nlist(lb, ub))
}

#' Additional response information for handling missing values
#'
#' Provide additional information on missing values in response variables
#' in \pkg{brms} models. Supported by all continuous families.
#' Required when using \code{\link{mi}} predictor terms.
#'
#' @param sdy Optional known measurement error of the response
#'   treated as standard deviation. If specified, handles
#'   measurement error and (completely) missing values
#'   at the same time using the plausible-values-technique.
#' @param idx Optional index variable assigning observations with missing responses
#'   to individual latent response values. This allows to express that multiple
#'   observations share the same latent missing value. It is also needed when
#'   using the \code{idx} argument in \code{\link{mi}} predictor terms.
#'
#' @return A list of additional response information to be processed further
#'   by \pkg{brms}.
#'
#' @details
#'   This function is almost solely useful when called in formulas passed to the
#'   \pkg{brms} package. Within formulas, the \code{resp_} prefix may be
#'   omitted. See the Details section of \code{\link{brmsformula}} (under
#'   "Additional response information") for more information.
#'   Please also see \code{\link{mi}} for more details and examples.
#'
#'   It is highly recommended to use single data variables as inputs
#'   for \code{sdy} and \code{idx} (instead of a more complicated expression)
#'   to make sure all post-processing functions work as expected.
#'
#' @seealso \code{\link{brmsformula}}, \code{\link{mi}}
#'
#' @examples
#' \dontrun{
#' data("nhanes", package = "mice")
#' N <- nrow(nhanes)
#'
#' # simple model with missing data
#' bform1 <- bf(bmi | mi() ~ age * mi(chl)) +
#'   bf(chl | mi() ~ age) +
#'   set_rescor(FALSE)
#'
#' fit1 <- brm(bform1, data = nhanes)
#'
#' summary(fit1)
#' plot(conditional_effects(fit1, resp = "bmi"), ask = FALSE)
#' loo(fit1, newdata = na.omit(fit1$data))
#' }
#'
#' @export
resp_mi <- function(sdy = NA, idx = NA) {
  sdy <- deparse0(substitute(sdy))
  idx <- deparse0(substitute(idx))
  class_resp_special("mi", call = match.call(), vars = nlist(sdy, idx))
}

#' (Deprecated) Additional response information for handling response indexes
#'
#' Provide additional information for handling response indexes
#' in \pkg{brms} models. This function is deprecated and replaced
#' by the \code{idx} argument of \code{\link{resp_mi}}.
#'
#' @param x A vector of index values that assigns each observation to an
#' index, which can then be used by other model terms.
#'
#' @return A list of additional response information to be processed further
#'   by \pkg{brms}.
#'
#' @details
#'   This function is almost solely useful when called in formulas passed to the
#'   \pkg{brms} package. Within formulas, the \code{resp_} prefix may be
#'   omitted. See the Details section of \code{\link{brmsformula}} (under
#'   "Additional response information") for more information.
#'
#' @seealso \code{\link{resp_mi}}, \code{\link{brmsformula}}
#'
#' @export
resp_index <- function(x) {
  index <- deparse0(substitute(x))
  class_resp_special("index", call = match.call(), vars = nlist(index))
}

#' @rdname addition-terms
#' @export
resp_rate <- function(denom) {
  denom <- deparse0(substitute(denom))
  class_resp_special("rate", call = match.call(), vars = nlist(denom))
}

#' Additional response information for subsetting data
#'
#' Provide additional information for subsetting data
#' in \pkg{brms} models. This is relevant only in multivariate models
#' where \code{resp_subset} allows to use different subsets of the data
#' in different univariate models.
#'
#' @param x A logical vector or one that can be transformed to it. Observations
#' for which \code{x} is \code{FALSE} will be excluded from the data used in
#' that specific univariate model.
#'
#' @return A list of additional response information to be processed further
#'   by \pkg{brms}.
#'
#' @details
#'   This function is almost solely useful when called in formulas passed to the
#'   \pkg{brms} package. Within formulas, the \code{resp_} prefix may be
#'   omitted. See the Details section of \code{\link{brmsformula}} (under
#'   "Additional response information") for more information.
#'
#' @examples
#' \dontrun{
#' # generate a simple dataset
#' dat <- data.frame(y = rnorm(10), z = rnorm(10),
#'                   sub = c(rep(TRUE, 5), rep(FALSE, 5)))
#'
#' # only use the first 5 observations to predict z
#' bform <- bf(y ~ 1) + bf(z | subset(sub) ~ y) + set_rescor(FALSE)
#' fit <- brm(bform, dat)
#' summary(fit)
#' }
#'
#' @export
resp_subset <- function(x) {
  subset <- deparse0(substitute(x))
  class_resp_special("subset", call = match.call(), vars = nlist(subset))
}

#' @rdname addition-terms
#' @export
resp_vreal <- function(...) {
  vars <- as.list(substitute(list(...)))[-1]
  class_resp_special("vreal", call = match.call(), vars = vars)
}

#' @rdname addition-terms
#' @export
resp_vint <- function(...) {
  vars <- as.list(substitute(list(...)))[-1]
  class_resp_special("vint", call = match.call(), vars = vars)
}

# class underlying response addition terms
# @param type type of the addition term
# @param call the call to the original addition term function
# @param vars named list of unevaluated variables
# @param flags named list of (evaluated) logical indicators
class_resp_special <- function(type, call, vars = list(), flags = list()) {
  type <- as_one_character(type)
  stopifnot(is.call(call), is.list(vars), is.list(flags))
  label <- deparse0(call)
  out <- nlist(type, call, label, vars, flags)
  class(out) <- c("resp_special")
  out
}

# computes data for addition arguments
eval_rhs <- function(formula, data = NULL) {
  formula <- as.formula(formula)
  eval(rhs(formula)[[2]], data, environment(formula))
}

# get expression for a variable of an addition term
# @param x list with potential $adforms elements
# @param ad name of the addition term
# @param target name of the element to extract
# @type type of the element to extract
# @return a character string or NULL
get_ad_expr <- function(x, ad, name, type = "vars") {
  ad <- as_one_character(ad)
  name <- as_one_character(name)
  type <- as_one_character(type)
  if (is.null(x$adforms[[ad]])) {
    return(NULL)
  }
  out <- eval_rhs(x$adforms[[ad]])[[type]][[name]]
  if (type == "vars" && is_equal(out, "NA")) {
    out <- NULL
  }
  out
}

# check if an addition term expression is present
has_ad_expr <- function(x, ad, name, type = "vars") {
  !is.null(get_ad_expr(x, ad, name, type))
}

# get values of a variable used in an addition term
# @return a vector of values or NULL
get_ad_values <- function(x, ad, name, data) {
  expr <- get_ad_expr(x, ad, name, type = "vars")
  eval2(expr, data)
}

# get a flag used in an addition term
# @return TRUE or FALSE
get_ad_flag <- function(x, ad, name) {
  expr <- get_ad_expr(x, ad, name, type = "flags")
  as_one_logical(eval2(expr))
}

# get variable names used in addition terms
get_ad_vars <- function(x, ...) {
  UseMethod("get_ad_vars")
}

#' @export
get_ad_vars.brmsterms <- function(x, ad, name = NULL, ...) {
  ad <- as_one_character(ad)
  if (is.null(name)) {
    out <- all_vars(x$adforms[[ad]])
  } else {
    out <- all_vars(get_ad_expr(x, ad = ad, name = name, type = "vars"))
  }
  out
}

#' @export
get_ad_vars.mvbrmsterms <- function(x, ad, ...) {
  unique(ulapply(x$terms, get_ad_vars, ad = ad, ...))
}

# coerce censored values into the right format
# @param x vector of censoring indicators
# @return transformed vector of censoring indicators
prepare_cens <- function(x) {
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
    return(x)
  }
  x <- unname(x)
  if (is.factor(x)) {
    x <- as.character(x)
  }
  ulapply(x, .prepare_cens)
}

# extract information on censoring of the response variable
# @return vector of censoring indicators or NULL in case of no censoring
get_cens <- function(bterms, data, resp = NULL) {
  if (!is.null(resp)) {
    bterms <- bterms$terms[[resp]]
  }
  out <- NULL
  if (is.formula(bterms$adforms$cens)) {
    out <- get_ad_values(bterms, "cens", "cens", data)
    out <- prepare_cens(out)
  }
  out
}

# indicates if the model may have interval censored observations
has_interval_cens <- function(bterms) {
  !is.null(get_ad_expr(bterms, "cens", "y2"))
}

# extract truncation boundaries
# @param bterms a brmsterms object
# @param data data.frame containing the truncation variables
# @param incl_family include the family in the derivation of the bounds?
# @param stan return bounds in form of Stan syntax?
# @return a list with elements 'lb' and 'ub' or corresponding Stan code
trunc_bounds <- function(bterms, data = NULL, incl_family = FALSE,
                         stan = FALSE, ...) {
  stopifnot(is.brmsterms(bterms))
  if (is.formula(bterms$adforms$trunc)) {
    trunc <- eval_rhs(bterms$adforms$trunc)
  } else {
    trunc <- resp_trunc()
  }
  out <- list(
    lb = eval2(trunc$vars$lb, data),
    ub = eval2(trunc$vars$ub, data)
  )
  if (incl_family) {
    family_bounds <- family_bounds(bterms)
    out$lb <- max(out$lb, family_bounds$lb)
    out$ub <- min(out$ub, family_bounds$ub)
  }
  if (stan) {
    if (any(out$lb > -Inf | out$ub < Inf)) {
      tmp <- c(
        if (out$lb > -Inf) paste0("lower=", out$lb),
        if (out$ub < Inf) paste0("upper=", out$ub)
      )
      out <- paste0("<", paste0(tmp, collapse = ","), ">")
    } else {
      out <- ""
    }
  }
  out
}

# check if addition argument 'subset' is used in the model
# works for both univariate and multivariate models
has_subset <- function(bterms) {
  if (is.brmsterms(bterms)) {
    out <- has_ad_terms(bterms, "subset")
  } else if (is.mvbrmsterms(bterms)) {
    out <- any(ulapply(bterms$terms, has_ad_terms, "subset"))
  } else {
    out <- FALSE
  }
  out
}

# check if a model has certain addition terms
has_ad_terms <- function(bterms, terms) {
  stopifnot(is.brmsterms(bterms), is.character(terms))
  any(ulapply(bterms$adforms[terms], is.formula))
}

# construct a list of indices for cross-formula referencing
# internal version of frame_index
frame_index <- function(x, ...) {
  UseMethod("frame_index")
}

# @param mv is this univariate model part of a multivariate model?
#' @export
frame_index.brmsterms <- function(x, data, mv = FALSE, ...) {
  out <- get_mi_index(x, data)
  if (!is.null(out)) {
    # index variable specified
    if (has_subset(x)) {
      len_old <- length(out)
      subset <- as.logical(get_ad_values(x, "subset", "subset", data))
      # same NA behavior as in subset_data
      subset[is.na(subset)] <- TRUE
      out <- out[subset]
      # if all observations are kept, it isn't really subsetting
      attr(out, "subset") <- length(out) < len_old
    }
  }
  if (!mv) {
    # ensure consistent format for both uni- and multivariate models
    out <- list(out)
    names(out)[1] <- terms_resp(x$respform)
  }
  out
}

#' @export
frame_index.mvbrmsterms <- function(x, data, ...) {
  lapply(x$terms, frame_index, data = data, mv = TRUE, ...)
}

# extract index values of response variables
# levels: shall the output be the (unique) levels of the latent values?
get_mi_index <- function(x, data, levels = FALSE) {
  stopifnot(is.brmsterms(x))
  levels <- as_one_logical(levels)
  if (!is.formula(x$adforms$mi)) {
    return(NULL)
  }
  idx <- get_ad_values(x, "mi", "idx", data)
  if (is.null(idx)) {
    idx <- get_ad_values(x, "index", "index", data)
    if (!is.null(idx)) {
      warning2(
        "Addition term 'index' is deprecated. ",
        "Use argument 'idx' of addition term 'mi' instead."
      )
    }
  }
  if (!levels) {
    return(idx)
  }
  y <- get_model_response(x, data)
  if (has_ad_expr(x, "mi", "sdy")) {
    # measurement error specified such that all observations
    # have corresponding latent values
    if (is.null(idx)) {
      idx <- seq_along(y)
    } # else idx is already in the right format
  } else {
    # idx does only apply to observations with missing responses
    # as such it has reduced length
    is_na_y <- is.na(y)
    if (is.null(idx)) {
      idx <- which(is_na_y)
    } else {
      idx <- idx[is_na_y]
    }
  }
  levels(factor(idx))
}

# check if cross-formula referencing is possible in subsetted models
check_cross_formula_indexing <- function(bterms) {
  sp_terms <- ulapply(get_effect(bterms, "sp"), all_terms)
  me_terms <- get_matches_expr(regex_sp("me"), sp_terms)
  if (length(me_terms)) {
    stop2("Cannot use me() terms in subsetted formulas.")
  }
  mi_terms <- get_matches_expr(regex_sp("mi"), sp_terms)
  idx_vars <- ulapply(mi_terms, function(x) eval2(x)$idx)
  if (any(idx_vars == "NA")) {
    stop2("'mi' predictor terms in subsetted formulas require ",
          "the 'idx' argument to be specified.")
  }
  invisible(TRUE)
}

# does an expression consist of a single variable?
is_single_variable <- function(x) {
  x <- as_one_character(x)
  is_equal(x, all_vars(x))
}
