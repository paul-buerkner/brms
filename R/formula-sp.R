# This file contains functions dealing with the extended
# formula syntax to specify special effects terms

#' Predictors with Measurement Error in \pkg{brms} Models
#'
#' (Soft deprecated) Specify predictors with measurement error. The function
#' does not evaluate its arguments -- it exists purely to help set up a model.
#'
#' @param x The variable measured with error.
#' @param sdx Known measurement error of \code{x}
#'   treated as standard deviation.
#' @param gr Optional grouping factor to specify which
#'   values of \code{x} correspond to the same value of the
#'   latent variable. If \code{NULL} (the default) each
#'   observation will have its own value of the latent variable.
#'
#' @details
#' For detailed documentation see \code{help(brmsformula)}.
#' \code{me} terms are soft deprecated in favor of the more
#' general and consistent \code{\link{mi}} terms.
#' By default, latent noise-free variables are assumed
#' to be correlated. To change that, add \code{set_mecor(FALSE)}
#' to your model formula object (see examples).
#'
#' @seealso
#' \code{\link{brmsformula}}, \code{\link{brmsformula-helpers}}
#'
#' @examples
#' \dontrun{
#' # sample some data
#' N <- 100
#' dat <- data.frame(
#'   y = rnorm(N), x1 = rnorm(N),
#'   x2 = rnorm(N), sdx = abs(rnorm(N, 1))
#'  )
#' # fit a simple error-in-variables model
#' fit1 <- brm(y ~ me(x1, sdx) + me(x2, sdx), data = dat,
#'             save_pars = save_pars(latent = TRUE))
#' summary(fit1)
#'
#' # turn off modeling of correlations
#' bform <- bf(y ~ me(x1, sdx) + me(x2, sdx)) + set_mecor(FALSE)
#' fit2 <- brm(bform, data = dat, save_pars = save_pars(latent = TRUE))
#' summary(fit2)
#' }
#'
#' @export
me <- function(x, sdx, gr = NULL) {
  # use 'term' for consistency with other special terms
  term <- deparse0(substitute(x))
  sdx <- deparse0(substitute(sdx))
  gr <- substitute(gr)
  if (!is.null(gr)) {
    gr <- deparse0(gr)
    stopif_illegal_group(gr)
  } else {
    gr <- ""
  }
  label <- deparse0(match.call())
  out <- nlist(term, sdx, gr, label)
  class(out) <- c("me_term", "sp_term")
  out
}

#' Predictors with Missing Values in \pkg{brms} Models
#'
#' Specify predictor term with missing values in \pkg{brms}. The function does
#' not evaluate its arguments -- it exists purely to help set up a model.
#' For documentation on how to specify missing values in response variables,
#' see \code{\link{resp_mi}}.
#'
#' @param x The variable containing missing values.
#' @param idx An optional variable containing indices of observations in \code{x}
#'   that are to be used in the model. This is mostly relevant in partially
#'   subsetted models (via \code{resp_subset}) but may also have other
#'   applications that I haven't thought of.
#'
#' @details For detailed documentation see \code{help(brmsformula)}.
#'
#' @seealso \code{\link{brmsformula}}, \code{\link{resp_mi}}
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
#'
#' # simulate some measurement noise
#' nhanes$se <- rexp(N, 2)
#'
#' # measurement noise can be handled within 'mi' terms
#' # with or without the presence of missing values
#' bform2 <- bf(bmi | mi() ~ age * mi(chl)) +
#'   bf(chl | mi(se) ~ age) +
#'   set_rescor(FALSE)
#'
#' fit2 <- brm(bform2, data = nhanes)
#'
#' summary(fit2)
#' plot(conditional_effects(fit2, resp = "bmi"), ask = FALSE)
#'
#' # 'mi' terms can also be used when some responses are subsetted
#' # indicates which observations in the chl model shall be kept
#' nhanes$sub <- TRUE
#' nhanes$sub[1:5] <- FALSE
#'
#' # unique IDs of the chl observations
#' nhanes$id <- 1:N
#'
#' # IDs of the chl observations as used in predicting bmi
#' nhanes$idx <- sample(6:N, N, TRUE)
#'
#' bform3 <- bf(bmi | mi() ~ age * mi(chl, idx = idx)) +
#'   bf(chl | mi(idx = id) + subset(sub) ~ age) +
#'   set_rescor(FALSE)
#'
#' fit3 <- brm(bform3, data = nhanes)
#'
#' summary(fit3)
#' plot(conditional_effects(fit3, resp = "bmi"), ask = FALSE)
#' }
#'
#' @export
mi <- function(x, idx = NA) {
  # use 'term' for consistency with other special terms
  term <- deparse0(substitute(x))
  term_vars <- all_vars(term)
  if (!is_equal(term, term_vars)) {
    stop2("'mi' only accepts single untransformed variables.")
  }
  idx <- deparse0(substitute(idx))
  if (idx != "NA") {
    idx_vars <- all_vars(idx)
    if (!is_equal(idx, idx_vars)) {
      stop2("'mi' only accepts single untransformed variables.")
    }
  }
  label <- deparse0(match.call())
  out <- nlist(term, idx, label)
  class(out) <- c("mi_term", "sp_term")
  out
}

#' Monotonic Predictors in \pkg{brms} Models
#'
#' Specify a monotonic predictor term in \pkg{brms}. The function does not
#' evaluate its arguments -- it exists purely to help set up a model.
#'
#' @param x An integer variable or an ordered factor to be modeled as monotonic.
#' @param id Optional character string. All monotonic terms
#'  with the same \code{id} within one formula will be modeled as
#'  having the same simplex (shape) parameter vector. If all monotonic terms
#'  of the same predictor have the same \code{id}, the resulting
#'  predictions will be conditionally monotonic for all values of
#'  interacting covariates (Bürkner & Charpentier, 2020).
#'
#' @details See Bürkner and Charpentier (2020) for the underlying theory. For
#'   detailed documentation of the formula syntax used for monotonic terms,
#'   see \code{help(brmsformula)} as well as \code{vignette("brms_monotonic")}.
#'
#' @seealso \code{\link{brmsformula}}
#'
#' @references
#' Bürkner P. C. & Charpentier E. (2020). Modeling Monotonic Effects of Ordinal
#' Predictors in Regression Models. British Journal of Mathematical and
#' Statistical Psychology. doi:10.1111/bmsp.12195
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
#' fit1 <- brm(ls ~ mo(income), data = dat)
#' summary(fit1)
#' plot(fit1, N = 6)
#' plot(conditional_effects(fit1), points = TRUE)
#'
#' # model interaction with other variables
#' dat$x <- sample(c("a", "b", "c"), 100, TRUE)
#' fit2 <- brm(ls ~ mo(income)*x, data = dat)
#' summary(fit2)
#' plot(conditional_effects(fit2), points = TRUE)
#'
#' # ensure conditional monotonicity
#' fit3 <- brm(ls ~ mo(income, id = "i")*x, data = dat)
#' summary(fit3)
#' plot(conditional_effects(fit3), points = TRUE)
#' }
#'
#' @export
mo <- function(x, id = NA) {
  # use 'term' for consistency with other special terms
  term <- deparse0(substitute(x))
  id <- as_one_character(id, allow_na = TRUE)
  label <- deparse0(match.call())
  out <- nlist(term, id, label)
  class(out) <- c("mo_term", "sp_term")
  out
}

#' Group-level effects as predictors in \pkg{brms} Models
#'
#' Specify a group-level predictor term in \pkg{brms}. That is,
#' use group-level effects defined somewhere in the model as
#' predictors in another part of the model. The function does not
#' evaluate its arguments -- it exists purely to help set up a model.
#'
#' @param gr Name of the grouping factor of the group-level effect
#'   to be used as predictor.
#' @param coef Optional name of the coefficient of the group-level effect.
#'   Defaults to \code{"Intercept"}.
#' @param resp Optional name of the response variable of the group-level effect.
#' @param dpar Optional name of the distributional parameter of the group-level effect.
#' @param nlpar Optional name of the non-linear parameter of the group-level effect.
#'
#' @seealso \code{\link{brmsformula}}
#'
#' @examples
#' \dontrun{
#' # use the group-level intercept of 'AY' for parameter 'ult'
#' # as predictor for the residual standard deviation 'sigma'
#' # multiplying by 1000 reduces the scale of 'ult' to roughly unity
#' bform <- bf(
#'   cum ~ 1000 * ult * (1 - exp(-(dev/theta)^omega)),
#'   ult ~ 1 + (1|AY), omega ~ 1, theta ~ 1,
#'   sigma ~ re(AY, nlpar = "ult"),
#'   nl = TRUE
#' )
#' bprior <- c(
#'   prior(normal(5, 1), nlpar = "ult"),
#'   prior(normal(1, 2), nlpar = "omega"),
#'   prior(normal(45, 10), nlpar = "theta"),
#'   prior(normal(0, 0.5), dpar = "sigma")
#' )
#'
#' fit <- brm(
#'   bform, data = loss,
#'   family = gaussian(),
#'   prior = bprior,
#'   control = list(adapt_delta = 0.9),
#'   chains = 2
#' )
#' summary(fit)
#'
#' # shows how sigma varies as a function of the AY levels
#' conditional_effects(fit, "AY", dpar = "sigma", re_formula = NULL)
#' }
#'
#' @export
re <- function(gr, coef = "Intercept", resp = "", dpar = "", nlpar = "") {
  term <- as_one_character(deparse_no_string(substitute(gr)))
  coef <- as_one_character(coef)
  resp <- as_one_character(resp)
  dpar <- as_one_character(dpar)
  nlpar <- as_one_character(nlpar)
  label <- deparse0(match.call())
  out <- nlist(term, coef, resp, dpar, nlpar, label)
  class(out) <- c("re_term", "sp_term")
  out
}

# find variable names for which to keep NAs
vars_keep_na <- function(x, ...) {
  UseMethod("vars_keep_na")
}

#' @export
vars_keep_na.mvbrmsterms <- function(x, ...) {
  resps <- get_element(x, "respform")
  resps <- ulapply(resps, terms_resp, check_names = FALSE)
  out <- lapply(x$terms, vars_keep_na, responses = resps, ...)
  vars_mi <- unique(ulapply(out, attr, "vars_mi"))
  out <- unique(unlist(out))
  miss_mi <- setdiff(vars_mi, out)
  if (length(miss_mi)) {
    stop2(
      "Response models of variables in 'mi' terms require " ,
      "specification of the addition argument 'mi'. See ?mi. ",
      "Error occurred for ", collapse_comma(miss_mi), "."
    )
  }
  out
}

#' @export
vars_keep_na.brmsterms <- function(x, responses = NULL, ...) {
  out <- character(0)
  if (is.formula(x$adforms$mi)) {
    mi_respcall <- terms_resp(x$respform, check_names = FALSE)
    mi_respvars <- all_vars(mi_respcall)
    mi_advars <- all_vars(x$adforms$mi)
    c(out) <- unique(c(mi_respcall, mi_respvars, mi_advars))
  }
  if (is.formula(x$adforms$cens)) {
    y2_expr <- get_ad_expr(x, "cens", "y2", type = "vars")
    c(out) <- all_vars(y2_expr)
  }
  uni_mi <- ulapply(get_effect(x, "sp"), attr, "uni_mi")
  if (length(uni_mi)) {
    vars_mi <- ulapply(uni_mi, function(term) eval2(term)$term)
    miss_mi <- setdiff(vars_mi, responses)
    if (length(miss_mi)) {
      stop2(
        "Variables in 'mi' terms should also be specified " ,
        "as response variables in the model. See ?mi. ",
        "Error occurred for ", collapse_comma(miss_mi), "."
      )
    }
    attr(out, "vars_mi") <- vars_mi
  }
  out
}

# extract unique names of noise-free terms
get_uni_me <- function(x) {
  uni_me <- ulapply(get_effect(x, "sp"), attr, "uni_me")
  if (!length(uni_me)) {
    return(NULL)
  }
  xname <- ulapply(uni_me, function(term) eval2(term)$term)
  df <- data.frame(xname, uni_me)
  df <- df[!duplicated(df), ]
  xdupl <- df$xname[duplicated(df$xname)]
  if (length(xdupl)) {
    calls <- df$uni_me[df$xname == xdupl[1]]
    stop2(
      "Variable '", xdupl[1], "' is used in different calls to 'me'.\n",
      "Associated calls are: ", collapse_comma(calls)
    )
  }
  unique(uni_me)
}

# save all me-terms within a tidy data.frame
frame_me <- function(x, data, old_levels = NULL) {
  uni_me <- get_uni_me(x)
  if (!length(uni_me)) {
    return(empty_meframe())
  }
  if (has_subset(x)) {
    # 'Xme' variables need to be the same across univariate models
    stop2("Argument 'subset' is not supported when using 'me' terms.")
  }
  out <- data.frame(
    term = uni_me, xname = "", grname = "",
    stringsAsFactors = FALSE
  )
  levels <- list()
  for (i in seq_rows(out)) {
    tmp <- eval2(out$term[i])
    out$xname[i] <- tmp$term
    if (isTRUE(nzchar(tmp$gr))) {
      out$grname[i] <- tmp$gr
      if (is.null(levels[[tmp$gr]])) {
        levels[[tmp$gr]] <- extract_levels(get(tmp$gr, data))
      }
    }
  }
  out$coef <- rename(paste0("me", out$xname))
  out$cor <- isTRUE(x$mecor)
  if (!is.null(old_levels)) {
    # for newdata numeration has to depend on the original levels
    set_levels(out) <- old_levels
    set_levels(out, "used") <- levels
  } else {
    set_levels(out) <- levels
  }
  class(out) <- meframe_class()
  out
}

empty_meframe <- function() {
  out <- data.frame(
    term = character(0), xname = character(0),
    grname = character(0), cor = logical(0),
    stringsAsFactors = FALSE
  )
  class(out) <- meframe_class()
  out
}

meframe_class <- function() {
  c("meframe", "data.frame")
}

is.meframe <- function(x) {
  inherits(x, "meframe")
}

# handle default of correlations between 'me' terms
default_mecor <- function(mecor = NULL) {
  if (is.null(mecor)) TRUE else as_one_logical(mecor)
}

# find names of all variables used in a special effects type
get_sp_vars <- function(x, type, name = NULL) {
  sp_terms <- ulapply(get_effect(x, "sp"), all_terms)
  sp_terms <- unique(get_matches_expr(regex_sp(type), sp_terms))
  if (is.null(name)) {
    # extract all variable names
    out <- all_vars(sp_terms)
  } else {
    # extract only variable names from a specific sp term argument
    out <- all_vars(ulapply(sp_terms, function(x) eval2(x)[[name]]))
  }
  out
}

# gather information of special effects terms
# @param x a formula, brmsterms, or brmsframe object
# @return a data.frame with one row per special term
# TODO: refactor to store in long format to avoid several list columns?
#   or go full out on list columns with one column per predictor type?
frame_sp <- function(x, data) {
  if (is.formula(x)) {
    x <- brmsterms(x, check_response = FALSE)$dpars$mu
  }
  form <- x[["sp"]]
  if (!is.formula(form)) {
    return(empty_data_frame())
  }
  mm <- sp_model_matrix(form, data, rename = FALSE)
  out <- data.frame(term = colnames(mm), stringsAsFactors = FALSE)
  out$coef <- rename(out$term)
  calls_cols <- c(paste0("calls_", all_sp_types()), "joint_call")
  list_cols <- c("vars_mi", "idx_mi", "idx2_mi", "ids_mo", "Imo", "reframe")
  for (col in c(calls_cols, list_cols)) {
    out[[col]] <- vector("list", nrow(out))
  }
  kmo <- 0
  terms_split <- strsplit(out$term, ":")
  for (i in seq_rows(out)) {
    # prepare mo terms
    take_mo <- grepl_expr(regex_sp("mo"), terms_split[[i]])
    if (any(take_mo)) {
      out$calls_mo[[i]] <- terms_split[[i]][take_mo]
      nmo <- length(out$calls_mo[[i]])
      out$Imo[[i]] <- (kmo + 1):(kmo + nmo)
      out$ids_mo[[i]] <- rep(NA, nmo)
      kmo <- kmo + nmo
      for (j in seq_along(out$calls_mo[[i]])) {
        mo_term <- out$calls_mo[[i]][[j]]
        mo_match <- get_matches_expr(regex_sp("mo"), mo_term)
        if (length(mo_match) > 1L || nchar(mo_match) < nchar(mo_term)) {
          stop2("The monotonic term '",  mo_term, "' is invalid.")
        }
        out$ids_mo[[i]][j] <- eval2(mo_term)$id
      }
    }
    # prepare me terms
    take_me <- grepl_expr(regex_sp("me"), terms_split[[i]])
    if (any(take_me)) {
      out$calls_me[[i]] <- terms_split[[i]][take_me]
      # remove 'I' (identity) function calls that
      # were used solely to separate formula terms
      out$calls_me[[i]] <- gsub("^I\\(", "(", out$calls_me[[i]])
    }
    # prepare mi terms
    take_mi <- grepl_expr(regex_sp("mi"), terms_split[[i]])
    if (any(take_mi)) {
      mi_parts <- terms_split[[i]][take_mi]
      out$calls_mi[[i]] <- get_matches_expr(regex_sp("mi"), mi_parts)
      out$vars_mi[[i]] <- out$idx_mi[[i]] <- rep(NA, length(out$calls_mi[[i]]))
      for (j in seq_along(out$calls_mi[[i]])) {
        mi_term <- eval2(out$calls_mi[[i]][[j]])
        out$vars_mi[[i]][j] <- mi_term$term
        if (mi_term$idx != "NA") {
          out$idx_mi[[i]][j] <- mi_term$idx
        }
      }
      # do it like terms_resp to ensure correct matching
      out$vars_mi[[i]] <- gsub("\\.|_", "", make.names(out$vars_mi[[i]]))
    }
    take_re <- grepl_expr(regex_sp("re"), terms_split[[i]])
    if (any(take_re)) {
      re_parts <- terms_split[[i]][take_re]
      out$calls_re[[i]] <- get_matches_expr(regex_sp("re"), re_parts)
      out$reframe[[i]] <- vector("list", length(out$calls_re[[i]]))
      for (j in seq_along(out$calls_re[[i]])) {
        re_call <- out$calls_re[[i]][[j]]
        re_term <- eval2(re_call)
        if (!is.null(x$frame$re)) {
          stopifnot(is.reframe(x$frame$re))
          cols <- c("coef", "resp", "dpar", "nlpar")
          rf <- subset2(x$frame$re, group = re_term$term, ls = re_term[cols])
          # Ideally we should check here if the required re term can be found.
          # However this will lead to errors in post-processing even if the
          # re terms are not actually evaluated. See prepare_predictions_sp
          # for more details. The necessary pre-processing validity check
          # is instead done in stan_sp.
          # if (!NROW(rf)) {
          #   stop2("Cannot find varying coefficients belonging to ", re_call, ".")
          # }
          # there should theoretically never be more than one matching row
          stopifnot(NROW(rf) <= 1L)
          if (isTRUE(rf$gtype == "mm")) {
            stop2("Multimembership terms are not yet supported by 're'.")
          }
          out$reframe[[i]][[j]] <- rf
        }
      }
      if (!isNULL(out$reframe[[i]])) {
        out$reframe[[i]] <- Reduce(rbind, out$reframe[[i]])
      } else {
        out$reframe[[i]] <- empty_reframe()
      }
    }
    has_sp_calls <- grepl_expr(regex_sp(all_sp_types()), terms_split[[i]])
    sp_calls <- sub("^I\\(", "(", terms_split[[i]][has_sp_calls])
    out$joint_call[[i]] <- paste0(sp_calls, collapse = " * ")
  }

  # extract data frame to track all required index variables
  uni_mi <- unique(data.frame(
    var = unlist(out$vars_mi),
    idx = unlist(out$idx_mi),
    stringsAsFactors = FALSE
  ))
  uni_mi$idx2 <- rep(NA, nrow(uni_mi))
  for (i in seq_rows(uni_mi)) {
    uni_mi_sub <- subset2(uni_mi, var = uni_mi$var[i])
    uni_mi$idx2[i] <- match(uni_mi$idx[i], na.omit(uni_mi_sub$idx))
  }
  attr(out, "uni_mi") <- uni_mi
  for (i in seq_rows(out)) {
    for (j in seq_along(out$idx_mi[[i]])) {
      sub <- subset2(
        uni_mi, var = out$vars_mi[[i]][j],
        idx = out$idx_mi[[i]][j]
      )
      out$idx2_mi[[i]][j] <- sub$idx2
    }
  }

  # extract information on covariates
  # only non-zero covariates are relevant to consider
  has_covars <- attr(mm, "covars")
  cumsum_covars <- cumsum(has_covars)
  out$Ic <- ifelse(has_covars, cumsum_covars, 0)
  class(out) <- spframe_class()
  out
}

spframe_class <- function() {
  c("spframe", "data.frame")
}

is.spframe <- function(x) {
  inherits(x, "spframe")
}

# extract names of monotonic simplex parameters
# @param spframe output of frame_sp
# @param use_id use the 'id' argument to construct simo labels?
# @return a character vector of length nrow(spframe)
get_simo_labels <- function(spframe, use_id = FALSE) {
  out <- named_list(spframe$term)
  I <- which(lengths(spframe$Imo) > 0)
  for (i in I) {
    # use the ID as label if specified
    out[[i]] <- ifelse(
      use_id & !is.na(spframe$ids_mo[[i]]), spframe$ids_mo[[i]],
      paste0(spframe$coef[i], seq_along(spframe$Imo[[i]]))
    )
  }
  unlist(out)
}

# standard errors of variables with missing values
get_sdy <- function(x, data = NULL) {
  stopifnot(is.brmsterms(x))
  miform <- x$adforms[["mi"]]
  sdy <- NULL
  if (is.formula(miform)) {
    mi <- eval_rhs(miform)
    if (mi$vars$sdy != "NA") {
      sdy <- eval2(mi$vars$sdy, data)
      if (!is.null(sdy) && !is.numeric(sdy)) {
        stop2("Measurement error should be numeric.")
      }
      if (isTRUE(any(sdy <= 0))) {
        stop2("Measurement error should be positive.")
      }
    }
  }
  sdy
}

# get names of grouping variables from me terms
get_me_group_vars <- function(x) {
  uni_me <- get_uni_me(x)
  out <- lapply(uni_me, eval2)
  out <- ufrom_list(out, "gr")
  out[nzchar(out)]
}

# extract reframe belonging to the unique re() predictor terms
get_uni_re_reframe <- function(x, reframe, unique_id = FALSE) {
  stopifnot(is.reframe(reframe))
  uni_re <- unique(ulapply(get_effect(x, "sp"), attr, "uni_re"))
  if (!length(uni_re)) {
    return(empty_reframe())
  }
  out <- vector("list", length(uni_re))
  for (i in seq_along(uni_re)) {
    re_term <- eval2(uni_re[i])
    cols <- c("coef", "resp", "dpar", "nlpar")
    out[[i]] <- subset2(reframe, group = re_term$term, ls = re_term[cols])
  }
  out <- Reduce(rbind, out)
  if (unique_id) {
    out <- out[!duplicated(out$id), ]
  }
  out
}

# get the design matrix of special effects terms
# @param formula a formula containing special effects terms
# @param data data.frame passed by the user
# @param types types of special terms to consider
# @param ... passed to get_model_matrix
# @details special terms will be evaluated to 1 so that columns
#   containing not only ones are those with covariates
# @return design matrix of special effects terms and their covariates
sp_model_matrix <- function(formula, data, types = all_sp_types(), ...) {
  attributes(data)$terms <- NULL
  terms_split <- strsplit(all_terms(formula), split = ":")
  terms_unique <- unique(unlist(terms_split))
  regex <- regex_sp(types)
  terms_replace <- terms_unique[grepl_expr(regex, terms_unique)]
  dummies <- paste0("dummy", seq_along(terms_replace), "__")
  data[dummies] <- list(1)
  terms_comb <- covars <- rep(NA, length(terms_split))
  # loop over terms and add dummy variables
  for (i in seq_along(terms_split)) {
    replace_i <- grepl_expr(regex, terms_split[[i]])
    terms_i_replace <- terms_split[[i]][replace_i]
    dummies_i <- dummies[match(terms_i_replace, terms_replace)]
    terms_split[[i]][replace_i] <- dummies_i
    terms_comb[i] <- paste0(terms_split[[i]], collapse = ":")
    # non-special covariates are part of the term
    covars[i] <- !all(replace_i)
  }
  new_formula <- str2formula(terms_comb)
  attributes(new_formula) <- attributes(formula)
  out <- get_model_matrix(new_formula, data, ...)
  # fixes issue #1504
  colnames(out) <- rm_wsp(colnames(out))
  # recover original column names
  colnames(out) <- rename(colnames(out), dummies, terms_replace)
  attr(out, "covars") <- covars
  out
}

# extract an me variable
get_me_values <- function(term, data) {
  term <- get_sp_term(term)
  stopifnot(is.me_term(term))
  x <- as.vector(eval2(term$term, data))
  if (!is.numeric(x)) {
    stop2("Noisy variables should be numeric.")
  }
  as.array(x)
}

# extract the measurement error of an me term
get_me_noise <- function(term, data) {
  term <- get_sp_term(term)
  stopifnot(is.me_term(term))
  sdx <- as.vector(eval2(term$sdx, data))
  if (length(sdx) == 0L) {
    stop2("Argument 'sdx' is missing in function 'me'.")
  } else if (length(sdx) == 1L) {
    sdx <- rep(sdx, nrow(data))
  }
  if (!is.numeric(sdx)) {
    stop2("Measurement error should be numeric.")
  }
  if (isTRUE(any(sdx <= 0))) {
    stop2("Measurement error should be positive.")
  }
  as.array(sdx)
}

# extract the grouping variable of an me term
get_me_group <- function(term, data) {
  term <- get_sp_term(term)
  stopifnot(is.me_term(term))
  as.array(eval2(term$gr, data))
}

# extract mo variables
get_mo_values <- function(term, data) {
  term <- get_sp_term(term)
  stopifnot(is.mo_term(term))
  x <- eval2(term$term, data)
  if (is.ordered(x)) {
    # counting starts at zero
    max_value <- length(levels(x)) - 1
    x <- as.numeric(x) - 1
  } else if (all(is_wholenumber(x))) {
    min_value <- attr(x, "min")
    if (is.null(min_value)) {
      min_value <- min(x)
    }
    x <- x - min_value
    max_value <- max(x)
  } else {
    stop2(
      "Monotonic predictors must be integers or ordered ",
      "factors. Error occurred for variable '", term$term, "'."
    )
  }
  x <- as.array(x)
  attr(x, "max") <- max_value
  x
}

# prepare 'sp_term' objects
get_sp_term <- function(term) {
  if (!is.sp_term(term)) {
    term <- eval2(as_one_character(term))
  }
  term
}

# all effects which fall under the 'sp' category of brms
all_sp_types <- function() {
  c("mo", "me", "mi", "re")
}

# classes used to set up special effects terms
is.sp_term <- function(x) {
  inherits(x, "sp_term")
}

is.mo_term <- function(x) {
  inherits(x, "mo_term")
}

is.me_term <- function(x) {
  inherits(x, "me_term")
}

is.mi_term <- function(x) {
  inherits(x, "mi_term")
}

is.re_term <- function(x) {
  inherits(x, "re_term")
}
