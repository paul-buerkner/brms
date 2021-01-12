# This file contains functions dealing with the extended 
# formula syntax to specify special effects terms

#' Predictors with Measurement Error in \pkg{brms} Models
#' 
#' Specify predictors with measurement error. The function does not evaluate its
#' arguments -- it exists purely to help set up a model.
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
#' 
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
#'            save_mevars = TRUE)
#' summary(fit1)
#' 
#' # turn off modeling of correlations
#' bform <- bf(y ~ me(x1, sdx) + me(x2, sdx)) + set_mecor(FALSE)
#' fit2 <- brm(bform, data = dat, save_mevars = TRUE)
#' summary(fit2)
#' } 
#' 
#' @export
me <- function(x, sdx, gr = NULL) {
  # use 'term' for consistency with other special terms
  term <- deparse(substitute(x))
  sdx <- deparse(substitute(sdx))
  gr <- substitute(gr)
  if (!is.null(gr)) {
    gr <- deparse_combine(gr)
    stopif_illegal_group(gr)
  } else {
    gr <- ""
  }
  label <- deparse(match.call())
  out <- nlist(term, sdx, gr, label)
  class(out) <- c("me_term", "sp_term")
  out
}

#' Predictors with Missing Values in \pkg{brms} Models
#' 
#' Specify predictor term with missing values in \pkg{brms}. The function does
#' not evaluate its arguments -- it exists purely to help set up a model.
#' 
#' @param x The variable containing missings.
#' 
#' @details For detailed documentation see \code{help(brmsformula)}. 
#' 
#' @seealso \code{\link{brmsformula}}
#'   
#' @examples 
#' \dontrun{
#' data("nhanes", package = "mice")
#' bform <- bf(bmi | mi() ~ age * mi(chl)) +
#'   bf(chl | mi() ~ age) + set_rescor(FALSE)
#' fit <- brm(bform, data = nhanes)
#' summary(fit)
#' plot(conditional_effects(fit, resp = "bmi"), ask = FALSE)
#' LOO(fit, newdata = na.omit(fit$data))
#' } 
#' 
#' @export
mi <- function(x) {
  # use 'term' for consistency with other special terms
  term <- substitute(x)
  vars <- all.vars(term)
  term <- deparse(term)
  if (!is_equal(term, vars)) {
    stop2("'mi' only accepts single untransformed variables.")
  }
  label <- deparse(match.call())
  out <- nlist(term, label)
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
#'  with the same \code{id} within one formula  will be modeled as
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
  term <- deparse(substitute(x))
  id <- as_one_character(id, allow_na = TRUE)
  label <- deparse(match.call())
  out <- nlist(term, id, label)
  class(out) <- c("mo_term", "sp_term")
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
  if (is.formula(x$adforms$mi)) {
    mi_respcall <- terms_resp(x$respform, check_names = FALSE)
    mi_respvars <- all_vars(mi_respcall)
    mi_advars <- all_vars(x$adforms$mi)
    out <- unique(c(mi_respcall, mi_respvars, mi_advars))
  } else {
    out <- character(0)
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
tidy_meef <- function(bterms, data, old_levels = NULL) {
  uni_me <- get_uni_me(bterms)
  if (!length(uni_me)) {
    return(empty_meef()) 
  }
  if (has_subset(bterms)) {
    # 'Xme' variables need to be the same across univariate models
    stop2("Argument 'subset' is not supported when using 'me' terms.")
  }
  out <- data.frame(
    term = uni_me, xname = "", grname = "", 
    stringsAsFactors = FALSE
  )
  levels <- vector("list", nrow(out))
  for (i in seq_rows(out)) {
    tmp <- eval2(out$term[i])
    out$xname[i] <- tmp$term
    if (isTRUE(nzchar(tmp$gr))) {
      out$grname[i] <- tmp$gr
      if (length(old_levels)) {
        levels[[i]] <- old_levels[[tmp$gr]]
      } else {
        levels[[i]] <- levels(factor(get(tmp$gr, data)))
      } 
    }
  }
  out$coef <- rename(paste0("me", out$xname))
  out$cor <- isTRUE(bterms$mecor)
  names(levels) <- out$grname
  levels <- levels[lengths(levels) > 0L]
  if (length(levels)) {
    levels <- levels[!duplicated(names(levels))]
    attr(out, "levels") <- levels
  }
  structure(out, class = c("meef_frame", "data.frame"))
}

empty_meef <- function() {
  out <- data.frame(
    term = character(0), xname = character(0),
    grname = character(0), cor = logical(0),
    stringsAsFactors = FALSE
  )
  structure(out, class = c("meef_frame", "data.frame"))
}

is.meef_frame <- function(x) {
  inherits(x, "meef_frame")
}

# handle default of correlations between 'me' terms
default_mecor <- function(mecor = NULL) {
  if (is.null(mecor)) TRUE else as_one_logical(mecor)
}

# find names of all variables used in a special effects type
get_sp_vars <- function(x, type) {
  sp_terms <- ulapply(get_effect(x, "sp"), all_terms)
  all_vars(str2formula(get_matches_expr(regex_sp(type), sp_terms)))
}

# gather information of special effects terms
# @param x either a formula or a list containing an element "sp"
# @param data data frame containing the monotonic variables
# @return a data.frame with one row per special term
tidy_spef <- function(x, data) {
  if (is.formula(x)) {
    x <- brmsterms(x, check_response = FALSE)$dpars$mu
  }
  form <- x[["sp"]]
  if (!is.formula(form)) {
    return(empty_data_frame())
  }
  mm <- sp_model_matrix(form, data, rename = FALSE)
  out <- data.frame(term = rm_wsp(colnames(mm)), stringsAsFactors = FALSE)
  out$coef <- rename(out$term)
  calls_cols <- paste0("calls_", all_sp_types())
  for (col in c(calls_cols, "joint_call", "vars_mi", "ids_mo", "Imo")) {
    out[[col]] <- vector("list", nrow(out))
  }
  kmo <- 0
  terms_split <- strsplit(out$term, ":")
  for (i in seq_rows(out)) {
    # prepare mo terms
    take_mo <- grepl_expr(regex_sp("mo"), terms_split[[i]])
    if (sum(take_mo)) {
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
        out$ids_mo[[i]][[j]] <- eval2(mo_term)[["id"]]
      }
    }
    # prepare me terms
    take_me <- grepl_expr(regex_sp("me"), terms_split[[i]])
    if (sum(take_me)) {
      out$calls_me[[i]] <- terms_split[[i]][take_me]
      # remove 'I' (identity) function calls that 
      # were used solely to separate formula terms
      out$calls_me[[i]] <- gsub("^I\\(", "(", out$calls_me[[i]])
    }
    # prepare mi terms 
    take_mi <- grepl_expr(regex_sp("mi"), terms_split[[i]])
    if (sum(take_mi)) {
      mi_parts <- terms_split[[i]][take_mi]
      out$calls_mi[[i]] <- get_matches_expr(regex_sp("mi"), mi_parts)
      out$vars_mi[[i]] <- all_vars(str2formula(out$calls_mi[[i]]))
      # do it like terms_resp to ensure correct matching
      out$vars_mi[[i]] <- gsub("\\.|_", "", make.names(out$vars_mi[[i]]))
    }
    has_sp_calls <- grepl_expr(regex_sp(all_sp_types()), terms_split[[i]])
    sp_calls <- sub("^I\\(", "(", terms_split[[i]][has_sp_calls])
    out$joint_call[[i]] <- paste0(sp_calls, collapse = " * ")
    out$Ic[i] <- any(!has_sp_calls)
  }
  not_one <- apply(mm, 2, function(x) any(x != 1))
  out$Ic <- cumsum(out$Ic | not_one)
  out
}

# extract names of monotonic simplex parameters 
# @param spef output of tidy_spef
# @param use_id use the 'id' argument to construct simo labels?
# @return a character vector of length nrow(spef)
get_simo_labels <- function(spef, use_id = FALSE) {
  out <- named_list(spef$term)
  I <- which(lengths(spef$Imo) > 0)
  for (i in I) {
    # use the ID as label if specified
    out[[i]] <- ifelse(
      use_id & !is.na(spef$ids_mo[[i]]), spef$ids_mo[[i]],
      paste0(spef$coef[i], seq_along(spef$Imo[[i]]))
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

# names of grouping variables used in measurement error terms
get_me_groups <- function(x) {
  uni_me <- get_uni_me(x)
  out <- lapply(uni_me, eval2) 
  out <- ulapply(out, "[[", "gr")
  out[nzchar(out)]
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
  terms_comb <- rep(NA, length(terms_split))
  # loop over terms and add dummy variables
  for (i in seq_along(terms_split)) {
    replace_i <- grepl_expr(regex, terms_split[[i]])
    terms_i_replace <- terms_split[[i]][replace_i]
    dummies_i <- dummies[match(terms_i_replace, terms_replace)]
    terms_split[[i]][replace_i] <- dummies_i
    terms_comb[i] <- paste0(terms_split[[i]], collapse = ":") 
  }
  new_formula <- str2formula(terms_comb)
  attributes(new_formula) <- attributes(formula)
  out <- get_model_matrix(new_formula, data, ...)
  # recover original column names
  colnames(out) <- rename(colnames(out), dummies, terms_replace)
  out
}

# formula of variables used in special effects terms
sp_fake_formula <- function(...) {
  dots <- c(...)
  out <- vector("list", length(dots))
  for (i in seq_along(dots)) {
    tmp <- eval2(dots[[i]])
    out[[i]] <- all_vars(c(tmp$term, tmp$sdx, tmp$gr))
  }
  str2formula(unique(unlist(out)))
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
    x <- as.numeric(x) - 1
  } else if (all(is_wholenumber(x))) {
    min_value <- attr(x, "min")
    if (is.null(min_value)) {
      min_value <- min(x)
    }
    x <- x - min_value
  } else {
    stop2(
      "Monotonic predictors must be integers or ordered ",
      "factors. Error occurred for variable '", term$term, "'."
    )
  }
  as.array(x)
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
  c("mo", "me", "mi")
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
