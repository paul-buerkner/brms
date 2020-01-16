#' Update \pkg{brms} models
#' 
#' This method allows to update an existing \code{brmsfit} object.
#' 
#' @param object An object of class \code{brmsfit}.
#' @param formula. Changes to the formula; for details see 
#'   \code{\link{update.formula}} and \code{\link{brmsformula}}.
#' @param newdata Optional \code{data.frame} 
#'   to update the model with new data.
#' @param recompile Logical, indicating whether the Stan model should 
#'   be recompiled. If \code{NULL} (the default), \code{update} tries
#'   to figure out internally, if recompilation is necessary. 
#'   Setting it to \code{FALSE} will cause all Stan code changing 
#'   arguments to be ignored. 
#' @param ... Other arguments passed to \code{\link{brm}}.
#'  
#' @details Sometimes, when updating the model formula, 
#'  it may happen that \R complains about a mismatch
#'  between \code{model frame} and \code{formula}.
#'  This error can be avoided by supplying your original data
#'  again via argument \code{newdata}.
#'  
#' @examples 
#' \dontrun{
#' fit1 <- brm(time | cens(censored) ~ age * sex + disease + (1|patient), 
#'             data = kidney, family = gaussian("log"))
#' summary(fit1)
#' 
#' ## remove effects of 'disease'
#' fit2 <- update(fit1, formula. = ~ . - disease)
#' summary(fit2)
#' 
#' ## remove the group specific term of 'patient' and
#' ## change the data (just take a subset in this example)
#' fit3 <- update(fit1, formula. = ~ . - (1|patient), 
#'                newdata = kidney[1:38, ])
#' summary(fit3)
#' 
#' ## use another family and add population-level priors
#' fit4 <- update(fit1, family = weibull(), inits = "0",
#'                prior = set_prior("normal(0,5)"))
#' summary(fit4)
#' }
#'
#' @export
update.brmsfit <- function(object, formula., newdata = NULL, 
                           recompile = NULL, ...) {
  dots <- list(...)
  testmode <- isTRUE(dots[["testmode"]])
  dots$testmode <- NULL
  object <- restructure(object)
  object$file <- NULL
  if (isTRUE(object$version$brms < "2.0.0")) {
    warning2("Updating models fitted with older versions of brms may fail.")
  }
  if ("data" %in% names(dots)) {
    # otherwise the data name cannot be found by substitute 
    stop2("Please use argument 'newdata' to update the data.")
  }
  if (!is.null(newdata)) {
    # TODO: update info stored in the families such as 'cats' or 'thres'
    dots$data <- newdata
    data.name <- substitute_name(newdata)
  } else {
    dots$data <- rm_attr(object$data, c("terms", "brmsframe"))
    data.name <- object$data.name
  }
  
  if (missing(formula.) || is.null(formula.)) {
    dots$formula <- object$formula
    if (!is.null(dots[["family"]])) {
      dots$formula <- bf(dots$formula, family = dots$family)
    } 
    if (!is.null(dots[["autocor"]])) {
      dots$formula <- bf(dots$formula, autocor = dots$autocor)
    }
  } else {
    # TODO: restructure updating of the model formula
    if (is.mvbrmsformula(formula.) || is.mvbrmsformula(object$formula)) {
      stop2("Updating formulas of multivariate models is not yet possible.")
    }
    if (is.brmsformula(formula.)) {
      nl <- get_nl(formula.)
    } else {
      formula. <- as.formula(formula.) 
      nl <- get_nl(formula(object))
    }
    family <- get_arg("family", formula., dots, object)
    autocor <- get_arg("autocor", formula., dots, object)
    dots$formula <- bf(formula., family = family, autocor = autocor, nl = nl)
    if (is_nonlinear(object)) {
      if (length(setdiff(all.vars(dots$formula$formula), ".")) == 0L) {
        dots$formula <- update(object$formula, dots$formula, mode = "keep")
      } else {
        dots$formula <- update(object$formula, dots$formula, mode = "replace")
        message("Argument 'formula.' will completely replace the ", 
                "original formula in non-linear models.")
      }
    } else {
      mvars <- all.vars(dots$formula$formula)
      mvars <- setdiff(mvars, c(names(object$data), "."))
      if (length(mvars) && is.null(newdata)) {
        stop2("New variables found: ", collapse_comma(mvars),
              "\nPlease supply your data again via argument 'newdata'.")
      }
      dots$formula <- update(formula(object), dots$formula)
    }
  }
  dots$formula <- validate_formula(dots$formula, data = dots$data)
  
  if (is.null(dots$prior)) {
    dots$prior <- object$prior
  } else {
    # update existing priors
    if (!is.brmsprior(dots$prior)) { 
      stop2("Invalid 'prior' argument.")
    }
    dots$prior <- rbind(dots$prior, object$prior)
    dupl_priors <- duplicated(dots$prior[, rcols_prior()])
    dots$prior <- dots$prior[!dupl_priors, ]
  }
  if (is.null(dots$sample_prior)) {
    dots$sample_prior <- attr(object$prior, "sample_prior")
    if (is.null(dots$sample_prior)) {
      has_prior_pars <- any(grepl("^prior_", parnames(object)))
      dots$sample_prior <- if (has_prior_pars) "yes" else "no"
    }
  }
  if (is.null(dots$save_ranef)) {
    dots$save_ranef <- isTRUE(attr(object$exclude, "save_ranef"))
  }
  if (is.null(dots$save_mevars)) {
    dots$save_mevars <- isTRUE(attr(object$exclude, "save_mevars"))
  }
  if (is.null(dots$save_all_pars)) {
    dots$save_all_pars <- isTRUE(attr(object$exclude, "save_all_pars"))
  }
  if (is.null(dots$knots)) {
    dots$knots <- attr(object$data, "knots")
  }
  arg_names <- c("data2", "cov_ranef", "stanvars", "stan_funs")
  old_args <- setdiff(arg_names, names(dots))
  dots[old_args] <- object[old_args]
  
  # update arguments controlling the sampling process
  if (is.null(dots$iter)) {
    # only keep old 'warmup' if also keeping old 'iter'
    dots$warmup <- first_not_null(dots$warmup, object$fit@sim$warmup)
  }
  dots$iter <- first_not_null(dots$iter, object$fit@sim$iter)
  dots$chains <- first_not_null(dots$chains, object$fit@sim$chains)
  dots$thin <- first_not_null(dots$thin, object$fit@sim$thin)
  control <- attr(object$fit@sim$samples[[1]], "args")$control
  control <- control[setdiff(names(control), names(dots$control))]
  dots$control[names(control)] <- control
  
  if (is.null(recompile)) {
    # only recompile if new and old stan code do not match
    new_stancode <- suppressMessages(do_call(make_stancode, dots))
    # stan code may differ just because of the version number (#288)
    new_stancode <- sub("^[^\n]+\n", "", new_stancode)
    old_stancode <- stancode(object, version = FALSE)
    recompile <- !is_equal(new_stancode, old_stancode)
    if (recompile) {
      message("The desired updates require recompiling the model") 
    }
  }
  recompile <- as_one_logical(recompile)
  if (recompile) {
    # recompliation is necessary
    dots$fit <- NA
    if (!testmode) {
      object <- do_call(brm, dots)
    }
  } else {
    # refit the model without compiling it again
    if (!is.null(dots$formula)) {
      object$formula <- dots$formula
      dots$formula <- NULL
    }
    bterms <- parse_bf(object$formula)
    object$data <- validate_data(dots$data, bterms = bterms)
    object$data2 <- validate_data2(dots$data2, bterms = bterms)
    object$family <- get_element(object$formula, "family")
    object$autocor <- get_element(object$formula, "autocor")
    object$ranef <- tidy_ranef(bterms, data = object$data)
    object$stanvars <- validate_stanvars(dots$stanvars)
    if (!is.null(dots$sample_prior)) {
      dots$sample_prior <- check_sample_prior(dots$sample_prior)
      attr(object$prior, "sample_prior") <- dots$sample_prior
    }
    object$exclude <- exclude_pars(
      bterms, data = object$data, ranef = object$ranef, 
      save_ranef = dots$save_ranef, save_mevars = dots$save_mevars,
      save_all_pars = dots$save_all_pars
    )
    if (!is.null(dots$algorithm)) {
      aopts <- c("sampling", "meanfield", "fullrank")
      algorithm <- match.arg(dots$algorithm, aopts)
      dots$algorithm <- object$algorithm <- algorithm
    } else if (!is.null(object$algorithm)) {
      dots$algorithm <- object$algorithm
    }
    if (!testmode) {
      dots$fit <- object
      object <- do_call(brm, dots)
    }
  }
  object$data.name <- data.name
  object
}

#' Update \pkg{brms} models based on multiple data sets
#' 
#' This method allows to update an existing \code{brmsfit_multiple} object.
#' 
#' @param object An object of class \code{brmsfit_multiple}.
#' @param formula. Changes to the formula; for details see 
#'   \code{\link{update.formula}} and \code{\link{brmsformula}}.
#' @param newdata List of \code{data.frames} to update the model with new data.
#'   Currently required even if the original data should be used.
#' @param ... Other arguments passed to \code{\link{update.brmsfit}}
#'   and \code{\link{brm_multiple}}.
#'  
#' @examples 
#' \dontrun{
#' library(mice)
#' imp <- mice(nhanes2)
#' 
#' # initially fit the model 
#' fit_imp1 <- brm_multiple(bmi ~ age + hyp + chl, data = imp, chains = 1)
#' summary(fit_imp1)
#' 
#' # update the model using fewer predictors
#' fit_imp2 <- update(fit_imp1, formula. = . ~ hyp + chl, newdata = imp)
#' summary(fit_imp2)
#' }
#'
#' @export
update.brmsfit_multiple <- function(object, formula., newdata = NULL, ...) {
  dots <- list(...)
  if ("data" %in% names(dots)) {
    # otherwise the data name cannot be found by substitute 
    stop2("Please use argument 'newdata' to update the data.")
  }
  if (is.null(newdata)) {
    stop2("'newdata' is required when updating a 'brmsfit_multiple' object.")
  }
  data.name <- substitute_name(newdata)
  if (inherits(newdata, "mids")) {
    require_package("mice", version = "3.0.0")
    newdata <- lapply(seq_len(newdata$m), mice::complete, data = newdata)
  } else if (!(is.list(newdata) && is.vector(newdata))) {
    stop2("'newdata' must be a list of data.frames.")
  }
  
  # update the template model using all arguments
  if (missing(formula.)) {
    formula. <- NULL
  }
  args <- c(nlist(object, formula., newdata = newdata[[1]]), dots)
  args$file <- NULL
  args$chains <- 0
  fit <- do_call(update.brmsfit, args)
  
  # arguments later passed to brm_multiple
  args <- c(nlist(fit, data = newdata), dots)
  # update arguments controlling the sampling process
  # they cannot be accessed directly from the template model 
  # as it does not contain any samples (chains = 0)
  if (is.null(args$iter)) {
    # only keep old 'warmup' if also keeping old 'iter'
    args$warmup <- first_not_null(args$warmup, object$fit@sim$warmup)
  }
  if (is.null(args$chains)) {
    # chains were combined across all submodels
    args$chains <- object$fit@sim$chains / max(NROW(object$rhats), 1)
  }
  args$iter <- first_not_null(args$iter, object$fit@sim$iter)
  args$thin <- first_not_null(args$thin, object$fit@sim$thin)
  control <- attr(object$fit@sim$samples[[1]], "args")$control
  control <- control[setdiff(names(control), names(args$control))]
  args$control[names(control)] <- control
  args$recompile <- NULL
  
  out <- do_call(brm_multiple, args)
  out$data.name <- data.name
  out
}
