#' Run the same \pkg{brms} model on multiple datasets
#' 
#' Run the same \pkg{brms} model on multiple datasets and then combine the
#' results into one fitted model object. This is useful in particular for
#' multiple missing value imputation, where the same model is fitted on multiple
#' imputed data sets. Models can be run in parallel using the \pkg{future}
#' package.
#' 
#' @inheritParams brm
#' @param data A list of data.frames each of which will be used to fit a
#'   separate model. Alternatively, a \code{mids} object from the \pkg{mice}
#'   package.
#' @param combine Logical; Indicates if the fitted models should be combined
#'   into a single fitted model object via \code{\link{combine_models}}.
#'   Defaults to \code{TRUE}.
#' @param fit An instance of S3 class \code{brmsfit_multiple} derived from a
#'   previous fit; defaults to \code{NA}. If \code{fit} is of class
#'   \code{brmsfit_multiple}, the compiled model associated with the fitted
#'   result is re-used and all arguments modifying the model code or data are
#'   ignored. It is not recommended to use this argument directly, but to call
#'   the \code{\link[brms:update.brmsfit_multiple]{update}} method, instead.
#' @param ... Further arguments passed to \code{\link{brm}}.
#' 
#' @details The combined model may issue false positive convergence warnings, as
#'   the MCMC chains corresponding to different datasets may not necessarily
#'   overlap, even if each of the original models did converge. To find out
#'   whether each of the original models converged, investigate
#'   \code{fit$rhats}, where \code{fit} denotes the output of
#'   \code{brm_multiple}.
#' 
#' @return If \code{combine = TRUE} a \code{brmsfit_multiple} object, which
#'   inherits from class \code{brmsfit} and behaves essentially the same. If
#'   \code{combine = FALSE} a list of \code{brmsfit} objects.
#' 
#' @examples
#' \dontrun{
#' library(mice)
#' imp <- mice(nhanes2)
#' 
#' # fit the model using mice and lm
#' fit_imp1 <- with(lm(bmi ~ age + hyp + chl), data = imp)
#' summary(pool(fit_imp1))
#' 
#' # fit the model using brms
#' fit_imp2 <- brm_multiple(bmi ~ age + hyp + chl, data = imp, chains = 1)
#' summary(fit_imp2)
#' plot(fit_imp2, pars = "^b_")
#' # investigate convergence of the original models
#' fit_imp2$rhats
#' 
#' # use the future package for parallelization
#' library(future)
#' plan(multiprocess)
#' fit_imp3 <- brm_multiple(bmi~age+hyp+chl, data = imp, chains = 1)
#' summary(fit_imp3)
#' }
#' 
#' @export
brm_multiple <- function(formula, data, family = gaussian(), prior = NULL, 
                         autocor = NULL, cov_ranef = NULL, 
                         sample_prior = c("no", "yes", "only"), 
                         sparse = NULL, knots = NULL, stanvars = NULL,
                         stan_funs = NULL, combine = TRUE, fit = NA,
                         seed = NA, file = NULL, ...) {
  combine <- as_one_logical(combine)
  if (!is.null(file)) {
    # optionally load saved model object
    if (!combine) {
      stop2("Cannot use 'file' if 'combine' is FALSE.")
    }
    file <- paste0(as_one_character(file), ".rds")
    x <- suppressWarnings(try(readRDS(file), silent = TRUE))
    if (!is(x, "try-error")) {
      if (!is.brmsfit_multiple(x)) {
        stop2("Object loaded via 'file' is not of class 'brmsfit_multiple'.")
      }
      return(x)
    }
  }
  
  data.name <- substitute_name(data)
  if (inherits(data, "mids")) {
    require_package("mice", version = "3.0.0")
    data <- lapply(seq_len(data$m), mice::complete, data = data)
  } else if (!is_data_list(data)) {
    stop2("'data' must be a list of data.frames.")
  }
  
  if (is.brmsfit(fit)) {
    # avoid complications when updating the model
    class(fit) <- setdiff(class(fit), "brmsfit_multiple")
  } else {
    args <- nlist(
      formula, data = data[[1]], family, prior, autocor, cov_ranef,
      sample_prior, sparse, knots, stanvars, stan_funs, seed, ...
    )
    args$chains <- 0
    fit <- do_call(brm, args)
  }
  
  fits <- futures <- rhats <- vector("list", length(data))
  for (i in seq_along(data)) {
    futures[[i]] <- future::future(
      update(fit, newdata = data[[i]], recompile = FALSE, ...),
      packages = "brms"
    )
  }
  for (i in seq_along(data)) {
    message("Fitting imputed model ", i)
    fits[[i]] <- future::value(futures[[i]]) 
    rhats[[i]] <- data.frame(as.list(rhat(fits[[i]])))
    if (any(rhats[[i]] > 1.1, na.rm = TRUE)) {
      warning2("Imputed model ", i, " did not converge.")
    }
  }
  if (combine) {
    fits <- combine_models(mlist = fits, check_data = FALSE)
    fits$data.name <- data.name
    fits$rhats <- do_call(rbind, rhats)
    class(fits) <- c("brmsfit_multiple", class(fits))
  }
  if (!is.null(file)) {
    saveRDS(fits, file = file)
  }
  fits
}

#' Combine Models fitted with \pkg{brms}
#' 
#' Combine multiple \code{brmsfit} objects, which fitted the same model.
#' This is usefuly for instance when having manually run models in parallel.
#' 
#' @param ... One or more \code{brmsfit} objects.
#' @param mlist Optional list of one or more \code{brmsfit} objects.
#' @param check_data Logical; indicates if the data should be checked
#'   for being the same across models (defaults to \code{TRUE}).
#'   Setting it to \code{FALSE} may be useful for instance
#'   when combining models fitted on multiple imputed data sets.
#'   
#' @details This function just takes the first model and replaces 
#'   its \code{stanfit} object (slot \code{fit}) by the combined 
#'   \code{stanfit} objects of all models.
#'   
#' @return A \code{brmsfit} object.
#' 
#' @export
combine_models <- function(..., mlist = NULL, check_data = TRUE) {
  models <- c(list(...), mlist)
  check_data <- as_one_logical(check_data)
  if (!length(models)) {
    stop2("No models supplied to 'combine_models'.")
  }
  for (i in seq_along(models)) {
    if (!is.brmsfit(models[[i]])) {
      stop2("Model ", i, " is no 'brmsfit' object.")
    }
    models[[i]] <- restructure(models[[i]])
  }
  ref_formula <- formula(models[[1]])
  ref_pars <- parnames(models[[1]])
  ref_mf <- model.frame(models[[1]]) 
  for (i in seq_along(models)[-1]) {
    if (!is_equal(formula(models[[i]]), ref_formula)) {
      stop2("Models 1 and ", i, " have different formulas.")
    }
    if (!is_equal(parnames(models[[i]]), ref_pars)) {
      stop2("Models 1 and ", i, " have different parameters.")
    }
    if (check_data && !is_equal(model.frame(models[[i]]), ref_mf)) {
      stop2(
        "Models 1 and ", i, " have different data. ", 
        "Set 'check_data' to FALSE to turn off checking of the data."
      )
    }
  }
  sflist <- lapply(models, "[[", "fit")
  models[[1]]$fit <- rstan::sflist2stanfit(sflist)
  models[[1]]
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

# validity check for data input of 'brm_multiple'
is_data_list <- function(x) {
  is.list(x) && is.vector(x)
}

warn_brmsfit_multiple <- function(x) {
  if (is.brmsfit_multiple(x)) {
    warning2(
      "Using only the first imputed data set. Please interpret the results ", 
      "with caution until a more principled approach has been implemented."
    )
  }
  invisible(x)
}
