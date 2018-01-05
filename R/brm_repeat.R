#' Run the same \pkg{brms} model on multiple datasets
#' 
#' Run the same \pkg{brms} model on multiple datasets and
#' then combine the results into one fitted model object. 
#' This is useful in particular for multiple missing value imputation,
#' where the same model is fitted on multiple imputed data sets.
#' 
#' @inheritParams brm
#' @param data A list of data.frames each of which will
#'   be used to fit a separate model. Alternatively, a 
#'   \code{mids} object from the \pkg{mice} package.
#' @param combine Logical; Indicates if the fitted models
#'   should be combined into a single fitted model object
#'   via \code{\link{combine_models}}. Defaults to \code{TRUE}.
#' @param ... Further arguments passed to \code{\link{brm}}.
#' 
#' @details The combined model may issue false positive 
#' convergence warnings, as the MCMC chains corresponding
#' to different datasets may not necessarily overlap, even 
#' if each of the original models did converge.
#' If you want to make sure that the convergence
#' warnings are not caused by real problems in some of the 
#' original models, run \code{brm_repeat} with \code{combine = FALSE},
#' check for convergence manually, and then combine the models
#' via \code{combine_models{mlist = fits, check_data = FALSE}},
#' where \code{fits} denotes the output of \code{brms_repeat}.
#' 
#' @return If \code{combine = TRUE} a single \code{brmsfit} object. 
#' If \code{combine = FALSE} a list of \code{brmsfit} objects. 
#' 
#' @examples
#' \dontrun{
#' library(mice)
#' imp <- mice(nhanes2)
#' 
#' # fit the model using mice and lm
#' fit_imp1 <- with(lm(bmi~age+hyp+chl), data = imp)
#' summary(pool(fit1))
#' 
#' # fit the model using brms
#' fit_imp2 <- brm_repeat(bmi~age+hyp+chl, data = imp, chains = 1)
#' summary(fit_imp2)
#' plot(fit_imp2, pars = "^b_")
#' }
#' 
#' @export
brm_repeat <- function(formula, data, combine = TRUE, ...) {
  dots <- list(...)
  combine <- as_one_logical(combine)
  data.name <- substr(collapse(deparse(substitute(data))), 1, 50)
  if (inherits(data, "mids")) {
    require_package("mice")
    data <- lapply(seq_len(data$m), mice::complete, x = data)
  } else if (!(is.list(data) && is.vector(data))) {
    stop2("'data' must be a list of data.frames.")
  }
  fits <- vector("list", length(data))
  args <- c(list(formula, data = data[[1]]), dots)
  fits[[1]] <- do.call(brm, args)
  fits[[1]]$data.name <- data.name
  for (i in seq_along(data)[-1]) {
    args <- c(list(fits[[1]], newdata = data[[i]]), dots)
    fits[[i]] <- do.call(update, args)
  }
  if (combine) {
    fits <- combine_models(mlist = fits, check_data = FALSE)
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
