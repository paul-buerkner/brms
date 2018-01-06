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
#' To find out whether each of the original models converged,
#' investigate \code{fit$rhats}, where \code{fit} 
#' denotes the output of \code{brm_multiple}.
#' 
#' @return If \code{combine = TRUE} a \code{brmsfit_multiple} object,
#' which inherits from class \code{brmsfit} and behaves essentially 
#' the same. If \code{combine = FALSE} a list of \code{brmsfit} objects. 
#' 
#' @examples
#' \dontrun{
#' library(mice)
#' imp <- mice(nhanes2)
#' 
#' # fit the model using mice and lm
#' fit_imp1 <- with(lm(bmi~age+hyp+chl), data = imp)
#' summary(pool(fit_imp1))
#' 
#' # fit the model using brms
#' fit_imp2 <- brm_multiple(bmi~age+hyp+chl, data = imp, chains = 1)
#' summary(fit_imp2)
#' plot(fit_imp2, pars = "^b_")
#' # investigate convergence of the original models
#' fit_imp2$rhats
#' }
#' 
#' @export
brm_multiple <- function(formula, data, combine = TRUE, ...) {
  combine <- as_one_logical(combine)
  data.name <- substr(collapse(deparse(substitute(data))), 1, 50)
  if (inherits(data, "mids")) {
    require_package("mice")
    data <- lapply(seq_len(data$m), mice::complete, x = data)
  } else if (!(is.list(data) && is.vector(data))) {
    stop2("'data' must be a list of data.frames.")
  }
  fits <- vector("list", length(data))
  message("Fitting imputed model 1")
  fits[[1]] <- brm(formula, data = data[[1]], ...)
  fits[[1]]$data.name <- data.name
  rhats <- data.frame(as.list(rhat(fits[[1]])))
  if (any(rhats > 1.1)) {
    warning2("Imputed model 1 did not converge.")
  }
  for (i in seq_along(data)[-1]) {
    message("Fitting imputed model ", i)
    fits[[i]] <- update(fits[[1]], newdata = data[[i]], ...)
    rhat_i <- data.frame(as.list(rhat(fits[[i]])))
    if (any(rhat_i > 1.1)) {
      warning2("Imputed model ", i, " did not converge.")
    }
    rhats <- rbind(rhats, rhat_i)
  }
  if (combine) {
    fits <- combine_models(mlist = fits, check_data = FALSE)
    fits$rhats <- rhats
    class(fits) <- c("brmsfit_multiple", class(fits))
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
