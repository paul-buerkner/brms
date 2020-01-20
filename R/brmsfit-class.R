#' Class \code{brmsfit} of models fitted with the \pkg{brms} package
#' 
#' Models fitted with the \code{\link[brms:brms]{brms}} package are 
#' represented as a \code{brmsfit} object, which contains the posterior 
#' samples, model formula, Stan code, relevant data, and other information.
#' 
#' @name brmsfit-class
#' @aliases brmsfit
#' @docType class
#' 
#' @details 
#' See \code{methods(class = "brmsfit")} for an overview of available methods.
#' 
#' @slot formula A \code{\link{brmsformula}} object
#' @slot family (Deprecated) A \code{\link{brmsfamily}} object
#' @slot data A \code{data.frame} containing all variables used in the model
#' @slot data.name The name of \code{data} as specified by the user
#' @slot data2 An optional \code{list} of data objects which cannot be passed
#'   via \code{data} 
#' @slot model The model code in \pkg{Stan} language
#' @slot prior A \code{\link{brmsprior}} object containing
#'   information on the priors used in the model
#' @slot autocor (Deprecated) An \code{\link{cor_brms}} object containing 
#'   the autocorrelation structure if specified
#' @slot ranef A \code{data.frame} containing the group-level structure
#' @slot cov_ranef A \code{list} of customized group-level covariance matrices
#' @slot stanvars A \code{\link{stanvars}} object or \code{NULL}
#' @slot stan_funs A character string of length one or \code{NULL}
#' @slot criteria An empty \code{list} for adding model fit criteria
#'   after estimation of the model.
#' @slot fit An object of class \code{\link[rstan:stanfit]{stanfit}}
#'   among others containing the posterior samples
#' @slot exclude The names of the parameters for which samples are not saved
#' @slot algorithm The name of the algorithm used to fit the model
#' @slot version The versions of \pkg{brms} and \pkg{rstan} with 
#'   which the model was fitted
#' @slot file Optional name of a file in which the model object was stored in
#'   or loaded from
#' 
#' @seealso 
#'   \code{\link{brms}}, 
#'   \code{\link{brm}}, 
#'   \code{\link{brmsformula}}, 
#'   \code{\link{brmsfamily}}
#' 
NULL

# brmsfit class
brmsfit <- function(formula = NULL, family = NULL, data = data.frame(), 
                    data.name = "", data2 = list(), model = "", 
                    prior = empty_prior(), autocor = NULL, 
                    ranef = empty_ranef(), cov_ranef = NULL, 
                    criteria = list(), stanvars = NULL, stan_funs = NULL, 
                    fit = NA, exclude = NULL, algorithm = "sampling",
                    file = NULL) {
  version <- list(
    brms = utils::packageVersion("brms"),
    rstan = utils::packageVersion("rstan")
  )
  x <- nlist(
    formula, family, data, data.name, data2, model, prior,
    autocor, ranef, cov_ranef, stanvars, stan_funs, 
    fit, exclude, algorithm, version, file
  )
  class(x) <- "brmsfit"
  x
}

#' Checks if argument is a \code{brmsfit} object
#' 
#' @param x An \R object
#' 
#' @export
is.brmsfit <- function(x) {
  inherits(x, "brmsfit")
}

#' Checks if argument is a \code{brmsfit_multiple} object
#' 
#' @param x An \R object
#' 
#' @export
is.brmsfit_multiple <- function(x) {
  inherits(x, "brmsfit_multiple")
}
