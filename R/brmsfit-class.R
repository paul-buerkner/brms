#' Class \code{brmsfit} of models fitted with the \pkg{brms} package
#'
#' Models fitted with the \code{\link[brms:brms-package]{brms}} package are
#' represented as a \code{brmsfit} object, which contains the posterior
#' draws (samples), model formula, Stan code, relevant data, and other information.
#'
#' @name brmsfit-class
#' @aliases brmsfit
#' @docType class
#'
#' @details
#' See \code{methods(class = "brmsfit")} for an overview of available methods.
#'
#' @slot formula A \code{\link{brmsformula}} object.
#' @slot data A \code{data.frame} containing all variables used in the model.
#' @slot data2 A \code{list} of data objects which cannot be passed
#'   via \code{data}.
#' @slot prior A \code{\link{brmsprior}} object containing
#'   information on the priors used in the model.
#' @slot stanvars A \code{\link{stanvars}} object.
#' @slot model The model code in \pkg{Stan} language.
#' @slot ranef A \code{data.frame} containing the group-level structure.
#' @slot exclude The names of the parameters for which draws are not saved.
#' @slot algorithm The name of the algorithm used to fit the model.
#' @slot backend The name of the backend used to fit the model.
#' @slot threads An object of class `brmsthreads` created by
#'   \code{\link{threading}}.
#' @slot opencl An object of class `brmsopencl` created by \code{\link{opencl}}.
#' @slot stan_args Named list of additional control arguments that were passed
#'   to the Stan backend directly.
#' @slot fit An object of class \code{\link[rstan:stanfit-class]{stanfit}}
#'   among others containing the posterior draws.
#' @slot criteria An empty \code{list} for adding model fit criteria
#'   after estimation of the model.
#' @slot file Optional name of a file in which the model object was stored in
#'   or loaded from.
#' @slot version The versions of \pkg{brms} and \pkg{rstan} with
#'   which the model was fitted.
#' @slot family (Deprecated) A \code{\link{brmsfamily}} object.
#' @slot autocor (Deprecated) An \code{\link{cor_brms}} object containing
#'   the autocorrelation structure if specified.
#' @slot cov_ranef (Deprecated) A \code{list} of customized group-level
#'   covariance matrices.
#' @slot stan_funs (Deprecated) A character string of length one or \code{NULL}.
#' @slot data.name (Deprecated) The name of \code{data} as specified by the user.
#'
#' @seealso
#'   \code{\link{brms}},
#'   \code{\link{brm}},
#'   \code{\link{brmsformula}},
#'   \code{\link{brmsfamily}}
#'
NULL

# brmsfit class
brmsfit <- function(formula = NULL, data = data.frame(), prior = empty_prior(),
                    data2 = list(), stanvars = NULL, model = "",
                    ranef = empty_ranef(), save_pars = NULL,
                    algorithm = "sampling", backend = "rstan",
                    threads = threading(), opencl = opencl(),
                    stan_args = list(), fit = NULL, criteria = list(),
                    file = NULL, family = NULL, autocor = NULL,
                    cov_ranef = NULL, stan_funs = NULL, data.name = "") {
  version <- list(
    brms = utils::packageVersion("brms"),
    rstan = utils::packageVersion("rstan"),
    stanHeaders = utils::packageVersion("StanHeaders")
  )
  if (backend == "cmdstanr") {
    require_package("cmdstanr")
    version$cmdstanr <- utils::packageVersion("cmdstanr")
    version$cmdstan <- as.package_version(cmdstanr::cmdstan_version())
  }
  x <- nlist(
    formula, data, prior, data2, stanvars, model, ranef,
    save_pars, algorithm, backend, threads, opencl, stan_args, fit, criteria,
    file, version, family, autocor, cov_ranef, stan_funs, data.name
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

is.stanfit <- function(x) {
  inherits(x, "stanfit")
}

