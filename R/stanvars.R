#' User-defined variables passed to Stan
#'
#' Prepare user-defined variables to be passed to one of Stan's
#' program blocks. This is primarily useful for defining more complex
#' priors, for refitting models without recompilation despite
#' changing priors, or for defining custom Stan functions.
#'
#' @aliases stanvars
#'
#' @param x An \R object containing data to be passed to Stan.
#'   Only required if \code{block = 'data'} and ignored otherwise.
#' @param name Optional character string providing the desired variable
#'  name of the object in \code{x}. If \code{NULL} (the default)
#'  the variable name is directly inferred from \code{x}.
#' @param scode Line of Stan code to define the variable
#'  in Stan language. If \code{block = 'data'}, the
#'  Stan code is inferred based on the class of \code{x} by default.
#' @param block Name of one of Stan's program blocks in
#'  which the variable should be defined. Can be \code{'data'},
#'  \code{'tdata'} (transformed data), \code{'parameters'},
#'  \code{'tparameters'} (transformed parameters), \code{'model'},
#'  \code{'likelihood'} (part of the model block where the likelihood is given),
#'  \code{'genquant'} (generated quantities) or \code{'functions'}.
#' @param position Name of the position within the block where the
#'  Stan code should be placed. Currently allowed are \code{'start'}
#'  (the default) and \code{'end'} of the block.
#' @param pll_args Optional Stan code to be put into the header
#'  of \code{partial_log_lik} functions. This ensures that the variables
#'  specified in \code{scode} can be used in the likelihood even when
#'  within-chain parallelization is activated via \code{\link{threading}}.
#'
#' @return An object of class \code{stanvars}.
#'
#' @details
#' The \code{stanvar} function is not vectorized. Instead, multiple
#' \code{stanvars} objects can be added together via \code{+} (see Examples).
#'
#' Special attention is necessary when using \code{stanvars} to inject
#' code into the \code{'likelihood'} block while having \code{\link{threading}}
#' activated. In this case, your custom Stan code may need adjustments to ensure
#' correct observation indexing. Please investigate the generated Stan code via
#' \code{\link[brms:make_stancode.default]{make_stancode}} to see which adjustments are necessary in your case.
#'
#' @examples
#' bprior <- prior(normal(mean_intercept, 10), class = "Intercept")
#' stanvars <- stanvar(5, name = "mean_intercept")
#' make_stancode(count ~ Trt, epilepsy, prior = bprior,
#'               stanvars = stanvars)
#'
#' # define a multi-normal prior with known covariance matrix
#' bprior <- prior(multi_normal(M, V), class = "b")
#' stanvars <- stanvar(rep(0, 2), "M", scode = "  vector[K] M;") +
#'   stanvar(diag(2), "V", scode = "  matrix[K, K] V;")
#' make_stancode(count ~ Trt + zBase, epilepsy,
#'               prior = bprior, stanvars = stanvars)
#'
#' # define a hierachical prior on the regression coefficients
#' bprior <- set_prior("normal(0, tau)", class = "b") +
#'   set_prior("target += normal_lpdf(tau | 0, 10)", check = FALSE)
#' stanvars <- stanvar(scode = "real<lower=0> tau;",
#'                     block = "parameters")
#' make_stancode(count ~ Trt + zBase, epilepsy,
#'               prior = bprior, stanvars = stanvars)
#'
#' # ensure that 'tau' is passed to the likelihood of a threaded model
#' # not necessary for this example but may be necessary in other cases
#' stanvars <- stanvar(scode = "real<lower=0> tau;",
#'                     block = "parameters", pll_args = "real tau")
#' make_stancode(count ~ Trt + zBase, epilepsy,
#'               stanvars = stanvars, threads = threading(2))
#'
#' @export
stanvar <- function(x = NULL, name = NULL, scode = NULL,
                    block = "data", position = "start",
                    pll_args = NULL) {
  vblocks <- c(
    "data", "tdata", "parameters", "tparameters",
    "model", "genquant", "functions", "likelihood"
  )
  block <- match.arg(block, vblocks)
  vpositions <- c("start", "end")
  position <- match.arg(position, vpositions)
  if (block == "data") {
    if (is.null(x)) {
      stop2("Argument 'x' is required if block = 'data'.")
    }
    if (is.null(name)) {
      name <- deparse0(substitute(x))
    }
    name <- as_one_character(name)
    if (!is_equal(name, make.names(name)) || grepl("\\.", name)) {
      stop2("'", limit_chars(name, 30), "' is not ",
            "a valid variable name in Stan.")
    }
    if (is.null(scode)) {
      # infer scode from x
      if (is.integer(x)) {
        if (length(x) == 1L) {
          scode <- paste0("int ", name)
        } else {
          scode <- paste0("int ", name, "[", length(x), "]")
        }
      } else if (is.vector(x)) {
        if (length(x) == 1L) {
          scode <- paste0("real ", name)
        } else {
          scode <- paste0("vector[", length(x), "] ", name)
        }
      } else if (is.array(x)) {
        if (length(dim(x)) == 1L) {
          scode <- paste0("vector[", length(x), "] ", name)
        } else if (is.matrix(x)) {
          scode <- paste0("matrix[", nrow(x), ", ", ncol(x), "] ", name)
        }
      }
      if (is.null(scode)) {
        stop2(
          "'stanvar' could not infer the Stan code for an object ",
          "of class '", class(x), "'. Please specify the Stan code ",
          "manually via argument 'scode'."
        )
      }
      scode <- paste0(scode, ";")
    }
    if (is.null(pll_args)) {
      # infer pll_args from x
      pll_type <- str_if(block %in% c("data", "tdata"), "data ")
      if (is.integer(x)) {
        if (length(x) == 1L) {
          pll_type <- paste0(pll_type, "int")
        } else {
          pll_type <- paste0(pll_type, "array[] int")
        }
      } else if (is.vector(x)) {
        if (length(x) == 1L) {
          pll_type <- paste0(pll_type, "real")
        } else {
          pll_type <- paste0(pll_type, "vector")
        }
      } else if (is.array(x)) {
        if (length(dim(x)) == 1L) {
          pll_type <- paste0(pll_type, "vector")
        } else if (is.matrix(x)) {
          pll_type <- paste0(pll_type, "matrix")
        }
      }
      if (!is.null(pll_type)) {
        pll_args <- paste0(pll_type, " ", name)
      } else {
        # don't throw an error because most people will not use threading
        pll_args <- character(0)
      }
    }
  } else {
    x <- NULL
    if (is.null(name)) {
      name <- ""
    }
    name <- as_one_character(name)
    if (is.null(scode)) {
      stop2("Argument 'scode' is required if block is not 'data'.")
    }
    scode <- as.character(scode)
    pll_args <- as.character(pll_args)
  }
  if (position == "end" && block %in% c("functions", "data")) {
    stop2("Position '", position, "' is not sensible for block '", block, "'.")
  }
  out <- nlist(name, sdata = x, scode, block, position, pll_args)
  structure(setNames(list(out), name), class = "stanvars")
}

# take a subset of a stanvars object
# @param x a stanvars object
# @param ... conditions defining the desired subset
subset_stanvars <- function(x, ...) {
  x <- validate_stanvars(x)
  structure_not_null(x[find_elements(x, ...)], class = "stanvars")
}

# collapse Stan code provided in a stanvars object
collapse_stanvars <- function(x, block = NULL, position = NULL) {
  x <- validate_stanvars(x)
  if (!length(x)) {
    return(character(0))
  }
  if (!is.null(block)) {
    x <- subset_stanvars(x, block = block)
  }
  if (!is.null(position)) {
    x <- subset_stanvars(x, position = position)
  }
  if (!length(x)) {
    return("")
  }
  collapse(wsp(nsp = 2), ufrom_list(x, "scode"), "\n")
}

# collapse partial log-lik code provided in a stanvars object
collapse_stanvars_pll_args <- function(x) {
  x <- validate_stanvars(x)
  if (!length(x)) {
    return(character(0))
  }
  out <- ufrom_list(x, "pll_args")
  if (!length(out)) {
    return("")
  }
  collapse(", ", out)
}

# validate 'stanvars' objects
validate_stanvars <- function(x, stan_funs = NULL) {
  if (is.null(x)) {
    x <- empty_stanvars()
  }
  if (!is.stanvars(x)) {
    stop2("Argument 'stanvars' is invalid. See ?stanvar for help.")
  }
  if (length(stan_funs) > 0) {
    warning2("Argument 'stan_funs' is deprecated. Please use argument ",
             "'stanvars' instead. See ?stanvar for more help.")
    stan_funs <- as_one_character(stan_funs)
    x <- x + stanvar(scode = stan_funs, block = "functions")
  }
  x
}

# add new data to stanvars
# @param x a 'stanvars' object
# @param newdata2 a list with new 'data2' objects
# @return a 'stanvars' object
add_newdata_stanvars <- function(x, newdata2) {
  stopifnot(is.stanvars(x))
  stanvars_data <- subset_stanvars(x, block = "data")
  for (name in names(stanvars_data)) {
    if (name %in% names(newdata2)) {
      x[[name]]$sdata <- newdata2[[name]]
    }
  }
  x
}

#' @export
c.stanvars <- function(x, ...) {
  dots <- lapply(list(...), validate_stanvars)
  class(x) <- "list"
  out <- unlist(c(list(x), dots), recursive = FALSE)
  svnames <- names(out)[nzchar(names(out))]
  if (any(duplicated(svnames))) {
    stop2("Duplicated names in 'stanvars' are not allowed.")
  }
  structure(out, class = "stanvars")
}

#' @export
"+.stanvars" <- function(e1, e2) {
  c(e1, e2)
}

is.stanvars <- function(x) {
  inherits(x, "stanvars")
}

empty_stanvars <- function() {
  structure(list(), class = "stanvars")
}
