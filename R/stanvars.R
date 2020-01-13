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
#'  the variable name is directly infered from \code{x}.
#' @param scode Line of Stan code to define the variable
#'  in Stan language. If \code{block = 'data'}, the
#'  Stan code is inferred based on the class of \code{x} by default.
#' @param block Name of one of Stan's program blocks in
#'  which the variable should be defined. Can be \code{'data'},
#'  \code{'tdata'} (transformed data), \code{'parameters'},
#'  \code{'tparameters'} (transformed parameters), \code{'model'},
#'  \code{'genquant'} (generated quantities) or \code{'functions'}.
#'  
#' @return An object of class \code{stanvars}.
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
#' @export
stanvar <- function(x = NULL, name = NULL, scode = NULL,
                    block = "data") {
  vblocks <- c(
    "data", "tdata", "parameters", "tparameters", 
    "model", "genquant", "functions"
  )
  block <- match.arg(block, vblocks)
  if (block == "data") {
    if (is.null(x)) {
      stop2("Argument 'x' is required if block = 'data'.")
    }
    if (is.null(name)) {
      name <- deparse(substitute(x))
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
      scode <- paste0("  ", scode, ";")
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
  }
  out <- nlist(name, sdata = x, scode, block)
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
collapse_stanvars <- function(x, block = NULL) {
  x <- validate_stanvars(x)
  if (!length(x)) {
    return(character(0))
  }
  if (!is.null(block)) {
    x <- subset_stanvars(x, block = block)
  }
  if (!length(x)) {
    return("")
  }
  collapse(ulapply(x, "[[", "scode"), "\n")
}

validate_stanvars <- function(x) {
  if (!is.null(x) && !is.stanvars(x)) {
    stop2("Argument 'stanvars' is invalid. See ?stanvar for help.")
  }
  x
}

# add new stanvars to a brmsfit object
add_new_stanvars <- function(x, newdata2) {
  stopifnot(is.brmsfit(x))
  newdata2 <- validate_data2(newdata2)
  stanvars_data <- subset_stanvars(x$stanvars, block = "data")
  for (name in names(stanvars_data)) {
    if (name %in% names(newdata2)) {
      x$stanvars[[name]]$sdata <- newdata2[[name]]
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
