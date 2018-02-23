#' User-defined variables passed to Stan
#' 
#' Prepare user-defined variables to be passed to Stan's data block.
#' This is primarily useful for defning more complex priors and 
#' for refitting models without recompilation despite changing priors.
#' 
#' @aliases stanvars
#'  
#' @param x An \R object containing data to be passed to Stan.
#' @param name Optinal character string providing the desired variable 
#'  name of the object in \code{x}. If \code{NULL} (the default)
#'  the variable name is directly infered from \code{x}.
#' @param scode Optional line of Stan code to define the variable
#'  in Stan's data block. If \code{NULL} (the default), the
#'  Stan code is inferred based on the class of \code{x}.
#'  
#' @return An object of class \code{stanvars}.
#' 
#' @examples 
#' bprior <- prior(normal(mean_intercept, 10), class = "Intercept")
#' stanvars <- stan_var(5, name = "mean_intercept")
#' make_stancode(count ~ Trt, epilepsy, prior = bprior, 
#'               stan_vars = stanvars)
#'               
#' # define a multi-normal prior with known covariance matrix
#' bprior <- prior(multi_normal(M, V), class = "b")
#' stanvars <- stan_var(2L, "KV") + 
#'   stan_var(rep(0, 2), "M", scode = "  vector[KV] M;") +
#'   stan_var(diag(2), "V", scode = "  matrix[KV, KV] V;") 
#' make_stancode(count ~ Trt + log_Base4_c, epilepsy,
#'               prior = bprior, stan_vars = stanvars)
#' 
#' @export
stan_var <- function(x, name = NULL, scode = NULL) {
  if (is.null(name)) {
    name <- deparse(substitute(x))
  }
  name <- as_one_character(name)
  if (is.null(scode)) {
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
        "'stan_vars' could not infer the Stan code for an object ",
        "of class '", class(x), "'. Please specify the Stan code ",
        "manually via argument 'scode'."
      )
    }
    scode <- paste0("  ", scode, ";")
  } else {
    scode <- as_one_character(scode)
  }
  out <- nlist(name, sdata = x, scode, block = "data")
  structure(setNames(list(out), name), class = "stanvars")
}

collapse_stanvars <- function(x) {
  x <- validate_stanvars(x)
  if (!length(x)) {
    return(character(0))
  }
  collapse(ulapply(x, "[[", "scode"), "\n")
}

validate_stanvars <- function(x) {
  if (!is.null(x) && !is.stanvars(x)) {
    stop2("Argument 'stan_vars' is invalid. See ?stanvars for help.")
  }
  x
}

#' @export
"c.stanvars" <- function(x, ...) {
  dots <- lapply(list(...), validate_stanvars)
  class(x) <- "list"
  out <- unlist(c(list(x), dots), recursive = FALSE)
  structure(out, class = "stanvars")
}

#' @export
"+.stanvars" <- function(e1, e2) {
  c(e1, e2)
}

is.stanvars <- function(x) {
  inherits(x, "stanvars")
}
