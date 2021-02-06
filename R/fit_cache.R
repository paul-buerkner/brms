# Normalizes Stan code to avoid triggering refit after whitespace and
# comment changes in the generated code.
# In some distant future, StanC3 may provide its own normalizing functions,
# until then this is a set of regex hacks.
# @param x a string containing the Stan code
normalize_stancode <- function(x) {
  x <- as_one_character(x)
  # Remove single-line comments
  x <- gsub("//[^\n\r]*[\n\r]", " ", x)
  x <- gsub("//[^\n\r]*$", " ", x)
  # Remove multi-line comments
  x <- gsub("/\\*([^*]*(\\*[^/])?)*\\*/", " ", x)
  # Standardize whitespace (including newlines)
  x <- gsub("[[:space:]]+"," ", x)
  
  trimws(x)
}


#' Check whether a given cached fit can be used without refitting.
#' Use with \code{verbose = TRUE} to get additional info on how the stored
#' fit differs from the given data and code.
#' @param cached_fit old fit (e.g. loaded from file)
#' @param sdata new Stan data (result of a call to \code{\link{make_standata}}).
#'   pass \code{NULL} to avoid data check.
#' @param scode new Stan code (result of a call to \code{\link{make_stancode}}).
#'   pass \code{NULL} to avoid code check.
#' @param algorithm new algorithm. Pass \code{NULL} to avoid algorithm check.
#' @param silent if \code{TRUE}, no messages will be given
#' @param verbose if \code{TRUE} detailed report of the differences is printed
#'   to the console.
#' @return a boolean indicating whether a refit is needed
#' @export
brmsfit_needs_refit <- function(cached_fit, sdata, scode, algorithm,
                                silent = FALSE,  verbose = FALSE) {
  stopifnot(is.brmsfit(cached_fit))
  if(!is.null(scode)) {
    stopifnot(!is.null(cached_fit$model))
  }
  if(!is.null(sdata)) {
    stopifnot(!is.null(standata(cached_fit)))
    stopifnot(is.list(sdata))
  }
  if(!is.null(algorithm)) {
    stopifnot(!is.null(cached_fit$algorithm))
  }
  

  refit <- FALSE
  if(!is.null(scode)) {
    if(normalize_stancode(scode) != normalize_stancode(cached_fit$model)) {
      if(!silent) {
        message("Stan code changed beyond whitespace/comments.")
        if(verbose) {
          if(!requireNamespace(package, quietly = TRUE)) {
            message("Install package 'diffobj' to show code differences.")
          } else {
            print(
              diffobj::diffChr(scode, cached_fit$model, format = "ansi8")
            )
          }
        }
      }
      refit <- TRUE
    }
  }
  
  if(!is.null(sdata)) {
    sdata_equality_res <- all.equal(sdata, standata(cached_fit), 
                                    check.attributes = FALSE)
    if(!isTRUE(sdata_equality_res)) {
      if(!silent) {
        message("The processed data for Stan changed.\n")
        if(verbose) {
          print(sdata_equality_res)
        }
      }
      refit <- TRUE
    }
  }
  
  if(!is.null(algorithm)) {
    if(as_one_character(algorithm) != cached_fit$algorithm) {
      if(!silent) {
        message(paste0("Algorthm changed from '", cached_fit$algorithm, 
                       "' to '", algorithm, "'\n"))
      }
      refit <- TRUE
    }
  }
    
  refit
}
