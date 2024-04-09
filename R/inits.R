set_inits <- function(distribution, class = "b", coef = "", group = "",
                      dpar = "", nlpar = "") {
  input <- nlist(distribution, class, coef, group, dpar, nlpar)
  input <- try(as.data.frame(input), silent = TRUE)
  if (is_try_error(input)) {
    stop2("Processing arguments of 'set_inits' has failed:\n", input)
  }
  out <- vector("list", nrow(input))
  for (i in seq_along(out)) {
    out[[i]] <- do_call(.set_inits, input[i, ])
  }
  Reduce("+", out)
}


# validate arguments passed to 'set_prior'
.set_inits <- function(distribution, class, coef, group,
                       dpar, nlpar) {
  distribution <- as_one_character(distribution)
  class <- as_one_character(class)
  group <- as_one_character(group)
  coef <- as_one_character(coef)
  dpar <- as_one_character(dpar)
  nlpar <- as_one_character(nlpar)
  if (dpar == "mu") {
    # distributional parameter 'mu' is currently implicit
    dpar <- ""
  }
  out <- data.frame(distribution, class, coef, group, dpar, nlpar)
  class(out) <- c("brmsinits", "data.frame")
  out
}


# combine multiple brmsinits objects into one brmsinits
#' @export
c.brmsinits <- function(x, ..., replace = FALSE) {
  dots <- list(...)
  if (all(sapply(dots, is.brmsinits))) {
    replace <- as_one_logical(replace)
    # don't use 'c()' here to avoid creating a recursion
    out <- do_call(rbind, list(x, ...))
    if (replace) {
      # update duplicated inits
      out <- unique(out, fromLast = TRUE)
    }
  } else {
    if (length(dots)) {
      stop2("Cannot add '", class(dots[[1]])[1], "' objects to the inits")
    }
    out <- c(as.data.frame(x))
  }
  out
}

#' @export
"+.brmsinits" <- function(e1, e2) {
  if (is.null(e2)) {
    return(e1)
  }
  if (!is.brmsinits(e2)) {
    stop2("Cannot add '", class(e2)[1], "' objects to the inits")
  }
  c(e1, e2)
}

is.brmsinits <- function(x) {
  inherits(x, "brmsinits")
}