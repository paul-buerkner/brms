#' Init definitions for **brms** models
#'
#' Define how initial values for specific parameters are generated.
#'
#' @inheritParams set_prior
#' @param distribution A character string specifying the distribution of the initial values
#'
#' @return An object of class `brmsinits` to be used in the `init` argument of [brm]
#' @export
#'
#' @examples
#' \dontrun{
#' inits <- set_inits("normal(0, 1)", class = "Intercept", coef = "mu") +
#'         set_inits("uniform(-1, 1)", class = "b", coef = "mu")
#' # use the inits in a brm call
#' fit <- brm(count ~ Trt + zAge, epilepsy, poisson(), init = inits)
#' }
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


# validate arguments passed to 'set_inits'
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

# Internal function for generating a list of inits to pass to stan from a
# brmsinits object created from set_inits()
# @param binits A brmsinits object
# @param bterms A brmsterms object
# @param data The data used in the model
# @param sdata The stan data list
# @return A list of inits to pass to stan
.inits_fun <- function(binits, bterms, data, sdata) {
  # TODO: check if inits are properly specified (similar to how the priors are checked)
  pars <- paste0(binits$dpar, binits$nlpar)
  sep <- ifelse(pars == "", "", "_")
  # temporary - works for Intercept and b, but not for sd, z, etc; needs to be generalized by using code from .stancode
  binits$stanpars <- paste0(binits$class, sep, pars)
  # get the information typically used in the parameters block of stancode
  info <- par_info(bterms, data)

  dims <- sdata[info$b_dim_name]
  dims <- ifelse(is.na(info$b_dim_name), 1, dims)
  prefixes <- ifelse(info$b_type == "real", "", "array(")
  suffixes <- ifelse(info$b_type == "real", "", ")") # here we would add dimensions as well

  # construct the call for generating inits for each row of binits
  out <- list()
  for (i in 1:nrow(binits)) {
    idx <- which(info$b_par == binits$stanpars[[i]])
    pinfo <- info[idx, ]
    dist <- parse_dist(binits$distribution[[i]])
    args <- paste0(dist$args, collapse = ", ")
    prefix <- prefixes[idx]
    suffix <- suffixes[idx]
    dim <- dims[[idx]]
    call <- glue('{prefix}{dist$fun}({dim}, {args}){suffix}')
    call <- parse(text = call)
    out[[binits$stanpars[[i]]]] <- eval(call)
  }

  out
}


# combine multiple brmsinits objects into one brmsinits (code almost identical to
# c.brmsprior)
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


# takes a character vector like 'normal(0, 1)' and returns a list with the
# r* function and its arguments
# to do - more careful checks of the passed format?
parse_dist <- function(x) {
  x <- as_one_character(x)
  x <- parse(text = x)[[1]]
  dist <- as.character(x[[1]])
  args <- as.list(x[-1])
  args <- lapply(args, function(x) {
    tmp <- as.character(x)
    as.numeric(collapse(tmp))
  })
  fun <- to_rfun(dist)
  nlist(fun, args)
}

# takes a character string and returns the corresponding r random generation
# function
to_rfun <- function(x) {
  x <- as_one_character(x)
  # TODO expandlist
  dists <- c(normal = 'norm', poisson = 'pois', binomial = 'binom',
             inv_gamma = 'invgamma', lognormal = 'lnorm', exponential = 'exp',
             uniform = 'unif')
  out <- dists[x]
  if (is.null(out) || is.na(out)) {
    out <- x
  }
  paste0("r", out)
}

par_info <- function(bterms, data, ...) {
  UseMethod("par_info")
}

#' @export
par_info.brmsterms <- function(bterms, data, ...) {
  out <- list()
  for (par in names(bterms$dpars)) {
    info <- par_info_fe(bterms$dpars[[par]], data)
    info <- as.data.frame(info)
    out <- rbind(out, info)
  }
  out
}

# internal function for extracting information about the population-effects parameters
# that is typically part of the parameters block of the stan code
# @param bterms A brmsterms object
# @param data The data used in the model
# @return A list with the following elements:
#   - b_type: the type of the parameter (real, vector, array)
#   - b_dim_name: the name of the dimension of the parameter (should match in standata)
#   - b_par: the name of the parameter in stan
# @details
#   if a parameter is described as vector[Kc_sigma] b_sigma, the output will be:
#   list(b_type = "vector", b_dim_name = "Kc_sigma", b_par = "b_sigma")
par_info_fe <- function(bterms, data) {
  out <- list()
  family <- bterms$family
  fixef <- colnames(data_fe(bterms, data)$X)
  center_X <- stan_center_X(bterms)
  ct <- str_if(center_X, "c")
  # remove the intercept from the design matrix?
  if (center_X) {
    fixef <- setdiff(fixef, "Intercept")
  }
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  resp <- usc(px$resp)

  out <- list()
  if (length(fixef)) {
    out$b_type <- "vector"
    out$b_dim_name <- glue("K{ct}{p}")
    out$b_par <- glue("b{p}")
  }

  if (center_X) {
    c(out$b_type) <- "real"
    c(out$b_dim_name) <- NA
    c(out$b_par) <- glue("Intercept{p}")
  }
  out
}