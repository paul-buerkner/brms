# unless otherwise specified, functions return a named list
# of Stan code snippets to be pasted together later on

# define custom Stan functions
# TODO: Figure out how to *uniquely* include Stan functions that are not stored
# in inst/chunks but rather defined programatically on the fly (e.g., in
# stan_ordinal_lpmf). Once this is done, stan_global_defs can be removed.
stan_global_defs <- function(bterms) {
  stopifnot(is.anybrmsterms(bterms))
  families <- family_names(bterms)
  links <- family_info(bterms, "link")
  unique_combs <- !duplicated(paste0(families, ":", links))
  families <- families[unique_combs]
  links <- links[unique_combs]
  out <- list()
  is_ordinal <- ulapply(families, is_ordinal)
  if (any(is_ordinal)) {
    ord_fams <- families[is_ordinal]
    ord_links <- links[is_ordinal]
    for (i in seq_along(ord_fams)) {
      if (has_extra_cat(ord_fams[i])) {
        str_add(out$fun) <- stan_hurdle_ordinal_lpmf(ord_fams[i], ord_links[i])
      } else {
        str_add(out$fun) <- stan_ordinal_lpmf(ord_fams[i], ord_links[i])
      }
    }
  }
  out
}

# link function in Stan language
# @param link name of the link function
# @param transform actually apply the link function?
stan_link <- function(link, transform = TRUE) {
  transform <- as_one_logical(transform %||% FALSE)
  if (!transform) {
    # we have a Stan lpdf that applies the link automatically
    # or we have a non-linear parameter that has no link function
    return("")
  }
  out <- switch(
    link,
    identity = "",
    log = "log",
    logm1 = "logm1",
    inverse = "inv",
    sqrt = "sqrt",
    "1/mu^2" = "inv_square",
    logit = "logit",
    probit = "inv_Phi",
    probit_approx = "inv_Phi",
    cloglog = "cloglog",
    cauchit = "cauchit",
    tan_half = "tan_half",
    log1p = "log1p",
    softplus = "log_expm1",
    squareplus = "inv_squareplus",
    softit = "softit"
  )
  out
}

# inverse link in Stan language
# @param link name of the link function
# @param transform actually apply the inv_link function?
stan_inv_link <- function(link, transform = TRUE) {
  transform <- as_one_logical(transform %||% FALSE)
  if (!transform) {
    # we have a Stan lpdf that applies the inv_link automatically
    # or we have a non-linear parameter that has no link function
    return("")
  }
  out <- switch(
    link,
    identity = "",
    log = "exp",
    logm1 = "expp1",
    inverse = "inv",
    sqrt = "square",
    "1/mu^2" = "inv_sqrt",
    logit = "inv_logit",
    probit = "Phi",
    probit_approx = "Phi_approx",
    cloglog = "inv_cloglog",
    cauchit = "inv_cauchit",
    tan_half = "inv_tan_half",
    log1p = "expm1",
    softplus = "log1p_exp",
    squareplus = "squareplus",
    softit = "inv_softit"
  )
  out
}

# define a vector in Stan language
stan_vector <- function(...) {
  paste0("transpose([", paste0(c(...), collapse = ", "), "])")
}

# prepare Stan code for correlations in the generated quantities block
# @param cor name of the correlation vector
# @param ncol number of columns of the correlation matrix
stan_cor_gen_comp <- function(cor, ncol) {
  Cor <- paste0(toupper(substring(cor, 1, 1)), substring(cor, 2))
  glue(
    "  // extract upper diagonal of correlation matrix\n",
    "  for (k in 1:{ncol}) {{\n",
    "    for (j in 1:(k - 1)) {{\n",
    "      {cor}[choose(k - 1, 2) + j] = {Cor}[j, k];\n",
    "    }}\n",
    "  }}\n"
  )
}

# indicates if a family-link combination has a built in
# function in Stan (such as binomial_logit)
# @param family a list with elements 'family' and 'link'
#   ideally a (brms)family object
# @param bterms brmsterms object of the univariate model
stan_has_built_in_fun <- function(family, bterms) {
  stopifnot(all(c("family", "link") %in% names(family)))
  stopifnot(is.brmsterms(bterms))
  cens_or_trunc <- stan_log_lik_adj(bterms$adforms, c("cens", "trunc"))
  link <- family[["link"]]
  dpar <- family[["dpar"]]
  if (cens_or_trunc) {
    # only few families have special lcdf and lccdf functions
    out <- has_built_in_fun(family, link, cdf = TRUE) ||
      has_built_in_fun(bterms, link, dpar = dpar, cdf = TRUE)
  } else {
    out <- has_built_in_fun(family, link) ||
      has_built_in_fun(bterms, link, dpar = dpar)
  }
  out
}

# get all variable names accepted in Stan
stan_all_vars <- function(x) {
  x <- gsub("\\.", "+", x)
  all_vars(x)
}

# transform names to be used as variable names in Stan
make_stan_names <- function(x) {
  gsub("\\.|_", "", make.names(x, unique = TRUE))
}

# functions to handle indexing when threading
stan_slice <- function(threads) {
  str_if(use_threading(threads), "[start:end]")
}

stan_nn <- function(threads) {
  str_if(use_threading(threads), "[nn]", "[n]")
}

stan_nn_def <- function(threads) {
  str_if(use_threading(threads), "    int nn = n + start - 1;\n")
}

stan_nn_regex <- function() {
  "\\[((n)|(nn))\\]"
}

# clean up arguments for partial_log_lik
# @param ... strings containing arguments of the form ', type identifier'
# @return named list of two elements:
#   typed: types + identifiers for use in the function header
#   plain: identifiers only for use in the function call
stan_clean_pll_args <- function(...) {
  args <- paste0(...)
  # split up header to remove duplicates
  typed <- unlist(strsplit(args, ", +"))[-1]
  typed <- unique(typed)
  plain <- rm_wsp(get_matches(" [^ ]+$", typed))
  typed <- collapse(", ", typed)
  plain <- collapse(", ", plain)
  nlist(typed, plain)
}

# prepare a string to be used as comment in Stan
stan_comment <- function(comment, wsp = 2) {
  comment <- as.character(comment)
  wsp <- wsp(nsp = wsp)
  if (!length(comment)) {
    return(character(0))
  }
  ifelse(nzchar(comment), paste0(wsp, "// ", comment), "")
}
