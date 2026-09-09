# unless otherwise specified, functions return a named list
# of Stan code snippets to be pasted together later on

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
# @param bterms brmsterms object of the univariate model
# @param family a list with elements 'family' and 'link'
#   ideally a (brms)family object
stan_has_built_in_fun <- function(bterms, family = NULL) {
  stopifnot(is.brmsterms(bterms))
  family <- family %||% bterms$family
  stopifnot(all(c("family", "link") %in% names(family)))
  link <- family[["link"]]
  dpar <- family[["dpar"]]
  # only few families have special lcdf and lccdf functions
  cdf <- has_ad_terms(bterms, c("cens", "trunc"))
  has_built_in_fun(family, link, cdf = cdf) ||
    has_built_in_fun(bterms, link, dpar = dpar, cdf = cdf)
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

stan_nn <- function(threads, subsample = NULL) {
  if (use_threading(threads)) return("[nn]")
  if (use_subsampling(subsample)) return("[nn]")
  "[n]"
}

stan_nn_def <- function(threads, subsample = NULL) {
  if (use_threading(threads)) return("    int nn = n + start - 1;\n")
  if (use_subsampling(subsample)) {
    return(glue("    int nn = {subsample$index_fn}(n);\n"))
  }
  ""
}

stan_nn_regex <- function() {
  "\\[((n)|(nn))\\]"
}

# expression for the number of observations, subsample-aware
# @param resp response suffix (possibly empty)
# @param subsample a brmssubsample object or NULL
stan_N_expr <- function(resp = "", subsample = NULL) {
  if (use_subsampling(subsample)) {
    glue("{subsample$size_fn}()")
  } else {
    glue("N{resp}")
  }
}

# wrap a Stan variable with a subsample getter function
# @param var_expr Stan variable expression (e.g. "Xc", "Y")
# @param subsample a brmssubsample object or NULL
# @return wrapped expression if a matching wrap entry exists, else unchanged
stan_subsample_wrap <- function(var_expr, subsample = NULL) {
  if (!use_subsampling(subsample) || length(subsample$wrap) == 0) {
    return(var_expr)
  }
  for (var_name in names(subsample$wrap)) {
    if (grepl(paste0("^", var_name, "(\\b|$)"), var_expr)) {
      return(glue("{subsample$wrap[[var_name]]}({var_expr})"))
    }
  }
  var_expr
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
