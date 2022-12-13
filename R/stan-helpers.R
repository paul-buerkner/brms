# unless otherwise specified, functions return a named list
# of Stan code snippets to be pasted together later on

# define Stan functions or globally used transformed data
# TODO: refactor to not require extraction of information from all model parts
#   'expand_include_statements' removes duplicates which opens the door
#   for adding Stan functions at better places rather than globally here
stan_global_defs <- function(bterms, prior, ranef, threads) {
  families <- family_names(bterms)
  links <- family_info(bterms, "link")
  unique_combs <- !duplicated(paste0(families, ":", links))
  families <- families[unique_combs]
  links <- links[unique_combs]
  out <- list()
  # TODO: detect these links in all dpars not just in 'mu'
  if (any(links == "cauchit")) {
    str_add(out$fun) <- "  #include 'fun_cauchit.stan'\n"
  } else if (any(links == "cloglog")) {
    str_add(out$fun) <- "  #include 'fun_cloglog.stan'\n"
  } else if (any(links == "softplus")) {
    str_add(out$fun) <- "  #include 'fun_softplus.stan'\n"
  } else if (any(links == "squareplus")) {
    str_add(out$fun) <- "  #include 'fun_squareplus.stan'\n"
  } else if (any(links == "softit")) {
    str_add(out$fun) <- "  #include 'fun_softit.stan'\n"
  }
  special <- get_special_prior(prior)
  if (!isNULL(lapply(special, "[[", "horseshoe"))) {
    str_add(out$fun) <- "  #include 'fun_horseshoe.stan'\n"
  }
  if (!isNULL(lapply(special, "[[", "R2D2"))) {
    str_add(out$fun) <- "  #include 'fun_r2d2.stan'\n"
  }
  if (nrow(ranef)) {
    r_funs <- NULL
    ids <- unique(ranef$id)
    for (id in ids) {
      r <- ranef[ranef$id == id, ]
      if (nrow(r) > 1L && r$cor[1]) {
        if (nzchar(r$by[1])) {
          if (nzchar(r$cov[1])) {
            c(r_funs) <- "  #include 'fun_scale_r_cor_by_cov.stan'\n"
          } else {
            c(r_funs) <- "  #include 'fun_scale_r_cor_by.stan'\n"
          }
        } else {
          if (nzchar(r$cov[1])) {
            c(r_funs) <- "  #include 'fun_scale_r_cor_cov.stan'\n"
          } else {
            c(r_funs) <- "  #include 'fun_scale_r_cor.stan'\n"
          }
        }
      }
    }
    str_add(out$fun) <- collapse(unique(r_funs))
  }
  family_files <- family_info(bterms, "include")
  if (length(family_files)) {
    str_add(out$fun) <- cglue("  #include '{family_files}'\n")
  }
  is_ordinal <- ulapply(families, is_ordinal)
  if (any(is_ordinal)) {
    ord_fams <- families[is_ordinal]
    ord_links <- links[is_ordinal]
    for (i in seq_along(ord_fams)) {
      str_add(out$fun) <- stan_ordinal_lpmf(ord_fams[i], ord_links[i])
    }
  }
  uni_mo <- ulapply(get_effect(bterms, "sp"), attr, "uni_mo")
  if (length(uni_mo)) {
    str_add(out$fun) <- "  #include 'fun_monotonic.stan'\n"
  }
  if (length(get_effect(bterms, "gp"))) {
    # TODO: include functions selectively
    str_add(out$fun) <- "  #include 'fun_gaussian_process.stan'\n"
    str_add(out$fun) <- "  #include 'fun_gaussian_process_approx.stan'\n"
    str_add(out$fun) <- "  #include 'fun_which_range.stan'\n"
  }
  acterms <- get_effect(bterms, "ac")
  acefs <- lapply(acterms, tidy_acef)
  if (any(ulapply(acefs, has_ac_subset, dim = "time", cov = TRUE))) {
    str_add(out$fun) <- glue(
      "  #include 'fun_is_equal.stan'\n"
    )
    if ("gaussian" %in% families) {
      str_add(out$fun) <- glue(
        "  #include 'fun_normal_time.stan'\n",
        "  #include 'fun_normal_time_se.stan'\n"
      )
    }
    if ("student" %in% families) {
      str_add(out$fun) <- glue(
        "  #include 'fun_student_t_time.stan'\n",
        "  #include 'fun_student_t_time_se.stan'\n"
      )
    }
    # TODO: include selectively once we have the 'latent' indicator
    str_add(out$fun) <- glue(
      "  #include 'fun_scale_time_err.stan'\n"
    )
    if (any(ulapply(acefs, has_ac_class, "arma"))) {
      str_add(out$fun) <- glue(
        "  #include 'fun_cholesky_cor_ar1.stan'\n",
        "  #include 'fun_cholesky_cor_ma1.stan'\n",
        "  #include 'fun_cholesky_cor_arma1.stan'\n"
      )
    }
    if (any(ulapply(acefs, has_ac_class, "cosy"))) {
      str_add(out$fun) <- glue(
        "  #include 'fun_cholesky_cor_cosy.stan'\n"
      )
    }
  }
  if (any(ulapply(acefs, has_ac_class, "sar"))) {
    if ("gaussian" %in% families) {
      str_add(out$fun) <- glue(
        "  #include 'fun_normal_lagsar.stan'\n",
        "  #include 'fun_normal_errorsar.stan'\n"
      )
    }
    if ("student" %in% families) {
      str_add(out$fun) <- glue(
        "  #include 'fun_student_t_lagsar.stan'\n",
        "  #include 'fun_student_t_errorsar.stan'\n"
      )
    }
  }
  if (any(ulapply(acefs, has_ac_class, "car"))) {
    str_add(out$fun) <- glue(
      "  #include 'fun_sparse_car_lpdf.stan'\n",
      "  #include 'fun_sparse_icar_lpdf.stan'\n"
    )
  }
  if (any(ulapply(acefs, has_ac_class, "fcor"))) {
    str_add(out$fun) <- glue(
      "  #include 'fun_normal_fcor.stan'\n",
      "  #include 'fun_student_t_fcor.stan'\n"
    )
  }
  if (use_threading(threads)) {
    str_add(out$fun) <- "  #include 'fun_sequence.stan'\n"
  }
  out
}

# link function in Stan language
# @param link name of the link function
stan_link <- function(link) {
  switch(link,
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
}

# inverse link in Stan language
# @param link name of the link function
stan_inv_link <- function(link) {
  switch(link,
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
