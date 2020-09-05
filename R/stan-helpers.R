# unless otherwise specifiedm functions return a named list 
# of Stan code snippets to be pasted together later on

# define Stan functions or globally used transformed data
stan_global_defs <- function(bterms, prior, ranef) {
  families <- family_names(bterms)
  links <- family_info(bterms, "link")
  unique_combs <- !duplicated(paste0(families, ":", links))
  families <- families[unique_combs]
  links <- links[unique_combs]
  out <- list()
  if (any(links == "cauchit")) {
    str_add(out$fun) <- "  #include 'fun_cauchit.stan'\n"
  } else if (any(links == "cloglog")) {
    str_add(out$fun) <- "  #include 'fun_cloglog.stan'\n"
  } else if (any(links == "softplus")) {
    str_add(out$fun) <- "  #include 'fun_softplus.stan'\n"
  }
  hs_dfs <- ulapply(attr(prior, "special"), "[[", "hs_df")
  if (any(nzchar(hs_dfs))) {
    str_add(out$fun) <- "  #include 'fun_horseshoe.stan'\n"
  }
  if (any(nzchar(ranef$by))) {
    str_add(out$fun) <- "  #include 'fun_scale_r_cor_by.stan'\n"
  }
  if (stan_needs_kronecker(ranef)) {
    str_add(out$fun) <- glue(
      "  #include 'fun_as_matrix.stan'\n",
      "  #include 'fun_kronecker.stan'\n"
    )
  }
  family_files <- family_info(bterms, "include")
  if (length(family_files)) {
    str_add(out$fun) <- cglue("  #include '{family_files}'\n")
  }
  const <- family_info(bterms, "const")
  if (length(const)) {
    str_add(out$tdata_def) <- cglue("  {const};\n")
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
  }
  acterms <- get_effect(bterms, "ac")
  acefs <- lapply(acterms, tidy_acef)
  if (any(ulapply(acefs, has_ac_subset, dim = "time", cov = TRUE))) {
    # TODO: include functions selectively
    str_add(out$fun) <- glue(
      "  #include 'fun_normal_time.stan'\n",
      "  #include 'fun_student_t_time.stan'\n",
      "  #include 'fun_scale_time_err.stan'\n",
      "  #include 'fun_cholesky_cor_ar1.stan'\n",
      "  #include 'fun_cholesky_cor_ma1.stan'\n",
      "  #include 'fun_cholesky_cor_arma1.stan'\n",
      "  #include 'fun_cholesky_cor_cosy.stan'\n"
    )
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
    softplus = "log_expm1"
  )
}

# inverse link in Stan language
# @param link name of the link function
stan_ilink <- function(link) {
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
    softplus = "log1p_exp"
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
# @param cens_or_trunc is the model censored or truncated?
stan_has_built_in_fun <- function(family, cens_or_trunc = FALSE) {
  stopifnot(all(c("family", "link") %in% names(family)))
  link <- family$link
  dpar <- family$dpar
  family <- family$family
  if (cens_or_trunc) {
    # only few families have special lcdf and lccdf functions
    log_families <- c("cox")
    logit_families <- character(0)
    logit_dpars <- character(0)
  } else {
    log_families <- c(
      "poisson", "negbinomial", "geometric", "com_poisson",
      "zero_inflated_poisson", "zero_inflated_negbinomial",
      "hurdle_poisson", "hurdle_negbinomial", "cox"
    )
    logit_families <- c(
      "binomial", "bernoulli", "cumulative", "categorical",
      "zero_inflated_binomial"
    )
    logit_dpars <- c("zi", "hu")
  }
  isTRUE(
    family %in% log_families && link == "log" ||
    family %in% logit_families && link == "logit" ||
    isTRUE(dpar %in% logit_dpars) && link == "logit"
  )
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

# checks if a model needs the kronecker product
# @param ranef output of tidy_ranef
# @return a single logical value
stan_needs_kronecker <- function(ranef) {
  ids <- unique(ranef$id)
  out <- FALSE
  for (id in ids) {
    r <- ranef[ranef$id == id, ]
    out <- out || nrow(r) > 1L && r$cor[1] && nzchar(r$cov[1])
  }
  out
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
