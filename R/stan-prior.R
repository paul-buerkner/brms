# unless otherwise specified, functions return a single character 
# string defining the likelihood of the model in Stan language

# Define priors for parameters in Stan language
# @param prior an object of class 'brmsprior'
# @param class the parameter class
# @param coef the coefficients of this class
# @param group the name of a grouping factor
# @param nlpar the name of a non-linear parameter
# @param prefix a prefix to put at the parameter class
# @param suffix a suffix to put at the parameter class
# @param matrix logical; corresponds the class to a parameter matrix?
# @param wsp a non-negative integer defining the number of spaces 
#   at the start of the output string
# @return
#   A character strings in stan language that defines priors 
#   for a given class of parameters. If a parameter has has 
#   no corresponding prior in prior, an empty string is returned.
stan_prior <- function(prior, class, coef = "", group = "", 
                       px = list(), prefix = "", suffix = "", 
                       wsp = 2, matrix = FALSE) {
  tp <- tp(wsp)
  wsp <- wsp(nsp = wsp)
  prior_only <- identical(attr(prior, "sample_prior"), "only")
  prior <- subset2(prior, 
    class = class, coef = c(coef, ""), group = c(group, "")
  )
  if (class %in% c("sd", "cor")) {
    # only sd and cor parameters have global priors
    px_tmp <- lapply(px, function(x) c(x, ""))
    prior <- subset2(prior, ls = px_tmp)
  } else {
    prior <- subset2(prior, ls = px)
  }
  if (!nchar(class) && nrow(prior)) {
    # unchecked prior statements are directly passed to Stan
    return(collapse(wsp, prior$prior, ";\n"))
  }
  # special priors cannot be passed literally to Stan
  is_special_prior <- is_special_prior(prior$prior)
  if (any(is_special_prior)) {
    special_prior <- prior$prior[is_special_prior]
    stop2("Prior ", collapse_comma(special_prior), " is used in an invalid ", 
          "context. See ?set_prior for details on how to use special priors.")
  }
  
  px <- as.data.frame(px)
  upx <- unique(px)
  if (nrow(upx) > 1L) {
    # can only happen for SD parameters of the same ID
    base_prior <- rep(NA, nrow(upx))
    for (i in seq_rows(upx)) {
      sub_upx <- lapply(upx[i, ], function(x) c(x, ""))
      sub_prior <- subset2(prior, ls = sub_upx) 
      base_prior[i] <- stan_base_prior(sub_prior)
    }
    if (length(unique(base_prior)) > 1L) {
      # define prior for single coefficients manually
      # as there is not single base_prior anymore
      prior_of_coefs <- prior[nzchar(prior$coef), vars_prefix()]
      take <- match_rows(prior_of_coefs, upx)
      prior[nzchar(prior$coef), "prior"] <- base_prior[take]
    }
    base_prior <- base_prior[1]
    bound <- ""
  } else {
    base_prior <- stan_base_prior(prior)
    bound <- prior[!nzchar(prior$coef), "bound"]
  }
  
  individual_prior <- function(i, prior, max_index) {
    # individual priors for each parameter of a class
    index <- ""
    if (max_index > 1L || matrix) {
      index <- glue("[{i}]")      
    }
    if (nrow(px) > 1L) {
      prior <- subset2(prior, ls = px[i, ])
    }
    uc_prior <- prior$prior[match(coef[i], prior$coef)]
    if (!is.na(uc_prior) & nchar(uc_prior)) { 
      # user defined prior for this parameter
      coef_prior <- uc_prior
    } else { 
      # base prior for this parameter
      coef_prior <- base_prior 
    }  
    if (nzchar(coef_prior)) {  
      # implies a proper prior
      pars <- paste0(class, index)
      out <- stan_target_prior(coef_prior, pars, bound = bound)
      out <- paste0(tp, out, ";\n")
    } else {
      # implies an improper flat prior
      out <- ""
    }
    return(out)
  }
  
  # generate stan prior statements
  class <- paste0(prefix, class, suffix)
  if (any(with(prior, nzchar(coef) & nzchar(prior)))) {
    # generate a prior for each coefficient
    out <- sapply(
      seq_along(coef), individual_prior, 
      prior = prior, max_index = length(coef)
    )
  } else if (nchar(base_prior) > 0) {
    if (matrix) {
      class <- glue("to_vector({class})")
    }
    out <- stan_target_prior(
      base_prior, class, ncoef = length(coef), bound = bound
    )
    out <- paste0(tp, out, ";\n")
  } else {
    out <- ""
  }
  out <- collapse(out)
  if (prior_only && nzchar(class) && !nchar(out)) {
    stop2("Sampling from priors is not possible as ", 
          "some parameters have no proper priors. ",
          "Error occured for class '", class, "'.")
  }
  out
}

# get the base prior for all coefficients
# this is the lowest level non-coefficient prior
# @param prior a brmsprior object
# @return a character string defining the base prior
stan_base_prior <- function(prior) {
  stopifnot(length(unique(prior$class)) <= 1)
  take <- with(prior, !nzchar(coef) & nzchar(prior))
  prior <- prior[take, ]
  if (!NROW(prior)) {
    return("")
  }
  vars <- c("group", "nlpar", "dpar", "resp", "class")
  for (v in vars) {
    take <- nzchar(prior[[v]])
    if (any(take)) {
      prior <- prior[take, ]
    }
  }
  stopifnot(NROW(prior) == 1)
  prior$prior
}

# Stan prior in target += notation
# @param prior character string defining the prior
# @param par name of the parameter on which to set the prior
# @param ncoef number of coefficients in the parameter
# @param bound bounds of the parameter in Stan language
# @return a character string defining the prior in Stan language
stan_target_prior <- function(prior, par, ncoef = 1, bound = "") {
  prior <- gsub("[[:space:]]+\\(", "(", prior)
  prior_name <- get_matches(
    "^[^\\(]+(?=\\()", prior, perl = TRUE, simplify = FALSE
  )
  for (i in seq_along(prior_name)) {
    if (length(prior_name[[i]]) != 1L) {
      stop2("The prior '", prior[i], "' is invalid.")
    }
  }
  prior_name <- unlist(prior_name)
  prior_args <- rep(NA, length(prior))
  for (i in seq_along(prior)) {
    prior_args[i] <- sub(glue("^{prior_name[i]}\\("), "", prior[i])
  }
  out <- glue("{prior_name}_lpdf({par} | {prior_args}")
  par_class <- unique(get_matches("^[^_]+", par))
  par_bound <- par_bounds(par_class, bound)
  prior_bound <- prior_bounds(prior_name)
  trunc_lb <- is.character(par_bound$lb) || par_bound$lb > prior_bound$lb
  trunc_ub <- is.character(par_bound$ub) || par_bound$ub < prior_bound$ub
  if (trunc_lb || trunc_ub) {
    wsp <- wsp(nsp = 4)
    if (trunc_lb && !trunc_ub) {
      str_add(out) <- glue(
        "\n{wsp}- {ncoef} * {prior_name}_lccdf({par_bound$lb} | {prior_args}"
      )
    } else if (!trunc_lb && trunc_ub) {
      str_add(out) <- glue(
        "\n{wsp}- {ncoef} * {prior_name}_lcdf({par_bound$ub} | {prior_args}"
      )
    } else if (trunc_lb && trunc_ub) {
      str_add(out) <- glue(
        "\n{wsp}- {ncoef} * log_diff_exp(", 
        "{prior_name}_lcdf({par_bound$ub} | {prior_args}, ",
        "{prior_name}_lcdf({par_bound$lb} | {prior_args})"
      )
    }
  }
  out
}

# Stan code for global parameters of special priors
# currently implemented are horseshoe and lasso
stan_special_prior_global <- function(bterms, data, prior, ...) {
  out <- list()
  tp <- tp()
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  prefix <- combine_prefix(px, keep_mu = TRUE)
  special <- attributes(prior)$special[[prefix]]
  if (!is.null(special[["hs_df"]])) {
    str_add(out$data) <- glue(
      "  // data for the horseshoe prior\n",
      "  real<lower=0> hs_df{p};  // local degrees of freedom\n",
      "  real<lower=0> hs_df_global{p};  // global degrees of freedom\n",
      "  real<lower=0> hs_df_slab{p};  // slab degrees of freedom\n",
      "  real<lower=0> hs_scale_global{p};  // global prior scale\n",
      "  real<lower=0> hs_scale_slab{p};  // slab prior scale\n"           
    )
    str_add(out$par) <- glue(
      "  // horseshoe shrinkage parameters\n",
      "  real<lower=0> hs_global{p}[2];  // global shrinkage parameters\n",
      "  real<lower=0> hs_c2{p};  // slab regularization parameter\n"
    )
    global_args <- glue("0.5 * hs_df_global{p}")
    global_args <- sargs(global_args, global_args)
    c2_args <- glue("0.5 * hs_df_slab{p}")
    c2_args <- sargs(c2_args, c2_args)
    str_add(out$prior) <- glue(
      "{tp}normal_lpdf(hs_global{p}[1] | 0, 1)\n",
      "    - 1 * log(0.5);\n",
      "{tp}inv_gamma_lpdf(hs_global{p}[2] | {global_args});\n",
      "{tp}inv_gamma_lpdf(hs_c2{p} | {c2_args});\n"
    )
  }
  if (!is.null(special[["lasso_df"]])) {
    str_add(out$data) <- glue(
      "  // data for the lasso prior\n",
      "  real<lower=0> lasso_df{p};  // prior degrees of freedom\n",
      "  real<lower=0> lasso_scale{p};  // prior scale\n"
    )
    str_add(out$par) <- glue(
      "  // lasso shrinkage parameter\n",
      "  real<lower=0> lasso_inv_lambda{p};\n"
    )
    str_add(out$prior) <- glue(
      "{tp}chi_square_lpdf(lasso_inv_lambda{p} | lasso_df{p});\n"
    )
  }
  out
}

# Stan code for local parameters of special priors
# currently implemented are 'horseshoe'
# @param class name of the parameter class
# @param prior a brmsprior object
# @param ncoef number of coefficients in the parameter
# @param px named list to subset 'prior'
# @param center_X is the design matrix centered?
# @param suffix optional suffix of the 'b' coefficient vector
stan_special_prior_local <- function(prior, class, ncoef, px, 
                                     center_X = FALSE, suffix = "") {
  class <- as_one_character(class)
  stopifnot(class %in% c("b", "bsp"))
  out <- list()
  p <- usc(combine_prefix(px))
  sp <- paste0(sub("^b", "", class), p)
  ct <- str_if(center_X, "c")
  tp <- tp()
  prefix <- combine_prefix(px, keep_mu = TRUE)
  special <- attributes(prior)$special[[prefix]]
  if (!is.null(special[["hs_df"]])) {
    str_add(out$par) <- glue(
      "  // local parameters for horseshoe prior\n",
      "  vector[K{ct}{sp}] zb{sp};\n",
      "  vector<lower=0>[K{ct}{sp}] hs_local{sp}[2];\n"
    )
    hs_scale_global <- glue("hs_scale_global{p}")
    if (isTRUE(special[["hs_autoscale"]])) {
      str_add(hs_scale_global) <- glue(" * sigma{usc(px$resp)}")
    }
    hs_args <- sargs(
      glue("zb{sp}"), glue("hs_local{sp}"), glue("hs_global{p}"), 
      hs_scale_global, glue("hs_scale_slab{p}^2 * hs_c2{p}")
    )
    str_add(out$tpar_comp) <- glue(
      "  // compute actual regression coefficients\n",
      "  b{sp}{suffix} = horseshoe({hs_args});\n"
    )
    local_args <- glue("0.5 * hs_df{p}")
    local_args <- sargs(local_args, local_args)
    str_add(out$prior) <- glue(
      "{tp}normal_lpdf(zb{sp} | 0, 1);\n",
      "{tp}normal_lpdf(hs_local{sp}[1] | 0, 1)\n", 
      "    - {ncoef} * log(0.5);\n",
      "{tp}inv_gamma_lpdf(hs_local{sp}[2] | {local_args});\n"
    )
  }
  out
}

# Stan code to sample separately from priors
# @param sample_prior take samples from priors?
# @param prior character string taken from stan_prior
# @param par_declars the parameters block of the Stan code
#     required to extract boundaries
# @param gen_quantities Stan code from the generated quantities block
# @param prior_special a list of values pertaining to special priors
#   such as horseshoe or lasso
stan_rngprior <- function(sample_prior, prior, par_declars,
                          gen_quantities, prior_special) {
  if (!sample_prior %in% "yes") {
    return(list())
  }
  prior <- strsplit(gsub(" |\\n", "", prior), ";")[[1]]
  # D will contain all relevant information about the priors
  D <- data.frame(prior = prior[nzchar(prior)])
  pars_regex <- "(?<=_lpdf\\()[^|]+" 
  D$par <- get_matches(pars_regex, D$prior, perl = TRUE, first = TRUE)
  # 'to_vector' should be removed from the parameter names
  has_tv <- grepl("^to_vector\\(", D$par)
  D$par[has_tv] <- gsub("^to_vector\\(|\\)$", "", D$par[has_tv])
  # do not sample from some auxiliary parameters
  excl_regex <- c("z", "zs", "zb", "zgp", "Xn", "Y", "hs", "temp")
  excl_regex <- paste0("(", excl_regex, ")", collapse = "|")
  excl_regex <- paste0("^(", excl_regex, ")(_|$)")
  D <- D[!grepl(excl_regex, D$par), ]
  if (!NROW(D)) return(list())
  
  # rename parameters containing indices
  has_ind <- grepl("\\[[[:digit:]]+\\]", D$par)
  D$par[has_ind] <- ulapply(D$par[has_ind], function(par) {
    ind_regex <- "(?<=\\[)[[:digit:]]+(?=\\])"
    ind <- get_matches(ind_regex, par, perl = TRUE)
    gsub("\\[[[:digit:]]+\\]", paste0("_", ind), par)
  })
  # cannot handle priors on variable transformations
  D <- D[D$par %in% stan_all_vars(D$par), ]
  if (!NROW(D)) return(list())
  
  class_old <- c("^L_", "^Lrescor")
  class_new <- c("cor_", "rescor")
  D$par <- rename(D$par, class_old, class_new, fixed = FALSE)
  dis_regex <- "(?<=target\\+=)[^\\(]+(?=_lpdf\\()"
  D$dist <- get_matches(dis_regex, D$prior, perl = TRUE, first = TRUE)
  D$dist <- sub("corr_cholesky$", "corr", D$dist)
  args_regex <- "(?<=\\|)[^$\\|]+(?=\\)($|-))"
  D$args <- get_matches(args_regex, D$prior, perl = TRUE, first = TRUE)
  
  # extract information from the initial parameter definition
  par_declars <- unlist(strsplit(par_declars, "\n", fixed = TRUE))
  par_declars <- gsub("^[[:blank:]]*", "", par_declars)
  par_declars <- par_declars[!grepl("^//", par_declars)]
  all_pars_regex <- "(?<= )[^[:blank:]]+(?=;)"
  all_pars <- get_matches(all_pars_regex, par_declars, perl = TRUE)
  all_pars <- rename(all_pars, class_old, class_new, fixed = FALSE)
  all_bounds <- get_matches("<.+>", par_declars, first = TRUE)
  all_types <- get_matches("^[^[:blank:]]+", par_declars)
  all_dims <- get_matches(
    "(?<=\\[)[^\\]]*", par_declars, first = TRUE, perl = TRUE
  )
  
  # define parameter types and boundaries
  D$dim <- D$bounds <- ""
  D$type <- "real"
  for (i in seq_along(all_pars)) {
    k <- which(grepl(paste0("^", all_pars[i]), D$par))
    D$dim[k] <- all_dims[i]
    D$bounds[k] <- all_bounds[i]
    if (grepl("^(simo_)|(theta)", all_pars[i])) {
      D$type[k] <- all_types[i]
    }
  }

  # exclude priors which depend on other priors
  # TODO: enable sampling from these priors as well
  found_vars <- lapply(D$args, find_vars, dot = FALSE, brackets = FALSE)
  contains_other_pars <- ulapply(found_vars, function(x) any(x %in% all_pars))
  D <- D[!contains_other_pars, ]
  if (!NROW(D)) return(list())
  
  out <- list()
  # sample priors in the generated quantities block
  D$lkj <- grepl("^lkj_corr$", D$dist)
  D$args <- paste0(ifelse(D$lkj, paste0(D$dim, ","), ""), D$args)
  D$lkj_index <- ifelse(D$lkj, "[1, 2]", "")
  D$prior_par <- glue("prior_{D$par}")
  str_add(out$gen_def) <- "  // additionally draw samples from priors\n"
  str_add(out$gen_def) <- cglue(
    "  {D$type} {D$prior_par} = {D$dist}_rng({D$args}){D$lkj_index};\n"
  )
  
  # sample from truncated priors using rejection sampling
  D$lb <- stan_extract_bounds(D$bounds, bound = "lower")
  D$ub <- stan_extract_bounds(D$bounds, bound = "upper")
  Ibounds <- which(nzchar(D$bounds))
  if (length(Ibounds)) {
    str_add(out$gen_comp) <- "  // use rejection sampling for truncated priors\n"
    for (i in Ibounds) {
      wl <- if (nzchar(D$lb[i])) glue("{D$prior_par[i]} < {D$lb[i]}")
      wu <- if (nzchar(D$ub[i])) glue("{D$prior_par[i]} > {D$ub[i]}")
      prior_while <- paste0(c(wl, wu), collapse = " || ")
      str_add(out$gen_comp) <- glue(
        "  while ({prior_while}) {{\n",
        "    {D$prior_par[i]} = {D$dist[i]}_rng({D$args[i]}){D$lkj_index[i]};\n",
        "  }}\n"
      )
    }
  }
  out
}

# indicate if the horseshoe prior is used in the predictor term
stan_use_horseshoe <- function(bterms, prior) {
  prefix <- combine_prefix(bterms, keep_mu = TRUE)
  special <- attr(prior, "special")[[prefix]]
  !is.null(special[["hs_df"]])
}

# extract Stan boundaries expression from a string
stan_extract_bounds <- function(x, bound = c("lower", "upper")) {
  bound <- match.arg(bound)
  x <- rm_wsp(x)
  regex <- glue("(?<={bound}=)[^,>]*")
  get_matches(regex, x, perl = TRUE, first = TRUE)
}
