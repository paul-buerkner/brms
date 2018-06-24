stan_prior <- function(prior, class, coef = "", group = "", 
                       px = list(), prefix = "", suffix = "", 
                       wsp = 2, matrix = FALSE) {
  # Define priors for parameters in Stan language
  # Args:
  #   prior: an object of class 'brmsprior'
  #   class: the parameter class
  #   coef: the coefficients of this class
  #   group: the name of a grouping factor
  #   nlpar: the name of a non-linear parameter
  #   prefix: a prefix to put at the parameter class
  #   suffix: a suffix to put at the parameter class
  #   matrix: logical; corresponds the class to a parameter matrix?
  #   wsp: an integer >= 0 defining the number of spaces 
  #        in front of the output string
  # Returns:
  #   A character strings in stan language that defines priors 
  #   for a given class of parameters. If a parameter has has 
  #   no corresponding prior in prior, an empty string is returned.
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
    return(collapse(wsp, prior$prior, "; \n"))
  } 
  
  px <- as.data.frame(px)
  upx <- unique(px)
  if (nrow(upx) > 1L) {
    # can only happen for SD parameters of the same ID
    base_prior <- rep(NA, nrow(upx))
    for (i in seq_len(nrow(upx))) {
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
    if (max_index > 1L || matrix) {
      index <- paste0("[", i, "]")      
    } else {
      index <- ""
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
      out <- paste0(tp, out, "; \n")
    } else {
      # implies an improper flat prior
      out <- ""
    }
    return(out)
  }
  
  # generate stan prior statements
  class <- paste0(prefix, class, suffix)
  if (any(with(prior, nchar(coef) & nchar(prior)))) {
    # generate a prior for each coefficient
    out <- sapply(
      seq_along(coef), individual_prior, 
      prior = prior, max_index = length(coef)
    )
  } else if (nchar(base_prior) > 0) {
    if (matrix) {
      class <- paste0("to_vector(", class, ")")
    }
    out <- stan_target_prior(
      base_prior, class, ncoef = length(coef), bound = bound
    )
    out <- paste0(tp, out, "; \n")
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

stan_base_prior <- function(prior) {
  # get base (highest level) prior of all priors
  # Args:
  #   prior: a prior.frame
  stopifnot(length(unique(prior$class)) <= 1L) 
  prior <- prior[with(prior, !nzchar(coef) & nzchar(prior)), ]
  vars <- c("group", "nlpar", "dpar", "resp", "class")
  i <- 1
  found <- FALSE
  base_prior <- ""
  take <- rep(FALSE, nrow(prior))
  while (!found && i <= length(vars)) {
    take <- nzchar(prior[[vars[i]]]) & !take
    if (any(take)) {
      base_prior <- prior[take, "prior"]
      found <- TRUE
    }
    i <- i + 1
  }
  stopifnot(length(base_prior) == 1L)
  base_prior
}

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
    prior_args[i] <- sub(
      paste0("^", prior_name[i], "\\("), "", prior[i]
    )
  }
  out <- paste0(prior_name, "_lpdf(", par, " | ", prior_args)
  par_class <- unique(get_matches("^[^_]+", par))
  par_bound <- par_bounds(par_class, bound)
  prior_bound <- prior_bounds(prior_name)
  trunc_lb <- is.character(par_bound$lb) || par_bound$lb > prior_bound$lb
  trunc_ub <- is.character(par_bound$ub) || par_bound$ub < prior_bound$ub
  if (trunc_lb || trunc_ub) {
    wsp <- wsp(nsp = 4)
    if (trunc_lb && !trunc_ub) {
      str_add(out) <- paste0(
        "\n", wsp, "- ", ncoef, " * ", prior_name, "_lccdf(", 
        par_bound$lb, " | ", prior_args
      )
    } else if (!trunc_lb && trunc_ub) {
      str_add(out) <- paste0(
        "\n", wsp, "- ", ncoef, " * ", prior_name, "_lcdf(", 
        par_bound$ub, " | ", prior_args
      )
    } else if (trunc_lb && trunc_ub) {
      str_add(out) <- paste0(
        "\n", wsp, "- ", ncoef, " * log_diff_exp(", 
        prior_name, "_lcdf(", par_bound$ub, " | ", prior_args, ", ",
        prior_name, "_lcdf(", par_bound$lb, " | ", prior_args, ")"
      )
    }
  }
  out
}

stan_special_prior_global <- function(bterms, data, prior) {
  # Stan code for global parameters of special priors
  # currently implemented are horseshoe and lasso
  out <- list()
  tp <- tp()
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  prefix <- combine_prefix(px, keep_mu = TRUE)
  special <- attributes(prior)$special[[prefix]]
  if (!is.null(special[["hs_df"]])) {
    str_add(out$data) <- paste0(
      "  real<lower=0> hs_df", p, "; \n",
      "  real<lower=0> hs_df_global", p, "; \n",
      "  real<lower=0> hs_df_slab", p, "; \n",
      "  real<lower=0> hs_scale_global", p, "; \n",
      "  real<lower=0> hs_scale_slab", p, "; \n"           
    )
    str_add(out$par) <- paste0(
      "  // horseshoe shrinkage parameters \n",
      "  real<lower=0> hs_global", p, "[2]; \n",
      "  real<lower=0> hs_c2", p, "; \n"
    )
    global_args <- paste0("0.5 * hs_df_global", p)
    global_args <- sargs(global_args, global_args)
    c2_args <- paste0("0.5 * hs_df_slab", p)
    c2_args <- sargs(c2_args, c2_args)
    str_add(out$prior) <- paste0(
      tp, "normal_lpdf(hs_global", p, "[1] | 0, 1)\n    - 1 * log(0.5);\n",
      tp, "inv_gamma_lpdf(hs_global", p, "[2] | ", global_args, ");\n",
      tp, "inv_gamma_lpdf(hs_c2", p, " | ", c2_args, ");\n"
    )
  }
  if (!is.null(special[["lasso_df"]])) {
    str_add(out$data) <- paste0(
      "  real<lower=0> lasso_df", p, "; \n",
      "  real<lower=0> lasso_scale", p, "; \n"
    )
    str_add(out$par) <- paste0(
      "  // lasso shrinkage parameter \n",
      "  real<lower=0> lasso_inv_lambda", p, "; \n"
    )
    str_add(out$prior) <- paste0(
      tp, "chi_square_lpdf(lasso_inv_lambda", p, " | lasso_df", p, ");\n"
    )
  }
  out
}

stan_special_prior_local <- function(class, prior, ncoef, px,
                                     center_X = FALSE)  {
  # Stan code for local parameters of special priors
  # currently implemented are horseshoe
  class <- as_one_character(class)
  stopifnot(class %in% c("b", "bsp"))
  out <- list()
  p <- usc(combine_prefix(px))
  sp <- paste0(sub("^b", "", class), p)
  ct <- ifelse(center_X, "c", "")
  tp <- tp()
  prefix <- combine_prefix(px, keep_mu = TRUE)
  special <- attributes(prior)$special[[prefix]]
  if (!is.null(special[["hs_df"]])) {
    str_add(out$par) <- paste0(
      "  // local parameters for horseshoe prior\n",
      "  vector[K", ct, sp, "] zb", sp, ";\n",
      "  vector<lower=0>[K", ct, sp, "] hs_local", sp, "[2];\n"
    )
    hs_scale_global <- paste0("hs_scale_global", p)
    if (isTRUE(special[["hs_autoscale"]])) {
      str_add(hs_scale_global) <- paste0(" * sigma", usc(px$resp))
    }
    hs_args <- sargs(
      paste0(c("zb", "hs_local"), sp), paste0("hs_global", p), 
      hs_scale_global, paste0("hs_scale_slab", p, "^2 * hs_c2", p)
    )
    str_add(out$tparD) <- paste0(
      "  vector[K", ct, sp, "] b", sp,
      " = horseshoe(", hs_args, "); \n"
    )
    local_args <- paste0("0.5 * hs_df", p)
    local_args <- sargs(local_args, local_args)
    str_add(out$prior) <- paste0(
      tp, "normal_lpdf(zb", sp, " | 0, 1);\n",
      tp, "normal_lpdf(hs_local", sp, "[1] | 0, 1)\n", 
      "    - ", ncoef, " * log(0.5);\n",
      tp, "inv_gamma_lpdf(hs_local", sp, "[2] | ", local_args, ");\n"
    )
  }
  out
}

stan_rngprior <- function(sample_prior, prior, par_declars,
                          gen_quantities, prior_special) {
  # stan code to sample from priors seperately
  # Args:
  #   sample_prior: take samples from priors?
  #   prior: character string taken from stan_prior
  #   par_declars: the parameters block of the Stan code
  #     required to extract boundaries
  #   gen_quantities: Stan code from the generated quantities block
  #   prior_special: a list of values pertaining to special priors
  #     such as horseshoe or lasso
  # Returns:
  #   a character string containing the priors to be sampled from in stan code
  out <- list()
  if (!sample_prior %in% "yes") {
    return(out)
  }
  prior <- strsplit(gsub(" |\\n", "", prior), ";")[[1]]
  prior <- prior[nzchar(prior)]
  pars_regex <- "(?<=_lpdf\\()[^|]+" 
  pars <- get_matches(pars_regex, prior, perl = TRUE, first = TRUE)
  pars <- gsub("to_vector\\(|\\)$", "", pars)
  excl_regex <- c("z", "zs", "zb", "zgp", "Xn", "Y", "hs")
  excl_regex <- paste0("(", excl_regex, ")", collapse = "|")
  excl_regex <- paste0("^(", excl_regex, ")(_|$)")
  take <- !grepl(excl_regex, pars)
  prior <- prior[take]
  pars <- pars[take]
  pars <- sub("^L_", "cor_", pars)
  pars <- sub("^Lrescor", "rescor", pars)
  dis_regex <- "(?<=target\\+=)[^\\(]+(?=_lpdf\\()"
  dis <- get_matches(dis_regex, prior, perl = TRUE, first = TRUE)
  dis <- sub("corr_cholesky$", "corr", dis)
  args_regex <- "(?<=\\|)[^$\\|]+(?=\\)($|-))"
  args <- get_matches(args_regex, prior, perl = TRUE, first = TRUE)
  args <- ifelse(grepl("^lkj_corr$", dis), paste0("2,", args), args)
  type <- rep("real", length(pars))
  
  # rename parameters containing indices
  has_ind <- grepl("\\[[[:digit:]]+\\]", pars)
  pars[has_ind] <- ulapply(pars[has_ind], function(par) {
    ind_regex <- "(?<=\\[)[[:digit:]]+(?=\\])"
    ind <- get_matches(ind_regex, par, perl = TRUE)
    gsub("\\[[[:digit:]]+\\]", paste0("_", ind), par)
  })
  
  # extract information from the initial parameter definition
  par_declars <- unlist(strsplit(par_declars, "\n", fixed = TRUE))
  par_declars <- gsub("^[[:blank:]]*", "", par_declars)
  par_declars <- par_declars[!grepl("^//", par_declars)]
  all_pars_regex <- "(?<= )[^[:blank:]]+(?=;)"
  all_pars <- get_matches(all_pars_regex, par_declars, perl = TRUE) 
  all_bounds <- get_matches("<.+>", par_declars, simplify = FALSE)
  all_bounds <- ulapply(all_bounds, function(x) if (length(x)) x else "")
  all_types <- get_matches("^[^[:blank:]]+", par_declars)
  
  # define parameter types and boundaries
  bounds <- rep("", length(pars))
  types <- rep("real", length(pars))
  for (i in seq_along(all_pars)) {
    k <- which(grepl(paste0("^", all_pars[i]), pars))
    bounds[k] <- all_bounds[i]
    if (grepl("^simo", all_pars[i])) {
      types[k] <- all_types[i]
    }
  }
  
  # distinguish between bounded and unbounded parameters
  has_bounds <- as.logical(nchar(bounds))
  if (any(has_bounds)) {  
    # bounded parameters have to be sampled in the model block
    str_add(out$par) <- paste0(
      "  // parameters to store prior samples\n",
      collapse(
        "  real", bounds[has_bounds], 
        " prior_", pars[has_bounds], ";\n"
      )
    )
    str_add(out$model) <- paste0(
      "  // additionally draw samples from priors\n",
      collapse(
        "  target += ", dis[has_bounds], "_lpdf(",
        "prior_", pars[has_bounds], " | ", args[has_bounds], ");\n"
      )
    )
  }
  no_bounds <- !has_bounds
  if (any(no_bounds)) {
    # use parameters sampled from priors for use in other priors
    spars <- character(0)
    # cannot sample from the horseshoe prior anymore as of brms 1.5.0
    lasso_prefix <- nzchar(ulapply(prior_special, "[[", "lasso_df"))
    lasso_prefix <- names(prior_special)[lasso_prefix]
    lasso_prefix <- usc(sub("^mu(_|$)", "", lasso_prefix))
    if (length(lasso_prefix)) {
      spars <- c(spars, paste0("lasso_inv_lambda", lasso_prefix))
    }
    if (length(spars)) {
      bpars <- grepl("^(b|(bsp)|(bcs))(_|$)", pars)
      args[bpars] <- rename(args[bpars], spars, paste0("prior_", spars))
    }
    lkj_index <- ifelse(grepl("^lkj_corr$", dis[no_bounds]), "[1, 2]", "")
    # unbounded parameters can be sampled in the generatated quantities block
    str_add(out$genD) <- paste0(
      "  // additionally draw samples from priors\n",
      collapse(
        "  ", types[no_bounds], " prior_", pars[no_bounds], 
        " = ", dis[no_bounds], "_rng(", args[no_bounds], ")",
        lkj_index, ";\n"
      )
    )
  }
  # compute priors for the actual population-level intercepts
  is_temp_intercept <- grepl("^temp.*_Intercept", pars)
  if (any(is_temp_intercept)) {
    temp_intercepts <- pars[is_temp_intercept]
    p <- gsub("^temp|_Intercept$", "", temp_intercepts)
    intercepts <- paste0("b", p, "_Intercept")
    regex <- paste0(" (\\+|-) dot_product\\(means_X", p, ", b", p, "\\)")
    X_means <- lapply(regex, get_matches, gen_quantities)
    X_means <- ulapply(X_means, function(x) if (length(x)) x[1] else "")
    str_add(out$genD) <- collapse(
      "  real prior_", intercepts,
      " = prior_", temp_intercepts, X_means, ";\n"
    )
  }
  out
}
