# unless otherwise specifiedm functions return a named list 
# of Stan code snippets to be pasted together later on

# Stan code for the response variables
stan_response <- function(bterms, data) {
  stopifnot(is.brmsterms(bterms))
  family <- bterms$family
  rtype <- str_if(use_int(family), "int", "real")
  multicol <- has_multicol(family)
  px <- check_prefix(bterms)
  resp <- usc(combine_prefix(px))
  out <- list()
  str_add(out$data) <- glue(
    "  int<lower=1> N{resp};  // number of observations\n"
  )
  
  if (has_cat(family)) {
    str_add(out$data) <- glue(
      "  int<lower=2> ncat{resp};  // number of categories\n"
    )
  }
  if (has_multicol(family)) {
    if (rtype == "real") {
      str_add(out$data) <- glue(
        "  vector[ncat{resp}] Y{resp}[N{resp}];  // response array\n"
      )
    } else if (rtype == "int") {
      str_add(out$data) <- glue(
        "  int Y{resp}[N{resp}, ncat{resp}];  // response array\n"
      )
    }
  } else {
    if (rtype == "real") {
      # don't use real Y[n]
      str_add(out$data) <- glue(
        "  vector[N{resp}] Y{resp};  // response variable\n"
      )
    } else if (rtype == "int") {
      str_add(out$data) <- glue(
        "  int Y{resp}[N{resp}];  // response variable\n"
      )
    }
  }
  if (has_ndt(family)) {
    str_add(out$tdata_def) <- glue(
      "  real min_Y{resp} = min(Y{resp});\n"
    )
  }
  if (has_trials(family) || is.formula(bterms$adforms$trials)) {
    str_add(out$data) <- glue(
      "  int trials{resp}[N{resp}];  // number of trials\n"
    )
  }
  if (is.formula(bterms$adforms$weights)) {
    str_add(out$data) <- glue(
      "  vector<lower=0>[N{resp}] weights{resp};  // model weights\n" 
    )
  }
  if (has_thres(family)) {
    groups <- get_thres_groups(family)
    if (any(nzchar(groups))) {
      str_add(out$data) <- glue(
        "  int<lower=1> ngrthres{resp};  // number of threshold groups\n",
        "  int<lower=1> nthres{resp}[ngrthres{resp}];  // number of thresholds\n",
        "  int<lower=1> Jthres{resp}[N{resp}, 2];  // threshold indices\n"
      )
      str_add(out$tdata_def) <- glue(
        "  int<lower=1> nmthres{resp} = prod(nthres{resp});", 
        "  // total number of thresholds\n",
        "  int<lower=1> Kthres_start{resp}[ngrthres{resp}];", 
        "  // start index per threshold group\n",
        "  int<lower=1> Kthres_end{resp}[ngrthres{resp}];",
        "  // end index per threshold group\n"
      )
      str_add(out$tdata_comp) <- glue(
        "  Kthres_start{resp}[1] = 1;\n",
        "  Kthres_end{resp}[1] = nthres{resp}[1];\n",
        "  for (i in 2:ngrthres{resp}) {{\n",
        "    Kthres_start{resp}[i] = Kthres_end{resp}[i-1] + 1;\n",
        "    Kthres_end{resp}[i] = Kthres_end{resp}[i-1] + nthres{resp}[i];\n",
        "  }}\n"
      )
    } else {
      str_add(out$data) <- glue(
        "  int<lower=2> nthres{resp};  // number of thresholds\n"
      )
    }
  }
  if (is.formula(bterms$adforms$se)) {
    str_add(out$data) <- glue(
      "  vector<lower=0>[N{resp}] se{resp};  // known sampling error\n"
    )
    str_add(out$tdata_def) <- glue(
      "  vector<lower=0>[N{resp}] se2{resp} = square(se{resp});\n"
    )
  }
  if (is.formula(bterms$adforms$dec)) {
    str_add(out$data) <- glue(
      "  int<lower=0,upper=1> dec{resp}[N{resp}];  // decisions\n"
    )
  }
  if (is.formula(bterms$adforms$rate)) {
    str_add(out$data) <- glue(
      "  vector<lower=0>[N{resp}] denom{resp};",
      "  // response denominator\n"
    )
    str_add(out$tdata_def) <- glue(
      "  // log response denominator\n",
      "  vector[N{resp}] log_denom{resp} = log(denom{resp});\n"
    )
  }
  if (is.formula(bterms$adforms$cens)) {
    cens <- eval_rhs(bterms$adforms$cens)
    str_add(out$data) <- glue(
      "  int<lower=-1,upper=2> cens{resp}[N{resp}];  // indicates censoring\n"
    )
    if (cens$vars$y2 != "NA") {
      # interval censoring is required
      rcens <- str_if(rtype == "int", 
                      glue("  int rcens{resp}[N{resp}];"), 
                      glue("  vector[N{resp}] rcens{resp};")
      )
      str_add(out$data) <- glue(
        rcens, "  // right censor points for interval censoring\n"
      )
    }
  }
  bounds <- trunc_bounds(bterms, data = data)
  if (any(bounds$lb > -Inf)) {
    str_add(out$data) <- glue(
      "  {rtype} lb{resp}[N{resp}];  // lower truncation bounds;\n"
    )
  }
  if (any(bounds$ub < Inf)) {
    str_add(out$data) <- glue(
      "  {rtype} ub{resp}[N{resp}];  // upper truncation bounds\n"
    )
  }
  if (is.formula(bterms$adforms$mi)) {
    Ybounds <- trunc_bounds(bterms, data, incl_family = TRUE, stan = TRUE)
    sdy <- get_sdy(bterms, data)
    if (is.null(sdy)) {
      # response is modeled without measurement error
      str_add(out$par) <- glue(
        "  vector{Ybounds}[Nmi{resp}] Ymi{resp};",
        "  // estimated missings\n"
      )
      str_add(out$data) <- glue(
        "  int<lower=0> Nmi{resp};  // number of missings\n",
        "  int<lower=1> Jmi{resp}[Nmi{resp}];",  
        "  // positions of missings\n"
      )
      str_add(out$model_def) <- glue(
        "  // vector combining observed and missing responses\n",
        "  vector[N{resp}] Yl{resp} = Y{resp};\n" 
      )
      str_add(out$model_comp_basic) <- glue(
        "  Yl{resp}[Jmi{resp}] = Ymi{resp};\n"
      )
    } else {
      str_add(out$data) <- glue(
        "  // data for measurement-error in the response\n",
        "  vector<lower=0>[N{resp}] noise{resp};\n",
        "  // information about non-missings\n",
        "  int<lower=0> Nme{resp};\n",
        "  int<lower=1> Jme{resp}[Nme{resp}];\n"
      )
      str_add(out$par) <- glue(
        "  vector{Ybounds}[N{resp}] Yl{resp};  // latent variable\n"
      )
      str_add(out$prior) <- glue(
        "  target += normal_lpdf(Y{resp}[Jme{resp}]",
        " | Yl{resp}[Jme{resp}],", 
        " noise{resp}[Jme{resp}]);\n"
      )
    }
  }
  if (is.formula(bterms$adforms$vreal)) {
    # vectors of real values for use in custom families
    vreal <- eval_rhs(bterms$adforms$vreal)
    k <- length(vreal$vars)
    str_add(out$data) <- cglue(
      "  // data for custom real vectors\n",
      "  vector[N{resp}] vreal{seq_len(k)}{resp};\n"
    )
  }
  if (is.formula(bterms$adforms$vint)) {
    # vectors of integer values for use in custom families
    vint <- eval_rhs(bterms$adforms$vint)
    k <- length(vint$vars)
    str_add(out$data) <- cglue(
      "  // data for custom integer vectors\n",
      "  int vint{seq_len(k)}{resp}[N{resp}];\n"
    )
  }
  out
}

# Stan code for ordinal thresholds
# intercepts in ordinal models require special treatment
# and must be present even when using non-linear predictors
# thus the relevant Stan code cannot be part of 'stan_fe'
stan_thres <- function(bterms, data, prior, ...) {
  stopifnot(is.btl(bterms) || is.btnl(bterms))
  out <- list()
  family <- bterms$family
  if (!is_ordinal(family)) {
    return(out)
  }
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  resp <- usc(px$resp)
  type <- str_if(has_ordered_thres(family), "ordered", "vector")
  gr <- grb <- ""
  groups <- get_thres_groups(bterms)
  if (has_thres_groups(bterms)) {
    # include one threshold vector per group
    gr <- usc(seq_along(groups))
    grb <- paste0("[", seq_along(groups), "]")
  }
  if (fix_intercepts(bterms)) {
    # identify ordinal mixtures by fixing their thresholds to the same values
    if (has_equidistant_thres(family)) {
      stop2("Cannot use equidistant and fixed thresholds at the same time.")
    }
    # separate definition from computation to support fixed parameters
    str_add(out$tpar_def) <- cglue(
      "  // ordinal thresholds\n",
      "  {type}[nthres{resp}{grb}] Intercept{p}{gr};\n"
    )
    str_add(out$tpar_comp) <- cglue(
      "  // fix thresholds across ordinal mixture components\n",
      "  Intercept{p}{gr} = fixed_Intercept{resp}{gr};\n"
    )
  } else {
    if (has_equidistant_thres(family)) {
      bound <- subset2(prior, class = "delta", group = "", ls = px)$bound
      for (i in seq_along(groups)) {
        str_add_list(out) <- stan_prior(
          prior, class = "Intercept", group = groups[i], 
          type = "real", prefix = "first_",
          suffix = glue("{p}{gr[i]}"), px = px, 
          comment = "first threshold"
        )
        str_add_list(out) <- stan_prior(
          prior, class = "delta", group = groups[i], 
          type = glue("real{bound}"), px = px, suffix = gr[i], 
          comment = "distance between thresholds"
        )
      }
      str_add(out$tpar_def) <- cglue(
        "  // temporary thresholds for centered predictors\n",
        "  {type}[nthres{resp}{grb}] Intercept{p}{gr};\n"
      )
      str_add(out$tpar_comp) <- cglue(
        "  // compute equidistant thresholds\n",
        "  for (k in 1:(nthres{resp}{grb})) {{\n",
        "    Intercept{p}{gr}[k] = first_Intercept{p}{gr}", 
        " + (k - 1.0) * delta{p}{gr};\n",
        "  }}\n"
      )
    } else {
      for (i in seq_along(groups)) {
        str_add_list(out) <- stan_prior(
          prior, class = "Intercept", group = groups[i], 
          coef = get_thres(bterms, group = groups[i]), 
          type = glue("{type}[nthres{resp}{grb[i]}]"),
          px = px, suffix = glue("{p}{gr[i]}"),
          comment = "temporary thresholds for centered predictors"
        )
      }
    }
  }
  if (has_thres_groups(bterms)) {
    # merge all group specific thresholds into one vector
    str_add(out$model_def) <- glue(
      "  vector[nmthres{resp}] merged_Intercept{p};  // merged thresholds\n"
    )
    grj <- seq_along(groups)
    grj <- glue("Kthres_start{resp}[{grj}]:Kthres_end{resp}[{grj}]")
    str_add(out$model_comp_basic) <- cglue(
      "  merged_Intercept{p}[{grj}] = Intercept{p}{gr};\n"
    )
  }
  sub_X_means <- ""
  if (stan_center_X(bterms) && length(all_terms(bterms$fe))) {
    # centering of the design matrix improves convergence
    # ordinal families either use thres - mu or mu - thres
    # both implies adding <mean_X, b> to the temporary intercept
    sub_X_means <- glue(" + dot_product(means_X{p}, b{p})")
  }
  str_add(out$gen_def) <- glue(
    "  // compute actual thresholds\n",
    "  vector[nthres{resp}{grb}] b{p}_Intercept{gr}",  
    " = Intercept{p}{gr}{sub_X_means};\n" 
  )
  out
}

# Stan code for the baseline functions of the Cox model
stan_bhaz <- function(bterms, prior, ...) {
  stopifnot(is.btl(bterms) || is.btnl(bterms))
  out <- list()
  if (!is_cox(bterms$family)) {
    return(out)
  }
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  resp <- usc(px$resp)
  str_add(out$data) <- glue(
    "  // data for flexible baseline functions\n",
    "  int Kbhaz{resp};  // number of basis functions\n",
    "  // design matrix of the baseline function\n",
    "  matrix[N{resp}, Kbhaz{resp}] Zbhaz{resp};\n",
    "  // design matrix of the cumulative baseline function\n",
    "  matrix[N{resp}, Kbhaz{resp}] Zcbhaz{resp};\n"
  )
  str_add_list(out) <- stan_prior(
    prior, class = "sbhaz", suffix = resp, px = px, 
    type = glue("vector<lower=0>[Kbhaz{resp}]"),
    comment = "baseline coefficients"
  )
  str_add(out$model_def) <- glue(
    "  // compute values of baseline function\n",
    "  vector[N{resp}] bhaz{resp} = Zbhaz{resp} * sbhaz{resp};\n",
    "  // compute values of cumulative baseline function\n",
    "  vector[N{resp}] cbhaz{resp} = Zcbhaz{resp} * sbhaz{resp};\n"
  )
  out
}

# Stan code specific to mixture families
stan_mixture <- function(bterms, data, prior) {
  out <- list()
  if (!is.mixfamily(bterms$family)) {
    return(out)
  }
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  nmix <- length(bterms$family$mix)
  theta_pred <- grepl("^theta", names(bterms$dpars))
  theta_pred <- bterms$dpars[theta_pred]
  theta_fix <- grepl("^theta", names(bterms$fdpars))
  theta_fix <- bterms$fdpars[theta_fix]
  def_thetas <- cglue(
    "  real<lower=0,upper=1> theta{1:nmix}{p};  // mixing proportion\n"
  )
  if (length(theta_pred)) {
    if (length(theta_pred) != nmix - 1) {
      stop2("Can only predict all but one mixing proportion.")
    }
    missing_id <- setdiff(1:nmix, dpar_id(names(theta_pred)))
    str_add(out$model_def) <- glue(
      "  vector[N{p}] theta{missing_id}{p} = rep_vector(0, N{p});\n",                   
      "  real log_sum_exp_theta;\n"      
    )
    sum_exp_theta <- glue("exp(theta{1:nmix}{p}[n])", collapse = " + ")
    str_add(out$model_comp_mix) <- glue(
      "  for (n in 1:N{p}) {{\n",
      "    // scale theta to become a probability vector\n",
      "    log_sum_exp_theta = log({sum_exp_theta});\n"
    )
    str_add(out$model_comp_mix) <- cglue(
      "    theta{1:nmix}{p}[n] = theta{1:nmix}{p}[n] - log_sum_exp_theta;\n"
    )
    str_add(out$model_comp_mix) <- "  }\n"
  } else if (length(theta_fix)) {
    # fix mixture proportions
    if (length(theta_fix) != nmix) {
      stop2("Can only fix no or all mixing proportions.")
    }
    str_add(out$data) <- "  // mixing proportions\n"
    str_add(out$data) <- cglue(
      "  real<lower=0,upper=1> theta{1:nmix}{p};\n"
    )
  } else {
    # estimate mixture proportions
    str_add(out$data) <- glue(
      "  vector[{nmix}] con_theta{p};  // prior concentration\n"                  
    )
    str_add(out$par) <- glue(
      "  simplex[{nmix}] theta{p};  // mixing proportions\n"
    )
    str_add(out$prior) <- glue(
      "  target += dirichlet_lpdf(theta{p} | con_theta{p});\n"                
    )
    # separate definition from computation to support fixed parameters
    str_add(out$tpar_def) <- "  // mixing proportions\n"
    str_add(out$tpar_def) <- cglue(
      "  real<lower=0,upper=1> theta{1:nmix}{p};\n"
    )
    str_add(out$tpar_comp) <- cglue(
      "  theta{1:nmix}{p} = theta{p}[{1:nmix}];\n"
    )
  }
  if (order_intercepts(bterms)) {
    # identify mixtures by ordering the intercepts of their components
    str_add(out$par) <- glue( 
      "  ordered[{nmix}] ordered_Intercept{p};  // to identify mixtures\n"
    )
  }
  if (fix_intercepts(bterms)) {
    # identify ordinal mixtures by fixing their thresholds to the same values
    stopifnot(is_ordinal(bterms))
    gr <- grb <- ""
    groups <- get_thres_groups(bterms)
    if (has_thres_groups(bterms)) {
      # include one threshold vector per group
      gr <- usc(seq_along(groups))
      grb <- paste0("[", seq_along(groups), "]")
    }
    type <- str_if(has_ordered_thres(bterms), "ordered", "vector")
    for (i in seq_along(groups)) {
      str_add_list(out) <- stan_prior(
        prior, class = "Intercept", px = px,  
        coef = get_thres(bterms, group = groups[i]), 
        type = glue("{type}[nthres{p}{grb[i]}]"),
        prefix = "fixed_", suffix = glue("{p}{gr[i]}"),
        comment = "thresholds fixed over mixture components"
      )
    }
  }
  out
}

# ordinal log-probability densitiy functions in Stan language
# @return a character string
stan_ordinal_lpmf <- function(family, link) {
  stopifnot(is.character(family), is.character(link))
  ilink <- stan_ilink(link)
  th <- function(k) {
    # helper function generating stan code inside ilink(.)
    if (family %in% c("cumulative", "sratio")) {
      out <- glue("thres[{k}] - mu")
    } else if (family %in% c("cratio", "acat")) {
      out <- glue("mu - thres[{k}]")
    }
    glue("disc * ({out})")
  }
  out <- glue(
    "  /* {family}-{link} log-PDF for a single response\n",
    "   * Args:\n",
    "   *   y: response category\n",
    "   *   mu: latent mean parameter\n",
    "   *   disc: discrimination parameter\n",
    "   *   thres: ordinal thresholds\n",
    "   * Returns:\n", 
    "   *   a scalar to be added to the log posterior\n",
    "   */\n",
    "   real {family}_{link}_lpmf(int y, real mu, real disc, vector thres) {{\n"
  )
  # define the function body
  if (family == "cumulative") {
    str_add(out) <- glue(
      "     int nthres = num_elements(thres);\n",
      "     real p;\n",
      "     if (y == 1) {{\n",
      "       p = {ilink}({th(1)});\n",
      "     }} else if (y == nthres + 1) {{\n",
      "       p = 1 - {ilink}({th('nthres')});\n",
      "     }} else {{\n",
      "       p = {ilink}({th('y')}) -\n",
      "           {ilink}({th('y - 1')});\n",
      "     }}\n",
      "     return log(p);\n",
      "   }}\n"
    )
  } else if (family %in% c("sratio", "cratio")) {
    sc <- str_if(family == "sratio", "1 - ")
    str_add(out) <- glue(
      "     int nthres = num_elements(thres);\n",
      "     vector[nthres + 1] p;\n",
      "     vector[nthres] q;\n",
      "     int k = 1;\n",
      "     while (k <= min(y, nthres)) {{\n",
      "       q[k] = {sc}{ilink}({th('k')});\n",
      "       p[k] = 1 - q[k];\n",
      "       for (kk in 1:(k - 1)) p[k] = p[k] * q[kk];\n", 
      "       k += 1;\n",
      "     }}\n",
      "     if (y == nthres + 1) {{\n",
      "       p[nthres + 1] = prod(q);\n",
      "     }}\n",
      "     return log(p[y]);\n",
      "   }}\n"
    )
  } else if (family == "acat") {
    if (ilink == "inv_logit") {
      str_add(out) <- glue(
        "     int nthres = num_elements(thres);\n",
        "     vector[nthres + 1] p;\n",
        "     p[1] = 0.0;\n",
        "     for (k in 1:(nthres)) {{\n",
        "       p[k + 1] = p[k] + {th('k')};\n",
        "     }}\n",
        "     p = exp(p);\n",
        "     return log(p[y] / sum(p));\n",
        "   }}\n"
      )
    } else {
      str_add(out) <- glue(   
        "     int nthres = num_elements(thres);\n",
        "     vector[nthres + 1] p;\n",
        "     vector[nthres] q;\n",
        "     for (k in 1:(nthres))\n",
        "       q[k] = {ilink}({th('k')});\n",
        "     for (k in 1:(nthres + 1)) {{\n",     
        "       p[k] = 1.0;\n",
        "       for (kk in 1:(k - 1)) p[k] = p[k] * q[kk];\n",
        "       for (kk in k:(nthres)) p[k] = p[k] * (1 - q[kk]);\n",   
        "     }}\n",
        "     return log(p[y] / sum(p));\n",
        "   }}\n"
      )
    }
  }
  # lpdf function for multiple merged thresholds
  str_add(out) <- glue(
    # TODO: include order_logistic_merged_lpdf
    "  /* {family}-{link} log-PDF for a single response and merged thresholds\n",
    "   * Args:\n",
    "   *   y: response category\n",
    "   *   mu: latent mean parameter\n",
    "   *   disc: discrimination parameter\n",
    "   *   thres: vector of merged ordinal thresholds\n",
    "   *   j: start and end index for the applid threshold within 'thres'\n",
    "   * Returns:\n", 
    "   *   a scalar to be added to the log posterior\n",
    "   */\n",
    "   real {family}_{link}_merged_lpmf(", 
    "int y, real mu, real disc, vector thres, int[] j) {{\n",
    "     return {family}_{link}_lpmf(y | mu, disc, thres[j[1]:j[2]]);\n",
    "   }}\n"
  )
  out
}
