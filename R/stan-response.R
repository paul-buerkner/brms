# unless otherwise specifiedm functions return a named list
# of Stan code snippets to be pasted together later on

# Stan code for the response variables
stan_response <- function(bterms, data, normalize) {
  stopifnot(is.brmsterms(bterms))
  lpdf <- stan_lpdf_name(normalize)
  family <- bterms$family
  rtype <- str_if(use_int(family), "int", "real")
  multicol <- has_multicol(family)
  px <- check_prefix(bterms)
  resp <- usc(combine_prefix(px))
  out <- list(resp_type = rtype)
  if (nzchar(resp)) {
    # global N is defined elsewhere
    str_add(out$data) <- glue(
      "  int<lower=1> N{resp};  // number of observations\n"
    )
    str_add(out$pll_def) <- glue(
      "  int N{resp} = end - start + 1;\n"
    )
  }
  if (has_cat(family)) {
    str_add(out$data) <- glue(
      "  int<lower=2> ncat{resp};  // number of categories\n"
    )
    str_add(out$pll_args) <- glue(", data int ncat{resp}")
  }
  if (has_multicol(family)) {
    if (rtype == "real") {
      str_add(out$data) <- glue(
        "  array[N{resp}] vector[ncat{resp}] Y{resp};  // response array\n"
      )
      str_add(out$pll_args) <- glue(", data array[] vector Y{resp}")
    } else if (rtype == "int") {
      str_add(out$data) <- glue(
        "  array[N{resp}, ncat{resp}] int Y{resp};  // response array\n"
      )
      str_add(out$pll_args) <- glue(", data array[,] int Y{resp}")
    }
  } else {
    if (rtype == "real") {
      # type vector (instead of array real) is required by some PDFs
      str_add(out$data) <- glue(
        "  vector[N{resp}] Y{resp};  // response variable\n"
      )
      str_add(out$pll_args) <- glue(", data vector Y{resp}")
    } else if (rtype == "int") {
      str_add(out$data) <- glue(
        "  array[N{resp}] int Y{resp};  // response variable\n"
      )
      str_add(out$pll_args) <- glue(", data array[] int Y{resp}")
    }
  }
  if (has_ndt(family)) {
    str_add(out$tdata_def) <- glue(
      "  real min_Y{resp} = min(Y{resp});\n"
    )
  }
  if (has_trials(family) || is.formula(bterms$adforms$trials)) {
    str_add(out$data) <- glue(
      "  array[N{resp}] int trials{resp};  // number of trials\n"
    )
    str_add(out$pll_args) <- glue(", data array[] int trials{resp}")
  }
  if (is.formula(bterms$adforms$weights)) {
    str_add(out$data) <- glue(
      "  vector<lower=0>[N{resp}] weights{resp};  // model weights\n"
    )
    str_add(out$pll_args) <- glue(", data vector weights{resp}")
  }
  if (has_thres(family)) {
    groups <- get_thres_groups(family)
    if (any(nzchar(groups))) {
      str_add(out$data) <- glue(
        "  int<lower=1> ngrthres{resp};  // number of threshold groups\n",
        "  array[ngrthres{resp}] int<lower=1> nthres{resp};  // number of thresholds\n",
        "  array[N{resp}, 2] int<lower=1> Jthres{resp};  // threshold indices\n"
      )
      str_add(out$tdata_def) <- glue(
        "  int<lower=1> nmthres{resp} = sum(nthres{resp});",
        "  // total number of thresholds\n",
        "  array[ngrthres{resp}] int<lower=1> Kthres_start{resp};",
        "  // start index per threshold group\n",
        "  array[ngrthres{resp}] int<lower=1> Kthres_end{resp};",
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
      str_add(out$pll_args) <- glue(
        ", data array[] int nthres{resp}, data array[,] int Jthres{resp}"
      )
    } else {
      str_add(out$data) <- glue(
        "  int<lower=2> nthres{resp};  // number of thresholds\n"
      )
      str_add(out$pll_args) <- glue(", data int nthres{resp}")
    }
  }
  if (is.formula(bterms$adforms$se)) {
    str_add(out$data) <- glue(
      "  vector<lower=0>[N{resp}] se{resp};  // known sampling error\n"
    )
    str_add(out$tdata_def) <- glue(
      "  vector<lower=0>[N{resp}] se2{resp} = square(se{resp});\n"
    )
    str_add(out$pll_args) <- glue(
      ", data vector se{resp}, data vector se2{resp}"
    )
  }
  if (is.formula(bterms$adforms$dec)) {
    str_add(out$data) <- glue(
      "  array[N{resp}] int<lower=0,upper=1> dec{resp};  // decisions\n"
    )
    str_add(out$pll_args) <- glue(", data array[] int dec{resp}")
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
    str_add(out$pll_args) <- glue(
      ", data vector denom{resp}, data vector log_denom{resp}"
    )
  }
  if (is.formula(bterms$adforms$cens)) {
    str_add(out$data) <- glue(
      "  array[N{resp}] int<lower=-1,upper=2> cens{resp};  // indicates censoring\n"
    )
    str_add(out$pll_args) <- glue(", data array[] int cens{resp}")
    y2_expr <- get_ad_expr(bterms, "cens", "y2")
    if (!is.null(y2_expr)) {
      # interval censoring is required
      if (rtype == "int") {
        str_add(out$data) <- glue(
          "  array[N{resp}] int rcens{resp};"
        )
        str_add(out$pll_args) <- glue(", data array[] int rcens{resp}")
      } else {
        str_add(out$data) <- glue(
          "  vector[N{resp}] rcens{resp};"
        )
        str_add(out$pll_args) <- glue(", data vector rcens{resp}")
      }
      str_add(out$data) <- "  // right censor points for interval censoring\n"
    }
  }
  bounds <- trunc_bounds(bterms, data = data)
  if (any(bounds$lb > -Inf)) {
    str_add(out$data) <- glue(
      "  array[N{resp}] {rtype} lb{resp};  // lower truncation bounds;\n"
    )
    str_add(out$pll_args) <- glue(", data array[] {rtype} lb{resp}")
  }
  if (any(bounds$ub < Inf)) {
    str_add(out$data) <- glue(
      "  array[N{resp}] {rtype} ub{resp};  // upper truncation bounds\n"
    )
    str_add(out$pll_args) <- glue(", data array[] {rtype} ub{resp}")
  }
  if (is.formula(bterms$adforms$mi)) {
    # TODO: pass 'Ybounds' via 'standata' instead of hardcoding them
    Ybounds <- trunc_bounds(bterms, data, incl_family = TRUE, stan = TRUE)
    sdy <- get_sdy(bterms, data)
    if (is.null(sdy)) {
      # response is modeled without measurement error
      str_add(out$data) <- glue(
        "  int<lower=0> Nmi{resp};  // number of missings\n",
        "  array[Nmi{resp}] int<lower=1> Jmi{resp};  // positions of missings\n"
      )
      str_add(out$par) <- glue(
        "  vector{Ybounds}[Nmi{resp}] Ymi{resp};  // estimated missings\n"
      )
      str_add(out$model_no_pll_def) <- glue(
        "  // vector combining observed and missing responses\n",
        "  vector[N{resp}] Yl{resp} = Y{resp};\n"
      )
      str_add(out$model_no_pll_comp_basic) <- glue(
        "  Yl{resp}[Jmi{resp}] = Ymi{resp};\n"
      )
      str_add(out$pll_args) <- glue(", vector Yl{resp}")
    } else {
      str_add(out$data) <- glue(
        "  // data for measurement-error in the response\n",
        "  vector<lower=0>[N{resp}] noise{resp};\n",
        "  // information about non-missings\n",
        "  int<lower=0> Nme{resp};\n",
        "  array[Nme{resp}] int<lower=1> Jme{resp};\n"
      )
      str_add(out$par) <- glue(
        "  vector{Ybounds}[N{resp}] Yl{resp};  // latent variable\n"
      )
      str_add(out$model_prior) <- glue(
        "  target += normal_{lpdf}(Y{resp}[Jme{resp}]",
        " | Yl{resp}[Jme{resp}], noise{resp}[Jme{resp}]);\n"
      )
      str_add(out$pll_args) <- glue(", vector Yl{resp}")
    }
  }
  if (is.formula(bterms$adforms$vreal)) {
    # vectors of real values for use in custom families
    vreal <- eval_rhs(bterms$adforms$vreal)
    k <- length(vreal$vars)
    str_add(out$data) <- cglue(
      "  // data for custom real vectors\n",
      "  array[N{resp}] real vreal{seq_len(k)}{resp};\n"
    )
    str_add(out$pll_args) <- cglue(", data array[] real vreal{seq_len(k)}{resp}")
  }
  if (is.formula(bterms$adforms$vint)) {
    # vectors of integer values for use in custom families
    vint <- eval_rhs(bterms$adforms$vint)
    k <- length(vint$vars)
    str_add(out$data) <- cglue(
      "  // data for custom integer vectors\n",
      "  array[N{resp}] int vint{seq_len(k)}{resp};\n"
    )
    str_add(out$pll_args) <- cglue(", data array[] int vint{seq_len(k)}{resp}")
  }
  out
}

# Stan code for ordinal thresholds
# intercepts in ordinal models require special treatment
# and must be present even when using non-linear predictors
# thus the relevant Stan code cannot be part of 'stan_fe'
stan_thres <- function(bterms, data, prior, normalize, ...) {
  stopifnot(is.btl(bterms) || is.btnl(bterms))
  out <- list()
  if (!is_ordinal(bterms)) {
    return(out)
  }
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  resp <- usc(px$resp)
  type <- str_if(has_ordered_thres(bterms), "ordered", "vector")
  coef_type <- str_if(has_ordered_thres(bterms), "", "real")
  gr <- grb <- ""
  groups <- get_thres_groups(bterms)
  if (has_thres_groups(bterms)) {
    # include one threshold vector per group
    gr <- usc(seq_along(groups))
    grb <- paste0("[", seq_along(groups), "]")
  }
  if (fix_intercepts(bterms)) {
    # identify ordinal mixtures by fixing their thresholds to the same values
    if (has_equidistant_thres(bterms)) {
      stop2("Cannot use equidistant and fixed thresholds at the same time.")
    }
    # separate definition from computation to support fixed parameters
    str_add(out$tpar_def) <- "  // ordinal thresholds\n"
    str_add(out$tpar_def) <- cglue(
      "  {type}[nthres{resp}{grb}] Intercept{p}{gr};\n"
    )
    str_add(out$tpar_comp) <-
      "  // fix thresholds across ordinal mixture components\n"
    str_add(out$tpar_comp) <- cglue(
      "  Intercept{p}{gr} = fixed_Intercept{resp}{gr};\n"
    )
  } else {
    if (has_equidistant_thres(bterms)) {
      bound <- subset2(prior, class = "delta", group = "", ls = px)$bound
      for (i in seq_along(groups)) {
        str_add_list(out) <- stan_prior(
          prior, class = "Intercept", group = groups[i],
          prefix = "first_", suffix = glue("{p}{gr[i]}"), px = px,
          comment = "first threshold", normalize = normalize
        )
        str_add_list(out) <- stan_prior(
          prior, class = "delta", group = groups[i], px = px, suffix = glue("{p}{gr[i]}"),
          comment = "distance between thresholds", normalize = normalize
        )
      }
      str_add(out$tpar_def) <-
        "  // temporary thresholds for centered predictors\n"
      str_add(out$tpar_def) <- cglue(
        "  {type}[nthres{resp}{grb}] Intercept{p}{gr};\n"
      )
      str_add(out$tpar_comp) <-
        "  // compute equidistant thresholds\n"
      str_add(out$tpar_comp) <- cglue(
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
          coef_type = coef_type, px = px, suffix = glue("{p}{gr[i]}"),
          comment = "temporary thresholds for centered predictors",
          normalize = normalize
        )
      }
    }
  }
  stz <- ""
  if (has_sum_to_zero_thres(bterms)) {
    stz <- "_stz"
    str_add(out$tpar_def) <- cglue(
      "  vector[nthres{resp}{grb}] Intercept{p}_stz{gr};",
      "  // sum-to-zero constraint thresholds\n"
    )
    str_add(out$tpar_comp) <- "  // compute sum-to-zero constraint thresholds\n"
    str_add(out$tpar_comp) <- cglue(
      "  Intercept{p}_stz{gr} = Intercept{p}{gr} - mean(Intercept{p}{gr});\n"
    )
  }
  if (has_thres_groups(bterms)) {
    # merge all group specific thresholds into one vector
    str_add(out$tpar_def) <- glue(
      "  vector[nmthres{resp}] merged_Intercept{p}{stz};  // merged thresholds\n"
    )
    str_add(out$tpar_comp) <- "  // merge thresholds\n"
    grj <- seq_along(groups)
    grj <- glue("Kthres_start{resp}[{grj}]:Kthres_end{resp}[{grj}]")
    str_add(out$tpar_comp) <- cglue(
      "  merged_Intercept{p}{stz}[{grj}] = Intercept{p}{stz}{gr};\n"
    )
    str_add(out$pll_args) <- cglue(", vector merged_Intercept{p}{stz}")
  } else {
    str_add(out$pll_args) <- glue(", vector Intercept{p}{stz}")
  }
  sub_X_means <- ""
  if (stan_center_X(bterms) && length(all_terms(bterms$fe))) {
    # centering of the design matrix improves convergence
    # ordinal families either use thres - mu or mu - thres
    # both implies adding <mean_X, b> to the temporary intercept
    sub_X_means <- glue(" + dot_product(means_X{p}, b{p})")
  }
  str_add(out$gen_def) <- "  // compute actual thresholds\n"
  str_add(out$gen_def) <- cglue(
    "  vector[nthres{resp}{grb}] b{p}_Intercept{gr}",
    " = Intercept{p}{stz}{gr}{sub_X_means};\n"
  )
  out
}

# Stan code for the baseline functions of the Cox model
stan_bhaz <- function(bterms, prior, threads, normalize, ...) {
  stopifnot(is.btl(bterms) || is.btnl(bterms))
  out <- list()
  if (!is_cox(bterms$family)) {
    return(out)
  }
  lpdf <- stan_lpdf_name(normalize)
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  resp <- usc(px$resp)
  slice <- stan_slice(threads)
  str_add(out$data) <- glue(
    "  // data for flexible baseline functions\n",
    "  int Kbhaz{resp};  // number of basis functions\n",
    "  // design matrix of the baseline function\n",
    "  matrix[N{resp}, Kbhaz{resp}] Zbhaz{resp};\n",
    "  // design matrix of the cumulative baseline function\n",
    "  matrix[N{resp}, Kbhaz{resp}] Zcbhaz{resp};\n",
    "  // a-priori concentration vector of baseline coefficients\n",
    "  vector<lower=0>[Kbhaz{resp}] con_sbhaz{resp};\n"
  )
  str_add(out$par) <- glue(
    "  simplex[Kbhaz{resp}] sbhaz{resp};  // baseline coefficients\n"
  )
  str_add(out$tpar_prior) <- glue(
    "  lprior += dirichlet_{lpdf}(sbhaz{resp} | con_sbhaz{resp});\n"
  )
  str_add(out$model_def) <- glue(
    "  // compute values of baseline function\n",
    "  vector[N{resp}] bhaz{resp} = Zbhaz{resp}{slice} * sbhaz{resp};\n",
    "  // compute values of cumulative baseline function\n",
    "  vector[N{resp}] cbhaz{resp} = Zcbhaz{resp}{slice} * sbhaz{resp};\n"
  )
  str_add(out$pll_args) <- glue(
    ", data matrix Zbhaz{resp}, data matrix Zcbhaz{resp}, vector sbhaz{resp}"
  )
  out
}

# Stan code specific to mixture families
stan_mixture <- function(bterms, data, prior, threads, normalize, ...) {
  out <- list()
  if (!is.mixfamily(bterms$family)) {
    return(out)
  }
  lpdf <- stan_lpdf_name(normalize)
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
      "  vector[N{p}] theta{missing_id}{p} = rep_vector(0.0, N{p});\n",
      "  real log_sum_exp_theta{p};\n"
    )
    sum_exp_theta <- glue("exp(theta{1:nmix}{p}[n])", collapse = " + ")
    str_add(out$model_comp_mix) <- glue(
      "  for (n in 1:N{p}) {{\n",
      "    // scale theta to become a probability vector\n",
      "    log_sum_exp_theta{p} = log({sum_exp_theta});\n"
    )
    str_add(out$model_comp_mix) <- cglue(
      "    theta{1:nmix}{p}[n] = theta{1:nmix}{p}[n] - log_sum_exp_theta{p};\n"
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
    str_add(out$pll_args) <- cglue(", real theta{1:nmix}{p}")
  } else {
    # estimate mixture proportions
    str_add(out$data) <- glue(
      "  vector[{nmix}] con_theta{p};  // prior concentration\n"
    )
    str_add(out$par) <- glue(
      "  simplex[{nmix}] theta{p};  // mixing proportions\n"
    )
    str_add(out$tpar_prior) <- glue(
      "  lprior += dirichlet_{lpdf}(theta{p} | con_theta{p});\n"
    )
    # separate definition from computation to support fixed parameters
    str_add(out$tpar_def) <- "  // mixing proportions\n"
    str_add(out$tpar_def) <- cglue(
      "  real<lower=0,upper=1> theta{1:nmix}{p};\n"
    )
    str_add(out$tpar_comp) <- cglue(
      "  theta{1:nmix}{p} = theta{p}[{1:nmix}];\n"
    )
    str_add(out$pll_args) <- cglue(", real theta{1:nmix}{p}")
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
    coef_type <- str_if(has_ordered_thres(bterms), "", "real")
    for (i in seq_along(groups)) {
      str_add_list(out) <- stan_prior(
        prior, class = "Intercept",
        coef = get_thres(bterms, group = groups[i]),
        type = glue("{type}[nthres{p}{grb[i]}]"),
        coef_type = coef_type, px = px,
        prefix = "fixed_", suffix = glue("{p}{gr[i]}"),
        comment = "thresholds fixed over mixture components",
        normalize = normalize
      )
    }
  }
  out
}

# ordinal log-probability densitiy functions in Stan language
# @return a character string
stan_ordinal_lpmf <- function(family, link) {
  stopifnot(is.character(family), is.character(link))
  inv_link <- stan_inv_link(link, vectorize = FALSE)
  th <- function(k) {
    # helper function generating stan code inside inv_link(.)
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
    if (inv_link == "inv_logit") {
      str_add(out) <- glue(
        "     int nthres = num_elements(thres);\n",
        "     if (y == 1) {{\n",
        "       return log_inv_logit({th(1)});\n",
        "     }} else if (y == nthres + 1) {{\n",
        "       return log1m_inv_logit({th('nthres')});\n",
        "     }} else {{\n",
        # TODO: replace with log_inv_logit_diff once rstan >= 2.25
        "       return log_diff_exp(\n",
        "         log_inv_logit({th('y')}), \n",
        "         log_inv_logit({th('y - 1')})\n",
        "       );\n",
        "     }}\n",
        "   }}\n"
      )
    } else {
      str_add(out) <- glue(
        "     int nthres = num_elements(thres);\n",
        "     real p;\n",
        "     if (y == 1) {{\n",
        "       p = {inv_link}({th(1)});\n",
        "     }} else if (y == nthres + 1) {{\n",
        "       p = 1 - {inv_link}({th('nthres')});\n",
        "     }} else {{\n",
        "       p = {inv_link}({th('y')}) -\n",
        "           {inv_link}({th('y - 1')});\n",
        "     }}\n",
        "     return log(p);\n",
        "   }}\n"
      )
    }
  } else if (family %in% c("sratio", "cratio")) {
    # TODO: support 'softit' link as well
    if (inv_link == "inv_cloglog") {
      qk <- str_if(
        family == "sratio",
        "-exp({th('k')})",
        "log1m_exp(-exp({th('k')}))"
      )
    } else if (inv_link == "inv_logit") {
      qk <- str_if(
        family == "sratio",
        "log1m_inv_logit({th('k')})",
        "log_inv_logit({th('k')})"
      )
    } else if (inv_link == "Phi") {
      # TODO: replace with more stable std_normal_lcdf once rstan >= 2.25
      qk <- str_if(
        family == "sratio",
        "normal_lccdf({th('k')}|0,1)",
        "normal_lcdf({th('k')}|0,1)"
      )
    } else if (inv_link == "Phi_approx") {
      qk <- str_if(
        family == "sratio",
        "log1m_inv_logit(0.07056 * pow({th('k')}, 3.0) + 1.5976 * {th('k')})",
        "log_inv_logit(0.07056 * pow({th('k')}, 3.0) + 1.5976 * {th('k')})"
      )
    } else if (inv_link == "inv_cauchit") {
      qk <- str_if(
        family == "sratio",
        "cauchy_lccdf({th('k')}|0,1)",
        "cauchy_lcdf({th('k')}|0,1)"
      )
    }
    qk <- glue(qk)
    str_add(out) <- glue(
      "     int nthres = num_elements(thres);\n",
      "     vector[nthres + 1] p;\n",
      "     vector[nthres] q;\n",
      "     int k = 1;\n",
      "     while (k <= min(y, nthres)) {{\n",
      "       q[k] = {qk};\n",
      "       p[k] = log1m_exp(q[k]);\n",
      "       for (kk in 1:(k - 1)) p[k] = p[k] + q[kk];\n",
      "       k += 1;\n",
      "     }}\n",
      "     if (y == nthres + 1) {{\n",
      "       p[nthres + 1] = sum(q);\n",
      "     }}\n",
      "     return p[y];\n",
      "   }}\n"
    )
  } else if (family == "acat") {
    if (inv_link == "inv_logit") {
      str_add(out) <- glue(
        "     int nthres = num_elements(thres);\n",
        "     vector[nthres + 1] p = append_row(0, cumulative_sum(disc * (mu - thres)));\n",
        "     return p[y] - log_sum_exp(p);\n",
        "   }}\n"
      )
    } else {
      str_add(out) <- glue(
        "     int nthres = num_elements(thres);\n",
        "     vector[nthres + 1] p;\n",
        "     vector[nthres] q;\n",
        "     for (k in 1:(nthres))\n",
        "       q[k] = {inv_link}({th('k')});\n",
        "     for (k in 1:(nthres + 1)) {{\n",
        "       p[k] = 1.0;\n",
        "       for (kk in 1:(k - 1)) p[k] = p[k] * q[kk];\n",
        "       for (kk in k:(nthres)) p[k] = p[k] * (1 - q[kk]);\n",
        "     }}\n",
        "     return log(p[y]) - log(sum(p));\n",
        "   }}\n"
      )
    }
  }
  # lpmf function for multiple merged thresholds
  str_add(out) <- glue(
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
    "int y, real mu, real disc, vector thres, array[] int j) {{\n",
    "     return {family}_{link}_lpmf(y | mu, disc, thres[j[1]:j[2]]);\n",
    "   }}\n"
  )
  if (family == "cumulative" && link == "logit") {
    # use the more efficient 'ordered_logistic' built-in function
    str_add(out) <- glue(
      "  /* ordered-logistic log-PDF for a single response and merged thresholds\n",
      "   * Args:\n",
      "   *   y: response category\n",
      "   *   mu: latent mean parameter\n",
      "   *   thres: vector of merged ordinal thresholds\n",
      "   *   j: start and end index for the applid threshold within 'thres'\n",
      "   * Returns:\n",
      "   *   a scalar to be added to the log posterior\n",
      "   */\n",
      "   real ordered_logistic_merged_lpmf(",
      "int y, real mu, vector thres, array[] int j) {{\n",
      "     return ordered_logistic_lpmf(y | mu, thres[j[1]:j[2]]);\n",
      "   }}\n"
    )
  }
  out
}

# log probability density for hurdle ordinal models
# @return a character string
stan_hurdle_ordinal_lpmf <- function(family, link) {
  stopifnot(is.character(family), is.character(link))
  # TODO: generalize to non-cumulative families?
  stopifnot(family == "hurdle_cumulative")
  inv_link <- stan_inv_link(link, vectorize = FALSE)
  th <- function(k) {
    out <- glue("thres[{k}] - mu")
    glue("disc * ({out})")
  }
  out <- glue(
    "  /* {family}-{link} log-PDF for a single response\n",
    "   * Args:\n",
    "   *   y: response category\n",
    "   *   mu: latent mean parameter\n",
    "   *   hu: hurdle probability\n",
    "   *   disc: discrimination parameter\n",
    "   *   thres: ordinal thresholds\n",
    "   * Returns:\n",
    "   *   a scalar to be added to the log posterior\n",
    "   */\n",
    "   real {family}_{link}_lpmf(int y, real mu, real hu, real disc, vector thres) {{\n",
    "\n"
  )
  # define the function body
  if (inv_link == "inv_logit") {
    str_add(out) <- glue(
      "     int nthres = num_elements(thres);\n",
      "     if (y == 0) {{\n",
      "       return bernoulli_lpmf(1 | hu);\n",
      "     }} else if (y == 1) {{\n",
      "       return log_inv_logit({th(1)}) +\n",
      "                bernoulli_lpmf(0 | hu);\n",
      "     }} else if (y == nthres + 2) {{\n",
      "       return log1m_inv_logit({th('nthres')}) +\n",
      "                bernoulli_lpmf(0 | hu);\n",
      "     }} else {{\n",
      # TODO: replace with log_inv_logit_diff once rstan >= 2.25
      "       return log_diff_exp(\n",
      "         log_inv_logit({th('y')}), \n",
      "         log_inv_logit({th('y - 1')})\n",
      "       ) + bernoulli_lpmf(0 | hu) ;\n",
      "     }}\n",
      "   }}\n"
    )
  } else {
    str_add(out) <- glue(
      "     int nthres = num_elements(thres);\n",
      "     real p;\n",
      "     if (y == 0){{\n",
      "       p = hu;\n",
      "     }} else if (y == 1) {{\n",
      "       p = {inv_link}({th(1)}) * (1 - hu);\n",
      "     }} else if (y == nthres + 1) {{\n",
      "       p = (1 - {inv_link}({th('nthres')})) * (1 - hu);\n",
      "     }} else {{\n",
      "       p = ({inv_link}({th('y')}) -\n",
      "           {inv_link}({th('y - 1')})) * (1 - hu);\n",
      "     }}\n",
      "     return log(p);\n",
      "   }}\n"
    )
  }

  # lpmf function for multiple merged thresholds
  str_add(out) <- glue(
    "  /* {family}-{link} log-PDF for a single response and merged thresholds\n",
    "   * Args:\n",
    "   *   y: response category\n",
    "   *   mu: latent mean parameter\n",
    "   *   hu: hurdle probability\n",
    "   *   disc: discrimination parameter\n",
    "   *   thres: vector of merged ordinal thresholds\n",
    "   *   j: start and end index for the applid threshold within 'thres'\n",
    "   * Returns:\n",
    "   *   a scalar to be added to the log posterior\n",
    "   */\n",
    "   real {family}_{link}_merged_lpmf(",
    "int y, real mu, real hu, real disc, vector thres, array[] int j) {{\n",
    "     return {family}_{link}_lpmf(y | mu, hu, disc, thres[j[1]:j[2]]);\n",
    "   }}\n"
  )

  if (link == "logit") {
    # use the more efficient ordered_logistic function when disc == 1
    str_add(out) <- glue(
      "\n",
      " // Use more efficient ordered_logistic function with disc == 1\n",
      "   real hurdle_cumulative_ordered_logistic_lpmf(int y, real mu, real hu, real disc, vector thres) {{\n",
      "     if (y == 0) {{\n",
      "       return bernoulli_lpmf(1 | hu);\n",
      "     }} else {{\n",
      "       return ordered_logistic_lpmf(y | mu, thres) +\n",
      "                bernoulli_lpmf(0 | hu);\n",
      "     }}\n",
      "   }}\n"
    )
    str_add(out) <- glue(
      "  /* use ordered-logistic log-PDF for a single response and merged thresholds\n",
      "   * Args:\n",
      "   *   y: response category\n",
      "   *   mu: latent mean parameter\n",
      "   *   hu: hurdle probability\n",
      "   *   thres: vector of merged ordinal thresholds\n",
      "   *   j: start and end index for the applid threshold within 'thres'\n",
      "   * Returns:\n",
      "   *   a scalar to be added to the log posterior\n",
      "   */\n",
      "   real hurdle_cumulative_ordered_logistic_merged_lpmf(",
      "int y, real mu, real hu, real disc, vector thres, array[] int j) {{\n",
      "     return hurdle_cumulative_ordered_logistic_lpmf(y | mu, hu, disc, thres[j[1]:j[2]]);\n",
      "   }}\n"
    )
  }
  out
}

