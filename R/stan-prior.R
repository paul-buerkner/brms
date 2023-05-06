# unless otherwise specified, functions return a single character
# string defining the likelihood of the model in Stan language

# Define priors for parameters in Stan language
# @param prior an object of class 'brmsprior'
# @param class the parameter class
# @param coef the coefficients of this class
# @param group the name of a grouping factor
# @param type Stan type used in the definition of the parameter
#   if type is empty the parameter is not initialized inside 'stan_prior'
# @param dim stan array dimension to be specified after the parameter name
#   cannot be expressed via 'suffix' as the latter should apply to
#   individual coefficients while 'dim' should not
#   TODO: decide whether to support arrays for parameters at all
#   an alternative would be to specify elements directly as parameters
# @param coef_type Stan type used in the definition of individual parameter
#   coefficients
# @param prefix a prefix to put at the parameter class
# @param suffix a suffix to put at the parameter class
# @param broadcast Stan type to which the prior should be broadcasted
#   in order to handle vectorized prior statements
#   supported values are 'vector' or 'matrix'
# @param comment character string containing a comment for the parameter
# @param px list or data.frame after which to subset 'prior'
# @return a named list of character strings in Stan language
stan_prior <- function(prior, class, coef = NULL, group = NULL,
                       type = "real", dim = "", coef_type = "real",
                       prefix = "", suffix = "", broadcast = "vector",
                       header_type = "", comment = "", px = list(),
                       normalize = TRUE) {
  prior_only <- isTRUE(attr(prior, "sample_prior") == "only")
  prior <- subset2(
    prior, class = class, coef = c(coef, ""),
    group = c(group, ""), ls = px
  )
  # special priors cannot be passed literally to Stan
  is_special_prior <- is_special_prior(prior$prior)
  if (any(is_special_prior)) {
    special_prior <- prior$prior[is_special_prior]
    stop2("Prior ", collapse_comma(special_prior), " is used in an invalid ",
          "context. See ?set_prior for details on how to use special priors.")
  }

  px <- as.data.frame(px, stringsAsFactors = FALSE)
  upx <- unique(px)
  if (nrow(upx) > 1L) {
    # TODO: find a better solution to handle this case
    # can only happen for SD parameters of the same ID
    base_prior <- lb <- ub <- rep(NA, nrow(upx))
    base_bounds <- data.frame(lb = lb, ub = ub)
    for (i in seq_rows(upx)) {
      sub_upx <- lapply(upx[i, ], function(x) c(x, ""))
      sub_prior <- subset2(prior, ls = sub_upx)
      base_prior[i] <- stan_base_prior(sub_prior)
      base_bounds[i, ] <- stan_base_prior(sub_prior, col = c("lb", "ub"))
    }
    if (length(unique(base_prior)) > 1L) {
      # define prior for single coefficients manually
      # as there is not single base_prior anymore
      take_coef_prior <- nzchar(prior$coef) & !nzchar(prior$prior)
      prior_of_coefs <- prior[take_coef_prior, vars_prefix()]
      take_base_prior <- match_rows(prior_of_coefs, upx)
      prior$prior[take_coef_prior] <- base_prior[take_base_prior]
    }
    base_prior <- base_prior[1]
    if (nrow(unique(base_bounds)) > 1L) {
      stop2("Conflicting boundary information for ",
            "coefficients of class '", class, "'.")
    }
    base_bounds <- base_bounds[1, ]
  } else {
    base_prior <- stan_base_prior(prior)
    # select both bounds together so that they come from the same base prior
    base_bounds <- stan_base_prior(prior, col = c("lb", "ub"))
  }
  bound <- convert_bounds2stan(base_bounds)

  # generate stan prior statements
  out <- list()
  par <- paste0(prefix, class, suffix)
  has_constant_priors <- FALSE
  has_coef_prior <- any(with(prior, nzchar(coef) & nzchar(prior)))
  if (has_coef_prior || nzchar(dim) && length(coef)) {
    # priors on individual coefficients are also individually set
    # priors are always set on individual coefficients for arrays
    index_two_dims <- is.matrix(coef)
    coef <- as.matrix(coef)
    prior <- subset2(prior, coef = coef)
    estimated_coef_indices <- list()
    used_base_prior <- FALSE
    for (i in seq_rows(coef)) {
      for (j in seq_cols(coef)) {
        index <- i
        if (index_two_dims) {
          c(index) <- j
        }
        prior_ij <- subset2(prior, coef = coef[i, j])
        if (NROW(px) > 1L) {
          # disambiguate priors of coefficients with the same name
          # coming from different model components
          stopifnot(NROW(px) == NROW(coef))
          prior_ij <- subset2(prior_ij, ls = px[i, ])
        }
        # zero rows can happen if only global priors present
        stopifnot(nrow(prior_ij) <= 1L)
        coef_prior <- prior_ij$prior
        if (!isTRUE(nzchar(coef_prior))) {
          used_base_prior <- TRUE
          coef_prior <- base_prior
        }
        if (!stan_is_constant_prior(coef_prior)) {
          # all parameters with non-constant priors are estimated
          c(estimated_coef_indices) <- list(index)
        }
        if (nzchar(coef_prior)) {
          # implies a proper prior or constant
          if (type == coef_type && !nzchar(dim)) {
            # the single coefficient of that parameter equals the parameter
            stopifnot(all(index == 1L))
            par_ij <- par
          } else {
            par_ij <- paste0(par, collapse("[", index, "]"))
          }
          if (stan_is_constant_prior(coef_prior)) {
            coef_prior <- stan_constant_prior(
              coef_prior, par_ij, broadcast = broadcast
            )
            str_add(out$tpar_prior_const) <- paste0(coef_prior, ";\n")
          } else {
            coef_prior <- stan_target_prior(
              coef_prior, par_ij, broadcast = broadcast,
              bound = bound, resp = px$resp[1], normalize = normalize
            )
            str_add(out$tpar_prior) <- paste0(lpp(), coef_prior, ";\n")
          }
        }
      }
    }
    # the base prior may be improper flat in which no Stan code is added
    # but we still have estimated coefficients if the base prior is used
    has_estimated_priors <- isTRUE(nzchar(out$tpar_prior)) ||
      used_base_prior && !stan_is_constant_prior(base_prior)
    has_constant_priors <- isTRUE(nzchar(out$tpar_prior_const))
    if (has_estimated_priors && has_constant_priors) {
      # need to mix definition in the parameters and transformed parameters block
      if (!nzchar(coef_type)) {
        stop2("Can either estimate or fix all values of parameter '", par, "'.")
      }
      coef_type <- stan_type_add_bounds(coef_type, bound)
      for (i in seq_along(estimated_coef_indices)) {
        index <- estimated_coef_indices[[i]]
        iu <- paste0(index, collapse = "_")
        str_add(out$par) <- glue(
          "  {coef_type} par_{par}_{iu};\n"
        )
        ib <- collapse("[", index, "]")
        str_add(out$tpar_prior_const) <- cglue(
          "  {par}{ib} = par_{par}_{iu};\n"
        )
      }
    }
  } else if (nzchar(base_prior)) {
    # only a global prior is present and will be broadcasted
    ncoef <- length(coef)
    has_constant_priors <- stan_is_constant_prior(base_prior)
    if (has_constant_priors) {
      constant_base_prior <- stan_constant_prior(
        base_prior, par = par, ncoef = ncoef, broadcast = broadcast
      )
      str_add(out$tpar_prior_const) <- paste0(constant_base_prior, ";\n")
    } else {
      target_base_prior <- stan_target_prior(
        base_prior, par = par, ncoef = ncoef, bound = bound,
        broadcast = broadcast, resp = px$resp[1], normalize = normalize
      )
      str_add(out$tpar_prior) <- paste0(lpp(), target_base_prior, ";\n")
    }
  }

  if (nzchar(type)) {
    # only define the parameter here if type is non-empty
    type <- stan_adjust_par_type(type, base_prior)
    type <- stan_type_add_bounds(type, bound)
    comment <- stan_comment(comment)
    par_definition <- glue("  {type} {par}{dim};{comment}\n")
    if (has_constant_priors) {
      # parameter must be defined in the transformed parameters block
      str_add(out$tpar_def) <- par_definition
    } else {
      # parameter can be defined in the parameters block
      str_add(out$par) <- par_definition
    }
    if (nzchar(header_type)) {
      str_add(out$pll_args) <- glue(", {header_type} {par}")
    }
  } else {
    if (has_constant_priors) {
      stop2("Cannot fix parameter '", par, "' in this model.")
    }
  }
  has_improper_prior <- !is.null(out$par) && is.null(out$tpar_prior)
  if (prior_only && has_improper_prior) {
    stop2("Sampling from priors is not possible as ",
          "some parameters have no proper priors. ",
          "Error occurred for parameter '", par, "'.")
  }
  out
}

# extract base prior information for a given set of priors
# the base prior is the lowest level, non-flat, non-coefficient prior
# @param prior a brmsprior object
# @param col columns for which base prior information is to be found
# @param sel_prior optional brmsprior object to subset 'prior' before
#   finding the base prior
# @return the 'col' columns of the identified base prior
stan_base_prior <- function(prior, col = "prior", sel_prior = NULL, ...) {
  stopifnot(all(col %in% c("prior", "lb", "ub")))
  if (!is.null(sel_prior)) {
    # find the base prior using sel_prior for subsetting
    stopifnot(is.brmsprior(sel_prior))
    prior <- subset2(
      prior, class = sel_prior$class, group = c(sel_prior$group, ""),
      dpar = sel_prior$dpar, nlpar = sel_prior$nlpar, resp = sel_prior$resp,
      ...
    )
  } else {
    prior <- subset2(prior, ...)
  }
  stopifnot(length(unique(prior$class)) <= 1)
  # take all rows with non-zero entries on any of the chosen columns
  take <- !nzchar(prior$coef) & Reduce("|", lapply(prior[col], nzchar))
  prior <- prior[take, ]
  if (!NROW(prior)) {
    if (length(col) == 1L) {
      return("")
    } else {
      return(brmsprior()[, col])
    }
  }
  vars <- c("group", "nlpar", "dpar", "resp", "class")
  for (v in vars) {
    take <- nzchar(prior[[v]])
    if (any(take)) {
      prior <- prior[take, ]
    }
  }
  stopifnot(NROW(prior) == 1L)
  prior[, col]
}

# Stan prior in target += notation
# @param prior character string defining the prior
# @param par name of the parameter on which to set the prior
# @param ncoef number of coefficients in the parameter
# @param bound bounds of the parameter in Stan language
# @param broadcast Stan type to which the prior should be broadcasted
# @param name of the response variable
# @return a character string defining the prior in Stan language
stan_target_prior <- function(prior, par, ncoef = 0, broadcast = "vector",
                              bound = "", resp = "", normalize = TRUE) {
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
    prior_args[i] <- sub(")$", "", prior_args[i])
  }
  if (broadcast == "matrix" && ncoef > 0) {
    # apply a scalar prior to all elements of a matrix
    par <- glue("to_vector({par})")
  }

  if (nzchar(prior_args)) {
    str_add(prior_args, start = TRUE) <- " | "
  }
  lpdf <- stan_lpdf_name(normalize)
  out <- glue("{prior_name}_{lpdf}({par}{prior_args})")
  par_class <- unique(get_matches("^[^_]+", par))
  par_bound <- convert_stan2bounds(bound)
  prior_bound <- prior_bounds(prior_name)
  trunc_lb <- is.character(par_bound$lb) || par_bound$lb > prior_bound$lb
  trunc_ub <- is.character(par_bound$ub) || par_bound$ub < prior_bound$ub
  if (normalize) {
    # obtain correct normalization constants for truncated priors
    if (trunc_lb || trunc_ub) {
      wsp <- wsp(nsp = 4)
      # scalar parameters are of length 1 but have no coefficients
      ncoef <- max(1, ncoef)
      if (trunc_lb && !trunc_ub) {
        str_add(out) <- glue(
          "\n{wsp}- {ncoef} * {prior_name}_lccdf({par_bound$lb}{prior_args})"
        )
      } else if (!trunc_lb && trunc_ub) {
        str_add(out) <- glue(
          "\n{wsp}- {ncoef} * {prior_name}_lcdf({par_bound$ub}{prior_args})"
        )
      } else if (trunc_lb && trunc_ub) {
        str_add(out) <- glue(
          "\n{wsp}- {ncoef} * log_diff_exp(",
          "{prior_name}_lcdf({par_bound$ub}{prior_args}), ",
          "{prior_name}_lcdf({par_bound$lb}{prior_args}))"
        )
      }
    }
  }
  out
}

# fix parameters to constants in Stan language
# @param prior character string defining the prior
# @param par name of the parameter on which to set the prior
# @param ncoef number of coefficients in the parameter
# @param broadcast Stan type to which the prior should be broadcasted
# @return a character string defining the prior in Stan language
stan_constant_prior <- function(prior, par, ncoef = 0, broadcast = "vector") {
  stopifnot(grepl("^constant\\(", prior))
  prior_args <- gsub("(^constant\\()|(\\)$)", "", prior)
  if (broadcast == "vector") {
    if (ncoef > 0) {
      # broadcast the scalar prior on the whole parameter vector
      prior_args <- glue("rep_vector({prior_args}, rows({par}))")
    }
    # no action required for individual coefficients of vectors
  } else if (broadcast == "matrix") {
    if (ncoef > 0) {
      # broadcast the scalar prior on the whole parameter matrix
      prior_args <- glue("rep_matrix({prior_args}, rows({par}), cols({par}))")
    } else {
      # single coefficient is a row in the parameter matrix
      prior_args <- glue("rep_row_vector({prior_args}, cols({par}))")
    }
  }
  glue("  {par} = {prior_args}")
}

# Stan code for global parameters of special shrinkage priors
stan_special_prior <- function(bterms, out, data, prior, ranef, normalize, ...) {
  stopifnot(is.list(out))
  tp <- tp()
  lpp <- lpp()
  lpdf <- stan_lpdf_name(normalize)
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  if (!has_special_prior(prior, px)) {
    return(out)
  }
  special <- get_special_prior(prior, px, main = TRUE)
  str_add(out$data) <- glue(
    "  int<lower=1> Kscales{p};  // number of local scale parameters\n"
  )
  if (special$name == "horseshoe") {
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
      "  real<lower=0> hs_global{p};  // global shrinkage parameter\n",
      "  real<lower=0> hs_slab{p};  // slab regularization parameter\n",
      "  vector<lower=0>[Kscales{p}] hs_local{p};  // local parameters for the horseshoe prior\n"
    )
    hs_scale_global <- glue("hs_scale_global{p}")
    if (isTRUE(special$autoscale)) {
      str_add(hs_scale_global) <- glue(" * sigma{usc(px$resp)}")
    }
    str_add(out$tpar_prior) <- glue(
      "{lpp}student_t_{lpdf}(hs_global{p} | hs_df_global{p}, 0, {hs_scale_global})",
      str_if(normalize, "\n    - 1 * log(0.5)"), ";\n",
      "{lpp}inv_gamma_{lpdf}(hs_slab{p} | 0.5 * hs_df_slab{p}, 0.5 * hs_df_slab{p});\n"
    )
    str_add(out$tpar_def) <- glue(
      "  vector<lower=0>[Kscales{p}] scales{p};  // local horseshoe scale parameters\n"
    )
    str_add(out$tpar_comp) <- glue(
      "  // compute horseshoe scale parameters\n",
      "  scales{p} = scales_horseshoe(hs_local{p}, hs_global{p}, hs_scale_slab{p}^2 * hs_slab{p});\n"
    )
    str_add(out$model_prior) <- glue(
      "{tp}student_t_{lpdf}(hs_local{p} | hs_df{p}, 0, 1)",
      str_if(normalize, "\n    - rows(hs_local{p}) * log(0.5)"), ";\n"
    )
  } else if (special$name == "R2D2") {
    str_add(out$data) <- glue(
      "  // data for the R2D2 prior\n",
      "  real<lower=0> R2D2_mean_R2{p};  // mean of the R2 prior\n",
      "  real<lower=0> R2D2_prec_R2{p};  // precision of the R2 prior\n",
      "  // concentration vector of the D2 prior\n",
      "  vector<lower=0>[Kscales{p}] R2D2_cons_D2{p};\n"
    )
    str_add(out$par) <- glue(
      "  // parameters of the R2D2 prior\n",
      "  real<lower=0,upper=1> R2D2_R2{p};\n",
      "  simplex[Kscales{p}] R2D2_phi{p};\n"
    )
    var_mult <- ""
    if (isTRUE(special$autoscale)) {
      var_mult <- glue("sigma{usc(px$resp)}^2 * ")
    }
    str_add(out$tpar_def) <- glue(
      "  real R2D2_tau2{p};  // global R2D2 scale parameter\n",
      "  vector<lower=0>[Kscales{p}] scales{p};  // local R2D2 scale parameters\n"
    )
    str_add(out$tpar_comp) <- glue(
      "  // compute R2D2 scale parameters\n",
      "  R2D2_tau2{p} = {var_mult}R2D2_R2{p} / (1 - R2D2_R2{p});\n",
      "  scales{p} = scales_R2D2(R2D2_phi{p}, R2D2_tau2{p});\n"
    )
    str_add(out$tpar_prior) <- glue(
      "{lpp}beta_{lpdf}(R2D2_R2{p} | R2D2_mean_R2{p} * R2D2_prec_R2{p}, ",
      "(1 - R2D2_mean_R2{p}) * R2D2_prec_R2{p});\n"
    )
    str_add(out$model_prior) <- glue(
      "{tp}dirichlet_{lpdf}(R2D2_phi{p} | R2D2_cons_D2{p});\n"
    )
  }

  if (has_special_prior(prior, px, class = "sd")) {
    # this has to be done here rather than in stan_re()
    # because the latter is not local to a linear predictor
    ids <- unique(subset2(ranef, ls = px)$id)
    str_add(out$prior_global_scales) <- glue(" sd_{ids}")
    str_add(out$prior_global_lengths) <- glue(" M_{ids}")
  }
  # split up scales into subsets belonging to different parameter classes
  # this connects the global to the local priors
  scales <- strsplit(trimws(out$prior_global_scales), " ")[[1]]
  lengths <- strsplit(trimws(out$prior_global_lengths), " ")[[1]]
  out$prior_global_scales <- out$prior_global_lengths <- NULL

  lengths <- c("1", lengths)
  for (i in seq_along(scales)) {
    lower <- paste0(lengths[1:i], collapse = "+")
    upper <- paste0(lengths[2:(i+1)], collapse = "+")
    # some scale parameters are a scalar not a vector
    bracket1 <- str_if(lengths[i+1] == "1", "[1]")
    str_add(out$tpar_comp) <- glue(
      "  {scales[i]} = scales{p}[({lower}):({upper})]{bracket1};\n"
    )
  }
  out
}

# Stan code of normal priors on regression coefficients
# in non-centered parameterization
# @param class name of the coefficient class
# @param suffix shared suffix of the involved variables
# @param suffix_class extra suffix of the class
# @param suffix_K extra suffix of K (number of coefficients)
stan_prior_non_centered <- function(class = "b", suffix = "", suffix_class = "",
                                    suffix_K = "", normalize = TRUE) {
  out <- list()
  csfx <- glue("{class}{suffix}")
  csfx2 <- glue("{class}{suffix_class}{suffix}")
  Ksfx <- glue("K{suffix_K}{suffix}")
  lpdf <- stan_lpdf_name(normalize)
  str_add(out$tpar_def) <- glue(
    "  vector[{Ksfx}] {csfx2};  // scaled coefficients\n"
  )
  str_add(out$par) <- glue(
    "  vector[{Ksfx}] z{csfx};  // unscaled coefficients\n"
  )
  str_add(out$tpar_def) <- glue(
    "  vector[{Ksfx}] sd{csfx};  // SDs of the coefficients\n"
  )
  str_add(out$tpar_special_prior) <- glue(
    "  {csfx2} = z{csfx} .* sd{csfx};  // scale coefficients\n"
  )
  str_add(out$model_prior) <- glue(
    "{tp()}std_normal_{lpdf}(z{csfx});\n"
  )
  str_add(out$prior_global_scales) <- glue(" sd{csfx}")
  str_add(out$prior_global_lengths) <- glue(" {Ksfx}")
  str_add(out$pll_args) <- glue(", vector {csfx2}")
  out
}

# combine unchecked priors for use in Stan
# @param prior a brmsprior object
# @return a single character string in Stan language
stan_unchecked_prior <- function(prior) {
  stopifnot(is.brmsprior(prior))
  if (all(nzchar(prior$class))) {
    return("")
  }
  prior <- subset2(prior, class = "")
  collapse("  ", prior$prior, ";\n")
}

# Stan code to sample separately from priors
# @param tpar_prior character string taken from stan_prior that contains
#   all priors that can potentially be sampled from separately
# @param par_declars the parameters block of the Stan code
#   required to extract boundaries
# @param gen_quantities Stan code from the generated quantities block
# @param special_prior a list of values pertaining to special priors
#   such as horseshoe or lasso
# @param sample_prior take draws from priors?
stan_rngprior <- function(tpar_prior, par_declars, gen_quantities,
                          special_prior, sample_prior = "yes") {
  if (!is_equal(sample_prior, "yes")) {
    return(list())
  }
  tpar_prior <- strsplit(gsub(" |\\n", "", tpar_prior), ";")[[1]]
  # D will contain all relevant information about the priors
  D <- data.frame(prior = tpar_prior[nzchar(tpar_prior)])
  pars_regex <- "(?<=(_lpdf\\())[^|]+"
  D$par <- get_matches(pars_regex, D$prior, perl = TRUE, first = TRUE)
  # 'std_normal' has no '|' and thus the above regex matches too much
  np <- !grepl("\\|", D$prior)
  np_regex <- ".+(?=\\)$)"
  D$par[np] <- get_matches(np_regex, D$par[np], perl = TRUE, first = TRUE)
  # 'to_vector' should be removed from the parameter names
  has_tv <- grepl("^to_vector\\(", D$par)
  tv_regex <- "(^to_vector\\()|(\\)(?=((\\[[[:digit:]]+\\])?)$))"
  D$par[has_tv] <- gsub(tv_regex, "", D$par[has_tv], perl = TRUE)
  # do not sample from some auxiliary parameters
  excl_regex <- c("tmp")
  excl_regex <- paste0("(", excl_regex, ")", collapse = "|")
  excl_regex <- paste0("^(", excl_regex, ")(_|$)")
  D <- D[!grepl(excl_regex, D$par), ]
  if (!NROW(D)) return(list())

  # rename parameters containing indices
  has_ind <- grepl("\\[[[:digit:]]+\\]", D$par)
  D$par[has_ind] <- ulapply(D$par[has_ind], function(par) {
    ind_regex <- "(?<=\\[)[[:digit:]]+(?=\\])"
    ind <- get_matches(ind_regex, par, perl = TRUE)
    gsub("\\[[[:digit:]]+\\]", paste0("__", ind), par)
  })
  # cannot handle priors on variable transformations
  D <- D[D$par %in% stan_all_vars(D$par), ]
  if (!NROW(D)) return(list())

  class_old <- c("^L_", "^Lrescor")
  class_new <- c("cor_", "rescor")
  D$par <- rename(D$par, class_old, class_new, fixed = FALSE)
  dis_regex <- "(?<=lprior\\+=)[^\\(]+(?=_lpdf\\()"
  D$dist <- get_matches(dis_regex, D$prior, perl = TRUE, first = TRUE)
  D$dist <- sub("corr_cholesky$", "corr", D$dist)
  args_regex <- "(?<=\\|)[^$\\|]+(?=\\)($|-))"
  D$args <- get_matches(args_regex, D$prior, perl = TRUE, first = TRUE)
  # 'std_normal_rng' does not exist in Stan
  has_std_normal <- D$dist == "std_normal"
  D$dist[has_std_normal] <- "normal"
  D$args[has_std_normal] <- "0,1"

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
    if (grepl("^((simo_)|(theta)|(R2D2_phi))", all_pars[i])) {
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
  str_add(out$gen_def) <- "  // additionally sample draws from priors\n"
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

# are multiple base priors supplied?
# px list of class, dpar, etc. elements used to infer parameter suffixes
stan_has_multiple_base_priors <- function(px) {
  px <- as.data.frame(px, stringsAsFactors = FALSE)
  nrow(unique(px)) > 1L
}

# check if any constant priors are present
# @param prior a vector of character strings
stan_is_constant_prior <- function(prior) {
  grepl("^constant\\(", prior)
}

# extract Stan boundaries expression from a string
stan_extract_bounds <- function(x, bound = c("lower", "upper")) {
  bound <- match.arg(bound)
  x <- rm_wsp(x)
  regex <- glue("(?<={bound}=)[^,>]*")
  get_matches(regex, x, perl = TRUE, first = TRUE)
}

# choose the right suffix for Stan probability densities
stan_lpdf_name <- function(normalize, int = FALSE) {
  if (normalize) {
    out <- ifelse(int, "lpmf", "lpdf")
  } else {
    out <- ifelse(int, "lupmf", "lupdf")
  }
  out
}

# add bounds to a Stan type specification which may include dimensions
stan_type_add_bounds <- function(type, bound) {
  regex_dim <- "\\[.*$"
  type_type <- sub(regex_dim, "", type)
  type_dim <- get_matches(regex_dim, type, first = TRUE)
  glue("{type_type}{bound}{type_dim}")
}

# adjust the type of a parameter based on the assigned prior
stan_adjust_par_type <- function(type, prior) {
  # TODO: add support for more type-prior combinations?
  combs <- data.frame(
    type = "vector",
    prior = "dirichlet",
    new_type = "simplex"
  )
  for (i in seq_rows(combs)) {
    regex_type <- paste0("^", combs$type[i], "\\[?")
    regex_prior <- paste0("^", combs$prior[i], "\\(")
    if (grepl(regex_type, type) && grepl(regex_prior, prior)) {
      brackets <- get_matches("\\[.*\\]$", type, first = TRUE)
      type <- paste0(combs$new_type[i], brackets)
      break
    }
  }
  type
}

# stops if a prior bound is given
stopif_prior_bound <- function(prior, class, ...) {
  lb <- stan_base_prior(prior, "lb", class = class, ...)
  ub <- stan_base_prior(prior, "ub", class = class, ...)
  if (nzchar(lb) || nzchar(ub)) {
    stop2("Cannot add bounds to class '", class, "' for this prior.")
  }
  return(invisible(NULL))
}

# lprior plus equal
lpp <- function(wsp = 2) {
  wsp <- collapse(rep(" ", wsp))
  paste0(wsp, "lprior += ")
}
