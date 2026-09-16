# Stan generation for physical sum-to-zero group effects

# Internal-output metadata is kept beside the generators that own these
# coordinates; ordinary save_pars handling only asks for the relevant names.
re_s2z_internal_fe <- function(plan) {
  p <- plan$p
  c(
    paste0(c("theta_s2z", "theta_s2z_active", "fixed_s2z",
             "zfixed_s2z", "sdfixed_s2z"), p),
    paste0("par_fixed_s2z", p, "_", seq_along(plan$inactive_q)),
    paste0("udf_b_s2z", p, "_", seq_along(plan$qnames))
  )
}

re_s2z_internal_re <- function(r) {
  classes <- c(
    "z_s2z", "z_mean_s2z", "r_s2z", "H_s2z",
    "q_explicit_s2z", "prior_mean_s2z", "prior_prec_s2z", "prior_scale_s2z",
    "W_matheron_s2z", "sqrt_W_matheron_s2z", "L_W_matheron_s2z",
    "theta_white_matheron_s2z", "group_scale_s2z", "group_prec_s2z",
    "group_info_s2z", "L_Sigma_s2z", "Q_Sigma_s2z", "P_group_s2z",
    "h_group_s2z", "P_s2z", "L_P_s2z", "H_joint_s2z", "h_joint_s2z",
    "D_s2z", "sqrt_D_s2z", "D_diag_s2z", "intercept_map_s2z",
    "rank1_info_s2z", "group_quad_s2z", "joint_quad_s2z", "mhat_s2z",
    "qhat_s2z", "white_s2z", "mean_r_s2z", "q_recovered_s2z"
  )
  rp <- usc(combine_prefix(r))
  c(paste0(classes, "_", r$id[1L]), paste0("r_s2z_", r$id, rp, "_", r$cn))
}

# Resolve each predictor's priors once, then select numerical work independently
# of public parameter and predictor ownership.
stan_re_s2z_systems <- function(bframe, prior) {
  lapply(Filter(has_re_s2z, all_bframel(bframe)), function(bfl) {
    plan <- re_s2z_plan(bfl, prior = prior)
    infos <- re_s2z_plan_infos(plan)
    system <- list(plan = plan, set_id = plan$set_id, ids = plan$ids,
                   bfl = bfl, infos = infos)
    if (.stan_re_s2z_uses_explicit_mean(system)) {
      kernel <- "explicit"
    } else if (length(plan$blocks) == 1L) {
      r <- plan$blocks[[1L]]$r
      kernel <- if (nrow(r) == 1L) "scalar" else
        if (isTRUE(r$cor[1])) "dense" else "independent"
    } else if (.stan_re_s2z_joint_uses_matheron(system)) {
      kernel <- "matheron"
    } else {
      kernel <- "joint"
    }
    system$kernel <- kernel
    system$small_matrix <- kernel %in% c("explicit", "dense", "joint") ||
      kernel == "matheron" && (
        length(.stan_re_s2z_joint_matheron_info(system)$P) > 1L ||
        any(vapply(plan$blocks, function(block) {
          nrow(block$r) > 1L && isTRUE(block$r$cor[1])
        }, logical(1)))
      )
    system
  })
}

stan_re_s2z_by_id <- function(systems) {
  out <- list()
  for (system in systems) {
    for (id in system$ids) out[[as.character(id)]] <- system
  }
  out
}

# The caller has already emitted the existing group data and scale priors.
stan_re_s2z_block <- function(id, system, bframe, prior, threads, normalize,
                              out = list()) {
  stopifnot(!is.null(system), id %in% system$ids)
  if (length(system$ids) > 1L) {
    return(.stan_re_s2z_joint_block(
      id, set = system, bframe = bframe, prior = prior, threads = threads,
      normalize = normalize, out = out
    ))
  }
  if (system$kernel == "explicit") {
    out <- .stan_re_s2z_explicit_block(
      id, set = system, prior = prior, normalize = normalize, out = out
    )
    return(.stan_re_s2z_explicit_system(system, normalize, out = out))
  }
  .stan_re_s2z(id, system = system, bframe = bframe, prior = prior,
               threads = threads, normalize = normalize, out = out)
}

# Singleton kernels emit their system work together with their one block to
# retain declaration and RNG ordering. Multiblock work follows all its blocks.
stan_re_s2z_system <- function(system, prior, normalize, ...) {
  if (length(system$ids) == 1L) return(list())
  .stan_re_s2z_joint(system, prior = prior, normalize = normalize, ...)
}

# The finite-population coefficient vector used by either likelihood path.
stan_re_s2z_coef <- function(bframe, prefix = NULL) {
  plan <- re_s2z_plan(bframe)
  stopifnot(!is.null(plan))
  p <- prefix %||% plan$p
  if (plan$center) {
    glue("tail(theta_s2z{p}, {length(plan$fixef)})")
  } else {
    glue("theta_s2z{p}")
  }
}

# Reconstruct ordinary group effects from physical deviations and a kernel's
# recovered mean. Splitting matrix/column work preserves existing output order.
stan_re_s2z_public_re_comp <- function(r, scalar_mean = FALSE, indent = "  ",
                                      part = "all") {
  id <- r$id[1L]
  J <- seq_rows(r)
  idp <- paste0(r$id, usc(combine_prefix(check_prefix(r))))
  r_public <- glue("r_{idp}_{r$cn}")
  r_s2z <- glue("r_s2z_{idp}_{r$cn}")
  is_cor <- nrow(r) > 1L && isTRUE(r$cor[1])
  if (is_cor) {
    out <- ""
    if (part %in% c("all", "matrix")) {
      str_add(out) <- glue(
        "{indent}r_{id} = r_s2z_{id};\n",
        "{indent}for (j in 1:N_{id}) r_{id}[j] += mean_r_s2z_{id}';\n"
      )
    }
    if (part %in% c("all", "columns")) {
      str_add(out) <- cglue("{indent}{r_public} = r_{id}[, {J}];\n")
    }
    return(out)
  }
  mean_index <- if (scalar_mean) "" else glue("[{J}]")
  cglue("{indent}{r_public} = {r_s2z} + mean_r_s2z_{id}{mean_index};\n")
}

# Public population coefficients have the same names and design-centering map
# for every numerical kernel. The kernels only supply the recovered vector.
stan_re_s2z_public_fe_def <- function(info) {
  p <- info$p
  q <- length(info$qnames)
  if (info$center) {
    glue(
      "  real Intercept{p};\n",
      str_if(length(info$fixef), glue("  vector[Kc{p}] b{p};\n")),
      "  real b{p}_Intercept;\n"
    )
  } else {
    glue("  vector[{q}] b{p};\n")
  }
}

stan_re_s2z_public_fe_comp <- function(info, set_id) {
  p <- info$p
  if (info$center) {
    glue(
      "  Intercept{p} = q_recovered_s2z_{set_id}[1];\n",
      str_if(
        length(info$fixef),
        glue(
          "  b{p} = tail(q_recovered_s2z_{set_id}, Kc{p});\n",
          "  b{p}_Intercept = Intercept{p} - dot_product(",
          "means_X{p}, b{p});\n"
        )
      ),
      str_if(!length(info$fixef), glue("  b{p}_Intercept = Intercept{p};\n"))
    )
  } else {
    glue("  b{p} = q_recovered_s2z_{set_id};\n")
  }
}

# Population coordinates for one physical S2Z predictor. Ordinary fixed-effect
# generation only dispatches here; the existing design-data setup is shared.
stan_fe_s2z <- function(bframe, prior, normalize) {
  out <- list()
  fixef <- bframe$frame$fe$vars_stan
  ct <- str_if(bframe$frame$fe$center, "c")
  px <- check_prefix(bframe)
  p <- usc(combine_prefix(px))
  plan <- re_s2z_plan(bframe)
  K_s2z <- length(plan$qnames)
  active_s2z <- plan$active_q
  inactive_s2z <- plan$inactive_q
  special_b_s2z <- has_special_prior(prior, bframe, class = "b")
  str_add(out$par) <- glue(
    "  vector[{length(active_s2z)}] theta_s2z_active{p};",
    "  // S2Z-active finite-population coefficients\n"
  )
  if (length(inactive_s2z)) {
    inactive_names <- bframe$frame$fe$vars[inactive_s2z]
    stopifnot(!"Intercept" %in% inactive_names)
    if (special_b_s2z) {
      # Validation permits a class-b shrinkage prior only when every
      # population slope is fixed-only. In that case the ordinary brms
      # non-centered representation remains exact, and Kscales is still
      # Kc (or K for an uncentered predictor).
      stopifnot(length(inactive_names) == length(fixef))
      inactive_prior <- stan_prior_non_centered(
        class = "fixed_s2z", suffix = p, suffix_K = ct,
        normalize = normalize
      )
      # The likelihood receives the assembled theta vector instead.
      inactive_prior$pll_args <- NULL
      str_add_list(out) <- inactive_prior
    } else {
      inactive_priors <- stan_re_s2z_inactive_priors(
        prior, coef = inactive_names, fixef = fixef, px = px
      )
      inactive_prior <- stan_prior(
        inactive_priors, class = "b", coef = inactive_names,
        type = glue("vector[{length(inactive_s2z)}]"),
        suffix = glue("_s2z_inactive{p}"), px = px,
        comment = "S2Z-inactive regression coefficients",
        normalize = normalize
      )
      # Keep this implementation detail out of public b_* discovery while
      # retaining stan_prior's coefficient, bound, and constant handling.
      inactive_prior <- lapply(inactive_prior, function(code) {
        gsub(
          glue("b_s2z_inactive{p}"), glue("fixed_s2z{p}"), code,
          fixed = TRUE
        )
      })
      str_add_list(out) <- inactive_prior
    }
  }
  str_add(out$tpar_def) <- glue(
    "  vector[{K_s2z}] theta_s2z{p};",
    "  // assembled finite-population coefficients\n"
  )
  theta_comp <- if (special_b_s2z) "tpar_special_prior" else "tpar_comp"
  for (a in seq_along(active_s2z)) {
    str_add(out[[theta_comp]]) <- glue(
      "  theta_s2z{p}[{active_s2z[a]}] = ",
      "theta_s2z_active{p}[{a}];\n"
    )
  }
  for (a in seq_along(inactive_s2z)) {
    str_add(out[[theta_comp]]) <- glue(
      "  theta_s2z{p}[{inactive_s2z[a]}] = ",
      "fixed_s2z{p}[{a}];\n"
    )
  }
  str_add(out$pll_args) <- glue(", vector theta_s2z{p}")
  out
}

# Keep global scalar-distribution priors on the original coefficient indices
# when the S2Z split shortens the fixed-only vector. Active-coordinate prior
# validation is unchanged; explicit coefficient priors retain their arguments.
stan_re_s2z_inactive_priors <- function(prior, coef, fixef, px) {
  local_prior <- subset2(
    prior, class = "b", coef = c(coef, ""), group = "", ls = px
  )
  base <- stan_base_prior(local_prior)
  name <- trimws(sub("\\(.*$", "", trimws(base)))
  if (!name %in% c("normal", "student_t", "cauchy", "logistic")) {
    return(prior)
  }
  # Split only top-level commas, retaining nested Stan expressions verbatim.
  args <- sub("^[^(]*\\(", "", base)
  args <- sub("\\)[[:space:]]*$", "", args)
  chars <- strsplit(args, "", fixed = TRUE)[[1L]]
  depth <- 0L
  quoted <- escaped <- FALSE
  commas <- integer()
  for (i in seq_along(chars)) {
    ch <- chars[i]
    if (quoted) {
      if (escaped) {
        escaped <- FALSE
      } else if (ch == "\\") {
        escaped <- TRUE
      } else if (ch == '"') {
        quoted <- FALSE
      }
    } else if (ch == '"') {
      quoted <- TRUE
    } else if (ch %in% c("(", "[", "{")) {
      depth <- depth + 1L
    } else if (ch %in% c(")", "]", "}")) {
      depth <- depth - 1L
    } else if (ch == "," && depth == 0L) {
      commas <- c(commas, i)
    }
  }
  args <- trimws(substring(
    args, c(1L, commas + 1L), c(commas - 1L, nchar(args))
  ))
  symbolic <- !is.finite(suppressWarnings(as.numeric(args)))
  if (!any(symbolic)) {
    return(prior)
  }
  indices <- match(coef, fixef)
  stopifnot(!anyNA(indices))
  base_tag <- stan_base_prior(local_prior, col = "tag")
  for (i in seq_along(coef)) {
    rows <- which(find_rows(prior, class = "b", coef = coef[i],
                            group = "", ls = px))
    if (any(nzchar(prior$prior[rows]))) {
      next
    }
    indexed <- args
    indexed[symbolic] <- glue(
      "s2z_prior_coordinate_brms({args[symbolic]}, ",
      "{indices[i]}, {length(fixef)})"
    )
    value <- paste0(name, "(", paste(indexed, collapse = ", "), ")")
    if (length(rows)) {
      prior$prior[rows] <- value
      prior$tag[rows] <- base_tag
    } else {
      prior <- rbind(prior, set_prior(
        value, class = "b", coef = coef[i], resp = px$resp,
        dpar = px$dpar, nlpar = px$nlpar, tag = base_tag
      ))
    }
  }
  prior
}

# The separated Gaussian mean/contrast law permits a Matheron update in the
# population-coordinate dimension instead of a complete square in the sum of
# all block dimensions. Conditional Student-t population priors remain
# eligible because their scale-mixture variables make them proper Gaussians.
# Student-t group effects do not separate their omitted means from their
# contrasts and therefore use the general joint system.
.stan_re_s2z_joint_matheron_info <- function(set) {
  infos <- set$infos
  if (length(infos) <= 1L) {
    return(NULL)
  }
  separated <- all(vapply(infos, function(info) {
    all(info$r$dist == "gaussian")
  }, logical(1)))
  if (!separated) {
    return(NULL)
  }
  # Only rows reached by an omitted group mean can couple blocks. A flat
  # population prior removes its row from that conditioning system, while a
  # proper row outside the group design remains an ordinary independent prior.
  active <- unique(unlist(lapply(infos, function(info) {
    rows <- info$match_q
    if (info$center && any(info$r$coef != "Intercept")) {
      rows <- c(rows, 1L)
    }
    rows
  }), use.names = FALSE))
  proper <- which(vapply(infos[[1L]]$prior, function(spec) {
    spec$dist %in% c("normal", "student")
  }, logical(1)))
  P <- sort(intersect(active, proper))
  total_M <- sum(vapply(infos, function(info) nrow(info$r), integer(1)))
  # The general complete square is already preferable at equal or lower
  # dimension. Strict inequality makes this dispatch a genuine fast path.
  if (length(P) >= total_M) {
    return(NULL)
  }
  list(P = P, inactive = setdiff(proper, P), total_M = total_M)
}

.stan_re_s2z_joint_uses_matheron <- function(set) {
  !is.null(.stan_re_s2z_joint_matheron_info(set))
}

# A logistic population prior cannot be integrated by the conditional
# Gaussian omitted-mean kernels. Such systems retain the omitted means
# explicitly in standardized coordinates and score the exact requested prior.
.stan_re_s2z_uses_explicit_mean <- function(set) {
  any(vapply(set$infos, function(info) {
    !is.null(info$prior) && any(vapply(
      info$prior, function(spec) identical(spec$dist, "logistic"), logical(1)
    ))
  }, logical(1)))
}

# Exact scalar population-prior statement for the explicit-mean path.
stan_re_s2z_prior_target <- function(spec, par, normalize) {
  stopifnot(is.list(spec), is.character(par), length(par) == 1L)
  if (identical(spec$dist, "flat")) {
    return("")
  }
  lpdf <- stan_lpdf_name(normalize)
  location <- stan_s2z_number(spec$location)
  scale <- stan_s2z_number(spec$scale)
  if (identical(spec$dist, "normal")) {
    glue("  lprior += normal_{lpdf}({par} | {location}, {scale});\n")
  } else if (identical(spec$dist, "student")) {
    df <- stan_s2z_number(spec$df)
    glue(
      "  lprior += student_t_{lpdf}({par} | {df}, {location}, {scale});\n"
    )
  } else if (identical(spec$dist, "logistic")) {
    glue("  lprior += logistic_{lpdf}({par} | {location}, {scale});\n")
  } else {
    stop2("Internal error: unsupported explicit S2Z population prior.")
  }
}

# One block in an exact explicit-mean system. The sampled mean coordinate is
# v = sqrt(N) L^-1 m, so m = L v / sqrt(N). With orthonormal S2Z contrasts,
# the Jacobian back to conventional group effects is |L|.
.stan_re_s2z_explicit_block <- function(id, set, prior, normalize,
                                        out = list()) {
  if (is.null(out[["tpar_prior"]])) {
    out[["tpar_prior"]] <- ""
  }
  take_info <- match(id, vapply(set$infos, `[[`, numeric(1), "id"))
  stopifnot(!is.na(take_info))
  info <- set$infos[[take_info]]
  r <- info$r
  q <- length(info$qnames)
  M <- nrow(r)
  J <- seq_rows(r)
  p <- info$p
  px <- check_prefix(r)
  idp <- paste0(r$id, usc(combine_prefix(px)))
  r_s2z <- glue("r_s2z_{idp}_{r$cn}")
  is_cor <- M > 1L && isTRUE(r$cor[1])
  is_student <- identical(r$dist[1], "student")

  if (is_cor) {
    str_add(out$data) <- glue(
      "  int<lower=1> NC_{id};  // number of group-level correlations\n"
    )
    str_add_list(out) <- stan_prior(
      prior, class = "L", group = r$group[1], suffix = usc(id),
      type = glue("cholesky_factor_corr[M_{id}]"),
      comment = "cholesky factor of correlation matrix",
      normalize = normalize
    )
  }
  if (identical(id, set$set_id)) {
    str_add(out$fun) <- "  #include 'fun_sum_to_zero.stan'\n"
  }

  str_add(out$par) <- glue(
    "  vector[M_{id} * (N_{id} - 1)] z_s2z_{id};",
    "  // physical orthonormal S2Z coordinates\n",
    "  vector[M_{id}] z_mean_s2z_{id};",
    "  // standardized omitted group mean\n"
  )
  str_add(out$tpar_def) <- glue(
    "  // exact explicit-mean S2Z block {id}\n",
    "  matrix[N_{id}, M_{id}] r_s2z_{id};\n",
    "  matrix[M_{id}, M_{id}] L_Sigma_s2z_{id};\n",
    "  vector[M_{id}] mean_r_s2z_{id};\n",
    "  real<lower=0> group_quad_s2z_{id};\n",
    str_if(
      is_student,
      glue(
        "  vector<lower=0>[N_{id}] group_scale_s2z_{id};\n",
        "  vector<lower=0>[N_{id}] group_prec_s2z_{id};\n"
      )
    ),
    cglue("  vector[N_{id}] {r_s2z};\n")
  )

  # H is fixed by the formula and centered predictor means. Build it once in
  # transformed data so it never enters the reverse-mode expression graph.
  str_add(out$tdata_def) <- glue(
    "  matrix[{q}, M_{id}] H_s2z_{id};\n"
  )
  str_add(out$tdata_comp) <- glue(
    "  H_s2z_{id} = rep_matrix(0.0, {q}, M_{id});\n"
  )
  for (j in seq_len(M)) {
    qi <- info$match_q[j]
    str_add(out$tdata_comp) <- glue(
      "  H_s2z_{id}[{qi}, {j}] = 1.0;\n"
    )
    if (info$center && info$r$coef[j] != "Intercept") {
      str_add(out$tdata_comp) <- glue(
        "  H_s2z_{id}[1, {j}] = means_X{p}[{qi - 1L}];\n"
      )
    }
  }
  if (is_cor) {
    str_add(out$tpar_comp) <- glue(
      "  L_Sigma_s2z_{id} = diag_pre_multiply(sd_{id}, L_{id});\n"
    )
  } else {
    str_add(out$tpar_comp) <- glue(
      "  L_Sigma_s2z_{id} = diag_matrix(sd_{id});\n"
    )
  }
  str_add(out$tpar_comp) <- glue(
    "  for (k in 1:M_{id}) {{\n",
    "    r_s2z_{id}[, k] = sum_to_zero_constrain_brms(",
    "segment(z_s2z_{id}, (k - 1) * (N_{id} - 1) + 1, ",
    "N_{id} - 1));\n",
    "  }}\n",
    "  mean_r_s2z_{id} = L_Sigma_s2z_{id} * z_mean_s2z_{id} / ",
    "sqrt(1.0 * N_{id});\n"
  )
  if (is_student) {
    tr <- subset_reframe_dist(r, "student")
    g <- usc(tr$ggn[1])
    str_add(out$tpar_comp) <- glue(
      "  group_scale_s2z_{id} = dfm{g};\n",
      "  group_prec_s2z_{id} = inv_square(group_scale_s2z_{id});\n"
    )
  }
  explicit_group_quad_code <- if (is_student) {
    glue(
      "    vector[M_{id}] mean_white_s2z = z_mean_s2z_{id} / ",
      "sqrt(1.0 * N_{id});\n",
      "    group_quad_s2z_{id} = 0.0;\n",
      "    for (j in 1:N_{id}) {{\n",
      "      vector[M_{id}] white_level_s2z = ",
      "white_group_s2z[, j] + mean_white_s2z;\n",
      "      group_quad_s2z_{id} += group_prec_s2z_{id}[j] * ",
      "dot_self(white_level_s2z);\n",
      "    }}\n"
    )
  } else {
    glue(
      "    group_quad_s2z_{id} = dot_self(to_vector(white_group_s2z)) ",
      "+ dot_self(z_mean_s2z_{id});\n"
    )
  }
  str_add(out$tpar_comp) <- glue(
    "  {{\n",
    "    matrix[M_{id}, N_{id}] white_group_s2z = ",
    "mdivide_left_tri_low_brms(L_Sigma_s2z_{id}, r_s2z_{id}');\n",
    "{explicit_group_quad_code}",
    "  }}\n",
    cglue("  {r_s2z} = r_s2z_{id}[, {J}];\n")
  )
  str_add(out$pll_args) <- cglue(", vector {r_s2z}")

  str_add(out$tpar_prior) <- glue(
    "  lprior += -0.5 * group_quad_s2z_{id}\n",
    "    - (N_{id} - 1) * ",
    "sum(log(diagonal(L_Sigma_s2z_{id})))\n",
    str_if(
      is_student,
      glue("    - M_{id} * sum(log(group_scale_s2z_{id}))\n")
    ),
    str_if(
      normalize,
      glue("    - 0.5 * N_{id} * M_{id} * log(2 * pi())\n")
    ),
    "  ;\n"
  )
  out
}

# Score and reconstruct one complete explicit-mean system after its blocks
# have been generated. This path is shared by singleton and multiblock models.
.stan_re_s2z_explicit_system <- function(set, normalize, out = list()) {
  infos <- set$infos
  info <- infos[[1L]]
  q <- length(info$qnames)
  set_id <- set$set_id
  p <- info$p
  active_q <- info$active_q
  if (is.null(active_q)) {
    active_q <- which(vapply(
      info$prior, function(spec) !identical(spec$dist, "flat"), logical(1)
    ))
  }

  str_add(out$tpar_def) <- glue(
    "  vector[{q}] q_explicit_s2z_{set_id};",
    "  // conventional coefficients for exact explicit S2Z means\n"
  )
  str_add(out$tpar_comp) <- glue(
    "  q_explicit_s2z_{set_id} = theta_s2z{p};\n"
  )
  for (block in infos) {
    str_add(out$tpar_comp) <- glue(
      "  q_explicit_s2z_{set_id} -= H_s2z_{block$id} * ",
      "mean_r_s2z_{block$id};\n"
    )
  }
  for (k in active_q) {
    str_add(out$tpar_prior) <- stan_re_s2z_prior_target(
      info$prior[[k]], glue("q_explicit_s2z_{set_id}[{k}]"), normalize
    )
  }

  str_add(out$gen_def) <- glue(
    "  vector[{q}] q_recovered_s2z_{set_id};\n"
  )
  str_add(out$gen_def) <- stan_re_s2z_public_fe_def(info)
  for (block in infos) {
    r <- block$r
    id <- block$id
    J <- seq_rows(r)
    idp <- paste0(r$id, usc(combine_prefix(check_prefix(r))))
    r_s2z <- glue("r_s2z_{idp}_{r$cn}")
    r_public <- glue("r_{idp}_{r$cn}")
    is_cor <- nrow(r) > 1L && isTRUE(r$cor[1])
    if (is_cor) {
      str_add(out$gen_def) <- glue(
        "  matrix[N_{id}, M_{id}] r_{id};\n",
        "  corr_matrix[M_{id}] Cor_{id}",
        " = multiply_lower_tri_self_transpose(L_{id});\n",
        "  vector<lower=-1,upper=1>[NC_{id}] cor_{id};\n"
      )
    }
    str_add(out$gen_def) <- cglue(
      "  vector[N_{id}] {r_public};\n"
    )
    str_add(out$gen_comp) <- stan_re_s2z_public_re_comp(
      r, indent = "  "
    )
  }
  str_add(out$gen_comp) <- glue(
    "  q_recovered_s2z_{set_id} = q_explicit_s2z_{set_id};\n"
  )
  for (block in infos) {
    r <- block$r
    id <- block$id
    if (nrow(r) > 1L && isTRUE(r$cor[1])) {
      str_add(out$gen_comp) <- stan_cor_gen_comp(
        cor = glue("cor_{id}"), ncol = glue("M_{id}")
      )
    }
  }
  str_add(out$gen_comp) <- stan_re_s2z_public_fe_comp(info, set_id)
  out
}

# Stan code local to one covariance block participating in a joint S2Z
# omitted-mean system. The block contributes its physical zero-sum effects and
# the Gaussian normal equations for its omitted mean in reference-whitened
# coordinates. Population-level priors and mean recovery are handled once by
# .stan_re_s2z_joint().
.stan_re_s2z_joint_block <- function(id, set, bframe, prior, threads,
                                     normalize, out = list()) {
  lpdf <- stan_lpdf_name(normalize)
  if (is.null(out[["tpar_prior"]])) {
    out[["tpar_prior"]] <- ""
  }
  if (.stan_re_s2z_uses_explicit_mean(set)) {
    return(.stan_re_s2z_explicit_block(
      id, set = set, prior = prior, normalize = normalize, out = out
    ))
  }
  r <- subset2(bframe$frame$re, id = id)
  stopifnot(is.reframe(r), has_rows(r), all(r$s2z), id %in% set$ids)
  bfl <- set$bfl
  info <- set$infos[[match(id, set$ids)]]
  q <- length(info$qnames)
  M <- nrow(r)
  J <- seq_rows(r)
  px <- check_prefix(r)
  idp <- paste0(r$id, usc(combine_prefix(px)))
  r_s2z <- glue("r_s2z_{idp}_{r$cn}")
  is_cor <- M > 1L && isTRUE(r$cor[1])
  is_student <- identical(r$dist[1], "student")
  use_matheron <- .stan_re_s2z_joint_uses_matheron(set)
  if (is_cor) {
    str_add(out$data) <- glue(
      "  int<lower=1> NC_{id};  // number of group-level correlations\n"
    )
    str_add_list(out) <- stan_prior(
      prior, class = "L", group = r$group[1], suffix = usc(id),
      type = glue("cholesky_factor_corr[M_{id}]"),
      comment = "cholesky factor of correlation matrix",
      normalize = normalize
    )
  }
  if (identical(id, set$set_id)) {
    str_add(out$fun) <- "  #include 'fun_sum_to_zero.stan'\n"
  }

  str_add(out$par) <- glue(
    "  vector[M_{id} * (N_{id} - 1)] z_s2z_{id};",
    "  // physical orthonormal S2Z coordinates\n"
  )

  str_add(out$tpar_def) <- glue(
    "  // S2Z block {id} in a joint omitted-mean system\n",
    "  matrix[N_{id}, M_{id}] r_s2z_{id};\n",
    "  matrix[M_{id}, M_{id}] L_Sigma_s2z_{id};\n",
    str_if(
      !use_matheron,
      glue(
        "  real<lower=0> group_info_s2z_{id};",
        "  // isotropic omitted-mean precision\n",
        "  vector[M_{id}] h_group_s2z_{id};\n"
      )
    ),
    "  real<lower=0> group_quad_s2z_{id};\n",
    str_if(
      is_student,
      glue(
        "  vector<lower=0>[N_{id}] group_scale_s2z_{id};\n",
        "  vector<lower=0>[N_{id}] group_prec_s2z_{id};\n"
      )
    ),
    "  // vectors used by the observation-level likelihood\n",
    cglue("  vector[N_{id}] {r_s2z};\n")
  )

  # This map is a deterministic function of formula structure and predictor
  # centering constants, so construct it outside the autodiff graph.
  str_add(out$tdata_def) <- glue(
    "  matrix[{q}, M_{id}] H_s2z_{id};\n"
  )
  str_add(out$tdata_comp) <- glue(
    "  H_s2z_{id} = rep_matrix(0.0, {q}, M_{id});\n"
  )
  for (j in seq_len(M)) {
    qi <- info$match_q[j]
    str_add(out$tdata_comp) <- glue(
      "  H_s2z_{id}[{qi}, {j}] = 1.0;\n"
    )
    if (info$center && info$r$coef[j] != "Intercept") {
      str_add(out$tdata_comp) <- glue(
        "  H_s2z_{id}[1, {j}] = means_X{info$p}[{qi - 1L}];\n"
      )
    }
  }

  if (is_student) {
    tr <- subset_reframe_dist(r, "student")
    g <- usc(tr$ggn[1])
    str_add(out$tpar_comp) <- glue(
      "  group_scale_s2z_{id} = dfm{g};\n",
      "  group_prec_s2z_{id} = inv_square(group_scale_s2z_{id});\n"
    )
  }

  scale <- glue("sd_{id}")
  if (is_cor) {
    str_add(out$tpar_comp) <- glue(
      "  L_Sigma_s2z_{id} = diag_pre_multiply({scale}, L_{id});\n"
    )
  } else {
    str_add(out$tpar_comp) <- glue(
      "  L_Sigma_s2z_{id} = diag_matrix({scale});\n"
    )
  }

  str_add(out$tpar_comp) <- glue(
    "  for (k in 1:M_{id}) {{\n",
    "    r_s2z_{id}[, k] = sum_to_zero_constrain_brms(",
    "segment(z_s2z_{id}, (k - 1) * (N_{id} - 1) + 1, ",
    "N_{id} - 1));\n",
    "  }}\n"
  )

  if (use_matheron && is_cor) {
    str_add(out$tpar_comp) <- glue(
      "  {{\n",
      "    matrix[M_{id}, N_{id}] white_group_s2z = ",
      "mdivide_left_tri_low_brms(L_Sigma_s2z_{id}, r_s2z_{id}');\n",
      "    group_quad_s2z_{id} = dot_self(to_vector(white_group_s2z));\n",
      "  }}\n"
    )
  } else if (use_matheron) {
    str_add(out$tpar_comp) <- glue(
      "  {{\n",
      "    group_quad_s2z_{id} = 0.0;\n",
      "    for (k in 1:M_{id}) {{\n",
      "      vector[N_{id}] white_group_s2z = ",
      "r_s2z_{id}[, k] / sd_{id}[k];\n",
      "      group_quad_s2z_{id} += dot_self(white_group_s2z);\n",
      "    }}\n",
      "  }}\n"
    )
  } else if (is_cor) {
    group_info <- if (is_student) {
      glue("sum(group_prec_s2z_{id})")
    } else {
      glue("1.0 * N_{id}")
    }
    group_score_code <- if (is_student) {
      glue(
        "    h_group_s2z_{id} = -white_group_s2z * ",
        "group_prec_s2z_{id};\n",
        "    group_quad_s2z_{id} = columns_dot_self(white_group_s2z) * ",
        "group_prec_s2z_{id};\n"
      )
    } else {
      glue(
        "    h_group_s2z_{id} = zeros_vector(M_{id});\n",
        "    group_quad_s2z_{id} = dot_self(to_vector(white_group_s2z));\n"
      )
    }
    str_add(out$tpar_comp) <- glue(
      "  {{\n",
      "    matrix[M_{id}, N_{id}] white_group_s2z = ",
      "mdivide_left_tri_low_brms(L_Sigma_s2z_{id}, r_s2z_{id}');\n",
      "    group_info_s2z_{id} = {group_info};\n",
      "{group_score_code}",
      "  }}\n"
    )
  } else {
    group_info <- if (is_student) {
      glue("sum(group_prec_s2z_{id})")
    } else {
      glue("1.0 * N_{id}")
    }
    group_score_code <- if (is_student) {
      glue(
        "    h_group_s2z_{id} = zeros_vector(M_{id});\n",
        "    group_quad_s2z_{id} = 0.0;\n",
        "    for (k in 1:M_{id}) {{\n",
        "      vector[N_{id}] white_group_s2z = r_s2z_{id}[, k] / ",
        "{scale}[k];\n",
        "      h_group_s2z_{id}[k] = -dot_product(",
        "white_group_s2z, group_prec_s2z_{id});\n",
        "      group_quad_s2z_{id} += dot_product(",
        "group_prec_s2z_{id}, square(white_group_s2z));\n",
        "    }}\n"
      )
    } else {
      glue(
        "    h_group_s2z_{id} = zeros_vector(M_{id});\n",
        "    group_quad_s2z_{id} = 0.0;\n",
        "    for (k in 1:M_{id}) {{\n",
        "      vector[N_{id}] white_group_s2z = r_s2z_{id}[, k] / ",
        "{scale}[k];\n",
        "      group_quad_s2z_{id} += dot_self(white_group_s2z);\n",
        "    }}\n"
      )
    }
    str_add(out$tpar_comp) <- glue(
      "  {{\n",
      "    group_info_s2z_{id} = {group_info};\n",
      "{group_score_code}",
      "  }}\n"
    )
  }
  str_add(out$tpar_comp) <- cglue(
    "  {r_s2z} = r_s2z_{id}[, {J}];\n"
  )
  str_add(out$pll_args) <- cglue(", vector {r_s2z}")
  out
}

# Fast separated Gaussian path for multiple S2Z blocks. Rather than factoring
# the covariance of every omitted block mean, it factors only the induced
# covariance of theta and uses one Matheron update to recover all conventional
# means and population coefficients jointly.
.stan_re_s2z_joint_matheron <- function(set, prior, threads, normalize, ...) {
  out <- list(tpar_prior = "")
  infos <- set$infos
  matheron <- .stan_re_s2z_joint_matheron_info(set)
  stopifnot(!is.null(matheron))
  info <- infos[[1L]]
  q <- length(info$qnames)
  P <- matheron$P
  inactive <- matheron$inactive
  rdim <- length(P)
  P_index <- paste0("{", paste(P, collapse = ", "), "}")
  set_id <- set$set_id
  p <- info$p
  lpdf <- stan_lpdf_name(normalize)

  # A Student-t population prior is conditionally Gaussian, so the same
  # update applies after drawing its usual inverse-chi-square mixture scale.
  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    if (spec$dist == "student") {
      str_add(out$par) <- glue(
        "  real<lower=0> udf_b_s2z{p}_{k};",
        "  // mixing variable for population coefficient {k}\n"
      )
      str_add(out$tpar_prior) <- glue(
        "  lprior += inv_chi_square_{lpdf}(udf_b_s2z{p}_{k} | ",
        "{stan_s2z_number(spec$df)});\n"
      )
    }
  }

  str_add(out$tpar_def) <- glue(
    "  // fast Gaussian Matheron system for S2Z blocks ",
    "{paste(set$ids, collapse = ', ')}\n",
    "  vector[{q}] prior_mean_s2z_{set_id};\n",
    "  vector<lower=0>[{q}] prior_scale_s2z_{set_id};\n",
    str_if(
      rdim == 1L,
      glue(
        "  real<lower=0> W_matheron_s2z_{set_id};\n",
        "  real<lower=0> sqrt_W_matheron_s2z_{set_id};\n",
        "  real theta_white_matheron_s2z_{set_id};\n"
      )
    ),
    str_if(
      rdim > 1L,
      glue(
        "  matrix[{rdim}, {rdim}] W_matheron_s2z_{set_id};\n",
        "  matrix[{rdim}, {rdim}] L_W_matheron_s2z_{set_id};\n",
        "  vector[{rdim}] theta_white_matheron_s2z_{set_id};\n"
      )
    ),
    "  real joint_quad_s2z_{set_id};\n"
  )
  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    str_add(out$tpar_comp) <- glue(
      "  prior_mean_s2z_{set_id}[{k}] = ",
      "{stan_s2z_number(spec$location)};\n"
    )
    cond_scale <- if (spec$dist == "flat") {
      "1.0"
    } else if (spec$dist == "normal") {
      stan_s2z_number(spec$scale)
    } else {
      glue(
        "{stan_s2z_number(spec$scale)} * sqrt(",
        "{stan_s2z_number(spec$df)} * udf_b_s2z{p}_{k})"
      )
    }
    str_add(out$tpar_comp) <- glue(
      "  prior_scale_s2z_{set_id}[{k}] = {cond_scale};\n"
    )
  }
  if (rdim == 1L) {
    str_add(out$tpar_comp) <- glue(
      "  W_matheron_s2z_{set_id} = ",
      "square(prior_scale_s2z_{set_id}[{P}]);\n"
    )
  }
  str_add(out$tpar_comp) <- glue(
    "  joint_quad_s2z_{set_id} = 0.0;\n"
  )
  for (b in seq_along(infos)) {
    id <- infos[[b]]$id
    if (rdim == 1L) {
      str_add(out$tpar_comp) <- glue(
        "  W_matheron_s2z_{set_id} += dot_self(",
        "H_s2z_{id}[{P}, ] * L_Sigma_s2z_{id}) / (1.0 * N_{id});\n"
      )
    } else if (rdim > 1L) {
      W_update <- if (b == 1L) {
        glue(
          "add_diag(\n",
          "      tcrossprod(H_active_s2z * L_Sigma_s2z_{id}) / ",
          "(1.0 * N_{id}),\n",
          "      square(prior_scale_s2z_{set_id}[{P_index}])\n",
          "    )"
        )
      } else {
        glue(
          "tcrossprod(",
          "H_active_s2z * L_Sigma_s2z_{id}) / (1.0 * N_{id})"
        )
      }
      W_operator <- if (b == 1L) "=" else "+="
      str_add(out$tpar_comp) <- glue(
        "  {{\n",
        "    matrix[{rdim}, M_{id}] H_active_s2z = ",
        "H_s2z_{id}[{P_index}, ];\n",
        "    W_matheron_s2z_{set_id} {W_operator} {W_update};\n",
        "  }}\n"
      )
    }
    str_add(out$tpar_comp) <- glue(
      "  joint_quad_s2z_{set_id} += group_quad_s2z_{id};\n"
    )
  }
  if (rdim == 1L) {
    str_add(out$tpar_comp) <- glue(
      "  sqrt_W_matheron_s2z_{set_id} = sqrt(W_matheron_s2z_{set_id});\n",
      "  theta_white_matheron_s2z_{set_id} = ",
      "(theta_s2z{p}[{P}] - prior_mean_s2z_{set_id}[{P}]) / ",
      "sqrt_W_matheron_s2z_{set_id};\n"
    )
  } else if (rdim > 1L) {
    str_add(out$tpar_comp) <- glue(
      "  L_W_matheron_s2z_{set_id} = ",
      "cholesky_decompose_brms(W_matheron_s2z_{set_id});\n",
      "  {{\n",
      "    vector[{rdim}] theta_difference_s2z = ",
      "theta_s2z{p}[{P_index}] - prior_mean_s2z_{set_id}[{P_index}];\n",
      "    theta_white_matheron_s2z_{set_id} = mdivide_left_tri_low_brms(",
      "L_W_matheron_s2z_{set_id}, theta_difference_s2z);\n",
      "  }}\n"
    )
  }

  for (k in inactive) {
    str_add(out$tpar_prior) <- glue(
      "  lprior += normal_{lpdf}(theta_s2z{p}[{k}] | ",
      "prior_mean_s2z_{set_id}[{k}], prior_scale_s2z_{set_id}[{k}]);\n"
    )
  }
  str_add(out$tpar_prior) <- glue(
    "  lprior += -0.5 * joint_quad_s2z_{set_id}\n"
  )
  if (rdim == 1L) {
    str_add(out$tpar_prior) <- glue(
      "    - 0.5 * square(theta_white_matheron_s2z_{set_id})\n",
      "    - log(sqrt_W_matheron_s2z_{set_id})\n"
    )
  } else if (rdim > 1L) {
    str_add(out$tpar_prior) <- glue(
      "    - 0.5 * dot_self(theta_white_matheron_s2z_{set_id})\n",
      "    - sum(log(diagonal(L_W_matheron_s2z_{set_id})))\n"
    )
  }
  for (b in seq_along(infos)) {
    id <- infos[[b]]$id
    str_add(out$tpar_prior) <- glue(
      "    - (N_{id} - 1) * ",
      "sum(log(diagonal(L_Sigma_s2z_{id})))\n"
    )
    if (normalize) {
      str_add(out$tpar_prior) <- glue(
        "    - 0.5 * (N_{id} - 1) * M_{id} * log(2 * pi())\n"
      )
    }
  }
  if (normalize && rdim > 0L) {
    str_add(out$tpar_prior) <- glue(
      "    - 0.5 * {rdim} * log(2 * pi())\n"
    )
  }
  str_add(out$tpar_prior) <- "  ;\n"

  str_add(out$gen_def) <- glue(
    "  vector[{q}] q_recovered_s2z_{set_id};\n"
  )
  str_add(out$gen_def) <- stan_re_s2z_public_fe_def(info)
  for (b in seq_along(infos)) {
    r <- infos[[b]]$r
    id <- infos[[b]]$id
    idp <- paste0(r$id, usc(combine_prefix(check_prefix(r))))
    r_public <- glue("r_{idp}_{r$cn}")
    is_cor <- nrow(r) > 1L && isTRUE(r$cor[1])
    str_add(out$gen_def) <- glue(
      "  vector[M_{id}] mean_r_s2z_{id};\n"
    )
    if (is_cor) {
      str_add(out$gen_def) <- glue(
        "  matrix[N_{id}, M_{id}] r_{id};\n"
      )
    }
    str_add(out$gen_def) <- cglue(
      "  vector[N_{id}] {r_public};\n"
    )
    if (is_cor) {
      str_add(out$gen_def) <- glue(
        "  // compute group-level correlations\n",
        "  corr_matrix[M_{id}] Cor_{id}",
        " = multiply_lower_tri_self_transpose(L_{id});\n",
        "  vector<lower=-1,upper=1>[NC_{id}] cor_{id};\n"
      )
    }
  }

  # Draw every m_f* and, when needed, the proper active part of beta*. One
  # shared innovation conditions them on theta. Flat active coordinates do not
  # enter W; beta is recovered as theta - sum(H_f m_f), which also enforces the
  # predictor identity to floating-point precision.
  str_add(out$gen_comp) <- glue(
    "  {{\n"
  )
  if (rdim > 0L) {
    str_add(out$gen_comp) <- glue(
      "    vector[{rdim}] theta_star_s2z;\n",
      "    vector[{rdim}] delta_s2z;\n"
    )
    for (a in seq_len(rdim)) {
      str_add(out$gen_comp) <- glue(
        "    theta_star_s2z[{a}] = prior_mean_s2z_{set_id}[{P[a]}] + ",
        "prior_scale_s2z_{set_id}[{P[a]}] * std_normal_rng();\n"
      )
    }
  }
  for (b in seq_along(infos)) {
    id <- infos[[b]]$id
    str_add(out$gen_comp) <- glue(
      "    {{\n",
      "      vector[M_{id}] z_mean_s2z;\n",
      "      for (k in 1:M_{id}) z_mean_s2z[k] = std_normal_rng();\n",
      "      mean_r_s2z_{id} = L_Sigma_s2z_{id} * z_mean_s2z / ",
      "sqrt(1.0 * N_{id});\n",
      "    }}\n"
    )
    if (rdim > 0L) {
      for (a in seq_len(rdim)) {
        str_add(out$gen_comp) <- glue(
          "    theta_star_s2z[{a}] += dot_product(",
          "H_s2z_{id}[{P[a]}, ], mean_r_s2z_{id}');\n"
        )
      }
    }
  }
  if (rdim == 1L) {
    str_add(out$gen_comp) <- glue(
      "    delta_s2z[1] = (theta_s2z{p}[{P}] - theta_star_s2z[1]) / ",
      "W_matheron_s2z_{set_id};\n"
    )
  } else if (rdim > 1L) {
    str_add(out$gen_comp) <- glue(
      "    {{\n",
      "      vector[{rdim}] theta_innovation_s2z;\n",
      "      vector[{rdim}] forward_solve_s2z;\n"
    )
    for (a in seq_len(rdim)) {
      str_add(out$gen_comp) <- glue(
        "      theta_innovation_s2z[{a}] = theta_s2z{p}[{P[a]}] - ",
        "theta_star_s2z[{a}];\n"
      )
    }
    str_add(out$gen_comp) <- glue(
      "      forward_solve_s2z = mdivide_left_tri_low_brms(",
      "L_W_matheron_s2z_{set_id}, theta_innovation_s2z);\n",
      "      delta_s2z = (mdivide_right_tri_low_brms(",
      "forward_solve_s2z', L_W_matheron_s2z_{set_id}))';\n",
      "    }}\n"
    )
  }
  str_add(out$gen_comp) <- glue(
    "    q_recovered_s2z_{set_id} = theta_s2z{p};\n"
  )
  for (b in seq_along(infos)) {
    r <- infos[[b]]$r
    id <- infos[[b]]$id
    J <- seq_rows(r)
    idp <- paste0(r$id, usc(combine_prefix(check_prefix(r))))
    r_s2z <- glue("r_s2z_{idp}_{r$cn}")
    r_public <- glue("r_{idp}_{r$cn}")
    is_cor <- nrow(r) > 1L && isTRUE(r$cor[1])
    if (rdim > 0L) {
      str_add(out$gen_comp) <- glue(
        "    {{\n",
        "      vector[M_{id}] active_score_s2z = ",
        "zeros_vector(M_{id});\n"
      )
      for (a in seq_len(rdim)) {
        str_add(out$gen_comp) <- glue(
          "      active_score_s2z += H_s2z_{id}[{P[a]}, ]' * ",
          "delta_s2z[{a}];\n"
        )
      }
      str_add(out$gen_comp) <- glue(
        "      mean_r_s2z_{id} += L_Sigma_s2z_{id} * ",
        "(L_Sigma_s2z_{id}' * active_score_s2z) / (1.0 * N_{id});\n",
        "    }}\n"
      )
    }
    str_add(out$gen_comp) <- glue(
      "    q_recovered_s2z_{set_id} -= H_s2z_{id} * mean_r_s2z_{id};\n"
    )
    str_add(out$gen_comp) <- stan_re_s2z_public_re_comp(
      r, indent = "    "
    )
  }
  str_add(out$gen_comp) <- "  }\n"
  for (b in seq_along(infos)) {
    r <- infos[[b]]$r
    id <- infos[[b]]$id
    if (nrow(r) > 1L && isTRUE(r$cor[1])) {
      str_add(out$gen_comp) <- stan_cor_gen_comp(
        cor = glue("cor_{id}"), ncol = glue("M_{id}")
      )
    }
  }
  str_add(out$gen_comp) <- stan_re_s2z_public_fe_comp(info, set_id)
  out
}

# Jointly integrate the omitted group-effect means of all physical S2Z blocks
# in one linear predictor. Block means are represented in their own reference-
# whitened coordinates, so the independent group hierarchies contribute a
# block-diagonal precision while population-level priors provide the only
# cross-block coupling.
.stan_re_s2z_joint <- function(set, prior, threads, normalize, ...) {
  out <- list(tpar_prior = "")
  infos <- set$infos
  stopifnot(length(infos) > 1L)
  if (.stan_re_s2z_uses_explicit_mean(set)) {
    return(.stan_re_s2z_explicit_system(
      set, normalize = normalize, out = out
    ))
  }
  if (.stan_re_s2z_joint_uses_matheron(set)) {
    return(.stan_re_s2z_joint_matheron(
      set, prior = prior, threads = threads, normalize = normalize, ...
    ))
  }
  info <- infos[[1L]]
  q <- length(info$qnames)
  Ms <- vapply(infos, function(x) nrow(x$r), integer(1))
  total_M <- sum(Ms)
  starts <- cumsum(c(1L, utils::head(Ms, -1L)))
  ends <- starts + Ms - 1L
  set_id <- set$set_id
  p <- info$p
  lpdf <- stan_lpdf_name(normalize)

  # Population Student-t and Cauchy priors use one mixing variable per fixed
  # coefficient for the entire joint system, rather than one copy per block.
  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    if (spec$dist == "student") {
      str_add(out$par) <- glue(
        "  real<lower=0> udf_b_s2z{p}_{k};",
        "  // mixing variable for population coefficient {k}\n"
      )
      str_add(out$tpar_prior) <- glue(
        "  lprior += inv_chi_square_{lpdf}(udf_b_s2z{p}_{k} | ",
        "{stan_s2z_number(spec$df)});\n"
      )
    }
  }

  str_add(out$tpar_def) <- glue(
    "  // joint omitted-mean system for S2Z blocks ",
    "{paste(set$ids, collapse = ', ')}\n",
    "  vector[{q}] prior_mean_s2z_{set_id};\n",
    "  vector<lower=0>[{q}] prior_prec_s2z_{set_id};\n",
    "  matrix[{total_M}, {total_M}] P_s2z_{set_id};\n",
    "  matrix[{total_M}, {total_M}] L_P_s2z_{set_id};\n",
    "  vector[{total_M}] h_joint_s2z_{set_id};\n",
    "  vector[{total_M}] mhat_s2z_{set_id};\n",
    "  real joint_quad_s2z_{set_id};\n"
  )
  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    loc <- stan_s2z_number(spec$location)
    str_add(out$tpar_comp) <- glue(
      "  prior_mean_s2z_{set_id}[{k}] = {loc};\n"
    )
    if (spec$dist == "flat") {
      prec <- "0.0"
    } else if (spec$dist == "normal") {
      prec <- glue("inv_square({stan_s2z_number(spec$scale)})")
    } else {
      prec <- glue(
        "inv_square({stan_s2z_number(spec$scale)} * sqrt(",
        "{stan_s2z_number(spec$df)} * udf_b_s2z{p}_{k}))"
      )
    }
    str_add(out$tpar_comp) <- glue(
      "  prior_prec_s2z_{set_id}[{k}] = {prec};\n"
    )
  }

  str_add(out$tpar_comp) <- glue(
    "  joint_quad_s2z_{set_id} = 0.0;\n",
    "  {{\n",
    "    matrix[{q}, {total_M}] prior_factor_s2z;\n",
    "    vector[{q}] prior_difference_s2z = ",
    "sqrt(prior_prec_s2z_{set_id}) .* (theta_s2z{p} - ",
    "prior_mean_s2z_{set_id});\n",
    "    vector[{total_M}] forward_solve_s2z;\n"
  )
  for (b in seq_along(infos)) {
    id <- infos[[b]]$id
    take <- glue("{starts[b]}:{ends[b]}")
    str_add(out$tpar_comp) <- glue(
      "    prior_factor_s2z[, {take}] = diag_pre_multiply(",
      "sqrt(prior_prec_s2z_{set_id}), H_s2z_{id} * ",
      "L_Sigma_s2z_{id});\n"
    )
  }
  str_add(out$tpar_comp) <- glue(
    "    P_s2z_{set_id} = crossprod(prior_factor_s2z);\n"
  )
  for (b in seq_along(infos)) {
    id <- infos[[b]]$id
    take <- glue("{starts[b]}:{ends[b]}")
    index <- if (starts[b] == 1L) "k" else glue("{starts[b] - 1L} + k")
    str_add(out$tpar_comp) <- glue(
      "    for (k in 1:M_{id}) {{\n",
      "      P_s2z_{set_id}[{index}, {index}] += group_info_s2z_{id};\n",
      "    }}\n",
      "    h_joint_s2z_{set_id}[{take}] = h_group_s2z_{id};\n",
      "    joint_quad_s2z_{set_id} += group_quad_s2z_{id};\n"
    )
  }
  str_add(out$tpar_comp) <- glue(
    "    h_joint_s2z_{set_id} += prior_factor_s2z' * ",
    "prior_difference_s2z;\n",
    "    L_P_s2z_{set_id} = cholesky_decompose_brms(P_s2z_{set_id});\n",
    "    forward_solve_s2z = mdivide_left_tri_low_brms(",
    "L_P_s2z_{set_id}, h_joint_s2z_{set_id});\n",
    "    mhat_s2z_{set_id} = (mdivide_right_tri_low_brms(",
    "forward_solve_s2z', L_P_s2z_{set_id}))';\n",
    "    joint_quad_s2z_{set_id} -= dot_self(forward_solve_s2z);\n",
    "  }}\n"
  )

  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    if (spec$dist == "flat") {
      next
    }
    if (spec$dist == "normal") {
      cond_scale <- stan_s2z_number(spec$scale)
    } else {
      cond_scale <- glue(
        "{stan_s2z_number(spec$scale)} * sqrt(",
        "{stan_s2z_number(spec$df)} * udf_b_s2z{p}_{k})"
      )
    }
    str_add(out$tpar_prior) <- glue(
      "  lprior += normal_{lpdf}(theta_s2z{p}[{k}] | ",
      "{stan_s2z_number(spec$location)}, {cond_scale});\n"
    )
  }
  str_add(out$tpar_prior) <- glue(
    "  lprior += -0.5 * joint_quad_s2z_{set_id}\n",
    "    - sum(log(diagonal(L_P_s2z_{set_id})))\n"
  )
  for (b in seq_along(infos)) {
    r <- infos[[b]]$r
    id <- infos[[b]]$id
    str_add(out$tpar_prior) <- glue(
      "    - (N_{id} - 1) * ",
      "sum(log(diagonal(L_Sigma_s2z_{id})))\n"
    )
    if (identical(r$dist[1], "student")) {
      str_add(out$tpar_prior) <- glue(
        "    - M_{id} * sum(log(group_scale_s2z_{id}))\n"
      )
    }
    if (normalize) {
      str_add(out$tpar_prior) <- glue(
        "    - 0.5 * N_{id} * M_{id} * log(2 * pi())\n",
        "    + 0.5 * M_{id} * log(1.0 * N_{id})\n"
      )
    }
  }
  if (normalize) {
    str_add(out$tpar_prior) <- glue(
      "    + 0.5 * {total_M} * log(2 * pi())\n"
    )
  }
  str_add(out$tpar_prior) <- "  ;\n"

  # One conditional draw preserves cross-block posterior covariance in the
  # conventional population/group parameterization.
  str_add(out$gen_def) <- glue(
    "  vector[{q}] q_recovered_s2z_{set_id};\n"
  )
  str_add(out$gen_def) <- stan_re_s2z_public_fe_def(info)
  for (b in seq_along(infos)) {
    r <- infos[[b]]$r
    id <- infos[[b]]$id
    idp <- paste0(r$id, usc(combine_prefix(check_prefix(r))))
    r_public <- glue("r_{idp}_{r$cn}")
    is_cor <- nrow(r) > 1L && isTRUE(r$cor[1])
    str_add(out$gen_def) <- glue(
      "  vector[M_{id}] mean_r_s2z_{id};\n"
    )
    if (is_cor) {
      str_add(out$gen_def) <- glue(
        "  matrix[N_{id}, M_{id}] r_{id};\n"
      )
    }
    str_add(out$gen_def) <- cglue(
      "  vector[N_{id}] {r_public};\n"
    )
    if (is_cor) {
      str_add(out$gen_def) <- glue(
        "  // compute group-level correlations\n",
        "  corr_matrix[M_{id}] Cor_{id}",
        " = multiply_lower_tri_self_transpose(L_{id});\n",
        "  vector<lower=-1,upper=1>[NC_{id}] cor_{id};\n"
      )
    }
  }

  str_add(out$gen_comp) <- glue(
    "  {{\n",
    "    vector[{total_M}] z_mean_s2z;\n",
    "    vector[{total_M}] mean_white_s2z;\n",
    "    for (k in 1:{total_M}) z_mean_s2z[k] = std_normal_rng();\n",
    "    mean_white_s2z = mhat_s2z_{set_id} + ",
    "(mdivide_right_tri_low_brms(z_mean_s2z', L_P_s2z_{set_id}))';\n",
    "    q_recovered_s2z_{set_id} = theta_s2z{p};\n"
  )
  for (b in seq_along(infos)) {
    r <- infos[[b]]$r
    id <- infos[[b]]$id
    M <- nrow(r)
    J <- seq_rows(r)
    take <- glue("{starts[b]}:{ends[b]}")
    idp <- paste0(r$id, usc(combine_prefix(check_prefix(r))))
    r_s2z <- glue("r_s2z_{idp}_{r$cn}")
    r_public <- glue("r_{idp}_{r$cn}")
    is_cor <- M > 1L && isTRUE(r$cor[1])
    mean_transform <- if (is_cor) {
      glue("L_Sigma_s2z_{id} * mean_white_s2z[{take}]")
    } else {
      glue("sd_{id} .* mean_white_s2z[{take}]")
    }
    str_add(out$gen_comp) <- glue(
      "    mean_r_s2z_{id} = {mean_transform};\n",
      "    q_recovered_s2z_{set_id} -= H_s2z_{id} * ",
      "mean_r_s2z_{id};\n"
    )
    str_add(out$gen_comp) <- stan_re_s2z_public_re_comp(
      r, indent = "    "
    )
  }
  str_add(out$gen_comp) <- "  }\n"
  for (b in seq_along(infos)) {
    r <- infos[[b]]$r
    id <- infos[[b]]$id
    if (nrow(r) > 1L && isTRUE(r$cor[1])) {
      str_add(out$gen_comp) <- stan_cor_gen_comp(
        cor = glue("cor_{id}"), ncol = glue("M_{id}")
      )
    }
  }
  str_add(out$gen_comp) <- stan_re_s2z_public_fe_comp(info, set_id)
  out
}

# Stan code for one physical sum-to-zero group-level block. The omitted common
# group-effect mean vector is integrated out analytically. Conditional
# Gaussian scale mixtures are handled by a group-specific scale, which makes
# this exact for both Gaussian and Student-t group effects.
.stan_re_s2z <- function(id, system, bframe, prior, threads, normalize,
                         out = list()) {
  lpdf <- stan_lpdf_name(normalize)
  # Avoid partial matching of $tpar_prior to $tpar_prior_const when a group
  # scale is fixed. Otherwise the constant assignment is appended to the
  # model block when normalize = FALSE.
  if (is.null(out[["tpar_prior"]])) {
    out[["tpar_prior"]] <- ""
  }
  r <- subset2(bframe$frame$re, id = id)
  stopifnot(is.reframe(r), has_rows(r), all(r$s2z))
  info <- system$infos[[1L]]
  q <- length(info$qnames)
  M <- nrow(r)
  J <- seq_rows(r)
  p <- info$p
  idp <- paste0(r$id, usc(combine_prefix(check_prefix(r))))
  is_cor <- M > 1L && isTRUE(r$cor[1])
  is_student <- identical(r$dist[1], "student")

  if (M == 1L) {
    return(.stan_re_s2z_scalar(
      id, r = r, info = info, normalize = normalize, out = out
    ))
  }

  if (!is_cor) {
    return(.stan_re_s2z_independent(
      id, r = r, info = info, normalize = normalize, out = out
    ))
  }

  # Keep brms's existing covariance parameters and priors. Only the
  # group-effect coordinates and their joint prior are replaced.
  if (is_cor) {
    str_add(out$data) <- glue(
      "  int<lower=1> NC_{id};  // number of group-level correlations\n"
    )
    str_add_list(out) <- stan_prior(
      prior, class = "L", group = r$group[1], suffix = usc(id),
      type = glue("cholesky_factor_corr[M_{id}]"),
      comment = "cholesky factor of correlation matrix",
      normalize = normalize
    )
  }

  str_add(out$fun) <- "  #include 'fun_sum_to_zero.stan'\n"
  str_add(out$par) <- glue(
    "  vector[M_{id} * (N_{id} - 1)] z_s2z_{id};",
    "  // physical orthonormal S2Z coordinates\n"
  )

  # Independent Student-t and Cauchy population priors are represented as
  # Gaussian scale mixtures, so the omitted group mean remains analytic.
  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    if (spec$dist == "student") {
      str_add(out$par) <- glue(
        "  real<lower=0> udf_b_s2z{p}_{k};",
        "  // mixing variable for population coefficient {k}\n"
      )
      str_add(out$tpar_prior) <- glue(
        "  lprior += inv_chi_square_{lpdf}(udf_b_s2z{p}_{k} | ",
        "{stan_s2z_number(spec$df)});\n"
      )
    }
  }

  str_add(out$tpar_def) <- glue(
    "  // physical sum-to-zero group-level effects of ID {id}\n",
    "  matrix[N_{id}, M_{id}] r_s2z_{id};\n",
    "  vector[{q}] prior_mean_s2z_{id};\n",
    "  vector<lower=0>[{q}] prior_prec_s2z_{id};\n",
    "  vector<lower=0>[N_{id}] group_scale_s2z_{id};\n",
    "  vector<lower=0>[N_{id}] group_prec_s2z_{id};\n",
    "  matrix[M_{id}, M_{id}] L_Sigma_s2z_{id};\n",
    "  matrix[M_{id}, M_{id}] Q_Sigma_s2z_{id};\n",
    "  matrix[M_{id}, M_{id}] P_s2z_{id};\n",
    "  matrix[M_{id}, M_{id}] L_P_s2z_{id};\n",
    "  vector[M_{id}] mhat_s2z_{id};\n",
    "  vector[{q}] qhat_s2z_{id};\n",
    "  real<lower=0> group_quad_s2z_{id};\n",
    "  // using vectors speeds up indexing in loops\n",
    cglue("  vector[N_{id}] r_s2z_{idp}_{r$cn};\n")
  )

  # The mapping H absorbs the omitted raw group mean into brms's population
  # coordinates. It is fixed by the formula and centered design-column means,
  # so construct it once outside the autodiff graph.
  str_add(out$tdata_def) <- glue(
    "  matrix[{q}, M_{id}] H_s2z_{id};\n"
  )
  str_add(out$tdata_comp) <- glue(
    "  H_s2z_{id} = rep_matrix(0.0, {q}, M_{id});\n"
  )
  for (j in seq_len(M)) {
    qi <- info$match_q[j]
    str_add(out$tdata_comp) <- glue(
      "  H_s2z_{id}[{qi}, {j}] = 1.0;\n"
    )
    if (info$center && info$r$coef[j] != "Intercept") {
      str_add(out$tdata_comp) <- glue(
        "  H_s2z_{id}[1, {j}] = means_X{p}[{qi - 1L}];\n"
      )
    }
  }
  if (is_cor) {
    str_add(out$tpar_comp) <- glue(
      "  L_Sigma_s2z_{id} = diag_pre_multiply(sd_{id}, L_{id});\n"
    )
  } else {
    str_add(out$tpar_comp) <- glue(
      "  L_Sigma_s2z_{id} = diag_matrix(sd_{id});\n"
    )
  }
  str_add(out$tpar_comp) <- glue(
    "  for (k in 1:M_{id}) {{\n",
    "    r_s2z_{id}[, k] = sum_to_zero_constrain_brms(",
    "segment(z_s2z_{id}, (k - 1) * (N_{id} - 1) + 1, ",
    "N_{id} - 1));\n",
    "  }}\n"
  )

  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    loc <- stan_s2z_number(spec$location)
    str_add(out$tpar_comp) <- glue(
      "  prior_mean_s2z_{id}[{k}] = {loc};\n"
    )
    if (spec$dist == "flat") {
      prec <- "0.0"
    } else if (spec$dist == "normal") {
      prec <- glue("inv_square({stan_s2z_number(spec$scale)})")
    } else {
      prec <- glue(
        "inv_square({stan_s2z_number(spec$scale)} * sqrt(",
        "{stan_s2z_number(spec$df)} * udf_b_s2z{p}_{k}))"
      )
    }
    str_add(out$tpar_comp) <- glue(
      "  prior_prec_s2z_{id}[{k}] = {prec};\n"
    )
  }

  if (is_student) {
    tr <- subset_reframe_dist(r, "student")
    g <- usc(tr$ggn[1])
    str_add(out$tpar_comp) <- glue(
      "  group_scale_s2z_{id} = dfm{g};\n"
    )
  } else {
    str_add(out$tpar_comp) <- glue(
      "  group_scale_s2z_{id} = rep_vector(1.0, N_{id});\n"
    )
  }
  direct_group_quad_code <- if (is_student) {
    glue(
      "    group_quad_s2z_{id} = 0.0;\n",
      "    for (j in 1:N_{id}) {{\n",
      "      vector[M_{id}] white_level_s2z = (white_s2z[, j] + ",
      "mean_white_s2z) / group_scale_s2z_{id}[j];\n",
      "      group_quad_s2z_{id} += dot_self(white_level_s2z);\n",
      "    }}\n"
    )
  } else {
    glue(
      "    group_quad_s2z_{id} = dot_self(to_vector(white_s2z)) + ",
      "N_{id} * dot_self(mean_white_s2z);\n"
    )
  }
  str_add(out$tpar_comp) <- glue(
    "  group_prec_s2z_{id} = inv_square(group_scale_s2z_{id});\n",
    "  {{\n",
    "    matrix[M_{id}, M_{id}] L_inv_s2z = mdivide_left_tri_low_brms(",
    "L_Sigma_s2z_{id}, diag_matrix(rep_vector(1.0, M_{id})));\n",
    "    vector[M_{id}] h_s2z;\n",
    "    Q_Sigma_s2z_{id} = crossprod(L_inv_s2z);\n",
    "    P_s2z_{id} = crossprod(diag_pre_multiply(",
    "sqrt(prior_prec_s2z_{id}), H_s2z_{id}))",
    " + sum(group_prec_s2z_{id}) * Q_Sigma_s2z_{id};\n",
    "    h_s2z = H_s2z_{id}' * (prior_prec_s2z_{id} .* ",
    "(theta_s2z{p} - prior_mean_s2z_{id})) - Q_Sigma_s2z_{id} * ",
    "(r_s2z_{id}' * group_prec_s2z_{id});\n",
    "    L_P_s2z_{id} = cholesky_decompose_brms(P_s2z_{id});\n",
    "    vector[M_{id}] forward_solve_s2z = mdivide_left_tri_low_brms(",
    "L_P_s2z_{id}, h_s2z);\n",
    "    mhat_s2z_{id} = (mdivide_right_tri_low_brms(",
    "forward_solve_s2z', L_P_s2z_{id}))';\n",
    "    qhat_s2z_{id} = theta_s2z{p} - H_s2z_{id} * mhat_s2z_{id};\n",
    "    matrix[M_{id}, N_{id}] white_s2z = mdivide_left_tri_low_brms(",
    "L_Sigma_s2z_{id}, r_s2z_{id}');\n",
    "    vector[M_{id}] mean_white_s2z = mdivide_left_tri_low_brms(",
    "L_Sigma_s2z_{id}, mhat_s2z_{id});\n",
    "{direct_group_quad_code}",
    "  }}\n",
    cglue("  r_s2z_{idp}_{r$cn} = r_s2z_{id}[, {J}];\n")
  )
  str_add(out$pll_args) <- cglue(
    ", vector r_s2z_{idp}_{r$cn}"
  )

  # Score the original conditional Gaussian hierarchy at its completed-square
  # mode, followed by the exact Gaussian integration factor. The final log(N)
  # is the measure correction for the orthonormal contrast coordinates.
  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    if (spec$dist == "flat") {
      next
    }
    if (spec$dist == "normal") {
      cond_scale <- stan_s2z_number(spec$scale)
    } else {
      cond_scale <- glue(
        "{stan_s2z_number(spec$scale)} * sqrt(",
        "{stan_s2z_number(spec$df)} * udf_b_s2z{p}_{k})"
      )
    }
    str_add(out$tpar_prior) <- glue(
      "  lprior += normal_{lpdf}(qhat_s2z_{id}[{k}] | ",
      "{stan_s2z_number(spec$location)}, {cond_scale});\n"
    )
  }
  str_add(out$tpar_prior) <- glue(
    "  lprior += -0.5 * group_quad_s2z_{id}\n",
    "    - N_{id} * sum(log(diagonal(L_Sigma_s2z_{id})))\n",
    "    - M_{id} * sum(log(group_scale_s2z_{id}))\n",
    "    - sum(log(diagonal(L_P_s2z_{id})))",
    str_if(
      normalize,
      glue(" - 0.5 * N_{id} * M_{id} * log(2 * pi())",
           " + 0.5 * M_{id} * log(2 * pi())",
           " + 0.5 * M_{id} * log(1.0 * N_{id})")
    ),
    ";\n"
  )

  # Recover conventional super-population coefficients and deviations using
  # one shared conditional draw of the analytically omitted group mean.
  str_add(out$gen_def) <- glue(
    "  vector[M_{id}] mean_r_s2z_{id};\n",
    "  vector[{q}] q_recovered_s2z_{id};\n",
    "  matrix[N_{id}, M_{id}] r_{id};\n"
  )
  str_add(out$gen_def) <- stan_re_s2z_public_fe_def(info)
  str_add(out$gen_def) <- cglue(
    "  vector[N_{id}] r_{idp}_{r$cn};\n"
  )
  if (is_cor) {
    str_add(out$gen_def) <- glue(
      "  // compute group-level correlations\n",
      "  corr_matrix[M_{id}] Cor_{id}",
      " = multiply_lower_tri_self_transpose(L_{id});\n",
      "  vector<lower=-1,upper=1>[NC_{id}] cor_{id};\n"
    )
  }

  str_add(out$gen_comp) <- glue(
    "  {{\n",
    "    vector[M_{id}] z_mean_s2z;\n",
    "    for (k in 1:M_{id}) z_mean_s2z[k] = std_normal_rng();\n",
    "    mean_r_s2z_{id} = mhat_s2z_{id} + ",
    "(mdivide_right_tri_low_brms(z_mean_s2z', L_P_s2z_{id}))';\n",
    "  }}\n",
    "  q_recovered_s2z_{id} = theta_s2z{p} - H_s2z_{id} * mean_r_s2z_{id};\n"
  )
  str_add(out$gen_comp) <- stan_re_s2z_public_re_comp(r, part = "matrix")
  str_add(out$gen_comp) <- stan_re_s2z_public_fe_comp(info, id)
  str_add(out$gen_comp) <- stan_re_s2z_public_re_comp(r, part = "columns")
  if (is_cor) {
    str_add(out$gen_comp) <- stan_cor_gen_comp(
      cor = glue("cor_{id}"), ncol = glue("M_{id}")
    )
  }
  out
}

# Component-wise specialization for multiple independent group-level effects.
# The omitted means are conditionally independent except for the single
# rank-one coupling induced by brms's centered population intercept.  A
# diagonal-plus-rank-one solve therefore replaces all M x M factorizations.
.stan_re_s2z_independent <- function(id, r, info, normalize, out = list()) {
  lpdf <- stan_lpdf_name(normalize)
  stopifnot(
    is.reframe(r), nrow(r) > 1L, !isTRUE(r$cor[1]), all(r$s2z)
  )
  q <- length(info$qnames)
  M <- nrow(r)
  J <- seq_rows(r)
  p <- info$p
  idp <- paste0(r$id, usc(combine_prefix(check_prefix(r))))
  is_student <- identical(r$dist[1], "student")
  r_s2z <- glue("r_s2z_{idp}_{r$cn}")
  r_public <- glue("r_{idp}_{r$cn}")

  str_add(out$fun) <- "  #include 'fun_sum_to_zero.stan'\n"
  str_add(out$par) <- glue(
    "  vector[M_{id} * (N_{id} - 1)] z_s2z_{id};",
    "  // physical orthonormal independent S2Z coordinates\n"
  )

  # Independent Student-t and Cauchy population priors remain conditionally
  # Gaussian after introducing one scalar mixing variable per coefficient.
  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    if (spec$dist == "student") {
      str_add(out$par) <- glue(
        "  real<lower=0> udf_b_s2z{p}_{k};",
        "  // mixing variable for population coefficient {k}\n"
      )
      str_add(out$tpar_prior) <- glue(
        "  lprior += inv_chi_square_{lpdf}(udf_b_s2z{p}_{k} | ",
        "{stan_s2z_number(spec$df)});\n"
      )
    }
  }

  str_add(out$tpar_def) <- glue(
    "  // component-wise physical S2Z effects of ID {id}\n",
    cglue("  vector[N_{id}] {r_s2z};\n"),
    "  vector[{q}] prior_mean_s2z_{id};\n",
    "  vector<lower=0>[{q}] prior_prec_s2z_{id};\n",
    str_if(
      is_student,
      glue(
        "  vector<lower=0>[N_{id}] group_scale_s2z_{id};\n",
        "  vector<lower=0>[N_{id}] group_prec_s2z_{id};\n"
      )
    ),
    "  vector<lower=0>[M_{id}] D_diag_s2z_{id};\n",
    "  real<lower=0> rank1_info_s2z_{id};\n",
    "  vector[M_{id}] mhat_s2z_{id};\n",
    "  vector[{q}] qhat_s2z_{id};\n",
    "  real<lower=0> group_quad_s2z_{id};\n"
  )

  # The centered-intercept coupling is determined entirely by the design.
  str_add(out$tdata_def) <- glue(
    "  vector[M_{id}] intercept_map_s2z_{id};\n"
  )
  str_add(out$tdata_comp) <- glue(
    "  intercept_map_s2z_{id} = zeros_vector(M_{id});\n"
  )
  if (info$center) {
    for (j in seq_len(M)) {
      qi <- info$match_q[j]
      value <- if (info$r$coef[j] == "Intercept") {
        "1.0"
      } else {
        glue("means_X{p}[{qi - 1L}]")
      }
      str_add(out$tdata_comp) <- glue(
        "  intercept_map_s2z_{id}[{j}] = {value};\n"
      )
    }
  }

  for (j in seq_len(M)) {
    str_add(out$tpar_comp) <- glue(
      "  {r_s2z[j]} = sum_to_zero_constrain_brms(",
      "segment(z_s2z_{id}, ({j} - 1) * (N_{id} - 1) + 1, ",
      "N_{id} - 1));\n"
    )
  }

  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    loc <- stan_s2z_number(spec$location)
    str_add(out$tpar_comp) <- glue(
      "  prior_mean_s2z_{id}[{k}] = {loc};\n"
    )
    if (spec$dist == "flat") {
      prec <- "0.0"
    } else if (spec$dist == "normal") {
      prec <- glue("inv_square({stan_s2z_number(spec$scale)})")
    } else {
      prec <- glue(
        "inv_square({stan_s2z_number(spec$scale)} * sqrt(",
        "{stan_s2z_number(spec$df)} * udf_b_s2z{p}_{k}))"
      )
    }
    str_add(out$tpar_comp) <- glue(
      "  prior_prec_s2z_{id}[{k}] = {prec};\n"
    )
  }

  if (is_student) {
    tr <- subset_reframe_dist(r, "student")
    g <- usc(tr$ggn[1])
    str_add(out$tpar_comp) <- glue(
      "  group_scale_s2z_{id} = dfm{g};\n",
      "  group_prec_s2z_{id} = inv_square(group_scale_s2z_{id});\n"
    )
  }

  # After multiplying each normal equation by its group SD squared, its
  # diagonal is D.  The only remaining coupling is lambda_Intercept * v v',
  # where v maps omitted means into the centered intercept.  Sherman-Morrison
  # gives both the mode and determinant in linear time and avoids unstable
  # inverse-square group SDs during warmup.
  str_add(out$tpar_comp) <- glue(
    "  {{\n",
    "    vector[M_{id}] base_info_s2z = zeros_vector(M_{id});\n",
    "    vector[M_{id}] base_score_s2z = zeros_vector(M_{id});\n",
    "    vector[M_{id}] scaled_score_s2z;\n",
    "    vector[M_{id}] independent_mode_s2z;\n",
    "    real group_info_s2z = ",
    str_if(is_student, glue("sum(group_prec_s2z_{id})"), glue("N_{id}")),
    ";\n"
  )
  for (j in seq_len(M)) {
    qi <- info$match_q[j]
    if (!info$center || qi != 1L) {
      str_add(out$tpar_comp) <- glue(
        "    base_info_s2z[{j}] = prior_prec_s2z_{id}[{qi}];\n",
        "    base_score_s2z[{j}] = prior_prec_s2z_{id}[{qi}] * ",
        "(theta_s2z{p}[{qi}] - prior_mean_s2z_{id}[{qi}]);\n"
      )
    }
  }
  str_add(out$tpar_comp) <- glue(
    "    D_diag_s2z_{id} = group_info_s2z + ",
    "square(sd_{id}) .* base_info_s2z;\n"
  )
  for (j in seq_len(M)) {
    contrast_score <- if (is_student) {
      glue(" - dot_product({r_s2z[j]}, group_prec_s2z_{id})")
    } else {
      ""
    }
    str_add(out$tpar_comp) <- glue(
      "    scaled_score_s2z[{j}] = square(sd_{id}[{j}]) * ",
      "base_score_s2z[{j}]{contrast_score};\n"
    )
  }
  coupling_prec <- str_if(info$center, glue("prior_prec_s2z_{id}[1]"), "0.0")
  coupling_resid <- str_if(
    info$center,
    glue("theta_s2z{p}[1] - prior_mean_s2z_{id}[1]"),
    "0.0"
  )
  str_add(out$tpar_comp) <- glue(
    "    independent_mode_s2z = scaled_score_s2z ./ D_diag_s2z_{id};\n",
    "    rank1_info_s2z_{id} = {coupling_prec} * dot_product(\n",
    "      square(sd_{id}) .* square(intercept_map_s2z_{id}),\n",
    "      1.0 ./ D_diag_s2z_{id}\n",
    "    );\n",
    "    mhat_s2z_{id} = independent_mode_s2z +\n",
    "      {coupling_prec} * square(sd_{id}) .* intercept_map_s2z_{id} ./\n",
    "      D_diag_s2z_{id} * ({coupling_resid} -\n",
    "      dot_product(intercept_map_s2z_{id}, independent_mode_s2z)) /\n",
    "      (1.0 + rank1_info_s2z_{id});\n",
    "  }}\n",
    "  qhat_s2z_{id} = theta_s2z{p};\n"
  )
  if (info$center) {
    str_add(out$tpar_comp) <- glue(
      "  qhat_s2z_{id}[1] -= dot_product(",
      "intercept_map_s2z_{id}, mhat_s2z_{id});\n"
    )
    for (j in seq_len(M)) {
      qi <- info$match_q[j]
      if (qi != 1L) {
        str_add(out$tpar_comp) <- glue(
          "  qhat_s2z_{id}[{qi}] -= mhat_s2z_{id}[{j}];\n"
        )
      }
    }
  } else {
    for (j in seq_len(M)) {
      str_add(out$tpar_comp) <- glue(
        "  qhat_s2z_{id}[{info$match_q[j]}] -= mhat_s2z_{id}[{j}];\n"
      )
    }
  }
  str_add(out$tpar_comp) <- glue(
    "  group_quad_s2z_{id} = 0.0;\n"
  )
  for (j in seq_len(M)) {
    standardized <- glue(
      "({r_s2z[j]} + mhat_s2z_{id}[{j}]) / sd_{id}[{j}]",
      str_if(is_student, glue(" ./ group_scale_s2z_{id}"))
    )
    str_add(out$tpar_comp) <- glue(
      "  group_quad_s2z_{id} += dot_self({standardized});\n"
    )
  }
  str_add(out$pll_args) <- cglue(", vector {r_s2z}")

  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    if (spec$dist == "flat") {
      next
    }
    if (spec$dist == "normal") {
      cond_scale <- stan_s2z_number(spec$scale)
    } else {
      cond_scale <- glue(
        "{stan_s2z_number(spec$scale)} * sqrt(",
        "{stan_s2z_number(spec$df)} * udf_b_s2z{p}_{k})"
      )
    }
    str_add(out$tpar_prior) <- glue(
      "  lprior += normal_{lpdf}(qhat_s2z_{id}[{k}] | ",
      "{stan_s2z_number(spec$location)}, {cond_scale});\n"
    )
  }
  str_add(out$tpar_prior) <- glue(
    "  lprior += -0.5 * group_quad_s2z_{id}\n",
    "    - (N_{id} - 1) * sum(log(sd_{id}))\n",
    str_if(
      is_student,
      glue("    - M_{id} * sum(log(group_scale_s2z_{id}))\n")
    ),
    "    - 0.5 * sum(log(D_diag_s2z_{id}))\n",
    "    - 0.5 * log1p(rank1_info_s2z_{id})",
    str_if(
      normalize,
      glue(" - 0.5 * N_{id} * M_{id} * log(2 * pi())",
           " + 0.5 * M_{id} * log(2 * pi())",
           " + 0.5 * M_{id} * log(1.0 * N_{id})")
    ),
    ";\n"
  )

  # Draw from the diagonal-minus-rank-one conditional covariance without a
  # dense Cholesky, then reconstruct conventional public b/r coordinates.
  str_add(out$gen_def) <- glue(
    "  vector[M_{id}] mean_r_s2z_{id};\n",
    "  vector[{q}] q_recovered_s2z_{id};\n"
  )
  str_add(out$gen_def) <- stan_re_s2z_public_fe_def(info)
  str_add(out$gen_def) <- cglue("  vector[N_{id}] {r_public};\n")

  str_add(out$gen_comp) <- glue(
    "  {{\n",
    "    vector[M_{id}] independent_noise_s2z;\n",
    "    real sqrt_rank1_s2z = sqrt(1.0 + rank1_info_s2z_{id});\n",
    "    real rank1_adjust_s2z;\n",
    "    for (k in 1:M_{id}) {{\n",
    "      independent_noise_s2z[k] = sd_{id}[k] * ",
    "std_normal_rng() / sqrt(D_diag_s2z_{id}[k]);\n",
    "    }}\n",
    "    rank1_adjust_s2z = {coupling_prec} / ",
    "(sqrt_rank1_s2z * (1.0 + sqrt_rank1_s2z)) *\n",
    "      dot_product(intercept_map_s2z_{id}, independent_noise_s2z);\n",
    "    mean_r_s2z_{id} = mhat_s2z_{id} + independent_noise_s2z -\n",
    "      rank1_adjust_s2z * square(sd_{id}) .* ",
    "intercept_map_s2z_{id} ./ D_diag_s2z_{id};\n",
    "  }}\n",
    "  q_recovered_s2z_{id} = theta_s2z{p};\n"
  )
  if (info$center) {
    str_add(out$gen_comp) <- glue(
      "  q_recovered_s2z_{id}[1] -= dot_product(",
      "intercept_map_s2z_{id}, mean_r_s2z_{id});\n"
    )
    for (j in seq_len(M)) {
      qi <- info$match_q[j]
      if (qi != 1L) {
        str_add(out$gen_comp) <- glue(
          "  q_recovered_s2z_{id}[{qi}] -= mean_r_s2z_{id}[{j}];\n"
        )
      }
    }
  } else {
    for (j in seq_len(M)) {
      str_add(out$gen_comp) <- glue(
        "  q_recovered_s2z_{id}[{info$match_q[j]}] -= ",
        "mean_r_s2z_{id}[{j}];\n"
      )
    }
  }
  str_add(out$gen_comp) <- stan_re_s2z_public_re_comp(r)
  str_add(out$gen_comp) <- stan_re_s2z_public_fe_comp(info, id)
  out
}

# Scalar specialization of the physical sum-to-zero block.
# In addition to avoiding all 1 x 1 matrix algebra, the Gaussian branch uses
# the exact zero sum of the orthonormal contrasts to remove a weighted dot
# product and the group-scale vectors altogether.
.stan_re_s2z_scalar <- function(id, r, info, normalize, out = list()) {
  lpdf <- stan_lpdf_name(normalize)
  stopifnot(is.reframe(r), nrow(r) == 1L, all(r$s2z))
  q <- length(info$qnames)
  p <- info$p
  idp <- paste0(r$id, usc(combine_prefix(check_prefix(r))))
  cn <- r$cn[1]
  is_student <- identical(r$dist[1], "student")
  r_s2z <- glue("r_s2z_{idp}_{cn}")
  r_public <- glue("r_{idp}_{cn}")

  str_add(out$fun) <- "  #include 'fun_sum_to_zero.stan'\n"
  str_add(out$par) <- glue(
    "  vector[N_{id} - 1] z_s2z_{id};",
    "  // physical orthonormal scalar S2Z coordinates\n"
  )

  # Independent Student-t and Cauchy population priors remain conditionally
  # Gaussian after introducing one scalar mixing variable per coefficient.
  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    if (spec$dist == "student") {
      str_add(out$par) <- glue(
        "  real<lower=0> udf_b_s2z{p}_{k};",
        "  // mixing variable for population coefficient {k}\n"
      )
      str_add(out$tpar_prior) <- glue(
        "  lprior += inv_chi_square_{lpdf}(udf_b_s2z{p}_{k} | ",
        "{stan_s2z_number(spec$df)});\n"
      )
    }
  }

  str_add(out$tpar_def) <- glue(
    "  // specialized scalar physical S2Z effects of ID {id}\n",
    "  vector[N_{id}] {r_s2z};\n",
    "  vector[{q}] prior_mean_s2z_{id};\n",
    "  vector<lower=0>[{q}] prior_prec_s2z_{id};\n",
    str_if(
      is_student,
      glue(
        "  vector<lower=0>[N_{id}] group_scale_s2z_{id};\n",
        "  vector<lower=0>[N_{id}] group_prec_s2z_{id};\n"
      )
    ),
    "  real<lower=0> D_s2z_{id};\n",
    "  real<lower=0> sqrt_D_s2z_{id};\n",
    "  real mhat_s2z_{id};\n",
    "  vector[{q}] qhat_s2z_{id};\n",
    "  real<lower=0> group_quad_s2z_{id};\n"
  )

  # H maps the omitted scalar group mean into all matching population
  # coordinates. A centered varying slope also shifts brms's temporary
  # centered intercept by the mean of its raw design column.
  qi <- info$match_q[1]
  str_add(out$tdata_def) <- glue(
    "  vector[{q}] H_s2z_{id};\n"
  )
  str_add(out$tdata_comp) <- glue(
    "  H_s2z_{id} = zeros_vector({q});\n",
    "  H_s2z_{id}[{qi}] = 1.0;\n",
    str_if(
      info$center && info$r$coef[1] != "Intercept",
      glue("  H_s2z_{id}[1] = means_X{p}[{qi - 1L}];\n")
    )
  )
  str_add(out$tpar_comp) <- glue(
    "  {r_s2z} = sum_to_zero_constrain_brms(z_s2z_{id});\n"
  )

  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    loc <- stan_s2z_number(spec$location)
    str_add(out$tpar_comp) <- glue(
      "  prior_mean_s2z_{id}[{k}] = {loc};\n"
    )
    if (spec$dist == "flat") {
      prec <- "0.0"
    } else if (spec$dist == "normal") {
      prec <- glue("inv_square({stan_s2z_number(spec$scale)})")
    } else {
      prec <- glue(
        "inv_square({stan_s2z_number(spec$scale)} * sqrt(",
        "{stan_s2z_number(spec$df)} * udf_b_s2z{p}_{k}))"
      )
    }
    str_add(out$tpar_comp) <- glue(
      "  prior_prec_s2z_{id}[{k}] = {prec};\n"
    )
  }

  # D = square(sd) * P is the scaled conditional precision of the omitted
  # mean. Working with D avoids forming inv_square(sd), which can overflow
  # near zero during warmup. For Gaussian effects, sum(r_s2z) is exactly zero
  # in the physical coordinates, so its contribution to the mode vanishes.
  if (is_student) {
    tr <- subset_reframe_dist(r, "student")
    g <- usc(tr$ggn[1])
    str_add(out$tpar_comp) <- glue(
      "  group_scale_s2z_{id} = dfm{g};\n",
      "  group_prec_s2z_{id} = inv_square(group_scale_s2z_{id});\n",
      "  {{\n",
      "    real tau_sq_s2z = square(sd_{id}[1]);\n",
      "    real prior_info_s2z = dot_product(prior_prec_s2z_{id}, ",
      "square(H_s2z_{id}));\n",
      "    real prior_score_s2z = dot_product(H_s2z_{id}, ",
      "prior_prec_s2z_{id} .* (theta_s2z{p} - ",
      "prior_mean_s2z_{id}));\n",
      "    D_s2z_{id} = tau_sq_s2z * prior_info_s2z + ",
      "sum(group_prec_s2z_{id});\n",
      "    mhat_s2z_{id} = (tau_sq_s2z * prior_score_s2z - ",
      "dot_product({r_s2z}, group_prec_s2z_{id})) / D_s2z_{id};\n",
      "  }}\n"
    )
  } else {
    str_add(out$tpar_comp) <- glue(
      "  {{\n",
      "    real tau_sq_s2z = square(sd_{id}[1]);\n",
      "    real prior_info_s2z = dot_product(prior_prec_s2z_{id}, ",
      "square(H_s2z_{id}));\n",
      "    real prior_score_s2z = dot_product(H_s2z_{id}, ",
      "prior_prec_s2z_{id} .* (theta_s2z{p} - ",
      "prior_mean_s2z_{id}));\n",
      "    D_s2z_{id} = tau_sq_s2z * prior_info_s2z + N_{id};\n",
      "    mhat_s2z_{id} = tau_sq_s2z * prior_score_s2z / D_s2z_{id};\n",
      "  }}\n"
    )
  }
  str_add(out$tpar_comp) <- glue(
    "  sqrt_D_s2z_{id} = sqrt(D_s2z_{id});\n",
    "  qhat_s2z_{id} = theta_s2z{p} - H_s2z_{id} * ",
    "mhat_s2z_{id};\n",
    "  {{\n",
    "    vector[N_{id}] white_s2z = ({r_s2z} + mhat_s2z_{id}) / ",
    "sd_{id}[1]",
    str_if(is_student, glue(" ./ group_scale_s2z_{id}")),
    ";\n",
    "    group_quad_s2z_{id} = dot_self(white_s2z);\n",
    "  }}\n"
  )
  str_add(out$pll_args) <- glue(", vector {r_s2z}")

  for (k in seq_len(q)) {
    spec <- info$prior[[k]]
    if (spec$dist == "flat") {
      next
    }
    if (spec$dist == "normal") {
      cond_scale <- stan_s2z_number(spec$scale)
    } else {
      cond_scale <- glue(
        "{stan_s2z_number(spec$scale)} * sqrt(",
        "{stan_s2z_number(spec$df)} * udf_b_s2z{p}_{k})"
      )
    }
    str_add(out$tpar_prior) <- glue(
      "  lprior += normal_{lpdf}(qhat_s2z_{id}[{k}] | ",
      "{stan_s2z_number(spec$location)}, {cond_scale});\n"
    )
  }
  str_add(out$tpar_prior) <- glue(
    "  lprior += -0.5 * group_quad_s2z_{id}\n",
    "    - (N_{id} - 1) * log(sd_{id}[1])\n",
    str_if(
      is_student,
      glue("    - sum(log(group_scale_s2z_{id}))\n")
    ),
    "    - 0.5 * log(D_s2z_{id})",
    str_if(
      normalize,
      glue(" - 0.5 * N_{id} * log(2 * pi())",
           " + 0.5 * log(2 * pi())",
           " + 0.5 * log(1.0 * N_{id})")
    ),
    ";\n"
  )

  # Reconstruct the conventional coefficient and deviations from one draw
  # of the analytically omitted scalar group mean.
  str_add(out$gen_def) <- glue(
    "  real mean_r_s2z_{id};\n",
    "  vector[{q}] q_recovered_s2z_{id};\n"
  )
  str_add(out$gen_def) <- stan_re_s2z_public_fe_def(info)
  str_add(out$gen_def) <- glue("  vector[N_{id}] {r_public};\n")

  str_add(out$gen_comp) <- glue(
    "  mean_r_s2z_{id} = mhat_s2z_{id} + ",
    "sd_{id}[1] * std_normal_rng() / sqrt_D_s2z_{id};\n",
    "  q_recovered_s2z_{id} = theta_s2z{p} - H_s2z_{id} * ",
    "mean_r_s2z_{id};\n"
  )
  str_add(out$gen_comp) <- stan_re_s2z_public_re_comp(r, scalar_mean = TRUE)
  str_add(out$gen_comp) <- stan_re_s2z_public_fe_comp(info, id)
  out
}

# Stable formatting for numeric constants inserted into generated Stan code.
stan_s2z_number <- function(x) {
  stopifnot(length(x) == 1L, is.finite(x))
  trimws(formatC(x, digits = 17, format = "g"))
}
