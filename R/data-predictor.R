#' Prepare Predictor Data
#'
#' Prepare data related to predictor variables in \pkg{brms}.
#' Only exported for use in package development.
#'
#' @param x An \R object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A named list of data related to predictor variables.
#'
#' @keywords internal
#' @export
data_predictor <- function(x, ...) {
  UseMethod("data_predictor")
}

#' @export
data_predictor.mvbrmsterms <- function(x, data, basis = NULL, ...) {
  out <- list(N = nrow(data))
  for (r in names(x$terms)) {
    bs <- basis$resps[[r]]
    c(out) <- data_predictor(x$terms[[r]], data = data, basis = bs, ...)
  }
  out
}

#' @export
data_predictor.brmsterms <- function(x, data, data2, prior, ranef,
                                     basis = NULL, ...) {
  out <- list()
  data <- subset_data(data, x)
  resp <- usc(combine_prefix(x))
  args_eff <- nlist(data, data2, ranef, prior, ...)
  for (dp in names(x$dpars)) {
    args_eff_spec <- list(x = x$dpars[[dp]], basis = basis$dpars[[dp]])
    c(out) <- do_call(data_predictor, c(args_eff_spec, args_eff))
  }
  for (dp in names(x$fdpars)) {
    if (is.numeric(x$fdpars[[dp]]$value)) {
      out[[paste0(dp, resp)]] <- x$fdpars[[dp]]$value
    }
  }
  for (nlp in names(x$nlpars)) {
    args_eff_spec <- list(x = x$nlpars[[nlp]], basis = basis$nlpars[[nlp]])
    c(out) <- do_call(data_predictor, c(args_eff_spec, args_eff))
  }
  c(out) <- data_gr_local(x, data = data, ranef = ranef)
  c(out) <- data_mixture(x, data2 = data2, prior = prior)
  out
}

# prepare data for all types of effects for use in Stan
# @param data the data passed by the user
# @param ranef object retuend by 'tidy_ranef'
# @param prior an object of class brmsprior
# @param basis information from original Stan data used to correctly
#   predict from new data. See 'standata_basis' for details.
# @param ... currently ignored
# @return a named list of data to be passed to Stan
#' @export
data_predictor.btl <- function(x, data, data2 = list(), ranef = empty_ranef(),
                               prior = brmsprior(), index = NULL, basis = NULL,
                               ...) {
  out <- c(
    data_fe(x, data),
    data_sp(x, data, data2 = data2, prior = prior, index = index, basis = basis$sp),
    data_re(x, data, ranef = ranef),
    data_cs(x, data),
    data_sm(x, data, basis = basis$sm),
    data_gp(x, data, basis = basis$gp),
    data_ac(x, data, data2 = data2, basis = basis$ac),
    data_offset(x, data),
    data_bhaz(x, data, data2 = data2, prior = prior, basis = basis$bhaz)
  )
  c(out) <- data_prior(x, data, prior = prior, sdata = out)
  out
}

# prepare data for non-linear parameters for use in Stan
#' @export
data_predictor.btnl <- function(x, data, data2 = list(), prior = brmsprior(),
                                basis = NULL, ...) {
  out <- list()
  c(out) <- data_cnl(x, data)
  c(out) <- data_ac(x, data, data2 = data2, basis = basis$ac)
  c(out) <- data_bhaz(x, data, data2 = data2, prior = prior, basis = basis$bhaz)
  out
}

# prepare data of fixed effects
data_fe <- function(bterms, data) {
  out <- list()
  p <- usc(combine_prefix(bterms))
  # the intercept is removed inside the Stan code for ordinal models
  cols2remove <- if (is_ordinal(bterms)) "(Intercept)"
  X <- get_model_matrix(rhs(bterms$fe), data, cols2remove = cols2remove)
  avoid_dpars(colnames(X), bterms = bterms)
  out[[paste0("K", p)]] <- ncol(X)
  out[[paste0("X", p)]] <- X
  out
}

# data preparation for splines
data_sm <- function(bterms, data, basis = NULL) {
  out <- list()
  smterms <- all_terms(bterms[["sm"]])
  if (!length(smterms)) {
    return(out)
  }
  p <- usc(combine_prefix(bterms))
  new <- length(basis) > 0L
  if (!new) {
    knots <- get_knots(data)
    basis <- named_list(smterms)
    for (i in seq_along(smterms)) {
      # the spline penalty has changed in 2.8.7 (#646)
      diagonal.penalty <- !require_old_default("2.8.7")
      basis[[i]] <- smoothCon(
        eval2(smterms[i]), data = data,
        knots = knots, absorb.cons = TRUE,
        diagonal.penalty = diagonal.penalty
      )
    }
  }
  bylevels <- named_list(smterms)
  ns <- 0
  lXs <- list()
  for (i in seq_along(basis)) {
    # may contain multiple terms when 'by' is a factor
    for (j in seq_along(basis[[i]])) {
      ns <- ns + 1
      sm <- basis[[i]][[j]]
      if (length(sm$by.level)) {
        bylevels[[i]][j] <- sm$by.level
      }
      if (new) {
        # prepare rasm for use with new data
        rasm <- s2rPred(sm, data)
      } else {
        rasm <- mgcv::smooth2random(sm, names(data), type = 2)
      }
      lXs[[ns]] <- rasm$Xf
      if (NCOL(lXs[[ns]])) {
        colnames(lXs[[ns]]) <- paste0(sm$label, "_", seq_cols(lXs[[ns]]))
      }
      Zs <- rasm$rand
      sfx <- paste0(p, "_", ns)
      out[[paste0("nb", sfx)]] <- length(Zs)
      if (length(Zs)) {
        names(Zs) <- paste0("Zs", sfx, "_", seq_along(Zs))
        c(out) <- Zs
        out[[paste0("knots", sfx)]] <- as.array(ulapply(Zs, ncol))
      } else {
        out[[paste0("knots", sfx)]] <- integer(0)
      }
    }
  }
  Xs <- do_call(cbind, lXs)
  avoid_dpars(colnames(Xs), bterms = bterms)
  smcols <- lapply(lXs, function(x) which(colnames(Xs) %in% colnames(x)))
  Xs <- structure(Xs, smcols = smcols, bylevels = bylevels)
  colnames(Xs) <- rename(colnames(Xs))
  out[[paste0("Ks", p)]] <- ncol(Xs)
  out[[paste0("Xs", p)]] <- Xs
  out
}

# prepare data for group-level effects for use in Stan
data_re <- function(bterms, data, ranef) {
  out <- list()
  px <- check_prefix(bterms)
  take <- find_rows(ranef, ls = px) & !find_rows(ranef, type = "sp")
  ranef <- ranef[take, ]
  if (!nrow(ranef)) {
    return(out)
  }
  gn <- unique(ranef$gn)
  for (i in seq_along(gn)) {
    r <- subset2(ranef, gn = gn[i])
    Z <- get_model_matrix(r$form[[1]], data = data, rename = FALSE)
    idp <- paste0(r$id[1], usc(combine_prefix(px)))
    Znames <- paste0("Z_", idp, "_", r$cn)
    if (r$gtype[1] == "mm") {
      ng <- length(r$gcall[[1]]$groups)
      if (r$type[1] == "cs") {
        stop2("'cs' is not supported in multi-membership terms.")
      }
      if (r$type[1] == "mmc") {
        # see issue #353 for the general idea
        mmc_expr <- "^mmc\\([^:]*\\)"
        mmc_terms <- get_matches_expr(mmc_expr, colnames(Z))
        for (t in mmc_terms) {
          pos <- which(grepl_expr(escape_all(t), colnames(Z)))
          if (length(pos) != ng) {
            stop2("Invalid term '", t, "': Expected ", ng,
                  " coefficients but found ", length(pos), ".")
          }
          for (j in seq_along(Znames)) {
            for (k in seq_len(ng)) {
              out[[paste0(Znames[j], "_", k)]] <- as.array(Z[, pos[k]])
            }
          }
        }
      } else {
        for (j in seq_along(Znames)) {
          out[paste0(Znames[j], "_", seq_len(ng))] <- list(as.array(Z[, j]))
        }
      }
    } else {
      if (r$type[1] == "cs") {
        ncatM1 <- nrow(r) / ncol(Z)
        Z_temp <- vector("list", ncol(Z))
        for (k in seq_along(Z_temp)) {
          Z_temp[[k]] <- replicate(ncatM1, Z[, k], simplify = FALSE)
        }
        Z <- do_call(cbind, unlist(Z_temp, recursive = FALSE))
      }
      if (r$type[1] == "mmc") {
        stop2("'mmc' is only supported in multi-membership terms.")
      }
      for (j in seq_cols(Z)) {
        out[[Znames[j]]] <- as.array(Z[, j])
      }
    }
  }
  out
}

# compute data for each group-level-ID per univariate model
data_gr_local <- function(bterms, data, ranef) {
  stopifnot(is.brmsterms(bterms))
  out <- list()
  ranef <- subset2(ranef, resp = bterms$resp)
  resp <- usc(bterms$resp)
  for (id in unique(ranef$id)) {
    id_ranef <- subset2(ranef, id = id)
    idresp <- paste0(id, resp)
    nranef <- nrow(id_ranef)
    group <- id_ranef$group[1]
    levels <- attr(ranef, "levels")[[group]]
    if (id_ranef$gtype[1] == "mm") {
      # multi-membership grouping term
      gs <- id_ranef$gcall[[1]]$groups
      ngs <- length(gs)
      weights <- id_ranef$gcall[[1]]$weights
      if (is.formula(weights)) {
        scale <- isTRUE(attr(weights, "scale"))
        weights <- as.matrix(eval_rhs(weights, data))
        if (!identical(dim(weights), c(nrow(data), ngs))) {
          stop2(
            "Grouping structure 'mm' expects 'weights' to be ",
            "a matrix with as many columns as grouping factors."
          )
        }
        if (scale) {
          if (isTRUE(any(weights < 0))) {
            stop2("Cannot scale negative weights.")
          }
          weights <- sweep(weights, 1, rowSums(weights), "/")
        }
      } else {
        # all members get equal weights by default
        weights <- matrix(1 / ngs, nrow = nrow(data), ncol = ngs)
      }
      for (i in seq_along(gs)) {
        gdata <- get(gs[i], data)
        J <- match(gdata, levels)
        if (anyNA(J)) {
          # occurs for new levels only
          new_gdata <- gdata[!gdata %in% levels]
          new_levels <- unique(new_gdata)
          J[is.na(J)] <- match(new_gdata, new_levels) + length(levels)
        }
        out[[paste0("J_", idresp, "_", i)]] <- as.array(J)
        out[[paste0("W_", idresp, "_", i)]] <- as.array(weights[, i])
      }
    } else {
      # ordinary grouping term
      g <- id_ranef$gcall[[1]]$groups
      gdata <- get(g, data)
      J <- match(gdata, levels)
      if (anyNA(J)) {
        # occurs for new levels only
        new_gdata <- gdata[!gdata %in% levels]
        new_levels <- unique(new_gdata)
        J[is.na(J)] <- match(new_gdata, new_levels) + length(levels)
      }
      out[[paste0("J_", idresp)]] <- as.array(J)
    }
  }
  out
}

# prepare global data for each group-level-ID
data_gr_global <- function(ranef, data2) {
  out <- list()
  for (id in unique(ranef$id)) {
    tmp <- list()
    id_ranef <- subset2(ranef, id = id)
    nranef <- nrow(id_ranef)
    group <- id_ranef$group[1]
    levels <- attr(ranef, "levels")[[group]]
    tmp$N <- length(levels)
    tmp$M <- nranef
    tmp$NC <- as.integer(nranef * (nranef - 1) / 2)
    # prepare number of levels of an optional 'by' variable
    if (nzchar(id_ranef$by[1])) {
      stopifnot(!nzchar(id_ranef$type[1]))
      bylevels <- id_ranef$bylevels[[1]]
      Jby <- match(attr(levels, "by"), bylevels)
      tmp$Nby <- length(bylevels)
      tmp$Jby <- as.array(Jby)
    }
    # prepare within-group covariance matrices
    cov <- id_ranef$cov[1]
    if (nzchar(cov)) {
      # validation is only necessary here for compatibility with 'cov_ranef'
      cov_mat <- validate_recov_matrix(data2[[cov]])
      found_levels <- rownames(cov_mat)
      found <- levels %in% found_levels
      if (any(!found)) {
        stop2("Levels of the within-group covariance matrix for '", group,
              "' do not match names of the grouping levels.")
      }
      cov_mat <- cov_mat[levels, levels, drop = FALSE]
      tmp$Lcov <- t(chol(cov_mat))
    }
    names(tmp) <- paste0(names(tmp), "_", id)
    c(out) <- tmp
  }
  out
}

# prepare data for special effects for use in Stan
data_sp <- function(bterms, data, data2, prior, index = NULL, basis = NULL) {
  out <- list()
  spef <- tidy_spef(bterms, data)
  if (!nrow(spef)) return(out)
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  # prepare general data
  out[[paste0("Ksp", p)]] <- nrow(spef)
  Csp <- sp_model_matrix(bterms$sp, data)
  avoid_dpars(colnames(Csp), bterms = bterms)
  Csp <- Csp[, spef$Ic > 0, drop = FALSE]
  Csp <- lapply(seq_cols(Csp), function(i) as.array(Csp[, i]))
  if (length(Csp)) {
    Csp_names <- paste0("Csp", p, "_", seq_along(Csp))
    out <- c(out, setNames(Csp, Csp_names))
  }
  if (any(lengths(spef$Imo) > 0)) {
    # prepare data specific to monotonic effects
    out[[paste0("Imo", p)]] <- max(unlist(spef$Imo))
    Xmo <- lapply(unlist(spef$calls_mo), get_mo_values, data = data)
    Xmo_names <- paste0("Xmo", p, "_", seq_along(Xmo))
    c(out) <- setNames(Xmo, Xmo_names)
    if (!is.null(basis$Jmo)) {
      # take information from original data
      Jmo <- basis$Jmo
    } else {
      Jmo <- as.array(ulapply(Xmo, max))
    }
    out[[paste0("Jmo", p)]] <- Jmo
    # prepare prior concentration of simplex parameters
    simo_coef <- get_simo_labels(spef, use_id = TRUE)
    ids <- unlist(spef$ids_mo)
    for (j in seq_along(simo_coef)) {
      # index of first ID appearance
      j_id <- match(ids[j], ids)
      if (is.na(ids[j]) || j_id == j) {
        # only evaluate priors without ID or first appearance of the ID
        # all other parameters will be copied over in the Stan code
        simo_prior <- subset2(prior,
          class = "simo", coef = simo_coef[j], ls = px
        )
        con_simo <- eval_dirichlet(simo_prior$prior, Jmo[j], data2)
        out[[paste0("con_simo", p, "_", j)]] <- as.array(con_simo)
      }
    }
  }
  uni_mi <- attr(spef, "uni_mi")
  for (j in seq_rows(uni_mi)) {
    if (!is.na(uni_mi$idx[j])) {
      idxl <- get(uni_mi$idx[j], data)
      if (is.null(index[[uni_mi$var[j]]])) {
        # the 'idx' argument needs to be mapped against 'index' addition terms
        stop2("Response '", uni_mi$var[j], "' needs to have an 'index' addition ",
              "term to compare with 'idx'. See ?mi for examples.")
      }
      idxl <- match(idxl, index[[uni_mi$var[j]]])
      if (anyNA(idxl)) {
        stop2("Could not match all indices in response '", uni_mi$var[j], "'.")
      }
      idxl_name <- paste0("idxl", p, "_", uni_mi$var[j], "_", uni_mi$idx2[j])
      out[[idxl_name]] <- as.array(idxl)
    } else if (isTRUE(attr(index[[uni_mi$var[j]]], "subset"))) {
      # cross-formula referencing is required for subsetted variables
      stop2("mi() terms of subsetted variables require ",
            "the 'idx' argument to be specified.")
    }
  }
  out
}

# prepare data for category specific effects
data_cs <- function(bterms, data) {
  out <- list()
  if (length(all_terms(bterms[["cs"]]))) {
    p <- usc(combine_prefix(bterms))
    Xcs <- get_model_matrix(bterms$cs, data)
    avoid_dpars(colnames(Xcs), bterms = bterms)
    out <- c(out, list(Kcs = ncol(Xcs), Xcs = Xcs))
    out <- setNames(out, paste0(names(out), p))
  }
  out
}

# prepare global data for noise free variables
data_Xme <- function(meef, data) {
  stopifnot(is.meef_frame(meef))
  out <- list()
  groups <- unique(meef$grname)
  for (i in seq_along(groups)) {
    g <- groups[i]
    K <- which(meef$grname %in% g)
    Mme <- length(K)
    out[[paste0("Mme_", i)]] <- Mme
    out[[paste0("NCme_", i)]] <- Mme * (Mme - 1) / 2
    if (nzchar(g)) {
      levels <- get_levels(meef)[[g]]
      gr <- get_me_group(meef$term[K[1]], data)
      Jme <- match(gr, levels)
      if (anyNA(Jme)) {
        # occurs for new levels only
        # replace NAs with unique values; fixes issue #706
        gr[is.na(gr)] <- paste0("new_", seq_len(sum(is.na(gr))), "__")
        new_gr <- gr[!gr %in% levels]
        new_levels <- unique(new_gr)
        Jme[is.na(Jme)] <- length(levels) + match(new_gr, new_levels)
      }
      ilevels <- unique(Jme)
      out[[paste0("Nme_", i)]] <- length(ilevels)
      out[[paste0("Jme_", i)]] <- Jme
    }
    for (k in K) {
      Xn <- get_me_values(meef$term[k], data)
      noise <- get_me_noise(meef$term[k], data)
      if (nzchar(g)) {
        for (l in ilevels) {
          # validate values of the same level
          take <- Jme %in% l
          if (length(unique(Xn[take])) > 1L ||
              length(unique(noise[take])) > 1L) {
            stop2(
              "Measured values and measurement error should be ",
              "unique for each group. Occured for level '",
              levels[l], "' of group '", g, "'."
            )
          }
        }
        Xn <- get_one_value_per_group(Xn, Jme)
        noise <- get_one_value_per_group(noise, Jme)
      }
      out[[paste0("Xn_", k)]] <- as.array(Xn)
      out[[paste0("noise_", k)]] <- as.array(noise)
    }
  }
  out
}

# prepare data for Gaussian process terms
# @param internal store some intermediate data for internal post-processing?
# @param ... passed to '.data_gp'
data_gp <- function(bterms, data, internal = FALSE, basis = NULL, ...) {
  out <- list()
  internal <- as_one_logical(internal)
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  gpef <- tidy_gpef(bterms, data)
  for (i in seq_rows(gpef)) {
    pi <- paste0(p, "_", i)
    Xgp <- lapply(gpef$covars[[i]], eval2, data)
    D <- length(Xgp)
    out[[paste0("Dgp", pi)]] <- D
    invalid <- ulapply(Xgp, function(x)
      !is.numeric(x) || isTRUE(length(dim(x)) > 1L)
    )
    if (any(invalid)) {
      stop2("Predictors of Gaussian processes should be numeric vectors.")
    }
    Xgp <- do_call(cbind, Xgp)
    cmc <- gpef$cmc[i]
    scale <- gpef$scale[i]
    gr <- gpef$gr[i]
    k <- gpef$k[i]
    c <- gpef$c[[i]]
    if (!isNA(k)) {
      out[[paste0("NBgp", pi)]] <- k ^ D
      Ks <- as.matrix(do_call(expand.grid, repl(seq_len(k), D)))
    }
    byvar <- gpef$byvars[[i]]
    byfac <- length(gpef$cons[[i]]) > 0L
    bynum <- !is.null(byvar) && !byfac
    if (byfac) {
      # for categorical 'by' variables prepare one GP per level
      # as.factor will keep unused levels needed for new data
      byval <- as.factor(get(byvar, data))
      byform <- str2formula(c(ifelse(cmc, "0", "1"), "byval"))
      con_mat <- model.matrix(byform)
      cons <- colnames(con_mat)
      out[[paste0("Kgp", pi)]] <- length(cons)
      Ngp <- Nsubgp <- vector("list", length(cons))
      for (j in seq_along(cons)) {
        # loop along contrasts of 'by'
        Cgp <- con_mat[, j]
        sfx <- paste0(pi, "_", j)
        tmp <- .data_gp(
          Xgp, k = k, gr = gr, sfx = sfx, Cgp = Cgp, c = c,
          scale = scale, internal = internal, basis = basis,
          ...
        )
        Ngp[[j]] <- attributes(tmp)[["Ngp"]]
        Nsubgp[[j]] <- attributes(tmp)[["Nsubgp"]]
        c(out) <- tmp
      }
      out[[paste0("Ngp", pi)]] <- unlist(Ngp)
      if (gr) {
        out[[paste0("Nsubgp", pi)]] <- unlist(Nsubgp)
      }
    } else {
      out[[paste0("Kgp", pi)]] <- 1L
      c(out) <- .data_gp(
        Xgp, k = k, gr = gr, sfx = pi, c = c,
        scale = scale, internal = internal, basis = basis,
        ...
      )
      if (bynum) {
        Cgp <- as.numeric(get(byvar, data))
        out[[paste0("Cgp", pi)]] <- as.array(Cgp)
      }
    }
  }
  if (length(basis)) {
    # original covariate values are required in new GP prediction
    Xgp_old <- basis[grepl("^Xgp", names(basis))]
    names(Xgp_old) <- paste0(names(Xgp_old), "_old")
    out[names(Xgp_old)] <- Xgp_old
  }
  out
}

# helper function to preparae GP related data
# @inheritParams data_gp
# @param Xgp matrix of covariate values
# @param k, gr, c see 'tidy_gpef'
# @param sfx suffix to put at the end of data names
# @param Cgp optional vector of values belonging to
#   a certain contrast of a factor 'by' variable
.data_gp <- function(Xgp, k, gr, sfx, Cgp = NULL, c = NULL,
                     scale = TRUE, internal = FALSE, basis = NULL) {
  out <- list()
  if (!is.null(Cgp)) {
    Cgp <- unname(Cgp)
    Igp <- which(Cgp != 0)
    Xgp <- Xgp[Igp, , drop = FALSE]
    out[[paste0("Igp", sfx)]] <- as.array(Igp)
    out[[paste0("Cgp", sfx)]] <- as.array(Cgp[Igp])
    attr(out, "Ngp") <- length(Igp)
  }
  if (gr) {
    groups <- factor(match_rows(Xgp, Xgp))
    ilevels <- levels(groups)
    Jgp <- match(groups, ilevels)
    Nsubgp <- length(ilevels)
    if (!is.null(Cgp)) {
      attr(out, "Nsubgp") <- Nsubgp
    } else {
      out[[paste0("Nsubgp", sfx)]]  <- Nsubgp
    }
    out[[paste0("Jgp", sfx)]] <- as.array(Jgp)
    not_dupl_Jgp <- !duplicated(Jgp)
    Xgp <- Xgp[not_dupl_Jgp, , drop = FALSE]
  }
  if (scale) {
    # scale predictor for easier specification of priors
    if (length(basis)) {
      # scale Xgp based on the original data
      dmax <- basis[[paste0("dmax", sfx)]]
    } else {
      dmax <- sqrt(max(diff_quad(Xgp)))
    }
    if (!isTRUE(dmax > 0)) {
      stop2("Could not scale GP covariates. Please set 'scale' to FALSE in 'gp'.")
    }
    if (internal) {
      # required for scaling of GPs with new data
      out[[paste0("dmax", sfx)]] <- dmax
    }
    Xgp <- Xgp / dmax
  }
  if (length(basis)) {
    # center Xgp based on the original data
    cmeans <- basis[[paste0("cmeans", sfx)]]
  } else {
    cmeans <- colMeans(Xgp)
  }
  if (internal) {
    # required for centering of approximate GPs with new data
    out[[paste0("cmeans", sfx)]] <- cmeans
    # required to compute inverse-gamma priors for length-scales
    out[[paste0("Xgp_prior", sfx)]] <- Xgp
  }
  if (!isNA(k)) {
    # basis function approach requires centered variables
    Xgp <- sweep(Xgp, 2, cmeans)
    D <- NCOL(Xgp)
    L <- choose_L(Xgp, c = c)
    Ks <- as.matrix(do_call(expand.grid, repl(seq_len(k), D)))
    XgpL <- matrix(nrow = NROW(Xgp), ncol = NROW(Ks))
    slambda <- matrix(nrow = NROW(Ks), ncol = D)
    for (m in seq_rows(Ks)) {
      XgpL[, m] <- eigen_fun_cov_exp_quad(Xgp, m = Ks[m, ], L = L)
      slambda[m, ] <- sqrt(eigen_val_cov_exp_quad(m = Ks[m, ], L = L))
    }
    out[[paste0("Xgp", sfx)]] <- XgpL
    out[[paste0("slambda", sfx)]] <- slambda
  } else {
    out[[paste0("Xgp", sfx)]] <- as.array(Xgp)
  }
  out
}

# data for autocorrelation variables
data_ac <- function(bterms, data, data2, basis = NULL, ...) {
  out <- list()
  N <- nrow(data)
  acef <- tidy_acef(bterms)
  if (has_ac_subset(bterms, dim = "time")) {
    gr <- get_ac_vars(acef, "gr", dim = "time")
    if (isTRUE(nzchar(gr))) {
      tgroup <- as.numeric(factor(data[[gr]]))
    } else {
      tgroup <- rep(1, N)
    }
  }
  if (has_ac_class(acef, "arma")) {
    # ARMA correlations
    acef_arma <- subset2(acef, class = "arma")
    out$Kar <- acef_arma$p
    out$Kma <- acef_arma$q
    if (!use_ac_cov_time(acef_arma)) {
      # data for the 'predictor' version of ARMA
      max_lag <- max(out$Kar, out$Kma)
      out$J_lag <- as.array(rep(0, N))
      for (n in seq_len(N)[-N]) {
        ind <- n:max(1, n + 1 - max_lag)
        # indexes errors to be used in the n+1th prediction
        out$J_lag[n] <- sum(tgroup[ind] %in% tgroup[n + 1])
      }
    }
  }
  if (use_ac_cov_time(acef)) {
    # data for the 'covariance' versions of time-series structures
    out$N_tg <- length(unique(tgroup))
    out$begin_tg <- as.array(ulapply(unique(tgroup), match, tgroup))
    out$nobs_tg <- as.array(with(out,
      c(if (N_tg > 1L) begin_tg[2:N_tg], N + 1) - begin_tg
    ))
    out$end_tg <- with(out, begin_tg + nobs_tg - 1)
    if (has_ac_class(acef, "unstr")) {
      time <- get_ac_vars(bterms, "time", dim = "time")
      time_data <- get(time, data)
      new_times <- extract_levels(time_data)
      if (length(basis)) {
        times <- basis$times
        # unstr estimates correlations only for given time point
        invalid_times <- setdiff(new_times, times)
        if (length(invalid_times)) {
          stop2("Cannot handle new time points in UNSTR models.")
        }
      } else {
        times <- new_times
      }
      out$n_unique_t <- length(times)
      out$n_unique_cortime <- out$n_unique_t * (out$n_unique_t - 1) / 2
      out$Iobs_tg <- matrix(0L, out$N_tg, max(out$nobs_tg))
      for (i in seq_len(out$N_tg)) {
        out$Iobs_tg[i, seq_len(out$nobs_tg[i])] <-
          time_data[out$begin_tg[i]:out$end_tg[i]]
      }
    }
  }
  if (has_ac_class(acef, "sar")) {
    acef_sar <- subset2(acef, class = "sar")
    M <- data2[[acef_sar$M]]
    rmd_rows <- attr(data, "na.action")
    if (!is.null(rmd_rows)) {
      class(rmd_rows) <- NULL
      M <- M[-rmd_rows, -rmd_rows, drop = FALSE]
    }
    if (!is_equal(dim(M), rep(N, 2))) {
      stop2("Dimensions of 'M' for SAR terms must be equal to ",
            "the number of observations.")
    }
    out$Msar <- as.matrix(M)
    out$eigenMsar <- eigen(M)$values
    # simplifies code of choose_N
    out$N_tg <- 1
  }
  if (has_ac_class(acef, "car")) {
    acef_car <- subset2(acef, class = "car")
    locations <- NULL
    if (length(basis)) {
      locations <- basis$locations
    }
    M <- data2[[acef_car$M]]
    if (acef_car$gr != "NA") {
      loc_data <- get(acef_car$gr, data)
      new_locations <- extract_levels(loc_data)
      if (is.null(locations)) {
        locations <- new_locations
      } else {
        invalid_locations <- setdiff(new_locations, locations)
        if (length(invalid_locations)) {
          stop2("Cannot handle new locations in CAR models.")
        }
      }
      Nloc <- length(locations)
      Jloc <- as.array(match(loc_data, locations))
      if (is.null(rownames(M))) {
        stop2("Row names are required for 'M' in CAR terms.")
      }
      found <- locations %in% rownames(M)
      if (any(!found)) {
        stop2("Row names of 'M' for CAR terms do not match ",
              "the names of the grouping levels.")
      }
      M <- M[locations, locations, drop = FALSE]
    } else {
      warning2(
        "Using CAR terms without a grouping factor is deprecated. ",
        "Please use argument 'gr' even if each observation ",
        "represents its own location."
      )
      Nloc <- N
      Jloc <- as.array(seq_len(Nloc))
      if (!is_equal(dim(M), rep(Nloc, 2))) {
        if (length(basis)) {
          stop2("Cannot handle new data in CAR terms ",
                "without a grouping factor.")
        } else {
          stop2("Dimensions of 'M' for CAR terms must be equal ",
                "to the number of observations.")
        }
      }
    }
    edges_rows <- (Matrix::tril(M)@i + 1)
    edges_cols <- sort(Matrix::triu(M)@i + 1) ## sort to make consistent with rows
    edges <- cbind("rows" = edges_rows, "cols" = edges_cols)
    c(out) <- nlist(
      Nloc, Jloc, Nedges = length(edges_rows),
      edges1 = as.array(edges_rows),
      edges2 = as.array(edges_cols)
    )
    if (acef_car$type %in% c("escar", "esicar")) {
      Nneigh <- Matrix::colSums(M)
      if (any(Nneigh == 0) && !length(basis)) {
        stop2(
          "For exact sparse CAR, all locations should have at ",
          "least one neighbor within the provided data set. ",
          "Consider using type = 'icar' instead."
        )
      }
      inv_sqrt_D <- diag(1 / sqrt(Nneigh))
      eigenMcar <- t(inv_sqrt_D) %*% M %*% inv_sqrt_D
      eigenMcar <- eigen(eigenMcar, TRUE, only.values = TRUE)$values
      c(out) <- nlist(Nneigh, eigenMcar)
    } else if (acef_car$type %in% "bym2") {
      c(out) <- list(car_scale = .car_scale(edges, Nloc))
    }
  }
  if (has_ac_class(acef, "fcor")) {
    acef_fcor <- subset2(acef, class = "fcor")
    M <- data2[[acef_fcor$M]]
    rmd_rows <- attr(data, "na.action")
    if (!is.null(rmd_rows)) {
      class(rmd_rows) <- NULL
      M <- M[-rmd_rows, -rmd_rows, drop = FALSE]
    }
    if (nrow(M) != N) {
      stop2("Dimensions of 'M' for FCOR terms must be equal ",
            "to the number of observations.")
    }
    out$Mfcor <- M
    # simplifies code of choose_N
    out$N_tg <- 1
  }
  if (length(out)) {
    resp <- usc(combine_prefix(bterms))
    out <- setNames(out, paste0(names(out), resp))
  }
  out
}

# prepare data of offsets for use in Stan
data_offset <- function(bterms, data) {
  out <- list()
  px <- check_prefix(bterms)
  if (is.formula(bterms$offset)) {
    p <- usc(combine_prefix(px))
    mf <- rm_attr(data, "terms")
    mf <- model.frame(bterms$offset, mf, na.action = na.pass)
    offset <- model.offset(mf)
    if (length(offset) == 1L) {
      offset <- rep(offset, nrow(data))
    }
    # use 'offsets' as 'offset' will be reserved in stanc3
    out[[paste0("offsets", p)]] <- as.array(offset)
  }
  out
}

# data for covariates in non-linear models
# @param x a btnl object
# @return a named list of data passed to Stan
data_cnl <- function(bterms, data) {
  stopifnot(is.btnl(bterms))
  out <- list()
  covars <- all.vars(bterms$covars)
  if (!length(covars)) {
    return(out)
  }
  p <- usc(combine_prefix(bterms))
  for (i in seq_along(covars)) {
    cvalues <- get(covars[i], data)
    if (is_like_factor(cvalues)) {
      # need to apply factor contrasts
      cform <- str2formula(covars[i])
      cvalues <- get_model_matrix(cform, data, cols2remove = "(Intercept)")
      if (NCOL(cvalues) > 1) {
        stop2("Factors with more than two levels are not allowed as covariates.")
      }
      cvalues <- cvalues[, 1]
    }
    out[[paste0("C", p, "_", i)]] <- as.array(cvalues)
  }
  out
}

# compute the spatial scaling factor of CAR models
# @param edges matrix with two columns defining the adjacency of the locations
# @param Nloc number of locations
# @return a scalar scaling factor
.car_scale <- function(edges, Nloc) {
  # amended from Imad Ali's code of CAR models in rstanarm
  stopifnot(is.matrix(edges), NCOL(edges) == 2)
  # Build the adjacency matrix
  adj_matrix <- Matrix::sparseMatrix(
    i = edges[, 1], j = edges[, 2], x = 1,
    symmetric = TRUE
  )
  # The ICAR precision matrix (which is singular)
  Q <- Matrix::Diagonal(Nloc, Matrix::rowSums(adj_matrix)) - adj_matrix
  # Add a small jitter to the diagonal for numerical stability
  Q_pert <- Q + Matrix::Diagonal(Nloc) *
    max(Matrix::diag(Q)) * sqrt(.Machine$double.eps)
  # Compute the diagonal elements of the covariance matrix subject to the
  # constraint that the entries of the ICAR sum to zero.
  .Q_inv <- function(Q) {
    Sigma <- Matrix::solve(Q)
    A <- matrix(1, 1, NROW(Sigma))
    W <- Sigma %*% t(A)
    Sigma <- Sigma - W %*% solve(A %*% W) %*% Matrix::t(W)
    return(Sigma)
  }
  Q_inv <- .Q_inv(Q_pert)
  # Compute the geometric mean of the variances (diagonal of Q_inv)
  exp(mean(log(Matrix::diag(Q_inv))))
}

# data for special priors such as horseshoe and lasso
data_prior <- function(bterms, data, prior, sdata = NULL) {
  out <- list()
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  special <- get_special_prior(prior, px)
  if (!is.null(special$horseshoe)) {
    # data for the horseshoe prior
    hs_names <- c("df", "df_global", "df_slab", "scale_global", "scale_slab")
    hs_data <- special$horseshoe[hs_names]
    if (!is.null(special$horseshoe$par_ratio)) {
      hs_data$scale_global <- special$horseshoe$par_ratio / sqrt(nrow(data))
    }
    names(hs_data) <- paste0("hs_", hs_names, p)
    out <- c(out, hs_data)
  }
  if (!is.null(special$R2D2)) {
    # data for the R2D2 prior
    R2D2_names <- c("mean_R2", "prec_R2", "cons_D2")
    R2D2_data <- special$R2D2[R2D2_names]
    # number of coefficients minus the intercept
    K <- sdata[[paste0("K", p)]] - ifelse(stan_center_X(bterms), 1, 0)
    if (length(R2D2_data$cons_D2) == 1L) {
      R2D2_data$cons_D2 <- rep(R2D2_data$cons_D2, K)
    }
    if (length(R2D2_data$cons_D2) != K) {
      stop2("Argument 'cons_D2' of the R2D2 prior must be of length 1 or ", K)
    }
    R2D2_data$cons_D2 <- as.array(R2D2_data$cons_D2)
    names(R2D2_data) <- paste0("R2D2_", R2D2_names, p)
    out <- c(out, R2D2_data)
  }
  if (!is.null(special$lasso)) {
    lasso_names <- c("df", "scale")
    lasso_data <- special$lasso[lasso_names]
    names(lasso_data) <- paste0("lasso_", lasso_names, p)
    out <- c(out, lasso_data)
  }
  out
}

# Construct design matrices for brms models
# @param formula a formula object
# @param data A data frame created with model.frame.
#   If another sort of object, model.frame is called first.
# @param cols2remove names of the columns to remove from
#   the model matrix; mainly used for intercepts
# @param rename rename column names via rename()?
# @param ... passed to stats::model.matrix
# @return
#   The design matrix for the given formula and data.
#   For details see ?stats::model.matrix
get_model_matrix <- function(formula, data = environment(formula),
                             cols2remove = NULL, rename = TRUE, ...) {
  stopifnot(is.atomic(cols2remove))
  terms <- validate_terms(formula)
  if (is.null(terms)) {
    return(NULL)
  }
  if (no_int(terms)) {
    cols2remove <- union(cols2remove, "(Intercept)")
  }
  X <- stats::model.matrix(terms, data, ...)
  cols2remove <- which(colnames(X) %in% cols2remove)
  if (length(cols2remove)) {
    X <- X[, -cols2remove, drop = FALSE]
  }
  if (rename) {
    colnames(X) <- rename(colnames(X), check_dup = TRUE)
  }
  X
}

# convenient wrapper around mgcv::PredictMat
PredictMat <- function(object, data, ...) {
  data <- rm_attr(data, "terms")
  out <- mgcv::PredictMat(object, data = data, ...)
  if (length(dim(out)) < 2L) {
    # fixes issue #494
    out <- matrix(out, nrow = 1)
  }
  out
}

# convenient wrapper around mgcv::smoothCon
smoothCon <- function(object, data, ...) {
  data <- rm_attr(data, "terms")
  vars <- setdiff(c(object$term, object$by), "NA")
  for (v in vars) {
    if (is_like_factor(data[[v]])) {
      # allow factor-like variables #562
      data[[v]] <- as.factor(data[[v]])
    } else if (inherits(data[[v]], "difftime")) {
      # mgcv cannot handle 'difftime' variables
      data[[v]] <- as.numeric(data[[v]])
    }
  }
  mgcv::smoothCon(object, data = data, ...)
}

# Aid prediction from smooths represented as 'type = 2'
# originally provided by Simon Wood
# @param sm output of mgcv::smoothCon
# @param data new data supplied for prediction
# @return A list of the same structure as returned by mgcv::smoothCon
s2rPred <- function(sm, data) {
  re <- mgcv::smooth2random(sm, names(data), type = 2)
  # prediction matrix for new data
  X <- PredictMat(sm, data)
  # transform to RE parameterization
  if (!is.null(re$trans.U)) {
    X <- X %*% re$trans.U
  }
  if (is.null(re$trans.D)) {
    # regression spline without penalization
    out <- list(Xf = X)
  } else {
    X <- t(t(X) * re$trans.D)
    # re-order columns according to random effect re-ordering
    X[, re$rind] <- X[, re$pen.ind != 0]
    # re-order penalization index in same way
    pen.ind <- re$pen.ind
    pen.ind[re$rind] <- pen.ind[pen.ind > 0]
    # start returning the object
    Xf <- X[, which(re$pen.ind == 0), drop = FALSE]
    out <- list(rand = list(), Xf = Xf)
    for (i in seq_along(re$rand)) {
      # loop over random effect matrices
      out$rand[[i]] <- X[, which(pen.ind == i), drop = FALSE]
      attr(out$rand[[i]], "s.label") <- attr(re$rand[[i]], "s.label")
    }
    names(out$rand) <- names(re$rand)
  }
  out
}
