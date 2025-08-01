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
data_predictor.mvbrmsterms <- function(x, data, sdata = NULL, ...) {
  out <- list(N = nrow(data))
  for (r in names(x$terms)) {
    c(out) <- data_predictor(x$terms[[r]], data = data, sdata = sdata, ...)
  }
  out
}

#' @export
data_predictor.brmsterms <- function(x, data, data2, prior, sdata = NULL, ...) {
  out <- list()
  data <- subset_data(data, x)
  resp <- usc(combine_prefix(x))
  args_eff <- nlist(data, data2, prior, sdata, ...)
  for (dp in names(x$dpars)) {
    args_eff_spec <- list(x = x$dpars[[dp]])
    c(out) <- do_call(data_predictor, c(args_eff_spec, args_eff))
  }
  for (dp in names(x$fdpars)) {
    if (is.numeric(x$fdpars[[dp]]$value)) {
      out[[paste0(dp, resp)]] <- x$fdpars[[dp]]$value
    }
  }
  for (nlp in names(x$nlpars)) {
    args_eff_spec <- list(x = x$nlpars[[nlp]])
    c(out) <- do_call(data_predictor, c(args_eff_spec, args_eff))
  }
  c(out) <- data_gr_local(x, data = data)
  c(out) <- data_mixture(x, data2 = data2, prior = prior)
  out
}

# prepare data for all types of effects for use in Stan
# @param data the data passed by the user
# @param prior an object of class brmsprior
# @param ... currently ignored
# @return a named list of data to be passed to Stan
#' @export
data_predictor.btl <- function(x, data, data2 = list(), prior = brmsprior(),
                               sdata = NULL, ...) {
  out <- c(
    data_fe(x, data),
    data_sp(x, data, data2 = data2, prior = prior),
    data_re(x, data),
    data_cs(x, data),
    data_sm(x, data),
    data_gp(x, data),
    data_ac(x, data, data2 = data2),
    data_offset(x, data),
    data_bhaz(x, data, data2 = data2, prior = prior)
  )
  c(out) <- data_special_prior(x, data, prior = prior, sdata = c(sdata, out))
  out
}

# prepare data for non-linear parameters for use in Stan
#' @export
data_predictor.btnl <- function(x, data, data2 = list(), prior = brmsprior(),
                                ...) {
  out <- list()
  c(out) <- data_cnl(x, data)
  c(out) <- data_ac(x, data, data2 = data2)
  c(out) <- data_bhaz(x, data, data2 = data2, prior = prior)
  out
}

# prepare data of fixed effects
data_fe <- function(bframe, data) {
  stopifnot(is.btl(bframe))
  if (!is.null(bframe$sdata$fe)) {
    # standata was already precomputed
    return(bframe$sdata$fe)
  }
  out <- list()
  p <- usc(combine_prefix(bframe))
  # the intercept is removed inside the Stan code for non-ordinal models
  is_ord <- is_ordinal(bframe)
  cols2remove <- if (is_ord) "(Intercept)"
  X <- get_model_matrix(rhs(bframe$fe), data, cols2remove = cols2remove)
  avoid_dpars(colnames(X), bframe)
  out[[paste0("K", p)]] <- ncol(X)
  if (stan_center_X(bframe)) {
    # relevant if the intercept is treated separately to enable centering
    out[[paste0("Kc", p)]] <- ncol(X) - ifelse(is_ord, 0, 1)
  }
  out[[paste0("X", p)]] <- X
  out
}

# data preparation for splines
data_sm <- function(bframe, data) {
  stopifnot(is.btl(bframe))
  if (!is.null(bframe$sdata$sm)) {
    # standata was already precomputed
    return(bframe$sdata$sm)
  }
  out <- list()
  smterms <- all_terms(bframe[["sm"]])
  if (!length(smterms)) {
    return(out)
  }
  p <- usc(combine_prefix(bframe))
  # basis contains information on the smooths from the original data
  basis <- bframe$basis$sm
  new <- length(basis) > 0L
  knots <- get_knots(data)
  diagonal.penalty <- !require_old_default("2.8.7")
  bylevels <- named_list(smterms)
  ns <- 0
  lXs <- list()
  for (i in seq_along(smterms)) {
    if (new) {
      sm <- basis[[i]]$sm
    } else {
      sm <- smoothCon(
        eval2(smterms[i]), data = data,
        knots = knots, absorb.cons = TRUE,
        diagonal.penalty = diagonal.penalty
      )
    }
    # may contain multiple terms when 'by' is a factor
    for (j in seq_along(sm)) {
      ns <- ns + 1
      if (length(sm[[j]]$by.level)) {
        bylevels[[i]][j] <- sm[[j]]$by.level
      }
      if (new) {
        # prepare smooths for use with new data
        # mgcv smooths are based on machine-specific SVD (#1465)
        re <- s2rPred(sm[[j]], re = basis[[i]]$re[[j]], data = data)
      } else {
        re <- mgcv::smooth2random(sm[[j]], names(data), type = 2)
      }
      lXs[[ns]] <- re$Xf
      if (NCOL(lXs[[ns]])) {
        colnames(lXs[[ns]]) <- paste0(sm[[j]]$label, "_", seq_cols(lXs[[ns]]))
      }
      Zs <- re$rand
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
  avoid_dpars(colnames(Xs), bframe)
  smcols <- lapply(lXs, function(x) which(colnames(Xs) %in% colnames(x)))
  Xs <- structure(Xs, smcols = smcols, bylevels = bylevels)
  colnames(Xs) <- rename(colnames(Xs))
  out[[paste0("Ks", p)]] <- ncol(Xs)
  out[[paste0("Xs", p)]] <- Xs
  out
}

# prepare data for group-level effects for use in Stan
data_re <- function(bframe, data) {
  stopifnot(is.bframel(bframe))
  out <- list()
  px <- check_prefix(bframe)
  reframe <- subset2(bframe$frame$re, type = "sp", fun = "%notin%")
  if (!has_rows(reframe)) {
    return(out)
  }
  gn <- unique(reframe$gn)
  for (i in seq_along(gn)) {
    r <- subset2(reframe, gn = gn[i])
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
data_gr_local <- function(bframe, data) {
  stopifnot(is.brmsframe(bframe))
  out <- list()
  reframe <- subset2(bframe$frame$re, resp = bframe$resp)
  resp <- usc(bframe$resp)
  for (id in unique(reframe$id)) {
    id_reframe <- subset2(reframe, id = id)
    idresp <- paste0(id, resp)
    nranef <- nrow(id_reframe)
    group <- id_reframe$group[1]
    levels <- get_levels(reframe)[[group]]
    if (id_reframe$gtype[1] == "mm") {
      # multi-membership grouping term
      gs <- id_reframe$gcall[[1]]$groups
      ngs <- length(gs)
      weights <- id_reframe$gcall[[1]]$weights
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
        # all members get equal membership weights by default
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
      g <- id_reframe$gcall[[1]]$groups
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
    # prepare data for group prior weights if specified
    if (nzchar(id_reframe$gcall[[1]]$pw)) {
      if (id_reframe$gtype[1] == "mm") {
        J <- unlist(out[paste0("J_", idresp, "_", seq_along(gs))])
      }
      # extract and validate prior weights
      group_prior_weights <- str2formula(id_reframe$gcall[[1]]$pw)
      group_prior_weights <- as.vector(eval_rhs(group_prior_weights, data))
      if (!is.numeric(group_prior_weights)) {
        stop2("Prior weights of grouping factors must be numeric.")
      }
      if (any(group_prior_weights < 0)) {
        warning2("Negative prior weights detected. Make sure this is intentional.")
      }
      # check that group-level weights do not vary within a group
      group_weights_consistent <- tapply(
        X = group_prior_weights, INDEX = J,
        FUN = function(x) length(unique(x)) == 1
      )
      if (!all(group_weights_consistent)) {
        stop2("Prior weights cannot vary within a group.")
      }
      # deduplicate weights vector (so length matches number of groups)
      # and order the weights vector to match groups' assigned indices
      distinct_J_indices <- !duplicated(J)
      group_prior_weights <- group_prior_weights[distinct_J_indices]
      group_prior_weights <- group_prior_weights[order(J[distinct_J_indices])]
      out[[paste0("PW_", id)]] <- as.array(group_prior_weights)
    }
  }
  out
}

# prepare global data for each group-level-ID
data_gr_global <- function(bframe, data2) {
  stopifnot(is.anybrmsframe(bframe))
  out <- list()
  reframe <- bframe$frame$re
  for (id in unique(reframe$id)) {
    tmp <- list()
    id_reframe <- subset2(reframe, id = id)
    nranef <- nrow(id_reframe)
    group <- id_reframe$group[1]
    levels <- attr(reframe, "levels")[[group]]
    tmp$N <- length(levels)
    tmp$M <- nranef
    tmp$NC <- as.integer(nranef * (nranef - 1) / 2)
    # prepare number of levels of an optional 'by' variable
    if (nzchar(id_reframe$by[1])) {
      stopifnot(!nzchar(id_reframe$type[1]))
      bylevels <- id_reframe$bylevels[[1]]
      Jby <- match(attr(levels, "by"), bylevels)
      tmp$Nby <- length(bylevels)
      tmp$Jby <- as.array(Jby)
    }
    # prepare within-group covariance matrices
    cov <- id_reframe$cov[1]
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
data_sp <- function(bframe, data, data2, prior) {
  stopifnot(is.bframel(bframe))
  if (!is.null(bframe$sdata$sp)) {
    # standata was already precomputed
    return(bframe$sdata$sp)
  }
  out <- list()
  spframe <- bframe$frame$sp
  if (!has_rows(spframe)) {
    return(out)
  }
  basis <- bframe$basis$sp
  px <- check_prefix(bframe)
  p <- usc(combine_prefix(px))
  # prepare general data
  out[[paste0("Ksp", p)]] <- nrow(spframe)
  Csp <- sp_model_matrix(bframe$sp, data)
  avoid_dpars(colnames(Csp), bframe)
  Csp <- Csp[, spframe$Ic > 0, drop = FALSE]
  Csp <- lapply(seq_cols(Csp), function(i) as.array(Csp[, i]))
  if (length(Csp)) {
    Csp_names <- paste0("Csp", p, "_", seq_along(Csp))
    out <- c(out, setNames(Csp, Csp_names))
  }
  if (any(lengths(spframe$Imo) > 0)) {
    # prepare data specific to monotonic effects
    out[[paste0("Imo", p)]] <- max(unlist(spframe$Imo))
    Xmo <- lapply(unlist(spframe$calls_mo), get_mo_values, data = data)
    Xmo_names <- paste0("Xmo", p, "_", seq_along(Xmo))
    c(out) <- setNames(Xmo, Xmo_names)
    if (!is.null(basis$Jmo)) {
      # take information from original data
      Jmo <- basis$Jmo
    } else {
      Jmo <- as.array(ulapply(Xmo, attr, "max"))
    }
    out[[paste0("Jmo", p)]] <- Jmo
    # prepare prior concentration of simplex parameters
    simo_coef <- get_simo_labels(spframe, use_id = TRUE)
    ids <- unlist(spframe$ids_mo)
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
  uni_mi <- attr(spframe, "uni_mi")
  index <- bframe$frame$index
  for (j in seq_rows(uni_mi)) {
    has_idx <- !is.na(uni_mi$idx[j])
    if (has_idx) {
      idxl <- get(uni_mi$idx[j], data)
      # a fully NA index variable will be treated as not having an index variable
      # this is useful in post-processing to turn off indexing entirely if desired
      has_idx <- !all(is.na(idxl))
    }
    has_subset <- isTRUE(attr(index[[uni_mi$var[j]]], "subset"))
    if (has_idx) {
      if (is.null(index[[uni_mi$var[j]]])) {
        # the 'idx' argument needs to be mapped against 'index' addition terms
        stop2("Response '", uni_mi$var[j], "' needs to have an 'mi' addition ",
              "term with an 'idx' variable specified to compare with 'idx' ",
              "variables in 'mi' predictor terms. See ?mi for examples.")
      }
      idxl <- match(idxl, index[[uni_mi$var[j]]])
      if (anyNA(idxl)) {
        stop2("Could not match all indices in response '", uni_mi$var[j], "'.")
      }
      idxl_name <- paste0("idxl", p, "_", uni_mi$var[j], "_", uni_mi$idx2[j])
      out[[idxl_name]] <- as.array(idxl)
    } else if (has_subset) {
      # cross-formula referencing is required for subsetted variables
      stop2("'mi' predictor terms of subsetted variables require ",
            "the 'idx' argument to be specified.")
    }
  }
  out
}

# prepare data for category specific effects
data_cs <- function(bframe, data) {
  stopifnot(is.btl(bframe))
  if (!is.null(bframe$sdata$cs)) {
    # standata was already precomputed
    return(bframe$sdata$cs)
  }
  out <- list()
  if (length(all_terms(bframe[["cs"]]))) {
    p <- usc(combine_prefix(bframe))
    Xcs <- get_model_matrix(bframe$cs, data)
    avoid_dpars(colnames(Xcs), bframe)
    out <- c(out, list(Kcs = ncol(Xcs), Xcs = Xcs))
    out <- setNames(out, paste0(names(out), p))
  }
  out
}

# prepare global data for noise free variables
data_Xme <- function(bframe, data) {
  stopifnot(is.anybrmsframe(bframe))
  meframe <- bframe$frame$me
  stopifnot(is.meframe(meframe))
  out <- list()
  groups <- unique(meframe$grname)
  for (i in seq_along(groups)) {
    g <- groups[i]
    K <- which(meframe$grname %in% g)
    Mme <- length(K)
    out[[paste0("Mme_", i)]] <- Mme
    out[[paste0("NCme_", i)]] <- Mme * (Mme - 1) / 2
    if (nzchar(g)) {
      levels <- get_levels(meframe)[[g]]
      gr <- get_me_group(meframe$term[K[1]], data)
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
      Xn <- get_me_values(meframe$term[k], data)
      noise <- get_me_noise(meframe$term[k], data)
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
data_gp <- function(bframe, data, internal = FALSE, ...) {
  stopifnot(is.bframel(bframe))
  if (!is.null(bframe$sdata$gp)) {
    # standata was already precomputed
    return(bframe$sdata$gp)
  }
  out <- list()
  internal <- as_one_logical(internal)
  px <- check_prefix(bframe)
  p <- usc(combine_prefix(px))
  basis <- bframe$basis$gp
  gpframe <- bframe$frame$gp
  for (i in seq_rows(gpframe)) {
    pi <- paste0(p, "_", i)
    Xgp <- lapply(gpframe$covars[[i]], eval2, data)
    D <- length(Xgp)
    out[[paste0("Dgp", pi)]] <- D
    invalid <- ulapply(Xgp, function(x)
      !is.numeric(x) || isTRUE(length(dim(x)) > 1L)
    )
    if (any(invalid)) {
      stop2("Predictors of Gaussian processes should be numeric vectors.")
    }
    Xgp <- do_call(cbind, Xgp)
    cmc <- gpframe$cmc[i]
    scale <- gpframe$scale[i]
    gr <- gpframe$gr[i]
    k <- gpframe$k[i]
    c <- gpframe$c[[i]]
    if (!isNA(k)) {
      out[[paste0("NBgp", pi)]] <- k ^ D
      Ks <- as.matrix(do_call(expand.grid, repl(seq_len(k), D)))
    }
    byvar <- gpframe$byvars[[i]]
    byfac <- length(gpframe$cons[[i]]) > 0L
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
# @param k, gr, c see 'frame_gp'
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

    if (length(basis)) {
      L <- basis[[paste0("Lgp", sfx)]]
    } else {
      # compute boundary factor L
      L <- choose_L(Xgp, c = c)
    }

    if (internal) {
      # required to compute eigenfunctions of approximate GPs with new data
      out[[paste0("Lgp", sfx)]] <- L
    }

    Ks <- as.matrix(do_call(expand.grid, repl(seq_len(k), D)))
    XgpL <- matrix(nrow = NROW(Xgp), ncol = NROW(Ks))
    slambda <- matrix(nrow = NROW(Ks), ncol = D)
    for (m in seq_rows(Ks)) {
      XgpL[, m] <- eigen_fun_laplacian(Xgp, m = Ks[m, ], L = L)
      slambda[m, ] <- sqrt(eigen_val_laplacian(m = Ks[m, ], L = L))
    }
    out[[paste0("Xgp", sfx)]] <- XgpL
    out[[paste0("slambda", sfx)]] <- slambda
  } else {
    out[[paste0("Xgp", sfx)]] <- as.array(Xgp)
  }
  out
}

# data for autocorrelation variables
data_ac <- function(bframe, data, data2, ...) {
  if (!is.null(bframe$sdata$ac)) {
    # standata was already precomputed
    return(bframe$sdata$ac)
  }
  out <- list()
  N <- nrow(data)
  basis <- bframe$basis$ac
  acframe <- bframe$frame$ac
  stopifnot(is.acframe(acframe))
  if (has_ac_subset(bframe, dim = "time")) {
    gr <- get_ac_vars(acframe, "gr", dim = "time")
    if (isTRUE(nzchar(gr))) {
      tgroup <- as.numeric(factor(data[[gr]]))
    } else {
      tgroup <- rep(1, N)
    }
  }
  if (has_ac_class(acframe, "arma")) {
    # ARMA correlations
    acframe_arma <- subset2(acframe, class = "arma")
    out$Kar <- acframe_arma$p
    out$Kma <- acframe_arma$q
    if (!use_ac_cov_time(acframe_arma)) {
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
  if (use_ac_cov_time(acframe)) {
    # data for the 'covariance' versions of time-series structures
    # TODO: change begin[i]:end[i] notation to slice[i]:(slice[i+1] - 1)
    #   see comment on PR #1435
    out$N_tg <- length(unique(tgroup))
    out$begin_tg <- as.array(ulapply(unique(tgroup), match, tgroup))
    out$nobs_tg <- as.array(with(out,
      c(if (N_tg > 1L) begin_tg[2:N_tg], N + 1) - begin_tg
    ))
    out$end_tg <- with(out, begin_tg + nobs_tg - 1)
    if (has_ac_class(acframe, "unstr")) {
      time <- get_ac_vars(bframe, "time", dim = "time")
      time_data <- get(time, data)
      new_times <- extract_levels(time_data)
      if (length(basis)) {
        times <- basis$times
        # unstr estimates correlations only for given time points
        invalid_times <- setdiff(new_times, times)
        if (length(invalid_times)) {
          stop2("Cannot handle new time points in UNSTR models.")
        }
      } else {
        times <- new_times
      }
      out$n_unique_t <- length(times)
      out$n_unique_cortime <- out$n_unique_t * (out$n_unique_t - 1) / 2
      Jtime <- match(time_data, times)
      out$Jtime_tg <- matrix(0L, out$N_tg, max(out$nobs_tg))
      for (i in seq_len(out$N_tg)) {
        out$Jtime_tg[i, seq_len(out$nobs_tg[i])] <-
          Jtime[out$begin_tg[i]:out$end_tg[i]]
      }
    }
  }
  if (has_ac_class(acframe, "sar")) {
    acframe_sar <- subset2(acframe, class = "sar")
    M <- data2[[acframe_sar$M]]
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
  if (has_ac_class(acframe, "car")) {
    acframe_car <- subset2(acframe, class = "car")
    locations <- NULL
    if (length(basis)) {
      locations <- basis$locations
    }
    M <- data2[[acframe_car$M]]
    if (acframe_car$gr != "NA") {
      loc_data <- get(acframe_car$gr, data)
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
    if (acframe_car$type %in% c("escar", "esicar")) {
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
    } else if (acframe_car$type %in% "bym2") {
      c(out) <- list(car_scale = .car_scale(edges, Nloc))
    }
  }
  if (has_ac_class(acframe, "fcor")) {
    acframe_fcor <- subset2(acframe, class = "fcor")
    M <- data2[[acframe_fcor$M]]
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
    resp <- usc(combine_prefix(bframe))
    out <- setNames(out, paste0(names(out), resp))
  }
  out
}

# prepare data of offsets for use in Stan
data_offset <- function(bframe, data) {
  stopifnot(is.btl(bframe))
  if (!is.null(bframe$sdata$offset)) {
    # standata was already precomputed
    return(bframe$sdata$offset)
  }
  out <- list()
  px <- check_prefix(bframe)
  if (is.formula(bframe$offset)) {
    p <- usc(combine_prefix(px))
    mf <- rm_attr(data, "terms")
    mf <- model.frame(bframe$offset, mf, na.action = na.pass)
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
data_cnl <- function(bframe, data) {
  stopifnot(is.btnl(bframe))
  if (!is.null(bframe$sdata$cnl)) {
    # standata was already precomputed
    return(bframe$sdata$cnl)
  }
  out <- list()
  covars <- all.vars(bframe$covars)
  if (!length(covars)) {
    return(out)
  }
  p <- usc(combine_prefix(bframe))
  for (i in seq_along(covars)) {
    cvalues <- get(covars[i], data)
    if (is_like_factor(cvalues)) {
      # need to apply factor contrasts
      cform <- str2formula(covars[i])
      cvalues <- get_model_matrix(cform, data, cols2remove = "(Intercept)")
      if (NCOL(cvalues) == 1L) {
        dim(cvalues) <- NULL
      }
    }
    if (isTRUE(dim(cvalues) > 2L)) {
      stop2("Non-linear covariates should be vectors or matrices.")
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

# data for special priors such as horseshoe and R2D2
data_special_prior <- function(bframe, data, prior, sdata = NULL) {
  out <- list()
  px <- check_prefix(bframe)
  p <- usc(combine_prefix(px))
  if (!has_special_prior(prior, px)) {
    return(out)
  }

  # number of coefficients affected by the shrinkage prior
  # fully compute this here to avoid having to pass the prior around
  # to all the individual data preparation functions
  # the order of adding things to Kscales doesn't matter but for consistency
  # it is still the same as the order in the Stan code
  Kscales <- 0
  if (has_special_prior(prior, px, class = "b")) {
    Kscales <- Kscales +
      sdata[[paste0("Kc", p)]] %||% sdata[[paste0("K", p)]] %||% 0 +
      sdata[[paste0("Ksp", p)]] %||% 0 +
      sdata[[paste0("Ks", p)]] %||% 0
  }
  if (has_special_prior(prior, px, class = "sds")) {
    take <- grepl(paste0("^nb", p, "_"), names(sdata))
    Kscales <- Kscales + sum(unlist(sdata[take]))
  }
  if (has_special_prior(prior, px, class = "sdgp")) {
    take <- grepl(paste0("^Kgp", p, "_"), names(sdata))
    Kscales <- Kscales + sum(unlist(sdata[take]))
  }
  if (has_special_prior(prior, px, class = "ar")) {
    Kscales <- Kscales + sdata[[paste0("Kar", p)]]
  }
  if (has_special_prior(prior, px, class = "ma")) {
    Kscales <- Kscales + sdata[[paste0("Kma", p)]]
  }
  if (has_special_prior(prior, px, class = "sderr")) {
    Kscales <- Kscales + 1
  }
  if (has_special_prior(prior, px, class = "sdcar")) {
    Kscales <- Kscales + 1
  }
  if (has_special_prior(prior, px, class = "sd")) {
    ids <- unique(bframe$frame$re$id)
    Kscales <- Kscales + sum(unlist(sdata[paste0("M_", ids)]))
  }
  out[[paste0("Kscales", p)]] <- Kscales

  special <- get_special_prior(prior, px, main = TRUE)
  if (special$name == "horseshoe") {
    # data for the horseshoe prior
    hs_names <- c("df", "df_global", "df_slab", "scale_global", "scale_slab")
    hs_data <- special[hs_names]
    if (!is.null(special$par_ratio)) {
      hs_data$scale_global <- special$par_ratio / sqrt(nrow(data))
    }
    names(hs_data) <- paste0("hs_", hs_names, p)
    c(out) <- hs_data
  } else if (special$name == "R2D2") {
    # data for the R2D2 prior
    R2D2_names <- c("mean_R2", "prec_R2", "cons_D2")
    R2D2_data <- special[R2D2_names]
    if (length(R2D2_data$cons_D2) == 1L) {
      R2D2_data$cons_D2 <- rep(R2D2_data$cons_D2, Kscales)
    }
    if (length(R2D2_data$cons_D2) != Kscales) {
      stop2("Argument 'cons_D2' of the R2D2 prior must be of length 1 or ", Kscales)
    }
    R2D2_data$cons_D2 <- as.array(R2D2_data$cons_D2)
    names(R2D2_data) <- paste0("R2D2_", R2D2_names, p)
    c(out) <- R2D2_data
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
  stopifnot(is_atomic_or_null(cols2remove))
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
  data <- sm_prepare_data(object, data)
  out <- mgcv::PredictMat(object, data = data, ...)
  if (length(dim(out)) < 2L) {
    # fixes issue #494
    out <- matrix(out, nrow = 1)
  }
  out
}

# convenient wrapper around mgcv::smoothCon
smoothCon <- function(object, data, ...) {
  data <- sm_prepare_data(object, data)
  mgcv::smoothCon(object, data = data, ...)
}

# mgcv doesn't handle a lot of special data types well
# need to prepare these variables manually beforehand
sm_prepare_data <- function(object, data) {
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
  data
}

# Aid prediction from smooths represented as 'type = 2'
# code obtained from the doc of ?mgcv::smooth2random
# @param sm output of mgcv::smoothCon
# @param re output of mgcv::smooth2random
# @param data new data supplied for prediction
# @return A list of the same structure as returned by mgcv::smooth2random
s2rPred <- function(sm, re, data) {
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
