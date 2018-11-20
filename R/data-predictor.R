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
data_predictor.mvbrmsterms <- function(x, old_sdata = NULL, ...) {
  out <- list()
  for (i in seq_along(x$terms)) {
    od <- old_sdata[[x$responses[i]]]
    out <- c(out, data_predictor(x$terms[[i]], old_sdata = od, ...))
  }
  out
}

#' @export
data_predictor.brmsterms <- function(x, data, prior, ranef, meef, 
                                     cov_ranef = NULL, knots = NULL, 
                                     not4stan = FALSE, old_sdata = NULL, ...) {
  out <- list()
  args_eff <- nlist(data, ranef, prior, knots, not4stan)
  for (dp in names(x$dpars)) {
    args_eff_spec <- list(x = x$dpars[[dp]], old_sdata = old_sdata[[dp]])
    c(out) <- run(data_predictor, c(args_eff_spec, args_eff))
  }
  for (dp in names(x$fdpars)) {
    resp <- usc(combine_prefix(x))
    out[[paste0(dp, resp)]] <- x$fdpars[[dp]]$value
  }
  for (nlp in names(x$nlpars)) {
    args_eff_spec <- list(x = x$nlpars[[nlp]], old_sdata = old_sdata[[nlp]])
    c(out) <- run(data_predictor, c(args_eff_spec, args_eff))
  }
  c(out,
    data_gr(ranef, data, cov_ranef = cov_ranef),
    data_Xme(meef, data),
    data_mixture(x, prior = prior)
  )
}

#' @export
data_predictor.btl <- function(x, data, ranef = empty_ranef(), 
                               prior = brmsprior(), knots = NULL, 
                               not4stan = FALSE, old_sdata = NULL, ...) {
  # prepare data for all types of effects for use in Stan
  # Args:
  #   data: the data passed by the user
  #   family: the model family
  #   prior: an object of class brmsprior
  #   autocor: object of class 'cor_brms'
  #   cov_ranef: name list of user-defined covariance matrices
  #   knots: optional knot values for smoothing terms
  #   nlpar: optional character string naming a non-linear parameter
  #   not4stan: is the data for use in S3 methods only?
  #   old_sdata: see 'extract_old_standata'
  # Returns:
  #   A named list of data to be passed to Stan
  c(data_fe(x, data, not4stan = not4stan),
    data_sp(x, data, prior = prior, Jmo = old_sdata$Jmo),
    data_re(x, data, ranef = ranef),
    data_cs(x, data),
    data_sm(x, data, knots = knots, smooths = old_sdata$smooths),
    data_gp(x, data, gps = old_sdata$gps),
    data_offset(x, data),
    data_prior(x, data, prior = prior)
  )
}

#' @export 
data_predictor.btnl <- function(x, data, not4stan = FALSE, ...) {
  # prepare data for non-linear parameters for use in Stan
  # matrix of covariates appearing in the non-linear formula
  out <- list()
  C <- get_model_matrix(x$covars, data = data)
  if (length(all.vars(x$covars)) != ncol(C)) {
    stop2("Factors with more than two levels are not allowed as covariates.")
  }
  # fixes issue #127 occuring for factorial covariates
  colnames(C) <- all.vars(x$covars)
  p <- usc(combine_prefix(x))
  if (not4stan) {
    out[[paste0("C", p)]] <- C
  } else {
    # use vectors as indexing matrices in Stan is slow
    if (ncol(C)) {
      Cnames <- paste0("C", p, "_", seq_cols(C))
      out <- c(out, setNames(as.list(as.data.frame(C)), Cnames))
    }
  }
  out
}

data_fe <- function(bterms, data, not4stan = FALSE) {
  # prepare data of fixed effects and smooth terms for use in Stan
  # handle smooth terms here as they also affect the FE design matrix 
  # Args: see data_predictor
  out <- list()
  p <- usc(combine_prefix(bterms))
  is_ordinal <- is_ordinal(bterms$family)
  is_bsts <- is.cor_bsts(bterms$autocor)
  # the intercept is removed inside the Stan code for ordinal models
  cols2remove <- if (is_ordinal && not4stan || is_bsts) "(Intercept)"
  X <- get_model_matrix(rhs(bterms$fe), data, cols2remove = cols2remove)
  avoid_dpars(colnames(X), bterms = bterms)
  out[[paste0("K", p)]] <- ncol(X)
  out[[paste0("X", p)]] <- X
  out
}

data_sm <- function(bterms, data, knots = NULL, smooths = NULL) {
  # data preparation for splines
  out <- list()
  smterms <- all_terms(bterms[["sm"]])
  if (!length(smterms)) {
    return(out)
  }
  p <- usc(combine_prefix(bterms))
  new <- length(smooths) > 0L
  if (!new) {
    smooths <- named_list(smterms)
    for (i in seq_along(smterms)) {
      smooths[[i]] <- mgcv::smoothCon(
        eval2(smterms[i]), data = data, 
        knots = knots, absorb.cons = TRUE
      )
    }
  }
  bylevels <- named_list(smterms)
  ns <- 0
  lXs <- list()
  for (i in seq_along(smooths)) {
    # may contain multiple terms when 'by' is a factor
    for (j in seq_along(smooths[[i]])) {
      ns <- ns + 1
      sm <- smooths[[i]][[j]]
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
      Zs <- setNames(Zs, paste0("Zs", p, "_", ns, "_", seq_along(Zs)))
      tmp <- list(length(Zs), as.array(ulapply(Zs, ncol)))
      tmp <- setNames(tmp, paste0(c("nb", "knots"), p, "_", ns))
      c(out) <- c(tmp, Zs)
    }
  }
  Xs <- run(cbind, lXs)
  avoid_dpars(colnames(Xs), bterms = bterms)
  smcols <- lapply(lXs, function(x) which(colnames(Xs) %in% colnames(x)))
  Xs <- structure(Xs, smcols = smcols, bylevels = bylevels)
  colnames(Xs) <- rename(colnames(Xs))
  out[[paste0("Ks", p)]] <- ncol(Xs)
  out[[paste0("Xs", p)]] <- Xs
  out
}

data_re <- function(bterms, data, ranef) {
  # prepare data for group-level effects for use in Stan
  # Args: see data_predictor
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
        Z <- run(cbind, unlist(Z_temp, recursive = FALSE))
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

data_gr <- function(ranef, data, cov_ranef = NULL) {
  # compute data specific for each group-level-ID
  # Args:
  #   ranef: data.frame returned by tidy_ranef
  #   cov_ranef: name list of user-defined covariance matrices
  out <- list()
  ids <- unique(ranef$id)
  for (id in ids) {
    id_ranef <- subset2(ranef, id = id)
    nranef <- nrow(id_ranef)
    group <- id_ranef$group[1]
    levels <- attr(ranef, "levels")[[group]]
    if (id_ranef$gtype[1] == "mm") {
      # multi-membership grouping term
      stopifnot(!nzchar(id_ranef$by[1]))
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
        J <- as.array(match(get(gs[i], data), levels))
        out[[paste0("J_", id, "_", i)]] <- J
        out[[paste0("W_", id, "_", i)]] <- as.array(weights[, i])
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
      out[[paste0("J_", id)]] <- as.array(J)
      if (nzchar(id_ranef$by[1])) {
        stopifnot(!nzchar(id_ranef$type[1]))
        bylevels <- id_ranef$bylevels[[1]]
        Jby <- match(attr(levels, "by"), bylevels)
        out[[paste0("Nby_", id)]] <- length(bylevels)
        out[[paste0("Jby_", id)]] <- as.array(Jby)
      }
    }
    temp <- list(length(levels), nranef, nranef * (nranef - 1) / 2)
    out <- c(out, setNames(temp, paste0(c("N_", "M_", "NC_"), id)))
    # prepare customized covariance matrices
    if (group %in% names(cov_ranef)) {
      cov_mat <- as.matrix(cov_ranef[[group]])
      if (!isSymmetric(unname(cov_mat))) {
        stop2("Covariance matrix of grouping factor '", 
              group, "' is not symmetric.")
      }
      found_levels <- rownames(cov_mat)
      if (is.null(found_levels)) {
        stop2("Row names are required for covariance matrix of '", group, "'.")
      }
      colnames(cov_mat) <- found_levels
      found <- levels %in% found_levels
      if (any(!found)) {
        stop2("Row names of covariance matrix of '", group, 
              "' do not match names of the grouping levels.")
      }
      cov_mat <- cov_mat[levels, levels, drop = FALSE]
      evs <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values
      if (min(evs) <= 0) {
        stop2("Covariance matrix of grouping factor '", 
              group, "' is not positive definite.")
      }
      out <- c(out, setNames(list(t(chol(cov_mat))), paste0("Lcov_", id)))
    }
  }
  out
}

data_sp <- function(bterms, data, prior = brmsprior(), Jmo = NULL) {
  # prepare data for special effects for use in Stan
  # Args: see data_predictor
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
    Xmo_fun <- function(x) as.array(attr(eval2(x, data), "var"))
    Xmo <- lapply(unlist(spef$calls_mo), Xmo_fun)
    Xmo_names <- paste0("Xmo", p, "_", seq_along(Xmo))
    out <- c(out, setNames(Xmo, Xmo_names))
    compute_Jmo <- is.null(Jmo)
    if (!length(Jmo)) {
      Jmo <- as.array(ulapply(Xmo, max))
    }
    out[[paste0("Jmo", p)]] <- Jmo
    # prepare prior concentration of simplex parameters
    simo_coef <- get_simo_labels(spef)
    for (i in seq_along(simo_coef)) {
      simo_prior <- subset2(prior, 
        class = "simo", coef = simo_coef[i], ls = px
      )
      simo_prior <- simo_prior$prior
      if (isTRUE(nzchar(simo_prior))) {
        simo_prior <- eval2(simo_prior)
        if (length(simo_prior) != Jmo[i]) {
          stop2("Invalid Dirichlet prior for the simplex of coefficient '",
                simo_coef[i], "'. Expected input of length ", Jmo[i], ".")
        }
      } else {
        simo_prior <- rep(1, Jmo[i])
      }
      out[[paste0("con_simo", p, "_", i)]] <- simo_prior
    }
  }
  out
}

data_cs <- function(bterms, data) {
  # prepare data for category specific effects
  # Args: see data_predictor
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

data_Xme <- function(meef, data) {
  # prepare global data for noise free variables
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
      gr <- attributes(eval2(meef$term[K[1]], data))[["gr"]]
      Jme <- match(gr, levels)
      if (anyNA(Jme)) {
        # occurs for new levels only
        new_gr <- gr[!gr %in% levels]
        new_levels <- unique(new_gr)
        Jme[is.na(Jme)] <- match(new_gr, new_levels) + length(levels)
        # represent all indices between 1 and length(unique(Jme))
        Jme <- as.numeric(factor(Jme))
      }
      ilevels <- unique(Jme)
      out[[paste0("Nme_", i)]] <- length(ilevels)
      out[[paste0("Jme_", i)]] <- Jme
    }
    for (k in K) {
      att <- attributes(eval2(meef$term[k], data))
      Xn <- as.array(att$var)
      noise <- as.array(att$sdx)
      if (nzchar(g)) {
        for (l in ilevels) {
          # validate values of the same level
          take <- Jme %in% l
          if (length(unique(Xn[take])) > 1L ||
              length(unique(noise[take])) > 1L ) {
            stop2(
              "Measured values and measurement error should be ", 
              "unique for each group. Occured for level '", 
              levels[l], "' of group '", g, "'."
            )
          }
        }
        not_dupl_Jme <- !duplicated(Jme)
        to_order <- order(Jme[not_dupl_Jme])
        Xn <- Xn[not_dupl_Jme][to_order]
        noise <- noise[not_dupl_Jme][to_order]
      }
      out[[paste0("Xn_", k)]] <- as.array(Xn)
      out[[paste0("noise_", k)]] <- as.array(noise)
    }
  }
  out
}

data_gp <- function(bterms, data, gps = NULL, ...) {
  # prepare data for Gaussian process terms
  # Args: see data_predictor
  #   ...: passed to '.data_gp'
  out <- list()
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
    Xgp <- run(cbind, Xgp)
    if (gpef$scale[i]) {
      # scale predictor for easier specification of priors
      if (length(gps) > 0L) {
        # scale Xgp based on the original data
        Xgp <- Xgp / gps[[i]]$dmax
      } else {
        dmax <- sqrt(max(diff_quad(Xgp)))
        Xgp <- Xgp / dmax
      }
    }
    cmc <- gpef$cmc[i]
    gr <- gpef$gr[i]
    k <- gpef$k[i]
    L <- gpef$L[[i]]
    if (!isNA(k)) {
      out[[paste0("NBgp", pi)]] <- k ^ D
      Ks <- as.matrix(do.call(expand.grid, repl(seq_len(k), D)))
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
        tmp <- .data_gp(Xgp, k = k, gr = gr, sfx = sfx, Cgp = Cgp, L = L, ...)
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
      c(out) <- .data_gp(Xgp, k = k, gr = gr, sfx = pi, L = L, ...)
      if (bynum) {
        Cgp <- as.numeric(get(byvar, data))
        out[[paste0("Cgp", pi)]] <- as.array(Cgp)
      }
    }
  }
  out
}

.data_gp <- function(Xgp, k, gr, sfx, Cgp = NULL, L = NULL, rawXgp = FALSE) {
  # helper function to preparae GP related data
  # Args:
  #   Xgp: matrix of covariate values
  #   k, gr, L: see tidy_gpef
  #   sfx: suffix to put at the end of data names
  #   Cgp: optional vector of values belonging to
  #     a certain contrast of a factor 'by' variable
  #   rawXgp: a flag to indicate if the covariate matrix should be returned 
  #     without further processing; required in 'def_lscale_prior'
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
    Xgp <-  Xgp[not_dupl_Jgp, , drop = FALSE]
  }
  if (rawXgp) {
    out[[paste0("Xgp", sfx)]] <- Xgp
    return(out)
  }
  if (!isNA(k)) {
    # basis function approach requires centered variables
    Xgp <- sweep(Xgp, 2, colMeans(Xgp))
    D <- NCOL(Xgp)
    L <- choose_L(Xgp, L = L)
    Ks <- as.matrix(do.call(expand.grid, repl(seq_len(k), D)))
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

data_offset <- function(bterms, data) {
  # prepare data of offsets for use in Stan
  # Args: see data_predictor
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
    out[[paste0("offset", p)]] <- as.array(offset)
  }
  out
}

data_autocor <- function(bterms, data, Y = NULL, new = FALSE,
                         old_locations = NULL) {
  # data for autocorrelation variables
  # Args:
  #   Y: vector of response values; only required in cor_arr
  #   new: does 'data' contain new data?
  #   old_locations: optional locations for CAR models 
  #     used when fitting the model
  stopifnot(is.brmsterms(bterms))
  autocor <- bterms$autocor
  N <- nrow(data)
  out <- list()
  if (is.cor_arma(autocor) || is.cor_bsts(autocor)) {
    if (length(bterms$time$group)) {
      tgroup <- as.numeric(factor(data[[bterms$time$group]]))
    } else {
      tgroup <- rep(1, N) 
    }
  }
  if (has_arma(autocor)) {
    Kar <- get_ar(autocor)
    Kma <- get_ma(autocor)
    Karr <- get_arr(autocor)
    if (Kar || Kma) {
      # ARMA correlations (of residuals)
      out$Kar <- Kar
      out$Kma <- Kma
      if (use_cov(autocor)) {
        # data for the 'covariance' version of ARMA 
        out$N_tg <- length(unique(tgroup))
        out$begin_tg <- as.array(ulapply(unique(tgroup), match, tgroup))
        out$nobs_tg <- as.array(with(out, 
          c(if (N_tg > 1L) begin_tg[2:N_tg], N + 1) - begin_tg
        ))
        out$end_tg <- with(out, begin_tg + nobs_tg - 1)
      } else {
        # data for the 'predictor' version of ARMA
        max_lag <- max(Kar, Kma)
        out$J_lag <- as.array(rep(0, N))
        for (n in seq_len(N)) {
          for (i in seq_len(max_lag)) {
            valid_lag <- n + 1 - i > 0 && n < N && 
              tgroup[n + 1] == tgroup[n + 1 - i]
            if (valid_lag) {
              out$J_lag[n] <- i
            }
          }
        }
      }
    }
    if (Karr) {
      # ARR effects (autoregressive effects of the response)
      out$Yarr <- arr_design_matrix(Y, Karr, tgroup)
      out$Karr <- Karr
    }
  } else if (is.cor_sar(autocor)) {
    if (!identical(dim(autocor$W), rep(N, 2))) {
      stop2("Dimensions of 'W' must be equal to the number of observations.")
    }
    out$W <- autocor$W
    # simplifies code of choose_N
    out$N_tg <- 1
  } else if (is.cor_car(autocor)) {
    if (isTRUE(nzchar(bterms$time$group))) {
      loc_data <- get(bterms$time$group, data)
      locations <- levels(factor(loc_data))
      if (!is.null(old_locations)) {
        new_locations <- setdiff(locations, old_locations)
        if (length(new_locations)) {
          stop2("Cannot handle new locations in CAR models.")
        }
      } else {
        old_locations <- locations
      }
      Nloc <- length(locations)
      Jloc <- as.array(match(loc_data, old_locations))
      found_locations <- rownames(autocor$W)
      if (is.null(found_locations)) {
        stop2("Row names are required for 'W'.")
      }
      colnames(autocor$W) <- found_locations
      found <- locations %in% found_locations
      if (any(!found)) {
        stop2("Row names of 'W' do not match ", 
              "the names of the grouping levels.")
      }
      autocor$W <- autocor$W[locations, locations, drop = FALSE]
    } else {
      Nloc <- N
      Jloc <- as.array(seq_len(Nloc))
      if (!identical(dim(autocor$W), rep(Nloc, 2))) {
        if (new) {
          stop2("Cannot handle new data in CAR models ",
                "without a grouping factor.")
        } else {
          stop2("Dimensions of 'W' must be equal ", 
                "to the number of observations.") 
        }
      }
    }
    W_tmp <- autocor$W
    W_tmp[upper.tri(W_tmp)] <- NA
    edges <- which(as.matrix(W_tmp == 1), arr.ind = TRUE)
    Nneigh <- Matrix::colSums(autocor$W)
    if (any(Nneigh == 0)) {
      stop2("All locations should have at least one neighbor.")
    }
    inv_sqrt_D <- diag(1 / sqrt(Nneigh))
    eigenW <- t(inv_sqrt_D) %*% autocor$W %*% inv_sqrt_D
    eigenW <- eigen(eigenW, TRUE, only.values = TRUE)$values
    out <- c(out, nlist(
      Nloc, Jloc, Nneigh, eigenW, Nedges = nrow(edges),  
      edges1 = as.array(edges[, 1]), edges2 = as.array(edges[, 2])
    ))
  } else if (is.cor_bsts(autocor)) {
    out$tg <- as.array(tgroup)
  } else if (is.cor_fixed(autocor)) {
    V <- autocor$V
    rmd_rows <- attr(data, "na.action")
    if (!is.null(rmd_rows)) {
      V <- V[-rmd_rows, -rmd_rows, drop = FALSE]
    }
    if (nrow(V) != N) {
      stop2("'V' must have the same number of rows as 'data'.")
    }
    if (min(eigen(V)$values <= 0)) {
      stop2("'V' must be positive definite.")
    }
    out$V <- V
    # simplifies code of choose_N
    out$N_tg <- 1
  }
  if (length(out)) {
    resp <- usc(combine_prefix(bterms))
    out <- setNames(out, paste0(names(out), resp))
  }
  out
}

#' Prepare Response Data
#' 
#' Prepare data related to response variables in \pkg{brms}. 
#' Only exported for use in package development.
#' 
#' @param x An \R object.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return A named list of data related to response variables.
#' 
#' @keywords internal
#' @export
data_response <- function(x, ...) {
  UseMethod("data_response")
}

#' @export
data_response.mvbrmsterms <- function(x, old_sdata = NULL, ...) {
  out <- list()
  for (i in seq_along(x$terms)) {
    od <- old_sdata[[x$responses[i]]]
    out <- c(out, data_response(x$terms[[i]], old_sdata = od, ...))
  }
  if (x$rescor) {
    out$nresp <- length(x$responses)
    out$nrescor <- out$nresp * (out$nresp - 1) / 2
  }
  out
}

#' @export
data_response.brmsterms <- function(x, data, check_response = TRUE,
                                    not4stan = FALSE, new = FALSE,
                                    old_sdata = NULL, ...) {
  # prepare data for the response variable
  N <- nrow(data)
  Y <- model.response(model.frame(x$respform, data, na.action = na.pass))
  out <- list(Y = unname(Y))
  if (is_binary(x$family) || is_categorical(x$family)) {
    out$Y <- as_factor(out$Y, levels = old_sdata$resp_levels)
    out$Y <- as.numeric(out$Y)
    if (is_binary(x$family)) {
      out$Y <- out$Y - 1
    }
  }
  if (is_ordinal(x$family) && is.ordered(out$Y)) {
    out$Y <- as.numeric(out$Y)
  }
  if (check_response) {
    family4error <- family_names(x$family)
    if (is.mixfamily(x$family)) {
      family4error <- paste0(family4error, collapse = ", ")
      family4error <- paste0("mixture(", family4error, ")")
    }
    if (!allow_factors(x$family) && !is.numeric(out$Y)) {
      stop2("Family '", family4error, "' requires numeric responses.")
    }
    if (is_binary(x$family)) {
      if (any(!out$Y %in% c(0, 1))) {
        stop2("Family '", family4error, "' requires responses ",
              "to contain only two different values.")
      }
    }
    if (is_ordinal(x$family)) {
      if (any(!is_wholenumber(out$Y)) || any(!out$Y > 0)) {
        stop2("Family '", family4error, "' requires either positive ",
              "integers or ordered factors as responses.")
      }
    }
    if (use_int(x$family)) {
      if (!all(is_wholenumber(out$Y))) {
        stop2("Family '", family4error, "' requires integer responses.")
      }
    }
    ybounds <- family_info(x$family, "ybounds")
    closed <- family_info(x$family, "closed")
    if (is.finite(ybounds[1])) {
      if (closed[1] && min(out$Y) < ybounds[1]) {
        stop2("Family '", family4error, "' requires response greater ",
              "than or equal to ", ybounds[1], ".")
      } else if (!closed[1] && min(out$Y) <= ybounds[1]) {
        stop2("Family '", family4error, "' requires response greater ",
              "than ", round(ybounds[1], 2), ".")
      }
    }
    if (is.finite(ybounds[2])) {
      if (closed[2] && max(out$Y) > ybounds[2]) {
        stop2("Family '", family4error, "' requires response smaller ",
              "than or equal to ", ybounds[2], ".")
      } else if (!closed[2] && max(out$Y) >= ybounds[2]) {
        stop2("Family '", family4error, "' requires response smaller ",
              "than ", round(ybounds[2], 2), ".")
      }
    }
    out$Y <- as.array(out$Y)
  }
  # data for addition arguments of the response
  if (has_trials(x$family) || is.formula(x$adforms$trials)) {
    if (!length(x$adforms$trials)) {
      out$trials <- max(out$Y, na.rm = TRUE)
      if (is.finite(out$trials)) {
        message("Using the maximum response value as the number of trials.")
        warning2(
          "Using 'binomial' families without specifying 'trials' ", 
          "on the left-hand side of the model formula is deprecated."
        )
      } else if (!is.null(old_sdata$trials)) {
        out$trials <- max(old_sdata$trials)
      } else {
        stop2("Could not compute the number of trials.")
      }
    } else if (is.formula(x$adforms$trials)) {
      out$trials <- eval_rhs(x$adforms$trials, data = data)
    } else {
      stop2("Argument 'trials' is misspecified.")
    }
    if (length(out$trials) == 1L) {
      out$trials <- rep(out$trials, nrow(data))
    }
    if (max(out$trials) == 1L && !not4stan) {
      message("Only 2 levels detected so that family 'bernoulli' ",
              "might be a more efficient choice.")
    }
    if (check_response && any(out$Y > out$trials)) {
      stop2("Number of trials is smaller than ",
            "the number of events.")
    }
    out$trials <- as.array(out$trials)
  }
  if (has_cat(x$family) || is.formula(x$adforms$cat)) {
    if (!length(x$adforms$cat)) {
      if (!is.null(old_sdata$ncat)) {
        out$ncat <- old_sdata$ncat
      } else {
        out$ncat <- max(out$Y)
      }
    } else if (is.formula(x$adforms$cat)) {
      out$ncat <- eval_rhs(x$adforms$cat, data = data)
    } else {
      stop2("Argument 'cat' is misspecified.")
    }
    if (out$ncat < 2L) {
      stop2("At least two response categories are required.")
    }
    if (out$ncat == 2L && !not4stan) {
      message("Only 2 levels detected so that family 'bernoulli' ",
              "might be a more efficient choice.")
    }
    if (check_response && any(out$Y > out$ncat)) {
      stop2("Number of categories is smaller than the response ",
            "variable would suggest.")
    }
  }
  if (is.formula(x$adforms$se)) {
    out$se <- as.array(eval_rhs(x$adforms$se, data = data))
  }
  if (is.formula(x$adforms$weights)) {
    out$weights <- as.array(eval_rhs(x$adforms$weights, data = data))
  }
  if (is.formula(x$adforms$dec)) {
    out$dec <- as.array(eval_rhs(x$adforms$dec, data = data))
  }
  if (is.formula(x$adforms$cens) && check_response) {
    cens <- eval_rhs(x$adforms$cens, data = data)
    out$cens <- rm_attr(cens, "y2")
    y2 <- attr(cens, "y2")
    if (!is.null(y2)) {
      icens <- cens %in% 2
      if (any(out$Y[icens] >= y2[icens])) {
        stop2("Left censor points must be smaller than right ",
              "censor points for interval censored data.")
      }
      y2[!icens] <- 0  # not used in Stan
      out$rcens <- as.array(y2)
    }
    out$cens <- as.array(out$cens)
  }
  if (is.formula(x$adforms$trunc)) {
    out <- c(out, eval_rhs(x$adforms$trunc, data = data))
    if (length(out$lb) == 1L) {
      out$lb <- rep(out$lb, N)
    }
    if (length(out$ub) == 1L) {
      out$ub <- rep(out$ub, N)
    }
    if (length(out$lb) != N || length(out$ub) != N) {
      stop2("Invalid truncation bounds.")
    }
    inv_bounds <- out$Y < out$lb | out$Y > out$ub
    if (check_response && isTRUE(any(inv_bounds))) {
      stop2("Some responses are outside of the truncation bounds.")
    }
  }
  if (is.formula(x$adforms$mi)) {
    sdy <- get_sdy(x, data)
    if (is.null(sdy)) {
      # missings only
      which_mi <- which(is.na(out$Y))
      out$Jmi <- as.array(which_mi)
      out$Nmi <- length(out$Jmi)
    } else {
      # measurement error in the response
      if (length(sdy) == 1L) {
        sdy <- rep(sdy, length(out$Y))
      }
      if (length(sdy) != length(out$Y)) {
        stop2("'sdy' must have the same length as the response.")
      }
      # all observations will have a latent score
      which_mi <- which(is.na(out$Y) | is.infinite(sdy))
      out$Jme <- as.array(setdiff(seq_along(out$Y), which_mi))
      out$Nme <- length(out$Jme)
      out$noise <- as.array(sdy)
    }
    if (!not4stan) {
      # Stan does not allow NAs in data
      # use Inf to that min(Y) is not affected
      out$Y[which_mi] <- out$noise[which_mi] <- Inf
    }
  } 
  resp <- usc(combine_prefix(x))
  c(setNames(out, paste0(names(out), resp)),
    # specify data for autocors here in order to pass Y
    data_autocor(
      x, data = data, Y = out$Y, new = new,
      old_locations = old_sdata$locations
    )
  )
}

data_mixture <- function(bterms, prior = brmsprior()) {
  # data specific for mixture models
  stopifnot(is.brmsterms(bterms))
  out <- list()
  if (is.mixfamily(bterms$family)) {
    families <- family_names(bterms$family)
    dp_classes <- dpar_class(names(c(bterms$dpars, bterms$fdpars)))
    if (!any(dp_classes %in% "theta")) {
      # estimate mixture probabilities directly
      take <- find_rows(prior, class = "theta", resp = bterms$resp)
      theta_prior <- prior$prior[take]
      if (isTRUE(nzchar(theta_prior))) {
        theta_prior <- eval2(theta_prior)
        if (length(theta_prior) != length(families)) {
          stop2("Invalid dirichlet prior for the ", 
                "mixture probabilities 'theta'.")
        }
        out$con_theta <- theta_prior
      } else {
        out$con_theta <- rep(1, length(families)) 
      }
      p <- usc(combine_prefix(bterms))
      names(out) <- paste0(names(out), p)
    }
  }
  out
}

data_prior <- function(bterms, data, prior) {
  # data for special priors such as horseshoe and lasso
  out <- list()
  px <- check_prefix(bterms)
  p <- usc(combine_prefix(px))
  prefix <- combine_prefix(px, keep_mu = TRUE)
  special <- attr(prior, "special")[[prefix]]
  if (!is.null(special[["hs_df"]])) {
    # data for the horseshoe prior
    hs_obj_names <- paste0("hs_", 
      c("df", "df_global", "df_slab", "scale_global", "scale_slab")
    )
    hs_data <- special[hs_obj_names]
    if (is.null(special[["hs_par_ratio"]])) {
      hs_data$hs_scale_global <- special$hs_scale_global
    } else {
      hs_data$hs_scale_global <- special$hs_par_ratio / sqrt(nrow(data))
    }
    names(hs_data) <- paste0(hs_obj_names, p) 
    out <- c(out, hs_data)
  }
  if (!is.null(special[["lasso_df"]])) {
    lasso_obj_names <- paste0("lasso_", c("df", "scale"))
    lasso_data <- special[lasso_obj_names]
    names(lasso_data) <- paste0(lasso_obj_names, p) 
    out <- c(out, lasso_data)
  }
  out
}
