#' @export
data_effects.btl <- function(x, data, ranef = empty_ranef(), 
                             prior = brmsprior(), knots = NULL, 
                             not4stan = FALSE, smooths = NULL, 
                             gps = NULL, Jmo = NULL) {
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
  #   smooths: optional list of smooth objects based on the original data
  #   gps: optional list of GP objects based on the original data
  #   Jmo: optional precomputed values of Jmo for monotonic effects
  # Returns:
  #   A named list of data to be passed to Stan
  px <- check_prefix(x)
  data_fe <- data_fe(
    x, data = data, knots = knots, px = px, 
    not4stan = not4stan, smooths = smooths
  )
  data_mo <- data_mo(
    x, data = data, ranef = ranef, 
    prior = prior, Jmo = Jmo, px = px
  )
  data_re <- data_re(
    ranef, data = data, px = px, not4stan = not4stan
  )
  data_me <- data_me(x, data = data, px = px)
  data_cs <- data_cs(x, data = data, px = px)
  data_gp <- data_gp(x, data = data, px = px, gps = gps)
  data_offset <- data_offset(x, data = data, px = px)
  data_prior <- data_prior(prior, data = data, px = px)
  c(data_fe, data_mo, data_re, data_me, data_cs, 
    data_gp, data_offset, data_prior)
}

#' @export 
data_effects.btnl <- function(x, data, ranef = empty_ranef(), 
                              prior = brmsprior(), knots = NULL, 
                              px = list(), not4stan = FALSE, 
                              smooths = NULL, gps = NULL, Jmo = NULL) {
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
      out <- c(out, setNames(
        as.list(as.data.frame(C)), 
        paste0("C", p, "_", seq_len(ncol(C)))
      ))
    }
  }
  for (nlp in names(x$nlpars)) {
    out <- c(out,
      data_effects(
        x$nlpars[[nlp]], data, ranef = ranef,
        prior = prior, knots = knots, not4stan = not4stan, 
        smooths = smooths[[nlp]], gps = gps[[nlp]], 
        Jmo = Jmo[[nlp]]
      )
    )
  }
  out
}

data_fe <- function(bterms, data, knots = NULL, px = list(), 
                    not4stan = FALSE, smooths = NULL) {
  # prepare data for fixed effects for use in Stan 
  # Args: see data_effects
  out <- list()
  p <- usc(combine_prefix(px))
  is_ordinal <- is_ordinal(bterms$family)
  is_bsts <- inherits(bterms$autocor, "cor_bsts")
  # the intercept is removed inside the Stan code for ordinal models
  cols2remove <- if (is_ordinal && not4stan || is_bsts) "(Intercept)"
  X <- get_model_matrix(rhs(bterms$fe), data, cols2remove = cols2remove)
  sm_labels <- get_sm_labels(bterms)
  if (length(sm_labels)) {
    stopifnot(is.null(smooths) || length(smooths) == length(sm_labels))
    Xs <- Zs <- list()
    new_smooths <- is.null(smooths)
    if (new_smooths) {
      smooths <- named_list(sm_labels)
      for (i in seq_along(sm_labels)) {
        smooths[[i]] <- mgcv::smoothCon(
          eval_smooth(sm_labels[i]), data = data, 
          knots = knots, absorb.cons = TRUE
        )
      }
    }
    by_levels <- named_list(sm_labels)
    ns <- 0
    for (i in seq_along(smooths)) {
      # may contain multiple terms when 'by' is a factor
      for (j in seq_along(smooths[[i]])) {
        ns <- ns + 1
        sm <- smooths[[i]][[j]]
        if (length(sm$by.level)) {
          by_levels[[i]][j] <- sm$by.level
        }
        if (!new_smooths) {
          sm$X <- mgcv::PredictMat(sm, rm_attr(data, "terms"))
        }
        rasm <- mgcv::smooth2random(sm, names(data))
        Xs[[ns]] <- rasm$Xf
        if (ncol(Xs[[ns]])) {
          colnames(Xs[[ns]]) <- paste0(sm$label, "_", seq_len(ncol(Xs[[ns]])))
        }
        Zs <- lapply(rasm$rand, attr, "Xr")
        Zs <- setNames(Zs, paste0("Zs", p, "_", ns, "_", seq_along(Zs)))
        knots <- list(length(Zs), as.array(ulapply(Zs, ncol)))
        knots <- setNames(knots, paste0(c("nb", "knots"), p, "_", ns))
        out <- c(out, knots, Zs)
      }
    }
    X <- cbind(X, do.call(cbind, Xs))
    scols <- lapply(Xs, function(x) which(colnames(X) %in% colnames(x)))
    X <- structure(X, smooth_cols = scols, by_levels = by_levels)
    colnames(X) <- rename(colnames(X))
  }
  avoid_dpars(colnames(X), bterms = bterms)
  c(out, setNames(list(ncol(X), X), paste0(c("K", "X"), p)))
}

data_mo <- function(bterms, data, ranef = empty_ranef(),
                    prior = brmsprior(), px = list(),
                    Jmo = NULL) {
  # prepare data for monotonic effects for use in Stan
  # Args: see data_effects
  p <- usc(combine_prefix(px))
  out <- list()
  if (is.formula(bterms[["mo"]])) {
    Xmo <- mo_design_matrix(bterms$mo, data, check = is.null(Jmo))
    avoid_dpars(colnames(Xmo), bterms = bterms)
    if (is.null(Jmo)) {
      Jmo <- as.array(apply(Xmo, 2, max))
    }
    out <- c(out, 
      setNames(
        list(ncol(Xmo), Xmo, Jmo), 
        paste0(c("Kmo", "Xmo", "Jmo"), p)
      )
    )
    # validate and assign vectors for dirichlet prior
    monef <- colnames(Xmo)
    for (i in seq_along(monef)) {
      simplex_prior <- subset2(prior,
        class = "simplex", coef = monef[i], ls = px
      )
      simplex_prior <- simplex_prior$prior
      if (isTRUE(nzchar(simplex_prior))) {
        simplex_prior <- eval2(simplex_prior)
        if (length(simplex_prior) != Jmo[i]) {
          stop2("Invalid Dirichlet prior for the simplex of coefficient '",
                monef[i], "'. Expected input of length ", Jmo[i], ".")
        }
        out[[paste0("con_simplex", p, "_", i)]] <- simplex_prior
      } else {
        out[[paste0("con_simplex", p, "_", i)]] <- rep(1, Jmo[i]) 
      }
    }
  }
  out
}

data_re <- function(ranef, data, px = list(), not4stan = FALSE) {
  # prepare data for group-level effects for use in Stan
  # Args: see data_effects
  out <- list()
  take <- find_rows(ranef, ls = px) &
    !find_rows(ranef, type = c("mo", "me"))
  # take <- ranef$nlpar == nlpar & !ranef$type %in% c("mo", "me")
  ranef <- ranef[take, ]
  if (nrow(ranef)) {
    Z <- lapply(ranef[!duplicated(ranef$gn), ]$form, 
                get_model_matrix, data = data)
    gn <- unique(ranef$gn)
    for (i in seq_along(gn)) {
      r <- subset2(ranef, gn = gn[i])
      idp <- paste0(r$id[1], usc(combine_prefix(px)))
      if (isTRUE(not4stan)) {
        # for internal use in S3 methods
        if (ncol(Z[[i]]) == 1L) {
          Z[[i]] <- as.vector(Z[[i]])
        }
        Zname <- paste0("Z_", gn[i])
        out <- c(out, setNames(Z[i], Zname))
      } else {
        if (r$type[1] == "cs") {
          ncatM1 <- nrow(r) / ncol(Z[[i]])
          Z_temp <- vector("list", ncol(Z[[i]]))
          for (k in seq_along(Z_temp)) {
            Z_temp[[k]] <- replicate(ncatM1, Z[[i]][, k])
          }
          Z[[i]] <- do.call(cbind, Z_temp)
        }
        Zname <- paste0("Z_", idp, "_", r$cn)
        for (j in seq_len(ncol(Z[[i]]))) {
          out <- c(out, setNames(list(as.array(Z[[i]][, j])), Zname[j]))
        }
      }
    }
  }
  out
}

data_gr <- function(ranef, data, cov_ranef = NULL) {
  # compute data specific for each group-level-ID
  # Args:
  #   ranef: data.frame returned by tidy_ranef
  #   data: the model.frame
  #   cov_ranef: name list of user-defined covariance matrices
  out <- list()
  ids <- unique(ranef$id)
  for (id in ids) {
    id_ranef <- ranef[ranef$id == id, ]
    nranef <- nrow(id_ranef)
    group <- id_ranef$group[1]
    levels <- attr(ranef, "levels")[[group]]
    if (id_ranef$gtype[1] == "mm") {
      gs <- id_ranef$gcall[[1]]$groups
      ngs <- length(gs)
      weights <- id_ranef$gcall[[1]]$weights
      if (is.formula(weights)) {
        scale <- isTRUE(attr(weights, "scale"))
        weights <- as.matrix(eval_rhs(weights, data))
        if (!identical(dim(weights), c(nrow(data), ngs))) {
          stop2("Grouping structure 'mm' expects 'weights' to be a matrix ", 
                "with as many columns as grouping factors.")
        }
        if (scale) {
          if (any(weights < 0)) {
            stop2("Cannot scale negative weights.")
          }         
          weights <- sweep(weights, 1, rowSums(weights), "/")
        }
      } else {
        # all members get equal weights by default
        weights <- matrix(1 / ngs, nrow = nrow(data), ncol = ngs)
      }
      for (i in seq_along(gs)) {
        temp <- list(as.array(match(get(gs[i], data), levels)), weights[, i])
        out <- c(out, setNames(temp, paste0(c("J_", "W_"), id, "_", i)))
      }
    } else {
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
    }
    temp <- list(length(levels), nranef, nranef * (nranef - 1) / 2)
    out <- c(out, setNames(temp, paste0(c("N_", "M_", "NC_"), id)))
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

data_cs <- function(bterms, data, px = list()) {
  # prepare data for category specific effects for use in Stan
  # Args:
  #   bterms: object of class brmsterms
  #   data: the data passed by the user
  out <- list()
  if (length(all_terms(bterms[["cs"]]))) {
    Xcs <- get_model_matrix(bterms$cs, data)
    avoid_dpars(colnames(Xcs), bterms = bterms)
    out <- c(out, list(Kcs = ncol(Xcs), Xcs = Xcs))
  }
  out
}

data_me <- function(bterms, data, px = list()) {
  # prepare data of meausurement error variables for use in Stan
  # Args:
  #   bterms: object of class brmsterms
  #   data: the data passed by the user
  out <- list()
  meef <- get_me_labels(bterms, data)
  if (length(meef)) {
    p <- usc(combine_prefix(px))
    Cme <- get_model_matrix(bterms$me, data)
    avoid_dpars(colnames(Cme), bterms = bterms)
    Cme <- Cme[, attr(meef, "not_one"), drop = FALSE]
    Cme <- lapply(seq_len(ncol(Cme)), function(i) Cme[, i])
    if (length(Cme)) {
      Cme <- setNames(Cme, paste0("Cme", p, "_", seq_along(Cme)))
    }
    uni_me <- attr(meef, "uni_me")
    Xn <- noise <- named_list(uni_me)
    for (i in seq_along(uni_me)) {
      temp <- eval2(uni_me[i], data)
      Xn[[i]] <- as.array(attr(temp, "var"))
      noise[[i]] <- as.array(attr(temp, "noise"))
    }
    names(Xn) <- paste0("Xn", p, "_", seq_along(Xn))
    names(noise) <- paste0("noise", p, "_", seq_along(Xn))
    Kme <- setNames(list(length(meef)), paste0("Kme", p))
    out <- c(out, Xn, noise, Cme, Kme)
  }
  out
}

data_gp <- function(bterms, data, px = list(), gps = NULL) {
  # prepare data for Gaussian process terms
  # Args:
  #   see data_effects.btl
  out <- list()
  gpef <- get_gp_labels(bterms)
  if (length(gpef)) {
    p <- usc(combine_prefix(px))
    for (i in seq_along(gpef)) {
      pi <- paste0(p, "_", i)
      gp <- eval2(gpef[i])
      Xgp <- lapply(gp$term, eval2, data)
      out[[paste0("Mgp", pi)]] <- length(Xgp)
      invalid <- ulapply(Xgp, function(x)
        !is.numeric(x) || isTRUE(length(dim(x)) > 1L)
      )
      if (any(invalid)) {
        stop2("Predictors of Gaussian processes should be numeric vectors.")
      }
      Xgp <- do.call(cbind, Xgp)
      if (gp$scale) {
        # scale predictor for easier specification of priors
        if (!is.null(gps)) {
          # scale Xgp based on the original data
          Xgp <- Xgp / gps[[i]]$dmax
        } else {
          dmax <- sqrt(max(diff_quad(Xgp)))
          Xgp <- Xgp / dmax
        }
      }
      out[[paste0("Xgp", pi)]] <- Xgp
      out[[paste0("Kgp", pi)]] <- 1L
      if (gp$by != "NA") {
        Cgp <- get(gp$by, data)
        if (is.numeric(Cgp)) {
          out[[paste0("Cgp", pi)]] <- Cgp
        } else {
          Cgp <- factor(Cgp)
          lCgp <- levels(Cgp)
          Jgp <- lapply(lCgp, function(x) which(Cgp == x))
          out[[paste0("Kgp", pi)]] <- length(Jgp)
          out[[paste0("Igp", pi)]] <- lengths(Jgp)
          Jgp_names <- paste0("Jgp", pi, "_", seq_along(Jgp))
          out <- c(out, setNames(Jgp, Jgp_names))
        }
      }
    }
  }
  out
} 

data_offset <- function(bterms, data, px = list()) {
  # prepare data of offsets for use in Stan
  out <- list()
  if (is.formula(bterms$offset)) {
    p <- usc(combine_prefix(px))
    mf <- model.frame(bterms$offset, rm_attr(data, "terms"))
    out[[paste0("offset", p)]] <- model.offset(mf)
  }
  out
}

data_mixture <- function(bterms, prior = brmsprior()) {
  # data specific for mixture models
  stopifnot(is.brmsterms(bterms))
  out <- list()
  if (is.mixfamily(bterms$family)) {
    families <- family_names(bterms$family)
    ap_classes <- dpar_class(
      names(c(bterms$dpars, bterms$fdpars))
    )
    if (!any(ap_classes %in% "theta")) {
      # estimate mixture probabilities directly
      take <- prior$class == "theta"
      theta_prior <- prior$prior[take]
      if (isTRUE(nzchar(theta_prior))) {
        theta_prior <- eval2(theta_prior)
        if (length(theta_prior) != length(families)) {
          stop2("Invalid dirichlet prior for the ", 
                "mixture probabilities 'theta'.")
        }
        out[["con_theta"]] <- theta_prior
      } else {
        out[["con_theta"]] <- rep(1, length(families)) 
      }
    }
  }
  out
}

data_prior <- function(prior, data, px = list()) {
  # data for special priors such as horseshoe and lasso
  out <- list()
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
