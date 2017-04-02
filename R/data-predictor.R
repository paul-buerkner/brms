#' @export
data_effects.btl <- function(x, data, ranef = empty_ranef(), 
                             prior = brmsprior(), knots = NULL, 
                             nlpar = "", not4stan = FALSE, 
                             smooth = NULL, Jmo = NULL) {
  # prepare data for all types of effects fpr use in Stan
  # Args:
  #   data: the data passed by the user
  #   family: the model family
  #   prior: an object of class brmsprior
  #   autocor: object of class 'cor_brms'
  #   cov_ranef: name list of user-defined covariance matrices
  #   knots: optional knot values for smoothing terms
  #   nlpar: optional character string naming a non-linear parameter
  #   not4stan: is the data for use in S3 methods only?
  #   smooth: optional list of smoothing objects based on 
  #           the original data
  #   Jmo: optional precomputed values of Jmo for monotonic effects
  # Returns:
  #   A named list of data to be passed to Stan
  nlpar <- check_nlpar(nlpar)
  data_fe <- data_fe(x, data = data, knots = knots,
                     nlpar = nlpar, not4stan = not4stan,
                     smooth = smooth)
  data_mo <- data_mo(x, data = data, ranef = ranef,
                     prior = prior, Jmo = Jmo, nlpar = nlpar)
  data_re <- data_re(ranef, data = data, nlpar = nlpar,
                     not4stan = not4stan)
  data_me <- data_me(x, data = data, nlpar = nlpar)
  data_cs <- data_cs(x, data = data, nlpar = nlpar)
  data_offset <- data_offset(x, data = data, nlpar = nlpar)
  c(data_fe, data_mo, data_re, data_me, data_cs, data_offset)
}

#' @export 
data_effects.btnl <- function(x, data, ranef = empty_ranef(), 
                              prior = brmsprior(), knots = NULL, 
                              nlpar = "", not4stan = FALSE, 
                              smooth = NULL, Jmo = NULL) {
  # prepare data for non-linear parameters for use in Stan
  # matrix of covariates appearing in the non-linear formula
  out <- list()
  C <- get_model_matrix(x$covars, data = data)
  if (length(all.vars(x$covars)) != ncol(C)) {
    stop2("Factors with more than two levels are not allowed as covariates.")
  }
  # fixes issue #127 occuring for factorial covariates
  colnames(C) <- all.vars(x$covars)
  if (not4stan) {
    out <- c(out, nlist(C))
  } else {
    # use vectors as indexing matrices in Stan is slow
    if (ncol(C)) {
      out <- c(out, setNames(
        as.list(as.data.frame(C)), 
        paste0("C_", seq_len(ncol(C)))
      ))
    }
  }
  for (nlp in names(x$nlpars)) {
    out <- c(out,
      data_effects(
        x$nlpars[[nlp]], data, ranef = ranef,
        prior = prior, knots = knots, nlpar = nlp,
        not4stan = not4stan, smooth = smooth, Jmo = Jmo
      )
    )
  }
  out
}

data_fe <- function(bterms, data, knots = NULL, nlpar = "", 
                    not4stan = FALSE, smooth = NULL) {
  # prepare data for fixed effects for use in Stan 
  # Args: see data_effects
  stopifnot(length(nlpar) == 1L)
  out <- list()
  p <- usc(nlpar, "prefix")
  is_ordinal <- is_ordinal(bterms$family)
  is_bsts <- inherits(bterms$autocor, "cor_bsts")
  # the intercept is removed inside the Stan code for ordinal models
  cols2remove <- if (is_ordinal && not4stan || is_bsts) "(Intercept)"
  X <- get_model_matrix(rhs(bterms$fe), data, cols2remove = cols2remove)
  sm_labels <- get_sm_labels(bterms)
  if (length(sm_labels)) {
    stopifnot(is.null(smooth) || length(smooth) == length(sm_labels))
    Xs <- Zs <- list()
    new_smooths <- is.null(smooth)
    if (new_smooths) {
      smooth <- named_list(sm_labels)
      for (i in seq_along(sm_labels)) {
        smooth[[i]] <- mgcv::smoothCon(
          eval_smooth(sm_labels[i]), data = data, 
          knots = knots, absorb.cons = TRUE
        )
      }
    }
    by_levels <- named_list(sm_labels)
    ns <- 0
    for (i in seq_along(smooth)) {
      # may contain multiple terms when 'by' is a factor
      for (j in seq_along(smooth[[i]])) {
        ns <- ns + 1
        sm <- smooth[[i]][[j]]
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
  avoid_auxpars(colnames(X), bterms = bterms)
  c(out, setNames(list(ncol(X), X), paste0(c("K", "X"), p)))
}

data_mo <- function(bterms, data, ranef = empty_ranef(),
                    prior = brmsprior(), nlpar = "",
                    Jmo = NULL) {
  # prepare data for monotonic effects for use in Stan
  # Args: see data_effects
  stopifnot(length(nlpar) == 1L)
  p <- if (nchar(nlpar)) paste0("_", nlpar) else ""
  out <- list()
  if (is.formula(bterms[["mo"]])) {
    Xmo <- prepare_mo_vars(bterms$mo, data, check = is.null(Jmo))
    avoid_auxpars(colnames(Xmo), bterms = bterms)
    if (is.null(Jmo)) {
      Jmo <- as.array(apply(Xmo, 2, max))
    }
    out <- c(out, 
      setNames(list(ncol(Xmo), Xmo, Jmo), paste0(c("Kmo", "Xmo", "Jmo"), p))
    )
    # validate and assign vectors for dirichlet prior
    monef <- colnames(Xmo)
    for (i in seq_along(monef)) {
      take <- prior$class == "simplex" & prior$coef == monef[i] & 
        prior$nlpar == nlpar  
      simplex_prior <- prior$prior[take]
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

data_re <- function(ranef, data, nlpar = "", not4stan = FALSE) {
  # prepare data for group-level effects for use in Stan
  # Args: see data_effects
  stopifnot(length(nlpar) == 1L)
  out <- list()
  take <- ranef$nlpar == nlpar & !ranef$type %in% c("mo", "me")
  ranef <- ranef[take, ]
  if (nrow(ranef)) {
    Z <- lapply(ranef[!duplicated(ranef$gn), ]$form, 
                get_model_matrix, data = data)
    gn <- unique(ranef$gn)
    for (i in seq_along(gn)) {
      r <- ranef[ranef$gn == gn[i], ]
      idp <- paste0(r$id[1], usc(nlpar, "prefix"))
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
        weights <- as.matrix(eval_rhs(weights, data))
        if (!identical(dim(weights), c(nrow(data), ngs))) {
          stop2("Grouping structure 'mm' expects 'weights' to be a matrix ", 
                "with as many columns as grouping factors.")
        }
        weights <- sweep(weights, 1, rowSums(weights), "/")
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

data_cs <- function(bterms, data, nlpar = "") {
  # prepare data for category specific effects for use in Stan
  # Args:
  #   bterms: object of class brmsterms
  #   data: the data passed by the user
  out <- list()
  if (length(all_terms(bterms[["cs"]]))) {
    stopifnot(!nzchar(nlpar))
    Xcs <- get_model_matrix(bterms$cs, data)
    avoid_auxpars(colnames(Xcs), bterms = bterms)
    out <- c(out, list(Kcs = ncol(Xcs), Xcs = Xcs))
  }
  out
}

data_me <- function(bterms, data, nlpar = "") {
  # prepare data of meausurement error variables for use in Stan
  # Args:
  #   bterms: object of class brmsterms
  #   data: the data passed by the user
  out <- list()
  meef <- get_me_labels(bterms, data)
  if (length(meef)) {
    p <- usc(nlpar, "prefix")
    Cme <- get_model_matrix(bterms$me, data)
    avoid_auxpars(colnames(Cme), bterms = bterms)
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

data_offset <- function(bterms, data, nlpar = "") {
  # prepare data of offsets for use in Stan
  out <- list()
  if (is.formula(bterms$offset)) {
    mf <- model.frame(bterms$offset, rm_attr(data, "terms"))
    out[[paste0("offset", usc(nlpar))]] <- model.offset(mf)
  }
  out
}

data_mixture <- function(bterms, prior = brmsprior()) {
  # data specific for mixture models
  stopifnot(is.brmsterms(bterms))
  out <- list()
  if (is.mixfamily(bterms$family)) {
    families <- family_names(bterms$family)
    ap_classes <- auxpar_class(
      names(c(bterms$auxpars, bterms$fauxpars))
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
