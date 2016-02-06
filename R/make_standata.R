#' Data for \pkg{brms} Models
#' 
#' Generate data for \pkg{brms} models to be passed to \pkg{Stan}
#'
#' @inheritParams brm
#' @param control A named list currently for internal usage only
#' @param ... Other potential arguments
#' 
#' @aliases brmdata brm.data
#' 
#' @return A named list of objects containing the required data 
#'   to fit a \pkg{brms} model with \pkg{Stan}. 
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' data1 <- make_standata(rating ~ treat + period + carry + (1|subject), 
#'                        data = inhaler, family = "cumulative")
#' names(data1)
#' 
#' data2 <- make_standata(count ~ log_Age_c + log_Base4_c * Trt_c 
#'                        + (1|patient) + (1|visit), 
#'                        data = epilepsy, family = "poisson")
#' names(data2)
#'          
#' @export
make_standata <- function(formula, data = NULL, family = "gaussian", 
                          autocor = NULL, nonlinear = NULL, partial = NULL, 
                          cov_ranef = NULL, control = NULL, ...) {
  # internal control arguments:
  #   is_newdata: logical; indicating if make_standata is called with new data
  #   keep_intercept: logical; indicating if the Intercept column
  #                   should be kept in the FE design matrix
  #   save_order: logical; should the initial order of the data be saved?
  dots <- list(...)
  # use deprecated arguments if specified
  cov_ranef <- use_alias(cov_ranef, dots$cov.ranef, warn = FALSE)
  # some input checks 
  formula <- update_formula(formula, data = data)
  autocor <- check_autocor(autocor)
  family <- check_family(family)
  is_linear <- is.linear(family)
  is_ordinal <- is.ordinal(family)
  is_count <- is.count(family)
  is_forked <- is.forked(family)
  et <- extract_time(autocor$formula)
  ee <- extract_effects(formula = formula, family = family, 
                        partial, et$all, nonlinear = nonlinear)
  na_action <- if (isTRUE(control$is_newdata)) na.pass else na.omit
  data <- update_data(data, family = family, effects = ee, et$group,
                      drop.unused.levels = !isTRUE(control$is_newdata),
                      na.action = na_action)
  
  # sort data in case of autocorrelation models
  if (has_arma(autocor)) {
    # amend if zero-inflated and hurdle models ever get 
    # autocorrelation structures as they are also using 'trait'
    if (is_forked) {
      stop("no autocorrelation allowed for this model", call. = FALSE)
    }
    if (is_linear && length(ee$response) > 1) {
      if (!grepl("^trait$|:trait$|^trait:|:trait:", et$group)) {
        stop(paste("autocorrelation structures for multiple responses must",
                   "contain 'trait' as grouping variable"), call. = FALSE)
      } else {
        to_order <- rmNULL(list(data[["trait"]], data[[et$group]], 
                                data[[et$time]]))
      }
    } else {
      to_order <- rmNULL(list(data[[et$group]], data[[et$time]]))
    }
    if (length(to_order)) {
      new_order <- do.call(order, to_order)
      data <- data[new_order, ]
      # old_order will allow to retrieve the initial order of the data
      attr(data, "old_order") <- order(new_order)
    }
  }
  
  # response variable
  standata <- list(N = nrow(data), Y = unname(model.response(data)))
  check_response <- !isTRUE(control$omit_response)
  if (check_response) {
    if (!(is_ordinal || family$family %in% c("bernoulli", "categorical")) && 
        !is.numeric(standata$Y)) {
      stop(paste("family", family$family, "expects numeric response variable"),
           call. = FALSE)
    }
    # transform and check response variable for different families
    regex_pos_int <- "(^|_)(binomial|poisson|negbinomial|geometric)$"
    if (grepl(regex_pos_int, family$family)) {
      if (!all(is.wholenumber(standata$Y)) || min(standata$Y) < 0) {
        stop(paste("family", family$family, "expects response variable", 
                   "of non-negative integers"), call. = FALSE)
      }
    } else if (family$family == "bernoulli") {
      standata$Y <- as.numeric(as.factor(standata$Y)) - 1
      if (any(!standata$Y %in% c(0,1))) {
        stop(paste("family", family$family, "expects response variable", 
                   "to contain only two different values"), call. = FALSE)
      }
    } else if (family$family %in% c("beta", "zero_inflated_beta")) {
      lower <- if (family$family == "beta") any(standata$Y <= 0)
      else any(standata$Y < 0)
      upper <- any(standata$Y >= 1)
      if (lower || upper) {
        stop("beta regression requires responses between 0 and 1", 
             call. = FALSE)
      }
    } else if (is.categorical(family)) { 
      standata$Y <- as.numeric(as.factor(standata$Y))
    } else if (is_ordinal) {
      if (is.factor(standata$Y)) {
        if (is.ordered(standata$Y)) standata$Y <- as.numeric(standata$Y)
        else stop(paste("family", family$family, "requires factored", 
                        "response variables to be ordered"), call. = FALSE)
      } else if (all(is.wholenumber(standata$Y))) {
        standata$Y <- standata$Y - min(standata$Y) + 1
      } else {
        stop(paste("family", family$family, "expects either integers or",
                   "ordered factors as response variables"), call. = FALSE)
      }
    } else if (is.skewed(family)) {
      if (min(standata$Y) <= 0) {
        stop(paste("family", family$family, "requires response variable", 
                   "to be positive"), call. = FALSE)
      }
    } else if (is.zero_inflated(family) || is.hurdle(family)) {
      if (min(standata$Y) < 0) {
        stop(paste("family", family$family, "requires response variable", 
                   "to be non-negative"), call. = FALSE)
      }
    }
  }
  # evaluate even if check_response is FALSE 
  # to ensure that N_trait is defined
  if (is_linear && length(ee$response) > 1) {
    standata$Y <- matrix(standata$Y, ncol = length(ee$response))
    NC_trait <- ncol(standata$Y) * (ncol(standata$Y) - 1) / 2
    standata <- c(standata, list(N_trait = nrow(standata$Y), 
                                 K_trait = ncol(standata$Y),
                                 NC_trait = NC_trait)) 
  } else if (is_forked) {
    # the second half of Y is not used because it is only dummy data
    # that was put into data to make melt_data work correctly
    standata$Y <- standata$Y[1:(nrow(data) / 2)] 
    standata$N_trait <- length(standata$Y)
  }
  
  # add an offset if present
  model_offset <- model.offset(data)
  if (!is.null(model_offset)) {
    standata$offset <- model_offset
  }
  
  # fixed effects data
  if (length(nonlinear)) {
    # fixed effects design matrices
    nlpars <- names(ee$nonlinear)
    fixed_list <- lapply(ee$nonlinear, function(par) par$fixed)
    X <- lapply(fixed_list, get_model_matrix, data = data)
    for (i in seq_along(nlpars)) {
      standata <- c(standata, setNames(list(ncol(X[[i]]), X[[i]]),
                                       paste0(c("K_", "X_"), nlpars[i])))
    }
    # matrix of covariances
    C <- get_model_matrix(ee$covars, data = data, rm_intercept = TRUE)
    if (length(all.vars(ee$covars)) != ncol(C)) {
      stop("Factors with more than two levels are not allowed as covariates",
           call. = FALSE)
    }
    standata <- c(standata, list(KC = ncol(C), C = C)) 
  } else {
    rm_intercept <- is_ordinal || !isTRUE(control$keep_intercept) ||
      isTRUE(attr(ee$fixed, "rsv_intercept"))
    X <- get_model_matrix(ee$fixed, data, rm_intercept = rm_intercept,
                          is_forked = is_forked)
    X_means <- colMeans(X)
    has_intercept <- attr(terms(formula), "intercept")
    if (!isTRUE(control$keep_intercept) && has_intercept) {
      # keep_intercept is TRUE when make_standata is called within S3 methods
      X <- sweep(X, 2, X_means, FUN = "-")
    }
    if (is.categorical(family)) {
      standata <- c(standata, 
                    list(Kp = ncol(X), Xp = X, Xp_means = as.array(X_means)))
    } else {
      standata <- c(standata, 
                    list(K = ncol(X), X = X, X_means = as.array(X_means)))
    }
  }
  # random effects data
  random <- get_random(ee)
  if (nrow(random)) {
    Z <- lapply(random$form, get_model_matrix, 
                data = data, is_forked = is_forked)
    r <- lapply(Z, colnames)
    ncolZ <- lapply(Z, ncol)
    # numeric levels passed to Stan
    expr <- expression(as.numeric(as.factor(get(g, data))), 
                       # number of levels
                       length(unique(get(g, data))),  
                       ncolZ[[i]],  # number of random effects
                       Z[[i]],  # random effects design matrix
                       #  number of correlations
                       ncolZ[[i]] * (ncolZ[[i]] - 1) / 2) 
    if (isTRUE(control$is_newdata)) {
      # for newdata only as levels are already defined in amend_newdata
      expr[1] <- expression(get(g, data)) 
    }
    for (i in 1:nrow(random)) {
      g <- random$group[[i]]
      if (length(nonlinear)) {
        pi <- paste0(rownames(random)[i], "_", i)
      } else pi <- i 
      name <- paste0(c("J_", "N_", "K_", "Z_", "NC_"), pi)
      if (ncolZ[[i]] == 1) {
        Z[[i]] <- as.vector(Z[[i]])
      }
      for (j in 1:length(name)) {
        standata <- c(standata, setNames(list(eval(expr[j])), name[j]))
      }
      if (g %in% names(cov_ranef)) {
        cov_mat <- as.matrix(cov_ranef[[g]])
        found_level_names <- rownames(cov_mat)
        colnames(cov_mat) <- found_level_names
        true_level_names <- sort(as.character(unique(data[[g]])))
        if (is.null(found_level_names)) 
          stop(paste("rownames are required for covariance matrix of", g),
               call. = FALSE)
        if (nrow(cov_mat) != length(true_level_names))
          stop(paste("dimension of covariance matrix of", g, "is incorrect"),
               call. = FALSE)
        if (any(sort(found_level_names) != true_level_names))
          stop(paste("rownames of covariance matrix of", g, 
                     "do not match names of the grouping levels"),
               call. = FALSE)
        if (!isSymmetric(unname(cov_mat)))
          stop(paste("covariance matrix of grouping factor", g, 
                     "is not symmetric"), call. = FALSE)
        if (min(eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values) <= 0)
          warning(paste("covariance matrix of grouping factor", g, 
                        "may not be positive definite"), call. = FALSE)
        cov_mat <- cov_mat[order(found_level_names), order(found_level_names)]
        if (length(r[[i]]) == 1 || !random$cor[[i]]) {
          # pivoting ensures that (numerically) semi-definite matrices can be used
          cov_mat <- suppressWarnings(chol(cov_mat, pivot = TRUE))
          cov_mat <- t(cov_mat[, order(attr(cov_mat, "pivot"))])
        } 
        standata <- c(standata, setNames(list(cov_mat), paste0("cov_",i)))
      }
    }
  }
  
  # addition and category specific variables
  if (is.formula(ee$se)) {
    standata <- c(standata, list(se = .addition(formula = ee$se, data = data)))
  }
  if (is.formula(ee$weights)) {
    standata <- c(standata, list(weights = .addition(formula = ee$weights, 
                                                     data = data)))
    if (is.linear(family) && length(ee$response) > 1 || is_forked) 
      standata$weights <- standata$weights[1:standata$N_trait]
  }
  if (is.formula(ee$cens) && check_response) {
    standata <- c(standata, list(cens = .addition(formula = ee$cens, 
                                                  data = data)))
    if (is.linear(family) && length(ee$response) > 1 || is_forked)
      standata$cens <- standata$cens[1:standata$N_trait]
  }
  if (is.formula(ee$trunc)) {
    standata <- c(standata, .addition(formula = ee$trunc))
    if (check_response && (min(standata$Y) < standata$lb || 
                           max(standata$Y) > standata$ub)) {
      stop("some responses are outside of the truncation boundaries",
           call. = FALSE)
    }
  }
  # data for specific families
  if (has_trials(family)) {
    if (!length(ee$trials)) {
      if (!is.null(control$trials)) {
        standata$trials <- control$trials
      } else {
        standata$trials <- max(standata$Y) 
      }
    } else if (is.wholenumber(ee$trials)) {
      standata$trials <- ee$trials
    } else if (is.formula(ee$trials)) {
      standata$trials <- .addition(formula = ee$trials, data = data)
    } else stop("Response part of formula is invalid.")
    standata$max_obs <- standata$trials  # for backwards compatibility
    if (max(standata$trials) == 1 && family$family == "binomial") 
      message(paste("Only 2 levels detected so that family bernoulli",
                    "might be a more efficient choice."))
    if (check_response && any(standata$Y > standata$trials))
      stop(paste("Number of trials is smaller than the response", 
                 "variable would suggest."), call. = FALSE)
  } else if (has_cat(family)) {
    if (!length(ee$cat)) {
      if (!is.null(control$ncat)) {
        standata$ncat <- control$ncat
      } else {
        standata$ncat <- max(standata$Y)
      }
    } else if (is.wholenumber(ee$cat)) { 
      standata$ncat <- ee$cat
    } else if (is.formula(ee$cat)) {
      warning("observations may no longer have different numbers of categories",
              call. = FALSE)
      standata$ncat <- max(.addition(formula = ee$cat, data = data))
    } else stop("Response part of formula is invalid.")
    standata$max_obs <- standata$ncat  # for backwards compatibility
    if (max(standata$ncat) == 2) {
      message(paste("Only 2 levels detected so that family bernoulli", 
                    "might be a more efficient choice."))
    }
    if (check_response && any(standata$Y > standata$ncat)) {
      stop(paste0("Number of categories is smaller than the response", 
                  "variable would suggest."), call. = FALSE)
    }
  } else if (family$family == "inverse.gaussian" && check_response) {
    # save as data to reduce computation time in Stan
    if (is.formula(ee[c("weights", "cens")])) {
      standata$log_Y <- log(standata$Y) 
    } else {
      standata$log_Y <- sum(log(standata$Y))
    }
    standata$sqrt_Y <- sqrt(standata$Y)
  } 
  
  # get data for category specific effects
  if (is.formula(partial)) {
    if (family$family %in% c("sratio","cratio","acat")) {
      Xp <- get_model_matrix(partial, data, rm_intercept = TRUE)
      standata <- c(standata, list(Kp = ncol(Xp), Xp = Xp))
      fp <- intersect(colnames(X), colnames(Xp))
      if (length(fp))
        stop(paste("Variables cannot be modeled as fixed and", 
                   "category specific effects at the same time.", 
                   "\nError occured for variables:", 
                   paste(fp, collapse = ", ")), call. = FALSE)
    } else {
      stop(paste("category specific effects are only meaningful for families", 
                 "'sratio', 'cratio', and 'acat'"), call. = FALSE)
    }
  } else if (!is.null(partial)) {
    stop("Argument partial must be a formula")
  }
  
  # autocorrelation variables
  if (has_arma(autocor)) {
    tgroup <- data[[et$group]]
    if (is.null(tgroup)) {
      tgroup <- rep(1, standata$N) 
    }
    Kar <- get_ar(autocor)
    Kma <- get_ma(autocor)
    Karr <- get_arr(autocor)
    if (Kar || Kma) {
      # ARMA effects (of residuals)
      standata$tgroup <- as.numeric(as.factor(tgroup))
      standata$E_pre <- matrix(0, nrow = standata$N, ncol = max(Kar, Kma))
      standata$Kar <- Kar
      standata$Kma <- Kma
      standata$Karma <- max(Kar, Kma)
      if (use_cov(autocor)) {
        # Modeling ARMA effects using a special covariance matrix
        # requires additional data
        standata$N_tg <- length(unique(standata$tgroup))
        standata$begin_tg <- as.array(with(standata, 
           ulapply(unique(tgroup), match, tgroup)))
        standata$nrows_tg <- as.array(with(standata, 
           c(if (N_tg > 1) begin_tg[2:N_tg], N + 1) - begin_tg))
        if (!is.null(standata$se)) {
          standata$squared_se <- standata$se^2
        } else {
          standata$squared_se <- rep(0, standata$N)
        }
      } 
    }
    if (Karr) {
      # ARR effects (autoregressive effects of the response)
      standata$Yarr <- arr_design_matrix(Y = standata$Y, r = Karr, 
                                         group = tgroup)
      standata$Karr <- Karr
    }
  } 
  if (isTRUE(control$save_order)) {
    attr(standata, "old_order") <- attr(data, "old_order")
  }
  standata
}  

#' @export
brmdata <- function(formula, data = NULL, family = "gaussian", 
                    autocor = NULL, partial = NULL, 
                    cov_ranef = NULL, ...)  {
  # deprectated alias of make_standata
  make_standata(formula = formula, data = data, 
                family = family, autocor = autocor,
                partial = partial, cov_ranef = cov_ranef, ...)
}

#' @export
brm.data <- function(formula, data = NULL, family = "gaussian", 
                     autocor = NULL, partial = NULL, 
                     cov_ranef = NULL, ...)  {
  # deprectated alias of make_standata
  make_standata(formula = formula, data = data, 
                family = family, autocor = autocor,
                partial = partial, cov_ranef = cov_ranef, ...)
}
