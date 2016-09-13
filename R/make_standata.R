#' Data for \pkg{brms} Models
#' 
#' Generate data for \pkg{brms} models to be passed to \pkg{Stan}
#'
#' @inheritParams brm
#' @param control A named list currently for internal usage only
#' @param ... Other potential arguments
#' 
#' @aliases brmdata
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
                          prior = NULL, autocor = NULL, nonlinear = NULL, 
                          partial = NULL, cov_ranef = NULL, 
                          sample_prior = FALSE, knots = NULL, 
                          control = list(), ...) {
  # internal control arguments:
  #   is_newdata: is make_standata is called with new data?
  #   not4stan: is make_standata called for use in S3 methods?
  #   save_order: should the initial order of the data be saved?
  #   omit_response: omit checking of the response?
  #   ntrials, ncat, Jm: standata based on the original data
  dots <- list(...)
  not4stan <- isTRUE(control$not4stan)
  is_newdata <- isTRUE(control$is_newdata)
  # use deprecated arguments if specified
  cov_ranef <- use_alias(cov_ranef, dots$cov.ranef, warn = FALSE)
  # some input checks 
  if (!(is.null(data) || is.list(data)))
    stop("argument 'data' must be a data.frame or list", call. = FALSE)
  family <- check_family(family)
  formula <- update_formula(formula, data = data, family = family,
                            partial = partial, nonlinear = nonlinear)
  old_mv <- isTRUE(attr(formula, "old_mv"))
  autocor <- check_autocor(autocor)
  is_linear <- is.linear(family)
  is_ordinal <- is.ordinal(family)
  is_count <- is.count(family)
  is_forked <- is.forked(family)
  is_categorical <- is.categorical(family)
  et <- extract_time(autocor$formula)
  ee <- extract_effects(formula, family = family, et$all)
  prior <- as.prior_frame(prior)
  check_prior_content(prior, family = family, warn = FALSE)
  na_action <- if (is_newdata) na.pass else na.omit
  data <- update_data(data, family = family, effects = ee, et$group,
                      drop.unused.levels = !is_newdata, 
                      na.action = na_action, knots = knots,
                      terms_attr = control$terms_attr)
  
  # sort data in case of autocorrelation models
  if (has_arma(autocor) || is(autocor, "cor_bsts")) {
    if (old_mv) {
      to_order <- rmNULL(list(data[["trait"]], data[[et$group]], 
                              data[[et$time]]))
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
    } else if (family$family %in% "bernoulli") {
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
        stop("The beta distribution requires responses between 0 and 1.", 
             call. = FALSE)
      }
    } else if (family$family %in% "von_mises") {
      if (any(standata$Y < -pi | standata$Y > pi)) {
        stop("The von_mises distribution requires ",
             "responses between -pi and pi.", call. = FALSE)
      }
    } else if (is_categorical) { 
      standata$Y <- as.numeric(as.factor(standata$Y))
      if (length(unique(standata$Y)) < 2L) {
        stop("At least two response categories are required.", call. = FALSE)
      }
    } else if (is_ordinal) {
      if (is.ordered(standata$Y)) {
        standata$Y <- as.numeric(standata$Y)
      } else if (all(is.wholenumber(standata$Y))) {
        standata$Y <- standata$Y - min(standata$Y) + 1
      } else {
        stop(paste("family", family$family, "expects either integers or",
                   "ordered factors as response variables"), call. = FALSE)
      }
      if (length(unique(standata$Y)) < 2L) {
        stop("At least two response categories are required.", call. = FALSE)
      }
    } else if (is.skewed(family) || is.lognormal(family)) {
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
    standata$Y <- as.array(standata$Y)
  }
  
  # data for various kinds of effects
  ranef <- tidy_ranef(ee, data)
  args_eff <- nlist(data, family, ranef, prior, autocor, knots, not4stan)
  if (length(ee$nonlinear)) {
    nlpars <- names(ee$nonlinear)
    # matrix of covariates appearing in the non-linear formula
    C <- get_model_matrix(ee$covars, data = data)
    if (length(all.vars(ee$covars)) != ncol(C)) {
      stop("Factors with more than two levels are not allowed as covariates",
           call. = FALSE)
    }
    standata <- c(standata, list(KC = ncol(C), C = C)) 
    for (nlp in nlpars) {
      args_eff_spec <- list(effects = ee$nonlinear[[nlp]], nlpar = nlp,
                            smooth = control$smooth[[nlp]],
                            Jm = control$Jm[[nlp]],
                            rm_intercept = FALSE)
      data_eff <- do.call(data_effects, c(args_eff_spec, args_eff))
      standata <- c(standata, data_eff)
    }
  } else {
    resp <- ee$response
    if (length(resp) > 1L && !old_mv) {
      args_eff_spec <- list(effects = ee, Jm = control$Jm[["mu"]],
                            smooth = control$smooth[["mu"]])
      for (r in resp) {
        data_eff <- do.call(data_effects, 
                            c(args_eff_spec, args_eff, nlpar = r))
        standata <- c(standata, data_eff)
        standata[[paste0("offset_", r)]] <- model.offset(data)
      }
      if (is.linear(family)) {
        standata$nresp <- length(resp) 
        standata$nrescor <- length(resp) * (length(resp) - 1) / 2 
      }
    } else {
      args_eff_spec <- list(effects = ee, Jm = control$Jm[["mu"]],
                            smooth = control$smooth[["mu"]])
      data_eff <- do.call(data_effects, c(args_eff_spec, args_eff))
      standata <- c(standata, data_eff, data_csef(ee, data = data))
      standata$offset <- model.offset(data)
    }
   
  }
  # data for predictors of scale / shape parameters
  for (ap in intersect(auxpars(), names(ee))) {
    args_eff_spec <- list(effects = ee[[ap]], nlpar = ap,
                          smooth = control$smooth[[ap]],
                          Jm = control$Jm[[ap]],
                          rm_intercept = FALSE)
    data_aux_eff <- do.call(data_effects, c(args_eff_spec, args_eff))
    standata <- c(standata, data_aux_eff)
  }
  # data for grouping factors separated after group-ID
  data_group <- data_group(ranef, data, cov_ranef = cov_ranef,
                           old_levels = control$old_levels)
  standata <- c(standata, data_group)
  
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
    if (max(standata$trials) == 1L && family$family == "binomial") 
      message(paste("Only 2 levels detected so that family bernoulli",
                    "might be a more efficient choice."))
    if (check_response && any(standata$Y > standata$trials))
      stop(paste("Number of trials is smaller than the response", 
                 "variable would suggest."), call. = FALSE)
  }
  if (has_cat(family)) {
    if (!length(ee$cat)) {
      if (!is.null(control$ncat)) {
        standata$ncat <- control$ncat
      } else {
        standata$ncat <- max(standata$Y)
      }
    } else if (is.wholenumber(ee$cat)) { 
      standata$ncat <- ee$cat
    } else stop("Addition argument 'cat' is misspecified.", call. = FALSE)
    standata$max_obs <- standata$ncat  # for backwards compatibility
    if (max(standata$ncat) == 2L) {
      message(paste("Only 2 levels detected so that family bernoulli", 
                    "might be a more efficient choice."))
    }
    if (check_response && any(standata$Y > standata$ncat)) {
      stop(paste0("Number of categories is smaller than the response", 
                  "variable would suggest."), call. = FALSE)
    }
  } 
  if (family$family == "inverse.gaussian" && check_response) {
    # save as data to reduce computation time in Stan
    if (is.formula(ee$weights) || is.formula(ee$cens)) {
      standata$log_Y <- log(standata$Y) 
    } else {
      standata$log_Y <- sum(log(standata$Y))
    }
    standata$sqrt_Y <- sqrt(standata$Y)
  }
  
  if (old_mv) {
    # deprecated as of brms 1.0.0
    # evaluate even if check_response is FALSE to ensure 
    # that N_trait is defined
    if (is_linear && length(ee$response) > 1L) {
      standata$Y <- matrix(standata$Y, ncol = length(ee$response))
      NC_trait <- ncol(standata$Y) * (ncol(standata$Y) - 1L) / 2L
      standata <- c(standata, list(N_trait = nrow(standata$Y), 
                                   K_trait = ncol(standata$Y),
                                   NC_trait = NC_trait)) 
      # for compatibility with the S3 methods of brms >= 1.0.0
      standata$nresp <- standata$K_trait
      standata$nrescor <- standata$NC_trait
    }
    if (is_forked) {
      # the second half of Y is only dummy data
      # that was put into data to make melt_data work correctly
      standata$N_trait <- nrow(data) / 2L
      standata$Y <- as.array(standata$Y[1L:standata$N_trait]) 
    }
    if (is_categorical && !isTRUE(control$old_cat == 1L)) {
      ncat1m <- standata$ncat - 1L
      standata$N_trait <- nrow(data) / ncat1m
      standata$Y <- as.array(standata$Y[1L:standata$N_trait])
      standata$J_trait <- as.array(matrix(1L:standata$N, ncol = ncat1m))
    }
  }
  
  # data for addition arguments
  if (is.formula(ee$se)) {
    standata <- c(standata, list(se = .addition(formula = ee$se, data = data)))
  }
  if (is.formula(ee$weights)) {
    standata <- c(standata, list(weights = .addition(ee$weights, data = data)))
    if (old_mv) {
      standata$weights <- standata$weights[1:standata$N_trait]
    }
  }
  if (is.formula(ee$disp)) {
    standata <- c(standata, list(disp = .addition(ee$disp, data = data)))
  }
  if (is.formula(ee$cens) && check_response) {
    standata <- c(standata, list(cens = .addition(ee$cens, data = data)))
    if (old_mv) {
      standata$cens <- standata$cens[1:standata$N_trait]
    }
  }
  if (is.formula(ee$trunc)) {
    standata <- c(standata, .addition(ee$trunc, data = data))
    if (length(standata$lb) == 1L) {
      standata$lb <- rep(standata$lb, standata$N)
    }
    if (length(standata$ub) == 1L) {
      standata$ub <- rep(standata$ub, standata$N)
    }
    if (length(standata$lb) != standata$N || 
        length(standata$ub) != standata$N) {
      stop("Invalid truncation bounds", call. = FALSE)
    }
    if (check_response && any(standata$Y < standata$lb | 
                              standata$Y > standata$ub)) {
      stop("Some responses are outside of the truncation bounds.",
           call. = FALSE)
    }
  }
  # autocorrelation variables
  if (has_arma(autocor)) {
    if (nchar(et$group)) {
      tgroup <- data[[et$group]]
    } else {
      tgroup <- rep(1, standata$N) 
    }
    Kar <- get_ar(autocor)
    Kma <- get_ma(autocor)
    Karr <- get_arr(autocor)
    if (Kar || Kma) {
      # ARMA effects (of residuals)
      standata$tg <- as.numeric(factor(tgroup))
      standata$Kar <- Kar
      standata$Kma <- Kma
      standata$Karma <- max(Kar, Kma)
      if (use_cov(autocor)) {
        # Modeling ARMA effects using a special covariance matrix
        # requires additional data
        standata$N_tg <- length(unique(standata$tg))
        standata$begin_tg <- as.array(with(standata, 
           ulapply(unique(tgroup), match, tgroup)))
        standata$nobs_tg <- as.array(with(standata, 
           c(if (N_tg > 1L) begin_tg[2:N_tg], N + 1) - begin_tg))
        standata$end_tg <- with(standata, begin_tg + nobs_tg - 1)
        if (!is.null(standata$se)) {
          standata$se2 <- standata$se^2
        } else {
          standata$se2 <- rep(0, standata$N)
        }
      } 
    }
    if (Karr) {
      if (length(ee$response) > 1L) {
        stop("ARR structure not yet implemented for multivariate models.",
             call. = FALSE)
      }
      # ARR effects (autoregressive effects of the response)
      standata$Yarr <- arr_design_matrix(Y = standata$Y, r = Karr, 
                                         group = tgroup)
      standata$Karr <- Karr
    }
  } 
  if (is(autocor, "cov_fixed")) {
    V <- autocor$V
    rmd_rows <- attr(data, "na.action")
    if (!is.null(rmd_rows)) {
      V <- V[-rmd_rows, -rmd_rows, drop = FALSE]
    }
    if (nrow(V) != nrow(data)) {
      stop("'V' must have the same number of rows as 'data'", call. = FALSE)
    }
    if (min(eigen(V)$values <= 0)) {
      stop("'V' must be positive definite", call. = FALSE)
    }
    standata$V <- V
  }
  if (is(autocor, "cor_bsts")) {
    if (length(ee$response) > 1L) {
      stop("BSTS structure not yet implemented for multivariate models.",
           call. = FALSE)
    }
    if (nchar(et$group)) {
      tgroup <- data[[et$group]]
    } else {
      tgroup <- rep(1, standata$N) 
    }
    standata$tg <- as.numeric(factor(tgroup))
  }
  standata$prior_only <- ifelse(identical(sample_prior, "only"), 1L, 0L)
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
