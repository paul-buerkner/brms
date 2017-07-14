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
make_standata <- function(formula, data, family = gaussian(), 
                          prior = NULL, autocor = NULL, nonlinear = NULL, 
                          cov_ranef = NULL, sample_prior = FALSE, 
                          knots = NULL, control = list(), ...) {
  # internal control arguments:
  #   is_newdata: is make_standata is called with new data?
  #   not4stan: is make_standata called for use in S3 methods?
  #   save_order: should the initial order of the data be saved?
  #   omit_response: omit checking of the response?
  #   only_response: only compute data related to the response?
  #   ntrials, ncat, Jmo: standata based on the original data
  dots <- list(...)
  not4stan <- isTRUE(control$not4stan)
  is_newdata <- isTRUE(control$is_newdata)
  # use deprecated arguments if specified
  cov_ranef <- use_alias(cov_ranef, dots$cov.ranef, warn = FALSE)
  # some input checks
  formula <- amend_formula(formula, data = data, family = family,
                           nonlinear = nonlinear)
  family <- formula$family
  old_mv <- isTRUE(formula[["old_mv"]])
  autocor <- check_autocor(autocor)
  is_linear <- is_linear(family)
  is_ordinal <- is_ordinal(family)
  is_count <- is_count(family)
  is_forked <- is_forked(family)
  is_categorical <- is_categorical(family)
  bterms <- parse_bf(formula, family = family, autocor = autocor)
  check_prior_content(prior, family = family, warn = FALSE)
  na_action <- if (is_newdata) na.pass else na.omit
  data <- update_data(
    data, bterms = bterms, na.action = na_action, 
    drop.unused.levels = !is_newdata, knots = knots,
    terms_attr = control$terms_attr
  )
  
  # sort data in case of autocorrelation models
  if (has_arma(autocor) || is(autocor, "cor_bsts")) {
    if (old_mv) {
      to_order <- rmNULL(list(
        data[["trait"]], 
        data[[bterms$time$group]], 
        data[[bterms$time$time]]
      ))
    } else {
      to_order <- rmNULL(list(
        data[[bterms$time$group]], 
        data[[bterms$time$time]]
      ))
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
  families <- family_names(family)
  if (is.mixfamily(family)) {
    family4error <- paste0(families, collapse = ", ")
    family4error <- paste0("mixture(", family4error, ")")
  } else {
    family4error <- families
  }
  if (check_response) {
    factors_allowed <- is_ordinal || 
      any(families %in% c("bernoulli", "categorical"))
    if (!factors_allowed && !is.numeric(standata$Y)) {
      stop2("Family '", family4error, "' requires numeric responses.")
    }
    # transform and check response variable for different families
    regex_pos_int <- "(^|_)(binomial|poisson|negbinomial|geometric)$"
    if (any(grepl(regex_pos_int, families))) {
      if (!all(is_wholenumber(standata$Y)) || min(standata$Y) < 0) {
        stop2("Family '", family4error, "' requires responses ", 
              "to be non-negative integers.")
      }
    } else if (any(families %in% "bernoulli")) {
      standata$Y <- as.numeric(as.factor(standata$Y)) - 1
      if (any(!standata$Y %in% c(0, 1))) {
        stop2("Family '", family4error, "' requires responses ", 
              "to contain only two different values.")
      }
    } else if (any(grepl("(^|_)beta$", families))) {
      if (any(families %in% "beta")) {
        lower <- any(standata$Y <= 0)
      } else {
        lower <- any(standata$Y < 0) 
      } 
      if (any(families %in% "zero_one_inflated_beta")) {
        upper <- any(standata$Y > 1) 
      } else {
        upper <- any(standata$Y >= 1) 
      }
      if (lower || upper) {
        stop2("Family '", family4error, "' requires responses ", 
              "between 0 and 1.")
      }
    } else if (any(families %in% "von_mises")) {
      if (any(standata$Y < -pi | standata$Y > pi)) {
        stop2("Family '", family4error, "' requires responses ",
              "between -pi and pi.")
      }
    } else if (is_categorical) { 
      standata$Y <- as.numeric(factor(standata$Y))
      if (length(unique(standata$Y)) < 3L) {
        stop2("At least three response categories are required.")
      }
    } else if (is_ordinal) {
      if (is.ordered(standata$Y)) {
        standata$Y <- as.numeric(standata$Y)
      } else if (all(is_wholenumber(standata$Y))) {
        standata$Y <- standata$Y - min(standata$Y) + 1
      } else {
        stop2("Family '", family4error, "' requires either integers or ",
              "ordered factors as responses.")
      }
      if (length(unique(standata$Y)) < 2L) {
        stop2("At least two response categories are required.")
      }
    } else if (is_skewed(family) || is_lognormal(family) || is_wiener(family)) {
      if (min(standata$Y) <= 0) {
        stop2("Family '", family4error, "' requires responses ", 
              "to be positive.")
      }
    } else if (is_zero_inflated(family) || is_hurdle(family)) {
      if (min(standata$Y) < 0) {
        stop2("Family '", family4error, "' requires responses ", 
              "to be non-negative.")
      }
    }
    standata$Y <- as.array(standata$Y)
  }
  
  # data for various kinds of effects
  only_response <- isTRUE(control$only_response)
  if (!only_response) {
    ranef <- tidy_ranef(bterms, data, ncat = control$ncat, 
                        old_levels = control$old_levels)
    args_eff <- nlist(data, ranef, prior, knots, not4stan)
    resp <- bterms$response
    if (length(resp) > 1L && !old_mv) {
      args_eff_spec <- list(
        x = bterms$auxpars[["mu"]], smooths = control$smooths[["mu"]],
        gps = control$gps[["mu"]], Jmo = control$Jmo[["mu"]]
      )
      bterms$auxpars[["mu"]] <- NULL
      for (r in resp) {
        data_eff <- do.call(
          data_effects, c(args_eff_spec, args_eff, nlpar = r)
        )
        standata <- c(standata, data_eff)
      }
      if (is_linear(family)) {
        standata$nresp <- length(resp)
        standata$nrescor <- length(resp) * (length(resp) - 1) / 2
        colnames(standata$Y) <- resp
      }
    }
    # data for predictors of auxiliary parameters
    for (ap in names(bterms$auxpars)) {
      args_eff_spec <- list(
        x = bterms$auxpars[[ap]], nlpar = ap, 
        smooths = control$smooths[[ap]],
        gps = control$gps[[ap]], Jmo = control$Jmo[[ap]]
      )
      data_aux_eff <- do.call(data_effects, c(args_eff_spec, args_eff))
      standata <- c(standata, data_aux_eff)
    }
    for (ap in names(bterms$fauxpars)) {
      standata[[ap]] <- bterms$fauxpars[[ap]]
    }
    standata <- c(standata,
      data_gr(ranef, data, cov_ranef = cov_ranef),
      data_mixture(bterms, prior = prior)
    )
  }
  
  # data for specific families
  if (has_trials(family)) {
    if (!length(bterms$adforms$trials)) {
      if (!is.null(control$trials)) {
        standata$trials <- control$trials
      } else {
        message("Using the maximum of the response ", 
                "variable as the number of trials.")
        standata$trials <- max(standata$Y) 
      }
    } else if (is.formula(bterms$adforms$trials)) {
      standata$trials <- eval_rhs(bterms$adforms$trials, data = data)
    } else {
      stop2("Argument 'trials' is misspecified.")
    }
    if (length(standata$trials) == 1L) {
      standata$trials <- rep(standata$trials, nrow(data))
    }
    if (max(standata$trials) == 1L && !not4stan) {
      message("Only 2 levels detected so that family 'bernoulli' ",
              "might be a more efficient choice.")
    }
    if (check_response && any(standata$Y > standata$trials)) {
      stop2("Number of trials is smaller than the response ", 
            "variable would suggest.")
    }
    standata$trials <- as.array(standata$trials)
  }
  if (has_cat(family)) {
    if (!length(bterms$adforms$cat)) {
      if (!is.null(control$ncat)) {
        standata$ncat <- control$ncat
      } else {
        standata$ncat <- length(unique(standata$Y))
      }
    } else if (is.formula(bterms$adforms$cat)) { 
      standata$ncat <- eval_rhs(bterms$adforms$cat, data = data)
    } else {
      stop2("Argument 'cat' is misspecified.")
    }
    if (max(standata$ncat) == 2L) {
      message("Only 2 levels detected so that family 'bernoulli' ",
              "might be a more efficient choice.")
    }
    if (check_response && any(standata$Y > standata$ncat)) {
      stop2("Number of categories is smaller than the response ", 
            "variable would suggest.")
    }
  }
  
  # data for addition arguments
  if (is.formula(bterms$adforms$se)) {
    standata[["se"]] <- as.array(eval_rhs(bterms$adforms$se, data = data))
  }
  if (is.formula(bterms$adforms$weights)) {
    standata[["weights"]] <- 
      as.array(eval_rhs(bterms$adforms$weights, data = data))
    if (old_mv) {
      standata$weights <- standata$weights[1:standata$N_trait]
    }
  }
  if (is.formula(bterms$adforms$disp)) {
    standata[["disp"]] <- as.array(eval_rhs(bterms$adforms$disp, data = data))
  }
  if (is.formula(bterms$adforms$dec)) {
    standata[["dec"]] <- as.array(eval_rhs(bterms$adforms$dec, data = data))
  }
  if (is.formula(bterms$adforms$cens) && check_response) {
    cens <- eval_rhs(bterms$adforms$cens, data = data)
    standata$cens <- rm_attr(cens, "y2")
    y2 <- attr(cens, "y2")
    if (!is.null(y2)) {
      icens <- cens %in% 2
      if (any(standata$Y[icens] >= y2[icens])) {
        stop2("Left censor points must be smaller than right ", 
              "censor points for interval censored data.")
      }
      y2[!icens] <- 0  # not used in Stan
      standata$rcens <- as.array(y2)
    }
    standata$cens <- as.array(standata$cens)
    if (old_mv) {
      standata$cens <- standata$cens[1:standata$N_trait]
    }
  }
  if (is.formula(bterms$adforms$trunc)) {
    standata <- c(standata, eval_rhs(bterms$adforms$trunc, data = data))
    if (length(standata$lb) == 1L) {
      standata$lb <- rep(standata$lb, standata$N)
    }
    if (length(standata$ub) == 1L) {
      standata$ub <- rep(standata$ub, standata$N)
    }
    if (length(standata$lb) != standata$N || 
        length(standata$ub) != standata$N) {
      stop2("Invalid truncation bounds.")
    }
    inv_bounds <- standata$Y < standata$lb | standata$Y > standata$ub
    if (check_response && any(inv_bounds)) {
      stop2("Some responses are outside of the truncation bounds.")
    }
  }
  
  # autocorrelation variables
  if (!only_response) {
    if (has_arma(autocor)) {
      if (nchar(bterms$time$group)) {
        tgroup <- data[[bterms$time$group]]
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
          # data for covariance matrices of ARMA effects 
          standata$N_tg <- length(unique(standata$tg))
          standata$begin_tg <- as.array(with(standata, 
            ulapply(unique(tgroup), match, tgroup)
          ))
          standata$nobs_tg <- as.array(with(standata, 
            c(if (N_tg > 1L) begin_tg[2:N_tg], N + 1) - begin_tg
          ))
          standata$end_tg <- with(standata, begin_tg + nobs_tg - 1)
        } 
      }
      if (Karr) {
        if (length(bterms$response) > 1L) {
          stop2("ARR structure not yet implemented for multivariate models.")
        }
        # ARR effects (autoregressive effects of the response)
        standata$Yarr <- arr_design_matrix(standata$Y, Karr, tgroup)
        standata$Karr <- Karr
      }
    } else if (is.cor_sar(autocor)) {
      if (!identical(dim(autocor$W), rep(standata$N, 2))) {
        stop2("Dimensions of 'W' must be equal to the number of observations.")
      }
      standata$W <- autocor$W
      # simplifies code of choose_N
      standata$N_tg <- 1
    } else if (is.cor_car(autocor)) {
      if (isTRUE(nzchar(bterms$time$group))) {
        loc_data <- get(bterms$time$group, data)
        locations <- levels(factor(loc_data))
        if (!is.null(control$old_locations)) {
          old_locations <- control$old_locations
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
        Nloc <- standata$N
        Jloc <- as.array(seq_len(Nloc))
        if (!identical(dim(autocor$W), rep(Nloc, 2))) {
          if (is_newdata) {
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
      standata <- c(standata, nlist(
        Nloc, Jloc, Nneigh, eigenW, Nedges = nrow(edges),  
        edges1 = as.array(edges[, 1]), edges2 = as.array(edges[, 2])
      ))
    } else if (is.cor_fixed(autocor)) {
      V <- autocor$V
      rmd_rows <- attr(data, "na.action")
      if (!is.null(rmd_rows)) {
        V <- V[-rmd_rows, -rmd_rows, drop = FALSE]
      }
      if (nrow(V) != nrow(data)) {
        stop2("'V' must have the same number of rows as 'data'.")
      }
      if (min(eigen(V)$values <= 0)) {
        stop2("'V' must be positive definite.")
      }
      standata$V <- V
      # simplifies code of choose_N
      standata$N_tg <- 1
    } else if (is.cor_bsts(autocor)) {
      if (length(bterms$response) > 1L) {
        stop2("BSTS structure not yet implemented for multivariate models.")
      }
      if (nchar(bterms$time$group)) {
        tgroup <- data[[bterms$time$group]]
      } else {
        tgroup <- rep(1, standata$N) 
      }
      standata$tg <- as.array(as.numeric(factor(tgroup)))
    }
  }
  
  if (old_mv) {
    # deprecated as of brms 1.0.0
    # evaluate even if check_response is FALSE to ensure 
    # that N_trait is defined
    if (is_linear && length(bterms$response) > 1L) {
      standata$Y <- matrix(standata$Y, ncol = length(bterms$response))
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
  
  standata$prior_only <- ifelse(identical(sample_prior, "only"), 1L, 0L)
  if (isTRUE(control$save_order)) {
    attr(standata, "old_order") <- attr(data, "old_order")
  }
  structure(standata, class = "standata")
}  

#' @export
brmdata <- function(...)  {
  # deprecated alias of make_standata
  warn_deprecated("make_standata")
  make_standata(...)
}
