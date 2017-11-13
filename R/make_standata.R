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
                          prior = NULL, autocor = NULL, 
                          nonlinear = NULL, cov_ranef = NULL, 
                          sample_prior = c("no", "yes", "only"), 
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
  formula <- amend_formula(
    formula, data = data, family = family, 
    autocor = autocor, nonlinear = nonlinear
  )
  family <- formula$family
  autocor <- formula$autocor
  bterms <- parse_bf(formula)
  old_mv <- isTRUE(formula[["old_mv"]])
  is_linear <- is_linear(family)
  is_ordinal <- is_ordinal(family)
  is_count <- is_count(family)
  is_forked <- is_forked(family)
  is_categorical <- is_categorical(family)
  sample_prior <- check_sample_prior(sample_prior)
  check_prior_content(prior, warn = FALSE)
  prior <- check_prior_special(bterms, prior = prior)
  na_action <- if (is_newdata) na.pass else na.omit
  data <- update_data(
    data, bterms = bterms, na.action = na_action, 
    drop.unused.levels = !is_newdata, knots = knots,
    terms_attr = control$terms_attr
  )
  if (has_arma(autocor) || is.cor_bsts(autocor)) {
    # order data in case of autocorrelation models
    data <- order_data(data, bterms = bterms, old_mv = old_mv)
  }
  
  # response variable
  out <- list(N = nrow(data), Y = unname(model.response(data)))
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
    if (!factors_allowed && !is.numeric(out$Y)) {
      stop2("Family '", family4error, "' requires numeric responses.")
    }
    # transform and check response variables for different families
    regex_pos_int <- "(^|_)(binomial|poisson|negbinomial|geometric)$"
    if (any(grepl(regex_pos_int, families))) {
      if (!all(is_wholenumber(out$Y)) || min(out$Y) < 0) {
        stop2("Family '", family4error, "' requires responses ", 
              "to be non-negative integers.")
      }
    } else if (any(families %in% "bernoulli")) {
      out$Y <- as.numeric(as.factor(out$Y)) - 1
      if (any(!out$Y %in% c(0, 1))) {
        stop2("Family '", family4error, "' requires responses ", 
              "to contain only two different values.")
      }
    } else if (any(grepl("(^|_)beta$", families))) {
      if (any(families %in% "beta")) {
        lower <- any(out$Y <= 0)
      } else {
        lower <- any(out$Y < 0) 
      } 
      if (any(families %in% "zero_one_inflated_beta")) {
        upper <- any(out$Y > 1) 
      } else {
        upper <- any(out$Y >= 1) 
      }
      if (lower || upper) {
        stop2("Family '", family4error, "' requires responses ", 
              "between 0 and 1.")
      }
    } else if (any(families %in% "von_mises")) {
      if (any(out$Y < -pi | out$Y > pi)) {
        stop2("Family '", family4error, "' requires responses ",
              "between -pi and pi.")
      }
    } else if (is_categorical) { 
      out$Y <- as.numeric(factor(out$Y))
      if (length(unique(out$Y)) < 3L) {
        stop2("At least three response categories are required.")
      }
    } else if (is_ordinal) {
      if (is.ordered(out$Y)) {
        out$Y <- as.numeric(out$Y)
      } else if (all(is_wholenumber(out$Y))) {
        out$Y <- out$Y - min(out$Y) + 1
      } else {
        stop2("Family '", family4error, "' requires either integers or ",
              "ordered factors as responses.")
      }
      if (length(unique(out$Y)) < 2L) {
        stop2("At least two response categories are required.")
      }
    } else if (is_skewed(family) || is_lognormal(family) || is_wiener(family)) {
      if (min(out$Y) <= 0) {
        stop2("Family '", family4error, "' requires responses ", 
              "to be positive.")
      }
    } else if (is_zero_inflated(family) || is_hurdle(family)) {
      if (min(out$Y) < 0) {
        stop2("Family '", family4error, "' requires responses ", 
              "to be non-negative.")
      }
    }
    out$Y <- as.array(out$Y)
  }
  
  # data for various kinds of effects
  only_response <- isTRUE(control$only_response)
  if (!only_response) {
    ranef <- tidy_ranef(
      bterms, data, ncat = control$ncat, 
      old_levels = control$old_levels
    )
    args_eff <- nlist(data, ranef, prior, knots, not4stan)
    resp <- bterms$response
    if (length(resp) > 1L && !old_mv) {
      args_eff_spec <- list(
        x = bterms$dpars[["mu"]], smooths = control$smooths[["mu"]],
        gps = control$gps[["mu"]], Jmo = control$Jmo[["mu"]]
      )
      bterms$dpars[["mu"]] <- NULL
      for (r in resp) {
        args_eff_spec$x$resp <- r
        data_eff <- do.call(data_effects, c(args_eff_spec, args_eff))
        out <- c(out, data_eff)
      }
      if (is_linear(family)) {
        out$nresp <- length(resp)
        out$nrescor <- length(resp) * (length(resp) - 1) / 2
        colnames(out$Y) <- resp
      }
    }
    # data for predictors of distributional parameters
    for (dp in names(bterms$dpars)) {
      args_eff_spec <- list(
        x = bterms$dpars[[dp]], 
        smooths = control$smooths[[dp]],
        gps = control$gps[[dp]], 
        Jmo = control$Jmo[[dp]]
      )
      data_aux_eff <- do.call(data_effects, c(args_eff_spec, args_eff))
      out <- c(out, data_aux_eff)
    }
    for (dp in names(bterms$fdpars)) {
      out[[dp]] <- bterms$fdpars[[dp]]$value
    }
    out <- c(out,
      data_gr(ranef, data, cov_ranef = cov_ranef),
      data_Xme(bterms, data),
      data_mixture(bterms, prior = prior)
    )
  }
  
  # data for specific families
  if (has_trials(family)) {
    if (!length(bterms$adforms$trials)) {
      if (!is.null(control$trials)) {
        out$trials <- control$trials
      } else {
        message("Using the maximum of the response ", 
                "variable as the number of trials.")
        out$trials <- max(out$Y) 
      }
    } else if (is.formula(bterms$adforms$trials)) {
      out$trials <- eval_rhs(bterms$adforms$trials, data = data)
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
  if (has_cat(family)) {
    if (!length(bterms$adforms$cat)) {
      if (!is.null(control$ncat)) {
        out$ncat <- control$ncat
      } else {
        out$ncat <- length(unique(out$Y))
      }
    } else if (is.formula(bterms$adforms$cat)) { 
      out$ncat <- eval_rhs(bterms$adforms$cat, data = data)
    } else {
      stop2("Argument 'cat' is misspecified.")
    }
    if (max(out$ncat) == 2L) {
      message("Only 2 levels detected so that family 'bernoulli' ",
              "might be a more efficient choice.")
    }
    if (check_response && any(out$Y > out$ncat)) {
      stop2("Number of categories is smaller than the response ", 
            "variable would suggest.")
    }
  }
  
  # data for addition arguments
  if (is.formula(bterms$adforms$se)) {
    out[["se"]] <- as.array(eval_rhs(bterms$adforms$se, data = data))
  }
  if (is.formula(bterms$adforms$weights)) {
    out[["weights"]] <- 
      as.array(eval_rhs(bterms$adforms$weights, data = data))
    if (old_mv) {
      out$weights <- out$weights[1:out$N_trait]
    }
  }
  if (is.formula(bterms$adforms$disp)) {
    out[["disp"]] <- as.array(eval_rhs(bterms$adforms$disp, data = data))
  }
  if (is.formula(bterms$adforms$dec)) {
    out[["dec"]] <- as.array(eval_rhs(bterms$adforms$dec, data = data))
  }
  if (is.formula(bterms$adforms$cens) && check_response) {
    cens <- eval_rhs(bterms$adforms$cens, data = data)
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
    if (old_mv) {
      out$cens <- out$cens[1:out$N_trait]
    }
  }
  if (is.formula(bterms$adforms$trunc)) {
    out <- c(out, eval_rhs(bterms$adforms$trunc, data = data))
    if (length(out$lb) == 1L) {
      out$lb <- rep(out$lb, out$N)
    }
    if (length(out$ub) == 1L) {
      out$ub <- rep(out$ub, out$N)
    }
    if (length(out$lb) != out$N || 
        length(out$ub) != out$N) {
      stop2("Invalid truncation bounds.")
    }
    inv_bounds <- out$Y < out$lb | out$Y > out$ub
    if (check_response && any(inv_bounds)) {
      stop2("Some responses are outside of the truncation bounds.")
    }
  }
  
  # autocorrelation variables
  if (!only_response) {
    if (is.cor_arma(autocor) || is.cor_bsts(autocor)) {
      if (nchar(bterms$time$group)) {
        tgroup <- as.numeric(factor(data[[bterms$time$group]]))
      } else {
        tgroup <- rep(1, out$N) 
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
          out$begin_tg <- as.array(
            ulapply(unique(tgroup), match, tgroup)
          )
          out$nobs_tg <- as.array(with(out, 
            c(if (N_tg > 1L) begin_tg[2:N_tg], N + 1) - begin_tg
          ))
          out$end_tg <- with(out, begin_tg + nobs_tg - 1)
        } else {
          # data for the 'predictor' version of ARMA
          max_lag <- max(Kar, Kma)
          out$J_lag <- rep(0, out$N)
          for (n in seq_len(out$N)) {
            for (i in seq_len(max_lag)) {
              valid_lag <- n + 1 - i > 0 && n < out$N && 
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
        out$Yarr <- arr_design_matrix(out$Y, Karr, tgroup)
        out$Karr <- Karr
      }
    } else if (is.cor_sar(autocor)) {
      if (!identical(dim(autocor$W), rep(out$N, 2))) {
        stop2("Dimensions of 'W' must be equal to the number of observations.")
      }
      out$W <- autocor$W
      # simplifies code of choose_N
      out$N_tg <- 1
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
        Nloc <- out$N
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
      if (nrow(V) != nrow(data)) {
        stop2("'V' must have the same number of rows as 'data'.")
      }
      if (min(eigen(V)$values <= 0)) {
        stop2("'V' must be positive definite.")
      }
      out$V <- V
      # simplifies code of choose_N
      out$N_tg <- 1
    }
  }
  
  if (old_mv) {
    # deprecated as of brms 1.0.0
    # evaluate even if check_response is FALSE to ensure 
    # that N_trait is defined
    if (is_linear && length(bterms$response) > 1L) {
      out$Y <- matrix(out$Y, ncol = length(bterms$response))
      NC_trait <- ncol(out$Y) * (ncol(out$Y) - 1L) / 2L
      out <- c(out, list(N_trait = nrow(out$Y), 
                         K_trait = ncol(out$Y),
                         NC_trait = NC_trait)) 
      # for compatibility with the S3 methods of brms >= 1.0.0
      out$nresp <- out$K_trait
      out$nrescor <- out$NC_trait
    }
    if (is_forked) {
      # the second half of Y is only dummy data
      # that was put into data to make melt_data work correctly
      out$N_trait <- nrow(data) / 2L
      out$Y <- as.array(out$Y[1L:out$N_trait]) 
    }
    if (is_categorical && !isTRUE(control$old_cat == 1L)) {
      ncat1m <- out$ncat - 1L
      out$N_trait <- nrow(data) / ncat1m
      out$Y <- as.array(out$Y[1L:out$N_trait])
      out$J_trait <- as.array(matrix(1L:out$N, ncol = ncat1m))
    }
  }
  
  out$prior_only <- as.integer(identical(sample_prior, "only"))
  if (isTRUE(control$save_order)) {
    attr(out, "old_order") <- attr(data, "old_order")
  }
  structure(out, class = "standata")
}  

#' @export
brmdata <- function(...)  {
  # deprecated alias of make_standata
  warn_deprecated("make_standata")
  make_standata(...)
}
