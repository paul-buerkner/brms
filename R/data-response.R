#' Extract response values
#' 
#' Extract response values from a \code{\link{brmsfit}} object.
#' 
#' @param x A \code{\link{brmsfit}} object.
#' @param resp Optional names of response variables for which to extract values.
#' @param warn For internal use only.
#' @param ... Further arguments passed to \code{\link{standata}}.
#' 
#' @return Returns a vector of response values for univariate models and a
#'   matrix of response values with one column per response variable for
#'   multivariate models.
#' 
#' @keywords internal
#' @export
get_y <- function(x, resp = NULL, warn = FALSE, ...) {
  stopifnot(is.brmsfit(x))
  resp <- validate_resp(resp, x)
  warn <- as_one_logical(warn)
  args <- list(x, resp = resp, ...)
  args$re_formula <- NA
  args$check_response <- TRUE
  args$only_response <- TRUE
  args$internal <- TRUE
  sdata <- do_call(standata, args)
  if (warn) {
    if (any(paste0("cens", usc(resp)) %in% names(sdata))) {
      warning2("Results may not be meaningful for censored models.")
    }
  }
  Ynames <- paste0("Y", usc(resp))
  if (length(Ynames) > 1L) {
    out <- do_call(cbind, sdata[Ynames])
    colnames(out) <- resp
  } else {
    out <- sdata[[Ynames]]
  }
  structure(out, old_order = attr(sdata, "old_order"))
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
    c(out) <- data_response(x$terms[[i]], old_sdata = od, ...)
  }
  if (x$rescor) {
    out$nresp <- length(x$responses)
    out$nrescor <- out$nresp * (out$nresp - 1) / 2
  }
  out
}

#' @export
data_response.brmsterms <- function(x, data, check_response = TRUE,
                                    not4stan = FALSE, old_sdata = NULL, ...) {
  # prepare data for the response variable
  data <- subset_data(data, x)
  N <- nrow(data)
  Y <- model.response(model.frame(x$respform, data, na.action = na.pass))
  out <- list(N = N, Y = unname(Y))
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
    if (has_multicol(x$family)) {
      if (!is.matrix(out$Y)) {
        stop2("This model requires a response matrix.")
      }
    }
    if (is_dirichlet(x$family)) {
      if (!is_equal(rowSums(out$Y), rep(1, nrow(out$Y)))) {
        stop2("Response values in dirichlet models must sum to 1.")
      }
    }
    ybounds <- family_info(x$family, "ybounds")
    closed <- family_info(x$family, "closed")
    if (is.finite(ybounds[1])) {
      y_min <- min(out$Y, na.rm = TRUE)
      if (closed[1] && y_min < ybounds[1]) {
        stop2("Family '", family4error, "' requires response greater ",
              "than or equal to ", ybounds[1], ".")
      } else if (!closed[1] && y_min <= ybounds[1]) {
        stop2("Family '", family4error, "' requires response greater ",
              "than ", round(ybounds[1], 2), ".")
      }
    }
    if (is.finite(ybounds[2])) {
      y_max <- max(out$Y, na.rm = TRUE)
      if (closed[2] && y_max > ybounds[2]) {
        stop2("Family '", family4error, "' requires response smaller ",
              "than or equal to ", ybounds[2], ".")
      } else if (!closed[2] && y_max >= ybounds[2]) {
        stop2("Family '", family4error, "' requires response smaller ",
              "than ", round(ybounds[2], 2), ".")
      }
    }
    out$Y <- as.array(out$Y)
  }
  # data for addition arguments of the response
  if (has_trials(x$family) || is.formula(x$adforms$trials)) {
    if (!length(x$adforms$trials)) {
      if (is_multinomial(x$family)) {
        stop2("Specifying 'trials' is required in multinomial models.")
      }
      out$trials <- round(max(out$Y, na.rm = TRUE))
      if (isTRUE(is.finite(out$trials))) {
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
      trials <- eval_rhs(x$adforms$trials)
      out$trials <- eval2(trials$vars$trials, data)
      if (!is.numeric(out$trials)) {
        stop2("Number of trials must be numeric.")
      }
      if (any(!is_wholenumber(out$trials) | out$trials < 1)) {
        stop2("Number of trials must be positive integers.")
      }
    } else {
      stop2("Argument 'trials' is misspecified.")
    }
    if (length(out$trials) == 1L) {
      out$trials <- rep(out$trials, nrow(data))
    }
    if (check_response) {
      if (is_multinomial(x$family)) {
        if (!is_equal(rowSums(out$Y), out$trials)) {
          stop2("Number of trials does not match the number of events.")
        }
      } else if (has_trials(x$family)) {
        if (max(out$trials) == 1L && !not4stan) {
          message("Only 2 levels detected so that family 'bernoulli' ",
                  "might be a more efficient choice.")
        }
        if (any(out$Y > out$trials)) {
          stop2("Number of trials is smaller than the number of events.")
        }
      }
    }
    out$trials <- as.array(out$trials)
  }
  if (has_cat(x$family) || is.formula(x$adforms$cat)) {
    if (!length(x$adforms$cat)) {
      if (!is.null(old_sdata$ncat)) {
        ncat <- old_sdata$ncat
      } else if (has_multicol(x$family)) {
        ncat <- NCOL(out$Y)
      } else {
        ncat <- max(out$Y)
      }
    } else if (is.formula(x$adforms$cat)) {
      cat <- eval_rhs(x$adforms$cat)
      ncat <- eval2(cat$vars$cat, data)
      grcat <- extract_grcat(x, data)
      if (!is.null(grcat)) {
        # TODO: add variable to get_group_vars()
        # TODO: handle new data
        # TODO: expect ncat to be a variable or named vector? #675
        grcat_levels <- levels(grcat)
        out$ngrcat <- length(grcat_levels)
        if (length(ncat) == 1L) {
          ncat <- rep(ncat, length(grcat))
        }
        if (length(ncat) != length(grcat)) {
          stop2("Variables passed to 'resp_cat' need to be of the same length.")
        }
        for (l in grcat_levels) {
          # validate values of the same level
          take <- grcat %in% l
          if (length(unique(ncat[take])) > 1L) {
            stop2(
              "Number of response categories should be ", 
              "unique for each group. Occured for level '", 
              grcat_levels[l], "' of group '", cat$vars$gr, "'."
            )
          }
        }
        # create an matrix of threshold indices per observation
        Jgrcat <- match(grcat, grcat_levels)
        not_dupl_Mcat <- !duplicated(Jgrcat)
        to_order <- order(Jgrcat[not_dupl_Mcat])
        ncat <- ncat[not_dupl_Mcat][to_order]
        Kthres_cumsum <- cumsum(ncat - 1)
        Kthres_start <- c(1, Kthres_cumsum[-length(ncat)] + 1)
        Kthres_end <- Kthres_cumsum
        Jthres <- cbind(Kthres_start, Kthres_end)[Jgrcat, ]
        out$Jthres <- Jthres
      } else {
        ncat <- as_one_numeric(ncat)
        if (!is_wholenumber(ncat) || ncat < 1L) {
          stop2("Number of categories must be a positive integer.")
        }
      }
    } else {
      stop2("Argument 'cat' is misspecified.")
    }
    if (min(ncat) < 2L) {
      stop2("At least two response categories are required.")
    }
    if (!has_multicol(x$family)) {
      if (ncat == 2L && !not4stan) {
        message("Only 2 levels detected so that family 'bernoulli' ",
                "might be a more efficient choice.")
      }
      if (check_response && any(out$Y > ncat)) {
        stop2("Number of categories is smaller than the response ",
              "variable would suggest.")
      }
    }
    out$ncat <- ncat
  }
  if (is.formula(x$adforms$se)) {
    se <- eval_rhs(x$adforms$se)
    out$se <- eval2(se$vars$se, data) 
    if (!is.numeric(out$se)) {
      stop2("Standard errors must be numeric.")
    }
    if (min(out$se) < 0) {
      stop2("Standard errors must be non-negative.")
    }
    out$se <- as.array(out$se)
  }
  if (is.formula(x$adforms$weights)) {
    weights <- eval_rhs(x$adforms$weights)
    out$weights <- eval2(weights$vars$weights, data)  
    if (!is.numeric(out$weights)) {
      stop2("Weights must be numeric.")
    }
    if (min(out$weights) < 0) {
      stop2("Weights must be non-negative.")
    }
    if (weights$flags$scale) {
      out$weights <- out$weights / sum(out$weights) * length(out$weights)
    }
    out$weights <- as.array(out$weights)
  }
  if (is.formula(x$adforms$dec)) {
    dec <- eval_rhs(x$adforms$dec)
    out$dec <- eval2(dec$vars$dec, data)
    if (is.character(out$dec) || is.factor(out$dec)) {
      if (!all(unique(out$dec) %in% c("lower", "upper"))) {
        stop2("Decisions should be 'lower' or 'upper' ",
              "when supplied as characters or factors.")
      }
      out$dec <- ifelse(out$dec == "lower", 0, 1)
    } else {
      out$dec <- as.numeric(as.logical(out$dec))
    }
    out$dec <- as.array(out$dec)
  }
  if (is.formula(x$adforms$rate)) {
    rate <- eval_rhs(x$adforms$rate)
    out$denom <- eval2(rate$vars$denom, data)
    if (!is.numeric(out$denom)) {
      stop2("Rate denomiators should be numeric.")
    }
    if (isTRUE(any(out$denom <= 0))) {
      stop2("Rate denomiators should be positive.")
    }
    out$denom <- as.array(out$denom)
  }
  if (is.formula(x$adforms$cens) && check_response) {
    cens <- eval_rhs(x$adforms$cens)
    out$cens <- eval2(cens$vars$cens, data)
    out$cens <- as.array(prepare_cens(out$cens))
    if (!all(is_wholenumber(out$cens) & out$cens %in% -1:2)) {
      stop2(
        "Invalid censoring data. Accepted values are ",
        "'left', 'none', 'right', and 'interval'\n",
        "(abbreviations are allowed) or -1, 0, 1, and 2.\n",
        "TRUE and FALSE are also accepted ",
        "and refer to 'right' and 'none' respectively."
      )
    }
    icens <- out$cens %in% 2
    if (any(icens)) {
      if (cens$vars$y2 == "NA") {
        stop2("Argument 'y2' is required for interval censored data.")
      }
      y2 <- unname(eval2(cens$vars$y2, data))
      if (any(out$Y[icens] >= y2[icens])) {
        stop2("Left censor points must be smaller than right ",
              "censor points for interval censored data.")
      }
      y2[!icens] <- 0  # not used in Stan
      out$rcens <- as.array(y2)
    }
  }
  if (is.formula(x$adforms$trunc)) {
    trunc <- eval_rhs(x$adforms$trunc)
    out$lb <- as.numeric(eval2(trunc$vars$lb, data))
    out$ub <- as.numeric(eval2(trunc$vars$ub, data))
    if (any(out$lb >= out$ub)) {
      stop2("Truncation bounds are invalid: lb >= ub")
    }
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
      if (!not4stan) {
        out$noise[which_mi] <- Inf
      }
    }
    if (!not4stan) {
      # Stan does not allow NAs in data
      # use Inf to that min(Y) is not affected
      out$Y[which_mi] <- Inf
    }
  }
  if (is.formula(x$adforms$vreal)) {
    # vectors of real values for use in custom families
    vreal <- eval_rhs(x$adforms$vreal)
    vreal <- lapply(vreal$vars, eval2, data)
    names(vreal) <- paste0("vreal", seq_along(vreal))
    for (i in seq_along(vreal)) {
      if (length(vreal[[i]]) == 1L) {
        vreal[[i]] <- rep(vreal[[i]], N)
      }
      vreal[[i]] <- as.array(as.numeric(vreal[[i]]))
    }
    c(out) <- vreal
  }
  if (is.formula(x$adforms$vint)) {
    # vectors of integer values for use in custom families
    vint <- eval_rhs(x$adforms$vint)
    vint <- lapply(vint$vars, eval2, data)
    names(vint) <- paste0("vint", seq_along(vint))
    for (i in seq_along(vint)) {
      if (length(vint[[i]]) == 1L) {
        vint[[i]] <- rep(vint[[i]], N)
      }
      if (!all(is_wholenumber(vint[[i]]))) {
        stop2("'vint' requires whole numbers as input.")
      }
      vint[[i]] <- as.array(vint[[i]])
    }
    c(out) <- vint
  }
  if (length(out)) {
    resp <- usc(combine_prefix(x))
    out <- setNames(out, paste0(names(out), resp))
  }
  out
}

# data specific for mixture models
data_mixture <- function(bterms, prior = brmsprior()) {
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
        theta_prior <- eval_dirichlet(theta_prior)
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

# data for the baseline functions of Cox models
data_bhaz <- function(bterms, data, basis = NULL) {
  out <- list()
  if (!is_cox(bterms$family)) {
    return(out) 
  }
  y <- model.response(model.frame(bterms$respform, data, na.action = na.pass))
  args <- bterms$family$bhaz 
  out$Zbhaz <- bhaz_basis_matrix(y, args, basis = basis)
  out$Zcbhaz <- bhaz_basis_matrix(y, args, integrate = TRUE, basis = basis)
  out$Kbhaz <- NCOL(out$Zbhaz)
  out
}

# Basis matrices for baseline hazard functions of the Cox model
# @param y vector of response values
# @param args arguments passed to the spline generating functions
# @param integrate compute the I-spline instead of the M-spline basis?
# @param basis optional precomputed basis matrix
# @return the design matrix of the baseline hazard function
bhaz_basis_matrix <- function(y, args = list(), integrate = FALSE, 
                              basis = NULL) {
  require_package("splines2")
  if (!is.null(basis)) {
    # perform predictions based on an existing basis matrix
    stopifnot(inherits(basis, "mSpline"))
    if (integrate) {
      # for predictions just the attibutes are required
      # which are the same of M-Splines and I-Splines
      class(basis) <- c("matrix", "iSpline")
    }
    return(predict(basis, y))
  }
  stopifnot(is.list(args))
  args$x <- y
  if (!is.null(args$intercept)) {
    args$intercept <- as_one_logical(args$intercept) 
  }
  if (is.null(args$Boundary.knots)) {
    if (isTRUE(args$intercept)) {
      lower_knot <- min(y)
      upper_knot <- max(y)
    } else {
      # we need a smaller lower boundary knot to avoid lp = -Inf 
      # the below choices are ad-hoc and may need further thought
      lower_knot <- max(min(y) - mad(y, na.rm = TRUE) / 10, 0)
      upper_knot <- max(y) + mad(y, na.rm = TRUE) / 10
    }
    args$Boundary.knots <- c(lower_knot, upper_knot)
  }
  if (integrate) {
    out <- do_call(splines2::iSpline, args)
  } else {
    out <- do_call(splines2::mSpline, args)
  }
  out
}

# extract names of response categories
extract_cat_names <- function(x, data) {
  stopifnot(is.brmsformula(x) || is.brmsterms(x))
  respform <- validate_resp_formula(x$formula)
  mr <- model.response(model.frame(respform, data))
  if (is_ordinal(x) && is.numeric(mr)) {
    out <- as.character(seq_len(max(mr)))
  } else if (has_multicol(x)) {
    mr <- as.matrix(mr)
    out <- as.character(colnames(mr))
    if (!length(out)) {
      out <- as.character(seq_cols(mr))
    }
  } else {
    out <- levels(factor(mr))
  }
  out
}

# coerce censored values into the right format
# @param x vector of censoring indicators
# @return transformed vector of censoring indicators
prepare_cens <- function(x) {
  .prepare_cens <- function(x) {  
    stopifnot(length(x) == 1L)
    regx <- paste0("^", x)
    if (grepl(regx, "left")) {
      x <- -1
    } else if (grepl(regx, "none") || isFALSE(x)) {
      x <- 0
    } else if (grepl(regx, "right") || isTRUE(x)) {
      x <- 1
    } else if (grepl(regx, "interval")) {
      x <- 2
    }
    return(x)
  }
  x <- unname(x)
  if (is.factor(x)) {
    x <- as.character(x)
  }
  ulapply(x, .prepare_cens)
}

# extract information on censoring of the response variable
# @param x a brmsfit object
# @param resp optional names of response variables for which to extract values
# @return vector of censoring indicators or NULL in case of no censoring
get_cens <- function(x, resp = NULL, newdata = NULL) {
  stopifnot(is.brmsfit(x))
  resp <- validate_resp(resp, x, multiple = FALSE)
  bterms <- parse_bf(x$formula)
  if (!is.null(resp)) {
    bterms <- bterms$terms[[resp]]
  }
  if (is.null(newdata)) {
    newdata <- model.frame(x)
  }
  out <- NULL
  if (is.formula(bterms$adforms$cens)) {
    cens <- eval_rhs(bterms$adforms$cens)
    out <- eval2(cens$vars$cens, newdata)
    out <- prepare_cens(out)
  }
  out
}

# extract truncation boundaries
trunc_bounds <- function(x, ...) {
  UseMethod("trunc_bounds")
}

# @return a named list with one element per response variable
#' @export
trunc_bounds.mvbrmsterms <- function(x, ...) {
  lapply(x$terms, trunc_bounds, ...)
}

# @param data data.frame containing the truncation variables
# @param incl_family include the family in the derivation of the bounds?
# @param stan return bounds in form of Stan syntax?
# @return a list with elements 'lb' and 'ub'
#' @export
trunc_bounds.brmsterms <- function(x, data = NULL, incl_family = FALSE, 
                                   stan = FALSE, ...) {
  if (is.formula(x$adforms$trunc)) {
    trunc <- eval_rhs(x$adforms$trunc)
  } else {
    trunc <- resp_trunc()
  }
  out <- list(
    lb = eval2(trunc$vars$lb, data),
    ub = eval2(trunc$vars$ub, data)
  )
  if (incl_family) {
    family_bounds <- family_bounds(x)
    out$lb <- max(out$lb, family_bounds$lb)
    out$ub <- min(out$ub, family_bounds$ub)
  }
  if (stan) {
    if (any(out$lb > -Inf | out$ub < Inf)) {
      tmp <- c(
        if (out$lb > -Inf) paste0("lower=", out$lb),
        if (out$ub < Inf) paste0("upper=", out$ub)
      )
      out <- paste0("<", paste0(tmp, collapse = ","), ">")
    } else {
      out <- ""
    }
  }
  out
}

# check if the model has group specific ordinal thresholds
# @param x list with potentail $adforms elements
# @return TRUE or FALSE
has_grcat <- function(x) {
  cat <- x$adforms$cat
  isTRUE(!is.null(cat) && eval_rhs(cat)$vars$gr != "NA")
}

# extract variabe for group specific ordinal thresholds
# @param x list with potentail $adforms elements
# @param data data passed by the user
# @return a vector of threshold groups or NULL
extract_grcat <- function(x, data) {
  if (!has_grcat(x)) {
    return(NULL)
  }
  cat <- eval_rhs(x$adforms$cat)
  factor(eval2(cat$vars$gr, data))
}

# check if addition argument 'subset' ist used in the model
has_subset <- function(bterms) {
  .has_subset <- function(x) {
    is.formula(x$adforms$subset)
  }
  if (is.brmsterms(bterms)) {
    out <- .has_subset(bterms)
  } else if (is.mvbrmsterms(bterms)) {
    out <- any(ulapply(bterms$terms, .has_subset))
  } else {
    out <- FALSE
  }
  out 
}
