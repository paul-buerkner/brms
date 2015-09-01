# melt data frame for multinormal models
#
# @param data a data.frame
# @param response names of the response variables
# @param family
#
# @return data in long format 
melt <- function(data, response, family) {
  if (length(response) > 1 && family != "gaussian")
    stop("multivariate models are currently only allowed for family 'gaussian'")
  else if (length(response) > 1 && family == "gaussian") {
    if (!is(data, "data.frame"))
      stop("data must be a data.frame in case of multiple responses")
    if ("trait" %in% names(data))
      stop("trait is a resevered variable name in case of multiple responses")
    new_columns <- structure(data.frame(unlist(lapply(response, rep, time = nrow(data))), 
                                        as.numeric(as.matrix(data[,response]))),
                             names = c("trait", response[1]))
    old_columns <- data[,which(!names(data) %in% response), drop = FALSE]
    old_columns <- do.call(rbind, lapply(response, function(i) old_columns))
    data <- cbind(old_columns, new_columns)
  }
  data
}  

# combine grouping factors
#
# @param data a data.frame
# @param ... the grouping factors to be combined. 
#
# @return a data.frame containing all old variables and the new combined grouping factors
combine_groups <- function(data, ...) {
  group <- c(...)
  if (length(group)) {
    for (i in 1:length(group)) {
      sgroup <- unlist(strsplit(group[[i]], "__"))
      if (length(sgroup) > 1) {
        new.var <- get(sgroup[1], data)
        for (j in 2:length(sgroup)) {
          new.var <- paste0(new.var, "_", get(sgroup[j], data))
        }
        data[[group[[i]]]] <- new.var
      }
    } 
  }
  data
}

# update data for use in brm
#
# @param data the original data.frame
# @param family
# @param effects output of extract_effects
#
# @return model.frame in long format with combined grouping variables if present
update_data <- function(data, family, effects, ...) {
  if (!"brms.frame" %in% class(data)) {
    data <- melt(data, response = effects$response, family = family)
    data <- stats::model.frame(effects$all, data = data, drop.unused.levels = TRUE)
    if (any(grepl("__", colnames(data))))
      stop("Variable names may not contain double underscores '__'")
    data <- combine_groups(data, effects$group, ...)
    class(data) <- c("brms.frame", "data.frame") 
  }
  data
}

#' Extract required data for \code{brms} models
#'
#' @inheritParams brm
#' 
#' @aliases brm.data
#' 
#' @return A named list of objects containing the required data to fit a \code{brms} model 
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' data1 <- brmdata(rating ~ treat + period + carry + (1|subject), 
#'                   data = inhaler, family = "cumulative")
#' names(data1)
#' 
#' data2 <- brmdata(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|visit), 
#'                   data = epilepsy, family = "poisson")
#' names(data2)
#'          
#' @export
brmdata <- function(formula, data = NULL, family = "gaussian", autocor = NULL, 
                    partial = NULL, cov.ranef = NULL) {
  family <- check_family(family[1])
  is_linear <- family %in% c("gaussian", "student", "cauchy")
  is_ordinal <- family %in% c("cumulative","cratio","sratio","acat")
  is_count <- family %in% c("poisson", "negbinomial", "geometric")
  is_skew <- family %in% c("gamma", "weibull", "exponential")
  if (is.null(autocor)) autocor <- cor_arma()
  if (!is(autocor,"cor_brms")) stop("cor must be of class cor_brms")
  
  et <- extract_time(autocor$formula)
  ee <- extract_effects(formula = formula, family = family, partial, et$all)
  data <- update_data(data, family = family, effects = ee, et$group)
  
  #sort data in case of autocorrelation models
  if (sum(autocor$p, autocor$q) > 0) {
    if (family == "gaussian" && length(ee$response) > 1) {
      if (!grepl("^trait$|__trait$|^trait__|__trait__", et$group))
        stop("autocorrelation structures for multiple responses must contain 'trait' as grouping variable")
      else to.order <- rmNULL(list(data[["trait"]], data[[et$group]], data[[et$time]]))
    }
    else to.order <- rmNULL(list(data[[et$group]], data[[et$time]]))
    if (length(to.order)) data <- data[do.call(order, to.order),]
  }
  
  #response variable
  standata <- list(N = nrow(data), Y = unname(model.response(data)))
  if (!is.numeric(standata$Y) && !(is_ordinal || family %in% c("bernoulli", "categorical"))) 
    stop(paste("family", family, "expects numeric response variable"))
  
  #transform and check response variable for different families
  if (is_count || family == "binomial") {
    if (!all(is.wholenumber(standata$Y)) || min(standata$Y) < 0)
      stop(paste("family", family, "expects response variable of non-negative integers"))
  }
  else if (family == "bernoulli") {
    standata$Y <- as.numeric(as.factor(standata$Y)) - 1
    if (any(!standata$Y %in% c(0,1)))
      stop("family bernoulli expects response variable to contain only two different values")
  }
  else if (family == "categorical") 
    standata$Y <- as.numeric(as.factor(standata$Y))
  else if (is_ordinal) {
    if (is.factor(standata$Y)) {
      if (is.ordered(standata$Y)) standata$Y <- as.numeric(standata$Y)
      else stop(paste("family", family, "requires factored response variables to be ordered"))
    }
    else if (all(is.wholenumber(standata$Y)))
      standata$Y <- standata$Y - min(standata$Y) + 1
    else stop(paste("family", family, "expects either integers or ordered factors as response variables"))
  }
  else if (is_skew) {
    if (min(standata$Y) < 0)
      stop(paste("family", family, "requires response variable to be non-negative"))
  }  
  else if (family == "gaussian" && length(ee$response) > 1) {
    standata$Y <- matrix(standata$Y, ncol = length(ee$response))
    standata <- c(standata, list(N_trait = nrow(standata$Y), K_trait = ncol(standata$Y)),
                   NC_trait = ncol(standata$Y) * (ncol(standata$Y)-1)/2) 
  }
  
  #fixed effects data
  X <- get_model_matrix(ee$fixed, data, rm.int = is_ordinal)
  if (family == "categorical") standata <- c(standata, list(Kp = ncol(X), Xp = X))
  else standata <- c(standata, list(K = ncol(X), X = X))
  
  #random effects data
  if (length(ee$random)) {
    Z <- lapply(ee$random, get_model_matrix, data = data)
    r <- lapply(Z, colnames)
    ncolZ <- lapply(Z, ncol)
    expr <- expression(as.numeric(as.factor(get(g, data))), length(unique(get(g, data))), 
                       ncolZ[[i]], Z[[i]], ncolZ[[i]]*(ncolZ[[i]]-1)/2)
    for (i in 1:length(ee$group)) {
      g <- ee$group[[i]]
      name <- paste0(c("lev_", "N_", "K_", "Z_", "NC_"), i)
      if (ncolZ[[i]] == 1) Z[[i]] <- as.vector(Z[[i]])
      for (j in 1:length(name)) standata <- c(standata, setNames(list(eval(expr[j])), name[j]))
      if (g %in% names(cov.ranef)) {
        cov_mat <- as.matrix(cov.ranef[[g]])
        found_level_names <- rownames(cov_mat)
        colnames(cov_mat) <- found_level_names
        true_level_names <- sort(as.character(unique(data[[g]])))
        if (is.null(found_level_names)) 
          stop(paste("Row names are required for covariance matrix of",g))
        if (nrow(cov_mat) != length(true_level_names))
          stop(paste("Dimension of covariance matrix of",g,"is incorrect"))
        if (any(sort(found_level_names) != true_level_names))
          stop(paste("Row names of covariance matrix of",g,"do not match names of the grouping levels"))
        if (!isSymmetric(unname(cov_mat)))
          stop(paste("Covariance matrix of grouping factor",g,"is not symmetric"))
        if (min(eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values) <= 0)
          warning(paste("Covariance matrix of grouping factor",g,"may not be positive definite"))
        cov_mat <- cov_mat[order(found_level_names), order(found_level_names)]
        if (length(r[[i]]) == 1) {
          cov_mat <- suppressWarnings(chol(cov_mat, pivot = TRUE))
          cov_mat <- t(cov_mat[, order(attr(cov_mat, "pivot"))])
        }  
        else if (length(r[[i]]) > 1 && !ee$cor[[i]])
          cov_mat <- t(suppressWarnings(chol(kronecker(cov_mat, diag(ncolZ[[i]])), pivot = TRUE)))
        standata <- c(standata, setNames(list(cov_mat), paste0("cov_",i)))
      }
    }
  }
  
  #addition and partial variables
  if (is.formula(ee$se)) {
    standata <- c(standata, list(sigma = .addition(formula = ee$se, data = data)))
  }
  if (is.formula(ee$weights)) {
    standata <- c(standata, list(weights = .addition(formula = ee$weights, data = data)))
    if (family == "gaussian" && length(ee$response) > 1) 
      standata$weights <- standata$weights[1:standata$N_trait]
  }
  if (is.formula(ee$cens)) {
    standata <- c(standata, list(cens = .addition(formula = ee$cens, data = data)))
  }
  if (family == "binomial") {
    standata$max_obs <- if (!length(ee$trials)) max(standata$Y)
                        else if (is.wholenumber(ee$trials)) ee$trials
                        else if (is.formula(ee$trials)) .addition(formula = ee$trials, data = data)
                        else stop("Response part of formula is invalid.")
    if (max(standata$max_obs) == 1) 
      message("Only 2 levels detected so that family 'bernoulli' might be a more efficient choice.")
    if (any(standata$Y > standata$max_obs))
      stop("Number of trials is smaller than the response variable would suggest.")
  }
  if (is_ordinal || family == "categorical") {
    standata$max_obs <- if (!length(ee$cat)) max(standata$Y)
                        else if (is.wholenumber(ee$cat)) ee$cat
                        else if (is.formula(ee$cat)) max(.addition(formula = ee$cat, data = data))
                        else stop("Response part of formula is invalid.")
    if (max(standata$max_obs) == 2) 
      message("Only 2 levels detected so that family 'bernoulli' might be a more efficient choice.")
    if (any(standata$Y > standata$max_obs))
      stop("Number of categories is smaller than the response variable would suggest.")
  }  
  
  #get data for partial effects
  if (is.formula(partial)) {
    if (family %in% c("sratio","cratio","acat")) {
      Xp <- get_model_matrix(partial, data, rm.int = TRUE)
      standata <- c(standata, list(Kp = ncol(Xp), Xp = Xp))
      fp <- intersect(colnames(X), colnames(Xp))
      if (length(fp))
        stop(paste("Variables cannot be modeled as fixed and partial effects at the same time.",
                   "Error occured for variables:", paste(fp, collapse = ", ")))
    } 
    else stop("partial is only meaningful for families 'sratio', 'cratio', and 'acat'")  
  }
  
  #autocorrelation variables
  if (is(autocor,"cor_arma") && autocor$p + autocor$q > 0) {
    tgroup <- data[[et$group]]
    if (is.null(tgroup)) tgroup <- rep(1, standata$N) 
    if (autocor$p > 0) {
      standata$Yar <- ar_design_matrix(Y = standata$Y, p = autocor$p, group = tgroup)
      standata$Kar <- autocor$p
    }
    if (autocor$q > 0 && is(autocor,"cor_arma")) {
      standata$Ema_pre <- matrix(0, nrow = standata$N, ncol = autocor$q)
      standata$Kma <- autocor$q
      standata$tgroup <- as.numeric(as.factor(tgroup))
    }
  } 
  standata
}  

# deprectated version of brmdata
#' @export
brm.data <- function(formula, data = NULL, family = "gaussian", autocor = NULL, 
                     partial = NULL, cov.ranef = NULL) 
  brmdata(formula = formula, data = data, family = family, autocor = autocor,
          partial = partial, cov.ranef = cov.ranef)

# Construct Design Matrices for \code{brms} models
# 
# @param formula An object of class "formula"
# @param data A data frame created with \code{model.frame}. If another sort of object, \code{model.frame} is called first.
# @param rm.int Flag indicating if the intercept column should be removed from the model.matrix. 
#   Primarily useful for ordinal models.
# 
# @return The design matrix for a regression-like model with the specified formula and data. 
#   For details see the documentation of \code{model.matrix}.
get_model_matrix <- function(formula, data = environment(formula), rm.int = FALSE) {
  if (!is(formula, "formula")) return(NULL) 
  X <- stats::model.matrix(formula, data)
  new_colnames <- rename(colnames(X), check_dup = TRUE)
  if (rm.int && "Intercept" %in% new_colnames) {
    X <- as.matrix(X[, -(1)])
    if (ncol(X)) colnames(X) <- new_colnames[2:length(new_colnames)]
  } 
  else colnames(X) <- new_colnames
  X   
}

#calculate design matrix for autoregressive effects
#
# @param Y a vector containing the response variable
# @param p autocor$p
# @param group vector containing the grouping variable for each observation
ar_design_matrix <- function(Y, p, group)  { 
  if (length(Y) != length(group)) 
    stop("Y and group must have the same length")
  if (p > 0) {
    U_group <- unique(group)
    N_group <- length(U_group)
    out <- matrix(0, nrow = length(Y), ncol = p)
    ptsum <- rep(0, N_group + 1)
    meanY <- mean(Y) 
    for (j in 1:N_group) {
      ptsum[j+1] <- ptsum[j] + sum(group == U_group[j])
      for (i in 1:p) {
        if (ptsum[j]+i+1 <= ptsum[j+1])
          out[(ptsum[j]+i+1):ptsum[j+1], i] <- Y[(ptsum[j]+1):(ptsum[j+1]-i)] - meanY
      }
    }
  }
  else out <- NULL
  out
}

# amend new data for predict method
# 
# @param new_data a data.frame containing new data for prediction 
# @param fit an object of class \code{brmsfit}
#
# @return updated data.frame to be compatible with fit$formula
amend_new_data <- function(new_data, fit) {
  ee <- extract_effects(fit$formula, family = fit$family)
  if (length(ee$group))
    stop("random effects models not yet supported for predicting new data")
  if (sum(fit$autocor$p, fit$autocor$q) > 0 && !all(ee$response %in% names(new_data))) 
    stop("response variables must be specified in new_data for autocorrelative models")
  else for (resp in ee$response) new_data[[resp]] <- 0 # add irrelevant response variables
  if (is.formula(ee$cens)) 
    for (cens in all.vars(ee$cens)) new_data[[cens]] <- 0 # add irrelevant censor variables
  data <- brmdata(fit$formula, data = new_data, family = fit$family,
                  autocor = fit$autocor, partial = fit$partial)
}

#computes data for addition arguments
.addition <- function(formula, data) {
  if (!is.formula(formula))
    formula <- as.formula(formula)
  eval(formula[[2]], data, environment(formula))
}

#standard errors for meta-analysis
.se <- function(x) {
  if (min(x) < 0) stop("standard errors must be non-negative")
  x  
}

#weights to be applied on any model
.weights <- function(x) {
  if (min(x) < 0) stop("weights must be non-negative")
  x
}

#trials for binomial models
.trials <- function(x) {
  if (any(!is.wholenumber(x) || x < 1))
    stop("number of trials must be positive integers")
  x
}

#number of categories for categorical and ordinal models
.cat <- function(x) {
  if (any(!is.wholenumber(x) || x < 1))
    stop("number of categories must be positive integers")
  x
}

#indicator for censoring
.cens <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  cens <- unname(sapply(x, function(x) {
    if (grepl(paste0("^",x), "right") || is.logical(x) && isTRUE(x)) x <- 1
    else if (grepl(paste0("^",x), "none") || is.logical(x) && !isTRUE(x)) x <- 0
    else if (grepl(paste0("^",x), "left")) x <- -1
    else x
  }))
  if (!all(unique(cens) %in% c(-1:1)))
    stop (paste0("Invalid censoring data. Accepted values are 'left', 'none', and 'right' \n",
                 "(abbreviations are allowed) or -1, 0, and 1. TRUE and FALSE are also accepted \n",
                 "and refer to 'right' and 'none' respectively."))
  cens
}
  