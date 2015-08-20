#' Parameters of interest for \code{brms} models (depricated)
#' 
#' @inheritParams brm
#' @param ranef logical; indicating if random effects estimates should be returned
#' 
#' @return A vector of character strings specifying parameters of interest for models produced by the \code{brms} package.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @details This function is depricated. Parameters of interest are now chosen by exclusion not by inclusion.
#' 
#' @examples 
#' brm.pars(rating ~ treat + period + carry + (1|subject),
#'          data = inhaler, family = "cumulative")
#
#  brm.pars(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|visit),
#         data = epilepsy, family = c("poisson", "log"))
#'               
#' @export
brm.pars = function(formula, data = NULL, family = "gaussian", autocor = NULL, partial = NULL,
                    threshold = "flexible", ranef = TRUE, predict = FALSE) {
  family <- family[1]
  if (is.null(autocor)) autocor <- cor.arma()
  if (!is(autocor,"cor.brms")) stop("cor must be of class cor.brms")
  ee <- extract.effects(formula = formula, family = family, partial = partial)
  data <- updateData(data, family = family, effects = ee)
    
  is.linear <- family %in% c("gaussian", "student", "cauchy")
  is.ordinal <- family  %in% c("cumulative","cratio","sratio","acat")
  is.skew <- family %in% c("gamma", "weibull", "exponential")
  if (!(is.linear || is.ordinal || is.skew || family %in% 
      c("poisson", "negbinomial", "geometric", "binomial","bernoulli", "categorical")))
    stop(paste(family,"is not a valid family"))
  
  f <- colnames(brm.model.matrix(ee$fixed, data, rm.int = is.ordinal))
  r <- lapply(lapply(ee$random, brm.model.matrix, data = data), colnames)
  p <- colnames(brm.model.matrix(partial, data, rm.int = TRUE))
  out <- NULL
  if (is.ordinal && threshold == "flexible") out <- c(out, "b_Intercept")
  if (is.ordinal && threshold == "equidistant") out <- c(out, "b_Intercept1", "delta")
  if (length(f) && family != "categorical") out <- c(out, "b")
  if (is.ordinal && length(p) || family == "categorical") out <- c(out, "bp")
  if (is.linear && !is(ee$se,"formula") && length(ee$response) == 1) out <- c(out, "sigma")
  if (family == "gaussian" && length(ee$response) > 1) out <- c(out, "sigma", "rescor")
  if (family == "student") out <- c(out,"nu")
  if (family %in% c("gamma", "weibull", "negbinomial")) out <- c(out,"shape")
  if (autocor$p > 0) out <- c(out,"ar")
  if (autocor$q > 0) out <- c(out,"ma")
  if (length(ee$group)) {
    out <- c(out, paste0("sd_",ee$group))
    out <- c(out, unlist(lapply(1:length(ee$group), function(i)
      if (length(r[[i]])>1 && ee$cor[[i]]) paste0("cor_",ee$group[[i]]))))
    if (ranef) out <- c(out, paste0("r_",ee$group))
  }
  if (predict) out <- c(out, "Y_pred")
  return(out)
}

#' Extract required data for \code{brms} models
#'
#' @inheritParams brm
#' 
#' @return A named list of objects containing the required data to fit a \code{brms} model 
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
#' 
#' @examples
#' data1 <- brm.data(rating ~ treat + period + carry + (1|subject), 
#'                   data = inhaler, family = "cumulative")
#' names(data1)
#' 
#' data2 <- brm.data(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|visit), 
#'                   data = epilepsy, family = "poisson")
#' names(data2)
#'          
#' @export
brm.data <- function(formula, data = NULL, family = "gaussian", prior = list(),
                     autocor = NULL, partial = NULL, cov.ranef = NULL) {
  family <- family[1]
  is.linear <- family %in% c("gaussian", "student", "cauchy")
  is.ordinal <- family  %in% c("cumulative","cratio","sratio","acat")
  is.count <- family %in% c("poisson", "negbinomial", "geometric")
  is.skew <- family %in% c("gamma", "weibull", "exponential")
  if (family == "multigaussian") 
    stop("family 'multigaussian' is depricated. Use family 'gaussian' instead")
  if (!(is.linear | is.ordinal | is.skew | is.count | family %in% 
        c("binomial", "bernoulli", "categorical")))
    stop(paste(family, "is not a valid family"))
  if (is.null(autocor)) autocor <- cor.arma()
  if (!is(autocor,"cor.brms")) stop("cor must be of class cor.brms")
  
  et <- extract.time(autocor$formula)
  ee <- extract.effects(formula = formula, family = family, partial, et$all)
  data <- updateData(data, family = family, effects = ee, et$group)
  group.names <- list()
  for (g in ee$group) { 
    group.names[[g]] <- sort(as.character(unique(data[[g]])))
    data[[g]] <- as.numeric(as.factor(data[[g]]))
  } 
  
  #sort data in case of autocorrelation models
  if (sum(autocor$p, autocor$q) > 0) {
    if (family == "gaussian" && length(ee$response) > 1 && !any(sapply(c("__trait","trait__"), grepl, x = et$group)))
      stop("autocorrelation structures for multiple responses must contain 'trait' as grouping variable")
    to.order <- rmNULL(list(data[["trait"]], data[[et$group]], data[[et$time]]))
    if (length(to.order)) 
      data <- data[do.call(order, to.order),]
  }
  
  #response variable
  standata <- list(N = nrow(data), Y = unname(model.response(data)))
  if (!is.numeric(standata$Y) && !(is.ordinal || family %in% c("bernoulli", "categorical"))) 
    stop(paste("family", family, "expects numeric response variable"))
  
  #transform and check response variable for different families
  if (is.count || family == "binomial") {
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
  else if (is.ordinal) {
    if (is.factor(standata$Y)) {
      if (is.ordered(standata$Y)) standata$Y <- as.numeric(standata$Y)
      else stop(paste("family", family, "requires factored response variables to be ordered"))
    }
    else if (all(is.wholenumber(standata$Y)))
      standata$Y <- standata$Y - min(standata$Y) + 1
    else stop(paste("family", family, "expects either integers or ordered factors as response variables"))
  }
  else if (is.skew) {
    if (min(standata$Y) < 0)
      stop(paste("family", family, "requires response variable to be non-negative"))
  }  
  else if (family == "gaussian" && length(ee$response) > 1) {
    standata$Y <- matrix(standata$Y, ncol = length(ee$response))
    standata <- c(standata, list(N_trait = nrow(standata$Y), K_trait = ncol(standata$Y)),
                   NC_trait = ncol(standata$Y) * (ncol(standata$Y)-1)/2) 
  }
  
  #fixed effects data
  X <- brm.model.matrix(ee$fixed, data, rm.int = is.ordinal)
  if (family == "categorical") standata <- c(standata, list(Kp = ncol(X), Xp = X))
  else standata <- c(standata, list(K = ncol(X), X = X))
  
  #random effects data
  if (length(ee$random)) {
    Z <- lapply(ee$random, brm.model.matrix, data = data)
    r <- lapply(Z, colnames)
    if (family != "categorical")
      to.zero <- unlist(lapply(unlist(lapply(r, intersect, y = colnames(X))), 
                               function(x) which(x == colnames(X)))) 
    else to.zero <- NULL
    ncolZ <- lapply(Z, ncol)
    expr <- expression(get(g, data), length(unique(get(g, data))), ncolZ[[i]], 
                       Z[[i]], ncolZ[[i]]*(ncolZ[[i]]-1)/2)
    for (i in 1:length(ee$group)) {
      g <- ee$group[[i]]
      name <- paste0(c("", "N_", "K_", "Z_", "NC_"), g)
      if (ncolZ[[i]] == 1) Z[[i]] <- as.vector(Z[[i]])
      for ( j in 1:length(name)) standata <- c(standata, setNames(list(eval(expr[j])), name[j]))
      if (g %in% names(cov.ranef)) {
        cov.ranef[[g]] <- as.matrix(cov.ranef[[g]])
        level.names <- rownames(cov.ranef[[g]])
        colnames(cov.ranef[[g]]) <- level.names
        if (is.null(level.names)) 
          stop(paste("Row names are required for covariance matrix of",g))
        if (nrow(cov.ranef[[g]]) != length(group.names[[g]]))
          stop(paste("Dimension of covariance matrix of",g,"is incorrect"))
        if (any(sort(level.names) != group.names[[g]]))
          stop(paste("Row names of covariance matrix of",g,"do not match names of the grouping levels"))
        if (!isSymmetric(unname(cov.ranef[[g]])))
          stop(paste("Covariance matrix of grouping factor",g,"is not symmetric"))
        if (min(eigen(cov.ranef[[g]], symmetric = TRUE, only.values = TRUE)$values) <= 0)
          warning(paste("Covariance matrix of grouping factor",g,"may not be positive definite"))
        cov.ranef[[g]] <- cov.ranef[[g]][order(level.names), order(level.names)]
        if (length(r[[i]]) == 1) 
          cov.ranef[[g]] <- t(suppressWarnings(chol(cov.ranef[[g]], pivot = TRUE)))
        else if (length(r[[i]]) > 1 && !ee$cor[[i]])
          cov.ranef[[g]] <- t(suppressWarnings(chol(kronecker(cov.ranef[[g]], diag(ncolZ[[i]])), pivot = TRUE)))
        standata <- c(standata, setNames(list(cov.ranef[[g]]), paste0("cov_",g)))
      }
    }
  }
  
  #addition and partial variables
  if (is.formula(ee$se)) {
    standata <- c(standata, list(sigma = unname(brm.model.matrix(ee$se, data, rm.int = TRUE)[,1])))
    if (min(standata$sigma) < 0) stop("standard errors must be non-negative")
  }
  if (is.formula(ee$weights)) {
    standata <- c(standata, list(weights = unname(brm.model.matrix(ee$weights, data, rm.int = TRUE)[,1])))
    if (family == "gaussian" && length(ee$response) > 1) standata$weights <- standata$weights[1:standata$N_trait]
    if (min(standata$weights) < 0) stop("weights must be non-negative")
  }
  if (is.formula(ee$cens)) {
    cens <- get(all.vars(ee$cens)[1], data)
    if (is.factor(cens)) cens <- as.character(cens)
    cens <- unname(sapply(cens, function(x) {
      if (grepl(paste0("^",x), "right") || is.logical(x) && isTRUE(x)) x <- 1
      else if (grepl(paste0("^",x), "none") || is.logical(x) && !isTRUE(x)) x <- 0
      else if (grepl(paste0("^",x), "left")) x <- -1
      else x
    }))
    if (!all(unique(cens) %in% c(-1:1)))
      stop (paste0("Invalid censoring data. Accepted values are 'left', 'none', and 'right' \n",
                   "(abbreviations are allowed) or -1, 0, and 1. TRUE and FALSE are also accepted \n",
                   "and refer to 'right' and 'none' respectively."))
    standata <- c(standata, list(cens = cens))
  }
  if (is.ordinal || family %in% c("binomial", "categorical")) {
    if (family == "binomial") add <- ee$trials
    else add <- ee$cat
    if (!length(add)) standata$max_obs <- max(standata$Y)
    else if (is.numeric(add)) standata$max_obs <- add
    else if (is.formula(add)) 
      standata$max_obs <- unname(brm.model.matrix(add, data, rm.int = TRUE)[,1])
    else stop("Response part of formula is invalid.")
    if (any(standata$Y > standata$max_obs))
      stop("The number of trials / categories is smaller the response variable would suggest.")
    if ((is.ordinal || family == "categorical") && max(standata$max_obs) == 2 ||
        family == "binomial" && max(standata$max_obs) == 1) 
      message("Only 2 levels detected so that family 'bernoulli' might be a more efficient choice.")
  } 
  
  #get data for partial effects
  if (is.formula(partial)) {
    if (family %in% c("sratio","cratio","acat")) {
      Xp <- brm.model.matrix(partial, data, rm.int = TRUE)
      standata <- c(standata, list(Kp = ncol(Xp), Xp = Xp))
      fp <- intersect(colnames(X), colnames(Xp))
      if (length(fp))
        stop(paste("Variables cannot be modeled as fixed and partial effects at the same time.",
                   "Error occured for variables:", paste(fp, collapse = ", ")))
    } 
    else stop("partial is only meaningful for families 'sratio', 'cratio', and 'acat'")  
  }
  
  #autocorrelation variables
  if (autocor$p + autocor$q > 0) {
    time <- data[[et$time]]
    if (is.null(time)) time <- 1:nrow(data)
    tgroup <- data[[et$group]]
    if (is.null(tgroup)) tgroup <- rep(1, standata$N) 
    U_tgroup <- unique(tgroup)
    N_tgroup <- length(U_tgroup)
    if (autocor$p > 0 && is(autocor,"cor.arma")) {
      standata$Yar <- matrix(0, nrow = standata$N, ncol = autocor$p)
      standata$Kar <- autocor$p
      ptsum <- rep(0, N_tgroup + 1)
      mean.Y <- mean(standata$Y) 
      for (j in 1:N_tgroup) {
        ptsum[j+1] <- ptsum[j] + sum(tgroup == U_tgroup[j])
        for (i in 1:autocor$p) {
          if (ptsum[j]+i+1 <= ptsum[j+1])
            standata$Yar[(ptsum[j]+i+1):ptsum[j+1], i] <- 
            standata$Y[(ptsum[j]+1):(ptsum[j+1]-i)] - mean.Y
        }
      }
    }
    if (autocor$q > 0 && is(autocor,"cor.arma")) {
      standata$Ema_pre <- matrix(0, nrow = standata$N, ncol = autocor$q)
      standata$Kma <- autocor$q
      standata$tgroup <- as.numeric(as.factor(tgroup))
    }
  } 
  standata
}  