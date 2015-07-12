#' Parameters of interest for \code{brms} models
#' 
#' @inheritParams brm
#' @param ranef logical; indicating if random effects estimates should be returned
#' 
#' @return A vector of character strings specifying parameters of interest for models produced by the \code{brms} package.
#' 
#' @author Paul-Christian Buerkner \email{paul.buerkner@@gmail.com}
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
                    threshold = "flexible", predict = FALSE, ranef = TRUE) {
  family <- family[1]
  if (is.null(autocor)) autocor <- cor.arma()
  if (!is(autocor,"cor.brms")) stop("cor must be of class cor.brms")
  ee <- extract.effects(formula = formula, family = family, partial = partial)
  data <- updateData(data, family = family, effects = ee)
    
  is.lin <- family %in% c("gaussian", "student", "cauchy")
  is.ord <- family  %in% c("cumulative","cratio","sratio","acat")
  is.skew <- family %in% c("gamma", "weibull", "exponential")
  if (!(is.lin | is.ord | is.skew | family %in% 
      c("poisson", "negbinomial", "geometric", "binomial","bernoulli", "categorical", "multigaussian")))
    stop(paste(family,"is not a valid family"))
  
  f <- colnames(brm.model.matrix(ee$fixed, data, rm.int = is.ord))
  r <- lapply(lapply(ee$random, brm.model.matrix, data = data), colnames)
  p <- colnames(brm.model.matrix(partial, data, rm.int = TRUE))
  out <- NULL
  if (is.ord & threshold == "flexible") out <- c(out, "b_Intercept")
  if (is.ord & threshold == "equidistant") out <- c(out, "b_Intercept1", "delta")
  if (length(f) & family != "categorical") out <- c(out, "b")
  if (is.ord & length(p) | family == "categorical") out <- c(out, "bp")
  if (is.lin & !is(ee$se,"formula")) out <- c(out, "sigma")
  if (family == "multigaussian") out <- c(out, "sigma", "rescor")
  if (family == "student") out <- c(out,"nu")
  if (family %in% c("gamma","weibull","negbinomial")) out <- c(out,"shape")
  if (autocor$p > 0) out <- c(out,"ar")
  if (autocor$q > 0) out <- c(out,"ma")
  if (length(ee$group)) {
    out <- c(out, paste0("sd_",ee$group))
    out <- c(out, unlist(lapply(1:length(ee$group), function(i)
      if (length(r[[i]])>1) paste0("cor_",ee$group[[i]]))))
    if (ranef) out <- c(out, paste0("r_",ee$group))
  }
  if (predict) out <- c(out,"Y_pred")
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
#'          data = inhaler, family = "cumulative")
#' names(data1)
#' 
#' data2 <- brm.data(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|visit), 
#'          data = epilepsy, family = c("poisson", "log"))
#' names(data2)
#'          
#' @export
brm.data <- function(formula, data = NULL, family = c("gaussian", "identity"), prior = list(),
                     autocor = NULL, partial = NULL, cov.ranef = NULL) {
  link <- brm.link(family)
  family <- family[1]
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
  if (is(autocor, "cor.brms")) {
    if (family == "multigaussian" & sum(autocor$p, autocor$q) > 0 & 
        !any(sapply(c("__trait","trait__"), grepl, x = et$group)))
      stop("autocorrelation structure for family 'multigaussian' must contain 'trait' as a grouping variable")
    to.order <- rmNULL(list(data[["trait"]], data[[et$group]], data[[et$time]]))
    if (length(to.order)) 
      data <- data[do.call(order, to.order),]
  }
  
  is.lin <- family %in% c("gaussian", "student", "cauchy")
  is.ord <- family  %in% c("cumulative","cratio","sratio","acat")
  is.skew <- family %in% c("gamma", "weibull", "exponential")
  is.count <- family %in% c("poisson", "negbinomial", "geometric")
  if (!(is.lin | is.ord | is.skew | is.count | family %in% 
      c("binomial", "bernoulli", "categorical", "multigaussian")))
    stop(paste(family,"is not a valid family"))
  
  supl.data <- list(N = nrow(data), Y = model.response(data))
  if (family == "multigaussian") {
    supl.data$Y <- matrix(supl.data$Y, ncol = length(ee$response))
    supl.data <- c(supl.data, list(N_trait = nrow(supl.data$Y), K_trait = ncol(supl.data$Y)),
                   NC_trait = ncol(supl.data$Y) * (ncol(supl.data$Y)-1)/2) 
  }
  X <- brm.model.matrix(ee$fixed, data, rm.int = is.ord)
  if (is.ord | family == "categorical") {
    if (is.factor(supl.data$Y)) supl.data$Y <- as.numeric(supl.data$Y)
    else supl.data$Y <- supl.data$Y - min(supl.data$Y) + 1
  }
  else if (is.factor(supl.data$Y)) 
    stop(paste("family", family, "expects numeric response variable")) 
  
  if (is.formula(ee$se))
    supl.data <- c(supl.data,list(sigma = brm.model.matrix(ee$se, data, rm.int = TRUE)[,1])) 
  if (is.formula(ee$weights)) {
    supl.data <- c(supl.data, list(weights = brm.model.matrix(ee$weights, data, rm.int = TRUE)[,1]))
    if (family == "multigaussian") supl.data$weights <- supl.data$weights[1:supl.data$N_trait]
  }
  if (is.formula(ee$cens)) {
    cens <- brm.model.matrix(ee$cens, data, rm.int = TRUE)[,1]
    cens <- sapply(cens, function(x) {
      if (grepl(paste("^",x), "right") | is.logical(x) & x) x <- 1
      else if (grepl(paste("^",x), "none") | is.logical(x) & !x) x <- 0
      else if (grepl(paste("^",x), "left")) x <- -1
      else x
    })
    if (!all(unique(cens) %in% c(-1:1)))
      stop (paste0("Invalid censoring data. Accepted values are 'left', 'none', and 'right' \n",
                   "(abbreviations are allowed) or -1, 0, and 1. TRUE and FALSE are also accepted \n",
                   "and refer to 'right' and 'none' respectively."))
    supl.data <- c(supl.data, list(cens = cens))
  }
  if (is.ord | family %in% c("binomial", "categorical")) {
    if (family == "binomial") add <- ee$trials
    else add <- ee$cat
    if (!length(add)) supl.data$max_obs <- max(supl.data$Y)
    else if (is.numeric(add)) supl.data$max_obs <- add
    else if (is(add, "formula")) 
      supl.data$max_obs <- brm.model.matrix(add, data, rm.int = TRUE)[,1]
    else stop("Response part of formula is invalid")
    
    if (is(partial,"formula")) {
      if (family %in% c("sratio","cratio","acat")) {
        Xp <- brm.model.matrix(partial, data, rm.int = TRUE)
        supl.data <- c(supl.data, list(Kp = ncol(Xp), Xp = Xp))
        fp <- intersect(colnames(X), colnames(Xp))
        if (length(fp))
          stop(paste("Variables cannot be modeled as fixed and partial effects at the same time.",
                     "Error occured for variables:", paste(fp, collapse = ", ")))
      } 
      else stop("partial is only meaningful for families 'sratio', 'cratio', and 'acat'")  
    } 
    two.cat <- (family == "categorical" | is.ord) & max(supl.data$max_obs) == 2
    if (two.cat) supl.data$Y <- supl.data$Y - 1
    family <- ifelse(family == "binomial" & max(supl.data$Y) == 1 | two.cat, 
                     "bernoulli", family)
  } 
  
  if (length(ee$random)) {
    Z <- lapply(ee$random, brm.model.matrix, data = data)
    r <- lapply(Z, colnames)
    if (family != "categorical")
      to.zero <- unlist(lapply(unlist(lapply(r, intersect, y = colnames(X))), 
                               function(x) which(x == colnames(X)))) 
    else to.zero <- NULL
    X[,to.zero] <- 0
    ncolZ <- lapply(Z, ncol)
    expr <- expression(get(g, data), length(unique(get(g, data))), ncolZ[[i]], 
                       Z[[i]], ncolZ[[i]]*(ncolZ[[i]]-1)/2)
    for (i in 1:length(ee$group)) {
      g <- ee$group[[i]]
      name <- paste0(c("", "N_", "K_", "Z_", "NC_"), g)
      if (ncolZ[[i]] == 1) Z[[i]] <- as.vector(Z[[i]])
      for ( j in 1:length(name)) supl.data <- c(supl.data, setNames(list(eval(expr[j])), name[j]))
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
        cov.ranef[[g]] <- nrow(cov.ranef[[g]])/sum(diag(cov.ranef[[g]])) * cov.ranef[[g]]
        cov.ranef[[g]] <- t(chol(cov.ranef[[g]][order(level.names), order(level.names)]))
        if (length(r[[i]]) > 1) 
          diag(cov.ranef[[g]]) <- diag(cov.ranef[[g]]) - 1
        supl.data <- c(supl.data, setNames(list(cov.ranef[[g]]), paste0("CF_cov_",g)))
      }
    }
  }
  if (family == "categorical") {
    supl.data <- c(supl.data, list(Kp = ncol(X), Xp = X))
    X <- data.frame()
  }
  if (ncol(X) > 0) supl.data <- c(supl.data, list(K = ncol(X), X = X))
  
  if (autocor$p + autocor$q > 0) {
    time <- data[[et$time]]
    if (is.null(time)) time <- 1:nrow(data)
    tgroup <- data[[et$group]]
    if (is.null(tgroup)) tgroup <- rep(1, supl.data$N) 
    U_tgroup <- unique(tgroup)
    N_tgroup <- length(U_tgroup)
    if (autocor$p > 0 & is(autocor,"cor.arma")) {
      supl.data$Yar <- matrix(0, nrow = supl.data$N, ncol = autocor$p)
      supl.data$Kar <- autocor$p
      ptsum <- rep(0, N_tgroup + 1)
      mean.Y <- mean(supl.data$Y) 
      for (j in 1:N_tgroup) {
        ptsum[j+1] <- ptsum[j] + sum(tgroup == U_tgroup[j])
        for (i in 1:autocor$p) {
          if (ptsum[j]+i+1 <= ptsum[j+1])
            supl.data$Yar[(ptsum[j]+i+1):ptsum[j+1], i] <- 
            supl.data$Y[(ptsum[j]+1):(ptsum[j+1]-i)] - mean.Y
        }
      }
    }
    if (autocor$q > 0 & is(autocor,"cor.arma")) {
      supl.data$Ema_pre <- matrix(0, nrow = supl.data$N, ncol = autocor$q)
      supl.data$Kma <- autocor$q
      supl.data$tgroup <- as.numeric(as.factor(tgroup))
    }
  } 
  supl.data
}  