# convert array to list of elements with reduced dimension
#
# @param x an arrary of dimension d
#
# @return a list of arrays of dimension d-1
array2list <- function(x) {
  if (is.null(dim(x))) stop("Argument x has no dimension")
  n.dim <- length(dim(x))
  l <- list(length = dim(x)[n.dim])
  ind <- collapse(rep(",",n.dim-1))
  for (i in 1:dim(x)[n.dim])
    l[[i]] <- eval(parse(text = paste0("x[",ind,i,"]")))
  names(l) <- dimnames(x)[[n.dim]]
  l
}

# convert list to array of increased dimension
#
# @param x a list containing arrary of dimension d
#
# @return an arrary of dimension d+1
list2array <- function(x) {
  if (!is.list(x) || length(x) == 0) 
    stop("x must be a non-empty list")
  x <- unlist(lapply(x, array2list), recursive = FALSE)
  dim_elements <- lapply(x, function(y) if (!is.null(dim(y))) dim(y) else length(y))
  dim_target <- dim_elements[[1]]
  if (!all(sapply(dim_elements, all.equal, current = dim_target)))
    stop("dimensions of list elements do not match")
  a <- array(NA, dim = c(dim_target, length(x)))
  ind <- collapse(rep(",",length(dim_target)))
  for (i in 1:length(x)) 
    eval(parse(text = paste0("a[",ind,i,"] <- x[[",i,"]]")))
  if (length(x) == 1) dimnames(a)[[length(dim_target)+1]] <- list(names(x))
  else dimnames(a)[[length(dim_target)+1]] <- names(x)
  a
}

# find the first element in A that is greater than target
#
# @param A a matrix
# @param target a vector of length nrow(A)
# @param i column of A being checked first
#
# @return a vector of the same length as target 
#   containing the column ids where A[,i] was first greater than target
first_greater <- function(A, target, i = 1) {
  ifelse(target <= A[,i] | ncol(A) == i, i, first_greater(A, target, i+1))
}

# apply the link function on x
link <- function(x, link) {
  if (link == "identity") x
  else if (link == "log") log(x)
  else if (link == "inverse") 1/x
  else if (link == "sqrt") sqrt(x)
  else if (link == "logit") logit(x)
  else if (link == "probit") qnorm(x)
  else if (link == "cloglog") log(-log(1-x))
  else if (link == "probit_approx") qnorm(x)
  else stop(paste("Link", link, "not supported"))
}

# apply the inverse link function on x
ilink <- function(x, link) {
  if (link == "identity") x
  else if (link == "log") exp(x)
  else if (link == "inverse") 1/x
  else if (link == "sqrt") x^2
  else if (link == "logit") ilogit(x)
  else if (link == "probit") pnorm(x)
  else if (link == "cloglog") 1 - exp(-exp(x))
  else if (link == "probit_approx") ilogit(0.07056*x^3 + 1.5976*x)
  else stop(paste("Link", link, "not supported"))
}

# get correlation names as combinations of variable names
# 
# @param names the variable names 
# @param type of the correlation to be put in front of the returned strings
# @param brackets should the correlation names contain brackets or underscores as seperators
# @param subset subset of correlation parameters to be returned. Currently only used in summary.brmsfit (s3.methods.R)
# @param the subtype of the correlation (e.g., g1 in cor_g1_x_y). Only used when subset is not NULL
get_cornames <- function(names, type = "cor", brackets = TRUE, subset = NULL, subtype = "") {
  cor_names <- NULL
  if (is.null(subset) && length(names) > 1) {
    for (i in 2:length(names)) {
      for (j in 1:(i-1)) {
        if (brackets) cor_names <- c(cor_names, paste0(type,"(",names[j],",",names[i],")"))
        else cor_names <- c(cor_names, paste0(type,"_",names[j],"_",names[i]))
      }
    }
  } else if (!is.null(subset)) {
    possible_values <- get_cornames(names = names, type = "", brackets = FALSE)
    subset <- rename(subset, paste0("^",type, if (nchar(subtype)) paste0("_",subtype)),
                     "", fixed = FALSE)
    matches <- which(possible_values %in% subset)
    cor_names <- get_cornames(names = names, type = type)[matches]
  }
  cor_names
}

# calculate estimates over posterior samples 
# 
# @param coef coefficient to be applied on the samples (e.g., "mean")
# @param samples the samples over which to apply coef
# @param margin see apply
# @param to.array logical; should theresult be transformed into an array of increased dimension?
# @param ... additional arguments passed to get(coef)
#
# @return typically a matrix with colnames(samples) as colnames
get_estimate <- function(coef, samples, margin = 2, to.array = FALSE, ...) {
  dots <- list(...)
  args <- list(X = samples, MARGIN = margin, FUN = coef)
  fun.args <- names(formals(coef))
  if (!"..." %in% fun.args)
    dots <- dots[names(dots) %in% fun.args]
  x <- do.call(apply, c(args, dots))
  if (is.null(dim(x))) 
    x <- matrix(x, dimnames = list(NULL, coef))
  else if (coef == "quantile") x <- aperm(x, length(dim(x)):1)
  if (to.array && length(dim(x)) == 2) 
    x <- array(x, dim = c(dim(x), 1), dimnames = list(NULL, NULL, coef))
  x 
}

#compute covariance and correlation matrices based on correlation and sd samples
#
# @param sd samples of standard deviations
# @param cor samples of correlations
#
# @details used in VarCorr.brmsfit
# 
# @return samples of covariance and correlation matrices
get_cov_matrix <- function(sd, cor = NULL) {
  if (any(sd < 0)) stop("standard deviations must be non negative")
  if (!is.null(cor)) {
    if (ncol(cor) != ncol(sd)*(ncol(sd)-1)/ 2 || nrow(sd) != nrow(cor))
      stop("dimensions of standard deviations and corrrelations do not match")
    if (any(cor < -1 || cor > 1)) 
      stop("correlations must be between -1 and 1")  
  }
  nsamples <- nrow(sd)
  nranef <- ncol(sd)
  cor_matrix <- cov_matrix <- aperm(array(diag(1, nranef), dim = c(nranef, nranef, nsamples)), c(3,1,2))
  for (i in 1:nranef) 
    cov_matrix[,i,i] <- sd[,i]^2 
  if (!is.null(cor)) {
    k <- 0 
    for (i in 2:nranef) {
      for (j in 1:(i-1)) {
        k = k + 1
        cor_matrix[,j,i] <- cor_matrix[,i,j] <- cor[,k]
        cov_matrix[,j,i] <- cov_matrix[,i,j] <- cor[,k] * sd[,i] * sd[,j]
      }
    }
  }
  list(cor = cor_matrix, cov = cov_matrix)
}

# calculate the evidence ratio between two disjunct hypotheses
# 
# @param x posterior samples 
# @param cut the cut point between the two hypotheses
# @param wsign direction of the hypothesis
# @param prior_samples optional prior samples for undirected hypothesis
# @param pow influence the accuracy of the density
# @param ... optional arguments passed to density.default
#
# @return the evidence ratio of the two hypothesis
evidence_ratio <- function(x, cut = 0, wsign = c("equal", "less", "greater"), 
                           prior_samples = NULL, pow = 12, ...) {
  wsign <- match.arg(wsign)
  if (wsign == "equal") {
    if (is.null(prior_samples)) out <- NA
    else {
      dots <- list(...)
      dots <- dots[names(dots) %in% names(formals("density.default"))]
      # compute prior and posterior densities
      prior_density <- do.call(density, c(list(x = prior_samples, n = 2^pow), dots))
      posterior_density <- do.call(density, c(list(x = x, n = 2^pow), dots))
      at_cut_prior <- match(min(abs(prior_density$x - cut)), abs(prior_density$x - cut))
      # evaluate densities at the cut point
      at_cut_posterior <- match(min(abs(posterior_density$x - cut)), abs(posterior_density$x - cut))
      out <- posterior_density$y[at_cut_posterior] / prior_density$y[at_cut_prior] 
    }
  }
  else if (wsign == "less") {
    out <- length(which(x < cut))
    out <- out / (length(x) - out)
  }  
  else if (wsign == "greater") {
    out <- length(which(x > cut))
    out <- out / (length(x) - out)  
  }
  out  
}

# compute eta for fixed effects
#
# @param X fixed effects design matrix
# @param b fixed effects samples
# 
# @return linear predictor for fixed effects
fixef_predictor <- function(X, b) {
  as.matrix(b) %*% t(as.matrix(X))
}

# compute eta for random effects
#  
# @param Z random effects design matrix
# @param gf levels of grouping factor for each observation
# @param r random effects samples
#
# @return linear predictor for random effects
ranef_predictor <- function(Z, gf, r) {
  Z <- expand_matrix(Z, gf)
  nlevels <- length(unique(gf))
  sort_levels <- unlist(lapply(1:nlevels, function(n) seq(n, ncol(r), nlevels)))
  as.matrix(r[, sort_levels]) %*% t(Z)
}

#compute eta for moving average effects
#
# @param data the data initially passed to stan
# @param ma moving average samples 
# @param eta previous linear predictor samples
# @param link
#
# @param new linear predictor samples updated by moving average effects
ma_predictor <- function(data, ma, eta, link = "identity") {
  ma <- as.matrix(ma)
  K <- ncol(ma)
  Ks <- 1:K
  Y <- link(data$Y, link)
  N <- length(Y)
  tg <- c(rep(0, K), data$tgroup)
  Ema <- array(0, dim = c(nrow(ma), K, N))
  e <- matrix(0, nrow = nrow(ma), ncol = N)
  for (n in 1:N) {
    eta[,n] <- eta[,n] + apply(ma * Ema[,,n], 1, sum)
    e[,n] <- Y[n] - eta[,n]
    if (n < N) {
      I <- which(n < N & tg[n+1+K] == tg[n+1+K-Ks])
      Ema[,I,n+1] <- e[,n+1-I]
    }
  }
  eta
}

# compute etap for partial and categorical effects
# 
# @param Xp partial design matrix 
# @param p partial effects samples
# @param max_obs number of categories
#
# @return linear predictor of partial effects as a 3D array (not as a matrix)
partial_predictor <- function(Xp, p, max_obs) {
  max_obs <- max(max_obs)
  etap <- array(0, dim = c(nrow(p), nrow(Xp), max_obs-1))
  for (k in 1:(max_obs-1)) {
    etap[,,k] <- as.matrix(p[, seq(k, (max_obs-1)*ncol(Xp), max_obs-1)]) %*% t(as.matrix(Xp))
  }
  etap
}

# expand a matrix into a sparse matrix of higher dimension
# 
# @param A matrix to be expanded
# @param x levels to expand the matrix
# 
# @details used in linear_predictor.brmsfit
#
# @return An expanded matrix of dimensions nrow(A) and ncol(A) * length(unique(x)) 
expand_matrix <- function(A, x) {
  A <- as.matrix(A)
  if (length(x) != nrow(A))
    stop("x must have nrow(A) elements")
  if (!all(is.wholenumber(x) & x > 0))
    stop("x must contain positive integers only")
  K <- ncol(A)
  v <- rep(0, K * max(x))
  do.call(rbind, lapply(1:nrow(A), function(n, v) {
    v[K*(x[n]-1) + 1:K] <- A[n, ] 
    return(v)}, v = v))
}

# amend new data for predict method
# 
# @param newdata a data.frame containing new data for prediction 
# @param formula a brms model formula
# @param family
# @param autocor
# @param partial
#
# @details used in predict.brmsfit and linear_predictor.brmsfit
#
# @return updated data.frame being compatible with formula
amend_newdata <- function(newdata, formula, family = "gaussian", autocor = cor_arma(),
                          partial = NULL) {
  ee <- extract_effects(formula, family = family)
  if (length(ee$group))
    stop("random effects models not yet supported for predicting new data")
  if (sum(autocor$p, autocor$q) > 0 && !all(ee$response %in% names(newdata))) 
    stop("response variables must be specified in newdata for autocorrelative models")
  else for (resp in ee$response) newdata[[resp]] <- 0 # add irrelevant response variables
  if (is.formula(ee$cens)) 
    for (cens in all.vars(ee$cens)) newdata[[cens]] <- 0 # add irrelevant censor variables
    data <- brmdata(formula, data = newdata, family = family, autocor = autocor, partial = partial)
}

# find all valid object names in a string (used in method hypothesis in s3.methods.R)
# 
# @x a character string
#
# @details currently used in hypothesis.brmsfit
#
# @return all valid variable names within the string
find_names <- function(x) {
  if (!is.character(x) || length(x) > 1) stop("x must be a character string of length 1")
  x <- gsub(" ", "", x)
  pos_fun <- gregexpr("([^([:digit:]|[:punct:])]|\\.)[[:alnum:]_\\.]*\\(", x)[[1]]
  pos_decnum <- gregexpr("\\.[[:digit:]]+", x)[[1]]
  pos_var <- list(rmMatch(gregexpr("([^([:digit:]|[:punct:])]|\\.)[[:alnum:]_\\.]*(\\[[[:digit:]]*\\])?", x)[[1]], 
                          pos_fun, pos_decnum))
  unlist(regmatches(x, pos_var))
}

# trace and density plots for one parameter
#
# @param par a single character string 
# @parm x a data.frame containing the samples
#
# @return trace and density plot for this parameter
td_plot <- function(par, x) {
  if (!is.character(par) || length(par) != 1)
    stop("par must be a character string")
  if (!is.data.frame(x))
    stop("x must be a data.frame")
  names(x)[match(par, names(x))] <- "value" 
  trace <- ggplot(x, aes_string(x = "iter", y = "value", group = "chains", colour = "chains")) +
    geom_line(alpha = 0.7) + 
    xlab("") + ylab("") + ggtitle(paste("Trace of", par)) + 
    theme(legend.position = "none",
          plot.title = element_text(size = 15, vjust = 1),
          plot.margin = grid::unit(c(0.2, 0, -0.8, -0.5), "lines"))
  density <- ggplot(x, aes_string(x = "value")) + 
    geom_density(aes_string(fill = "chains"), alpha = 0.5) + 
    xlab("") + ylab("") + ggtitle(paste("Density of", par)) + 
    theme(plot.title = element_text(size = 15, vjust = 1),
          plot.margin = grid::unit(c(0.2, 0, -0.8, -0.5), "lines"))
  return(gridExtra::arrangeGrob(trace, density, ncol = 2, nrow = 1))
}

#' @export
print.brmssummary <- function(x, digits = 2, ...) {
  cat(paste0(" Family: ", x$family, " (", x$link, ") \n"))
  cat(paste("Formula:", gsub(" {1,}", " ", Reduce(paste, deparse(x$formula))), "\n"))
  cat(paste0("   Data: ", x$data.name, " (Number of observations: ",x$nobs,") \n"))
  if (x$sampler == "") {
    cat(paste("\nThe model does not contain posterior samples."))
  }
  else {
    cat(paste0("Samples: ", x$n.chains, " chains, each with n.iter = ", x$n.iter, 
               "; n.warmup = ", x$n.warmup, "; n.thin = ", x$n.thin, "; \n",
               "         total post-warmup samples = ", (x$n.iter-x$n.warmup)/x$n.thin*x$n.chains, "\n"))
    cat(paste0("   WAIC: ", ifelse(is.numeric(x$WAIC), round(x$WAIC, digits = digits), x$WAIC), "\n \n"))
    
    if (length(x$group)) {
      cat("Random Effects: \n")
      for (i in 1:length(x$group)) {
        cat(paste0("~",x$group[i], " (Number of levels: ",x$ngrps[[x$group[i]]],") \n"))
        x$random[[x$group[i]]][,"Eff.Sample"] <- round(x$random[[x$group[i]]][,"Eff.Sample"], digits = 0)
        print(round(x$random[[x$group[i]]], digits = digits))
        cat("\n")
      }
    }
    
    if (nrow(x$cor_pars)) {
      cat("Correlation Structure: "); print(x$autocor); cat("\n")
      x$cor_pars[,"Eff.Sample"] <- round(x$cor_pars[,"Eff.Sample"], digits = 0)
      print(round(x$cor_pars, digits = digits))
      cat("\n")
    }
    
    cat("Fixed Effects: \n")
    x$fixed[,"Eff.Sample"] <- round(x$fixed[,"Eff.Sample"], digits = 0)
    print(round(x$fixed, digits = digits)) 
    cat("\n")
    
    if (nrow(x$spec_pars)) {
      cat("Family Specific Parameters: \n")
      x$spec_pars[,"Eff.Sample"] <- round(x$spec_pars[,"Eff.Sample"], digits = 0)
      print(round(x$spec_pars, digits = digits))
      cat("\n")
    }
    
    cat(paste0("Samples were drawn using ",x$sampler,". For each parameter, Eff.Sample is a \n",
               "crude measure of effective sample size, and Rhat is the potential scale \n",
               "reduction factor on split chains (at convergence, Rhat = 1)."))
  }
}  

#' @export
print.brmshypothesis <- function(x, digits = 2, ...) {
  cat(paste0("Hypothesis Tests for class ", x$class, ":\n"))
  x$hypothesis[,1:5] <- round(x$hypothesis[,1:5], digits = digits)
  print(x$hypothesis, quote = FALSE)
  cat(paste0("---\n'*': The expected value under the hypothesis lies outside the ",(1-x$alpha)*100,"% CI."))
}

#' @export
print.brmsmodel <- function(x, ...) cat(x)

#' @export
print.ic <- function(x, digits = 2, ...) {
  ic <- names(x)[3]
  mat <- matrix(c(x[[ic]], x[[paste0("se_",ic)]]), ncol = 2, 
                dimnames = list("", c(toupper(ic), "SE")))
  print(round(mat, digits = digits))
}

#' @export
print.iclist <- function(x, digits = 2, ...) {
  ic <- names(x[[1]])[3]
  mat <- matrix(0, nrow = length(x), ncol = 2, 
                dimnames = list(names(x), c(toupper(ic), "SE")))
  for (i in 1:length(x))
    mat[i, ] <- c(x[[i]][[ic]], x[[i]][[paste0("se_",ic)]])
  if (is.matrix(attr(x, "compare")))
    mat <- rbind(mat, attr(x, "compare"))
  print(round(mat, digits = digits))
}