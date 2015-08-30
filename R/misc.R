# convert array to list of elements with reduced dimension
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

#convert list to array of increased dimension
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

isNULL <- function(x) is.null(x) || ifelse(is.vector(x), all(sapply(x, is.null)), FALSE)

rmNULL <- function(x) {
  x <- Filter(Negate(isNULL), x)
  lapply(x, function(x) if (is.list(x)) rmNULL(x) else x)
}

rmNum <- function(x) x[sapply(x, Negate(is.numeric))]

#remove all elements in x that also appear in ... while keeping all attributes
rmMatch <- function(x, ...) {
  att <- attributes(x)
  keep <- which(!(x %in% c(...)))
  x <- x[keep]
  attributes(x) <- att
  attr(x, "match.length") <- att$match.length[keep] 
  x
} 

#check if x is a whole number (integer)
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {  
  if (!is.numeric(x)) return(FALSE)
  else return(abs(x - round(x)) < tol)
}  

#wrapper for paste with collapse
collapse <- function(..., sep = "")
  paste(..., sep = sep, collapse = "")


#find the first element in row that is greater than target (function is vectorized)
first_greater <- function(A, target, i = 1) {
  ifelse(target <= A[,i] | ncol(A) == i, i, first_greater(A, target, i+1))
}

#compute the logit
logit <- function(p) {
  log(p/(1-p))
}

#compute the inverse of logit
ilogit <- function(x) { 
  exp(x) / (1 + exp(x))
}

#apply the link function on x
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

#apply the inverse link function on x
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

#calculate estimates over posterior samples 
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
cov_matrix <- function(sd, cor = NULL) {
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

#calculate the evidence ratio between two disjunct hypotheses
eratio <- function(x, cut = 0, wsign = c("equal", "less", "greater"), prior_samples = NULL, pow = 12, ...) {
  wsign <- match.arg(wsign)
  if (wsign == "equal") 
    if (is.null(prior_samples)) out <- NA
    else {
      dots <- list(...)
      dots <- dots[names(dots) %in% names(formals("density.default"))]
      prior_density <- do.call(density, c(list(x = prior_samples, n = 2^pow), dots))
      posterior_density <- do.call(density, c(list(x = x, n = 2^pow), dots))
      at_cut_prior <- match(min(abs(prior_density$x - cut)), abs(prior_density$x - cut))
      at_cut_posterior <- match(min(abs(posterior_density$x - cut)), abs(posterior_density$x - cut))
      out <- posterior_density$y[at_cut_posterior] / prior_density$y[at_cut_prior] 
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

# Links for \code{brms} families
# 
# @param family A vector of one or two character strings. The first string indicates the distribution of the dependent variable (the 'family'). Currently, the following distributions are supported:
#  \code{"gaussian"}, \code{"student"}, \code{"cauchy"}, \code{"poisson"}, \code{"binomial"}, \code{"categorical"}, 
#  \code{"gamma"}, \code{"exponential"}, \code{"weibull"}, \code{"cumulative"}, \cr
#  \code{"cratio"}, \code{"sratio"}, and \code{"acat"}.
#  The second string indicates the link function, which must supported by the distribution of the dependent variable. If not specified, default link functions are used (see 'Details').
# @return The second element of \code{family} (if present) or else the default link of the specified family
#   
# @details The families \code{gaussian}, \code{student}, and \code{cauchy} accept the links (as names) \code{identity}, \code{log}, and \code{inverse};
# the \code{poisson} family the links \code{log}, \code{identity}, and \code{sqrt}; 
# families \code{binomial}, \code{cumulative}, \code{cratio}, \code{sratio}, and \code{acat} the links \code{logit}, \code{probit}, \code{probit_approx}, and \code{cloglog};
# family  \code{categorical} the link \code{logit}; families \code{gamma}, \code{weibull}, and \code{exponential} the links \code{log}, \code{identity}, and \code{inverse}. 
# The first link mentioned for each family is the default.     
# 
# @examples brm.link("gaussian")
# brm.link(c("gaussian","log"))
link4family <- function(family) {
  link <- family[2]
  family <- family[1]
  is_linear <- family %in% c("gaussian", "student", "cauchy")
  is_skew <- family %in% c("gamma", "weibull", "exponential")
  is_cat <- family %in% c("cumulative", "cratio", "sratio", "acat", "binomial", "bernoulli")                    
  is_count <- family %in% c("poisson", "negbinomial", "geometric")

  if (is.na(link)) {
    if (is_linear) link <- "identity"
    else if (is_skew || is_count) link <- "log"
    else if (is_cat || family == "categorical") link <- "logit"
  }
  else if (is_linear && !link %in% c("identity", "log", "inverse") ||
           is_count && !link %in% c("log", "identity", "sqrt") ||
           is_cat && !link %in% c("logit", "probit", "probit_approx", "cloglog") ||
           family == "categorical" && link != "logit" ||
           is_skew && !link %in% c("log", "identity", "inverse"))
    stop(paste(link, "is not a valid link for family", family))
  else if (is_count && link == "sqrt") 
    warning(paste(family, "model with sqrt link may not be uniquely identified"))
  link
}

#check validity of family
check_family <- function(family) {
  if (family == "normal") family <- "gaussian"
  if (family == "multigaussian") 
    stop("family 'multigaussian' is deprecated. Use family 'gaussian' instead")
  if (!family %in% c("gaussian", "student", "cauchy", "binomial", "bernoulli", "categorical",
                     "poisson", "negbinomial", "geometric", "gamma", "weibull", "exponential",
                     "cumulative", "cratio", "sratio", "acat"))
    stop(paste(family, "is not a valid family"))
  family
}

#list irrelevant parameters not to be saved by Stan
exclude_pars <- function(formula, ranef = TRUE) {
  ee <- extract_effects(formula = formula, add.ignore = TRUE)
  out <- c("eta", "etam", "etap", "b_Intercept1", "Lrescor", "Rescor",
           "p", "q", "e", "Ema", "lp_pre")
  if (length(ee$group)) {
    for (i in 1:length(ee$group)) {
      out <- c(out, paste0("pre_",i), paste0("L_",i), paste0("Cor_",i))
      if (!ranef) out <- c(out, paste0("r_",i))
    }
  }
  out
}