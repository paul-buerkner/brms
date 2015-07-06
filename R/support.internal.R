# Construct Design Matrices for \code{brms} models
# 
# @param formula An object of class "formula"
# @param data A data frame created with \code{model.frame}. If another sort of object, \code{model.frame} is called first.
# @param rm.int Flag indicating if the intercept column should be removed from the model.matrix. 
#   Primarily useful for ordinal models.
# 
# @return The design matrix for a regression-like model with the specified formula and data. 
#   For details see the documentation of \code{model.matrix}.
brm.model.matrix = function(formula, data = environment(formula), rm.int = FALSE) {
  if (!is(formula, "formula")) return(NULL) 
  X <- model.matrix(formula,data)
  #cn.new <- brm.replace(colnames(X))
  cn.new <- colnames(X)
  if (rm.int & "Intercept" %in% cn.new) {
    X <- as.matrix(X[,-(1)])
    if (ncol(X)) colnames(X) <- cn.new[2:length(cn.new)]
  } 
  else colnames(X) <- cn.new
  X   
}

#  replace certain symbols in a character vector
brm.replace <- function(names, symbols = NULL, subs = NULL) {
  if (is.null(symbols))
    symbols <- c(" |\\.|\\(|\\)|\\[|\\]|,", ":", "\\+", "-", "\\*", "/", "\\^", "=", "!=")
  if (is.null(subs))
    subs <- c("", "__", "P", "M", "T", "D", "E", "EQ", "NEQ")
  if (length(symbols) != length(subs)) 
    stop("length(symbols) != length(subs)")
  for (i in 1:length(symbols)) 
    names <- gsub(symbols[i], subs[i], names)
  names
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
brm.link <- function(family) {
  link <- family[2]
  family <- family[1]
  is.lin <- family %in% c("gaussian", "student", "cauchy", "multigaussian")
  is.skew <- family %in% c("gamma", "weibull", "exponential")
  is.bin <- family %in% c("cumulative", "cratio", "sratio", "acat","binomial", "bernoulli")                    
  is.count <- family %in% c("poisson", "negbinomial", "geometric")
  if (is.na(link)) {
    if (is.lin) link <- "identity"
    else if (is.skew | is.count) link <- "log"
    else if (is.bin | family == "categorical") link <- "logit"
  }
  else if (is.lin & !is.element(link, c("identity", "log", "inverse")) |
             is.count & !link %in% c("log", "identity", "sqrt") |
             is.bin & !link %in% c("logit", "probit", "probit_approx", "cloglog") |
             family == "categorical" & link != "logit" |
             is.skew & !link %in% c("log", "identity", "inverse"))
    stop(paste(link, "is not a valid link for family", family))
  else if (is.count & link == "sqrt") 
    warning(paste(family, "model with sqrt link may not be uniquely identified"))
  link
}

# convert array to list of elements with reduced dimension
array2list <- function(x) {
  if (is.null(dim(x))) stop("Argument x has no dimension")
  n.dim <- length(dim(x))
  l <- list(length = dim(x)[n.dim])
  ind <- paste0(rep(",",n.dim-1),collapse="")
  for (i in 1:dim(x)[n.dim])
    l[[i]] <- eval(parse(text = paste0("x[",ind,i,"]")))
  names(l) <- dimnames(x)[[n.dim]]
  l
}

#calculate estimates over posterior samples 
get.estimate <- function(coef, samples, margin = 1, to.array = FALSE, ...) {
  dots <- list(...)
  args <- list(X = samples, MARGIN = margin, FUN = coef)
  fun.args <- names(formals(coef))
  if (!"..." %in% fun.args) 
    dots <- dots[fun.args %in% names(dots)] 
  x <- do.call(apply, c(args, dots))
  if (is.null(dim(x))) 
    x <- matrix(x, dimnames = list(NULL, coef))
  else if (coef == "quantile") x <- aperm(x, length(dim(x)):1)
  if (to.array & length(dim(x)) == 2) 
    x <- array(x, dim = c(dim(x), 1), dimnames = list(NULL, NULL, coef))
  x 
}

#get correlation names of random effects
get.cor.names <- function(names) {
  cor.names <- NULL
  if (length(names) > 1)
    for (i in 2:length(names)) 
      for (j in 1:(i-1)) 
        cor.names <- c(cor.names, paste0("cor(",names[j],",",names[i],")"))
  cor.names
}

isNULL <- function(x) is.null(x) | all(sapply(x, is.null))

rmNULL <- function(x) {
  x <- Filter(Negate(isNULL), x)
  lapply(x, function(x) if (is.list(x)) rmNULL(x) else x)
}

rmNum <- function(x) x[sapply(x, Negate(is.numeric))]

#removes all elements in x appearing also in y
rmMatch <- function(x, y) {
  att <- attributes(x)
  keep <- which(!(x %in% y))
  x <- x[keep]
  attributes(x) <- att
  attr(x, "match.length") <- att$match.length[keep] 
  x
} 

#melt data frame for family = "multigaussian"
brm.melt <- function(data, response, family) {
  if (length(response) > 1 & family != "multigaussian")
    stop("multivariate models are currently only allowed for family 'multigaussian'")
  else if (length(response) == 1 & family == "multigaussian")
    stop("Only one response variable detected. Use family 'gaussian' instead of 'multigaussian'")
  else if (!is(data, "data.frame"))
    stop("data must be a data.frame if family 'multigaussian' is used")
  else if (length(response) > 1 & family == "multigaussian") {
    if ("trait" %in% names(data))
      stop("trait is a resevered variable name for family 'multigaussian'")
    data <- reshape2::melt(data, measure.vars = response)
    names(data)[(ncol(data)-1):ncol(data)] <- c("trait", response[1])
  }
  data
}  

#rename parameters
#' @export
rename.pars <- function(x, ...) {
  chains <- length(x$fit@sim$samples) 
  f <- colnames(x$data$X)
  ee <- extract.effects(x$formula, family = x$family)
  r <- lapply(lapply(ee$group, function(g) get(paste0("Z_",g), x$data)), colnames)
  pars <- dimnames(x$fit)$parameters
  
  #rename fixed effects
  bs <- grepl("^b\\[", pars)
  x$fit@sim$fnames_oi[bs] <- paste0("b_",f)
  for (i in 1:chains) names(x$fit@sim$samples[[i]])[bs] <- paste0("b_",f)
  
  #rename random effects
  for (j in 1:length(r)) {
    sds <- grepl(paste0("^sd_",ee$group[[j]]), pars)
    sds_names <- paste0("sd_",ee$group[[j]],"_",r[[j]])
    cors <- grepl(paste0("^cor_",ee$group[[j]]), pars)
    cors_names <- unlist(lapply(1:length(ee$group), function(i)
      if (length(r[[i]])>1) paste0("cor_",ee$group[[i]],"_", unlist(lapply(2:length(r[[i]]), function(j) 
        lapply(1:(j-1), function(k) paste0(r[[i]][k],"_",r[[i]][j]))))))) 
    x$fit@sim$fnames_oi[sds] <- sds_names
    x$fit@sim$fnames_oi[cors] <- cors_names
    for (i in 1:chains) {
      names(x$fit@sim$samples[[i]])[sds] <- sds_names
      names(x$fit@sim$samples[[i]])[cors] <- cors_names
    }  
  }
  
  #rename residuals for family "multigaussian"
  
  
  x
}