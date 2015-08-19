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
  else if (link == "probit") logit(0.07056*x^3 + 1.5976*x)
  else if (link == "cloglog") log(-log(1-x))
  else if (link == "probit_approx") logit(0.07056*x^3 + 1.5976*x)
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
