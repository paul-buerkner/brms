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

isNULL <- function(x) is.null(x) || all(sapply(x, is.null))

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

#check if x is a whole number (integer)
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {  
  if (!is.numeric(x)) return(FALSE)
  else return(abs(x - round(x)) < tol)
}  
