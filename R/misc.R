# check if an object is NULL
isNULL <- function(x) is.null(x) || ifelse(is.vector(x), all(sapply(x, is.null)), FALSE)

# recursively removes NULL entries from an object
rmNULL <- function(x) {
  x <- Filter(Negate(isNULL), x)
  lapply(x, function(x) if (is.list(x)) rmNULL(x) else x)
}

# remove all numeric elements from an object
rmNum <- function(x) x[sapply(x, Negate(is.numeric))]

# remove all elements in x that also appear in ... while keeping all attributes
rmMatch <- function(x, ...) {
  att <- attributes(x)
  keep <- which(!(x %in% c(...)))
  x <- x[keep]
  attributes(x) <- att
  attr(x, "match.length") <- att$match.length[keep] 
  x
} 

# check if x is a whole number (integer)
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {  
  if (!is.numeric(x)) return(FALSE)
  else return(abs(x - round(x)) < tol)
}  

# wrapper for paste with collapse
collapse <- function(..., sep = "")
  paste(..., sep = sep, collapse = "")

# compute the logit
logit <- function(p) {
  log(p/(1-p))
}

# compute the inverse of logit
ilogit <- function(x) { 
  exp(x) / (1 + exp(x))
}

#checks if x is formula (or list of formulas)
#
# @param x s formula or a list of formulas
# @param or logical; indicates if any element must be a formula (or = TRUE) or if all elements must be formulas
# 
# @return TRUE or FALSE 
is.formula <- function(x, or = TRUE) {
  if (!is.list(x)) x <- list(x)
  out <- sapply(x, function(y) is(y, "formula"))
  if (or) out <- any(out)
  else out <- all(out)
  out
}

# converts formula to string
#
# @param formula a model formula
# @param rm a vector of to elements indicating how many characters should be removed at the beginning
#  and end of the string respectively
#
# @return the formula as string 
formula2string <- function(formula, rm = c(0, 0)) {
  if (!is.formula(formula))
    stop(paste(deparse(substitute(formula)),"must be of class formula"))
  if (is.na(rm[2])) rm[2] <- 0
  x <- gsub(" ","", Reduce(paste, deparse(formula)))
  x <- substr(x, 1 + rm[1], nchar(x)-rm[2])
  x
} 