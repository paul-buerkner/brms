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
  cn.new <- brm.replace(colnames(X))
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
    symbols <- c(" |\\.|\\(|\\)|\\[|\\]", ":", "\\+", "-", "\\*", "/", "\\^", "=", "!=")
  if (is.null(subs))
    subs <- c("", "__", "P", "M", "T", "D", "E", "EQ", "NEQ")
  if (length(symbols) != length(subs)) 
    stop("length(symbols) != length(subs)")
  for (i in 1:length(symbols)) 
    names <- gsub(symbols[i], subs[i], names)
  names
}

# Extract fixed and random effects from a formula
# 
# @param formula An object of class "formula" using mostly the syntax of the \code{lme4} package
# @param ... Additional objects of class "formula"
# 
# @return A named list of the following six objects: \cr
#   \code{fixed}:    An object of class "formula" that contains the fixed effects including the dependent variable. \cr
#   \code{random}:   A list of formulas containing the random effects per grouping variable. \cr
#   \code{group}:    A list of names of the grouping variables. \cr
#   \code{add}:      A one sided formula containing the \code{add} part of \code{formula = y | add ~ predictors} if present. \cr
#   \code{add2}:  A one sided formula containing the \code{add2} part of \code{formula = y || add2 ~ predictors} if present. \cr
#   \code{all}:      A formula that contains every variable mentioned in \code{formula} and \code{...} 
# 
# @examples
# \dontrun{ 
# # fixed effects model
# extract.effects(response ~ I(1/a) + b)
# 
# # mixed effects model
# extract.effects(response ~ I(1/a) + b + (1 + c | d))
# 
# # mixed effects model with additional information on the response variable 
# # (e.g., standard errors in a gaussian linear model)
# extract.effects(response | se:sei ~ I(1/a) + b + (1 + c | d))
# }
extract.effects <- function(formula, ..., family = "none", add.ignore = FALSE) {
  formula <- gsub(" ","",Reduce(paste, deparse(formula)))  
  fixed <- gsub(paste0("\\([^(\\||~)]*\\|[^\\)]*\\)\\+|\\+\\([^(\\||~)]*\\|[^\\)]*\\)",
                       "|\\([^(\\||~)]*\\|[^\\)]*\\)"),"",formula)
  fixed <- gsub("\\|+[^~]*~", "~", fixed)
  if (substr(fixed, nchar(fixed), nchar(fixed)) == "~") fixed <- paste0(fixed, "1")
  fixed <- formula(fixed)
  if (length(fixed) < 3) stop("invalid formula: response variable is missing")
  
  rg <- unlist(regmatches(formula, gregexpr("\\([^\\|\\)]*\\|[^\\)]*\\)", formula)))
  random <- lapply(regmatches(rg, gregexpr("\\([^\\|]*", rg)), function(r) 
    formula(paste0("~ ",substr(r, 2, nchar(r)))))
  group <- lapply(regmatches(rg, gregexpr("\\|[^\\)]*", rg)), function(g) 
    substr(g, 2, nchar(g)))
  x <- list(fixed = fixed, random = random, group = group)
  
  fun <- c("se", "weights", "trials", "cat", "cens")
  if (!add.ignore) {
    add <- unlist(regmatches(formula, gregexpr("\\|[^~]*~", formula)))[1]
    add <- substr(add, 2, nchar(add)-1)
    families <- list(se = c("gaussian","student","cauchy"), weights = c("gaussian","student","cauchy"),
      trials = c("binomial"), cat = c("categorical", "cumulative", "cratio", "sratio", "acat"), 
      cens = c("gaussian","student","cauchy","binomial","poisson","geometric","negbinomial","exponential",
               "weibull","gamma"))
    for (f in fun) {
      x[[f]] <- unlist(regmatches(add, gregexpr(paste0(f,"\\([^\\|]*\\)"), add)))[1]
      add <- gsub(paste0(f,"\\([^~|\\|]*\\)\\|*"), "", add)
      if (is.na(x[[f]])) x[[f]] <- NULL
      else if (family %in% families[[f]] | families[[f]][1] == "all") {
        x[[f]] <- substr(x[[f]], nchar(f) + 2, nchar(x[[f]]) -1)
        if (is.na(suppressWarnings(as.numeric(x[[f]])))) {
          x[[f]] <- as.formula(paste0("~", x[[f]]))
          if (length(all.vars(x[[f]])) > 1) 
            stop(paste("Argument",f,"in formula contains more than one variable"))
        }  
        else x[[f]] <- as.numeric(x[[f]])
      }  
      else stop(paste("Argument",f,"in formula is not supported by family",family))
    }
    if (nchar(gsub("\\|", "", add)) > 0 & !is.na(add))
      stop(paste0("Invalid addition part of formula. Please see the 'Details' section of help(brm) ",
        "for further information. \nNote that the syntax of addition has changed in brms 0.2.1 as ",
        "the old one was not flexible enough."))
  }
  
  if (length(group)) group <- lapply(paste("~",group),"formula") 
  up.formula <- unlist(lapply(c(random, group, rmNULL(rmNum(x[fun])), ...), 
                       function(x) paste0("+", Reduce(paste, deparse(x[[2]])))))
  up.formula <- paste0("update(",Reduce(paste, deparse(fixed)),", . ~ .",paste0(up.formula, collapse=""),")")
  all <- eval(parse(text = up.formula))
  environment(all) <- globalenv()
  return(c(x, all = all))
} 

# extract time and grouping variabels for correlation structure
extract.time <- function(formula) {
  if (is.null(formula)) return(NULL)
  formula <- gsub(" ","",Reduce(paste, deparse(formula))) 
  time <- gsub("~|\\|[[:print:]]*", "", formula)
  group <- gsub("~[^\\|]*|\\|", "", formula)
  x <- list(time = time, group = group)
  if (!nchar(group)) group <- NULL
  x$all <- formula(paste("~",paste(c(time, group), collapse = "+")))
  x
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
  is.lin <- family %in% c("gaussian", "student", "cauchy")
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
  x <- apply(samples, margin, coef, ...)
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

is.formula <- function(x, or = TRUE) {
  if (!is.list(x)) x <- list(x)
  out <- sapply(x, function(y) is(y, "formula"))
  if (or) out <- any(out)
  else out <- all(out)
  out
}
