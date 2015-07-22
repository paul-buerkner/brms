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
# # in this case standard errors in a gaussian linear model
# extract.effects(response | se(sei) ~ I(1/a) + b + (1 + c | d), family = "gaussian")
# }
extract.effects <- function(formula, ..., family = "none", add.ignore = FALSE) {
  formula <- formula2string(formula)  
  fixed <- gsub(paste0("\\([^(\\||~)]*\\|[^\\)]*\\)\\+|\\+\\([^(\\||~)]*\\|[^\\)]*\\)",
                       "|\\([^(\\||~)]*\\|[^\\)]*\\)"),"",formula)
  fixed <- gsub("\\|+[^~]*~", "~", fixed)
  if (substr(fixed, nchar(fixed), nchar(fixed)) == "~") fixed <- paste0(fixed, "1")
  fixed <- formula(fixed)
  if (length(fixed) < 3) stop("invalid formula: response variable is missing")
  
  rg <- unlist(regmatches(formula, gregexpr("\\([^\\|\\)]*\\|[^\\)]*\\)", formula)))
  random <- lapply(regmatches(rg, gregexpr("\\([^\\|]*", rg)), function(r) 
    formula(paste0("~ ",substr(r, 2, nchar(r)))))
  cor <- lapply(regmatches(rg, gregexpr("\\|[^\\)]*", rg)), function(g) cor <- substr(g, 1, 2) != "||")
  group <- lapply(regmatches(rg, gregexpr("\\|[^\\)]*", rg)), function(g) {
    g <- ifelse(substr(g, 1, 2) == "||", substr(g, 3, nchar(g)), substr(g, 2, nchar(g)))
    if (nchar(gsub(":", "", gsub("[^([:digit:]|[:punct:])][[:alnum:]_\\.]*", "", g))))
      stop(paste("Illegal grouping term:",g,"\nGrouping terms may contain only variable names",
                 "combined by the interaction symbol ':'"))
    return(formula(paste("~",g)))})
  x <- list(fixed = fixed, random = random, cor = cor,
            group = lapply(group, function(g) paste0(all.vars(g), collapse = "__")))
  
  fun <- c("se", "weights", "trials", "cat", "cens")
  if (!add.ignore) {
    add <- unlist(regmatches(formula, gregexpr("\\|[^~]*~", formula)))[1]
    add <- substr(add, 2, nchar(add)-1)
    families <- list(se = c("gaussian","student","cauchy"), weights = c("all"),
                     trials = c("binomial"), cat = c("categorical", "cumulative", "cratio", "sratio", "acat"), 
                     cens = c("gaussian","student","cauchy","binomial","poisson","geometric","negbinomial","exponential",
                              "weibull","gamma"))
    for (f in fun) {
      x[[f]] <- unlist(regmatches(add, gregexpr(paste0(f,"\\([^\\|]*\\)"), add)))[1]
      add <- gsub(paste0(f,"\\([^~|\\|]*\\)\\|*"), "", add)
      if (is.na(x[[f]])) x[[f]] <- NULL
      else if (family %in% families[[f]] || families[[f]][1] == "all") {
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
    if (nchar(gsub("\\|", "", add)) > 0 && !is.na(add))
      stop(paste0("Invalid addition part of formula. Please see the 'Details' section of help(brm) ",
                  "for further information. \nNote that the syntax of addition has changed in brms 0.2.1 as ",
                  "the old one was not flexible enough."))
  }
  
  up.formula <- unlist(lapply(c(random, group, rmNULL(rmNum(x[fun])), ...), 
                              function(x) paste0("+", Reduce(paste, deparse(x[[2]])))))
  up.formula <- paste0("update(",Reduce(paste, deparse(fixed)),", . ~ .",paste0(up.formula, collapse=""),")")
  x$all <- eval(parse(text = up.formula))
  environment(x$all) <- globalenv()
  x$response = all.vars(x$all[[2]])
  if (length(x$response) > 1) {
    x$fixed <- eval(parse(text = paste0("update(x$fixed, ", x$response[1], " ~ .)"))) 
    x$all <- eval(parse(text = paste0("update(x$all, ", x$response[1], " ~ .)"))) 
  }  
  x
} 

# extract time and grouping variabels for correlation structure
extract.time <- function(formula) {
  if (is.null(formula)) return(NULL)
  formula <- gsub(" ","",Reduce(paste, deparse(formula))) 
  time <- all.vars(as.formula(paste("~", gsub("~|\\|[[:print:]]*", "", formula))))
  if (length(time) > 1) stop("Autocorrelation structures may only contain 1 time variable")
  x <- list(time = ifelse(length(time), time, ""), group = "")

  group <- gsub("~[^\\|]*|\\|", "", formula)
  if (nchar(group)) {
    if (nchar(gsub(":", "", gsub("[^([:digit:]|[:punct:])][[:alnum:]_\\.]*", "", group))))
      stop(paste("Illegal grouping term:",group,"\nGrouping terms may contain only variable names",
                 "combined by the interaction symbol ':'"))
    group <- formula(paste("~", group))
    x$group <- paste0(all.vars(group), collapse = "__")
  }
  x$all <- formula(paste("~",paste(c("1", time, all.vars(group)), collapse = "+")))
  x
}

# incorporate addition and partial arguments into formula 
brm.update.formula <- function(formula, addition = NULL, partial = NULL) {
  var.name <- names(addition)
  addition <- lapply(addition, formula2string, rm = 1)
  fnew <- "."
  if (length(addition)) 
    for (i in 1:length(addition))
      fnew <- paste0(fnew, " | ", var.name[i], "(", addition[[i]], ")")
  fnew <- paste(fnew, "~ .")
  if (is.formula(partial)) {
    partial <- formula2string(partial, rm=1)
    fnew <- paste(fnew, "+ partial(", partial, ")")
  }
  update.formula(formula, formula(fnew))
}

#find all valid object names in a string
find.names <- function(x) {
  if (!is.character(x)) stop("x must be of class character")
  fun.pos <- gregexpr("([^([:digit:]|[:punct:])]|\\.|_)[[:alnum:]_\\.]*\\(", x)[[1]]
  decnum.pos <- gregexpr("\\.[[:digit:]]+", x)[[1]]
  var.pos <- list(rmMatch(gregexpr("([^([:digit:]|[:punct:])]|\\.|_)[[:alnum:]_\\.]*(\\[[[:digit:]]*\\])?", x)[[1]], 
                          fun.pos, decnum.pos))
  unlist(regmatches(x, var.pos))
}

#checks if x is formula (or list of formulas)
is.formula <- function(x, or = TRUE) {
  if (!is.list(x)) x <- list(x)
  out <- sapply(x, function(y) is(y, "formula"))
  if (or) out <- any(out)
  else out <- all(out)
  out
}

#converts formula to string
formula2string <- function(formula, rm = c(0,0)) {
  if (!is.formula(formula))
    stop(paste(deparse(substitute(formula)),"must be of class formula"))
  if (is.na(rm[2])) rm[2] <- 0
  x <- gsub(" ","", Reduce(paste, deparse(formula)))
  x <- substr(x, 1 + rm[1], nchar(x)-rm[2])
  x
} 