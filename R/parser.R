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
# extract_effects(response ~ I(1/a) + b)
# 
# # mixed effects model
# extract_effects(response ~ I(1/a) + b + (1 + c | d))
# 
# # mixed effects model with additional information on the response variable 
# # in this case standard errors in a gaussian linear model
# extract_effects(response | se(sei) ~ I(1/a) + b + (1 + c | d), family = "gaussian")
# }
extract_effects <- function(formula, ..., family = "none") {
  formula <- formula2string(formula)  
  fixed <- gsub(paste0("\\([^(\\||~)]*\\|[^\\)]*\\)\\+|\\+\\([^(\\||~)]*\\|[^\\)]*\\)",
                       "|\\([^(\\||~)]*\\|[^\\)]*\\)"),"",formula)
  fixed <- gsub("\\|+[^~]*~", "~", fixed)
  if (substr(fixed, nchar(fixed), nchar(fixed)) == "~") fixed <- paste0(fixed, "1")
  fixed <- formula(fixed)
  if (family %in% c("cumulative", "sratio", "cratio", "acat"))
    fixed <- update.formula(fixed, . ~ . +1)
  if (length(fixed) < 3) stop("invalid formula: response variable is missing")
  
  rg <- unlist(regmatches(formula, gregexpr("\\([^\\|\\)]*\\|[^\\)]*\\)", formula)))
  random <- lapply(regmatches(rg, gregexpr("\\([^\\|]*", rg)), function(r) 
    formula(paste0("~ ",substr(r, 2, nchar(r)))))
  cor <- unlist(lapply(regmatches(rg, gregexpr("\\|[^\\)]*", rg)), function(g) substr(g, 1, 2) != "||"))
  group_formulas <- lapply(regmatches(rg, gregexpr("\\|[^\\)]*", rg)), function(g) {
    g <- ifelse(substr(g, 1, 2) == "||", substr(g, 3, nchar(g)), substr(g, 2, nchar(g)))
    if (nchar(gsub(":", "", gsub("[^([:digit:]|[:punct:])][[:alnum:]_\\.]*", "", g))))
      stop(paste("Illegal grouping term:",g,"\nGrouping terms may contain only variable names",
                 "combined by the interaction symbol ':'"))
    return(formula(paste("~",g)))})
  group <- unlist(lapply(group_formulas, function(g) paste0(all.vars(g), collapse = ":")))
  # ordering is to ensure that all REs of the same grouping factor are next to each other
  x <- list(fixed = fixed, 
            random = if (length(group)) random[order(group)] else random, 
            cor = if (length(group)) cor[order(group)] else cor,
            group = if (length(group)) group[order(group)] else group)
  
  fun <- c("se", "weights", "trials", "cat", "cens")
  add_vars <- list()
  if (family != "none") {
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
        args <- substr(x[[f]], nchar(f) + 2, nchar(x[[f]]) -1)
        if (is.na(suppressWarnings(as.numeric(args)))) {
          x[[f]] <- as.formula(paste0("~ .", x[[f]]))
          add_vars[[f]] <- as.formula(paste("~", paste(all.vars(x[[f]]), collapse = "+")))
        }  
        else x[[f]] <- as.numeric(args)
      }  
      else stop(paste("Argument",f,"in formula is not supported by family",family))
    }
    if (nchar(gsub("\\|", "", add)) > 0 && !is.na(add))
      stop(paste0("Invalid addition part of formula. Please see the 'Details' section of help(brm) ",
                  "for further information. \nNote that the syntax of addition has changed in brms 0.2.1 as ",
                  "the old one was not flexible enough."))
  }
  
  new_formula <- unlist(lapply(c(random, group_formulas, add_vars, ...), 
                              function(x) paste0("+", Reduce(paste, deparse(x[[2]])))))
  new_formula <- paste0("update(",Reduce(paste, deparse(fixed)),", . ~ .",paste0(new_formula, collapse=""),")")
  x$all <- eval(parse(text = new_formula))
  environment(x$all) <- globalenv()
  x$response = all.vars(x$all[[2]])
  if (length(x$response) > 1) {
    if (!is.null(x$cens) || !is.null(x$se))
      stop("multivariate models currently allow only weights as addition arguments")
    x$fixed <- eval(parse(text = paste0("update(x$fixed, ", x$response[1], " ~ .)"))) 
    x$all <- eval(parse(text = paste0("update(x$all, ", x$response[1], " ~ .)"))) 
  }  
  x
} 

# extract time and grouping variabels for correlation structure
# 
# @formula a one sided formula of the form ~ time|group typically taken from a cor_brms object
# 
# @return a list with elements time, group, and all, where all contains a formula with all variables in formula
extract_time <- function(formula) {
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
    x$group <- paste0(all.vars(group), collapse = ":")
  }
  x$all <- formula(paste("~",paste(c("1", time, all.vars(group)), collapse = "+")))
  x
}

# incorporate addition and partial arguments into formula 
# 
# @param formula a model formula 
# @param addition a list with one sided formulas taken from the addition arguments in brm
# @param partial a one sided formula containing partial effects
#
# @return an updated formula containing the addition and partial effects
update_formula <- function(formula, addition = NULL, partial = NULL) {
  var_names <- names(addition)
  addition <- lapply(addition, formula2string, rm = 1)
  fnew <- "."
  if (length(addition)) 
    for (i in 1:length(addition))
      fnew <- paste0(fnew, " | ", var_names[i], "(", addition[[i]], ")")
  fnew <- paste(fnew, "~ .")
  if (is.formula(partial)) {
    partial <- formula2string(partial, rm=1)
    fnew <- paste(fnew, "+ partial(", partial, ")")
  }
  update.formula(formula, formula(fnew))
}

# find all valid object names in a string (used in method hypothesis in s3.methods.R)
# 
# @x a character string
#
# @return all valid variable names within the string
find_names <- function(x) {
  if (!is.character(x) || length(x) > 1) stop("x must be a character string of length 1")
  pos_fun <- gregexpr("([^([:digit:]|[:punct:])]|\\.|_)[[:alnum:]_\\.]*\\(", x)[[1]]
  pos_decnum <- gregexpr("\\.[[:digit:]]+", x)[[1]]
  pos_var <- list(rmMatch(gregexpr("([^([:digit:]|[:punct:])]|\\.|_)[[:alnum:]_\\.]*(\\[[[:digit:]]*\\])?", x)[[1]], 
                          pos_fun, pos_decnum))
  unlist(regmatches(x, pos_var))
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