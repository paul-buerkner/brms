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
  cn.new <- rename(colnames(X), check_dup = TRUE)
  if (rm.int && "Intercept" %in% cn.new) {
    X <- as.matrix(X[,-(1)])
    if (ncol(X)) colnames(X) <- cn.new[2:length(cn.new)]
  } 
  else colnames(X) <- cn.new
  X   
}

#  rename certain symbols in a character vector
rename <- function(names, symbols = NULL, subs = NULL, fixed = TRUE, check_dup = FALSE) {
  if (is.null(symbols))
    symbols <- c(" ", "(", ")", "[", "]", ",", "+", "-", "*", "/", "^", "=", "!=")
  if (is.null(subs))
    subs <- c(rep("", 6), "P", "M", "MU", "D", "E", "EQ", "NEQ")
  if (length(symbols) != length(subs)) 
    stop("length(symbols) != length(subs)")
  new.names <- names
  for (i in 1:length(symbols)) 
    new.names <- gsub(symbols[i], subs[i], new.names, fixed = fixed)
  dup <- duplicated(new.names)
  if (check_dup && any(dup)) 
    stop(paste0("Internal renaming of variables led to duplicated names. \n",
      "Occured for variables: ", paste(names[which(new.names %in% new.names[dup])], collapse = ", ")))
  new.names
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
    else if (is.skew || is.count) link <- "log"
    else if (is.bin || family == "categorical") link <- "logit"
  }
  else if (is.lin && !is.element(link, c("identity", "log", "inverse")) ||
           is.count && !link %in% c("log", "identity", "sqrt") ||
           is.bin && !link %in% c("logit", "probit", "probit_approx", "cloglog") ||
           family == "categorical" && link != "logit" ||
           is.skew && !link %in% c("log", "identity", "inverse"))
    stop(paste(link, "is not a valid link for family", family))
  else if (is.count && link == "sqrt") 
    warning(paste(family, "model with sqrt link may not be uniquely identified"))
  link
}

#melt data frame for family = "multigaussian"
brm.melt <- function(data, response, family) {
  if (length(response) > 1 && family != "gaussian")
    stop("multivariate models are currently only allowed for family 'gaussian'")
  #else if (length(response) == 1 && family == "multigaussian")
  #  stop("Only one response variable detected. Use family 'gaussian' instead of 'multigaussian'")
  #else if (!is(data, "data.frame") && family == "gaussian")
  #  stop("data must be a data.frame if family 'multigaussian' is used")
  else if (length(response) > 1 && family == "gaussian") {
    if (!is(data, "data.frame"))
      stop("data must be a data.frame in case of multiple responses")
    if ("trait" %in% names(data))
      stop("trait is a resevered variable name in case of multiple responses")
    data <- reshape2::melt(data, measure.vars = response)
    names(data)[(ncol(data)-1):ncol(data)] <- c("trait", response[1])
  }
  data
}  

#combine grouping factors
combine.groups <- function(data, ...) {
  group <- c(...)
  if (length(group)) {
    for (i in 1:length(group)) {
      sgroup <- unlist(strsplit(group[[i]], "__"))
      if (length(sgroup) > 1) {
        new.var <- get(sgroup[1], data)
        for (j in 2:length(sgroup)) {
          new.var <- paste0(new.var, "_", get(sgroup[j], data))
        }
        data[[group[[i]]]] <- new.var
      }
    } 
  }
  data
}

#update data for use in brm
updateData <- function(data, family, effects, ...) {
  if (!"brms.frame" %in% class(data)) {
    data <- brm.melt(data, response = effects$response, family = family)
    data <- stats::model.frame(effects$all, data = data, drop.unused.levels = TRUE)
    if (any(grepl("__", colnames(data))))
      stop("Variable names may not contain double underscores '__'")
    data <- combine.groups(data, effects$group, ...)
    class(data) <- c("brms.frame", "data.frame") 
  }
  data
}

#list irrelevant parameters not to be saved by Stan
exclude_pars <- function(formula, ranef = TRUE) {
  ee <- extract.effects(formula = formula, add.ignore = TRUE)
  out <- c("eta", "b_Intercept1")
  if (length(ee$group)) {
    out <- c(out, paste0("pre_",ee$group))
    if (!ranef) out <- c(out, paste0("r_",ee$group))
  }
  out
}