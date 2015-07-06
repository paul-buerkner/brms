# Regression models in Bugs
# 
# @inheritParams brm
# @return A character string containing the model in bugs language
brm.bugs <- function(formula, data = NULL, family = "gaussian", link = "identiy", 
                       prior = list(), partial = NULL, threshold = "flexible", 
                       predict = FALSE, save.model = FALSE, ...) {  
  ef <- extract.effects(formula = formula, family = family, partial = partial)
  #data <- brm.melt(data, response = ef$response, family = family[1])
  #data <- model.frame(ef$all, data=data, drop.unused.levels=TRUE)
  
  is.ord <- family %in% c("cumulative", "cratio", "sratio", "acat")
  X = brm.model.matrix(ef$fixed, data, rm.int = is.ord)
  f = colnames(X)
  r = sapply(lapply(ef$random, brm.model.matrix, data = data, rm.int = is.ord), colnames)
  n = ifelse(is(ef$trials,"formula") | is(ef$cat,"formula"), "[n]", "")
  ilink <- c(identity = "", log = "exp", inverse = "inv", sqrt = "square", logit = "ilogit", 
             probit = "phi", probit.approx = "Phi_approx", cloglog = "icloglog")[[link]]
  llh <- bugs.llh(family, link = link, se = is(ef$se, "formula"))
  fe.only <- setdiff(f, unlist(lapply(r, intersect, y = f)))
  if (length(f)) eta.fe <- unlist(sapply(eval(1:length(f)), function(x)
    if (is.element(f[x],fe.only)) paste0("b_",f[x],"*X[n,",x,"]")))
  else eta.fe <- NULL
  eta.re <- unlist(lapply(ef$group, function(g)
    if (length(r[[match(g, ef$group)]])) paste0("inprod(r_",g,"[",g,"[n],],Z_",g,"[n,])")))
  eta <- paste0(c(eta.fe, eta.re), collapse=" + ")
  if (!length(f) & all(!sapply(r,length))) eta <- "0"

  ord.lh.jags <- function(family) {
    th = function(k,r=ef$random,g=ef$group) {
      out <- ""
      if (length(g)) {
        r <- lapply(r,function(r) colnames(brm.model.matrix(r,data)))
        out <- paste0(unlist(sapply(eval(1:length(g)), function(i)
          if (is.element("Intercept",r[[i]])) paste0("rI_",g[[i]],"[",g[[i]],"[n],",k,"]"))),collapse=" + ")
      }  
      if (out=="") out <- paste0("b_Intercept[",k,"]")
      out
    } 
    ord <- ""
    if (family=="cumulative")
      ord <- paste0(
        "  p[n,1] <- ",ilink,"(",th(1)," - eta[n]) \n",
        "  for ( k in 2:(max_obs",n,"-1)) { \n", 
        "    p[n,k] <- ",ilink,"(",th("k")," - eta[n]) - ",ilink,"(",th("k-1")," - eta[n]) \n", 
        "  } \n",  
        "  p[n,max_obs",n,"] <- 1 - ",ilink,"(",th(paste0("max_obs",n,"-1"))," - eta[n]) \n")
    else if (family=="cratio") 
      ord <- paste0(
        "  q[n,1] <- ",ilink,"(eta[n] - ",th(1),") \n",
        "  p[n,1] <- 1-q[n,1] \n", 
        "  for ( k in 2:(max_obs",n,"-1)) { \n",
        "    q[n,k] <- ",ilink,"(eta[n] - ",th("k"),") \n",
        "    p[n,k] <- (1-q[n,k]) * prod(q[n,1:(k-1)]) \n",
        "  } \n",
        "  p[n,max_obs",n,"] <- prod(q[n,1:(max_obs",n,"-1)]) \n") 
    else if (family=="sratio")
      ord <- paste0(
        "  q[n,1] <- 1-",ilink,"(",th(1)," - eta[n]) \n",
        "  p[n,1] <- 1-q[n,1] \n", 
        "  for ( k in 2:(max_obs",n,"-1)) { \n",
        "    q[n,k] <- 1-",ilink,"(",th("k")," - eta[n]) \n",
        "    p[n,k] <- (1-q[n,k]) * prod(q[n,1:(k-1)]) \n",
        "  } \n",
        "  p[n,max_obs",n,"] <- prod(q[n,1:(max_obs",n,"-1)]) \n") 
    else if (family=="acat") {
      if (link=="logit") 
        ord <- paste0(
          "  q[n,1] <- 1 \n",
          "  for ( k in 1:(max_obs",n,"-1)) { \n",
          "    qq[n,k] <- eta[n] - ",th("k")," \n",
          "    q[n,k+1] <- exp(sum(qq[n,1:k])) \n",
          "  } \n",
          "  p[n,1:max_obs",n,"] <- q[n,1:max_obs",n,"]/sum(q[n,1:max_obs",n,"]) \n")
      else
        ord <- paste0(
          "  qq[n,1] <- ",ilink,"(eta[n] -",th(1),") \n",
          "  q[n,1] <- prod(1-qq[n,1:(max_obs",n,"-1)]) \n",
          "  for ( k in 2:(max_obs",n,"-1)) { \n",
          "    qq[n,k] <- ",ilink,"(eta[n] - ",th("k"),") \n",
          "    q[n,k] <- prod(qq[n,1:(k-1)]) * prod(1-qq[n,k:(max_obs",n,"-1)]) \n",
          "  } \n",
          "  q[n,max_obs",n,"] <- prod(qq[n,1:(max_obs",n,"-1)]) \n",
          "  p[n,1:max_obs",n,"] <- q[n,1:max_obs",n,"]/sum(q[n,1:max_obs",n,"]) \n")
    }
  }

  b.priors <- unlist(lapply(f, function(x) paste0("b_",x)))
  prior[c(b.priors, "sigma", "shape")] <- bugs.prior(c(b.priors, "sigma", "shape"), 
                                            prior = prior, engine = "jags", s = 0) 
  fe.prior.jags <- paste0(prior[b.priors], collapse = "")
  if (family=="gaussian" & !is(ef$se,"formula")) 
    fe.prior.jags <- paste0(fe.prior.jags, prior["sigma"])
  if (is.element(family,c("gamma","weibull"))) 
    fe.prior.jags <- paste0(fe.prior.jags, prior["shape"])
  
  re.prior.jags <- function(rgroup) {
    random <- rgroup[[1]]
    group <- rgroup[[2]]
    r <- colnames(brm.model.matrix(random, data, rm.int = is.ord))
    if (!length(r)) return("")   
    prior[paste0("P_", group)] <- bugs.prior(paste0("P_", group), prior = prior, s = 0,
                                            df = length(r), engine = "jags")
    rprior = ifelse(length(r) > 1, " ~ dmnorm"," ~ dnorm")
    out = paste0("for ( g in 1:N_",group,") { \n",
                 "  r_",group,"[g,1:K_",group,"]", rprior, "(M_",group,",P_",group,") \n", 
                 "} \n", 
                 prior[paste0("P_",group)], "V_",group," <- inverse(P_",group,") \n")  
    for ( z in 1:length(r)) {
      if ( is.element(r[z],f)) {
        out <- paste0(out,"M_",group,"[",z,"] <- b_",r[z],"\n") 
      }
      else out <- paste0(out,"M_",group,"[",z,"] <- 0 \n")
    }
    out
  }
  re.prior.jags <- paste0(sapply(mapply(c, ef$random, ef$group, SIMPLIFY=FALSE), re.prior.jags),
                          collapse = "")
  
  ord.fe.prior.jags <- c(bugs.prior("b_Intercept", prior = prior, engine = "jags", ind = "k"), "")
  if (threshold == "equidistant" ) { 
    ord.fe.prior.jags <- c(paste0("  b_Intercept[k] <- b_Intercept1 + (k-1)*delta \n"),
                           paste0(bugs.prior(c("b_Intercept1","delta"), prior=prior,engine="jags",s=0),collapse=""))
  }
  ord.re.prior.jags <- function(rgroup) {
    random <- rgroup[[1]]
    group <- rgroup[[2]]
    r = colnames(brm.model.matrix(random,data))
    if ( !is.element("Intercept",r)) return(rep("",2))
    prior[paste0("PI_",group)] = bugs.prior(paste0("PI_",group), prior=prior, engine="jags",s=0)                                 
    out=rep("",2)
    out[1] <- paste0(out[1], "  for ( g in 1:N_",group,") { \n",
                     "    rI_",group,"[g,k] ~ dnorm(b_Intercept[k],PI_",group,") \n", "  } \n")
    out[2] <- paste0(out[2],prior[paste0("PI_",group)], "VI_",group," <- inverse(PI_",group,") \n")  
    return(out)  
  }
  ord.re.prior.jags <- lapply(mapply(c, ef$random, ef$group, SIMPLIFY=FALSE), ord.re.prior.jags)
  ord.re.prior.jags <- c(paste0("",sapply(ord.re.prior.jags, function(x) return(x[1])),collapse=""),
                         paste0("",sapply(ord.re.prior.jags, function(x) return(x[2])),collapse=""))
  
  model <- paste0( 
    "model { \n",
    "### likelihood \n",
    "for ( n in 1:N) { \n",
    "  Y[n] ~ ",llh, "\n",
    if (predict) 
      paste0("  Y_pred[n] ~ ",llh, "\n"),
    "  eta[n] <- ",eta," \n",
       ord.lh.jags(family), 
    "} \n \n",
    "### prior for fixed effects \n",
    fe.prior.jags, "\n",
    "### prior for random effects \n",
    re.prior.jags, "\n",
    "### prior for parameters specific to ordinal models \n",
    if (is.ord) paste0(
      "for ( k in 1:(max(max_obs)-1)) { \n",
        ord.fe.prior.jags[1], ord.re.prior.jags[1], "} \n",
        ord.fe.prior.jags[2], ord.re.prior.jags[2]), "\n",
    "}")

  class(model) <- c("character", "brmsmodel")
  if (is.character(save.model)) {
    sink(save.model)
    cat(model)
    sink()
  }
  model
}


# Priors in Bugs
bugs.prior = function(par, prior = list(), ind = rep("", length(par)), s = 0, ...) { 
  if (length(par) != length(ind)) 
    stop("The lengths of par and ind must be the same")
  fun <- function(par_ind) {
    par <- par_ind[[1]]
    ind <- par_ind[[2]]
    if (ind != "") ind <- paste0("[",ind,"]")
    par.type <- unlist(regmatches(par, gregexpr("[^_]*", par)))[1]
    if (is.na(par.type)) par.type <- ""
    if (is.null(prior[[par]])) {
      if (!is.null(prior[[par.type]])) prior[[par]] <- prior[[par.type]]   
      else {
        df <- ifelse(is.null(list(...)$df), 1, list(...)$df)
        if (par.type == "b" | par == "delta") prior[[par]] <- "dnorm(0,1e-06)"
        else if (par.type == "PI") prior[[par]] <- "dgamma(0.001,0.001)"
        else if (par == "sigma") prior[[par]] <- "dunif(0,100)"
        else if (par == "shape") prior[[par]] <- "dgamma(1,0.001)"
        else if (par.type == "P")
          if (df == 1) prior[[par]] <- "dgamma(0.001,0.001)"
        else prior[[par]] <- paste0("dwish(Sigma_",substr(par,3,nchar(par)),",",df,")")
        else return("")
      }
      prior[[par]] <- paste("~", prior[[par]])
    }
    else {
      p <- gsub(" ","",prior[[par]])
      if (!substr(p,1,1) == "~" & !substr(p,1,2) == "<-") 
        prior[[par]] <- paste("~",prior[[par]])
    }  
    return(paste0(paste0(rep(" ", s), collapse = ""), par, ind, " ", prior[[par]], " \n")) 
  }
  return(sapply(mapply(list, par, ind, SIMPLIFY = FALSE), fun))
}

# Likelihoods in stan language
# 
# Define the likelihood of the dependent variable in stan language
# 
# @inheritParams brm
# @param add A flag inicating if the model contains additional information of the response variable
#   (e.g., standard errors in a gaussian linear model)
# @param weights A flag indicating if the response variable should have unequal weights 
#   Only used if \code{family} is either \code{"gaussian", "student"}, or \code{"cauchy"}.
#    
# @return A character string defining a line of code in stan language 
#   that contains the likelihood of the dependent variable. 
# @examples 
# \dontrun{
# brm.llh(family = "gaussian")
# brm.llh(family = "cumulative", link = "logit")
# }
bugs.llh <- function(family, link = "identity", predict = FALSE, se = FALSE,
                     weights = FALSE, cens = FALSE) {
  is.ord <- is.element(family, c("cumulative", "cratio", "sratio", "acat"))
  n <- ifelse(predict | is.element(link, c("inverse","sqrt")), "[n]", "")
  ilink <- c(identity = "", log = "exp", inverse = "inv", sqrt = "square", logit = "ilogit", 
             probit = "phi", probit_approx = "Phi_approx", cloglog = "icloglog")[[link]]
  llh <- list(gaussian = paste0("dnorm(eta[n],pow(sigma",if (se)"[n]",",-2))"), 
              poisson = "dpois(exp(eta[n]))", 
              binomial = paste0("dbin(",ilink,"(eta[n]), max_obs",n,")"), 
              gamma = "dgamma(shape, shape/exp(eta[n]))",
              weibull = paste0("dweib(shape, exp(-eta[n]/shape))"),
              exponential = paste0("dexp(exp(-eta[n]))"),
              ordinal = paste0("dcat(p[n,1:max_obs",n,"])"))[[ifelse(is.ord,"ordinal",family)]]
  llh
}

# Initial values for Bugs models in \code{brms}
# 
# Generate initial values for Bugs models in \code{brms}
# 
# @param range A positive number. Defines the interval around zero (from \code{-range} to \code{range}) 
#   from which the regression parameters are generated randomly (default is 10).
# @inheritParams brm
# 
# @return A named list that specifies the initial values for parameters
# @examples 
# bugs.inits(value~treat+period+carry+(1|subject), data = inhaler, family = "gaussian")
# 
# bugs.inits(value~treat+period+carry+(1|subject), data = inhaler, family = "sratio",
#           threshold = "equidistant")
bugs.inits <- function(formula, data = NULL, range = 2, family="gaussian", partial = NULL, 
                       threshold = "flexible", engine="stan", seed=NULL) {
  if (is.null(range)) range <- 2
  ef <- extract.effects(formula = formula, family = family, partial = partial) 
  #data <- brm.melt(data, response = ef$response, family = family[1])
  #data <- model.frame(ef$all, data = data, drop.unused.levels = TRUE)
  
  if(length(as.integer(seed)) == 1) set.seed(as.integer(seed))
  is.lin <- is.element(family, c("gaussian", "student", "cauchy"))
  is.ord <- is.element(family,c("cumulative","cratio","sratio","acat"))
  if (is.ord | family == "binomial") {
    if (family == "binomial") add <- ef$trials
    else add <- ef$cat
    if (!length(add)) max_obs <- max(as.numeric(model.response(data)))
    else if (is.numeric(add)) max_obs <- add
    else if (is(add, "formula")) 
      max_obs <- brm.model.matrix(add, data, rm.int = TRUE)
    else stop("response part of argument formula is invalid")
  }
  f <- colnames(brm.model.matrix(ef$fixed,data, rm.int = is.ord))
  r <- lapply(ef$random, function(r,data) colnames(brm.model.matrix(r, data)), data = data)
  rb <- sapply(r,function(r) length(setdiff(r, ifelse(is.ord, "Intercept", NA))))
  rth <- unlist(lapply(r, function(r) is.element("Intercept",r)))
  p <- colnames(brm.model.matrix(partial,data, rm.int = TRUE))
  
  inits = list()
  if (length(f)) inits <- c(inits, setNames(as.list(runif(length(f), -range, range)), paste0("b_",f)))
  if (is.lin & !is(ef$se,"formula")) inits <- c(inits,sigma=runif(1,2,10))
  if (family == "student") inits <- c(inits, nu = runif(1,1,60))
  else if (is.element(family,c("gamma","weibull"))) inits <- c(inits, shape = runif(1,0.2,5))
  if (length(ef$random))
    for ( i in (1:length(ef$random))[rb>0]) {
      lg <- length(unique(get(ef$group[[i]],data)))
      mat <- matrix(runif(lg*rb[i],-range,range),nrow=lg,ncol=rb[i])
      Pmat <- runif(rb[i],0.1,0.5)
      if (rb[i]>1) Pmat <- diag(Pmat)
      inits <- c(inits, setNames(list(mat), paste0("r_",ef$group[[i]])))
      inits <- c(inits, setNames(list(Pmat), paste0("P_",ef$group[[i]])))
    }
  if (is.ord) {
    if (threshold == "flexible") 
      inits <- c(inits, list(b_Intercept = sort(runif(max(max_obs)-1,-range,range))))
    else if (threshold == "equidistant") 
      inits <- c(inits, list(b_Intercept1 = runif(1, -range, range), delta=runif(1,0.2,1)))
    if (length(p)) 
      inits <- c(inits, setNames(replicate(length(p), runif(max(max_obs)-1, -range, range),
                                           simplify = FALSE), paste0("b_",p)))
    for (i in (1:length(ef$random))[rth]) {
      lg <- length(unique(get(ef$group[[i]],data)))
      mat <- t(apply(matrix(runif(lg*(max(max_obs) - 1), -range, range),
                            nrow = lg, ncol = max(max_obs) - 1), 1, sort))
      inits <- c(inits, setNames(list(mat, runif(1, 0.1, 10)),
                                 paste0(c("rI_", "PI_"), ef$group[[i]])))
    } 
  }
  return(inits)
}
