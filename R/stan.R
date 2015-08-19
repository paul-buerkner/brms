# Regression models in Stan
# 
# @inheritParams brm
# @return A character string containing the model in stan language
# @examples
# \dontrun{
# stan.model <- brm.stan(y ~ log_Age + log_Base4 * Trt, data = epilepsy, 
#                   family = "poisson", link = "log", prior = list(b_ = "normal(0,5)"))
# }
stan.model <- function(formula, data = NULL, family = "gaussian", link = "identity",
                       prior = list(), partial = NULL, threshold = "flexible", cov.ranef = NULL,
                       WAIC = FALSE, predict = FALSE, sample.prior = FALSE, autocor = cor.arma(), 
                       save.model = NULL) {
  ee <- extract.effects(formula = formula, family = family, partial = partial) 
  family <- ifelse(family == "gaussian" && length(ee$response) > 1, "multinormal", family) 
  is.linear <- family %in% c("gaussian", "student", "cauchy")
  is.ordinal <- family %in% c("cumulative", "cratio", "sratio", "acat") 
  is.skew <- family %in% c("gamma", "weibull", "exponential")
  is.count <- family %in% c("poisson", "negbinomial", "geometric")
  is.multi <- family == "multinormal"

  if (family == "categorical") {
    X <- data.frame()
    Xp <- brm.model.matrix(ee$fixed, data, rm.int = is.ordinal)
  }
  else {
    X <- brm.model.matrix(ee$fixed, data, rm.int = is.ordinal)
    Xp <- brm.model.matrix(partial, data, rm.int = TRUE)
  }  
  f <- colnames(X)
  p <- colnames(Xp)
  Z <- lapply(ee$random, brm.model.matrix, data = data)
  r <- lapply(Z,colnames)
  n <- ifelse(is.formula(ee[c("trials","cat")]), "[n]", "")
  trait <- ifelse(is.multi, "_trait", "")
  
  ranef <- unlist(lapply(mapply(list, r, ee$group, ee$cor, SIMPLIFY = FALSE), stan.ranef, 
                         prior = prior, cov.ranef = cov.ranef))
  names.ranef <- unique(names(ranef))
  if (length(ranef)) ranef <- sapply(1:length(names.ranef), function(x) 
    collapse(ranef[seq(x, length(ranef), length(names.ranef))]))
  ranef <- setNames(as.list(ranef), names.ranef)
  
  max_obs <- ifelse(rep(is.formula(ee[c("trials", "cat")]), 3), 
    c("MAX_obs", "  int MAX_obs; \n", "  MAX_obs <- max(max_obs); \n"), c("max_obs", rep("", 2)))
  eta <- stan.eta(family = family, link = link, f = f, p = p, group = ee$group, 
                  autocor = autocor, max_obs = max_obs)
  ma <- stan.ma(family = family, link = link, autocor = autocor, group = ee$group, N = nrow(data),
                levels = unlist(lapply(ee$group, function(g) length(unique(get(g, data))))))
  ordinal <- stan.ordinal(family = family, link = link, partial = length(p), max_obs = max_obs[1], n = n, 
                  threshold = threshold, simplify = !(predict || WAIC))  
  multi <- stan.multi(family, response = ee$response)
  llh <- stan.llh(family, link = link, add = is.formula(ee[c("se", "trials", "cat")]), 
                  weights = is.formula(ee$weights), cens = is.formula(ee$cens))
  genquant <- stan.genquant(family, link = link, add = is.formula(ee[c("se", "trials", "cat")]),
                            predict = predict, logllh = WAIC)
  priors <- paste0(
    if (length(f)) stan.prior(paste0("b_",f), prior, ind = 1:length(f)),
    if (length(p)) stan.prior(paste0("b_",p), prior, ind = 1:length(p), partial = TRUE), 
    if (autocor$p) stan.prior("ar", prior),
    if (autocor$q) stan.prior("ma", prior),
    if (is.ordinal && threshold == "flexible") 
      stan.prior("b_Intercept", prior, add.type = "Intercept")
    else if (is.ordinal && threshold == "equidistant") 
      paste0(stan.prior("b_Intercept1", prior, add.type = "Intercept1"), stan.prior("delta", prior)),
    if (family %in% c("gamma", "weibull", "negbinomial")) stan.prior("shape", prior),
    if (family == "student") stan.prior("nu", prior),
    if (is.linear && !is.formula(ee$se)) stan.prior(paste0("sigma_",ee$response), prior), 
    if (is.multi) paste0(stan.prior(paste0("sigma_",ee$response), prior, ind = 1:length(ee$response)), 
                      stan.prior("Lrescor", prior)),
    ranef$model)
  rngprior <- stan.rngprior(sample.prior = sample.prior, priors = priors, family = family, fixed = f,
                            partial = p, response = ee$response, random = r, group = ee$group)
  
  kronecker <- any(sapply(mapply(list, r, ee$group, SIMPLIFY = FALSE), function(x, cov.ranef = "")
    length(x[[1]]) > 1 && x[[2]] %in% cov.ranef, cov.ranef = cov.ranef))
  functions.text <- stan.function(kronecker = kronecker)
  
  data.text <- paste0(
    "data { \n",
    "  int<lower=1> N; \n", 
    if (is.linear || is.skew) 
      "  real Y[N]; \n"
    else if (is.count || is.ordinal || family %in% c("binomial", "bernoulli", "categorical")) 
      "  int Y[N]; \n"
    else if (is.multi) paste0(
      "  int<lower=1> N_trait; \n  int<lower=1> K_trait; \n",  
      "  int NC_trait; \n  vector[K_trait] Y[N_trait]; \n"),
    if (ncol(X)) "  int<lower=1> K; \n  matrix[N,K] X; \n",
    if (length(p)) "  int<lower=1> Kp; \n  matrix[N,Kp] Xp; \n",  
    if (autocor$p && is(autocor, "cor.arma")) 
      "  int<lower=1> Kar; \n  matrix[N,Kar] Yar; \n",
    if (autocor$q && is(autocor, "cor.arma")) 
      "  int<lower=1> Kma; \n  row_vector[Kma] Ema_pre[N]; \n  vector[N] tgroup; \n",
    if (is.linear && is.formula(ee$se))
      "  real<lower=0> sigma[N]; \n",
    if (is.formula(ee$weights))
      paste0("  vector<lower=0>[N",trait,"] weights; \n"),
    if (is.ordinal || family %in% c("binomial", "categorical")) 
      paste0("  int max_obs",toupper(n),"; \n"),
    if (is.formula(ee$cens) && !(is.ordinal || family == "categorical"))
      "  vector[N] cens; \n",
    ranef$data,
    "} \n")
  
  zero <- ifelse(rep(family == "categorical", 2), c("  row_vector[1] zero; \n", "  zero[1] <- 0; \n"), "")
  trans.data.text <- paste0(
    "transformed data { \n",
      max_obs[2], zero[1], ranef$tdataD, 
      max_obs[3], zero[2], ranef$tdataC,
    "} \n")
  
  par.text <- paste0(
    "parameters { \n",
    if (length(f)) "  vector[K] b; \n",
    if (length(p)) paste0("  matrix[Kp,",max_obs[1],"-1] bp; \n"),
    ordinal$par, ranef$par,
    if (autocor$p && is(autocor, "cor.arma")) "  vector[Kar] ar; \n",
    if (autocor$q && is(autocor, "cor.arma")) "  vector[Kma] ma; \n",
    if (is.linear && !is.formula(ee$se)) 
      "  real<lower=0> sigma; \n",
    if (family == "student") 
      "  real<lower=0> nu; \n",
    if (is.multi) 
      "  vector<lower=0>[K_trait] sigma; \n  cholesky_factor_corr[K_trait] Lrescor; \n",
    if (family %in% c("gamma", "weibull", "negbinomial")) 
      "  real<lower=0> shape; \n",
    rngprior$par,
    "} \n")
  
  vectorize <- c(!(length(ee$group) || autocor$q || eta$transform ||
    (is.ordinal && !(family == "cumulative" && link == "logit" && !predict && !WAIC && !is.formula(ee$cat)))),            
    !(is.formula(ee$cens) || is.formula(ee$weights) || is.ordinal || family == "categorical")) 
  if (!vectorize[1] && !is.multi)
    loop.trans <- c("  for (n in 1:N) { \n", "  } \n")
  else if (is.multi)
    loop.trans <- c(paste0("  for (m in 1:N_trait) { \n  for (k in 1:K_trait) { \n",    
                           "    int n; \n    n <- (k-1)*N_trait + m; \n"), "  }} \n")
  else loop.trans <- rep("", 2)

  model <- paste0(functions.text,
  data.text, trans.data.text, par.text,
  "transformed parameters { \n",
    eta$transD, ma$transD, ordinal$transD, ranef$transD, 
    eta$transC1, ma$transC1, ordinal$transC1, ranef$transC, 
    loop.trans[1],
      eta$transC2, ma$transC2, ordinal$transC2, eta$transC3, 
    loop.trans[2],
  "} \n",
  "model { \n",
    if (is.formula(ee$weights) && !is.formula(ee$cens)) 
      paste0("  vector[N",trait,"] lp_pre; \n"),
    priors, 
    ifelse(vectorize[2], llh, paste0("  for (n in 1:N",trait,") { \n  ", llh,"  } \n")),
    if (is.formula(ee$weights) && !is.formula(ee$cens)) 
    "  increment_log_prob(dot_product(weights,lp_pre)); \n",
    rngprior$model,
  "} \n",
  "generated quantities { \n",
    multi$genD, ranef$genD, genquant$genD, rngprior$genD, 
    multi$genC, ranef$genC, genquant$genC, rngprior$genC,
  "} \n")
  
  class(model) <- c("character", "brmsmodel")
  if (is.character(save.model)) {
    sink(save.model)
    cat(model)
    sink()
  }
  model
}

# Random effects in Stan 
# 
# @return A vector of strings containing the random effects in stan language
stan.ranef <- function(rg, prior = list(), cov.ranef = "") {
  r <- rg[[1]]
  g <- rg[[2]]
  cor <- rg[[3]]
  c.cov <- g %in% cov.ranef
  out <- setNames(as.list(rep("", 9)), c("data", "tdataD", "tdataC", "par", "model", "tranD", "transC", "genD", "genC"))
  out$data <- paste0("  int<lower=1> ",g,"[N]; \n",
                     "  int<lower=1> N_",g,"; \n",
                     "  int<lower=1> K_",g,"; \n",
                     if (c.cov && (cor || length(r) == 1)) 
                       paste0("  matrix[N_",g,", N_",g,"] cov_",g,"; \n"),
                     if (c.cov && !cor && length(r) > 1) 
                       paste0("  matrix[N_",g,"*K_",g,", N_",g,"*K_",g,"] cov_",g,"; \n"))
  out$model <- paste0(stan.prior(paste0("sd_",g,"_",r), add.type = g, prior = prior ,
                      ind = ifelse(length(r) == 1, "", list(1:length(r)))[[1]]))
  
  if (length(r) == 1) {
    out$data <- paste0(out$data, "  real Z_",g,"[N]; \n")
    out$par <- paste0("  vector[N_",g,"] pre_",g,"; \n",
                      "  real<lower=0> sd_",g,"; \n")
    out$model <- paste0(out$model,"  pre_",g," ~ normal(0,1); \n")
    out$transD <- paste0("  vector[N_",g,"] r_",g,"; \n")
    out$transC <- paste0("  r_",g, " <- sd_",g," * (", 
                         if (c.cov) paste0("cov_",g," * "), "pre_",g,"); \n")
  }  
  else if (length(r) > 1) {
    out$data <- paste0(out$data,  "  row_vector[K_",g,"] Z_",g,"[N]; \n  int NC_",g,"; \n")
    out$par <- paste0("  matrix[N_",g,",K_",g,"] pre_",g,"; \n",
                      "  vector<lower=0>[K_",g,"] sd_",g,"; \n",
                      if (cor) paste0("  cholesky_factor_corr[K_",g,"] L_",g,"; \n"))
    out$model <- paste0(out$model, ifelse(cor, stan.prior(paste0("L_",g), prior = prior, add.type = g), ""),
                     "  to_vector(pre_",g,") ~ normal(0,1); \n")
    out$transD <- paste0("  vector[K_",g,"] r_",g,"[N_",g,"]; \n")
    if (c.cov) {
      if (cor) out$transC <- paste0("  r_",g," <- to_array(kronecker_cholesky(cov_",g,", L_",g,", sd_",g,") * ",
                                    "to_vector(pre_",g,"), N_",g,", K_",g,"); \n")
      else out$transC <- paste0("  r_",g," <- to_array(to_vector(rep_matrix(sd_",g,", N_",g,")) .* ",
                                "(cov_",g," * to_vector(pre_",g,")), N_",g,", K_",g,"); \n")
    }
    else out$transC <- paste0("  for (i in 1:N_",g,") { \n",
                              "    r_",g, "[i] <- sd_",g," .* (", if (cor) paste0("L_",g," * "), 
                              "to_vector(pre_",g,"[i])); \n  } \n")
    if (cor) {
      out$genD <- paste0("  corr_matrix[K_",g,"] Cor_",g,"; \n",
                         "  vector<lower=-1,upper=1>[NC_",g,"] cor_",g,"; \n")
      out$genC <- paste0("  Cor_",g," <- multiply_lower_tri_self_transpose(L_",g,"); \n",
                         collapse(unlist(lapply(2:length(r),function(i) lapply(1:(i-1), function(j)
                          paste0("  cor_",g,"[",(i-1)*(i-2)/2+j,"] <- Cor_",g,"[",j,",",i,"]; \n")))))) 
    }  
  }
  out
}

# Likelihoods in stan language
stan.llh <- function(family, link, predict = FALSE, add = FALSE, weights = FALSE, 
                     cens = FALSE, logllh = FALSE) {
  is.cat <- family %in% c("cumulative", "cratio", "sratio", "acat", "categorical")
  is.count <- family %in% c("poisson","negbinomial", "geometric")
  is.skew <- family %in% c("gamma","exponential","weibull")
  is.dich <- family %in% c("binomial", "bernoulli")
  
  simplify <- !(cens || predict || logllh) && (is.dich && link == "logit" || is.count && link == "log" ||
                family %in% c("cumulative", "categorical") && link == "logit" && !add) 
  n <- ifelse(predict || logllh || cens || weights || is.cat, "[n]", "")
  ns <- ifelse(add && (predict || logllh || cens || weights), "[n]", "")
  ilink <- ifelse((cens || predict || logllh) && (is.dich && link == "logit" || is.count && link == "log"), 
                  stan.ilink(link), "")
  
  lin.args <- paste0("eta",n,",sigma",ns)
  if (simplify) llh.pre <- switch(family,
     poisson = c("poisson_log", paste0("eta",n)), 
     negbinomial = c("neg_binomial_2_log", paste0("eta",n,",shape")),
     geometric = c("neg_binomial_2_log", paste0("eta",n,",1")),
     cumulative = c("ordered_logistic", "eta[n],b_Intercept"),
     categorical = c("categorical_logit", "to_vector(append_col(zero, eta[n] + etap[n]))"), 
     binomial = c("binomial_logit", paste0("max_obs",ns,",eta",n)), 
     bernoulli = c("bernoulli_logit", paste0("eta",n)))
  else llh.pre <- switch(ifelse(is.cat, "categorical", 
                         ifelse(family == "gaussian" && link == "log", "lognormal", family)),
     gaussian = c("normal", lin.args),
     student = c("student_t", paste0("nu,",lin.args)),
     cauchy = c("cauchy", lin.args),
     multinormal = c("multi_normal_cholesky", paste0("etam",n,",diag_pre_multiply(sigma,Lrescor)")),
     lognormal = c("lognormal", lin.args),
     poisson = c("poisson", paste0(ilink,"(eta",n,")")),
     negbinomial = c("neg_binomial_2", paste0(ilink,"(eta",n,"),shape")),
     geometric = c("neg_binomial_2", paste0(ilink,"(eta",n,"),1")),
     binomial = c("binomial", paste0("max_obs",ns,",",ilink,"(eta",n,")")),
     bernoulli = c("bernoulli", paste0(ilink,"(eta",n,")")), 
     gamma = c("gamma", paste0("shape,eta",n)), 
     exponential = c("exponential", paste0("eta",n)),
     weibull = c("weibull", paste0("shape,eta",n)), 
     categorical = c("categorical","p[n]"))
  
  type <- c("predict", "logllh", "cens", "weights")[match(TRUE, c(predict, logllh, cens, weights))]
  if (is.na(type)) type <- "general"
  addW <- ifelse(weights, "weights[n] * ", "")
  llh <- switch(type, 
    predict = paste0("    Y_pred[n] <- ", llh.pre[1],"_rng(",llh.pre[2],"); \n"),
    logllh = paste0("    log_llh[n] <- ", llh.pre[1], "_log(Y[n],",llh.pre[2],"); \n"),
    cens = paste0("if (cens[n] == 0) ", 
      ifelse(!weights, paste0("Y[n] ~ ", llh.pre[1],"(",llh.pre[2],"); \n"),
             paste0("increment_log_prob(", addW, llh.pre[1], "_log(Y[n],",llh.pre[2],")); \n")),
      "    else { \n",         
      "      if (cens[n] == 1) increment_log_prob(", addW, llh.pre[1], "_ccdf_log(Y[n],", llh.pre[2],")); \n",
      "      else increment_log_prob(", addW, llh.pre[1], "_cdf_log(Y[n],", llh.pre[2],")); \n",
      "    } \n"),
    weights = paste0("  lp_pre[n] <- ", llh.pre[1], "_log(Y[n],",llh.pre[2],"); \n"),
    general = paste0("  Y", n, " ~ ", llh.pre[1],"(",llh.pre[2],"); \n")) 
  llh
}

# prediction part in Stan
stan.eta <- function(family, link, f, p = NULL, group = list(), 
                     autocor = cor.arma(), max_obs = "max_obs") {
  is.linear <- family %in% c("gaussian", "student", "cauchy")
  is.ordinal <- family %in% c("cumulative", "cratio", "sratio", "acat") 
  is.skew <- family %in% c("gamma", "weibull", "exponential")
  is.count <- family %in% c("poisson", "negbinomial", "geometric")
  is.dich <- family %in% c("binomial", "bernoulli")
  is.multi <- family == "multinormal"
  
  eta <- list()
  # initialize eta
  eta$transD <- paste0("  vector[N] eta; \n", 
                       ifelse(length(p), paste0("  matrix[N,",max_obs[1],"-1] etap; \n"), ""),
                       ifelse(is.multi, "  vector[K_trait] etam[N_trait]; \n", ""))
  eta.multi <- ifelse(is.multi, "etam[m,k]", "eta[n]")
  
  # transform eta before it is passed to the likelihood
  ilink <- stan.ilink(link)
  eta$transform <- !(link == "identity" || family == "gaussian" && link == "log" ||
                     is.ordinal || family == "categorical" || is.count && link == "log" ||
                     is.dich && link == "logit")
  eta.ilink <- rep("", 2)
  if (eta$transform) {
    eta.ilink <- switch(family, c(paste0(ilink,"("), ")"),
                   gamma = c(paste0("shape/(",ilink,"("), "))"), 
                   exponential = c(paste0(ilink,"(-("), "))"), 
                   weibull = c(paste0("inv(",ilink,"(-("), ")/shape))"))
    if (autocor$q > 0) {
      eta$transC3 <- paste0("    ",eta.multi," <- ",eta.ilink[1], eta.multi, eta.ilink[2],"; \n")
      eta.ilink <- rep("", 2)  
    }
  }
  
  #define fixed, random and autocorrelation effects
  eta$transC1 <- paste0("  eta <- ", ifelse(length(f), "X*b", "rep_vector(0,N)"), 
                        if (autocor$p && is(autocor, "cor.arma")) " + Yar*ar", "; \n", if (length(p)) "  etap <- Xp * bp; \n")
  eta.re <- ifelse(length(group), collapse(" + Z_",group,"[n]*r_",group,"[",group,"[n]]"), "")
  eta.ma <- ifelse(autocor$q && is(autocor, "cor.arma"), " + Ema[n]*ma", "")
  if (nchar(eta.re) || nchar(eta.ma) || is.multi || nchar(eta.ilink[1])) {
    eta$transC2 <- paste0("    ",eta.multi," <- ",eta.ilink[1],"eta[n]", eta.ma, eta.re, eta.ilink[2],"; \n")
  }
  eta
}

# moving average autocorrelation in Stan
stan.ma <- function(family, link, autocor, group, levels, N) {
  is.linear <- family %in% c("gaussian", "student", "cauchy")
  is.multi <- family == "multinormal"
  ma <- list()
  if (is(autocor, "cor.arma") && autocor$q) {
    link.fun <- c(identity = "", log = "log", inverse = "inv")[link]
    if (!(is.linear || is.multi) && suppressWarnings(max(levels)) < N) 
      stop(paste0("moving-average models for family ",family," require a random effect with the same number \n",
                  "of levels as observations in the data"))
    e.ranef <- ifelse(!(is.linear || is.multi), paste0("r_", group[levels == max(levels)][[1]]), "")
    index <- ifelse(is.multi, "m,k", "n")
    ma$transD <- paste0("  row_vector[Kma] Ema[N]; \n", if(is.linear || is.multi) "  vector[N] e; \n") 
    ma$transC1 <- "  Ema <- Ema_pre; \n" 
    ma$transC2 <- paste0(ifelse(is.linear || is.multi, paste0("    e[n] <- ",link.fun,"(Y[",index,"]) - eta[n]", "; \n"), ""), 
                         "    for (i in 1:Kma) if (n+1-i > 0 && n < N && tgroup[n+1] == tgroup[n+1-i]) \n",
                         "      Ema[n+1,i] <- ",ifelse(is.linear || is.multi, "e[n+1-i]", paste0(e.ranef,"[n+1-i]")), "; \n")
  }
  ma
}

#stan code for posterior predictives and log likelihood
stan.genquant <- function(family, link, predict = FALSE, logllh = FALSE, add = FALSE) {
  if (predict || logllh) {
    trait <- ifelse(family == "multinormal", "_trait", "")
    out <- list(genC = paste0("  for (n in 1:N",trait,") { \n"))
    if (predict) {
      llh.pred <- stan.llh(family = family, link = link, add = add, predict = TRUE) 
      if (family %in% c("gaussian", "student", "cauchy", "gamma", "weibull", "exponential"))
        out$genD <- "  real Y_pred[N]; \n"  
      else if (family %in% c("binomial", "bernoulli", "poisson", "negbinomial", "geometric",
                             "categorical", "cumulative", "cratio", "sratio", "acat"))
        out$genD <- "  int Y_pred[N]; \n"
      else if (family == "multinormal")
        out$genD <-"  vector[K_trait] Y_pred[N_trait]; \n"
      out$genC <- paste0(out$genC, stan.llh(family = family, link = link, add = add, predict = TRUE))
    }
    if (logllh) {
      out$genD <- paste0(out$genD, "  vector[N",trait,"] log_llh; \n")
      out$genC <- paste0(out$genC, stan.llh(family = family, link = link, add = add, logllh = TRUE))
    }
    out$genC <- paste0(out$genC, "  } \n")
  }
  else out <- list()
  out
}

#stan code for user defined functions
stan.function <- function(kronecker = FALSE) {
  out <- NULL
  if (kronecker) out <- paste0(
    "  // calculate the cholesky factor of the kronecker covariance matrix \n",
    "  matrix kronecker_cholesky(matrix X, matrix L, vector sd) { \n",
    "    matrix[rows(X)*rows(L), cols(X)*cols(L)] kron; \n",
    "    matrix[rows(L), cols(L)] C; \n",
    "    int rX; \n",
    "    int rC; \n",
    "    C <- multiply_lower_tri_self_transpose(L); \n",
    "    rX <- rows(X); \n",
    "    rC <- rows(C); \n",
    "    for (i in 1:rX) { \n",
    "      for (j in 1:rC) { \n",
    "        for (k in 1:rX) { \n",
    "          for (l in 1:rC) { \n",
    "            kron[(k-1)*rC+l, (i-1)*rC+j] <- sd[l] * sd[j] * X[k,i] * C[l,j]; \n",
    "          } \n",
    "        } \n",
    "      } \n",
    "    } \n",
    "    return cholesky_decompose(kron); \n",
    "  } \n",
    "  // turn a vector into a 2 dimensional array \n",
    "  vector[] to_array(vector X, int N, int K) { \n",
    "    vector[K] Y[N]; \n",
    "    for (i in 1:N) \n",
    "      Y[i] <- segment(X, (i-1)*K+1, K); \n",
    "    return Y; \n",
    "  } \n")
  return(paste0("functions { \n", out, "} \n"))
}

# multinormal effects in Stan
stan.multi <- function(family, response) {
  out <- list()
  if (family == "multinormal") {
   out$genD <- paste0("  matrix[K_trait,K_trait] Rescor; \n",
    "  vector<lower=-1,upper=1>[NC_trait] rescor; \n")
   out$genC <- paste0("  Rescor <- multiply_lower_tri_self_transpose(Lrescor); \n",
        collapse(unlist(lapply(2:length(response),function(i) lapply(1:(i-1), function(j)
        paste0("  rescor[",(i-1)*(i-2)/2+j,"] <- Rescor[",j,",",i,"]; \n"))))))
  }
  out
}

# Ordinal effects in Stan
# 
# @return A vector of strings containing the ordinal effects in stan language
stan.ordinal <- function(family, link, partial = FALSE, max_obs = "max_obs", 
                     n = "", threshold = "flexible", simplify = TRUE) {
  is.ordinal <- family %in% c("cumulative", "cratio", "sratio", "acat")
  if (!(is.ordinal || family == "categorical")) return(list())
  ilink <- c(identity = "", log = "exp", inverse = "inv", sqrt = "square", logit = "inv_logit", 
             probit = "Phi", probit_approx = "Phi_approx", cloglog = "inv_cloglog")[link]
  th <- function(k) {
    sign <- ifelse(is.element(family, c("cumulative", "sratio"))," - ", " + ")
    ptl <- ifelse(partial, paste0(sign, "etap[n,k]"), "") 
    if (sign == " - ") out <- paste0("b_Intercept[",k,"]", ptl, " - eta[n]")
    else out <- paste0("eta[n]", ptl, " - b_Intercept[",k,"]")
  }  
  add.loop <- ifelse(n == "[n]", paste0("    for (k in (max_obs[n]+1):MAX_obs) p[n,k] <- 0.0; \n"), "")
  hd <- ifelse(rep(n == "[n]" && is.element(family, c("sratio","cratio")), 2), 
               c("head(", paste0(",max_obs",n,"-1)")), "")
  sc <- ifelse(family == "sratio", "1-", "")
  intercept <- paste0("  ", ifelse(family == "cumulative", "ordered", "vector"), "[",max_obs,"-1] b_Intercept; \n")
  
  out <- list()
  if (is.ordinal) {
    out$par <-  ifelse(threshold == "flexible", intercept,
      paste0("  real b_Intercept1; \n  real", ifelse(family == "cumulative", "<lower=0>", "")," delta; \n")) 
    out$transC1 <- ifelse(threshold == "equidistant", 
      paste0("  for (k in 1:(",max_obs,"-1)) { \n    b_Intercept[k] <- b_Intercept1 + (k-1.0)*delta; \n  } \n"), "")
    out$transD <- ifelse(threshold == "equidistant", intercept, "")
  }
  if (!(family %in% c("cumulative", "categorical") && ilink == "inv_logit" && n != "[n]" && simplify)) {
    out$transD <- paste0(out$transD, "  vector[",max_obs,"] p[N]; \n", 
      if (!family %in% c("cumulative", "categorical")) paste0("  vector[",max_obs,"-1] q[N]; \n"))
    if (family == "categorical" && ilink == "inv_logit") out$transC <- paste0(
      "    p[n,1] <- 1.0; \n",
      "    for (k in 2:max_obs",n,") { \n",
      "      p[n,k] <- exp(eta[n,k-1]); \n",
      "    } \n", add.loop,
      "    p[n] <- p[n]/sum(p[n]); \n")
    else if (family == "cumulative") out$transC2 <- paste0(
      "    p[n,1] <- ",ilink,"(",th(1),"); \n",
      "    for (k in 2:(max_obs",n,"-1)) { \n", 
      "      p[n,k] <- ",ilink,"(",th("k"),") - ",ilink,"(",th("k-1"),"); \n", 
      "    } \n",
      "    p[n,max_obs",n,"] <- 1 - ",ilink,"(",th(paste0("max_obs",n,"-1")),"); \n", 
      add.loop)
    else if (is.element(family,c("sratio", "cratio"))) out$transC2 <- paste0(
      "    for (k in 1:(max_obs",n,"-1)) { \n",
      "      q[n,k] <- ",sc, ilink,"(",th("k"),"); \n",
      "      p[n,k] <- 1-q[n,k]; \n",
      "      for (kk in 1:(k-1)) p[n,k] <- p[n,k] * q[n,kk]; \n", 
      "    } \n",
      "    p[n,max_obs",n,"] <- prod(",hd[1],"q[n]",hd[2],"); \n",
      add.loop)
    else if (family == "acat") 
      if (ilink == "inv_logit") out$transC2 <- paste0(
        "    p[n,1] <- 1.0; \n",
        "    for (k in 1:(max_obs",n,"-1)) { \n",
        "      q[n,k] <- ",th("k"),"; \n",
        "      p[n,k+1] <- q[n,1]; \n",
        "      for (kk in 2:k) p[n,k+1] <- p[n,k+1] + q[n,kk]; \n",
        "      p[n,k+1] <- exp(p[n,k+1]); \n",
        "    } \n", add.loop,
        "    p[n] <- p[n]/sum(p[n]); \n")
    else out$transC2 <- paste0(                   
      "    for (k in 1:(max_obs",n,"-1)) \n",
      "      q[n,k] <- ",ilink,"(",th("k"),"); \n",
      "    for (k in 1:max_obs",n,") { \n",     
      "      p[n,k] <- 1.0; \n",
      "      for (kk in 1:(k-1)) p[n,k] <- p[n,k] * q[n,kk]; \n",
      "      for (kk in k:(max_obs",n,"-1)) p[n,k] <- p[n,k] * (1-q[n,kk]); \n",      
      "    } \n", add.loop,
      "    p[n] <- p[n]/sum(p[n]); \n")
  }
  out
}

# Priors in Stan
# 
# Define priors for parameters in Stan language
# 
# @param par A vector of parameter names
# @param prior A named list of strings containing the priors for parameters in \code{par}.
# @param ind An optional index to allow for different priors for different elements of a parameter vector (see 'Examples')
# @param s An integer >= 0 defining the number of spaces in front of the output string
# @param ... Other potential arguments
# @inheritParams brm
# 
# @return A vector of character strings each defining a line of code in stan language that
#   defines the prior of a parameter in \code{par}. If a parameter has has no corresponding prior in \code{prior} 
#   and also no internal default in \code{stan.prior} (see 'Details'), an empty string is returned.
#      
# @examples 
# stan.prior(c("b_x1","sd_Site","sd_obs"), 
#           prior = list(b_ = "normal(0,1)", sd = "cauchy(0,2.5)"))
#  
# # returns a cauchy prior for sd_Site and a uniform prior for sd_obs                
# stan.prior(c("sd_Site","sd_obs"), 
#           prior = list(sd = "cauchy(0,5)", sd_obs = "uniform(0,100)"))            
# 
# # returns a uniform prior for the first element of b                      
# stan.prior("b", prior = list(beta = "uniform(0,10)"), ind = "[1]")
# 
# # returns the default prior for nu
# stan.prior("nu")
stan.prior = function(par, prior = list(), add.type = NULL, ind = rep("", length(par)), 
                      partial = FALSE, s = 2) { 
  if (length(par) != length(ind)) 
    stop("Arguments par and ind must have the same length")
  if (length(par) > 1 && all(ind == ""))
    warning("Indices are all zero")
  type <- unique(unlist(lapply(par, function(par) 
    unlist(regmatches(par, gregexpr("[^_]*", par)))[1])))
  if (length(type) > 1) stop("Only one parameter type can be handled at once")
  if (!is.null(add.type)) {
    type <- c(type, paste0(type,"_",add.type))
    if (!all(grepl(type[2], par))) 
      stop("Additional parameter type not present in all parameters")
  }  
  default.prior <- list(b = "", bp = "", sigma = "cauchy(0,5)", delta = "", ar = "", 
    ma = "", L = "lkj_corr_cholesky(1.0)", sd = "cauchy(0,5)", Lrescor = "lkj_corr_cholesky(1.0)",  
    nu = "uniform(1,100)", shape = "gamma(0.01,0.01)") 
  if (type[2] %in% names(prior)) base.prior <- prior[[type[2]]]
  else if (type[1] %in% names(prior)) base.prior <- prior[[type[1]]]
  else base.prior <- default.prior[[type[1]]]
  
  if (type[1] == "b" && partial) type[1] <- "bp"
  type <- ifelse(is.na(type[2]), type[1], type[2])  
  if (any(par %in% names(prior)))
    out <- sapply(1:length(par), function(i, par, ind) {
      if (ind[i] != "") 
        ind[i] <- paste0("[",ind[i],"]")
      if (!par[i] %in% names(prior)) 
        prior[[par[i]]] <- base.prior
      if (nchar(paste0(prior[[par[i]]], "")) > 0)  
        return(paste0(collapse(rep(" ", s)), type, ind[i], " ~ ", prior[[par[i]]], "; \n"))
      else return("") }, 
    par = par, ind = ind)
  else if (nchar(paste0(base.prior, "")) > 0)
    out <- paste0(collapse(rep(" ", s)), 
      ifelse(identical(type,"bp"), "to_vector(bp)", type), " ~ ", base.prior, "; \n")
  else out <- ""
  return(collapse(out))
}

# stan code to sample from priors seperately
stan.rngprior <- function(sample.prior, priors = "", family = "gaussian", fixed = NULL,
                          partial = NULL, response = NULL, random = list(), group = list()) {
  out <- list()
  if (sample.prior) {
    priors <- gsub(" ", "", paste0("\n",priors))
    pars <- gsub("\\\n|to_vector\\(|\\)", "", regmatches(priors, gregexpr("\\\n[^~]+", priors))[[1]])
    take <- !grepl("^pre_", pars)
    pars <- rename(pars[take], symbols = c("^L_", "^rescorL_"), subs = c("cor_", "rescor_"), 
                   fixed = FALSE)
    dis <- gsub("~", "", regmatches(priors, gregexpr("~[^\\(]+", priors))[[1]])[take]
    args <- regmatches(priors, gregexpr("\\([^;]+\\);", priors))[[1]][take]
    
    #rename parameters containing indices
    has_ind <- grepl("\\[[[:digit:]]+\\]", pars)
    pars[has_ind] <- sapply(pars[has_ind], function(par) {
      ind <- regmatches(par, gregexpr("\\[[[:digit:]]+\\]", par))
      ind <- as.numeric(substr(ind, 2, nchar(ind) - 1))
      if (grepl("^b\\[", par)) par <- paste0("b_",fixed[ind])
      else if (grepl("^bp\\[", par)) par <- paste0("b_",partial[ind])
      else if (grepl("^sigma\\[", par)) par <- paste0("sigma_",response[ind])
      else if (grepl("^sd_", par)) {
        ind_group <- which(grepl(paste0("^sd_",group), par))
        par <- gsub("\\[[[:digit:]]+\\]", paste0("_",random[[ind_group]][ind]), par)
      }
      return(par)})
    
    #special treatment of lkj_corr_cholesky priors
    args <- ifelse(grepl("corr_cholesky$", dis), paste0("(2,", substr(args, 2, nchar(args)-1), "[1,2];"), args)
    dis <- sub("corr_cholesky$", "corr", dis)
    
    #distinguish between bounded and unbounded parameters
    bound <- grepl("^sd|^sigma|^shape$|^nu$", pars) | family == "cumulative" & grepl("^delta$", pars)
    if (any(bound)) {
      out$par <- collapse("  real<lower=0> prior_",pars[bound],"; \n")
      out$model <- collapse("  prior_",pars[bound]," ~ ",dis[bound],args[bound]," \n")
    }
    if (any(!bound)) {
      out$genD <- collapse("  real prior_",pars[!bound],"; \n")
      out$genC <- collapse("  prior_",pars[!bound]," <- ",dis[!bound],"_rng",args[!bound]," \n")
    }
  }
  out
}

# find the inverse link to a given link function
stan.ilink <- function(link) {
  switch(link, identity = "", log = "exp", inverse = "inv", sqrt = "square", 
         logit = "inv_logit", probit = "Phi", probit_approx = "Phi_approx", cloglog = "inv_cloglog")
}
