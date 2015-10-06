fitted_response <- function(x, eta, data) {
  # comnpute fitted values on the response scale
  # Args:
  #   x: A brmsfit object
  #   eta: linear predictor
  #   data: data initially passed to Stan
  family <- x$family
  nresp <- length(extract_effects(x$formula, family = family)$response)
  is_catordinal <- indicate_ordinal(family) || family == "categorical"
  if (family == "gaussian" && x$link == "log" && nresp == 1) {
    family <- "lognormal"
  }
  
  # compute (mean) fitted values
  if (family == "binomial") {
    max_obs <- matrix(rep(data$max_obs, nrow(eta)), nrow = nrow(eta), byrow = TRUE)
    # scale eta from [0,1] to [0,max_obs]
    mu <- ilink(eta, x$link) * max_obs 
  } else if (family == "lognormal") {
    sigma <- get_sigma(x, data = data, method = "fitted", n = nrow(eta))
    mu <- ilink(eta + sigma^2 / 2, x$link)  
  } else if (family == "weibull") {
    shape <- posterior_samples(x, "^shape$")$shape
    mu <- 1 / (ilink(-eta / shape, x$link)) * gamma(1 + 1 / shape)  
  } else if (is_catordinal) {
    mu <- fitted_catordinal(eta, max_obs = data$max_obs, family = x$family,
                            link = x$link)
  } else if (indicate_hurdle(family)) {
    shape <- posterior_samples(x, "^shape$")$shape 
    mu <- fitted_hurdle(eta, shape = shape, N_trait = data$N_trait,
                        family = x$family, link = x$link)
  } else if (indicate_zero_inflated(family)) {
    mu <- fitted_zero_inflated(eta, N_trait = data$N_trait, link = x$link)
  } else {
    # for any other distribution, ilink(eta) is already the mean fitted value
    mu <- ilink(eta, x$link)
  }
  
  # fitted values for truncated models
  if (!(is.null(data$lb) && is.null(data$ub))) {
    if (family != "gaussian") {
      stop(paste("fitted values on the respone scale not yet implemented",
                 "for non-gaussian truncated models"))
    } else {
      # fitted values for truncated normal models
      sigma <- get_sigma(x, data = data, method = "fitted", n = nrow(eta))
      lb <- ifelse(is.null(data$lb), -Inf, data$lb)
      ub <- ifelse(is.null(data$ub), Inf, data$ub)
      zlb <- (lb - mu) / sigma
      zub <- (ub - mu) / sigma
      scale <- (dnorm(zlb) - dnorm(zub)) / (pnorm(zub) - pnorm(zlb))  
      mu <- mu + scale * sigma  
    }
  } 
  mu
}

fitted_catordinal <- function(eta, max_obs, family, link) {
  # compute fitted values for categorical and ordinal families
  ncat <- max(max_obs)
  # get probabilities of each category
  get_density <- function(n) {
    do.call(paste0("d", family), 
            list(1:ncat, eta = eta[, n, ], ncat = ncat, link = link))
  }
  aperm(abind(lapply(1:ncol(eta), get_density), along = 3), perm = c(1, 3, 2))
}

fitted_hurdle <- function(eta, shape, N_trait, family, link) {
  n_base <- 1:N_trait
  n_hu <- n_base + N_trait
  pre_mu <- ilink(eta[, n_base], link)
  # adjust pre_mu as it is no longer the mean of the truncated distributions
  if (family == "hurdle_poisson") {
    adjusted_mu <- pre_mu / (1 - exp(-pre_mu))
  } else if (family == "hurdle_negbinomial") {
    adjusted_mu <- pre_mu / (1 - (shape / (pre_mu + shape))^shape)
  } else {
    adjusted_mu <- pre_mu
  }
  # incorporate hurdle part
  pre_mu * (1 - ilink(eta[, n_hu], "logit")) 
}

fitted_zero_inflated <- function(eta, N_trait, link) {
  n_base <- 1:N_trait
  n_zi <- n_base + N_trait
  # incorporate zero-inflation part
  ilink(eta[, n_base], link) * (1 - ilink(eta[, n_zi], "logit")) 
}