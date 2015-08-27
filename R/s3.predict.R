#' @export 
predict.brmsfit <- function(object, ...) {
  if (!is(object$fit, "stanfit") || !length(object$fit@sim)) 
    stop("The model does not contain posterior samples")
  ee <- extract.effects(object$formula, family = object$family)
  if (object$link == "log" && object$family == "gaussian" && length(ee$response) == 1) 
    object$family <- "lognormal"
  if (object$family == "gaussian" && length(ee$response) > 1)
    object$family <- "multinormal"
  
  #compute all necessary samples
  samples <- list(eta = linear.predictor(object))
  if (object$family %in% c("gaussian", "student", "cauchy", "lognormal", "multinormal") && !is.formula(ee$se))
    samples$sigma <- as.matrix(posterior.samples(object, parameters = "^sigma_"))
  if (object$family == "student") 
    samples$nu <- as.matrix(posterior.samples(object, parameters = "^nu$"))
  if (object$family %in% c("gamma", "weibull","negbinomial")) 
    samples$shape <- as.matrix(posterior.samples(object, parameters = "^shape$"))
  if (object$family == "multinormal") {
    samples$rescor <- as.matrix(posterior.samples(object, parameters = "^rescor_"))
    samples$Sigma <- cov_matrix(sd = samples$sigma, cor = samples$rescor)$cov
    message(paste("Computing posterior predictive samples of multinormal distribution. \n", 
            "This may take a while."))
  }
  
  #call predict functions
  predict_fun <- get(paste0("predict_",object$family))
  samples <- do.call(cbind, lapply(1:nrow(as.matrix(object$data$Y)), function(n) 
    do.call(predict_fun, list(n = n, data = object$data, samples = samples, link = object$link))))
  out <- do.call(cbind, lapply(c("mean", "sd", "quantile"), get.estimate, 
                               samples = samples, probs = c(0.025, 0.975)))
  colnames(out) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI")
  
  #sort rows in case of multinormal models
  if (object$family == "multinormal") {
    nobs <- object$data$N_trait * object$data$K_trait
    out <- out[unlist(lapply(1:object$data$K_trait, function(k) seq(k, nobs, object$data$K_trait))),]
  }
  out
}

predict_gaussian <- function(n, data, samples, link) {
  sigma <- if (!is.null(samples$sigma)) samples$sigma
           else data$sigma
  rnorm(nrow(samples$eta), mean = ilink(samples$eta[,n], link), sd = sigma)
}

predict_student <- function(n, data, samples, link) {
  sigma <- if (!is.null(samples$sigma)) samples$sigma
           else data$sigma
  rstudent(nrow(samples$eta), df = samples$nu, mu = ilink(samples$eta[,n], link), 
           sigma = sigma)
}

predict_cauchy <- function(n, data, samples, link) {
  sigma <- if (!is.null(samples$sigma)) samples$sigma
           else data$sigma
  rstudent(nrow(samples$eta), df = 1, mu = ilink(samples$eta[,n], link), sigma = sigma)
}

predict_lognormal <- function(n, data, samples, link) {
  sigma <- if (!is.null(samples$sigma)) samples$sigma
           else data$sigma
  rlnorm(nrow(samples$eta), meanlog = samples$eta[,n], sdlog = sigma)
}

predict_multinormal <- function(n, data, samples, link) {
  nobs <- data$N_trait * data$K_trait
  do.call(rbind, lapply(1:nrow(samples$eta), function(i) 
    rmultinormal(1, Sigma = samples$Sigma[i,,],
                 mu = samples$eta[i, seq(n, nobs, data$N_trait)])))
}

predict_binomial <- function(n, data, samples, link) {
  rbinom(nrow(samples$eta), size = data$max_obs, 
         prob = ilink(samples$eta[,n], link))
}  

predict_bernoulli <- function(n, data, samples, link) {
  rbinom(nrow(samples$eta), size = 1, 
         prob = ilink(samples$eta[,n], link))
}

predict_poisson <- function(n, data, samples, link) {
  rpois(nrow(samples$eta), lambda = ilink(samples$eta[,n], link))
}

predict_negbinomial <- function(n, data, samples, link) {
  rnbinom(nrow(samples$eta), mu = ilink(samples$eta[,n], link), 
          size = samples$shape)
}

predict_geometric <- function(n, data, samples, link) {
  rnbinom(nrow(samples$eta), mu = ilink(samples$eta[,n], link), size = 1)
}

predict_exponential <-  function(n, data, samples, link) {
  rexp(nrow(samples$eta), rate = ilink(-samples$eta[,n], link))
}

predict_gamma <- function(n, data, samples, link) {
  rgamma(nrow(samples$eta), shape = samples$shape,
         scale = ilink(samples$eta[,n], link) / samples$shape)
}

predict_weibull <- function(n, data, samples, link) {
  rweibull(nrow(samples$eta), shape = samples$shape,
           scale = 1/(ilink(-samples$eta[,n]/samples$shape, link)))
}

predict_categorical <- function(n, data, samples, link) {
  max_obs <- ifelse(length(data$max_obs) > 1, data$max_obs[n], data$max_obs) 
  if (link == "logit") {
    p <- exp(cbind(rep(0, nrow(samples$eta)), samples$eta[,n,1:(max_obs-1)]))
    p <- do.call(cbind, lapply(1:max_obs, function(j) rowSums(as.matrix(p[,1:j]))))
  }
  else stop(paste("Link",link,"not supported"))
  first_greater(p, target = runif(nrow(samples$eta), min = 0, max = p[,max_obs]))
}

predict_cumulative <- function(n, data, samples, link) {
  max_obs <- ifelse(length(data$max_obs) > 1, data$max_obs[n], data$max_obs) 
  p <- cbind(ilink(samples$eta[,n,1], link), 
             if (max_obs > 2) sapply(2:(max_obs-1), function(k) 
               ilink(samples$eta[,n,k], link) - ilink(samples$eta[,n,k-1], link)),
             1 - ilink(samples$eta[,n,max_obs-1], link))
  p <- do.call(cbind, lapply(1:max_obs, function(j) rowSums(as.matrix(p[,1:j]))))
  first_greater(p, target = runif(nrow(samples$eta), min = 0, max = p[,max_obs]))
}

predict_sratio <- function(n, data, samples, link) {
  max_obs <- ifelse(length(data$max_obs) > 1, data$max_obs[n], data$max_obs) 
  y <- data$Y[n]
  q <- sapply(1:(max_obs-1), function(k) 
    1 - ilink(samples$eta[,n,k], link))
  p <- cbind(1 - q[,1], 
             if (max_obs > 2) sapply(2:(max_obs-1), function(k)
               (1 - q[,k]) * apply(as.matrix(q[,1:(k-1)]), 1, prod)),
             apply(q, 1, prod))
  p <- do.call(cbind, lapply(1:max_obs, function(j) rowSums(as.matrix(p[,1:j]))))
  first_greater(p, target = runif(nrow(samples$eta), min = 0, max = p[,max_obs]))
}

predict_cratio <- function(n, data, samples, link) {
  max_obs <- ifelse(length(data$max_obs) > 1, data$max_obs[n], data$max_obs) 
  y <- data$Y[n]
  q <- sapply(1:(max_obs-1), function(k) 
    ilink(samples$eta[,n,k], link))
  p <- cbind(1 - q[,1], 
             if (max_obs > 2) sapply(2:(max_obs-1), function(k)
               (1 - q[,k]) * apply(as.matrix(q[,1:(k-1)]), 1, prod)),
             apply(q, 1, prod))
  p <- do.call(cbind, lapply(1:max_obs, function(j) rowSums(as.matrix(p[,1:j]))))
  first_greater(p, target = runif(nrow(samples$eta), min = 0, max = p[,max_obs]))
}

predict_acat <- function(n, data, samples, link) {
  max_obs <- ifelse(length(data$max_obs) > 1, data$max_obs[n], data$max_obs) 
  y <- data$Y[n]
  if (link == "logit") {
    q <- sapply(1:(max_obs-1), function(k) samples$eta[,n,k])
    p <- cbind(rep(0, nrow(samples$eta)), q[,1], 
               matrix(0, nrow = nrow(samples$eta), ncol = max_obs - 2))
    if (max_obs > 2) 
      p[,3:max_obs] <- sapply(3:max_obs, function(k) rowSums(q[,1:(k-1)]))
  }
  else {
    q <- sapply(1:(max_obs-1), function(k) 
      ilink(samples$eta[,n,k], link))
    p <- cbind(apply(1 - q[,1:(max_obs-1)], 1, prod), 
               matrix(0, nrow = nrow(samples$eta), ncol = max_obs - 1))
    if (max_obs > 2)
      p[,2:(max_obs-1)] <- sapply(2:(max_obs-1), function(k) 
        apply(as.matrix(q[,1:(k-1)]), 1, prod) * apply(as.matrix(1 - q[,k:(max_obs-1)]), 1, prod))
    p[,max_obs] <- apply(q[,1:(max_obs-1)], 1, prod)
  }
  p <- do.call(cbind, lapply(1:max_obs, function(j) rowSums(as.matrix(p[,1:j]))))
  first_greater(p, target = runif(nrow(samples$eta), min = 0, max = p[,max_obs]))
}  