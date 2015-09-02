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
  max_obs <- ifelse(length(data$max_obs) > 1, data$max_obs[n], data$max_obs) 
  rbinom(nrow(samples$eta), size = max_obs, 
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