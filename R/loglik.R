#compute WAIC and LOO using the 'loo' package
calculate_ic <- function(x, ic = c("waic", "loo")) {
  ic <- match.arg(ic)
  ee <- extract_effects(x$formula)
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples") 
  IC <- do.call(eval(parse(text = paste0("loo::",ic))), list(loglik(x)))
  class(IC) <- c("ic", "loo")
  return(IC)
}

#compare information criteria of different models
compare_ic <- function(x, ic = c("waic", "loo")) {
  ic <- match.arg(ic)
  n_models <- length(x)
  compare_matrix <- matrix(0, nrow = n_models * (n_models-1) / 2, ncol = 2)
  rnames <- rep("", nrow(compare_matrix))
  n <- 1
  for (i in 1:(n_models-1)) {
    for (j in (i+1):n_models) {
      temp <- loo::compare(x[[j]], x[[i]])
      compare_matrix[n,] <- c(-2 * temp$elpd_diff, 2 * temp$se) 
      rnames[n] <- paste(names(x)[i], "-", names(x)[j])
      n <- n + 1
    }
  }
  rownames(compare_matrix) <- rnames
  colnames(compare_matrix) <- c(toupper(ic), "SE")
  compare_matrix
}

loglik_gaussian <- function(n, data, samples, link) {
  sigma <- if (!is.null(samples$sigma)) samples$sigma
           else data$sigma
  out <- if (is.null(data$cens) || data$cens[n] == 0)
    dnorm(data$Y[n], mean = ilink(samples$eta[,n], link), 
          sd = sigma, log = TRUE)
  else if (data$cens[n] == 1)
    pnorm(data$Y[n], mean = ilink(samples$eta[,n], link), 
          sd = sigma, lower.tail = FALSE, log.p = TRUE)
  else if (data$cens[n] == -1)
    pnorm(data$Y[n], mean = ilink(samples$eta[,n], link), 
          sd = sigma, log.p = TRUE)
  if ("weights" %in% names(data)) out <- out * data$weights[n]
  out
}

loglik_student <- function(n, data, samples, link) {
  sigma <- if (!is.null(samples$sigma)) samples$sigma
           else data$sigma
  out <- if (is.null(data$cens) || data$cens[n] == 0)
    dstudent(data$Y[n], df = samples$nu, mu = ilink(samples$eta[,n], link), 
             sigma = sigma, log = TRUE)
  else if (data$cens[n] == 1)
    pstudent(data$Y[n], df = samples$nu, mu = ilink(samples$eta[,n], link), 
             sigma = sigma, lower.tail = FALSE, log.p = TRUE)
  else if (data$cens[n] == -1)
    pstudent(data$Y[n], df = samples$nu, mu = ilink(samples$eta[,n], link), 
             sigma = sigma, log.p = TRUE)
  if ("weights" %in% names(data)) out <- out * data$weights[n]
  out
}

loglik_cauchy <- function(n, data, samples, link) {
  sigma <- if (!is.null(samples$sigma)) samples$sigma
           else data$sigma
  out <- if (is.null(data$cens) || data$cens[n] == 0)
    dstudent(data$Y[n], df = 1, mu = ilink(samples$eta[,n], link), 
             sigma = sigma, log = TRUE)
  else if (data$cens[n] == 1)
    pstudent(data$Y[n], df = 1, mu = ilink(samples$eta[,n], link), 
             sigma = sigma, lower.tail = FALSE, log.p = TRUE)
  else if (data$cens[n] == -1)
    pstudent(data$Y[n], df = 1, mu = ilink(samples$eta[,n], link), 
             sigma = sigma, log.p = TRUE)
  if ("weights" %in% names(data)) out <- out * data$weights[n]
  out
}

loglik_lognormal <- function(n, data, samples, link) {
  sigma <- if (!is.null(samples$sigma)) samples$sigma
           else data$sigma
  out <- if (is.null(data$cens) || data$cens[n] == 0)
    dlnorm(data$Y[n], meanlog = samples$eta[,n], 
           sdlog = sigma, log = TRUE)
  else if (data$cens[n] == 1)
    plnorm(data$Y[n], meanlog = ilink(samples$eta[,n], link), 
           sdlog = sigma, lower.tail = FALSE, log.p = TRUE)
  else if (data$cens[n] == -1)
    plnorm(data$Y[n], meanlog = ilink(samples$eta[,n], link), 
           sdlog = sigma, log.p = TRUE)
  if ("weights" %in% names(data)) out <- out * data$weights[n]
  out
}

loglik_multinormal <- function(n, data, samples, link) {
  nobs <- data$N_trait * data$K_trait
  out <- sapply(1:nrow(samples$eta), function(i) 
    dmultinormal(data$Y[n,], Sigma = samples$Sigma[i,,], log = TRUE,
                 mu = samples$eta[i, seq(n, nobs, data$N_trait)]))
  if ("weights" %in% names(data)) out <- out * data$weights[n]
  out
}

loglik_binomial <- function(n, data, samples, link) {
  max_obs <- ifelse(length(data$max_obs) > 1, data$max_obs[n], data$max_obs) 
  out <- if (is.null(data$cens) || data$cens[n] == 0)
    dbinom(data$Y[n], size = max_obs, 
           prob = ilink(samples$eta[,n], link), log = TRUE)
  else if (data$cens[n] == 1)
    pbinom(data$Y[n], size = max_obs, 
           prob = ilink(samples$eta[,n], link), 
           lower.tail = FALSE, log.p = TRUE)
  else if (data$cens[n] == -1)
    pbinom(data$Y[n], size = max_obs, 
           prob = ilink(samples$eta[,n], link), log.p = TRUE)
  if ("weights" %in% names(data)) out <- out * data$weights[n]
  out
}  

loglik_bernoulli <- function(n, data, samples, link) {
  out <- if (is.null(data$cens) || data$cens[n] == 0)
    dbinom(data$Y[n], size = 1, 
           prob = ilink(samples$eta[,n], link), log = TRUE)
  else if (data$cens[n] == 1)
    pbinom(data$Y[n], size = 1, 
           prob = ilink(samples$eta[,n], link), 
           lower.tail = FALSE, log.p = TRUE)
  else if (data$cens[n] == -1)
    pbinom(data$Y[n], size = 1, 
           prob = ilink(samples$eta[,n], link), 
           log.p = TRUE)
  if ("weights" %in% names(data)) out <- out * data$weights[n]
  out
}

loglik_poisson <- function(n, data, samples, link) {
  out <- if (is.null(data$cens) || data$cens[n] == 0)
    dpois(data$Y[n], lambda = ilink(samples$eta[,n], link), log = TRUE)
  else if (data$cens[n] == 1)
    ppois(data$Y[n], lambda = ilink(samples$eta[,n], link), 
          lower.tail = FALSE, log.p = TRUE)
  else if (data$cens[n] == -1)
    ppois(data$Y[n], lambda = ilink(samples$eta[,n], link), 
          log.p = TRUE)
  if ("weights" %in% names(data)) out <- out * data$weights[n]
  out
}

loglik_negbinomial <- function(n, data, samples, link) {
  out <- if (is.null(data$cens) || data$cens[n] == 0)
    dnbinom(data$Y[n], mu = ilink(samples$eta[,n], link), 
            size = samples$shape, log = TRUE)
  else if (data$cens[n] == 1)
    pnbinom(data$Y[n], mu = ilink(samples$eta[,n], link), 
            size = samples$shape, lower.tail = FALSE, log.p = TRUE)
  else if (data$cens[n] == -1)
    pnbinom(data$Y[n], mu = ilink(samples$eta[,n], link), 
            size = samples$shape, log.p = TRUE)
  if ("weights" %in% names(data)) out <- out * data$weights[n]
  out
}

loglik_geometric <- function(n, data, samples, link) {
  out <- if (is.null(data$cens) || data$cens[n] == 0)
    dnbinom(data$Y[n], mu = ilink(samples$eta[,n], link), 
            size = 1, log = TRUE)
  else if (data$cens[n] == 1)
    pnbinom(data$Y[n], mu = ilink(samples$eta[,n], link), 
            size = 1, lower.tail = FALSE, log.p = TRUE)
  else if (data$cens[n] == -1)
    pnbinom(data$Y[n], mu = ilink(samples$eta[,n], link), 
            size = 1, log.p = TRUE)
  if ("weights" %in% names(data)) out <- out * data$weights[n]
  out
}

loglik_exponential <-  function(n, data, samples, link) {
  out <- if (is.null(data$cens) || data$cens[n] == 0)
    dexp(data$Y[n], rate = ilink(-samples$eta[,n], link), log = TRUE)
  else if (data$cens[n] == 1)
    pexp(data$Y[n], rate = ilink(-samples$eta[,n], link), 
         lower.tail = FALSE, log.p = TRUE)
  else if (data$cens[n] == -1)
    pexp(data$Y[n], rate = ilink(-samples$eta[,n], link), 
         log.p = TRUE)
  if ("weights" %in% names(data)) out <- out * data$weights[n]
  out
}

loglik_gamma <- function(n, data, samples, link) {
  out <- if (is.null(data$cens) || data$cens[n] == 0)
    dgamma(data$Y[n], shape = samples$shape, log = TRUE,
           scale = ilink(samples$eta[,n], link) / samples$shape)
  else if (data$cens[n] == 1)
    pgamma(data$Y[n], shape = samples$shape, 
           scale = ilink(samples$eta[,n], link) / samples$shape,
           lower.tail = FALSE, log.p = TRUE)
  else if (data$cens[n] == -1)
    pgamma(data$Y[n], shape = samples$shape, 
           scale = ilink(samples$eta[,n], link) / samples$shape,
           log.p = TRUE)
  if ("weights" %in% names(data)) out <- out * data$weights[n]
  out
}

loglik_weibull <- function(n, data, samples, link) {
  out <- if (is.null(data$cens) || data$cens[n] == 0)
    dweibull(data$Y[n], shape = samples$shape, log = TRUE,
             scale = 1/(ilink(-samples$eta[,n]/samples$shape, link)))
  else if (data$cens[n] == 1)
    pweibull(data$Y[n], shape = samples$shape, 
             scale = 1/(ilink(-samples$eta[,n]/samples$shape, link)),
             log.p = TRUE, lower.tail = FALSE)
  else if (data$cens[n] == -1)
    pweibull(data$Y[n], shape = samples$shape, 
             scale = 1/(ilink(-samples$eta[,n]/samples$shape, link)),
             log.p = TRUE)
  if ("weights" %in% names(data)) out <- out * data$weights[n]
  out
}

loglik_categorical <- function(n, data, samples, link) {
  max_obs <- ifelse(length(data$max_obs) > 1, data$max_obs[n], data$max_obs) 
  if (link == "logit") {
    p <- cbind(rep(0, nrow(samples$eta)), samples$eta[,n,1:(max_obs-1)])
    out <- p[,data$Y[n]] - log(rowSums(exp(p)))
  }
  else stop(paste("Link",link,"not supported"))
  if ("weights" %in% names(data)) out <- out * data$weights[n]
  out
}

loglik_cumulative <- function(n, data, samples, link) {
  max_obs <- ifelse(length(data$max_obs) > 1, data$max_obs[n], data$max_obs) 
  y <- data$Y[n]
  if (y == 1) out <- log(ilink(samples$eta[,n,1], link))
  else if (y == max_obs) out <- log(1 - ilink(samples$eta[,n,y-1], link))
  else out <- log(ilink(samples$eta[,n,y], link) - 
                  ilink(samples$eta[,n,y-1], link))
  if ("weights" %in% names(data)) out <- out * data$weights[n]
  out
}

loglik_sratio <- function(n, data, samples, link) {
  max_obs <- ifelse(length(data$max_obs) > 1, data$max_obs[n], data$max_obs) 
  y <- data$Y[n]
  q <- sapply(1:min(y, max_obs-1), function(k) 
    1 - ilink(samples$eta[,n,k], link))
  if (y == 1) out <- log(1 - q[,1])
  else if (y == 2) out <- log(1 - q[,2]) + log(q[,1])
  else if (y == max_obs) out <- rowSums(log(q))
  else out <- log(1 - q[,y]) + rowSums(log(q[,1:(y-1)]))
  if ("weights" %in% names(data)) out <- out * data$weights[n]
  out
}

loglik_cratio <- function(n, data, samples, link) {
  max_obs <- ifelse(length(data$max_obs) > 1, data$max_obs[n], data$max_obs) 
  y <- data$Y[n]
  q <- sapply(1:min(y, max_obs-1), function(k) 
    ilink(samples$eta[,n,k], link))
  if (y == 1) out <- log(1 - q[,1])
  else if (y == 2) out <- log(1 - q[,2]) + log(q[,1])
  else if (y == max_obs) out <- rowSums(log(q))
  else out <- log(1 - q[,y]) + rowSums(log(q[,1:(y-1)]))
  if ("weights" %in% names(data)) out <- out * data$weights[n]
  out
}

loglik_acat <- function(n, data, samples, link) {
  max_obs <- ifelse(length(data$max_obs) > 1, data$max_obs[n], data$max_obs) 
  y <- data$Y[n]
  if (link == "logit") {
    q <- sapply(1:(max_obs-1), function(k) samples$eta[,n,k])
    p <- cbind(rep(0, nrow(samples$eta)), q[,1], 
               matrix(0, nrow = nrow(samples$eta), ncol = max_obs - 2))
    if (max_obs > 2) 
      p[,3:max_obs] <- sapply(3:max_obs, function(k) rowSums(q[,1:(k-1)]))
    out <- p[,y] - log(rowSums(exp(p)))
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
    out <- log(p[,y]) - log(apply(p, 1, sum))
  }
  if ("weights" %in% names(data)) out <- out * data$weights[n]
  out
}  

