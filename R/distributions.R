dstudent <- function(x, df = stop("df is required"), mu = 0, sigma = 1, log = FALSE) {
  # density of student's distribution 
  #
  # Args:
  #  x: the value(s) at which the density should be evaluated
  #  df: degrees of freedom
  #  mu: the mean
  #  sigma: the scale parameter
  #  log: logical; return on log scale?
  if (log) {
    dt((x - mu) / sigma, df = df, log = TRUE) - log(sigma)
  } else {
    dt((x - mu) / sigma, df = df) / sigma
  }
}


pstudent <- function(q, df = stop("df is required"), mu = 0, sigma = 1, 
                     lower.tail = TRUE, log.p = FALSE) {
  # distribution function of student's distribution
  #
  # Args:
  #  q: the value(s) at which the distribution should be evaluated
  #  df: degrees of freedom
  #  mu: the mean
  #  sigma: the scale parameter
  #  lower.tail: same as for every pdist function
  #  log.p: logical; return on log scale?
  pt((q - mu)/sigma, df = df, lower.tail = lower.tail, log.p = log.p)
}

qstudent <-  function(p, df = stop("df is required"), mu = 0, sigma = 1) {
  # quantiles of student's distribution
  #
  # Args:
  #  p: the probabilities to find quantiles for
  #  df: degrees of freedom
  #  mu: the mean
  #  sigma: the scale parameter
  mu + sigma * qt(p, df = df)
}

rstudent <-  function(n, df = stop("df is required"), mu = 0, sigma = 1) {
  # random values of student's distribution 
  #
  # Args:
  #  n: number of random values
  #  df: degrees of freedom
  #  mu: the mean
  #  sigma: the scale parameter
  mu + sigma * rt(n, df = df)
}

dmultinormal <- function(x, mu, Sigma, log = TRUE) {
  # density of the multinormal distribution
  #
  # Args:
  #   x: the value(s) at which the density should be evaluated
  #   mu: mean vector
  #   sigma: covariance matrix
  #   log: return on log scale?
  k <- length(x)
  rooti <- backsolve(chol(Sigma), diag(k))
  quads <- colSums((crossprod(rooti, (x - mu)))^2)
  out <- -(k / 2) * log(2 * pi) + sum(log(diag(rooti))) - .5 * quads
  if (!log) 
    out <- exp(out)
  out
}

rmultinormal <- function(n, mu, Sigma, check = FALSE) {
  # random values of the multinormal distribution 
  #
  # Args:
  #   n: number of random values
  #   mu: mean vector
  #   sigma: covariance matrix
  #   check: check sigma for symmetry?
  p <- length(mu)
  if (check) {
    if (!all(dim(Sigma) == c(p, p))) 
      stop("incompatible arguments")
    if (!isSymmetric(unname(Sigma)))
      stop("Sigma is not symmetric")
  }
  samples <- matrix(rnorm(n * p), nrow = n, ncol = p)
  mu + samples %*% chol(Sigma)
}

dinv_gaussian <- function (x, mu, lambda, log = FALSE) {
  # density of the inverse gaussian distribution
  #
  # Args:
  #   x: the value(s) at which the density should be evaluated
  #   mu: mean parameter
  #   shape: shape parameter
  #   log: return on log scale?
  if (any(mu <= 0)) 
    stop("mu must be positive")
  if (any(lambda <= 0)) 
    stop("lambda must be positive")
  out <- 0.5 * log(lambda / (2 * pi * x^3)) -
    lambda * (x - mu)^2 / (2 * mu^2 * x)
  if (!log) out <- exp(out)
  out
}

pinv_gaussian <- function (q, mu, lambda, lower.tail = TRUE, log.p = FALSE) {
  # quantiles of the inverse gaussian distribution
  #
  # Args:
  #  p: the probabilities to find quantiles for
  #  mu: mean oarameter
  #  lambda: shape parameter
  if (any(mu <= 0)) 
    stop("mu must be positive")
  if (any(lambda <= 0)) 
    stop("lambda must be positive")
  out <- pnorm(sqrt(lambda / q) * (q / mu - 1)) + 
    exp(2 * lambda / mu) * pnorm(-sqrt(lambda / q) * (q / mu + 1))
  if (!lower.tail)
    out <- 1 - out
  if (log.p)
    out <- log(out)
  out
}

rinv_gaussian <- function(n, mu, lambda) {
  # random values of the inverse gaussian distribution 
  #
  # Args:
  #  n: number of random values
  #  mu: mean oarameter
  #  lambda: shape parameter
  u <- runif(n)
  z <- rnorm(n)^2
  phi <- lambda / mu
  y1 <- 1 - 0.5 * (sqrt(z^2 + 4 * phi * z) - z) / phi
  mu * ifelse((1 + y1) * u > 1, 1 / y1, y1)
}

dcategorical <- function(x, eta, ncat, link = "logit") {
  # density of the categorical distribution
  # 
  # Args:
  #   x: positive integers not greater than ncat
  #   mu: the linear predictor (of length or ncol ncat-1)  
  #   ncat: the number of categories
  #   link: the link function
  if (is.null(dim(eta))) 
    eta <- matrix(eta, nrow = 1)
  if (length(dim(eta)) != 2 || !is.numeric(eta)) 
    stop("eta must be a numeric vector or matrix")
  if (missing(ncat)) 
    ncat <- ncol(eta)
  if (link == "logit") {
    p <- exp(cbind(rep(0, nrow(eta)), eta[, 1:(ncat - 1)]))
  } else {
    stop(paste("Link", link, "not supported"))
  }
  p <- p / rowSums(p)
  p[, x]
}

pcategorical <- function(q, eta, ncat, link = "logit") {
  # distribution functions for the categorical family
  #
  # Args:
  #   q: positive integers not greater than ncat
  #   mu: the linear predictor (of length or ncol ncat-1)  
  #   ncat: the number of categories
  #
  # Retruns: 
  #   probabilites P(x <= q)
  p <- dcategorical(1:max(q), eta = eta, ncat = ncat, link = link)
  do.call(cbind, lapply(q, function(j) rowSums(as.matrix(p[, 1:j]))))
}

dcumulative <- function(x, eta, ncat, link = "logit") {
  # density of the cumulative distribution
  #
  # Args: same as dcategorical
  if (is.null(dim(eta))) eta <- matrix(eta, nrow = 1)
  if (length(dim(eta)) != 2 || !is.numeric(eta)) 
    stop("eta must be a numeric vector or matrix")
  if (missing(ncat)) ncat <- ncol(eta)
  mu <- ilink(eta, link)
  p <- cbind(mu[, 1], 
             if (ncat > 2) 
               sapply(2:(ncat - 1), function(k) mu[, k] - mu[, k - 1]), 
             1 - mu[, ncat - 1])
  p[, x]
}

dsratio <- function(x, eta, ncat, link = "logit") {
  # density of the sratio distribution
  #
  # Args: same as dcategorical
  if (is.null(dim(eta))) eta <- matrix(eta, nrow = 1)
  if (length(dim(eta)) != 2 || !is.numeric(eta)) 
    stop("eta must be a numeric vector or matrix")
  if (missing(ncat)) ncat <- ncol(eta)
  mu <- ilink(eta, link)
  p <- cbind(mu[, 1], 
        if (ncat > 2) sapply(2:(ncat - 1), function(k)
          (mu[, k]) * apply(as.matrix(1 - mu[, 1:(k - 1)]), 1, prod)),
        apply(1 - mu, 1, prod))
  p[, x]
}

dcratio <- function(x, eta, ncat, link = "logit") {
  # density of the cratio distribution
  #
  # Args: same as dcategorical
  if (is.null(dim(eta))) 
    eta <- matrix(eta, nrow = 1)
  if (length(dim(eta)) != 2 || !is.numeric(eta)) 
    stop("eta must be a numeric vector or matrix")
  if (missing(ncat)) 
    ncat <- ncol(eta)
  mu <- ilink(eta, link)
  p <- cbind(1 - mu[, 1], 
        if (ncat > 2) sapply(2:(ncat - 1), function(k)
          (1 - mu[, k]) * apply(as.matrix(mu[, 1:(k - 1)]), 1, prod)),
        apply(mu, 1, prod))
  p[, x]
}

dacat <- function(x, eta, ncat, link = "logit") {
  # density of the acat distribution
  #
  # Args: same as dcategorical
  if (is.null(dim(eta))) 
    eta <- matrix(eta, nrow = 1)
  if (length(dim(eta)) != 2 || !is.numeric(eta)) 
    stop("eta must be a numeric vector or matrix")
  if (missing(ncat)) 
    ncat <- ncol(eta)
  
  if (link == "logit") { # faster evaluation in this case
    p <- cbind(rep(1, nrow(eta)), exp(eta[,1]), 
               matrix(NA, nrow = nrow(eta), ncol = ncat - 2))
    if (ncat > 2) 
      p[, 3:ncat] <- exp(sapply(3:ncat, function(k) rowSums(eta[, 1:(k-1)])))
  } else {
    mu <- ilink(eta, link)
    p <- cbind(apply(1 - mu[,1:(ncat - 1)], 1, prod), 
               matrix(0, nrow = nrow(eta), ncol = ncat - 1))
    if (ncat > 2)
      p[, 2:(ncat - 1)] <- sapply(2:(ncat - 1), function(k) 
        apply(as.matrix(mu[, 1:(k - 1)]), 1, prod) * 
          apply(as.matrix(1 - mu[, k:(ncat - 1)]), 1, prod))
    p[, ncat] <- apply(mu[, 1:(ncat - 1)], 1, prod)
  }
  p <- p / rowSums(p)
  p[, x]
}

pordinal <- function(q, eta, ncat, family, link = "logit") {
  # distribution functions for ordinal families
  #
  # Args:
  #   q: positive integers not greater than ncat
  #   mu: the linear predictor (of length or ncol ncat-1)  
  #   ncat: the number of categories
  #
  # Returns: 
  #   probabilites P(x <= q)
  p <- do.call(paste0("d", family), list(1:max(q), eta = eta, 
                                         ncat = ncat, link = link))
  do.call(cbind, lapply(q, function(j) rowSums(as.matrix(p[, 1:j]))))
}