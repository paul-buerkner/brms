dstudent <- function(x, df = stop("df is required"), mu = 0, sigma = 1, log = FALSE) {
  # density of student's distribution 
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

dmulti_normal <- function(x, mu, Sigma, log = TRUE,
                          check = FALSE) {
  # density of the multivariate normal distribution 
  # not vectorized to increase speed when x is only a vector not a matrix
  # Args:
  #   x: the value(s) at which the density should be evaluated
  #   mu: mean vector
  #   sigma: covariance matrix
  #   log: return on log scale?
  #   check: check arguments for validity?
  # Returns:
  #   density of the multi_normal distribution a values x
  p <- length(x)
  if (check) {
    if (length(mu) != p) {
      stop("dimension of mu is incompatible")
    }
    if (!all(dim(Sigma) == c(p, p))) {
      stop("dimension of Sigma is incompatible")
    }
    if (!isSymmetric(Sigma, tol = sqrt(.Machine$double.eps), 
                     check.attributes = FALSE)) {
      stop("Sigma must be a symmetric matrix")
    }
  }
  rooti <- backsolve(chol(Sigma), diag(p))
  quads <- colSums((crossprod(rooti, (x - mu)))^2)
  out <- -(p / 2) * log(2 * pi) + sum(log(diag(rooti))) - .5 * quads
  if (!log) 
    out <- exp(out)
  out
}

rmulti_normal <- function(n, mu, Sigma, check = FALSE) {
  # random values of the multivariate normal distribution 
  # Args:
  #   n: number of random values
  #   mu: mean vector
  #   sigma: covariance matrix
  #   check: check arguments for validity?
  # Returns:
  #   n samples of multi_normal distribution of dimension length(mu) 
  p <- length(mu)
  if (check) {
    if (!(is.wholenumber(n) && n > 0)) {
      stop("n must be a positive integer")
    }
    if (!all(dim(Sigma) == c(p, p))) {
      stop("dimension of Sigma is incompatible")
    }
    if (!isSymmetric(Sigma, tol = sqrt(.Machine$double.eps), 
                     check.attributes = FALSE)) {
      stop("Sigma must be a symmetric matrix")
    }
  }
  samples <- matrix(rnorm(n * p), nrow = n, ncol = p)
  mu + samples %*% chol(Sigma)
}

dmulti_student <- function(x, df, mu, Sigma, log = TRUE,
                           check = FALSE) {
  # density of the multivariate student-t distribution 
  # Args:
  #   x: the value(s) at which the density should be evaluated
  #   df: degrees of freedom
  #   mu: mean vector
  #   sigma: covariance matrix
  #   log: return on log scale?
  #   check: check arguments for validity?
  # Returns:
  #   density of the multi_student distribution a values x
  if (is.vector(x)) {
    x <- matrix(x, ncol = length(x))
  }
  p <- ncol(x)
  if (check) {
    if (df <= 0) {
      stop("df must be greater zero")
    }
    if (length(mu) != p) {
      stop("dimension of mu is incompatible")
    }
    if (!all(dim(Sigma) == c(p, p))) {
      stop("dimension of Sigma is incompatible")
    }
    if (!isSymmetric(Sigma, tol = sqrt(.Machine$double.eps), 
                     check.attributes = FALSE)) {
      stop("Sigma must be a symmetric matrix")
    }
  }
  chol_Sigma <- chol(Sigma)
  rooti <- backsolve(chol_Sigma, t(x) - mu, transpose = TRUE)
  quads <- colSums(rooti^2)
  out <- lgamma((p + df)/2) - (lgamma(df / 2) + sum(log(diag(chol_Sigma))) + 
         p / 2 * log(pi * df)) - 0.5 * (df + p) * log1p(quads / df)
  if (!log) 
    out <- exp(out)
  out
}

rmulti_student <- function(n, df, mu, Sigma, log = TRUE, 
                           check = FALSE) {
  # random values of the multivariate student-t distribution 
  # Args:
  #   n: number of random values
  #   df: degrees of freedom
  #   mu: mean vector
  #   sigma: covariance matrix
  #   check: check arguments for validity?
  # Returns:
  #   n samples of multi_student distribution of dimension length(mu) 
  p <- length(mu)
  if (check) {
    if (df <= 0) {
      stop("df must be greater zero")
    }
  }
  samples <- rmulti_normal(n, mu = rep(0, p), Sigma = Sigma, check = check) / 
               sqrt(rchisq(n, df = df) / df)
  sweep(samples, 2, mu, "+")
}

dvon_mises <- function(x, mu, kappa, log = FALSE) {
  # density function of the von Mises distribution
  # CircStats::dvm has support within [0, 2*pi], 
  # but in brms we use [-pi, pi]
  # Args:
  #    mu: location parameter
  #    kappa: precision parameter
  out <- CircStats::dvm(x + base::pi, mu + base::pi, kappa)
  if (log) {
    out <- log(out)
  }
  out
}

pvon_mises <- function(q, mu, kappa, lower.tail = TRUE, 
                       log.p = FALSE, ...) {
  # distribution function of the von Mises distribution
  q <- q + base::pi
  mu <- mu + base::pi
  out <- .pvon_mises(q, mu, kappa, ...)
  if (!lower.tail) {
    out <- 1 - out
  }
  if (log.p) {
    out <- log(out)
  }
  out
}

# vectorized version of CircStats::pvm
.pvon_mises <- Vectorize(CircStats::pvm, c("theta", "mu", "kappa"))

rvon_mises <- function(n, mu, kappa) {
  # sample random numbers from the von Mises distribution
  stopifnot(n %in% c(1, max(length(mu), length(kappa))))
  mu <- mu + base::pi
  .rvon_mises(1, mu, kappa) - base::pi
}

# vectorized version of CircStats::rvm
.rvon_mises <- Vectorize(CircStats::rvm, c("mean", "k"))

dexgauss <- function (x, mean, sigma, beta, log = FALSE) {
  # PDF of the exponentially modified gaussian distribution
  # Args:
  #   mean: mean of the distribution; mean = mu + beta
  #   sigma: SD of the gaussian comoponent
  #   beta: scale / inverse rate of the exponential component
  if (any(sigma < 0)) {
    stop2("sigma must be greater than 0.")
  }
  if (any(beta < 0)) {
    stop2("beta must be greater than 0.")
  }
  mu <- mean - beta
  z <- x - mu - sigma^2 / beta
  out <- ifelse(beta > 0.05 * sigma, 
    -log(beta) - (z + sigma^2 / (2 * beta)) / beta + log(pnorm(z / sigma)), 
    dnorm(x, mean = mu, sd = sigma, log = TRUE))
  if (!log) {
    out <- exp(out)
  }
  out
}

pexgauss <- function(q, mean, sigma, beta, 
                     lower.tail = TRUE, log.p = FALSE) {
  # CDF of the exponentially modified gaussian distribution
  # Args:
  #   see dexgauss
  if (any(sigma < 0)) {
    stop2("sigma must be greater than 0.")
  }
  if (any(beta < 0)) {
    stop2("beta must be greater than 0.")
  } 
  mu <- mean - beta
  z <- q - mu - sigma^2 / beta
  out <- ifelse(beta > 0.05 * sigma, 
    pnorm((q - mu) / sigma) - pnorm(z / sigma) * 
      exp(((mu + sigma^2 / beta)^2 - mu^2 - 2 * q * sigma^2 / beta) / 
            (2 * sigma^2)), 
    pnorm(q, mean = mu, sd = sigma))
  if (!lower.tail) {
    out <- 1 - out
  } 
  if (log.p) {
    out <- log(out) 
  } 
  out
}

rexgauss <- function(n, mean, sigma, beta) {
  # create random numbers of the exgaussian distribution
  # Args:
  #   see dexgauss
  if (any(sigma < 0)) {
    stop2("sigma must be greater than 0.")
  }
  if (any(beta < 0)) {
    stop2("beta must be greater than 0.")
  } 
  mu <- mean - beta
  rnorm(n, mean = mu, sd = sigma) + rexp(n, rate = 1 / beta)
}

dcategorical <- function(x, eta, ncat, link = "logit") {
  # density of the categorical distribution
  # 
  # Args:
  #   x: positive integers not greater than ncat
  #   eta: the linear predictor (of length or ncol ncat-1)  
  #   ncat: the number of categories
  #   link: the link function
  #
  # Returns:
  #   probabilities P(X = x)
  if (is.null(dim(eta))) 
    eta <- matrix(eta, nrow = 1)
  if (length(dim(eta)) != 2) 
    stop("eta must be a numeric vector or matrix")
  if (missing(ncat))
    ncat <- ncol(eta) + 1
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
  #   eta: the linear predictor (of length or ncol ncat-1)  
  #   ncat: the number of categories
  #   link: a character string naming the link
  #
  # Retruns: 
  #   probabilities P(x <= q)
  p <- dcategorical(1:max(q), eta = eta, ncat = ncat, link = link)
  do.call(cbind, lapply(q, function(j) rowSums(as.matrix(p[, 1:j]))))
}

dcumulative <- function(x, eta, ncat, link = "logit") {
  # density of the cumulative distribution
  #
  # Args: same as dcategorical
  if (is.null(dim(eta))) 
    eta <- matrix(eta, nrow = 1)
  if (length(dim(eta)) != 2) 
    stop("eta must be a numeric vector or matrix")
  if (missing(ncat)) 
    ncat <- ncol(eta) + 1
  mu <- ilink(eta, link)
  rows <- list(mu[, 1])
  if (ncat > 2) {
    .fun <- function(k) {
      mu[, k] - mu[, k - 1]
    }
    rows <- c(rows, lapply(2:(ncat - 1), .fun))
  }
  rows <- c(rows, list(1 - mu[, ncat - 1]))
  p <- do.call(cbind, rows)
  p[, x]
}

dsratio <- function(x, eta, ncat, link = "logit") {
  # density of the sratio distribution
  #
  # Args: same as dcategorical
  if (is.null(dim(eta))) 
    eta <- matrix(eta, nrow = 1)
  if (length(dim(eta)) != 2) 
    stop("eta must be a numeric vector or matrix")
  if (missing(ncat)) 
    ncat <- ncol(eta) + 1
  mu <- ilink(eta, link)
  rows <- list(mu[, 1])
  if (ncat > 2) {
    .fun <- function(k) {
      (mu[, k]) * apply(as.matrix(1 - mu[, 1:(k - 1)]), 1, prod)
    }
    rows <- c(rows, lapply(2:(ncat - 1), .fun))
  }
  rows <- c(rows, list(apply(1 - mu, 1, prod)))
  p <- do.call(cbind, rows)
  p[, x]
}

dcratio <- function(x, eta, ncat, link = "logit") {
  # density of the cratio distribution
  #
  # Args: same as dcategorical
  if (is.null(dim(eta))) 
    eta <- matrix(eta, nrow = 1)
  if (length(dim(eta)) != 2) 
    stop("eta must be a numeric vector or matrix")
  if (missing(ncat)) 
    ncat <- ncol(eta) + 1
  mu <- ilink(eta, link)
  rows <- list(1 - mu[, 1])
  if (ncat > 2) {
    .fun <- function(k) {
      (1 - mu[, k]) * apply(as.matrix(mu[, 1:(k - 1)]), 1, prod)
    }
    rows <- c(rows, lapply(2:(ncat - 1), .fun))
  }
  rows <- c(rows, list(apply(mu, 1, prod)))
  p <- do.call(cbind, rows)
  p[, x]
}

dacat <- function(x, eta, ncat, link = "logit") {
  # density of the acat distribution
  #
  # Args: same as dcategorical
  if (is.null(dim(eta))) 
    eta <- matrix(eta, nrow = 1)
  if (length(dim(eta)) != 2) 
    stop("eta must be a numeric vector or matrix")
  if (missing(ncat)) 
    ncat <- ncol(eta) + 1
  if (link == "logit") { # faster evaluation in this case
    p <- cbind(rep(1, nrow(eta)), exp(eta[,1]), 
               matrix(NA, nrow = nrow(eta), ncol = ncat - 2))
    if (ncat > 2) {
      .fun <- function(k) {
        rowSums(eta[, 1:(k-1)])
      }
      p[, 3:ncat] <- exp(sapply(3:ncat, .fun))
    }
  } else {
    mu <- ilink(eta, link)
    p <- cbind(apply(1 - mu[,1:(ncat - 1)], 1, prod), 
               matrix(0, nrow = nrow(eta), ncol = ncat - 1))
    if (ncat > 2) {
      .fun <- function(k) {
        apply(as.matrix(mu[, 1:(k - 1)]), 1, prod) * 
          apply(as.matrix(1 - mu[, k:(ncat - 1)]), 1, prod)
      }
      p[, 2:(ncat - 1)] <- sapply(2:(ncat - 1), .fun)
    }
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
  #   eta: the linear predictor (of length or ncol ncat-1)  
  #   ncat: the number of categories
  #   family: a character string naming the family
  #   link: a character string naming the link
  #
  # Returns: 
  #   probabilites P(x <= q)
  args <- list(1:max(q), eta = eta, ncat = ncat, link = link)
  p <- do.call(paste0("d", family), args)
  .fun <- function(j) {
    rowSums(as.matrix(p[, 1:j]))
  }
  do.call(cbind, lapply(q, .fun))
}