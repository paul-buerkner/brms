#' The Student-t Distribution
#' 
#' Density, distribution function, quantile function and random generation 
#' for the Student-t distribution with location \code{mu}, scale \code{sigma},
#' and degrees of freedom \code{df}.
#' 
#' @name StudentT 
#' 
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of samples to draw from the distribution.
#' @param mu Vector of location values.
#' @param sigma Vector of scale values.
#' @param df Vector of degrees of freedom.
#' @param log,log.p Logical; If \code{TRUE}, values are returned on the log scale.
#' @param lower.tail Logical; If \code{TRUE} (default), return P(X <= x). 
#'   Else, return P(X > x) .
#'   
#' @details See \code{vignette("brms_families")} for details
#' on the parameterization.
#' 
#' @seealso \code{\link[stats:TDist]{TDist}}
#' 
#' @export
dstudent_t <- function(x, df, mu = 0, sigma = 1, log = FALSE) {
  if (any(sigma <= 0)) {
    stop2("sigma must be greater than 0.")
  }
  if (log) {
    dt((x - mu) / sigma, df = df, log = TRUE) - log(sigma)
  } else {
    dt((x - mu) / sigma, df = df) / sigma
  }
}

#' @rdname StudentT
#' @export
pstudent_t <- function(q, df, mu = 0, sigma = 1, 
                       lower.tail = TRUE, log.p = FALSE) {
  if (any(sigma <= 0)) {
    stop2("sigma must be greater than 0.")
  }
  pt((q - mu) / sigma, df = df, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname StudentT
#' @export
qstudent_t <-  function(p, df, mu = 0, sigma = 1) {
  if (any(sigma <= 0)) {
    stop2("sigma must be greater than 0.")
  }
  mu + sigma * qt(p, df = df)
}

#' @rdname StudentT
#' @export
rstudent_t <- function(n, df, mu = 0, sigma = 1) {
  if (any(sigma <= 0)) {
    stop2("sigma must be greater than 0.")
  }
  mu + sigma * rt(n, df = df)
}

#' The Multivariate Normal Distribution
#' 
#' Density function and random generation for the multivariate normal 
#' distribution with mean vector \code{mu} and covariance matrix \code{Sigma}.
#' 
#' @name MultiNormal
#' 
#' @inheritParams StudentT
#' @param x Vector or matrix of quantiles. If \code{x} is a matrix, 
#'   each row is taken to be a quantile.
#' @param mu Mean vector with length equal to the number of dimensions.
#' @param Sigma Covariance matrix.
#' @param check Logical; Indicates whether several input checks
#'   should be performed. Defaults to \code{FALSE} to improve
#'   efficiency.
#'   
#' @details See the Stan user's manual \url{http://mc-stan.org/documentation/}
#' for details on the parameterization
#'   
#' @export
dmulti_normal <- function(x, mu, Sigma, log = FALSE, check = FALSE) {
  if (is.vector(x) || length(dim(x)) == 1L) {
    x <- matrix(x, ncol = length(x))
  }
  p <- ncol(x)
  if (check) {
    if (length(mu) != p) {
      stop2("Dimension of mu is incorrect.")
    }
    if (!all(dim(Sigma) == c(p, p))) {
      stop2("Dimension of Sigma is incorrect.")
    }
    if (!is_symmetric(Sigma)) {
      stop2("Sigma must be a symmetric matrix.")
    }
  }
  chol_Sigma <- chol(Sigma)
  rooti <- backsolve(chol_Sigma, t(x) - mu, transpose = TRUE)
  quads <- colSums(rooti^2)
  out <- -(p / 2) * log(2 * pi) - sum(log(diag(chol_Sigma))) - .5 * quads
  if (!log) {
    out <- exp(out)
  }
  out
}

#' @rdname MultiNormal
#' @export
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
    if (!(is_wholenumber(n) && n > 0)) {
      stop2("n must be a positive integer.")
    }
    if (!all(dim(Sigma) == c(p, p))) {
      stop2("Dimension of Sigma is incorrect.")
    }
    if (!is_symmetric(Sigma)) {
      stop2("Sigma must be a symmetric matrix.")
    }
  }
  samples <- matrix(rnorm(n * p), nrow = n, ncol = p)
  mu + samples %*% chol(Sigma)
}

#' The Multivariate Student-t Distribution
#' 
#' Density function and random generation for the multivariate Student-t
#' distribution with location vector \code{mu}, covariance matrix \code{Sigma},
#' and degrees of freedom \code{df}.
#' 
#' @name MultiStudentT
#' 
#' @inheritParams StudentT
#' @param x Vector or matrix of quantiles. If \code{x} is a matrix, 
#'   each row is taken to be a quantile.
#' @param mu Location vector with length equal to the number of dimensions.
#' @param Sigma Covariance matrix.
#' @param check Logical; Indicates whether several input checks
#'   should be performed. Defaults to \code{FALSE} to improve
#'   efficiency.
#'   
#' @details See the Stan user's manual \url{http://mc-stan.org/documentation/}
#'   for details on the parameterization
#'   
#' @export
dmulti_student_t <- function(x, df, mu, Sigma, log = FALSE, check = FALSE) {
  if (is.vector(x) || length(dim(x)) == 1L) {
    x <- matrix(x, ncol = length(x))
  }
  p <- ncol(x)
  if (check) {
    if (any(df <= 0)) {
      stop2("df must be greater than 0.")
    }
    if (length(mu) != p) {
      stop2("Dimension of mu is incorrect.")
    }
    if (!all(dim(Sigma) == c(p, p))) {
      stop2("Dimension of Sigma is incorrect.")
    }
    if (!is_symmetric(Sigma)) {
      stop2("Sigma must be a symmetric matrix.")
    }
  }
  chol_Sigma <- chol(Sigma)
  rooti <- backsolve(chol_Sigma, t(x) - mu, transpose = TRUE)
  quads <- colSums(rooti^2)
  out <- lgamma((p + df)/2) - (lgamma(df / 2) + sum(log(diag(chol_Sigma))) + 
         p / 2 * log(pi * df)) - 0.5 * (df + p) * log1p(quads / df)
  if (!log) {
    out <- exp(out)
  }
  out
}

#' @rdname MultiStudentT
#' @export
rmulti_student_t <- function(n, df, mu, Sigma, check = FALSE) {
  p <- length(mu)
  if (any(df <= 0)) {
    stop2("df must be greater than 0.")
  }
  samples <- rmulti_normal(n, mu = rep(0, p), Sigma = Sigma, check = check)
  samples <- samples / sqrt(rchisq(n, df = df) / df)
  sweep(samples, 2, mu, "+")
}

#' The von Mises Distribution
#' 
#' Density, distribution function, and random generation for the 
#' von Mises distribution with location \code{mu}, and precision \code{kappa}.
#' 
#' @name VonMises
#' 
#' @inheritParams StudentT
#' @param x,q Vector of quantiles.
#' @param kappa Vector of precision values.
#' @param acc Accuracy of numerical approximations.
#'   
#' @details See \code{vignette("brms_families")} for details
#' on the parameterization.
#' 
#' @export
dvon_mises <- function(x, mu, kappa, log = FALSE) {
  # CircStats::dvm has support within [0, 2*pi], 
  # but in brms we use [-pi, pi]
  if (any(kappa < 0)) {
    stop2("kappa must be non-negative")
  }
  be <- besselI(kappa, nu = 0, expon.scaled = TRUE)
  out <- - log(2 * pi * be) + kappa * (cos(x - mu) - 1)
  if (!log) {
    out <- exp(out)
  }
  out
}

#' @rdname VonMises
#' @export
pvon_mises <- function(q, mu, kappa, lower.tail = TRUE, 
                       log.p = FALSE, acc = 1e-20) {
  # code basis taken from CircStats::pvm but improved 
  # considerably with respect to speed and stability
  if (any(kappa < 0)) {
    stop2("kappa must be non-negative")
  }
  pi <- base::pi
  pi2 <- 2 * pi
  q <- (q + pi) %% pi2
  mu <- (mu + pi) %% pi2
  args <- expand(q = q, mu = mu, kappa = kappa)
  q <- args$q
  mu <- args$mu
  kappa <- args$kappa
  rm(args)
  
  rec_sum <- function(q, kappa, acc, sum = 0, i = 1) {
    # compute the sum of of besselI functions recursively
    term <- (besselI(kappa, nu = i) * sin(i * q)) / i
    sum <- sum + term
    rd <- abs(term) >= acc
    if (sum(rd)) {
      sum[rd] <- rec_sum(
        q[rd], kappa[rd], acc, sum = sum[rd], i = i + 1
      ) 
    }
    sum
  }
  
  .pvon_mises <- function(q, kappa, acc) {
    sum <- rec_sum(q, kappa, acc)
    q / pi2 + sum / (pi * besselI(kappa, nu = 0))
  }
  
  out <- rep(NA, length(mu))
  zero_mu <- mu == 0
  if (sum(zero_mu)) {
    out[zero_mu] <- .pvon_mises(q[zero_mu], kappa[zero_mu], acc)
  }
  lq_mu <- q <= mu
  if (sum(lq_mu)) {
    upper <- (q[lq_mu] - mu[lq_mu]) %% pi2
    upper[upper == 0] <- pi2
    lower <- (-mu[lq_mu]) %% pi2
    out[lq_mu] <- 
      .pvon_mises(upper, kappa[lq_mu], acc) -
      .pvon_mises(lower, kappa[lq_mu], acc)
  }
  uq_mu <- q > mu
  if (sum(uq_mu)) {
    upper <- q[uq_mu] - mu[uq_mu]
    lower <- mu[uq_mu] %% pi2
    out[uq_mu] <- 
      .pvon_mises(upper, kappa[uq_mu], acc) +
      .pvon_mises(lower, kappa[uq_mu], acc)
  }
  if (!lower.tail) {
    out <- 1 - out
  }
  if (log.p) {
    out <- log(out)
  }
  out
}

#' @rdname VonMises
#' @export
rvon_mises <- function(n, mu, kappa) {
  # code basis taken from CircStats::rvm but improved 
  # considerably with respect to speed and stability
  if (any(kappa < 0)) {
    stop2("kappa must be non-negative")
  }
  args <- expand(mu = mu, kappa = kappa, length = n)
  mu <- args$mu
  kappa <- args$kappa
  rm(args)
  pi <- base::pi
  mu <- mu + pi
  
  rvon_mises_outer <- function(r, mu, kappa) {
    n <- length(r)
    U1 <- runif(n, 0, 1)
    z <- cos(pi * U1)
    f <- (1 + r * z) / (r + z)
    c <- kappa * (r - f)
    U2 <- runif(n, 0, 1)
    outer <- is.na(f) | is.infinite(f) |
      !(c * (2 - c) - U2 > 0 | log(c / U2) + 1 - c >= 0)
    inner <- !outer
    out <- rep(NA, n)
    if (sum(inner)) {
      out[inner] <- rvon_mises_inner(f[inner], mu[inner]) 
    }
    if (sum(outer)) {
      # evaluate recursively until a valid sample is found
      out[outer] <- rvon_mises_outer(r[outer], mu[outer], kappa[outer])      
    }
    out
  }
  
  rvon_mises_inner <- function(f, mu) {
    n <- length(f)
    U3 <- runif(n, 0, 1)
    (sign(U3 - 0.5) * acos(f) + mu) %% (2 * pi)
  }
  
  a <- 1 + (1 + 4 * (kappa^2))^0.5
  b <- (a - (2 * a)^0.5) / (2 * kappa)
  r <- (1 + b^2) / (2 * b)
  # indicates underflow due to kappa being close to zero
  is_uf <- is.na(r) | is.infinite(r) 
  not_uf <- !is_uf
  out <- rep(NA, n)
  if (sum(is_uf)) {
    out[is_uf] <- runif(sum(is_uf), 0, 2 * pi) 
  } 
  if (sum(not_uf)) {
    out[not_uf] <- rvon_mises_outer(r[not_uf], mu[not_uf], kappa[not_uf])
  }
  out - pi
}

#' The Exponentially Modified Gaussian Distribution
#' 
#' Density, distribution function, and random generation 
#' for the exponentially modified Gaussian distribution with
#' mean \code{mu} and standard deviation \code{sigma} of the gaussian 
#' component, as well as scale \code{beta} of the exponential
#' component.
#' 
#' @name ExGaussian 
#' 
#' @inheritParams StudentT
#' @param x,q Vector of quantiles.
#' @param mu Vector of means of the gaussian component.
#' @param sigma Vector of standard deviations of the gaussian component.
#' @param beta Vector of scales of the exponential component.
#'   
#' @details See \code{vignette("brms_families")} for details
#' on the parameterization.
#' 
#' @export
dexgaussian <- function(x, mu, sigma, beta, log = FALSE) {
  if (any(sigma <= 0)) {
    stop2("sigma must be greater than 0.")
  }
  if (any(beta <= 0)) {
    stop2("beta must be greater than 0.")
  }
  args <- nlist(x, mu, sigma, beta)
  args <- do.call(expand, args)
  args$z <- with(args, x - mu - sigma^2 / beta)
  
  out <- with(args, ifelse(
    beta > 0.05 * sigma, 
    -log(beta) - (z + sigma^2 / (2 * beta)) / beta + log(pnorm(z / sigma)),
    dnorm(x, mean = mu, sd = sigma, log = TRUE)
  ))
  if (!log) {
    out <- exp(out)
  }
  out
}

#' @rdname ExGaussian
#' @export
pexgaussian <- function(q, mu, sigma, beta, 
                        lower.tail = TRUE, log.p = FALSE) {
  # CDF of the exponentially modified gaussian distribution
  # Args:
  #   see dexgauss
  if (any(sigma <= 0)) {
    stop2("sigma must be greater than 0.")
  }
  if (any(beta <= 0)) {
    stop2("beta must be greater than 0.")
  }
  args <- nlist(q, mu, sigma, beta)
  args <- do.call(expand, args)
  args$z <- with(args, q - mu - sigma^2 / beta)
  
  out <- with(args, ifelse(
    beta > 0.05 * sigma, 
    pnorm((q - mu) / sigma) - pnorm(z / sigma) * 
      exp(((mu + sigma^2 / beta)^2 - mu^2 - 2 * q * sigma^2 / beta) / 
            (2 * sigma^2)), 
    pnorm(q, mean = mu, sd = sigma)
  ))
  if (!lower.tail) {
    out <- 1 - out
  } 
  if (log.p) {
    out <- log(out) 
  } 
  out
}

#' @rdname ExGaussian
#' @export
rexgaussian <- function(n, mu, sigma, beta) {
  # create random numbers of the exgaussian distribution
  # Args:
  #   see dexgauss
  if (any(sigma <= 0)) {
    stop2("sigma must be greater than 0.")
  }
  if (any(beta <= 0)) {
    stop2("beta must be greater than 0.")
  }
  rnorm(n, mean = mu, sd = sigma) + rexp(n, rate = 1 / beta)
}

#' The Frechet Distribution
#' 
#' Density, distribution function, quantile function and random generation 
#' for the Frechet distribution with location \code{loc}, scale \code{scale},
#' and shape \code{shape}.
#' 
#' @name Frechet
#' 
#' @inheritParams StudentT
#' @param x,q Vector of quantiles.
#' @param loc Vector of locations.
#' @param scale Vector of scales.
#' @param shape Vector of shapes.
#'   
#' @details See \code{vignette("brms_families")} for details
#' on the parameterization.
#' 
#' @export
dfrechet <- function (x, loc = 0, scale = 1, shape = 1, log = FALSE) {
  if (isTRUE(any(scale <= 0))) {
    stop2("Argument 'scale' must be positive.")
  }
  if (isTRUE(any(shape <= 0))) {
    stop2("Argument 'shape' must be positive.")
  }
  x <- (x - loc) / scale
  args <- nlist(x, loc, scale, shape)
  args <- do.call(expand, args)
  out <- with(args, 
    log(shape / scale) - (1 + shape) * log(x) - x^(-shape)
  )
  if (!log) { 
    out <- exp(out)
  }
  out
}

#' @rdname Frechet
#' @export
pfrechet <- function(q, loc = 0, scale = 1, shape = 1, 
                     lower.tail = TRUE, log.p = FALSE) {
  if (isTRUE(any(scale <= 0))) {
    stop2("Argument 'scale' must be positive.")
  }
  if (isTRUE(any(shape <= 0))) {
    stop2("Argument 'shape' must be positive.")
  }
  q <- pmax((q - loc) / scale, 0)
  out <- exp(-q^(-shape))
  if (!lower.tail) {
    out <- 1 - out
  }
  if (log.p) {
    out <- log(out)
  }
  out
}

#' @rdname Frechet
#' @export
qfrechet <- function (p, loc = 0, scale = 1, shape = 1, 
                      lower.tail = TRUE, log.p = FALSE) {
  if (isTRUE(any(p <= 0)) || isTRUE(any(p >= 1))) {
    stop("'p' must contain probabilities in (0,1)")
  }
  if (isTRUE(any(scale <= 0))) {
    stop2("Argument 'scale' must be positive.")
  }
  if (isTRUE(any(shape <= 0))) {
    stop2("Argument 'shape' must be positive.")
  }
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  loc + scale * (-log(p))^(-1/shape)
}

#' @rdname Frechet
#' @export
rfrechet <- function(n, loc = 0, scale = 1, shape = 1) {
  if (isTRUE(any(scale <= 0))) {
    stop2("Argument 'scale' must be positive.")
  }
  if (isTRUE(any(shape <= 0))) {
    stop2("Argument 'shape' must be positive.")
  }
  loc + scale * rexp(n)^(-1 / shape)
}

#' The Inverse Gaussian Distribution
#' 
#' Density, distribution function, and random generation 
#' for the inverse Gaussian distribution with location \code{mu}, 
#' and shape \code{shape}.
#' 
#' @name InvGaussian
#' 
#' @inheritParams StudentT
#' @param x,q Vector of quantiles.
#' @param mu Vector of locations.
#' @param shape Vector of shapes.
#'   
#' @details See \code{vignette("brms_families")} for details
#' on the parameterization.
#' 
#' @export
dinv_gaussian <- function(x, mu = 1, shape = 1, log = FALSE) {
  if (isTRUE(any(mu <= 0))) {
    stop2("Argument 'mu' must be positive.")
  }
  if (isTRUE(any(shape <= 0))) {
    stop2("Argument 'shape' must be positive.")
  }
  args <- nlist(x, mu, shape)
  args <- do.call(expand, args)
  out <- with(args,
    0.5 * log(shape / (2 * pi)) -  
    1.5 * log(x) - 0.5 * shape * (x - mu)^2 / (x * mu^2)
  )
  if (!log) {
    out <- exp(out)
  }
  out
}

#' @rdname InvGaussian
#' @export
pinv_gaussian <- function(q, mu = 1, shape = 1, lower.tail = TRUE,
                          log.p = FALSE) {
  if (isTRUE(any(mu <= 0))) {
    stop2("Argument 'mu' must be positive.")
  }
  if (isTRUE(any(shape <= 0))) {
    stop2("Argument 'shape' must be positive.")
  }
  args <- nlist(q, mu, shape)
  args <- do.call(expand, args)
  out <- with(args,
    pnorm(sqrt(shape / q) * (q / mu - 1)) + 
      exp(2 * shape / mu) * pnorm(- sqrt(shape / q) * (q / mu + 1))
  )
  if (!lower.tail) {
    out <- 1 - out
  }
  if (log.p) {
    out <- log(out)
  }
  out
}

#' @rdname InvGaussian
#' @export
rinv_gaussian <- function(n, mu = 1, shape = 1) {
  # create random numbers for the inverse gaussian distribution
  # Args:
  #   Args: see dinv_gaussian
  if (isTRUE(any(mu <= 0))) {
    stop2("Argument 'mu' must be positive.")
  }
  if (isTRUE(any(shape <= 0))) {
    stop2("Argument 'shape' must be positive.")
  }
  args <- nlist(mu, shape, length = n)
  args <- do.call(expand, args)
  # algorithm from wikipedia
  args$y <- rnorm(n)^2
  args$x <- with(args, 
    mu + (mu^2 * y) / (2 * shape) - mu / (2 * shape) * 
      sqrt(4 * mu * shape * y + mu^2 * y^2) 
  )
  args$z <- runif(n)
  with(args, ifelse(z <= mu / (mu + x), x, mu^2 / x))
}

#' The Generalized Extreme Value Distribution
#' 
#' Density, distribution function, and random generation 
#' for the generalized extreme value distribution with 
#' location \code{mu}, scale \code{sigma} and shape \code{xi}.
#' 
#' @name GenExtremeValue
#' 
#' @inheritParams StudentT
#' @param x,q Vector of quantiles.
#' @param mu Vector of locations.
#' @param sigma Vector of scales.
#' @param xi Vector of shapes.
#'   
#' @details See \code{vignette("brms_families")} for details
#' on the parameterization.
#' 
#' @export
dgen_extreme_value <- function(x, mu = 0, sigma = 1, 
                               xi = 0, log = FALSE) {
  if (any(sigma <= 0)) {
    stop2("sigma bust be greater than 0.")
  }
  x <- (x - mu) / sigma
  args <- nlist(x, mu, sigma, xi)
  args <- do.call(expand, args)
  args$t <- with(args, 1 + xi * x)
  out <- with(args, ifelse(
    xi == 0, 
    - log(sigma) - x - exp(-x),
    - log(sigma) - (1 + 1 / xi) * log(t) - t^(-1 / xi)
  ))
  if (!log) {
    out <- exp(out)
  } 
  out
}

#' @rdname GenExtremeValue
#' @export
pgen_extreme_value <- function(q, mu = 0, sigma = 1, xi = 0,
                               lower.tail = TRUE, log.p = FALSE) {
  if (any(sigma <= 0)) {
    stop2("sigma bust be greater than 0.")
  }
  q <- (q - mu) / sigma
  args <- nlist(q, mu, sigma, xi)
  args <- do.call(expand, args)
  out <- with(args, ifelse(
    xi == 0, 
    exp(-exp(-q)),
    exp(-(1 + xi * q)^(-1 / xi))
  ))
  if (!lower.tail) {
    out <- 1 - out
  }
  if (log.p) {
    out <- log(out)
  }
  out
}

#' @rdname GenExtremeValue
#' @export
rgen_extreme_value <- function(n, mu = 0, sigma = 1, xi = 0) {
  if (any(sigma <= 0)) {
    stop2("sigma bust be greater than 0.")
  }
  args <- nlist(mu, sigma, xi, length = n)
  args <- do.call(expand, args)
  with(args, ifelse(
    xi == 0,
    mu - sigma * log(rexp(n)),
    mu + sigma * (rexp(n)^(-xi) - 1) / xi
  ))
}

#' The Asymmetric Laplace Distribution
#' 
#' Density, distribution function, quantile function and random generation 
#' for the asymmetric Laplace distribution with location \code{mu}, 
#' scale \code{sigma} and asymmetry parameter \code{quantile}.
#' 
#' @name AsymLaplace
#' 
#' @inheritParams StudentT
#' @param x,q Vector of quantiles.
#' @param mu Vector of locations.
#' @param sigma Vector of scales.
#' @param quantile Asymmetry parameter corresponding to quantiles
#'   in quantile regression (hence the name).
#'   
#' @details See \code{vignette("brms_families")} for details
#' on the parameterization.
#' 
#' @export
dasym_laplace <- function(x, mu = 0, sigma = 1, quantile = 0.5, 
                          log = FALSE) {
  out <- ifelse(x < mu, 
    yes = (quantile * (1 - quantile) / sigma) * 
           exp((1 - quantile) * (x - mu) / sigma),
    no = (quantile * (1 - quantile) / sigma) * 
          exp(-quantile * (x - mu) / sigma)
  )
  if (log) {
    out <- log(out)
  }
  out
}

#' @rdname AsymLaplace
#' @export
pasym_laplace <- function(q, mu = 0, sigma = 1, quantile = 0.5,
                          lower.tail = TRUE, log.p = FALSE) {
  out <- ifelse(q < mu, 
    yes = quantile * exp((1 - quantile) * (q - mu) / sigma), 
    no = 1 - (1 - quantile) * exp(-quantile * (q - mu) / sigma)
  )
  if (!lower.tail) {
    out <- 1 - out
  }
  if (log.p) {
    out <- log(out) 
  }
  out
}

#' @rdname AsymLaplace
#' @export
qasym_laplace <- function(p, mu = 0, sigma = 1, quantile = 0.5,
                          lower.tail = TRUE, log.p = FALSE) {
  # quantile function of the asymmetric laplace distribution
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  if (length(quantile) == 1L) {
    quantile <- rep(quantile, length(mu))
  }
  ifelse(p < quantile, 
    yes = mu + ((sigma * log(p / quantile)) / (1 - quantile)), 
    no = mu - ((sigma * log((1 - p) / (1 - quantile))) / quantile)
  )
}

#' @rdname AsymLaplace
#' @export
rasym_laplace <- function(n, mu = 0, sigma = 1, quantile = 0.5) {
  # random numbers of the asymmetric laplace distribution
  u <- runif(n)
  qasym_laplace(u, mu = mu, sigma = sigma, quantile = quantile)
}

#' The Wiener Diffusion Model Distribution
#' 
#' Density function and random generation for the Wiener
#' diffusion model distribution with boundary separation \code{alpha},
#' non-decision time \code{tau}, bias  \code{beta} and 
#' drift rate \code{delta}.
#' 
#' @name Wiener
#' 
#' @inheritParams StudentT
#' @param alpha Boundary separation parameter.
#' @param tau Non-decision time parameter.
#' @param beta Bias parameter.
#' @param delta Drift rate parameter.
#' @param resp Response: \code{"upper"} or \code{"lower"}. 
#'   If no character vector, it is coerced to logical 
#'   where \code{TRUE} indicates \code{"upper"} and 
#'   \code{FALSE} indicates \code{"lower"}.
#' @param types Which types of responses to return? By default,
#'   return both the response times \code{"q"} and the dichotomous 
#'   responses \code{"resp"}. If either \code{"q"} or \code{"resp"}, 
#'   return only one of the two types.
#'   
#' @details These are wrappers around functions of the \pkg{RWiener} package.  
#' See \code{vignette("brms_families")} for details on the parameterization.
#' 
#' @seealso \code{\link[RWiener:wienerdist]{wienerdist}}
#' 
#' @export
dwiener <- function(x, alpha, tau, beta, delta, resp = 1, log = FALSE) {
  alpha <- as.numeric(alpha)
  tau <- as.numeric(tau)
  beta <- as.numeric(beta)
  delta <- as.numeric(delta)
  if (!is.character(resp)) {
    resp <- ifelse(resp, "upper", "lower") 
  }
  # vectorized version of RWiener::dwiener
  .dwiener <- Vectorize(
    RWiener::dwiener, 
    c("alpha", "tau", "beta", "delta")
  )
  args <- nlist(q = x, alpha, tau, beta, delta, resp, give_log = log)
  do.call(.dwiener, args)
}

#' @rdname Wiener
#' @export
rwiener <- function(n, alpha, tau, beta, delta, types = c("q", "resp")) {
  stopifnot(all(types %in% c("q", "resp")))
  alpha <- as.numeric(alpha)
  tau <- as.numeric(tau)
  beta <- as.numeric(beta)
  delta <- as.numeric(delta)
  max_len <- max(lengths(list(alpha, tau, beta, delta)))
  n <- n[1]
  if (max_len > 1L) {
    if (!n %in% c(1, max_len)) {
      stop2("Can only sample exactly once for each condition.")
    }
    n <- 1
  }
  .rwiener <- function(...) {
    # vectorized version of RWiener::rwiener
    # returns a numeric vector
    fun <- Vectorize(
      rwiener_num, 
      c("alpha", "tau", "beta", "delta"),
      SIMPLIFY = FALSE
    )
    do.call(rbind, fun(...))
  }
  args <- nlist(n, alpha, tau, beta, delta, types)
  do.call(.rwiener, args)
}

rwiener_num <- function(n, alpha, tau, beta, delta, types) {
  # helper function to return a numeric vector instead
  # of a data.frame with two columns as for RWiener::rwiener
  out <- RWiener::rwiener(n, alpha, tau, beta, delta)
  out[["resp"]] <- ifelse(out[["resp"]] == "upper", 1, 0)
  if (length(types) == 1L) {
    out <- out[[types]]
  }
  out
}

dcategorical <- function(x, eta, ncat, link = "logit") {
  # density of the categorical distribution
  # Args:
  #   x: positive integers not greater than ncat
  #   eta: the linear predictor (of length or ncol ncat-1)  
  #   ncat: the number of categories
  #   link: the link function
  # Returns:
  #   probabilities P(X = x)
  if (is.null(dim(eta))) {
    eta <- matrix(eta, nrow = 1)
  }
  if (length(dim(eta)) != 2L) {
    stop2("eta must be a numeric vector or matrix.")
  }
  if (missing(ncat)) {
    ncat <- ncol(eta) + 1
  }
  if (link == "logit") {
    p <- exp(cbind(rep(0, nrow(eta)), eta[, 1:(ncat - 1)]))
  } else {
    stop2("Link '", link, "' not supported.")
  }
  p <- p / rowSums(p)
  p[, x]
}

pcategorical <- function(q, eta, ncat, link = "logit") {
  # distribution functions for the categorical family
  # Args:
  #   q: positive integers not greater than ncat
  #   eta: the linear predictor (of length or ncol ncat-1)  
  #   ncat: the number of categories
  #   link: a character string naming the link
  # Retruns: 
  #   probabilities P(x <= q)
  p <- dcategorical(1:max(q), eta = eta, ncat = ncat, link = link)
  do.call(cbind, lapply(q, function(j) rowSums(as.matrix(p[, 1:j]))))
}

dcumulative <- function(x, eta, ncat, link = "logit") {
  # density of the cumulative distribution
  # Args: same as dcategorical
  if (is.null(dim(eta))) {
    eta <- matrix(eta, nrow = 1)
  }
  if (length(dim(eta)) != 2) {
    stop2("eta must be a numeric vector or matrix.")
  }
  if (missing(ncat)) {
    ncat <- ncol(eta) + 1
  }
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
  # Args: same as dcategorical
  if (is.null(dim(eta))) {
    eta <- matrix(eta, nrow = 1)
  }
  if (length(dim(eta)) != 2) {
    stop2("eta must be a numeric vector or matrix.")
  }
  if (missing(ncat)) {
    ncat <- ncol(eta) + 1
  }
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
  # Args: same as dcategorical
  if (is.null(dim(eta))) {
    eta <- matrix(eta, nrow = 1)
  }
  if (length(dim(eta)) != 2) {
    stop2("eta must be a numeric vector or matrix.")
  }
  if (missing(ncat)) {
    ncat <- ncol(eta) + 1
  }
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
  # Args: same as dcategorical
  if (is.null(dim(eta))) {
    eta <- matrix(eta, nrow = 1)
  }
  if (length(dim(eta)) != 2) {
    stop2("eta must be a numeric vector or matrix.")
  }
  if (missing(ncat)) {
    ncat <- ncol(eta) + 1
  }
  if (link == "logit") { 
    # faster evaluation in this case
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
  # Args:
  #   q: positive integers not greater than ncat
  #   eta: the linear predictor (of length or ncol ncat-1)  
  #   ncat: the number of categories
  #   family: a character string naming the family
  #   link: a character string naming the link
  # Returns: 
  #   probabilites P(x <= q)
  args <- list(1:max(q), eta = eta, ncat = ncat, link = link)
  p <- do.call(paste0("d", family), args)
  .fun <- function(j) {
    rowSums(as.matrix(p[, 1:j]))
  }
  do.call(cbind, lapply(q, .fun))
}
