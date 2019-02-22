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

#' The Skew-Normal Distribution
#' 
#' Density, distribution function, and random generation for the 
#' skew-normal distribution with mean \code{mu}, 
#' standard deviation \code{sigma}, and skewness \code{alpha}.
#' 
#' @name SkewNormal
#' 
#' @inheritParams StudentT
#' @param x,q Vector of quantiles.
#' @param mu Vector of mean values.
#' @param sigma Vector of standard deviation values.
#' @param alpha Vector of skewness values.
#' @param xi Optional vector of location values.
#'   If \code{NULL} (the default), will be computed internally.
#' @param omega Optional vector of scale values.
#'   If \code{NULL} (the default), will be computed internally. 
#' @param tol Tolerance of the approximation used in the 
#'   computation of quantiles.
#'   
#' @details See \code{vignette("brms_families")} for details
#' on the parameterization.
#' 
#' @export
dskew_normal <- function(x, mu = 0, sigma = 1, alpha = 0, 
                         xi = NULL, omega = NULL, log = FALSE) {
  if (any(sigma <= 0)) {
    stop2("sigma must be greater than 0.")
  }
  args <- cp2dp(mu, sigma, alpha, xi = xi, omega = omega, x = x)
  out <- with(args, {
    # do it like sn::dsn
    z <- (x - xi) / omega 
    if (length(alpha) == 1L) {
      alpha <- rep(alpha, length(z))
    }
    logN <- -log(sqrt(2 * pi)) - log(omega) - z^2 / 2
    logS <- ifelse(
      abs(alpha) < Inf, 
      pnorm(alpha * z, log.p = TRUE),
      log(as.numeric(sign(alpha) * z > 0))
    )
    out <- logN + logS - pnorm(0, log.p = TRUE)
    ifelse(abs(z) == Inf, -Inf, out)
  })
  if (!log) {
    out <- exp(out)
  }
  out
}

#' @rdname SkewNormal
#' @export 
pskew_normal <- function(q, mu = 0, sigma = 1, alpha = 0, 
                         xi = NULL, omega = NULL,
                         lower.tail = TRUE, log.p = FALSE) {
  require_package("mnormt")
  if (any(sigma <= 0)) {
    stop2("sigma must be greater than 0.")
  }
  args <- cp2dp(mu, sigma, alpha, xi = xi, omega = omega, q = q)
  out <- with(args, {
    # do it like sn::psn
    z <- (q - xi) / omega
    nz <- length(z)
    is_alpha_inf <- abs(alpha) == Inf
    delta[is_alpha_inf] <- sign(alpha[is_alpha_inf])
    out <- numeric(nz)
    for (k in seq_len(nz)) {
      if (is_alpha_inf[k]) {
        if (alpha[k] > 0) {
          out[k] <- 2 * (pnorm(pmax(z[k], 0)) - 0.5)
        } else {
          out[k] <- 1 - 2 * (0.5 - pnorm(pmin(z[k], 0)))
        }
      } else {
        S <- matrix(c(1, -delta[k], -delta[k], 1), 2, 2)
        out[k] <- 2 * mnormt::biv.nt.prob(
          0, lower = rep(-Inf, 2), upper = c(z[k], 0),
          mean = c(0, 0), S = S
        )
      }
    }
    pmin(1, pmax(0, out))
  })
  if (!lower.tail) {
    out <- 1 - out
  }
  if (log.p) {
    out <- log(out)
  }
  out
}

#' @rdname SkewNormal
#' @export
qskew_normal <- function(p, mu = 0, sigma = 1, alpha = 0, 
                         xi = NULL, omega = NULL,
                         lower.tail = TRUE, log.p = FALSE, 
                         tol = 1e-8) {
  if (any(sigma <= 0)) {
    stop2("sigma must be greater than 0.")
  }
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  args <- cp2dp(mu, sigma, alpha, xi = xi, omega = omega, p = p)
  out <- with(args, {
    # do it like sn::qsn
    na <- is.na(p) | (p < 0) | (p > 1)
    zero <- (p == 0)
    one <- (p == 1)
    p <- replace(p, (na | zero | one), 0.5)
    cum <- skew_normal_cumulants(0, 1, alpha, n = 4)
    g1 <- cum[, 3] / cum[, 2]^(3 / 2)
    g2 <- cum[, 4] / cum[, 2]^2
    x <- qnorm(p)
    x <- x + (x^2 - 1) * g1 / 6 + 
      x * (x^2 - 3) * g2 / 24 - 
      x * (2 * x^2 - 5) * g1^2 / 36
    x <- cum[, 1] + sqrt(cum[, 2]) * x
    px <- pskew_normal(x, xi = 0, omega = 1, alpha = alpha)
    max_err <- 1
    while (max_err > tol) {
      x1 <- x - (px - p) / 
        dskew_normal(x, xi = 0, omega = 1, alpha = alpha)
      x <- x1
      px <- pskew_normal(x, xi = 0, omega = 1, alpha = alpha)
      max_err <- max(abs(px - p))
      if (is.na(max_err)) {
        warning2("Approximation in 'qskew_normal' might have failed.")
      }
    }
    x <- replace(x, na, NA)
    x <- replace(x, zero, -Inf)
    x <- replace(x, one, Inf)
    as.numeric(xi + omega * x)
  })
  out
}

#' @rdname SkewNormal
#' @export
rskew_normal <- function(n, mu = 0, sigma = 1, alpha = 0,
                         xi = NULL, omega = NULL) {
  if (any(sigma <= 0)) {
    stop2("sigma must be greater than 0.")
  }
  args <- cp2dp(mu, sigma, alpha, xi = xi, omega = omega)
  with(args, {
    # do it like sn::rsn
    z1 <- rnorm(n)
    z2 <- rnorm(n)
    id <- z2 > args$alpha * z1
    z1[id] <- -z1[id]
    xi + omega * z1 
  })
}

cp2dp <- function(mu = 0, sigma = 1, alpha = 0, 
                  xi = NULL, omega = NULL, ...) {
  # convert skew-normal mixed-CP to DP parameterization
  # Returns:
  #   A data.frame containing all relevant parameters 
  delta <- alpha / sqrt(1 + alpha^2)
  if (is.null(omega)) {
    omega <- sigma / sqrt(1 - 2 / pi * delta^2)
  } 
  if (is.null(xi)) {
    xi <- mu - omega * delta * sqrt(2 / pi)
  }
  expand(dots = nlist(mu, sigma, alpha, xi, omega, delta, ...))
}

skew_normal_cumulants <- function(xi = 0, omega = 1, alpha = 0, n = 4) {
  # helper function for qskew_normal 
  # code basis taken from sn::sn.cumulants
  # this function uses xi and omega rather than mu and sigma!
  cumulants_half_norm <- function(n) {
    n <- max(n, 2)
    n <- as.integer(2 * ceiling(n/2))
    half.n <- as.integer(n/2)
    m <- 0:(half.n - 1)
    a <- sqrt(2/pi)/(gamma(m + 1) * 2^m * (2 * m + 1))
    signs <- rep(c(1, -1), half.n)[seq_len(half.n)]
    a <- as.vector(rbind(signs * a, rep(0, half.n)))
    coeff <- rep(a[1], n)
    for (k in 2:n) {
      ind <- seq_len(k - 1)
      coeff[k] <- a[k] - sum(ind * coeff[ind] * a[rev(ind)]/k)
    }
    kappa <- coeff * gamma(seq_len(n) + 1)
    kappa[2] <- 1 + kappa[2]
    return(kappa)
  }
  
  args <- expand(dots = nlist(xi, omega, alpha))
  with(args, {
    # do it like sn::sn.cumulants
    delta <- alpha / sqrt(1 + alpha^2)
    kv <- cumulants_half_norm(n)
    if (length(kv) > n)  {
      kv <- kv[-(n + 1)] 
    }
    kv[2] <- kv[2] - 1
    kappa <- outer(delta, 1:n, "^") * 
      matrix(rep(kv, length(xi)), ncol = n, byrow = TRUE)
    kappa[, 2] <- kappa[, 2] + 1
    kappa <- kappa * outer(omega, 1:n, "^")
    kappa[, 1] <- kappa[, 1] + xi
    kappa
  })
}

pinvgamma <- function(q, shape, rate, lower.tail = TRUE, log.p = FALSE) {
  # CDF of the inverse gamma function
  pgamma(1/q, shape, rate = rate, lower.tail = !lower.tail, log.p = log.p)
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
  out <- -log(2 * pi * be) + kappa * (cos(x - mu) - 1)
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
#' @param mu Vector of means of the combined distribution.
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
  args <- do_call(expand, args)
  args$mu <- with(args, mu - beta)
  args$z <- with(args, x - mu - sigma^2 / beta)
  
  out <- with(args, 
    -log(beta) - (z + sigma^2 / (2 * beta)) / beta + 
      pnorm(z / sigma, log.p = TRUE)
  )
  if (!log) {
    out <- exp(out)
  }
  out
}

#' @rdname ExGaussian
#' @export
pexgaussian <- function(q, mu, sigma, beta, 
                        lower.tail = TRUE, log.p = FALSE) {
  if (any(sigma <= 0)) {
    stop2("sigma must be greater than 0.")
  }
  if (any(beta <= 0)) {
    stop2("beta must be greater than 0.")
  }
  args <- nlist(q, mu, sigma, beta)
  args <- do_call(expand, args)
  args$mu <- with(args, mu - beta)
  args$z <- with(args, q - mu - sigma^2 / beta)
  
  out <- with(args, 
    pnorm((q - mu) / sigma) - pnorm(z / sigma) * 
      exp(((mu + sigma^2 / beta)^2 - mu^2 - 2 * q * sigma^2 / beta) / 
            (2 * sigma^2))
  )
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
  if (any(sigma <= 0)) {
    stop2("sigma must be greater than 0.")
  }
  if (any(beta <= 0)) {
    stop2("beta must be greater than 0.")
  }
  mu <- mu - beta
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
  args <- do_call(expand, args)
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

#' The Shifted Log Normal Distribution
#' 
#' Density, distribution function, quantile function and random generation 
#' for the shifted log normal distribution with mean \code{meanlog}, 
#' standard deviation \code{sdlog}, and shift parameter \code{shift}.
#' 
#' @name Shifted_Lognormal
#' 
#' @inheritParams StudentT
#' @param x,q Vector of quantiles.
#' @param meanlog Vector of means.
#' @param sdlog Vector of standard deviations.
#' @param shift Vector of shifts.
#'   
#' @details See \code{vignette("brms_families")} for details
#' on the parameterization.
#' 
#' @export
dshifted_lnorm <- function(x, meanlog = 0, sdlog = 1, shift = 0, log = FALSE) {
  args <- nlist(dist = "lnorm", x, shift, meanlog, sdlog, log)
  do_call(dshifted, args)
}

#' @rdname Shifted_Lognormal
#' @export
pshifted_lnorm <- function(q, meanlog = 0, sdlog = 1, shift = 0, 
                           lower.tail = TRUE, log.p = FALSE) {
  args <- nlist(dist = "lnorm", q, shift, meanlog, sdlog, lower.tail, log.p)
  do_call(pshifted, args)
}

#' @rdname Shifted_Lognormal
#' @export
qshifted_lnorm <- function(p, meanlog = 0, sdlog = 1, shift = 0, 
                           lower.tail = TRUE, log.p = FALSE) {
  args <- nlist(dist = "lnorm", p, shift, meanlog, sdlog, lower.tail, log.p)
  do_call(qshifted, args)
}

#' @rdname Shifted_Lognormal
#' @export
rshifted_lnorm <- function(n, meanlog = 0, sdlog = 1, shift = 0) {
  args <- nlist(dist = "lnorm", n, shift, meanlog, sdlog)
  do_call(rshifted, args)
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
  args <- do_call(expand, args)
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
  args <- do_call(expand, args)
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
  args <- do_call(expand, args)
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
  args <- do_call(expand, args)
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
  args <- do_call(expand, args)
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
  args <- do_call(expand, args)
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
  u <- runif(n)
  qasym_laplace(u, mu = mu, sigma = sigma, quantile = quantile)
}

#' The Discrete Weibull Distribution
#' 
#' Density, distribution function, quantile function and random generation 
#' for the discrete Weibull distribution with location \code{mu} and
#' shape \code{shape}.
#' 
#' @name DiscreteWeibull
#' 
#' @inheritParams StudentT
#' @param mu Location parameter in the unit interval.
#' @param shape Positive shape parameter.
#' 
#' @details See \code{vignette("brms_families")} for details
#' on the parameterization.
#' 
#' @export
ddiscrete_weibull <- function(x, mu, shape, log = FALSE) {
  if (any(mu < 0 | mu > 1)) {
    stop2("mu bust be between 0 and 1.")
  }
  if (any(shape <= 0)) {
    stop2("shape bust be greater than 0.")
  }
  x <- round(x)
  out <- mu^x^shape - mu^(x + 1)^shape
  out[x < 0] <- 0
  if (log) {
    out <- log(out)
  }
  out
}

#' @rdname DiscreteWeibull
#' @export
pdiscrete_weibull <- function(x, mu, shape, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu < 0 | mu > 1)) {
    stop2("mu bust be between 0 and 1.")
  }
  if (any(shape <= 0)) {
    stop2("shape bust be greater than 0.")
  }
  x <- round(x)
  if (lower.tail) {
    out <- 1 - mu^(x + 1)^shape
    out[x < 0] <- 0
  } else {
    out <- mu^(x + 1)^shape
    out[x < 0] <- 1
  }
  if (log.p) {
    out <- log(out)
  }
  out
}

#' @rdname DiscreteWeibull
#' @export
qdiscrete_weibull <- function(p, mu, shape, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu < 0 | mu > 1)) {
    stop2("mu bust be between 0 and 1.")
  }
  if (any(shape <= 0)) {
    stop2("shape bust be greater than 0.")
  }
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  ceiling((log(1 - p) / log(mu))^(1 / shape) - 1)
}

#' @rdname DiscreteWeibull
#' @export
rdiscrete_weibull <- function(n, mu, shape) {
  r <- runif(n, 0, 1)
  qdiscrete_weibull(r, mu, shape)
}

#' The Dirichlet Distribution
#' 
#' Density function and random number generation for the dirichlet
#' distribution with shape parameter vector \code{alpha}.
#' 
#' @name Dirichlet
#' 
#' @inheritParams StudentT
#' @param x Matrix of quantiles. Each row corresponds to one probability vector.
#' @param alpha Matrix of positive shape parameters. Each row corresponds to one
#'   probability vector.
#'
#' @details See \code{vignette("brms_families")} for details on the 
#' parameterization.
#' 
#' @export
ddirichlet <- function(x, alpha, log = FALSE) {
  log <- as_one_logical(log)
  if (!is.matrix(x)) {
    x <- matrix(x, nrow = 1)
  }
  if (!is.matrix(alpha)) {
    alpha <- matrix(alpha, nrow(x), length(alpha), byrow = TRUE)
  }
  if (nrow(x) == 1L && nrow(alpha) > 1L) {
    x <- repl(x, nrow(alpha))
    x <- do_call(rbind, x)
  } else if (nrow(x) > 1L && nrow(alpha) == 1L) {
    alpha <- repl(alpha, nrow(x))
    alpha <- do_call(rbind, alpha)
  }
  if (any(x < 0)) {
    stop2("x must be non-negative.")
  }
  if (!is_equal(rowSums(x), rep(1, nrow(x)))) {
    stop2("x must sum to 1 per row.")
  }
  if (any(alpha <= 0)) {
    stop2("alpha must be positive.")
  }
  out <- lgamma(rowSums(alpha)) - rowSums(lgamma(alpha)) + 
    rowSums((alpha - 1) * log(x))
  if (!log) {
    out <- exp(out)
  } 
  return(out)
}

#' @rdname Dirichlet
#' @export
rdirichlet <- function(n, alpha) {
  n <- as_one_numeric(n)
  if (!is.matrix(alpha)) {
    alpha <- matrix(alpha, nrow = 1)
  }
  if (prod(dim(alpha)) == 0) {
    stop2("alpha should be non-empty.")
  }
  if (any(alpha <= 0)) {
    stop2("alpha must be positive.")
  }
  if (n == 1) {
    n <- nrow(alpha)
  }
  if (n > nrow(alpha)) {
    alpha <- matrix(alpha, nrow = n, ncol = ncol(alpha), byrow = TRUE)
  }
  x <- matrix(rgamma(ncol(alpha) * n, alpha), ncol = ncol(alpha))
  x / rowSums(x)
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
  require_package("RWiener")
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
  do_call(.dwiener, args)
}

#' @rdname Wiener
#' @export
rwiener <- function(n, alpha, tau, beta, delta, types = c("q", "resp")) {
  require_package("RWiener")
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
    do_call(rbind, fun(...))
  }
  args <- nlist(n, alpha, tau, beta, delta, types)
  do_call(.rwiener, args)
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

#' Zero-Inflated Distributions
#' 
#' Density and distribution functions for zero-inflated distributions.
#' 
#' @name ZeroInflated
#' 
#' @inheritParams StudentT
#' @param zi zero-inflation propability
#' @param mu,lambda location parameter
#' @param shape,shape1,shape2 shape parameter
#' @param size number of trials
#' @param prob probability of success on each trial
#' 
#' @details 
#' The density of a zero-inflated distribution can be specified as follows.
#' If \eqn{x = 0} set \eqn{f(x) = \theta + (1 - \theta) * g(0)}.
#' Else set \eqn{f(x) = (1 - \theta) * g(x)},
#' where \eqn{g(x)} is the density of the non-zero-inflated part.
NULL

#' @rdname ZeroInflated
#' @export
dzero_inflated_poisson <- function(x, lambda, zi, log = FALSE) {
  pars <- nlist(lambda)
  .dzero_inflated(x, "pois", zi, pars, log)
}

#' @rdname ZeroInflated
#' @export
pzero_inflated_poisson <- function(q, lambda, zi, lower.tail = TRUE, 
                                   log.p = FALSE) {
  pars <- nlist(lambda)
  .pzero_inflated(q, "pois", zi, pars, lower.tail, log.p)
}

#' @rdname ZeroInflated
#' @export
dzero_inflated_negbinomial <- function(x, mu, shape, zi, log = FALSE) {
  pars <- nlist(mu, size = shape)
  .dzero_inflated(x, "nbinom", zi, pars, log)
}

#' @rdname ZeroInflated
#' @export
pzero_inflated_negbinomial <- function(q, mu, shape, zi, lower.tail = TRUE,
                                       log.p = FALSE) {
  pars <- nlist(mu, size = shape)
  .pzero_inflated(q, "nbinom", zi, pars, lower.tail, log.p)
}

#' @rdname ZeroInflated
#' @export
dzero_inflated_binomial <- function(x, size, prob, zi, log = FALSE) {
  pars <- nlist(size, prob)
  .dzero_inflated(x, "binom", zi, pars, log)
}

#' @rdname ZeroInflated
#' @export
pzero_inflated_binomial <- function(q, size, prob, zi, lower.tail = TRUE,
                                    log.p = FALSE) {
  pars <- nlist(size, prob)
  .pzero_inflated(q, "binom", zi, pars, lower.tail, log.p)
}

#' @rdname ZeroInflated
#' @export
dzero_inflated_beta <- function(x, shape1, shape2, zi, log = FALSE) {
  pars <- nlist(shape1, shape2)
  # zi_beta is technically a hurdle model
  .dhurdle(x, "beta", zi, pars, log, type = "real")
}

#' @rdname ZeroInflated
#' @export
pzero_inflated_beta <- function(q, shape1, shape2, zi, lower.tail = TRUE, 
                                log.p = FALSE) {
  pars <- nlist(shape1, shape2)
  # zi_beta is technically a hurdle model
  .phurdle(q, "beta", zi, pars, lower.tail, log.p, type = "real")
}

.dzero_inflated <- function(x, dist, zi, pars, log) {
  # density function of zero-inflated models
  # Args:
  #   dist: name of the distribution
  #   zi: bernoulli zero-inflated parameter
  #   pars: list of parameters passed to pdf
  stopifnot(is.list(pars))
  dist <- as_one_character(dist)
  log <- as_one_logical(log)
  args <- expand(dots = c(nlist(x, zi), pars))
  x <- args$x
  zi <- args$zi
  pars <- args[names(pars)]
  pdf <- paste0("d", dist) 
  out <- ifelse(x == 0, 
    log(zi + (1 - zi) * do_call(pdf, c(0, pars))),
    log(1 - zi) + do_call(pdf, c(list(x), pars, log = TRUE))
  )
  if (!log) {
    out <- exp(out)
  }
  out
}

.pzero_inflated <- function(q, dist, zi, pars, lower.tail, log.p) {
  # distribution function of zero-inflated models
  # Args:
  #   dist: name of the distribution
  #   zi: bernoulli zero-inflated parameter
  #   pars: list of parameters passed to pdf
  stopifnot(is.list(pars))
  dist <- as_one_character(dist)
  lower.tail <- as_one_logical(lower.tail)
  log.p <- as_one_logical(log.p)
  args <- expand(dots = c(nlist(q, zi), pars))
  q <- args$q
  zi <- args$zi
  pars <- args[names(pars)]
  cdf <- paste0("p", dist)
  out <- log(1 - zi) +
    do_call(cdf, c(list(q), pars, lower.tail = FALSE, log = TRUE))
  if (lower.tail) {
    out <- 1 - exp(out)
    if (log.p) {
      out <- log(out)
    }
  } else {
    if (!log.p) {
      out <- exp(out)
    }
  }
  out
}

#' Hurdle Distributions
#' 
#' Density and distribution functions for hurdle distributions.
#' 
#' @name Hurdle
#' 
#' @inheritParams StudentT
#' @param hu hurdle propability
#' @param mu,lambda location parameter
#' @param shape shape parameter
#' @param sigma,scale scale parameter
#' 
#' @details 
#' The density of a hurdle distribution can be specified as follows.
#' If \eqn{x = 0} set \eqn{f(x) = \theta}. Else set
#' \eqn{f(x) = (1 - \theta) * g(x) / (1 - G(0))}
#' where \eqn{g(x)} and \eqn{G(x)} are the density and distribution
#' function of the non-hurdle part, respectively.
NULL

#' @rdname Hurdle
#' @export
dhurdle_poisson <- function(x, lambda, hu, log = FALSE) {
  pars <- nlist(lambda)
  .dhurdle(x, "pois", hu, pars, log, type = "int")
}

#' @rdname Hurdle
#' @export
phurdle_poisson <- function(q, lambda, hu, lower.tail = TRUE, 
                            log.p = FALSE) {
  pars <- nlist(lambda)
  .phurdle(q, "pois", hu, pars, lower.tail, log.p, type = "int")
}

#' @rdname Hurdle
#' @export
dhurdle_negbinomial <- function(x, mu, shape, hu, log = FALSE) {
  pars <- nlist(mu, size = shape)
  .dhurdle(x, "nbinom", hu, pars, log, type = "int")
}

#' @rdname Hurdle
#' @export
phurdle_negbinomial <- function(q, mu, shape, hu, lower.tail = TRUE, 
                                log.p = FALSE) {
  pars <- nlist(mu, size = shape)
  .phurdle(q, "nbinom", hu, pars, lower.tail, log.p, type = "int")
}

#' @rdname Hurdle
#' @export
dhurdle_gamma <- function(x, shape, scale, hu, log = FALSE) {
  pars <- nlist(shape, scale)
  .dhurdle(x, "gamma", hu, pars, log, type = "real")
}

#' @rdname Hurdle
#' @export
phurdle_gamma <- function(q, shape, scale, hu, lower.tail = TRUE, 
                          log.p = FALSE) {
  pars <- nlist(shape, scale)
  .phurdle(q, "gamma", hu, pars, lower.tail, log.p, type = "real")
}

#' @rdname Hurdle
#' @export
dhurdle_lognormal <- function(x, mu, sigma, hu, log = FALSE) {
  pars <- list(meanlog = mu, sdlog = sigma)
  .dhurdle(x, "lnorm", hu, pars, log, type = "real")
}

#' @rdname Hurdle
#' @export
phurdle_lognormal <- function(q, mu, sigma, hu, lower.tail = TRUE,
                              log.p = FALSE) {
  pars <- list(meanlog = mu, sdlog = sigma)
  .phurdle(q, "lnorm", hu, pars, lower.tail, log.p, type = "real")
}

.dhurdle <- function(x, dist, hu, pars, log, type) {
  # density function of hurdle models
  # Args:
  #   dist: name of the distribution
  #   hu: bernoulli hurdle parameter
  #   pars: list of parameters passed to pdf
  #   type: support of distribution (int or real)
  stopifnot(is.list(pars))
  dist <- as_one_character(dist)
  log <- as_one_logical(log)
  type <- match.arg(type, c("int", "real"))
  args <- expand(dots = c(nlist(x, hu), pars))
  x <- args$x
  hu <- args$hu
  pars <- args[names(pars)]
  pdf <- paste0("d", dist)
  if (type == "int") {
    lccdf0 <- log(1 - do_call(pdf, c(0, pars)))
  } else {
    lccdf0 <- 0
  }
  out <- ifelse(x == 0, 
    log(hu),
    log(1 - hu) + do_call(pdf, c(list(x), pars, log = TRUE)) - lccdf0
  )
  if (!log) {
    out <- exp(out)
  }
  out
}

.phurdle <- function(q, dist, hu, pars, lower.tail, log.p, type) {
  # distribution function of hurdle models
  # Args:
  #   dist: name of the distribution
  #   hu: bernoulli hurdle parameter
  #   pars: list of parameters passed to pdf
  #   type: support of distribution (int or real)
  stopifnot(is.list(pars))
  dist <- as_one_character(dist)
  lower.tail <- as_one_logical(lower.tail)
  log.p <- as_one_logical(log.p)
  type <- match.arg(type, c("int", "real"))
  args <- expand(dots = c(nlist(q, hu), pars))
  q <- args$q
  hu <- args$hu
  pars <- args[names(pars)]
  cdf <- paste0("p", dist)
  out <- log(1 - hu) +
    do_call(cdf, c(list(q), pars, lower.tail = FALSE, log = TRUE))
  if (type == "int") {
    pdf <- paste0("d", dist)
    out <- out - log(1 - do_call(pdf, c(0, pars)))
  }
  if (lower.tail) {
    out <- 1 - exp(out)
    if (log.p) {
      out <- log(out)
    }
  } else {
    if (!log.p) {
      out <- exp(out)
    }
  }
  out
}

dcategorical <- function(x, eta, log = FALSE) {
  # density of the categorical distribution
  # with the softmax response function
  # Args:
  #   x: positive integers not greater than ncat
  #   eta: the linear predictor (of length or ncol ncat-1)
  #   log: return values on the log scale?
  if (is.null(dim(eta))) {
    eta <- matrix(eta, nrow = 1)
  }
  if (length(dim(eta)) != 2L) {
    stop2("eta must be a numeric vector or matrix.")
  }
  if (log) {
    out <- log_softmax(eta)
  } else {
    out <- softmax(eta)
  }
  out[, x]
}

pcategorical <- function(q, eta, log = FALSE) {
  # distribution function of the categorical distribution
  # with the softmax response function
  # Args:
  #   q: positive integers not greater than ncat
  #   eta: the linear predictor (of length or ncol ncat-1)  
  #   log: return values on the log scale?
  p <- dcategorical(seq_len(max(q)), eta = eta)
  out <- do_call(cbind, lapply(q, function(j) rowSums(as.matrix(p[, 1:j]))))
  if (log) {
    out <- log(out)
  }
  out
}

dmultinomial <- function(x, eta, log = FALSE) {
  # density of the multinomial distribution
  # with the softmax response function
  # Args:
  #   x: positive integers not greater than ncat
  #   eta: the linear predictor (of length or ncol ncat-1)
  #   log: return values on the log scale?
  if (is.null(dim(eta))) {
    eta <- matrix(eta, nrow = 1)
  }
  if (length(dim(eta)) != 2L) {
    stop2("eta must be a numeric vector or matrix.")
  }
  log_prob <- log_softmax(eta)
  size <- sum(x)
  x <- as_draws_matrix(x, dim = dim(eta))
  out <- lgamma(size + 1) + rowSums(x * log_prob - lgamma(x + 1))
  if (!log) {
    out <- exp(out)
  }
  out
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
  p <- do_call(cbind, rows)
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
  p <- do_call(cbind, rows)
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
  p <- do_call(cbind, rows)
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
        rowSums(eta[, 1:(k - 1)])
      }
      p[, 3:ncat] <- exp(sapply(3:ncat, .fun))
    }
  } else {
    mu <- ilink(eta, link)
    p <- cbind(apply(1 - mu[, 1:(ncat - 1)], 1, prod), 
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
  p <- do_call(paste0("d", family), args)
  .fun <- function(j) {
    rowSums(as.matrix(p[, 1:j]))
  }
  do_call(cbind, lapply(q, .fun))
}

# helper functions to shift arbitrary distributions
dshifted <- function(dist, x, shift = 0, ...) {
  do_call(paste0("d", dist), list(x - shift, ...))
}

pshifted <- function(dist, q, shift = 0, ...) {
  do_call(paste0("p", dist), list(q - shift, ...))
}

qshifted <- function(dist, p, shift = 0, ...) {
  do_call(paste0("q", dist), list(p, ...)) + shift
}

rshifted <- function(dist, n, shift = 0, ...) {
  do_call(paste0("r", dist), list(n, ...)) + shift
}
