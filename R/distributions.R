#' The Student-t Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Student-t distribution with location \code{mu}, scale \code{sigma},
#' and degrees of freedom \code{df}.
#'
#' @name StudentT
#'
#' @param x Vector of quantiles.
#' @param q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of draws to sample from the distribution.
#' @param mu Vector of location values.
#' @param sigma Vector of scale values.
#' @param df Vector of degrees of freedom.
#' @param log Logical; If \code{TRUE}, values are returned on the log scale.
#' @param log.p Logical; If \code{TRUE}, values are returned on the log scale.
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
  if (isTRUE(any(sigma < 0))) {
    stop2("sigma must be non-negative.")
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
  if (isTRUE(any(sigma < 0))) {
    stop2("sigma must be non-negative.")
  }
  pt((q - mu) / sigma, df = df, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname StudentT
#' @export
qstudent_t <-  function(p, df, mu = 0, sigma = 1,
                        lower.tail = TRUE, log.p = FALSE) {
  if (isTRUE(any(sigma < 0))) {
    stop2("sigma must be non-negative.")
  }
  p <- validate_p_dist(p, lower.tail = lower.tail, log.p = log.p)
  mu + sigma * qt(p, df = df)
}

#' @rdname StudentT
#' @export
rstudent_t <- function(n, df, mu = 0, sigma = 1) {
  if (isTRUE(any(sigma < 0))) {
    stop2("sigma must be non-negative.")
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
#' @details See the Stan user's manual \url{https://mc-stan.org/documentation/}
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
  draws <- matrix(rnorm(n * p), nrow = n, ncol = p)
  mu + draws %*% chol(Sigma)
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
#' @details See the Stan user's manual \url{https://mc-stan.org/documentation/}
#'   for details on the parameterization
#'
#' @export
dmulti_student_t <- function(x, df, mu, Sigma, log = FALSE, check = FALSE) {
  if (is.vector(x) || length(dim(x)) == 1L) {
    x <- matrix(x, ncol = length(x))
  }
  p <- ncol(x)
  if (check) {
    if (isTRUE(any(df <= 0))) {
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
  if (isTRUE(any(df <= 0))) {
    stop2("df must be greater than 0.")
  }
  draws <- rmulti_normal(n, mu = rep(0, p), Sigma = Sigma, check = check)
  draws <- draws / sqrt(rchisq(n, df = df) / df)
  sweep(draws, 2, mu, "+")
}

#' The (Multivariate) Logistic Normal Distribution
#'
#' Density function and random generation for the (multivariate) logistic normal
#' distribution with latent mean vector \code{mu} and covariance matrix \code{Sigma}.
#'
#' @name LogisticNormal
#'
#' @inheritParams StudentT
#' @param x Vector or matrix of quantiles. If \code{x} is a matrix,
#'   each row is taken to be a quantile.
#' @param mu Mean vector with length equal to the number of dimensions.
#' @param Sigma Covariance matrix.
#' @param refcat A single integer indicating the reference category.
#'   Defaults to \code{1}.
#' @param check Logical; Indicates whether several input checks
#'   should be performed. Defaults to \code{FALSE} to improve
#'   efficiency.
#'
#' @export
dlogistic_normal <- function(x, mu, Sigma, refcat = 1, log = FALSE,
                             check = FALSE) {
  if (is.vector(x) || length(dim(x)) == 1L) {
    x <- matrix(x, ncol = length(x))
  }
  lx <- link_categorical(x, refcat)
  out <- dmulti_normal(lx, mu, Sigma, log = TRUE) - rowSums(log(x))
  if (!log) {
    out <- exp(out)
  }
  out
}

#' @rdname LogisticNormal
#' @export
rlogistic_normal <- function(n, mu, Sigma, refcat = 1, check = FALSE) {
  out <- rmulti_normal(n, mu, Sigma, check = check)
  inv_link_categorical(out, refcat = refcat)
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
  if (isTRUE(any(sigma < 0))) {
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
  if (isTRUE(any(sigma < 0))) {
    stop2("sigma must be non-negative.")
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
  if (isTRUE(any(sigma < 0))) {
    stop2("sigma must be non-negative.")
  }
  p <- validate_p_dist(p, lower.tail = lower.tail, log.p = log.p)
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
  if (isTRUE(any(sigma < 0))) {
    stop2("sigma must be non-negative.")
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

# convert skew-normal mixed-CP to DP parameterization
# @return a data.frame containing all relevant parameters
cp2dp <- function(mu = 0, sigma = 1, alpha = 0,
                  xi = NULL, omega = NULL, ...) {
  delta <- alpha / sqrt(1 + alpha^2)
  if (is.null(omega)) {
    omega <- sigma / sqrt(1 - 2 / pi * delta^2)
  }
  if (is.null(xi)) {
    xi <- mu - omega * delta * sqrt(2 / pi)
  }
  expand(dots = nlist(mu, sigma, alpha, xi, omega, delta, ...))
}

# helper function for qskew_normal
# code basis taken from sn::sn.cumulants
# uses xi and omega rather than mu and sigma
skew_normal_cumulants <- function(xi = 0, omega = 1, alpha = 0, n = 4) {
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

# CDF of the inverse gamma function
pinvgamma <- function(q, shape, rate, lower.tail = TRUE, log.p = FALSE) {
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
  if (isTRUE(any(kappa < 0))) {
    stop2("kappa must be non-negative")
  }
  # expects x in [-pi, pi] rather than [0, 2*pi] as CircStats::dvm
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
  if (isTRUE(any(kappa < 0))) {
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

  # code basis taken from CircStats::pvm but improved
  # considerably with respect to speed and stability
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
  if (isTRUE(any(kappa < 0))) {
    stop2("kappa must be non-negative")
  }
  args <- expand(mu = mu, kappa = kappa, length = n)
  mu <- args$mu
  kappa <- args$kappa
  rm(args)
  pi <- base::pi
  mu <- mu + pi

  # code basis taken from CircStats::rvm but improved
  # considerably with respect to speed and stability
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
  if (isTRUE(any(sigma < 0))) {
    stop2("sigma must be non-negative.")
  }
  if (isTRUE(any(beta < 0))) {
    stop2("beta must be non-negative.")
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
  if (isTRUE(any(sigma < 0))) {
    stop2("sigma must be non-negative.")
  }
  if (isTRUE(any(beta < 0))) {
    stop2("beta must be non-negative.")
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
  if (isTRUE(any(sigma < 0))) {
    stop2("sigma must be non-negative.")
  }
  if (isTRUE(any(beta < 0))) {
    stop2("beta must be non-negative.")
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
dfrechet <- function(x, loc = 0, scale = 1, shape = 1, log = FALSE) {
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
qfrechet <- function(p, loc = 0, scale = 1, shape = 1,
                     lower.tail = TRUE, log.p = FALSE) {
  if (isTRUE(any(scale <= 0))) {
    stop2("Argument 'scale' must be positive.")
  }
  if (isTRUE(any(shape <= 0))) {
    stop2("Argument 'shape' must be positive.")
  }
  p <- validate_p_dist(p, lower.tail = lower.tail, log.p = log.p)
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
      exp(2 * shape / mu) * pnorm(-sqrt(shape / q) * (q / mu + 1))
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

#' The Beta-binomial Distribution
#'
#' Cumulative density & mass functions, and random number generation for the
#' Beta-binomial distribution using the following re-parameterisation of the
#' \href{https://mc-stan.org/docs/2_29/functions-reference/beta-binomial-distribution.html}{Stan
#' Beta-binomial definition}:
#' \itemize{
#'  \item{\code{mu = alpha * beta}} mean probability of trial success.
#'  \item{\code{phi = (1 - mu) * beta}} precision or over-dispersion, component.
#' }
#'
#' @name BetaBinomial
#'
#' @inheritParams StudentT
#' @param x,q Vector of quantiles.
#' @param size Vector of number of trials (zero or more).
#' @param mu Vector of means.
#' @param phi Vector of precisions.
#'
#' @export
dbeta_binomial <- function(x, size, mu, phi, log = FALSE) {
  require_package("extraDistr")
  alpha <- mu * phi
  beta <- (1 - mu) * phi
  extraDistr::dbbinom(x, size, alpha = alpha, beta = beta, log = log)
}

#' @rdname BetaBinomial
#' @export
pbeta_binomial <- function(q, size, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  require_package("extraDistr")
  alpha <- mu * phi
  beta <- (1 - mu) * phi
  extraDistr::pbbinom(q, size, alpha = alpha, beta = beta,
                      lower.tail = lower.tail, log.p = log.p)
}

#' @rdname BetaBinomial
#' @export
rbeta_binomial <- function(n, size, mu, phi) {
  # beta location-scale probabilities
  probs <- rbeta(n, mu * phi, (1 - mu) * phi)
  # binomial draws
  rbinom(n, size = size, prob = probs)
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
dgen_extreme_value <- function(x, mu = 0, sigma = 1, xi = 0, log = FALSE) {
  if (isTRUE(any(sigma <= 0))) {
    stop2("sigma bust be positive.")
  }
  x <- (x - mu) / sigma
  args <- nlist(x, mu, sigma, xi)
  args <- do_call(expand, args)
  args$t <- with(args, 1 + xi * x)
  out <- with(args, ifelse(
    xi == 0,
    -log(sigma) - x - exp(-x),
    -log(sigma) - (1 + 1 / xi) * log(t) - t^(-1 / xi)
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
  if (isTRUE(any(sigma <= 0))) {
    stop2("sigma bust be positive.")
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
qgen_extreme_value <- function(p, mu = 0, sigma = 1, xi = 0,
                               lower.tail = TRUE, log.p = FALSE) {
  if (isTRUE(any(sigma <= 0))) {
    stop2("sigma bust be positive.")
  }
  p <- validate_p_dist(p, lower.tail = lower.tail, log.p = log.p)
  args <- nlist(p, mu, sigma, xi)
  args <- do_call(expand, args)
  out <- with(args, ifelse(
    xi == 0,
    mu - sigma * log(-log(p)),
    mu + (sigma * (1 - (-log(p))^xi)) / xi
  ))
  out
}

#' @rdname GenExtremeValue
#' @export
rgen_extreme_value <- function(n, mu = 0, sigma = 1, xi = 0) {
  if (isTRUE(any(sigma <= 0))) {
    stop2("sigma bust be positive.")
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
  p <- validate_p_dist(p, lower.tail = lower.tail, log.p = log.p)
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

# The Discrete Weibull Distribution
#
# Density, distribution function, quantile function and random generation
# for the discrete Weibull distribution with location \code{mu} and
# shape \code{shape}.
#
# @name DiscreteWeibull
#
# @inheritParams StudentT
# @param mu Location parameter in the unit interval.
# @param shape Positive shape parameter.
#
# @details See \code{vignette("brms_families")} for details
# on the parameterization.
#
# @export
ddiscrete_weibull <- function(x, mu, shape, log = FALSE) {
  if (isTRUE(any(mu < 0 | mu > 1))) {
    stop2("mu bust be between 0 and 1.")
  }
  if (isTRUE(any(shape <= 0))) {
    stop2("shape bust be positive.")
  }
  x <- round(x)
  out <- mu^x^shape - mu^(x + 1)^shape
  out[x < 0] <- 0
  if (log) {
    out <- log(out)
  }
  out
}

# @rdname DiscreteWeibull
# @export
pdiscrete_weibull <- function(x, mu, shape, lower.tail = TRUE, log.p = FALSE) {
  if (isTRUE(any(mu < 0 | mu > 1))) {
    stop2("mu bust be between 0 and 1.")
  }
  if (isTRUE(any(shape <= 0))) {
    stop2("shape bust be positive.")
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

# @rdname DiscreteWeibull
# @export
qdiscrete_weibull <- function(p, mu, shape, lower.tail = TRUE, log.p = FALSE) {
  if (isTRUE(any(mu < 0 | mu > 1))) {
    stop2("mu bust be between 0 and 1.")
  }
  if (isTRUE(any(shape <= 0))) {
    stop2("shape bust be positive.")
  }
  p <- validate_p_dist(p, lower.tail = lower.tail, log.p = log.p)
  ceiling((log(1 - p) / log(mu))^(1 / shape) - 1)
}

# @rdname DiscreteWeibull
# @export
rdiscrete_weibull <- function(n, mu, shape) {
  u <- runif(n, 0, 1)
  qdiscrete_weibull(u, mu, shape)
}

# mean of the discrete weibull distribution
# @param mu location parameter
# @param shape shape parameter
# @param M maximal evaluated element of the series
# @param thres threshold for new elements at which to stop evaluation
mean_discrete_weibull <- function(mu, shape, M = 1000, thres = 0.001) {
  opt_M <- ceiling(max((log(thres) / log(mu))^(1 / shape)))
  if (opt_M <= M) {
    M <- opt_M
  } else {
    # avoid the loop below running too slow
    warning2(
      "Approximating the mean of the 'discrete_weibull' ",
      "distribution failed and results be inaccurate."
    )
  }
  out <- 0
  for (y in seq_len(M)) {
    out <- out + mu^y^shape
  }
  # approximation of the residual series (see Englehart & Li, 2011)
  # returns unreasonably large values presumably due to numerical issues
  out
}

# PDF of the COM-Poisson distribution
# com_poisson in brms uses the mode parameterization
dcom_poisson <- function(x, mu, shape, log = FALSE) {
  x <- round(x)
  log_mu <- log(mu)
  log_Z <- log_Z_com_poisson(log_mu, shape)
  out <- shape * (x * log_mu - lgamma(x + 1)) - log_Z
  if (!log) {
    out <- exp(out)
  }
  out
}

# random numbers from the COM-Poisson distribution
rcom_poisson <- function(n, mu, shape, M = 10000) {
  n <- check_n_rdist(n, mu, shape)
  M <- as.integer(as_one_numeric(M))
  log_mu <- log(mu)
  # approximating log_Z may yield too large random draws
  log_Z <- log_Z_com_poisson(log_mu, shape, approx = FALSE)
  u <- runif(n, 0, 1)
  cdf <- exp(-log_Z)
  lfac <- 0
  y <- 0
  out <- rep(0, n)
  not_found <- cdf < u
  while (any(not_found) && y <= M) {
    y <- y + 1
    out[not_found] <- y
    lfac <- lfac + log(y)
    cdf <- cdf + exp(shape * (y * log_mu - lfac) - log_Z)
    not_found <- cdf < u
  }
  if (any(not_found)) {
    out[not_found] <- NA
    nfailed <- sum(not_found)
    warning2(
      "Drawing random numbers from the 'com_poisson' ",
      "distribution failed in ", nfailed, " cases."
    )
  }
  out
}

# CDF of the COM-Poisson distribution
pcom_poisson <- function(x, mu, shape, lower.tail = TRUE, log.p = FALSE) {
  x <- round(x)
  args <- expand(x = x, mu = mu, shape = shape)
  x <- args$x
  mu <- args$mu
  shape <- args$shape

  log_mu <- log(mu)
  log_Z <- log_Z_com_poisson(log_mu, shape)
  out <- rep(0, length(x))
  dim(out) <- attributes(args)$max_dim
  out[x > 0] <- log1p_exp(shape * log_mu)
  k <- 2
  lfac <- 0
  while (any(x >= k)) {
    lfac <- lfac + log(k)
    term <- shape * (k * log_mu - lfac)
    out[x >= k] <- log_sum_exp(out[x >= k], term)
    k <- k + 1
  }
  out <- out - log_Z
  out[out > 0] <- 0
  if (!lower.tail) {
    out <- log1m_exp(out)
  }
  if (!log.p) {
    out <- exp(out)
  }
  out
}

# log normalizing constant of the COM Poisson distribution
# @param log_mu log location parameter
# @param shape shape parameter
# @param M maximal evaluated element of the series
# @param thres threshold for new elements at which to stop evaluation
# @param approx use a closed form approximation of the mean if appropriate?
log_Z_com_poisson <- function(log_mu, shape, M = 10000, thres = 1e-16,
                              approx = TRUE) {
  if (isTRUE(any(shape <= 0))) {
    stop2("'shape' must be positive.")
  }
  if (isTRUE(any(shape == Inf))) {
    stop2("'shape' must be finite.")
  }
  approx <- as_one_logical(approx)
  args <- expand(log_mu = log_mu, shape = shape)
  log_mu <- args$log_mu
  shape <- args$shape

  out <- rep(NA, length(log_mu))
  dim(out) <- attributes(args)$max_dim
  use_poisson <- shape == 1
  if (any(use_poisson)) {
    # shape == 1 implies the poisson distribution
    out[use_poisson] <- exp(log_mu[use_poisson])
  }
  if (approx) {
    # use a closed form approximation if appropriate
    use_approx <- log_mu * shape >= log(1.5) & log_mu >= log(1.5)
    if (any(use_approx)) {
      out[use_approx] <- log_Z_com_poisson_approx(
        log_mu[use_approx], shape[use_approx]
      )
    }
  }
  use_exact <- is.na(out)
  if (any(use_exact)) {
    # direct computation of the truncated series
    M <- as.integer(as_one_numeric(M))
    thres <- as_one_numeric(thres)
    log_thres <- log(thres)
    log_mu <- log_mu[use_exact]
    shape <- shape[use_exact]
    # first 2 terms of the series
    out_exact <- log1p_exp(shape * log_mu)
    lfac <- 0
    k <- 2
    converged <- FALSE
    while (!converged && k <= M) {
      lfac <- lfac + log(k)
      term <- shape * (k * log_mu - lfac)
      out_exact <- log_sum_exp(out_exact, term)
      converged <- all(term <= log_thres)
      k <- k + 1
    }
    out[use_exact] <- out_exact
    if (!converged) {
      warning2(
        "Approximating the normalizing constant of the 'com_poisson' ",
        "distribution failed and results may be inaccurate."
      )
    }
  }
  out
}

# approximate the log normalizing constant of the COM Poisson distribution
# based on doi:10.1007/s10463-017-0629-6
log_Z_com_poisson_approx <- function(log_mu, shape) {
  shape_mu <- shape * exp(log_mu)
  shape2 <- shape^2
  # first 4 terms of the residual series
  log_sum_resid <- log(
    1 + shape_mu^(-1) * (shape2 - 1) / 24 +
      shape_mu^(-2) * (shape2 - 1) / 1152 * (shape2 + 23) +
      shape_mu^(-3) * (shape2 - 1) / 414720 *
        (5 * shape2^2 - 298 * shape2 + 11237)
  )
  shape_mu + log_sum_resid -
    ((log(2 * pi) + log_mu) * (shape - 1) / 2 + log(shape) / 2)
}

# compute the log mean of the COM Poisson distribution
# @param mu location parameter
# @param shape shape parameter
# @param M maximal evaluated element of the series
# @param thres threshold for new elements at which to stop evaluation
# @param approx use a closed form approximation of the mean if appropriate?
mean_com_poisson <- function(mu, shape, M = 10000, thres = 1e-16,
                             approx = TRUE) {
  if (isTRUE(any(shape <= 0))) {
    stop2("'shape' must be positive.")
  }
  if (isTRUE(any(shape == Inf))) {
    stop2("'shape' must be finite.")
  }
  approx <- as_one_logical(approx)
  args <- expand(mu = mu, shape = shape)
  mu <- args$mu
  shape <- args$shape

  out <- rep(NA, length(mu))
  dim(out) <- attributes(args)$max_dim
  use_poisson <- shape == 1
  if (any(use_poisson)) {
    # shape == 1 implies the poisson distribution
    out[use_poisson] <- mu[use_poisson]
  }
  if (approx) {
    # use a closed form approximation if appropriate
    use_approx <- mu^shape >= 1.5 & mu >= 1.5
    if (any(use_approx)) {
      out[use_approx] <- mean_com_poisson_approx(
        mu[use_approx], shape[use_approx]
      )
    }
  }
  use_exact <- is.na(out)
  if (any(use_exact)) {
    # direct computation of the truncated series
    M <- as.integer(as_one_numeric(M))
    thres <- as_one_numeric(thres)
    log_thres <- log(thres)
    mu <- mu[use_exact]
    shape <- shape[use_exact]
    log_mu <- log(mu)
    # first 2 terms of the series
    log_num <- shape * log_mu  # numerator
    log_Z <- log1p_exp(shape * log_mu)  # denominator
    lfac <- 0
    k <- 2
    converged <- FALSE
    while (!converged && k <= M) {
      log_k <- log(k)
      lfac <- lfac + log_k
      term <- shape * (k * log_mu - lfac)
      log_num <- log_sum_exp(log_num, log_k + term)
      log_Z <- log_sum_exp(log_Z, term)
      converged <- all(term <= log_thres)
      k <- k + 1
    }
    if (!converged) {
      warning2(
        "Approximating the mean of the 'com_poisson' ",
        "distribution failed and results may be inaccurate."
      )
    }
    out[use_exact] <- exp(log_num - log_Z)
  }
  out
}

# approximate the mean of COM-Poisson distribution
# based on doi:10.1007/s10463-017-0629-6
mean_com_poisson_approx <- function(mu, shape) {
  term <- 1 - (shape - 1) / (2 * shape) * mu^(-1) -
    (shape^2 - 1) / (24 * shape^2) * mu^(-2) -
    (shape^2 - 1) / (24 * shape^3) * mu^(-3)
  mu * term
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
  if (isTRUE(any(x < 0))) {
    stop2("x must be non-negative.")
  }
  if (!is_equal(rowSums(x), rep(1, nrow(x)))) {
    stop2("x must sum to 1 per row.")
  }
  if (isTRUE(any(alpha <= 0))) {
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
  if (isTRUE(any(alpha <= 0))) {
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
#' @param backend Name of the package to use as backend for the computations.
#'   Either \code{"Rwiener"} (the default) or \code{"rtdists"}.
#'   Can be set globally for the current \R session via the
#'   \code{"wiener_backend"} option (see \code{\link{options}}).
#'
#' @details
#' These are wrappers around functions of the \pkg{RWiener} or \pkg{rtdists}
#' package (depending on the chosen \code{backend}). See
#' \code{vignette("brms_families")} for details on the parameterization.
#'
#' @seealso \code{\link[RWiener:wienerdist]{wienerdist}},
#'   \code{\link[rtdists:Diffusion]{Diffusion}}
#'
#' @export
dwiener <- function(x, alpha, tau, beta, delta, resp = 1, log = FALSE,
                    backend = getOption("wiener_backend", "Rwiener")) {
  alpha <- as.numeric(alpha)
  tau <- as.numeric(tau)
  beta <- as.numeric(beta)
  delta <- as.numeric(delta)
  if (!is.character(resp)) {
    resp <- ifelse(resp, "upper", "lower")
  }
  log <- as_one_logical(log)
  backend <- match.arg(backend, c("Rwiener", "rtdists"))
  .dwiener <- paste0(".dwiener_", backend)
  args <- nlist(x, alpha, tau, beta, delta, resp)
  args <- as.list(do_call(expand, args))
  args$log <- log
  do_call(.dwiener, args)
}

# dwiener using Rwiener as backend
.dwiener_Rwiener <- function(x, alpha, tau, beta, delta, resp, log) {
  require_package("RWiener")
  .dwiener <- Vectorize(
    RWiener::dwiener,
    c("q", "alpha", "tau", "beta", "delta", "resp")
  )
  args <- nlist(q = x, alpha, tau, beta, delta, resp, give_log = log)
  do_call(.dwiener, args)
}

# dwiener using rtdists as backend
.dwiener_rtdists <- function(x, alpha, tau, beta, delta, resp, log) {
  require_package("rtdists")
  args <- list(
    rt = x, response = resp, a = alpha,
    t0 = tau, z = beta * alpha, v = delta
  )
  out <- do_call(rtdists::ddiffusion, args)
  if (log) {
    out <- log(out)
  }
  out
}

#' @rdname Wiener
#' @export
rwiener <- function(n, alpha, tau, beta, delta, types = c("q", "resp"),
                    backend = getOption("wiener_backend", "Rwiener")) {
  n <- as_one_numeric(n)
  alpha <- as.numeric(alpha)
  tau <- as.numeric(tau)
  beta <- as.numeric(beta)
  delta <- as.numeric(delta)
  types <- match.arg(types, several.ok = TRUE)
  backend <- match.arg(backend, c("Rwiener", "rtdists"))
  .rwiener <- paste0(".rwiener_", backend)
  args <- nlist(n, alpha, tau, beta, delta, types)
  do_call(.rwiener, args)
}

# rwiener using Rwiener as backend
.rwiener_Rwiener <- function(n, alpha, tau, beta, delta, types) {
  require_package("RWiener")
  max_len <- max(lengths(list(alpha, tau, beta, delta)))
  if (max_len > 1L) {
    if (!n %in% c(1, max_len)) {
      stop2("Can only sample exactly once for each condition.")
    }
    n <- 1
  }
  # helper function to return a numeric vector instead
  # of a data.frame with two columns as for RWiener::rwiener
  .rwiener_num <- function(n, alpha, tau, beta, delta, types) {
    out <- RWiener::rwiener(n, alpha, tau, beta, delta)
    out$resp <- ifelse(out$resp == "upper", 1, 0)
    if (length(types) == 1L) {
      out <- out[[types]]
    }
    out
  }
  # vectorized version of .rwiener_num
  .rwiener <- function(...) {
    fun <- Vectorize(
      .rwiener_num,
      c("alpha", "tau", "beta", "delta"),
      SIMPLIFY = FALSE
    )
    do_call(rbind, fun(...))
  }
  args <- nlist(n, alpha, tau, beta, delta, types)
  do_call(.rwiener, args)
}

# rwiener using rtdists as backend
.rwiener_rtdists <- function(n, alpha, tau, beta, delta, types) {
  require_package("rtdists")
  max_len <- max(lengths(list(alpha, tau, beta, delta)))
  if (max_len > 1L) {
    if (!n %in% c(1, max_len)) {
      stop2("Can only sample exactly once for each condition.")
    }
    n <- max_len
  }
  out <- rtdists::rdiffusion(
    n, a = alpha, t0 = tau, z = beta * alpha, v = delta
  )
  # TODO: use column names of rtdists in the output?
  names(out)[names(out) == "rt"] <- "q"
  names(out)[names(out) == "response"] <- "resp"
  out$resp <- ifelse(out$resp == "upper", 1, 0)
  if (length(types) == 1L) {
    out <- out[[types]]
  }
  out
}

# density of the cox proportional hazards model
# @param x currently ignored as the information is passed
#   via 'bhaz' and 'cbhaz'. Before exporting the cox distribution
#   functions, this needs to be refactored so that x is actually used
# @param mu positive location parameter
# @param bhaz baseline hazard
# @param cbhaz cumulative baseline hazard
dcox <- function(x, mu, bhaz, cbhaz, log = FALSE) {
  out <- hcox(x, mu, bhaz, cbhaz, log = TRUE) +
    pcox(x, mu, bhaz, cbhaz, lower.tail = FALSE, log.p = TRUE)
  if (!log) {
    out <- exp(out)
  }
  out
}

# hazard function of the cox model
hcox <- function(x, mu, bhaz, cbhaz, log = FALSE) {
  out <- log(bhaz) + log(mu)
  if (!log) {
    out <- exp(out)
  }
  out
}

# distribution function of the cox model
pcox <- function(q, mu, bhaz, cbhaz, lower.tail = TRUE, log.p = FALSE) {
  log_surv <- -cbhaz * mu
  if (lower.tail) {
    if (log.p) {
      out <- log1m_exp(log_surv)
    } else {
      out <- 1 - exp(log_surv)
    }
  } else {
    if (log.p) {
      out <- log_surv
    } else {
      out <- exp(log_surv)
    }
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
#' @param zi zero-inflation probability
#' @param mu,lambda location parameter
#' @param shape,shape1,shape2 shape parameter
#' @param phi precision parameter
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
dzero_inflated_beta_binomial <- function(x, size, mu, phi, zi, log = FALSE) {
  pars <- nlist(size, mu, phi)
  .dzero_inflated(x, "beta_binomial", zi, pars, log)
}

#' @rdname ZeroInflated
#' @export
pzero_inflated_beta_binomial <- function(q, size, mu, phi, zi,
                                         lower.tail = TRUE, log.p = FALSE) {
  pars <- nlist(size, mu, phi)
  .pzero_inflated(q, "beta_binomial", zi, pars, lower.tail, log.p)
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

# @rdname ZeroInflated
# @export
dzero_inflated_asym_laplace <- function(x, mu, sigma, quantile, zi,
                                        log = FALSE) {
  pars <- nlist(mu, sigma, quantile)
  # zi_asym_laplace is technically a hurdle model
  .dhurdle(x, "asym_laplace", zi, pars, log, type = "real")
}

# @rdname ZeroInflated
# @export
pzero_inflated_asym_laplace <- function(q, mu, sigma, quantile, zi,
                                        lower.tail = TRUE, log.p = FALSE) {
  pars <- nlist(mu, sigma, quantile)
  # zi_asym_laplace is technically a hurdle model
  .phurdle(q, "asym_laplace", zi, pars, lower.tail, log.p,
           type = "real", lb = -Inf, ub = Inf)
}

# density of a zero-inflated distribution
# @param dist name of the distribution
# @param zi bernoulli zero-inflated parameter
# @param pars list of parameters passed to pdf
.dzero_inflated <- function(x, dist, zi, pars, log) {
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

# CDF of a zero-inflated distribution
# @param dist name of the distribution
# @param zi bernoulli zero-inflated parameter
# @param pars list of parameters passed to pdf
# @param lb lower bound of the conditional distribution
# @param ub upper bound of the conditional distribution
.pzero_inflated <- function(q, dist, zi, pars, lower.tail, log.p,
                            lb = 0, ub = Inf) {
  stopifnot(is.list(pars))
  dist <- as_one_character(dist)
  lower.tail <- as_one_logical(lower.tail)
  log.p <- as_one_logical(log.p)
  lb <- as_one_numeric(lb)
  ub <- as_one_numeric(ub)
  args <- expand(dots = c(nlist(q, zi), pars))
  q <- args$q
  zi <- args$zi
  pars <- args[names(pars)]
  cdf <- paste0("p", dist)
  # compute log CCDF values
  out <- log(1 - zi) +
    do_call(cdf, c(list(q), pars, lower.tail = FALSE, log.p = TRUE))
  # take the limits of the distribution into account
  out <- ifelse(q < lb, 0, out)
  out <- ifelse(q > ub, -Inf, out)
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
#' @param hu hurdle probability
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

# density of a hurdle distribution
# @param dist name of the distribution
# @param hu bernoulli hurdle parameter
# @param pars list of parameters passed to pdf
# @param type support of distribution (int or real)
.dhurdle <- function(x, dist, hu, pars, log, type) {
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

# CDF of a hurdle distribution
# @param dist name of the distribution
# @param hu bernoulli hurdle parameter
# @param pars list of parameters passed to pdf
# @param type support of distribution (int or real)
# @param lb lower bound of the conditional distribution
# @param ub upper bound of the conditional distribution
.phurdle <- function(q, dist, hu, pars, lower.tail, log.p, type,
                     lb = 0, ub = Inf) {
  stopifnot(is.list(pars))
  dist <- as_one_character(dist)
  lower.tail <- as_one_logical(lower.tail)
  log.p <- as_one_logical(log.p)
  type <- match.arg(type, c("int", "real"))
  lb <- as_one_numeric(lb)
  ub <- as_one_numeric(ub)
  args <- expand(dots = c(nlist(q, hu), pars))
  q <- args$q
  hu <- args$hu
  pars <- args[names(pars)]
  cdf <- paste0("p", dist)
  # compute log CCDF values
  out <- log(1 - hu) +
    do_call(cdf, c(list(q), pars, lower.tail = FALSE, log.p = TRUE))
  if (type == "int") {
    pdf <- paste0("d", dist)
    out <- out - log(1 - do_call(pdf, c(0, pars)))
  }
  out <- ifelse(q < 0, log_sum_exp(out, log(hu)), out)
  # take the limits of the distribution into account
  out <- ifelse(q < lb, 0, out)
  out <- ifelse(q > ub, -Inf, out)
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

#' @rdname Mixcure
#' @export
dmixcure_lognormal <- function(x, mu, sigma, inc, log = FALSE) {
  pars <- list(meanlog = mu, sdlog = sigma)
  .dmixcure(x, "lnorm", inc, pars, log)
}

#' @rdname Mixcure
#' @export
pmixcure_lognormal <- function(q, mu, sigma, inc, lower.tail = TRUE, log.p = FALSE) {
  pars <- list(meanlog = mu, sdlog = sigma)
  .pmixcure(q, "lnorm", inc, pars, lower.tail, log.p)
}

#' @rdname Mixcure
#' @export
dmixcure_weibull <- function(x, shape, scale, inc, log = FALSE) {
  pars <- list(shape = shape, scale = scale)
  .dmixcure(x, "weibull", inc, pars, log)
}

#' @rdname Mixcure
#' @export
pmixcure_weibull <- function(q, shape, scale, inc, lower.tail = TRUE, log.p = FALSE) {
  pars <- list(shape = shape, scale = scale)
  .pmixcure(q, "weibull", inc, pars, lower.tail, log.p)
}

# density of a mixcure distribution
# @param dist name of the distribution
# @param inc bernoulli incidence parameter
# @param pars list of parameters passed to pdf
.dmixcure <- function(x, dist, inc, pars, log) {
  stopifnot(is.list(pars))
  dist <- as_one_character(dist)
  log <- as_one_logical(log)
  args <- expand(dots = c(nlist(x, inc), pars))
  x <- args$x
  inc <- args$inc
  pars <- args[names(pars)]
  pdf <- paste0("d", dist)
  # incidence part (not censored): pi(z) * f(t | x)
  out <- log(inc) + do_call(pdf, c(list(x), pars, log = TRUE))
  if (!log) {
    out <- exp(out)
  }
  out
}

# CDF of a mixcure distribution
# @param dist name of the distribution
# @param inc bernoulli incidence parameter
# @param pars list of parameters passed to pdf
# @param lb lower bound of the conditional distribution
# @param ub upper bound of the conditional distribution
.pmixcure <- function(q, dist, inc, pars, lower.tail, log.p, lb = 0, ub = Inf) {
  stopifnot(is.list(pars))
  dist <- as_one_character(dist)
  lower.tail <- as_one_logical(lower.tail)
  log.p <- as_one_logical(log.p)
  args <- expand(dots = c(nlist(q, inc), pars))
  q <- args$q
  inc <- args$inc
  pars <- args[names(pars)]
  cdf <- paste0("p", dist)
  # compute log CCDF values
  # latency part (right-censored): [1 - pi(z)] + pi(z) * S(t | x)
  out <- matrixStats::logSumExp(c(
    1, -inc,
    inc + do_call(cdf, c(list(q), pars, lower.tail = FALSE, log.p = TRUE))
  ))
  # # take the limits of the distribution into account
  out <- ifelse(q < lb, 0, out)
  out <- ifelse(q > ub, -Inf, out)
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

# density of the categorical distribution with the softmax transform
# @param x positive integers not greater than ncat
# @param eta the linear predictor (of length or ncol ncat)
# @param log return values on the log scale?
dcategorical <- function(x, eta, log = FALSE) {
  if (is.null(dim(eta))) {
    eta <- matrix(eta, nrow = 1)
  }
  if (length(dim(eta)) != 2L) {
    stop2("eta must be a numeric vector or matrix.")
  }
  out <- inv_link_categorical(eta, log = log, refcat = NULL)
  out[, x, drop = FALSE]
}

# generic inverse link function for the categorical family
#
# @param x Matrix (S x `ncat` or S x `ncat - 1` (depending on `refcat_obj`),
#   with S denoting the number of posterior draws and `ncat` denoting the number
#   of response categories) with values of `eta` for one observation (see
#   dcategorical()) or an array (S x N x `ncat` or S x N x `ncat - 1` (depending
#   on `refcat_obj`)) containing the same values as the matrix just described,
#   but for N observations.
# @param refcat Integer indicating the reference category to be inserted in 'x'.
#   If NULL, `x` is not modified at all.
# @param log Logical (length 1) indicating whether to log the return value.
#
# @return If `x` is a matrix, then a matrix (S x `ncat`, with S denoting the
#   number of posterior draws and `ncat` denoting the number of response
#   categories) containing the values of the inverse-link function applied to
#   `x`. If `x` is an array, then an array (S x N x `ncat`) containing the same
#   values as the matrix just described, but for N observations.
inv_link_categorical <- function(x, refcat = 1, log = FALSE) {
  if (!is.null(refcat)) {
    x <- insert_refcat(x, refcat = refcat)
  }
  out <- log_softmax(x)
  if (!log) {
    out <- exp(out)
  }
  out
}

# generic link function for the categorical family
#
# @param x Matrix (S x `ncat`, with S denoting the number of posterior draws and
#   `ncat` denoting the number of response categories) of probabilities for the
#   response categories or an array (S x N x `ncat`) containing the same values
#   as the matrix just described, but for N observations.
# @param refcat Numeric (length 1) giving the index of the reference category.
# @param return_refcat Logical (length 1) indicating whether to include the
#   reference category in the return value.
#
# @return If `x` is a matrix, then a matrix (S x `ncat` or S x `ncat - 1`
#   (depending on `return_refcat`), with S denoting the number of posterior
#   draws and `ncat` denoting the number of response categories) containing the
#   values of the link function applied to `x`. If `x` is an array, then an
#   array (S x N x `ncat` or S x N x `ncat - 1` (depending on `return_refcat`))
#   containing the same values as the matrix just described, but for N
#   observations.
link_categorical <- function(x, refcat = 1, return_refcat = FALSE) {
  ndim <- length(dim(x))
  marg_noncat <- seq_along(dim(x))[-ndim]
  if (return_refcat) {
    x_tosweep <- x
  } else {
    x_tosweep <- slice(x, ndim, -refcat, drop = FALSE)
  }
  log(sweep(
    x_tosweep,
    MARGIN = marg_noncat,
    STATS = slice(x, ndim, refcat),
    FUN = "/"
  ))
}

# CDF of the categorical distribution with the softmax transform
# @param q positive integers not greater than ncat
# @param eta the linear predictor (of length or ncol ncat)
# @param log.p return values on the log scale?
pcategorical <- function(q, eta, log.p = FALSE) {
  p <- dcategorical(seq_len(max(q)), eta = eta)
  out <- cblapply(q, function(j) rowSums(p[, 1:j, drop = FALSE]))
  if (log.p) {
    out <- log(out)
  }
  out
}

# density of the multinomial distribution with the softmax transform
# @param x positive integers not greater than ncat
# @param eta the linear predictor (of length or ncol ncat)
# @param log return values on the log scale?
dmultinomial <- function(x, eta, log = FALSE) {
  if (is.null(dim(eta))) {
    eta <- matrix(eta, nrow = 1)
  }
  if (length(dim(eta)) != 2L) {
    stop2("eta must be a numeric vector or matrix.")
  }
  log_prob <- log_softmax(eta)
  size <- sum(x)
  x <- data2draws(x, dim = dim(eta))
  out <- lgamma(size + 1) + rowSums(x * log_prob - lgamma(x + 1))
  if (!log) {
    out <- exp(out)
  }
  out
}

# density of the cumulative distribution
#
# @param x Integer vector containing response category indices to return the
#   "densities" (probability masses) for.
# @param eta Vector (length S, with S denoting the number of posterior draws) of
#   linear predictor draws.
# @param thres Matrix (S x `ncat - 1`, with S denoting the number of posterior
#   draws and `ncat` denoting the number of response categories) of threshold
#   draws.
# @param disc Vector (length S, with S denoting the number of posterior draws,
#   or length 1 for recycling) of discrimination parameter draws.
# @param link Character vector (length 1) giving the name of the link function.
#
# @return A matrix (S x `length(x)`) containing the values of the inverse-link
#   function applied to `disc * (thres - eta)`.
dcumulative <- function(x, eta, thres, disc = 1, link = "logit") {
  eta <- disc * (thres - eta)
  if (link == "identity") {
    out <- eta
  } else {
    out <- inv_link_cumulative(eta, link = link)
  }
  out[, x, drop = FALSE]
}

# generic inverse link function for the cumulative family
#
# @param x Matrix (S x `ncat - 1`, with S denoting the number of posterior draws
#   and `ncat` denoting the number of response categories) with values of
#   `disc * (thres - eta)` for one observation (see dcumulative()) or an array
#   (S x N x `ncat - 1`) containing the same values as the matrix just
#   described, but for N observations.
# @param link Character vector (length 1) giving the name of the link function.
#
# @return If `x` is a matrix, then a matrix (S x `ncat`, with S denoting the
#   number of posterior draws and `ncat` denoting the number of response
#   categories) containing the values of the inverse-link function applied to
#   `x`. If `x` is an array, then an array (S x N x `ncat`) containing the same
#   values as the matrix just described, but for N observations.
inv_link_cumulative <- function(x, link) {
  x <- inv_link(x, link)
  ndim <- length(dim(x))
  dim_noncat <- dim(x)[-ndim]
  ones_arr <- array(1, dim = c(dim_noncat, 1))
  zeros_arr <- array(0, dim = c(dim_noncat, 1))
  abind::abind(x, ones_arr) - abind::abind(zeros_arr, x)
}

# generic link function for the cumulative family
#
# @param x Matrix (S x `ncat`, with S denoting the number of posterior draws and
#   `ncat` denoting the number of response categories) of probabilities for the
#   response categories or an array (S x N x `ncat`) containing the same values
#   as the matrix just described, but for N observations.
# @param link Character string (length 1) giving the name of the link function.
#
# @return If `x` is a matrix, then a matrix (S x `ncat - 1`, with S denoting the
#   number of posterior draws and `ncat` denoting the number of response
#   categories) containing the values of the link function applied to `x`. If
#   `x` is an array, then an array (S x N x `ncat - 1`) containing the same
#   values as the matrix just described, but for N observations.
link_cumulative <- function(x, link) {
  ndim <- length(dim(x))
  ncat <- dim(x)[ndim]
  dim_noncat <- dim(x)[-ndim]
  nthres <- dim(x)[ndim] - 1
  marg_noncat <- seq_along(dim(x))[-ndim]
  dim_t <- c(nthres, dim_noncat)
  x <- apply(slice(x, ndim, -ncat, drop = FALSE), marg_noncat, cumsum)
  x <- aperm(array(x, dim = dim_t), perm = c(marg_noncat + 1, 1))
  link(x, link)
}

# density of the sratio distribution
#
# @param x Integer vector containing response category indices to return the
#   "densities" (probability masses) for.
# @param eta Vector (length S, with S denoting the number of posterior draws) of
#   linear predictor draws.
# @param thres Matrix (S x `ncat - 1`, with S denoting the number of posterior
#   draws and `ncat` denoting the number of response categories) of threshold
#   draws.
# @param disc Vector (length S, with S denoting the number of posterior draws,
#   or length 1 for recycling) of discrimination parameter draws.
# @param link Character vector (length 1) giving the name of the link function.
#
# @return A matrix (S x `length(x)`) containing the values of the inverse-link
#   function applied to `disc * (thres - eta)`.
dsratio <- function(x, eta, thres, disc = 1, link = "logit") {
  eta <- disc * (thres - eta)
  if (link == "identity") {
    out <- eta
  } else {
    out <- inv_link_sratio(eta, link = link)
  }
  out[, x, drop = FALSE]
}

# generic inverse link function for the sratio family
#
# @param x Matrix (S x `ncat - 1`, with S denoting the number of posterior draws
#   and `ncat` denoting the number of response categories) with values of
#   `disc * (thres - eta)` for one observation (see dsratio()) or an array
#   (S x N x `ncat - 1`) containing the same values as the matrix just
#   described, but for N observations.
# @param link Character vector (length 1) giving the name of the link function.
#
# @return If `x` is a matrix, then a matrix (S x `ncat`, with S denoting the
#   number of posterior draws and `ncat` denoting the number of response
#   categories) containing the values of the inverse-link function applied to
#   `x`. If `x` is an array, then an array (S x N x `ncat`) containing the same
#   values as the matrix just described, but for N observations.
inv_link_sratio <- function(x, link) {
  x <- inv_link(x, link)
  ndim <- length(dim(x))
  dim_noncat <- dim(x)[-ndim]
  nthres <- dim(x)[ndim]
  marg_noncat <- seq_along(dim(x))[-ndim]
  ones_arr <- array(1, dim = c(dim_noncat, 1))
  dim_t <- c(nthres, dim_noncat)
  Sx_cumprod <- aperm(
    array(apply(1 - x, marg_noncat, cumprod), dim = dim_t),
    perm = c(marg_noncat + 1, 1)
  )
  abind::abind(x, ones_arr) * abind::abind(ones_arr, Sx_cumprod)
}

# generic link function for the sratio family
#
# @param x Matrix (S x `ncat`, with S denoting the number of posterior draws and
#   `ncat` denoting the number of response categories) of probabilities for the
#   response categories or an array (S x N x `ncat`) containing the same values
#   as the matrix just described, but for N observations.
# @param link Character string (length 1) giving the name of the link function.
#
# @return If `x` is a matrix, then a matrix (S x `ncat - 1`, with S denoting the
#   number of posterior draws and `ncat` denoting the number of response
#   categories) containing the values of the link function applied to `x`. If
#   `x` is an array, then an array (S x N x `ncat - 1`) containing the same
#   values as the matrix just described, but for N observations.
link_sratio <- function(x, link) {
  ndim <- length(dim(x))
  .F_k <- function(k) {
    if (k == 1) {
      prev_res <- list(F_k = NULL, S_km1_prod = 1)
    } else {
      prev_res <- .F_k(k - 1)
    }
    F_k <- slice(x, ndim, k, drop = FALSE) / prev_res$S_km1_prod
    .out <- list(
      F_k = abind::abind(prev_res$F_k, F_k),
      S_km1_prod = prev_res$S_km1_prod * (1 - F_k)
    )
    return(.out)
  }
  x <- .F_k(dim(x)[ndim] - 1)$F_k
  link(x, link)
}

# density of the cratio distribution
#
# @param x Integer vector containing response category indices to return the
#   "densities" (probability masses) for.
# @param eta Vector (length S, with S denoting the number of posterior draws) of
#   linear predictor draws.
# @param thres Matrix (S x `ncat - 1`, with S denoting the number of posterior
#   draws and `ncat` denoting the number of response categories) of threshold
#   draws.
# @param disc Vector (length S, with S denoting the number of posterior draws,
#   or length 1 for recycling) of discrimination parameter draws.
# @param link Character vector (length 1) giving the name of the link function.
#
# @return A matrix (S x `length(x)`) containing the values of the inverse-link
#   function applied to `disc * (thres - eta)`.
dcratio <- function(x, eta, thres, disc = 1, link = "logit") {
  eta <- disc * (eta - thres)
  if (link == "identity") {
    out <- eta
  } else {
    out <- inv_link_cratio(eta, link = link)
  }
  out[, x, drop = FALSE]
}

# generic inverse link function for the cratio family
#
# @param x Matrix (S x `ncat - 1`, with S denoting the number of posterior draws
#   and `ncat` denoting the number of response categories) with values of
#   `disc * (thres - eta)` for one observation (see dcratio()) or an array
#   (S x N x `ncat - 1`) containing the same values as the matrix just
#   described, but for N observations.
# @param link Character vector (length 1) giving the name of the link function.
#
# @return If `x` is a matrix, then a matrix (S x `ncat`, with S denoting the
#   number of posterior draws and `ncat` denoting the number of response
#   categories) containing the values of the inverse-link function applied to
#   `x`. If `x` is an array, then an array (S x N x `ncat`) containing the same
#   values as the matrix just described, but for N observations.
inv_link_cratio <- function(x, link) {
  x <- inv_link(x, link)
  ndim <- length(dim(x))
  dim_noncat <- dim(x)[-ndim]
  nthres <- dim(x)[ndim]
  marg_noncat <- seq_along(dim(x))[-ndim]
  ones_arr <- array(1, dim = c(dim_noncat, 1))
  dim_t <- c(nthres, dim_noncat)
  x_cumprod <- aperm(
    array(apply(x, marg_noncat, cumprod), dim = dim_t),
    perm = c(marg_noncat + 1, 1)
  )
  abind::abind(1 - x, ones_arr) * abind::abind(ones_arr, x_cumprod)
}

# generic link function for the cratio family
#
# @param x Matrix (S x `ncat`, with S denoting the number of posterior draws and
#   `ncat` denoting the number of response categories) of probabilities for the
#   response categories or an array (S x N x `ncat`) containing the same values
#   as the matrix just described, but for N observations.
# @param link Character string (length 1) giving the name of the link function.
#
# @return If `x` is a matrix, then a matrix (S x `ncat - 1`, with S denoting the
#   number of posterior draws and `ncat` denoting the number of response
#   categories) containing the values of the link function applied to `x`. If
#   `x` is an array, then an array (S x N x `ncat - 1`) containing the same
#   values as the matrix just described, but for N observations.
link_cratio <- function(x, link) {
  ndim <- length(dim(x))
  .F_k <- function(k) {
    if (k == 1) {
      prev_res <- list(F_k = NULL, F_km1_prod = 1)
    } else {
      prev_res <- .F_k(k - 1)
    }
    F_k <- 1 - slice(x, ndim, k, drop = FALSE) / prev_res$F_km1_prod
    .out <- list(
      F_k = abind::abind(prev_res$F_k, F_k),
      F_km1_prod = prev_res$F_km1_prod * F_k
    )
    return(.out)
  }
  x <- .F_k(dim(x)[ndim] - 1)$F_k
  link(x, link)
}

# density of the acat distribution
#
# @param x Integer vector containing response category indices to return the
#   "densities" (probability masses) for.
# @param eta Vector (length S, with S denoting the number of posterior draws) of
#   linear predictor draws.
# @param thres Matrix (S x `ncat - 1`, with S denoting the number of posterior
#   draws and `ncat` denoting the number of response categories) of threshold
#   draws.
# @param disc Vector (length S, with S denoting the number of posterior draws,
#   or length 1 for recycling) of discrimination parameter draws.
# @param link Character vector (length 1) giving the name of the link function.
#
# @return A matrix (S x `length(x)`) containing the values of the inverse-link
#   function applied to `disc * (thres - eta)`.
dacat <- function(x, eta, thres, disc = 1, link = "logit") {
  eta <- disc * (eta - thres)
  if (link == "identity") {
    out <- eta
  } else {
    out <- inv_link_acat(eta, link = link)
  }
  out[, x, drop = FALSE]
}

# generic inverse link function for the acat family
#
# @param x Matrix (S x `ncat - 1`, with S denoting the number of posterior draws
#   and `ncat` denoting the number of response categories) with values of
#   `disc * (thres - eta)` (see dacat()).
# @param link Character vector (length 1) giving the name of the link function.
#
# @return A matrix (S x `ncat`, with S denoting the number of posterior draws
#   and `ncat` denoting the number of response categories) containing the values
#   of the inverse-link function applied to `x`.
inv_link_acat <- function(x, link) {
  ndim <- length(dim(x))
  dim_noncat <- dim(x)[-ndim]
  nthres <- dim(x)[ndim]
  marg_noncat <- seq_along(dim(x))[-ndim]
  ones_arr <- array(1, dim = c(dim_noncat, 1))
  dim_t <- c(nthres, dim_noncat)
  if (link == "logit") {
    # faster evaluation in this case
    exp_x_cumprod <- aperm(
      array(apply(exp(x), marg_noncat, cumprod), dim = dim_t),
      perm = c(marg_noncat + 1, 1)
    )
    out <- abind::abind(ones_arr, exp_x_cumprod)
  } else {
    x <- inv_link(x, link)
    x_cumprod <- aperm(
      array(apply(x, marg_noncat, cumprod), dim = dim_t),
      perm = c(marg_noncat + 1, 1)
    )
    Sx_cumprod_rev <- apply(
      1 - slice(x, ndim, rev(seq_len(nthres)), drop = FALSE),
      marg_noncat, cumprod
    )
    Sx_cumprod_rev <- aperm(
      array(Sx_cumprod_rev, dim = dim_t),
      perm = c(marg_noncat + 1, 1)
    )
    Sx_cumprod_rev <- slice(
      Sx_cumprod_rev, ndim, rev(seq_len(nthres)), drop = FALSE
    )
    out <- abind::abind(ones_arr, x_cumprod) *
      abind::abind(Sx_cumprod_rev, ones_arr)
  }
  catsum <- array(apply(out, marg_noncat, sum), dim = dim_noncat)
  sweep(out, marg_noncat, catsum, "/")
}

# generic link function for the acat family
#
# @param x Matrix (S x `ncat`, with S denoting the number of posterior draws and
#   `ncat` denoting the number of response categories) of probabilities for the
#   response categories or an array (S x N x `ncat`) containing the same values
#   as the matrix just described, but for N observations.
# @param link Character string (length 1) giving the name of the link function.
#
# @return If `x` is a matrix, then a matrix (S x `ncat - 1`, with S denoting the
#   number of posterior draws and `ncat` denoting the number of response
#   categories) containing the values of the link function applied to `x`. If
#   `x` is an array, then an array (S x N x `ncat - 1`) containing the same
#   values as the matrix just described, but for N observations.
link_acat <- function(x, link) {
  ndim <- length(dim(x))
  ncat <- dim(x)[ndim]
  x <- slice(x, ndim, -1, drop = FALSE) / slice(x, ndim, -ncat, drop = FALSE)
  if (link == "logit") {
    # faster evaluation in this case
    out <- log(x)
  } else {
    x <- inv_odds(x)
    out <- link(x, link)
  }
  out
}

# CDF for ordinal distributions
# @param q positive integers not greater than ncat
# @param eta draws of the linear predictor
# @param thres draws of threshold parameters
# @param disc draws of the discrimination parameter
# @param family a character string naming the family
# @param link a character string naming the link
# @return a matrix of probabilities P(x <= q)
pordinal <- function(q, eta, thres, disc = 1, family = NULL, link = "logit") {
  family <- as_one_character(family)
  link <- as_one_character(link)
  args <- nlist(x = seq_len(max(q)), eta, thres, disc, link)
  p <- do_call(paste0("d", family), args)
  .fun <- function(j) rowSums(as.matrix(p[, 1:j, drop = FALSE]))
  cblapply(q, .fun)
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

# validate argument p in q<dist> functions
validate_p_dist <- function(p, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  if (isTRUE(any(p <= 0)) || isTRUE(any(p >= 1))) {
    stop2("'p' must contain probabilities in (0,1)")
  }
  p
}

# check if 'n' in r<dist> functions is valid
# @param n number of desired random draws
# @param .. parameter vectors
# @return validated 'n'
check_n_rdist <- function(n, ...) {
  n <- as.integer(as_one_numeric(n))
  max_len <- max(lengths(list(...)))
  if (max_len > 1L) {
    if (!n %in% c(1, max_len)) {
      stop2("'n' must match the maximum length of the parameter vectors.")
    }
    n <- max_len
  }
  n
}
