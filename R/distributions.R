#density of student's distribution with parameters df, mu, and sigma
dstudent <- function(x, df, mu = 0, sigma = 1, log = FALSE) {
  log.const = log(gamma((df+1)/2)) - log(gamma(df/2)) - log(df*pi)/2 - log(sigma)
  out <- log.const - log(1 + 1/df * ((x - mu)/sigma)^2) * (df + 1)/2
  if (!log) out <- exp(out)
  out
}

#density of the multinormal distribution
dmultinormal <- function(x, mu, Sigma, log = TRUE) {
  k <- length(x)
  rooti <- backsolve(chol(Sigma),diag(k))
  quads <- colSums((crossprod(rooti, (x-mu)))^2)
  out <- -(k/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads
  if (!log) out <- exp(out)
  out
}