#density of student's distribution with parameters df, mu, and sigma
dstudent <- function(x, df = stop("df is required"), mu = 0, sigma = 1, log = FALSE) {
  if (log) dt((x - mu)/sigma, df = df, log = TRUE) - log(sigma)
  else dt((x - mu)/sigma, df = df)/sigma
}

#distribution function of student's distribution with parameters df, mu, and sigma
pstudent <- function(q, df = stop("df is required"), mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  pt((q - mu)/sigma, df = df, lower.tail = lower.tail, log.p = log.p)
}

#rquantiles of student's distribution with parameters df, mu, and sigma
qstudent <-  function(p, df = stop("df is required"), mu = 0, sigma = 1) {
  mu + sigma * qt(p, df = df)
}

#random values of student's distribution with parameters df, mu, and sigma
rstudent <-  function(n, df = stop("df is required"), mu = 0, sigma = 1) {
  mu + sigma * rt(n, df = df)
}

#density of the multinormal distribution with parameters mu and Sigma
dmultinormal <- function(x, mu, Sigma, log = TRUE) {
  k <- length(x)
  rooti <- backsolve(chol(Sigma),diag(k))
  quads <- colSums((crossprod(rooti, (x-mu)))^2)
  out <- -(k/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads
  if (!log) out <- exp(out)
  out
}