# density of student's distribution with parameters df, mu, and sigma
dstudent <- function(x, df = stop("df is required"), mu = 0, sigma = 1, log = FALSE) {
  if (log) dt((x - mu)/sigma, df = df, log = TRUE) - log(sigma)
  else dt((x - mu)/sigma, df = df)/sigma
}

# distribution function of student's distribution with parameters df, mu, and sigma
pstudent <- function(q, df = stop("df is required"), mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  pt((q - mu)/sigma, df = df, lower.tail = lower.tail, log.p = log.p)
}

# rquantiles of student's distribution with parameters df, mu, and sigma
qstudent <-  function(p, df = stop("df is required"), mu = 0, sigma = 1) {
  mu + sigma * qt(p, df = df)
}

# random values of student's distribution with parameters df, mu, and sigma
rstudent <-  function(n, df = stop("df is required"), mu = 0, sigma = 1) {
  mu + sigma * rt(n, df = df)
}

# density of the multinormal distribution with parameters mu and Sigma (not vectorized)
dmultinormal <- function(x, mu, Sigma, log = TRUE) {
  k <- length(x)
  rooti <- backsolve(chol(Sigma),diag(k))
  quads <- colSums((crossprod(rooti, (x-mu)))^2)
  out <- -(k/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads
  if (!log) out <- exp(out)
  out
}

# random values of the multinormal distribution with parameters mu and Sigma
rmultinormal <- function(n, mu, Sigma, check = FALSE) {
  p <- length(mu)
  if (check) {
    if (!all(dim(Sigma) == c(p, p))) 
      stop("incompatible arguments")
    if (!isSymmetric(unname(Sigma)))
      stop("Sigma is not symmetric")
  }
  cholSigma <- chol(Sigma)
  samples <- matrix(rnorm(n*p), nrow = n, ncol = p)
  mu + samples %*% cholSigma
}

# density of the categorical distribution
# 
# @param same arguments as dcumulative
dcategorical <- function(x, eta, cat, link = "logit") {
  if (is.null(dim(eta))) eta <- matrix(eta, nrow = 1)
  if (length(dim(eta)) != 2 || !is.numeric(eta)) 
    stop("eta must be a numeric vector or matrix")
  if (missing(cat)) cat <- ncol(eta)
  if (link == "logit")
    p <- exp(cbind(rep(0, nrow(eta)), eta[,1:(cat-1)]))
  else stop(paste("Link", link, "not supported"))
  p <- p / rowSums(p)
  p[,x]
}

# distribution functions for the categorical family
#
# @param q positive integers not greater than cat
# @param mu the linear predictor (of length or ncol cat-1)  
# @param cat the number of categories
#
# @return probabilites P(x <= q)
pcategorical <- function(q, eta, cat, link = "logit") {
  p <- dcategorical(1:max(q), eta = eta, cat = cat, link = link)
  do.call(cbind, lapply(q, function(j) rowSums(as.matrix(p[,1:j]))))
}

# density of the cumulative distribution
# 
# @param x positive integers not greater than cat
# @param mu the linear predictor (of length or ncol cat-1)  
# @param cat the number of categories
#
# @return the probabilities of the values in x
dcumulative <- function(x, eta, cat, link = "logit") {
  if (is.null(dim(eta))) eta <- matrix(eta, nrow = 1)
  if (length(dim(eta)) != 2 || !is.numeric(eta)) 
    stop("eta must be a numeric vector or matrix")
  if (missing(cat)) cat <- ncol(eta)
  mu <- ilink(eta, link)
  p <- cbind(mu[,1], 
             if (cat > 2) sapply(2:(cat-1), function(k) mu[,k] - mu[,k-1]), 
             1 - mu[,cat-1])
  p[,x]
}

# density of the sratio distribution
# 
# @param same arguments as dcumulative
dsratio <- function(x, eta, cat, link = "logit") {
  if (is.null(dim(eta))) eta <- matrix(eta, nrow = 1)
  if (length(dim(eta)) != 2 || !is.numeric(eta)) 
    stop("eta must be a numeric vector or matrix")
  if (missing(cat)) cat <- ncol(eta)
  mu <- ilink(eta, link)
  p <- cbind(mu[,1], 
        if (cat > 2) sapply(2:(cat-1), function(k)
          (mu[,k]) * apply(as.matrix(1-mu[,1:(k-1)]), 1, prod)),
        apply(1-mu, 1, prod))
  p[,x]
}

# density of the cratio distribution
# 
# @param same arguments as dcumulative
dcratio <- function(x, eta, cat, link = "logit") {
  if (is.null(dim(eta))) eta <- matrix(eta, nrow = 1)
  if (length(dim(eta)) != 2 || !is.numeric(eta)) 
    stop("eta must be a numeric vector or matrix")
  if (missing(cat)) cat <- ncol(eta)
  mu <- ilink(eta, link)
  p <- cbind(1 - mu[,1], 
        if (cat > 2) sapply(2:(cat-1), function(k)
          (1 - mu[,k]) * apply(as.matrix(mu[,1:(k-1)]), 1, prod)),
        apply(mu, 1, prod))
  p[,x]
}

# density of the acat family
# 
# @param same arguments as dcumulative
dacat <- function(x, eta, cat, link = "logit") {
  if (is.null(dim(eta))) eta <- matrix(eta, nrow = 1)
  if (length(dim(eta)) != 2 || !is.numeric(eta)) 
    stop("eta must be a numeric vector or matrix")
  if (missing(cat)) cat <- ncol(eta)
  
  if (link == "logit") {
    p <- cbind(rep(1, nrow(eta)), exp(eta[,1]), 
               matrix(NA, nrow = nrow(eta), ncol = cat - 2))
    if (cat > 2) 
      p[,3:cat] <- exp(sapply(3:cat, function(k) rowSums(eta[,1:(k-1)])))
  } 
  else {
    mu <- ilink(eta, link)
    p <- cbind(apply(1 - mu[,1:(cat-1)], 1, prod), 
               matrix(0, nrow = nrow(eta), ncol = cat - 1))
    if (cat > 2)
      p[,2:(cat-1)] <- sapply(2:(cat-1), function(k) 
        apply(as.matrix(mu[,1:(k-1)]), 1, prod) * apply(as.matrix(1 - mu[,k:(cat-1)]), 1, prod))
    p[,cat] <- apply(mu[,1:(cat-1)], 1, prod)
  }
  p <- p / rowSums(p)
  p[,x]
}

# distribution functions for ordinal families
#
# @param q positive integers not greater than cat
# @param mu the linear predictor (of length or ncol cat-1)  
# @param cat the number of categories
#
# @return probabilites P(x <= q)
pordinal <- function(q, eta, cat, family, link = "logit") {
  p <- do.call(paste0("d", family), list(1:max(q), eta = eta, cat = cat, link = link))
  do.call(cbind, lapply(q, function(j) rowSums(as.matrix(p[,1:j]))))
}