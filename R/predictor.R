#' @export
linear.predictor.brmsfit <- function(x, ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  n.samples <- nrow(posterior.samples(x, parameters = "^lp__$"))
  Y <- x$data$Y
  eta <- matrix(0, nrow = n.samples, ncol = length(Y))
  X <- x$data$X
  if (!is.null(X)) {
    b <- posterior.samples(x, parameters = "^b_[^\\[]+$")
    eta <- eta + fixef_predictor(X = X, b = b)  
  }
  group <- names(x$ranef)
  if (length(group)) {
    for (i in 1:length(group)) {
      Z <- get(paste0("Z_",group[i]), x$data)
      gf <- get(group[i], x$data)
      r <- posterior.samples(x, parameters = paste0("^r_",group[i],"\\["))
      eta <- eta + ranef_predictor(Z = Z, gf = gf, r = r) 
    }
  }
  if (x$autocor$p > 0) {
    Yar <- as.matrix(x$data$Yar)
    ar <- posterior.samples(x, parameters = "^ar\\[")
    eta <- eta + fixef_predictor(X = Yar, b = ar)
  }
  if (x$autocor$q > 0) {
    ma <- posterior.samples(x, parameters = "^ma\\[")
    eta <- ma_predictor(data = x$data, ma = ma, eta = eta, link = x$link)
  }
  eta
}

#compute eta for fixed effects
fixef_predictor <- function(X, b) {
  as.matrix(b) %*% t(as.matrix(X))
}
  
#compute eta for random effects
ranef_predictor <- function(Z, gf, r) {
  Z <- expand_matrix(Z, gf)
  nlevels <- length(unique(gf))
  sort_levels <- unlist(lapply(1:nlevels, function(n) seq(n, ncol(r), nlevels)))
  as.matrix(r[,sort_levels]) %*% t(Z)
}

#compute eta for moving average effects
ma_predictor <- function(data, ma, eta, link = "identity") {
  ma <- as.matrix(ma)
  K <- 1:ncol(ma)
  Ks <- 1:K
  Y <- link(data$Y, link)
  N <- length(Y)
  tg <- c(rep(0, K), data$tgroup)
  Ema <- array(0, dim = c(nrow(ma), K, N))
  e <- matrix(0, nrow = nrow(ma), ncol = N)
  for (n in 1:N) {
    eta[,n] <- eta[,n] + apply(ma * Ema[,,n], 1, sum)
    e[,n] <- Y[n] - eta[,n]
    if (n < N) {
      I <- which(n < N & tg[n+1+K] == tg[n+1+K-Ks])
      Ema[,I,n+1] <- e[,n+1-I]
    }
  }
  eta
}

#expand a matrix into a sparse matrix of higher dimension
expand_matrix <- function(A, x) {
  A <- as.matrix(A)
  if (length(x) != nrow(A))
    stop("x must have nrow(A) elements")
  if (!all(is.wholenumber(x) & x > 0))
    stop("x must contain positive integers only")
  K <- ncol(A)
  v <- rep(0, K * max(x))
  do.call(rbind, lapply(1:nrow(A), function(n, v) {
    v[K*(x[n]-1) + 1:K] <- A[n, ] 
    return(v)}, v = v))
}
