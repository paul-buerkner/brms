#' @export
linear_predictor.brmsfit <- function(x, ...) {
  if (!is(x$fit, "stanfit") || !length(x$fit@sim)) 
    stop("The model does not contain posterior samples")
  n.samples <- nrow(posterior.samples(x, parameters = "^lp__$"))
  Y <- x$data$Y
  eta <- matrix(0, nrow = n.samples, ncol = length(Y))
  X <- x$data$X
  if (!is.null(X) && ncol(X) && x$family != "categorical") {
    b <- posterior.samples(x, parameters = "^b_[^\\[]+$")
    eta <- eta + fixef_predictor(X = X, b = b)  
  }

  group <- names(x$ranef)
  all_groups <- extract_effects(x$formula)$group #may contain the same group more than ones
  if (length(group)) {
    for (i in 1:length(group)) {
      if (any(grepl(paste0("^lev_"), names(x$data)))) { # implies brms > 0.4.2
        #create a single RE design matrix for every grouping factor
        Z <- do.call(cbind, lapply(which(all_groups == group[i]), function(k) 
                     get(paste0("Z_",k), x$data)))
        gf <- get(paste0("lev_",match(group[i], all_groups)), x$data)
      } else { # implies brms < 0.4.2
        Z <- get(paste0("Z_",group[i]), x$data)
        gf <- get(group[i], x$data)
      }
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
  if (x$family %in% c("cumulative", "cratio", "sratio", "acat")) {
    Intercept <- posterior.samples(x, "^b_Intercept\\[")
    if (!is.null(x$data$Xp) && ncol(x$data$Xp)) {
      p <- posterior.samples(x, paste0("^b_",colnames(x$data$Xp),"\\["))
      etap <- partial_predictor(Xp = x$data$Xp, p = p, max_obs = x$data$max_obs)
    }  
    else etap <- array(0, dim = c(dim(eta), x$data$max_obs-1))
    for (k in 1:(x$data$max_obs-1)) {
      etap[,,k] <- etap[,,k] + eta
      if (x$family %in% c("cumulative", "sratio")) etap[,,k] <-  Intercept[,k] - etap[,,k]
      else etap[,,k] <- etap[,,k] - Intercept[,k]
    }
    eta <- etap
  }
  else if (x$family == "categorical") {
    if (!is.null(x$data$X)) {
      p <- posterior.samples(x, parameters = "^b_")
      etap <- partial_predictor(x$data$X, p, x$data$max_obs)
    }
    else etap <- array(0, dim = c(dim(eta), x$data$max_obs-1))
    for (k in 1:(x$data$max_obs-1)) etap[,,k] <- etap[,,k] + eta
    eta <- etap
  }
  eta
}

# compute eta for fixed effects
#
# @param X fixed effects design matrix
# @param b fixed effects samples
# 
# @return linear predictor for fixed effects
fixef_predictor <- function(X, b) {
  as.matrix(b) %*% t(as.matrix(X))
}
  
# compute eta for random effects
#  
# @param Z random effects design matrix
# @param gf levels of grouping factor for each observation
# @param r random effects samples
#
# @return linear predictor for random effects
ranef_predictor <- function(Z, gf, r) {
  Z <- expand_matrix(Z, gf)
  nlevels <- length(unique(gf))
  sort_levels <- unlist(lapply(1:nlevels, function(n) seq(n, ncol(r), nlevels)))
  as.matrix(r[,sort_levels]) %*% t(Z)
}

#compute eta for moving average effects
#
# @param data the data initially passed to stan
# @param ma moving average samples 
# @param eta previous linear predictor samples
# @param link
#
# @param new linear predictor samples updated by moving average effects
ma_predictor <- function(data, ma, eta, link = "identity") {
  ma <- as.matrix(ma)
  K <- ncol(ma)
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

# compute etap for partial and categorical effects
# 
# @param Xp partial design matrix 
# @param p partial effects samples
# @param max_obs number of categories
#
# @return linear predictor of partial effects as a 3D array (not as a matrix)
partial_predictor <- function(Xp, p, max_obs) {
  max_obs <- max(max_obs)
  etap <- array(0, dim = c(nrow(p), nrow(Xp), max_obs-1))
  for (k in 1:(max_obs-1)) {
    etap[,,k] <- as.matrix(p[,seq(k, (max_obs-1)*ncol(Xp), max_obs-1)]) %*% t(as.matrix(Xp))
  }
  etap
}

# expand a matrix into a sparse matrix of higher dimension
# 
# @param A matrix to be expanded
# @param x levels to expand the matrix
#
# @return An expanded matrix of dimensions nrow(A) and ncol(A) * length(unique(x)) 
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
