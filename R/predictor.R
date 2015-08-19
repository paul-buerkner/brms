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
  eta
}

#calculate eta for fixed effects
fixef_predictor <- function(X, b) {
  X <- as.matrix(X)
  do.call(rbind, lapply(1:nrow(b), function(i) t(X %*% as.numeric(b[i,]))))
}
  
#calculate eta for random effects
ranef_predictor <- function(Z, gf, r) {
  Z <- as.matrix(Z)
  N <- nrow(Z)
  K <- ncol(Z)
  nlevels <- length(unique(gf))
  sort_levels <- unlist(lapply(1:nlevels, function(n) seq(n, ncol(r), nlevels)))
  r <- r[,sort_levels]
  index <- lapply(1:N, function(n) K*(gf[n]-1) + 1:K)
  do.call(rbind, lapply(1:nrow(r), function(i) {
    r_row <- as.numeric(r[i, ])
    sapply(1:N, function(n) Z[n,] %*% r_row[index[[n]]])}))
}
  
