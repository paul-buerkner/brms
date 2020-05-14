### Rudimentary support for brms objects in emmeans package. 
### Obviously this is way less than is needed, but it does support simpler models

# These are NOT exported here. They are dynamically exported if emmeans is installed.
# See code in zzz.R

# Note from RVL: It may be desirable to add optional arguments, especially to
# the emm_basis method, to allow special options -- e.g., for hurdle
# models, which part of the model to estimate. If you do that, obviously
# you will want to document it in some way.


# Recover the predictors used in the fixed-effects part of the model
recover_data.brmsfit <- function(object, data, ...) {
    bt <- brmsterms(formula(object))
    if (class(bt) != "brmsterms")
        stop("This model is currently not supported.")
    mt <- attr(model.frame(bt$dpars$mu$fe, data = object$data), "terms")
    ### we don't have a call component so I'll just put in a dummy
    emmeans::recover_data(call("brms"), mt, "na.omit", data = object$data, ...)
}

# Calculate the basis for making predictions. This is essentially the
# inside of the predict() function with new data on the linkm scale. 
# Transforming to response scale, if desired, is handled by emmeans
emm_basis.brmsfit <- function(object, trms, xlev, grid, vcov., ...) {
    m <- model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    contr <- lapply(object$data, function(.) attr(., "contrasts"))
    contr <- contr[!sapply(contr, is.null)]
    X <- model.matrix(trms, m, contrasts.arg = contr)
    nm <- rename(colnames(X))
    V <- vcov(object)[nm, nm, drop = FALSE]
    nbasis <- estimability::all.estble  # we're assuming no rank deficiency
    dfargs <- list()
    dffun <- function(k, dfargs) Inf
    misc <- .std.link.labels(brmsterms(formula(object))$dpars$mu$family, list())
    post.beta <- as.matrix(object, pars = paste0("b_", nm), fixed = TRUE)
    bhat <- apply(post.beta, 2, mean)
    
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, dfargs = dfargs, 
         misc = misc, post.beta = post.beta)
}
