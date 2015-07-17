#calculate estimates over posterior samples 
get.estimate <- function(coef, samples, margin = 2, to.array = FALSE, ...) {
  dots <- list(...)
  args <- list(X = samples, MARGIN = margin, FUN = coef)
  fun.args <- names(formals(coef))
  if (!"..." %in% fun.args) 
    dots <- dots[fun.args %in% names(dots)] 
  x <- do.call(apply, c(args, dots))
  if (is.null(dim(x))) 
    x <- matrix(x, dimnames = list(NULL, coef))
  else if (coef == "quantile") x <- aperm(x, length(dim(x)):1)
  if (to.array & length(dim(x)) == 2) 
    x <- array(x, dim = c(dim(x), 1), dimnames = list(NULL, NULL, coef))
  x 
}

#get correlation names
get.cor.names <- function(names, type = "cor", eval = TRUE) {
  cor.names <- NULL
  if (length(names) > 1 & eval)
    for (i in 2:length(names)) 
      for (j in 1:(i-1)) 
        cor.names <- c(cor.names, paste0(type,"(",names[j],",",names[i],")"))
  cor.names
}

#rename parameters
rename.pars <- function(x, ...) {
  if (!length(x$fit@sim)) return(x)
  chains <- length(x$fit@sim$samples) 
  n.pars <- length(x$fit@sim$fnames_oi)
  x$fit@sim$fnames_oi[1:(n.pars-1)] <- gsub("__", ":", x$fit@sim$fnames_oi[1:(n.pars-1)])
  for (i in 1:chains) names(x$fit@sim$samples[[i]]) <- x$fit@sim$fnames_oi
  pars <- dimnames(x$fit)$parameters
  change <- list()
  
  #find positions of parameters and define new names
  f <- colnames(x$data$X)
  if (length(f)) 
    change[[length(change)+1]] <- list(pos = grepl("^b\\[", pars), names = paste0("b_",f))
  if (is.formula(x$partial) | x$family == "categorical") {
    if (x$family == "categorical") p <- colnames(x$data$X)
    else p <- colnames(x$data$Xp)
   change[[length(change)+1]] <- list(pos = grepl("^bp\\[", pars), sort = TRUE,
      names = paste0("b_",sapply(1:(max(x$data$max_obs) - 1), function(i) 
                     sapply(p, function(p) paste0(p,"[",i,"]")))))
  }  
  group <- names(x$ranef)
  if (length(x$ranef)) {
    for (j in 1:length(x$ranef)) {
     change[[length(change)+1]] <- list(pos = grepl(paste0("^sd_",group[j],"(\\[|$)"), pars),
                                        names = paste0("sd_",group[j],"_", x$ranef[[j]]))
      cor_names <- unlist(lapply(1:length(group), function(i)
        if (length(x$ranef[[i]])>1) paste0("cor_",group[j],"_", unlist(lapply(2:length(x$ranef[[i]]), 
           function(j) lapply(1:(j-1), function(k) paste0(x$ranef[[i]][k],"_",x$ranef[[i]][j]))))))) 
     change[[length(change)+1]] <- list(pos = grepl(paste0("^cor_",group[j],"(\\[|$)"), pars),
                                        names = cor_names) 
    }
  }
  ee <- extract.effects(x$formula, family = x$family)
  if (x$family %in% c("gaussian", "student", "cauchy", "multigaussian")) {
   change[[length(change)+1]] <- list(pos = grepl("^sigma", pars), names = paste0("sigma_",ee$response))
    if (x$family == "multigaussian") {
      rescor_names <- paste0("rescor_",unlist(lapply(2:length(ee$response), function(j) 
          lapply(1:(j-1), function(k) paste0(ee$response[k],"_",ee$response[j])))))
     change[[length(change)+1]] <- list(pos = grepl("^rescor\\[", pars), names = rescor_names)
    }
  } 
  
  #rename parameters
  for (c in 1:length(change)) {
    sort <- !is.null(change[[c]]$sort)
    if (sort) x$fit@sim$fnames_oi[change[[c]]$pos] <- sort(change[[c]]$names)
    else x$fit@sim$fnames_oi[change[[c]]$pos] <- change[[c]]$names
    for (i in 1:chains) {
      names(x$fit@sim$samples[[i]])[change[[c]]$pos] <- change[[c]]$names
      if (sort) x$fit@sim$samples[[i]][change[[c]]$pos] <- 
          x$fit@sim$samples[[i]][change[[c]]$pos][order(change[[c]]$names)]
    }  
  }
  x
}