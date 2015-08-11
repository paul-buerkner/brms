#calculate estimates over posterior samples 
get.estimate <- function(coef, samples, margin = 2, to.array = FALSE, ...) {
  dots <- list(...)
  args <- list(X = samples, MARGIN = margin, FUN = coef)
  fun.args <- names(formals(coef))
  if (!"..." %in% fun.args)
    dots <- dots[names(dots) %in% fun.args]
  x <- do.call(apply, c(args, dots))
  if (is.null(dim(x))) 
    x <- matrix(x, dimnames = list(NULL, coef))
  else if (coef == "quantile") x <- aperm(x, length(dim(x)):1)
  if (to.array && length(dim(x)) == 2) 
    x <- array(x, dim = c(dim(x), 1), dimnames = list(NULL, NULL, coef))
  x 
}

#calculate the evidence ratio between two disjunct hypotheses
eratio <- function(x, cut = 0, wsign = c("equal", "less", "greater"), prior_samples = NULL, pow = 12, ...) {
  wsign <- match.arg(wsign)
  if (wsign == "equal") 
    if (is.null(prior_samples)) out <- NA
    else {
      dots <- list(...)
      dots <- dots[names(dots) %in% names(formals("density.default"))]
      prior_density <- do.call(density, c(list(x = prior_samples, n = 2^pow), dots))
      posterior_density <- do.call(density, c(list(x = x, n = 2^pow), dots))
      at_cut_prior <- which(abs(prior_density$x - cut) == min(abs(prior_density$x - cut)))
      at_cut_posterior <- which(abs(posterior_density$x - cut) == min(abs(posterior_density$x - cut)))
      out <- posterior_density$y[at_cut_posterior] / prior_density$y[at_cut_prior] 
    }
  else if (wsign == "less") {
    out <- length(which(x < cut))
    out <- out/(length(x) - out)
  }  
  else if (wsign == "greater") {
    out <- length(which(x > cut))
    out <- out/(length(x) - out)
  }
  out  
}

#get correlation names
get.cor.names <- function(names, type = "cor", eval = TRUE) {
  cor.names <- NULL
  if (length(names) > 1 && eval)
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
  if (is.formula(x$partial) || x$family == "categorical") {
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
  if (x$family %in% c("gaussian", "student", "cauchy")) {
   change[[length(change)+1]] <- list(pos = grepl("^sigma", pars), names = paste0("sigma_",ee$response))
    if (x$family == "gaussian" && length(ee$response) > 1) {
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

#get appropriate prior samples for hypothesis testing
get_prior_samples <- function(x, pars) {
  if (!is(x, "brmsfit")) stop("x must be of class brmsfit")
  if (!is.character(pars)) stop("pars must be a character vector")
  par_names <- par.names(x)
  prior_names <- par_names[grepl("^prior_", par_names)]
  prior_samples <- posterior.samples(x, parameters = prior_names, fixed = TRUE)
  matches <- lapply(sub("^prior_", "", prior_names), regexpr, text = pars, fixed = TRUE)
  matches <- matrix(unlist(matches), ncol = length(pars), byrow = TRUE)
  matches <- apply(matches, 2, function(table) match(1, table))
  if (!anyNA(matches)) prior_samples <- as.data.frame(prior_samples[,matches])
  else prior_samples <- NULL
  prior_samples
}