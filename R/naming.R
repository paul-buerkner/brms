#  rename certain symbols in a character vector
rename <- function(names, symbols = NULL, subs = NULL, fixed = TRUE, check_dup = FALSE) {
  if (is.null(symbols))
    symbols <- c(" ", "(", ")", "[", "]", ",", "+", "-", "*", "/", "^", "=", "!=")
  if (is.null(subs))
    subs <- c(rep("", 6), "P", "M", "MU", "D", "E", "EQ", "NEQ")
  if (length(symbols) != length(subs)) 
    stop("length(symbols) != length(subs)")
  new.names <- names
  for (i in 1:length(symbols)) 
    new.names <- gsub(symbols[i], subs[i], new.names, fixed = fixed)
  dup <- duplicated(new.names)
  if (check_dup && any(dup)) 
    stop(paste0("Internal renaming of variables led to duplicated names. \n",
                "Occured for variables: ", paste(names[which(new.names %in% new.names[dup])], collapse = ", ")))
  new.names
}

#get correlation names
get.cor.names <- function(names, type = "cor", eval = TRUE, brackets = TRUE) {
  cor.names <- NULL
  if (length(names) > 1 && eval) {
    for (i in 2:length(names)) {
      for (j in 1:(i-1)) {
        if (brackets) cor.names <- c(cor.names, paste0(type,"(",names[j],",",names[i],")"))
        else cor.names <- c(cor.names, paste0(type,"_",names[j],"_",names[i]))
      }
    }
  }
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
  ee <- extract.effects(x$formula, family = x$family)
  change <- list()
  
  #find positions of parameters and define new names
  f <- colnames(x$data$X)
  if (length(f) && x$family != "categorical") 
    change[[length(change)+1]] <- list(pos = grepl("^b\\[", pars), 
                                       oldname = "b", 
                                       pnames = paste0("b_",f), 
                                       fnames = paste0("b_",f))
  if (is.formula(x$partial) || x$family == "categorical") {
    if (x$family == "categorical") p <- colnames(x$data$X)
    else p <- colnames(x$data$Xp)
    thres <- (max(x$data$max_obs) - 1)
    change[[length(change)+1]] <- list(pos = grepl("^bp\\[", pars), 
                                       oldname = "bp", 
                                       pnames = paste0("b_",p), 
                                       fnames = paste0("b_", sapply(p, function(p) 
                                         sapply(1:thres, function(i) paste0(p,"[",i,"]")))),
                                       dim = thres,
                                       sort = unlist(lapply(1:length(p), function(k) 
                                         seq(k, thres*length(p), length(p)))))
  }  
  group <- names(x$ranef)
  if (length(x$ranef)) {
    for (j in 1:length(x$ranef)) {
      change[[length(change)+1]] <- list(pos = grepl(paste0("^sd_",group[j],"(\\[|$)"), pars),
                                         oldname = paste0("sd_",group[j]),
                                         pnames = paste0("sd_",group[j],"_", x$ranef[[j]]),
                                         fnames = paste0("sd_",group[j],"_", x$ranef[[j]]))
      if (length(x$ranef[[j]]) > 1 && ee$cor[[j]]) {
        cor_names <- get.cor.names(x$ranef[[j]], type = paste0("cor_",group[j]), brackets = FALSE)
        change[[length(change)+1]] <- list(pos = grepl(paste0("^cor_",group[j],"(\\[|$)"), pars),
                                           oldname = paste0("cor_",group[j]),
                                           pnames = cor_names,
                                           fnames = cor_names) 
      }
    }
  }
  if (x$family %in% c("gaussian", "student", "cauchy") && !is.formula(ee$se)) {
   change[[length(change)+1]] <- list(pos = grepl("^sigma", pars), 
                                      oldname = "sigma",
                                      pnames = paste0("sigma_",ee$response),
                                      fnames = paste0("sigma_",ee$response))
    if (x$family == "gaussian" && length(ee$response) > 1) {
      rescor_names <- paste0("rescor_",unlist(lapply(2:length(ee$response), function(j) 
          lapply(1:(j-1), function(k) paste0(ee$response[k],"_",ee$response[j])))))
     change[[length(change)+1]] <- list(pos = grepl("^rescor\\[", pars), 
                                        oldname = "rescor",
                                        pnames = rescor_names,
                                        fnames = rescor_names)
    }
  } 
  
  #rename parameters
  if (length(change)) {
    for (c in 1:length(change)) {
      x$fit@sim$fnames_oi[change[[c]]$pos] <- change[[c]]$fnames
      for (i in 1:chains) {
        names(x$fit@sim$samples[[i]])[change[[c]]$pos] <- change[[c]]$fnames
        if (!is.null(change[[c]]$sort)) x$fit@sim$samples[[i]][change[[c]]$pos] <- 
            x$fit@sim$samples[[i]][change[[c]]$pos][change[[c]]$sort]
      }
      onp <- match(change[[c]]$oldname, names(x$fit@sim$dims_oi))
      x$fit@sim$dims_oi <- c(if (onp > 1) x$fit@sim$dims_oi[1:(onp-1)], 
                             setNames(lapply(change[[c]]$pnames, function(x) 
                               if (is.null(change[[c]]$dim)) numeric(0)
                               else change[[c]]$dim), 
                               change[[c]]$pnames),
                             x$fit@sim$dims_oi[(onp+1):length(x$fit@sim$dims_oi)])
    }
  }
  x$fit@sim$pars_oi <- names(x$fit@sim$dims_oi)
  x
}

#' Extract parameter names
#' 
#' Extract all parameter names for which priors may be specified
#' 
#' @param x An object of class \code{formula}
#' @inheritParams brm
#' @param internal A flag indicating if the names of additional internal parameters should be displayed. 
#'   Setting priors on these parameters is not recommended
#' @param ... Currently ignored
#' 
#' @return A list of character vectors containing the parameter names for which priors may be specified
#' 
#' @examples 
#' par.names(rating ~ treat + period + carry + (1+carry|subject), 
#'           data = inhaler, family = "student")
#'           
#' par.names(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|visit),
#'           data = epilepsy, family = "poisson")          
#' 
#' @export
par.names.formula <- function(x, data = NULL, family = "gaussian", autocor = NULL, 
                              partial = NULL, threshold = "flexible", internal = FALSE, ...) {
  if (is.null(autocor)) autocor <- cor.arma()
  if (!is(autocor, "cor.brms")) stop("cor must be of class cor.brms")
  if (!threshold %in% c("flexible","equidistant")) 
    stop("threshold must be either flexible or equidistant")
  family <- family[1]
  ee <- extract.effects(x, family = family)
  data <- update_data(data, family = family, effects = ee)
  out <- list(fixef = paste0("b_",colnames(brm.model.matrix(ee$fixed, data = data))),
              ranef = list(), other = NULL)
  if (is.formula(partial)) 
    out$fixef <- c(out$fixef, colnames(brm.model.matrix(partial, data = data, rm.int = TRUE)))
  if (length(ee$group)) {
    gs <- unlist(ee$group)
    for (i in 1:length(gs)) {
      ranef <- colnames(brm.model.matrix(ee$random[[i]], data = data))
      out$ranef[[gs[i]]] <- c(paste0("sd_",gs[i],"_",ranef),
                              if (ee$cor[[i]] && length(ranef) > 1) 
                                c(paste0("cor_",gs[i]), if(internal) paste0("L_",gs[i]))) 
    }
  }
  if (is(autocor, "cor.arma") && autocor$p) out$other <- c(out$other, "ar")
  if (is(autocor, "cor.arma") && autocor$q) out$other <- c(out$other, "ma")
  if (family %in% c("gaussian", "student", "cauchy") && !is.formula(ee$se))
    out$other <- c(out$other, paste0("sigma_",ee$response))
  if (family == "gaussian" && length(ee$response) > 1)
    out$other <- c(out$other, "rescor", if (internal) "Lrescor")
  if (family == "student") out$other <- c(out$other, "nu")
  if (family %in% c("gamma", "weibull", "negbinomial")) 
    out$other <- c(out$other, "shape")
  if (family %in% c("cumulative", "sratio", "cratio", "acat") && threshold == "equidistant")
    out$other <- c(out$other, "delta")
  out
}