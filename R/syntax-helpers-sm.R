# This file contains functions dealing with the extended 
# formula syntax to specify smooth terms via mgcv

# extract information about smooth terms
# @param x either a formula or a list containing an element "sm"
# @param data data.frame containing the covariates
tidy_smef <- function(x, data) {
  if (is.formula(x)) {
    x <- parse_bf(x, check_response = FALSE)$dpars$mu
  }
  form <- x[["sm"]] 
  if (!is.formula(form)) {
    return(empty_data_frame())
  }
  out <- data.frame(term = all_terms(form), stringsAsFactors = FALSE)
  nterms <- nrow(out)
  out$sfun <- get_matches("^[^\\(]+", out$term)
  out$vars <- out$byvars <- out$covars <- vector("list", nterms)
  for (i in seq_len(nterms)) {
    sm <- eval2(out$term[i])
    out$covars[[i]] <- sm$term
    if (sm$by != "NA") {
      out$byvars[[i]] <- sm$by
    }
    out$vars[[i]] <- c(out$covars[[i]], out$byvars[[i]])
  }
  out$label <- paste0(out$sfun, rename(ulapply(out$vars, collapse)))
  # prepare information inferred from the data
  sdata <- data_sm(x, data, knots = attr(data, "knots"))
  bylevels <- attr(sdata$Xs, "bylevels")
  nby <- lengths(bylevels)
  tmp <- vector("list", nterms)
  for (i in seq_len(nterms)) {
    tmp[[i]] <- out[i, , drop = FALSE]
    tmp[[i]]$termnum <- i
    if (nby[i] > 0L) {
      tmp[[i]] <- do_call(rbind, repl(tmp[[i]], nby[i]))
      tmp[[i]]$bylevel <- rm_wsp(bylevels[[i]])
      tmp[[i]]$byterm <- paste0(tmp[[i]]$term, tmp[[i]]$bylevel)
      str_add(tmp[[i]]$label) <- rename(tmp[[i]]$bylevel)
    } else {
      tmp[[i]]$bylevel <- NA
      tmp[[i]]$byterm <- tmp[[i]]$term
    }
  }
  out <- do_call(rbind, tmp)
  out$knots <- sdata[grepl("^knots_", names(sdata))]
  out$nbases <- lengths(out$knots)
  attr(out, "Xs_names") <- colnames(sdata$Xs)
  rownames(out) <- NULL
  out
}

# check if smooths are present in the model
has_smooths <- function(bterms) {
  length(get_effect(bterms, target = "sm")) > 0L
}
