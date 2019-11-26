# has the model group specific ordinal thresholds
# @param x list with potentail $adforms elements
# @return TRUE or FALSE
has_grcat <- function(x) {
  cat <- x$adforms$cat
  isTRUE(!is.null(cat) && eval_rhs(cat)$vars$gr != "NA")
}

# extract variabe for group specific ordinal thresholds
# @param x list with potentail $adforms elements
# @param data data passed by the user
# @return a vector of threshold groups or NULL
extract_grcat <- function(x, data) {
  if (!has_grcat(x)) {
    return(NULL)
  }
  cat <- eval_rhs(x$adforms$cat)
  factor(eval2(cat$vars$gr, data))
}

