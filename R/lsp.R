# find all namespace entries of a package, which are of
# a particular type for instance all exported objects
# retrieved from https://github.com/raredd/rawr
# @param package the package name
# @param what type of the objects to retrieve ("all" for all objects)
# @param pattern regex that must be matches by the object names
# @return a character vector of object names
lsp <- function(package, what = "all", pattern = ".*") {
  if (!is.character(substitute(package)))
    package <- deparse0(substitute(package))
  ns <- asNamespace(package)

  ## base package does not have NAMESPACE
  if (isBaseNamespace(ns)) {
    res <- ls(.BaseNamespaceEnv, all.names = TRUE)
    return(res[grep(pattern, res, perl = TRUE, ignore.case = TRUE)])
  } else {
    ## for non base packages
    if (exists('.__NAMESPACE__.', envir = ns, inherits = FALSE)) {
      wh <- get('.__NAMESPACE__.', inherits = FALSE,
                envir = asNamespace(package, base.OK = FALSE))
      what <- if (missing(what)) 'all'
      else if ('?' %in% what) return(ls(wh))
      else ls(wh)[pmatch(what[1], ls(wh))]
      if (!is.null(what) && !any(what %in% c('all', ls(wh))))
        stop('\'what\' should be one of ',
             paste0(shQuote(ls(wh)), collapse = ', '),
             ', or \'all\'', domain = NA)
      res <- sapply(ls(wh), function(x) getNamespaceInfo(ns, x))
      res <- rapply(res, ls, classes = 'environment',
                    how = 'replace', all.names = TRUE)
      if (is.null(what))
        return(res[grep(pattern, res, perl = TRUE, ignore.case = TRUE)])
      if (what %in% 'all') {
        res <- ls(getNamespace(package), all.names = TRUE)
        return(res[grep(pattern, res, perl = TRUE, ignore.case = TRUE)])
      }
      if (any(what %in% ls(wh))) {
        res <- res[[what]]
        return(res[grep(pattern, res, perl = TRUE, ignore.case = TRUE)])
      }
    } else stop(sprintf('no NAMESPACE file found for package %s', package))
  }
}
