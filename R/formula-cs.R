#' Category Specific Predictors in \pkg{brms} Models
#'
#' @aliases cse
#'
#' @param expr Expression containing predictors,
#'  for which category specific effects should be estimated.
#'  For evaluation, \R formula syntax is applied.
#'
#' @details For detailed documentation see \code{help(brmsformula)}
#'   as well as \code{vignette("brms_overview")}.
#'
#' This function is almost solely useful when
#' called in formulas passed to the \pkg{brms} package.
#'
#' @seealso \code{\link{brmsformula}}
#'
#' @examples
#' \dontrun{
#' fit <- brm(rating ~ period + carry + cs(treat),
#'   data = inhaler, family = sratio("cloglog"),
#'   prior = set_prior("normal(0,5)"), chains = 2
#' )
#' summary(fit)
#' plot(fit, ask = FALSE)
#' }
#'
#' @export
cs <- function(expr) {
  deparse_no_string(substitute(expr))
}

# alias of function 'cs' used in the JSS paper of brms
#' @export
cse <- function(expr) {
  deparse_no_string(substitute(expr))
}
