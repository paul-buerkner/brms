#' Infections in kidney patients
#'
#' @description This dataset, originally discussed in
#'   McGilchrist and Aisbett (1991), describes the first and second
#'   (possibly right censored) recurrence time of
#'   infection in kidney patients using portable dialysis equipment.
#'   In addition, information on the risk variables age, sex and disease
#'   type is provided.
#'
#' @format A data frame of 76 observations containing
#'   information on the following 7 variables.
#' \describe{
#'  \item{time}{The time to first or second recurrence of the infection,
#'    or the time of censoring}
#'  \item{recur}{A factor of levels \code{1} or \code{2}
#'    indicating if the infection recurred for the first
#'    or second time for this patient}
#'  \item{censored}{Either \code{0} or \code{1}, where \code{0} indicates
#'    no censoring of recurrence time and \code{1} indicates right censoring}
#'  \item{patient}{The patient number}
#'  \item{age}{The age of the patient}
#'  \item{sex}{The sex of the patient}
#'  \item{disease}{A factor of levels \code{other, GN, AN},
#'    and \code{PKD} specifying the type of disease}
#' }
#'
#' @source McGilchrist, C. A., & Aisbett, C. W. (1991).
#'   Regression with frailty in survival analysis.
#'   \emph{Biometrics}, 47(2), 461-466.
#'
#' @examples
#' \dontrun{
#' ## performing surivival analysis using the "weibull" family
#' fit1 <- brm(time | cens(censored) ~ age + sex + disease,
#'   data = kidney, family = weibull, init = "0"
#' )
#' summary(fit1)
#' plot(fit1)
#'
#' ## adding random intercepts over patients
#' fit2 <- brm(time | cens(censored) ~ age + sex + disease + (1 | patient),
#'   data = kidney, family = weibull(), init = "0",
#'   prior = set_prior("cauchy(0,2)", class = "sd")
#' )
#' summary(fit2)
#' plot(fit2)
#' }
#'
"kidney"


#' Clarity of inhaler instructions
#'
#' @description Ezzet and Whitehead (1991) analyze data from a two-treatment,
#'   two-period crossover trial to compare 2 inhalation devices for
#'   delivering the drug salbutamol in 286 asthma patients.
#'   Patients were asked to rate the clarity of leaflet instructions
#'   accompanying each device, using a 4-point ordinal scale.
#'
#' @format A data frame of 572 observations containing
#'   information on the following 5 variables.
#' \describe{
#'  \item{subject}{The subject number}
#'  \item{rating}{The rating of the inhaler instructions
#'    on a scale ranging from 1 to 4}
#'  \item{treat}{A contrast to indicate which of
#'    the two inhaler devices was used}
#'  \item{period}{A contrast to indicate the time of administration}
#'  \item{carry}{A contrast to indicate possible carry over effects}
#' }
#'
#' @source Ezzet, F., & Whitehead, J. (1991).
#'   A random effects model for ordinal responses from a crossover trial.
#'   \emph{Statistics in Medicine}, 10(6), 901-907.
#'
#' @examples
#' \dontrun{
#' ## ordinal regression with family "sratio"
#' fit1 <- brm(rating ~ treat + period + carry,
#'   data = inhaler, family = sratio(),
#'   prior = set_prior("normal(0,5)")
#' )
#' summary(fit1)
#' plot(fit1)
#'
#' ## ordinal regression with family "cumulative"
#' ## and random intercept over subjects
#' fit2 <- brm(rating ~ treat + period + carry + (1 | subject),
#'   data = inhaler, family = cumulative(),
#'   prior = set_prior("normal(0,5)")
#' )
#' summary(fit2)
#' plot(fit2)
#' }
#'
"inhaler"


#' Epileptic seizure counts
#'
#' @description Breslow and Clayton (1993) analyze data initially
#'   provided by Thall and Vail (1990) concerning
#'   seizure counts in a randomized trial of anti-convulsant
#'   therapy in epilepsy. Covariates are treatment,
#'   8-week baseline seizure counts, and age of the patients in years.
#'
#' @format A data frame of 236 observations containing information
#'   on the following 9 variables.
#' \describe{
#'  \item{Age}{The age of the patients in years}
#'  \item{Base}{The seizure count at 8-weeks baseline}
#'  \item{Trt}{Either \code{0} or \code{1} indicating
#'    if the patient received anti-convulsant therapy}
#'  \item{patient}{The patient number}
#'  \item{visit}{The session number from \code{1} (first visit)
#'    to \code{4} (last visit)}
#'  \item{count}{The seizure count between two visits}
#'  \item{obs}{The observation number, that is
#'    a unique identifier for each observation}
#'  \item{zAge}{Standardized \code{Age}}
#'  \item{zBase}{Standardized \code{Base}}
#' }
#'
#' @source Thall, P. F., & Vail, S. C. (1990).
#'  Some covariance models for longitudinal count data with overdispersion.
#'  \emph{Biometrics, 46(2)}, 657-671. \cr
#'
#' Breslow, N. E., & Clayton, D. G. (1993).
#'  Approximate inference in generalized linear mixed models.
#'  \emph{Journal of the American Statistical Association}, 88(421), 9-25.
#'
#' @examples
#' \dontrun{
#' ## poisson regression without random effects.
#' fit1 <- brm(count ~ zAge + zBase * Trt,
#'   data = epilepsy, family = poisson()
#' )
#' summary(fit1)
#' plot(fit1)
#'
#' ## poisson regression with varying intercepts of patients
#' ## as well as normal priors for overall effects parameters.
#' fit2 <- brm(count ~ zAge + zBase * Trt + (1 | patient),
#'   data = epilepsy, family = poisson(),
#'   prior = set_prior("normal(0,5)")
#' )
#' summary(fit2)
#' plot(fit2)
#' }
#'
"epilepsy"

#' Cumulative Insurance Loss Payments
#'
#' @description This dataset, discussed in Gesmann & Morris (2020), contains
#'   cumulative insurance loss payments over the course of ten years.
#'
#' @format A data frame of 55 observations containing information
#'   on the following 4 variables.
#' \describe{
#'  \item{AY}{Origin year of the insurance (1991 to 2000)}
#'  \item{dev}{Deviation from the origin year in months}
#'  \item{cum}{Cumulative loss payments}
#'  \item{premium}{Achieved premiums for the given origin year}
#' }
#'
#' @source Gesmann M. & Morris J. (2020). Hierarchical Compartmental Reserving
#'   Models. \emph{CAS Research Papers}.
#'
#' @examples
#' \dontrun{
#' # non-linear model to predict cumulative loss payments
#' fit_loss <- brm(
#'   bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
#'     ult ~ 1 + (1 | AY), omega ~ 1, theta ~ 1,
#'     nl = TRUE
#'   ),
#'   data = loss, family = gaussian(),
#'   prior = c(
#'     prior(normal(5000, 1000), nlpar = "ult"),
#'     prior(normal(1, 2), nlpar = "omega"),
#'     prior(normal(45, 10), nlpar = "theta")
#'   ),
#'   control = list(adapt_delta = 0.9)
#' )
#'
#' # basic summaries
#' summary(fit_loss)
#' conditional_effects(fit_loss)
#'
#' # plot predictions per origin year
#' conditions <- data.frame(AY = unique(loss$AY))
#' rownames(conditions) <- unique(loss$AY)
#' me_loss <- conditional_effects(
#'   fit_loss,
#'   conditions = conditions,
#'   re_formula = NULL, method = "predict"
#' )
#' plot(me_loss, ncol = 5, points = TRUE)
#' }
#'
"loss"
