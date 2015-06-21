#' Infections in kidney patients
#' 
#' @description This dataset, originally discussed in McGilchrist and Aisbett (1991), describes the first and second (possibly right censored) recurrence time of 
#' infection in kideny patients using portable dialysis equipment. In addition, information on the risk variables age, sex and disease type is provided.
#' 
#' @format A dataframe of 76 observations containing information on the following 7 variables.
#' \describe{
#'  \item{time}{The time to first or second recurrence of the infection, or the time of censoring}
#'  \item{recur}{A factor of levels \code{1} or \code{2} indicating if the infection recurred for the first or second time for this patient}
#'  \item{censored}{Either \code{0} or \code{1}, where \code{0} indicates no censoring of recurrence time and \code{1} indicates right censoring} 
#'  \item{patient}{The patient number}
#'  \item{age}{The age of the patient}
#'  \item{sex}{The sex of the patient}
#'  \item{disease}{A factor of levels \code{other, GN, AN}, and \code{PKD} specifiying the type of disease}
#' }
#' 
#' @examples 
#' \dontrun{
#' ## performing surivival analysis using the "weibull" family
#' ## time | cens indicates which values in variable time are right censored
#' fit_k1 <- brm(time | cens(censored) ~ age + sex + disease, data = kidney, 
#'               family = "weibull")
#' summary(fit_k1) 
#' plot(fit_k1)
#' 
#' ## adding random intercepts over patients and using weakly informative priors 
#' ## for regression parameters and standard deviations of random effects 
#' fit_k2 <- brm(time | cens(censored) ~ age + sex + disease + (1|patient), 
#'               data = kidney, prior = list(sd = "uniform(0,20)"), 
#'               family = "weibull", n.iter = 3000, silent = TRUE)
#' summary(fit_k2) 
#' plot(fit_k2)         
#' }
#' 
#' @source McGilchrist, C. A., & Aisbett, C. W. (1991). Regression with frailty in survival analysis. 
#'   \emph{Biometrics, 47(2)}, 461-466.
"kidney"


#' Clarity of inhaler instructions
#' 
#' @description Ezzet and Whitehead (1991) analyse data from a two-treatment, two-period crossover trial to compare 2 inhalation devices for 
#'   delivering the drug salbutamol in 286 asthma patients. Patients were asked to rate the clarity of leaflet instructions accompanying each device, 
#'   using a 4-point ordinal scale.
#'   
#' @format A dataframe of 572 observations containing information on the following 5 variables. 
#' \describe{
#'  \item{subject}{The subject number}
#'  \item{rating}{The rating of the inhaler instructions on a scale ranging from 1 to 4}
#'  \item{treat}{A contrast to indicate which of the two inhaler devices was used} 
#'  \item{period}{A contrast to indicate the time of administration}
#'  \item{carry}{A contrast to indicate possible carry over effects}
#' } 
#' 
#' @examples
#' \dontrun{
#' ## ordinal regression with family "sratio"
#' fit_i1 <- brm(rating ~ treat + period + carry, data = inhaler, 
#'               family = "sratio", prior = list(b = "normal(0,5)"))
#' summary(fit_i1) 
#' plot(fit_i1)
#'        
#' ## ordinal regression with family "cumulative" and random intercept over subjects             
#' fit_i2 <- brm(rating ~ treat + period + carry + (1|subject), data = inhaler, 
#'               family = "cumulative", prior = list(b = "normal(0,5)"))
#' summary(fit_i2) 
#' plot(fit_i2)
#' }
#' 
#' @source Ezzet, F., & Whitehead, J. (1991). A random effects model for ordinal responses from a crossover trial. 
#'   \emph{Statistics in Medicine, 10(6)}, 901-907.
"inhaler"


#' Epileptic seizure counts
#' 
#' @description Breslow and Clayton (1993) analyse data initially provided by Thall and Vail (1990) concerning 
#'   seizure counts in a randomised trial of anti-convulsant therapy in epilepsy. Covariates are treatment, 
#'   8-week baseline seizure counts, and age of the patients in years. 
#'   
#' @format A dataframe of 236 observations containing information on the following 9 variables. 
#' \describe{
#'  \item{Age}{The age of the patients in years}
#'  \item{Base}{The seizure count at 8-weeks baseline}
#'  \item{Trt}{Either \code{0} or \code{1} indicating if the patient recieved anti-convulsant therapy} 
#'  \item{log_Age_c}{The logarithm of Age centered arounds its mean}
#'  \item{log_Base4_c}{The logarithm of Base divided by 4 (i.e. log(Base/4)) centered around its mean}
#'  \item{Trt_c}{Trt centered around its mean}
#'  \item{visit}{The session number from \code{1} (first visit) to \code{4} (last visit)}
#'  \item{count}{The seizure count between two visits}
#'  \item{patient}{The patient number}
#' } 
#' 
#' @examples
#' \dontrun{
#' ## poisson regression without random effects. 
#' ## family = c("poisson", "log") is equivalent to family = "poisson"
#' fit_e1 <- brm(count ~ log_Age_c + log_Base4_c * Trt_c, 
#'             data = epilepsy, family = c("poisson", "log"))
#' summary(fit_e1) 
#' plot(fit_e1)             
#'     
#' ## poisson regression with random intercepts over patients and visits
#' ## as well as normal priors for fixed effects parameters.    
#' fit_e2 <- brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|visit), 
#'             data = epilepsy, family = "poisson", prior = list(b = "normal(0,5)"))
#' summary(fit_e2) 
#' plot(fit_e2)
#' }
#'  
#' @source Thall, P. F., & Vail, S. C. (1990). Some covariance models for longitudinal count data with overdispersion. 
#'    \emph{Biometrics, 46(2)}, 657-671. \cr
#'    
#' Breslow, N. E., & Clayton, D. G. (1993). Approximate inference in generalized linear mixed models. 
#'    \emph{Journal of the American Statistical Association, 88(421)}, 9-25.
"epilepsy"