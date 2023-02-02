/* mixcure Weibull (AFT) log-PDF of a single response
 * identity parameterization of the incidence part
 */
real mixcure_weibull_lpdf(real y, real mu, real shape, real inc) {
    return bernoulli_lpmf(1 | inc) + weibull_lpdf(y | shape, mu);
}
real mixcure_weibull_lccdf(real y, real mu, real shape, real inc) {
    return log_sum_exp(
        bernoulli_lpmf(0 | inc),
        bernoulli_lpmf(1 | inc) + weibull_lccdf(y | shape, mu)
    );
}
real mixcure_weibull_lcdf(real y, real mu, real shape, real inc) {
    return log1m_exp(mixcure_weibull_lccdf(y | mu, shape, inc));
}
/* mixcure Weibull (AFT) log-PDF of a single response
 * logit parameterization of the incidence part
 */
real mixcure_weibull_logit_lpdf(real y, real mu, real shape, real inc) {
    return bernoulli_logit_lpmf(1 | inc) + weibull_lpdf(y | shape, mu);
}
real mixcure_weibull_logit_lccdf(real y, real mu, real shape, real inc) {
    return log_sum_exp(
        bernoulli_logit_lpmf(0 | inc),
        bernoulli_logit_lpmf(1 | inc) + weibull_lccdf(y | shape, mu)
    );
}
real mixcure_weibull_logit_lcdf(real y, real mu, real shape, real inc) {
    return log1m_exp(mixcure_weibull_logit_lccdf(y | mu, shape, inc));
}