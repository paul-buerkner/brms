/* Weibull AFT mixture cure model
 * identity parameterization of the incidence part
 */
real weibull_mixcure_lpdf(real y, real mu, real shape, real inc) {
    return bernoulli_lpmf(1 | inc) + weibull_lpdf(y | shape, mu);
}
real weibull_mixcure_lccdf(real y, real mu, real shape, real inc) {
    return log_sum_exp(
        bernoulli_lpmf(0 | inc),
        bernoulli_lpmf(1 | inc) + weibull_lccdf(y | shape, mu)
    );
}
real weibull_mixcure_lcdf(real y, real mu, real shape, real inc) {
    return log1m_exp(weibull_mixcure_lccdf(y | mu, shape, inc));
}
/* Weibull AFT mixture cure model
 * logit parameterization of the incidence part
 */
real weibull_mixcure_logit_lpdf(real y, real mu, real shape, real inc) {
    return bernoulli_logit_lpmf(1 | inc) + weibull_lpdf(y | shape, mu);
}
real weibull_mixcure_logit_lccdf(real y, real mu, real shape, real inc) {
    return log_sum_exp(
        bernoulli_logit_lpmf(0 | inc),
        bernoulli_logit_lpmf(1 | inc) + weibull_lccdf(y | shape, mu)
    );
}
real weibull_mixcure_logit_lcdf(real y, real mu, real shape, real inc) {
    return log1m_exp(weibull_mixcure_logit_lccdf(y | mu, shape, inc));
}