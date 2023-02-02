/* mixcure lognormal (AFT) log-PDF of a single response
 * identity parameterization of the incidence part
 */
real mixcure_lognormal_lpdf(real y, real mu, real sigma, real inc) {
    return bernoulli_lpmf(1 | inc) + lognormal_lpdf(y | mu, sigma);
}
real mixcure_lognormal_lccdf(real y, real mu, real sigma, real inc) {
    return log_sum_exp(
        bernoulli_lpmf(0 | inc),
        bernoulli_lpmf(1 | inc) + lognormal_lccdf(y | mu, sigma)
    );
}
real mixcure_lognormal_lcdf(real y, real mu, real sigma, real inc) {
    return log1m_exp(mixcure_lognormal_lccdf(y | mu, sigma, inc));
}
/* mixcure lognormal (AFT) log-PDF of a single response
 * logit parameterization of the incidence part
 */
real mixcure_lognormal_logit_lpdf(real y, real mu, real sigma, real inc) {
    return bernoulli_logit_lpmf(1 | inc) + lognormal_lpdf(y | mu, sigma);
}
real mixcure_lognormal_logit_lccdf(real y, real mu, real sigma, real inc) {
    return log_sum_exp(
        bernoulli_logit_lpmf(0 | inc),
        bernoulli_logit_lpmf(1 | inc) + lognormal_lccdf(y | mu, sigma)
    );
}
real mixcure_lognormal_logit_lcdf(real y, real mu, real sigma, real inc) {
    return log1m_exp(mixcure_lognormal_logit_lccdf(y | mu, sigma, inc));
}
