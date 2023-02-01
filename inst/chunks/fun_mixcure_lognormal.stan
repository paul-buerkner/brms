/* Lognormal AFT mixture cure model
 * identity parameterization of the incidence part
 */
real lognormal_mixcure_lpdf(real y, real mu, real sigma, real inc) {
    return bernoulli_lpmf(1 | inc) + lognormal_lpdf(y | mu, sigma);
}
real lognormal_mixcure_lccdf(real y, real mu, real sigma, real inc) {
    return log_sum_exp(
        bernoulli_lpmf(0 | inc),
        bernoulli_lpmf(1 | inc) + lognormal_lccdf(y | mu, sigma)
    );
}
real lognormal_mixcure_lcdf(real y, real mu, real sigma, real inc) {
    return log1m_exp(lognormal_mixcure_lccdf(y | mu, sigma, inc));
}
/* Lognormal AFT mixture cure model
 * logit parameterization of the incidence part
 */
real lognormal_mixcure_logit_lpdf(real y, real mu, real sigma, real inc) {
    return bernoulli_logit_lpmf(1 | inc) + lognormal_lpdf(y | mu, sigma);
}
real lognormal_mixcure_logit_lccdf(real y, real mu, real sigma, real inc) {
    return log_sum_exp(
        bernoulli_logit_lpmf(0 | inc),
        bernoulli_logit_lpmf(1 | inc) + lognormal_lccdf(y | mu, sigma)
    );
}
real lognormal_mixcure_logit_lcdf(real y, real mu, real sigma, real inc) {
    return log1m_exp(lognormal_mixcure_logit_lccdf(y | mu, sigma, inc));
}
