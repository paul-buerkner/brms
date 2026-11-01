  /* compute log(Phi(y)) with exact values and gradients in both tails, via
   * erfc and, below y = -37, the asymptotic series of the Mills ratio; the
   * gradient of std_normal_lcdf is off by up to 6e-5 even in the bulk
   * Args:
   *   y: a scalar in (-Inf, Inf)
   * Returns:
   *   a scalar in (-Inf, 0)
   */
   real log_Phi(real y) {
     if (y < -37) {
       // 1 - z + 3 z^2 - 15 z^3 + ... - 1, with z = 1 / y^2
       real z = inv_square(y);
       real s = z * (-15 + z * (105 + z * (-945 + z * 10395)));
       s = z * (-1 + z * (3 + s));
       return -0.5 * square(y) - log(-y) - 0.5 * log(2 * pi()) + log1p(s);
     } else if (y < 0) {
       return log(erfc(-y / sqrt2())) - log2();
     }
     return log1m(0.5 * erfc(y / sqrt2()));
   }
