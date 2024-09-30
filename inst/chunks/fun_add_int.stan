  /* add a single integer to an array of integers
   * Args:
   *   x: array of integers
   *   y: a single integer
   * Returns:
   *   an array of intergers of the same length as x
   */
  array[] int add_int(array[] int x, int y) {
    array[num_elements(x)] int out;
    for (n in 1:num_elements(x)) {
      out[n] = x[n] + y;
    }
    return out;
  }
