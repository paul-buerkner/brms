  /* how many elements are in a range of integers?
   * Args:
   *   x: an integer array
   *   start: start of the range (inclusive)
   *   end: end of the range (inclusive)
   * Returns:
   *   a scalar integer
   */
  int size_range(int[] x, int start, int end) {
    int out = 0;
    for (i in 1:size(x)) {
      out += (x[i] >= start && x[i] <= end);
    }
    return out;
  }
  /* which elements are in a range of integers?
   * Args:
   *   x: an integer array
   *   start: start of the range (inclusive)
   *   end: end of the range (inclusive)
   * Returns:
   *   an integer array
   */
  int[] which_range(int[] x, int start, int end) {
    array[size_range(x, start, end)] int out;
    int j = 1;
    for (i in 1:size(x)) {
      if (x[i] >= start && x[i] <= end) {
        out[j] = i;
        j += 1;
      }
    }
    return out;
  }
  /* adjust array values to x - start + 1
   * Args:
   *   x: an integer array
   *   start: start of the range of values in x (inclusive)
   * Returns:
   *   an integer array
   */
  int[] start_at_one(int[] x, int start) {
    array[size(x)] int out;
    for (i in 1:size(x)) {
      out[i] = x[i] - start + 1;
    }
    return out;
  }
