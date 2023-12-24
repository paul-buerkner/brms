  /* grouped data stored linearly in "data" as indexed by begin and end
   * is repacked to be stacked into an array of vectors.
   */
  vector[] stack_vectors(vector long_data, int n, array[] int stack,
                         array[] int begin, array[] int end) {
    int S = sum(stack);
    int G = size(stack);
    array[S] vector[n] stacked;
    int j = 1;
    for (i in 1:G) {
      if (stack[i] == 1) {
        stacked[j] = long_data[begin[i]:end[i]];
        j += 1;
      }
    }
    return stacked;
  }
