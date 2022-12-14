  /* grouped data stored linearly in "data" as indexed by begin and end
   * is repacked to be stacked into an array of vectors.
   */
  vector[] stack_vectors(vector long_data, int n, int[] stack,
                         int[] begin, int[] end) {
    int S = sum(stack);
    int G = size(stack);
    vector[n] stacked[S];
    int j = 1;
    for (i in 1:G) {
      if (stack[i] == 1) {
        stacked[j] = long_data[begin[i]:end[i]];
        j += 1;
      }
    }
    return stacked;
  }
