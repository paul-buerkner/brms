  real monotonous(vector scale, int i) {
    if (i == 0) {
      return 0;
    } else {
      return(sum(scale[1:i]));
    }
  }
