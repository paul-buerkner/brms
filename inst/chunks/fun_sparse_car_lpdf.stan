 /* Return the log probability of a proper conditional autoregressive (CAR)
  * prior with a sparse representation for the adjacency matrix
  * Full credit to Max Joseph (https://github.com/mbjoseph/CARstan)
  * Args:
  *   phi: Vector containing the CAR parameters for each location
  *   car: Dependence (usually spatial) parameter for the CAR prior
  *   sdcar: Standard deviation parameter for the CAR prior
  *   Nloc: Number of locations
  *   Nedges: Number of edges (adjacency pairs)
  *   Nneigh: Number of neighbors for each location
  *   eigenW: Eigenvalues of D^(-1/2) * W * D^(-1/2)
  *   edges1, edges2: Sparse representation of adjacency matrix
  * Details:
  *   D = Diag(Nneigh)
  * Returns:
  *   Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector phi, real car, real sdcar, int Nloc,
                       int Nedges, data vector Nneigh, data vector eigenW,
                       array[] int edges1, array[] int edges2) {
    real tau;  // precision parameter
    row_vector[Nloc] phit_D;  // phi' * D
    row_vector[Nloc] phit_W;  // phi' * W
    vector[Nloc] ldet;
    tau = inv_square(sdcar);
    phit_D = (phi .* Nneigh)';
    phit_W = rep_row_vector(0, Nloc);
    for (i in 1:Nedges) {
      phit_W[edges1[i]] = phit_W[edges1[i]] + phi[edges2[i]];
      phit_W[edges2[i]] = phit_W[edges2[i]] + phi[edges1[i]];
    }
    for (i in 1:Nloc) {
      ldet[i] = log1m(car * eigenW[i]);
    }
    return 0.5 * (Nloc * log(tau) + sum(ldet) -
           tau * (phit_D * phi - car * (phit_W * phi)));
  }
