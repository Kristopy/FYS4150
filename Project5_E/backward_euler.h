#ifndef BACKWARD_EULER_H
#define BACKWARD_EULER_H
#include <armadillo>
using namespace arma;
void diffusion_BE(mat &u, int N, int n, vec v, rowvec w, rowvec z, bool save = 0);
#endif // BACKWARD_EULER_H
