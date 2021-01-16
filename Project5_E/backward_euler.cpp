#include <iostream>
#include <armadillo>
#include "sym_tridiag_solver.cpp"

using namespace std;
using namespace arma;


void diffusion_BE(mat &u, int N, int n, vec v, rowvec w, rowvec z, bool save = 0)
{

  double dx    = (double) 1/(N - 1) ;
  double dt    =          dx*dx/2   ;
  double alpha =          dt/(dx*dx);


  u.zeros();

  u(span::all, 0)     = v;
  u(0, span::all)     = w;
  u(N - 1, span::all) = z;

  vec x(N);

  for (int j = 1; j < n; ++j) {

    sym_tridiag_solver(x, u(span::all, j - 1), N, -alpha, 1 + 2*alpha);
    u(span::all, j) = x;

  }

  if (save) { u.save("data_BE.dat", raw_ascii) ;}
  return;
}
