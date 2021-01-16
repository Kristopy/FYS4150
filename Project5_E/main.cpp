#include <iostream>
#include <armadillo>
#include "backward_euler.h"
#include "sym_tridiag_solver.h"
using namespace std;
using namespace arma;

void diffusion_FE(mat &u, int N, int n, vec v, rowvec w, rowvec z, bool save);

int main()
{

  int N = 11;
  int n = 11;

  mat u(N, n); vec v(N); rowvec w(n); rowvec z(n);



  double dx    = (double) 1/(N - 1) ;
  double dt    =          dx*dx/2   ;
  double alpha =          dt/(dx*dx);
  vec x(N);
  for (int j = 1; j < n; ++j) {
  sym_tridiag_solver(x, u(span::all, j -1), N, -alpha, 1 + 2*alpha);
  cout << x << endl;
  }

  diffusion_FE(u, N, n, v.zeros(), w.zeros(), z.ones(), 1);
  diffusion_BE(u, N, n, v.zeros(), w.zeros(), z.ones(), 1);
  cout << "Done!" << endl;

  return 0;

}

void diffusion_FE(mat &u, int N, int n, vec v, rowvec w, rowvec z, bool save = 0)
{
  /*

    Solves the one-dimensional diffusion equation
    on the interval 0 <= x <= 1 for t >= 0
    using Forward Euler

    mat u   : N x n - matrix
    int N   : # of grid points
    int n   : # of time steps
    vec v   : Initial conditions u(x, t=0)
    rowvec w: Boundary condition u(x=0, t)
    rowvec z: Boundary condition u(x=1, t)

  */

  double dx    = (double) 1/(N - 1) ;
  double dt    =          dx*dx/2   ;
  double alpha =          dt/(dx*dx);

  u.zeros();

  u(span::all, 0)     = v;
  u(0, span::all)     = w;
  u(N - 1, span::all) = z;

  for (int j = 1; j < n - 1; ++j) {
  for (int i = 1; i < N - 1; ++i) {

    u(i, j+1) = alpha*(u(i+1, j) + u(i-1, j)) + (1 - 2*alpha)*u(i, j);

  }}

  if (save) { u.save("data_FE.dat", raw_ascii) ;}

  return;

}
