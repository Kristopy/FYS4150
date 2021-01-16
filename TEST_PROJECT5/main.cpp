#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>

using namespace  std;
using namespace arma;



void diffusion_FE(mat &u, int N, int T, double dx, rowvec v, vec w, vec z, bool save );
void diffusion_BE(mat &u, int N, int n, double dx, rowvec v, vec w, vec z, bool save );
void diffusion_CN(mat &u, int N, int n, double dx, rowvec v, vec w, vec z, bool save );

void tridiag(rowvec &x, rowvec y, int N, double b, double d);

int main(int argc, char* argv[]){

  // Step length in position x
  double dx;
  cout << "Choose step-length dx:" << endl;
  cin >> dx;

  // Step length in time
  double dt = 0.5*dx*dx;

  // Defining alpha
  double alpha = dt/dx/dx;

  // Number of integration points along x-axis (inner points only)
  int N = int(1.0/(dx));

  // Number of time steps till final time
  int T;
  cout << "Choose Total time iteration:" << endl;
  cin >> T;

  cout << "Parameters" << endl;
  cout << "dx = " << dx << endl;
  cout << "dt = " << dt << endl;
  cout << "alpha = " << alpha << endl;
  cout << "N = " << N << endl;

  // Defining u
  mat u = zeros<mat>(T,N);

  //Initial condition
  rowvec v(N); //Initial conditions u(x, t=0)
  rowvec x(N); //Goes in tridiag

  vec w(T); //Boundary condition u(x=0, t)
  vec z(T); //Boundary condition u(x=1, t)

  diffusion_FE(u, N, T, dx, v.zeros(),w.zeros(), z.ones(),1);
  cout << u << endl;
  cout << "-------------------------------------"<< endl;
  diffusion_BE(u, N, T, dx, v.zeros(),w.zeros(), z.ones(),1);
  cout << u << endl;
  cout << "-------------------------------------"<< endl;
  diffusion_CN(u, N, T, dx, v.zeros(),w.zeros(), z.ones(),1);
  cout << u << endl;
  return 0;
}

void tridiag(rowvec &x, rowvec y, int N, double b, double d)
{

  /*

    | 1            |    |x0  |   |y0  |
    | b d b        |    |x1  |   |y1  |
    |     ...      | *  |... | = |... |
    |         b d b|    |    |   |    |
    |             1|    |xN-1|   |yN-1|

  */
  x.set_size(N);

  double bb = b*b;

  vec diag(N); diag.fill(d);

  y(1) += -b*y(0);
  for (int i = 2; i < N - 1; ++i) {
    diag(i) += -bb/diag(i-1);
       y(i) += -b*y(i-1)/diag(i-1);
  }

  x(0) = y(0); x(N - 1) = y(N - 1);
  for (int i = N - 2; i > 0; --i) {
    x(i) = (y(i) - b*x(i + 1))/diag(i);
  }

  return;
}

void diffusion_FE(mat &u, int N, int T, double dx, rowvec v, vec w, vec z, bool save = 0)
{
  /*
    Solves the one-dimensional diffusion equation
    on the interval 0 <= x <= 1 for t >= 0
    using Forward Euler

    mat u   : T x N - matrix-- Time depicted as rows
    int N   : # of grid points
    int T   : # of time steps

    rowvec v   : Initial conditions u(x, t=0)
    vec w: Boundary condition u(x=0, t)
    vec z: Boundary condition u(x=1, t)
  */

  double dt    = dx*dx/2;
  double alpha = dt/(dx*dx);

  //Matrix u(T,N) is all zeros to start with.
  u.zeros();

  //
  u(0, span::all)    = v;
  u(span::all, 0)    = w;
  u(span::all,N - 1) = z;

  for (int j = 1; j < T - 1; ++j) {
    for (int i = 1; i < N - 1; ++i) {
        u(j+1, i) = alpha*(u(j, i+1) + u(j, i-1)) + (1 - 2*alpha)*u(j, i);
    }
  }
  if (save) { u.save("data_FE.dat", raw_ascii) ;}
  return;
}

void diffusion_BE(mat &u, int N, int T, double dx, rowvec v, vec w, vec z, bool save = 0)
{

  double dt    = dx*dx/2;
  double alpha = dt/(dx*dx);


  u.zeros();

  u(0, span::all)    = v;
  u(span::all, 0)    = w;
  u(span::all,N - 1) = z;

  rowvec x(N);

  for (int j = 1; j < T; ++j) {

    tridiag(x, u(j - 1,span::all), N, -alpha, 1 + 2*alpha);
    u(j,span::all) = x;

  }
  if (save) { u.save("data_BE.dat", raw_ascii) ;}
  return;
}

void diffusion_CN(mat &u, int N, int T, double dx, rowvec v, vec w, vec z, bool save = 0){

    double dt    = dx*dx/2;
    double alpha = dt/(dx*dx);


    u.zeros();

    u(0, span::all)    = v;
    u(span::all, 0)    = w;
    u(span::all,N - 1) = z;


    rowvec x(N);
    for (int j = 1; j < T; j++){
        tridiag(x, u(j - 1,span::all), N, -alpha, 2 + 2*alpha);
        u(j,span::all) = x;

      }
    if (save) {u.save("data_CN.dat", raw_ascii) ;}
    return;
}

