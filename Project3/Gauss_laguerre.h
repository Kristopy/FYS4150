#ifndef GAUSS_LAGUERRE_H
#define GAUSS_LAGUERRE_H

#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;

void gauss_laguerre(double *x, double *w, int n, double alf);
double gammln( double xx);

#endif // GAUSS_LAGUERRE_H
