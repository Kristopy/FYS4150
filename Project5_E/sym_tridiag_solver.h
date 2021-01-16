#ifndef FUNCTION_H
#define FUNCTION_H

#include <armadillo>

using namespace arma;
using namespace  std;

void sym_tridiag_solver(vec &x, vec y, int N, double b, double d);
#endif
