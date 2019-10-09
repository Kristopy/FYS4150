#ifndef VARIOUS_FUNCTIONS_H
#define VARIOUS_FUNCTIONS_H
#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;

double int_function(double x);
void gauleg(double x1, double x2, double x[], double w[], int n);
double trapezoidal_rule(double a, double b, int n, double (*func)(double));
double simpson(double a, double b, int n, double (*func)(double));

#endif // VARIOUS_FUNCTIONS_H
