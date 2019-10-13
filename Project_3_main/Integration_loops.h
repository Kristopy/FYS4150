#ifndef INTEGRATION_LOOPS_H
#define INTEGRATION_LOOPS_H

double Gauss_legendre(double a, double b, int N,  double& int_legendre, double& Legendre_time);
double Gauss_laguerre(int N, double PI, double& int_laguerre, double& Laguerre_time);
double integration_limit();

#endif // INTEGRATION_LOOPS_H
