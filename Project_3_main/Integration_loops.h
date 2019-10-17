#ifndef INTEGRATION_LOOPS_H
#define INTEGRATION_LOOPS_H

double GQ_legendre(double a, double b, int N,  double& int_legendre, double& Legendre_time);
double GQ_laguerre(int N, double PI, double& int_laguerre, double& Laguerre_time);
double integration_limit_main();

#endif // INTEGRATION_LOOPS_H
