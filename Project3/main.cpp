#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#define EPS 3.0e-14
#define MAXIT 10
#define ZERO 1.0E-10

#include "Gauss_laguerre.h"              /* Contains: Jacobi() */
#include "Various_functions.h"    /* Contains: Identity() and Toeplitz() */

// Preferences/envirement/system
// Users/kristoffervarslott/Qt/Qt Creator.app/Contents/MacOS/../Resources/scripts/openTerminal.py

using namespace std;
using namespace arma;


//     Here we define various functions called by the main program
/*
double int_function(double x);
void gauss_laguerre(double *, double *, int, double);
double trapezoidal_rule ( double, double, int, double (*func)(double) );
double simpson ( double, double, int, double (*func)(double) );
void gauleg(double, double, double *, double *, int);
double gammln(double);
*/
//   Main function begins here
int main()
{
     int n;
     double a, b;
     cout << "Read in the number of integration points" << endl;
     cin >> n;
     cout << "Read in integration limits" << endl;
     cin >> a >> b;
//   reserve space in memory for vectors containing the mesh points
//   weights and function values for the use of the gauss-legendre
//   method
     double *x = new double [n];
     double *w = new double [n];
     // Gauss-Laguerre is old-fashioned translation of F77 --> C++
     // arrays start at 1 and end at n
     double *xgl = new double [n+1];
     double *wgl = new double [n+1];
     // These arrays are used for improved Gauss-Legendre, mapping of
     // x \in [-1,1] to x \in [0, infinity)
     double *r = new double [n];
     double *s = new double [n];

     for ( int i = 0;  i < n; i++){
        cout << x[i] << endl;
     }

//   set up the mesh points and weights
     gauleg(a, b,x,w, n);
//   set up the mesh points and weights and the power of x^alf
     double alf = 1.0;
     gauss_laguerre(xgl,wgl, n, alf);
//   evaluate the integral with the Gauss-Legendre method
//   Note that we initialize the sum. Here brute force gauleg
     double int_gauss = 0.;
     for ( int i = 0;  i < n; i++){
        int_gauss+=w[i]*int_function(x[i]);
     }
//   evaluate the integral with the Gauss-Laguerre method
//   Note that we initialize the sum
     double int_gausslag = 0.;
     for ( int i = 1;  i <= n; i++){
       int_gausslag += wgl[i];//*sin(xgl[i]);
     }
//   evaluate the integral with the Gauss-Laguerre method
//   Here we change the mesh points with a tangent mapping.
//   Need to call gauleg from -1 to + 1
     gauleg(-1.0, 1.0,x,w, n);
     double pi_4 = acos(-1.0)*0.25;
     for ( int i = 0;  i < n; i++){
       double xx=pi_4*(x[i]+1.0);
       r[i]= tan(xx);
       s[i]=pi_4/(cos(xx)*cos(xx))*w[i];
     }
     double int_gausslegimproved = 0.;
     for ( int i = 0;  i < n; i++){
       int_gausslegimproved += s[i]*int_function(r[i]);
     }

        //cout << int_function(x[i]) << endl;

//    final output
      cout  << setiosflags(ios::showpoint | ios::uppercase);
      cout  << "Trapez-rule = " << setw(20) << setprecision(15) << trapezoidal_rule(a, b,n, &int_function)
           << endl;
      cout << "Simpson's rule = " << setw(20) << setprecision(15) << simpson(a, b,n, &int_function)
           << endl;
      cout << "Gaussian-Legendre quad = "<< setw(20) << setprecision(15)  << int_gauss << endl;
      cout << "Gaussian-Laguerre quad = " << setw(20) << setprecision(15) << int_gausslag << endl;
      cout << "Gaussian-Legendre improved quad = " << setw(20) << setprecision(15) << int_gausslegimproved << endl;

      delete [] x;
      delete [] w;
      delete [] xgl;
      delete [] wgl;
      delete [] s;
      delete [] r;
      return 0;
}  // end of main program



