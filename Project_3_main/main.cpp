#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>  //setiosflags, setw, setprecision
//#include <stdio.h>
//#include <stdlib.h>
#define ZERO 1.0E-12
#include "legendre.h"
#include "repulsion_function.h"

using namespace std;


int main()
{
    int c;
    //Checking for valid integration limits [Lambda]
    double int_limit = 0.0;
    for (int l = 1;  l < 100; l++){
        c = l;
        for ( int i = 0;  i < c; i++){
            for (int j = 0; j < c; j++){
                for (int k = 0; k < c; k++){
                    int_limit = integration_limit(i,j,k);
                }
            }
        }
        //If the exponetial value is lower than the defined zero
        //we set the integration limit, as well as breaking the loop
        if (int_limit < ZERO){
            cout <<  "Integration limits can be set to:" << "["<<-c<<","<< c<< "]" << endl;
            break;
        }
    }



    int n;
    double a, b;
    cout << "Read in the number of integration points" << endl;
    cin >> n;
    cout << "Read in integration limits" << endl;
    cin >> a >> b;

    double *x = new double [n];
    double *w = new double [n];

    gauleg(a, b, x, w, n);
    double int_gauss = 0.0;
    double iterations = 0.0;
    for ( int i = 0;  i < n; i++){
        for (int j = 0; j < n; j++){
            for (int k = 0; k < n; k++){
                for (int l = 0; l < n; l++){
                    for (int p = 0; p < n; p++){
                        for (int q = 0; q < n; q++){
                            int_gauss+=w[i]*w[j]*w[k]*w[l]*w[p]*w[q]*Cartesian_Function(x[i],x[j],x[k],x[l],x[p],x[q]);
                            iterations++;
                        }
                    }
                }
            }
        }
    }


    //Limits of [-3, 3] yields 9.59929050878060E-07
    //Limits of [-4, 4] yields 9.40499789667483E-10
    //Limits of [-5, 5] yields 9.21463782719653E-13
    //Limits of [-6, 6] yields 9.02813070446519E-16
    const double PI = std::atan(1.0)*4;
    double Exact_value = (5*pow(PI,2))/(256);
    cout << "Iterations: " << iterations << endl;
    cout  << setiosflags(ios::showpoint | ios::uppercase);
    cout << "Gaussian-Legendre quad = "<< setw(20) << setprecision(15)  << int_gauss << endl;
    cout << "Value of exponential = "<< setw(20) << setprecision(15)  << int_limit << endl;
    cout << "Excact result = "<< Exact_value << endl;
}  // end of main program










//   reserve space in memory for vectors containing the mesh points
//   weights and function values for the use of the gauss-legendre
//   method
// Gauss-Laguerre is old-fashioned translation of F77 --> C++
// arrays start at 1 and end at n
//double *xgl = new double [n+1];
//double *wgl = new double [n+1];
// These arrays are used for improved Gauss-Legendre, mapping of
// x \in [-1,1] to x \in [0, infinity)
//double *r = new double [n];
//double *s = new double [n];
//   set up the mesh points and weights
