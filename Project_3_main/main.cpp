#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>  //setiosflags, setw, setprecision
#include <string>
//#include <stdio.h>
//#include <stdlib.h>
#define ZERO 1.0E-5
#include "legendre.h"
#include "Laguerre.h"
#include "repulsion_function.h"


using namespace std;
ofstream ofile;

int main(int argc, char* argv[])
{
    char *outfilename;
    outfilename = argv[1];

    const double PI = std::atan(1.0)*4;
    double Exact_value = (5*pow(PI,2))/(256);

    //----------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------------
    int c;
    //Checking for valid integration limits [Lambda]
    double int_limit = 0;
    for (int l = 1;  l < 100; l++){
        c = l;
        for ( int i = 0;  i < c; i++){
            for (int j = 0; j < c; j++){
                for (int k = 0; k < c; k++){
                    int_limit = integration_limit(i,j,k);
                    //cout << "Det her er i: " << i << endl;
                    //cout << "Det her er j: " << j << endl;
                    //cout << "Det her er k: " << k << endl;
                    //cout << "-----------------------" << endl;
                    }
                }
            }
    //If the exponetial value is lower than the defined zero
    //we set the integration limit, as well as breaking the loop
    if (int_limit < ZERO){
        double h = sqrt(c*c + c*c+ c*c);
        cout <<  "Integration limits can be set to:" << "["<<-h<<","<< h<< "]" << endl;
        break;}
    }
    cout << "Value of exponential = "<< setw(20) << setprecision(15)  << int_limit << endl;

    //----------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------------

    //Reading in Number of iterations and integration limits.
    int n;
    double a, b;
    cout << "Read in the number of integration points" << endl;
    cin >> n;
    cout << "Read in integration limits" << endl;
    cin >> a >> b;

    //----------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------------

    //Starting writing in the values we are interested in:
    ofile.open(outfilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "|   N:  |  Legendre:    |   Laguerre:  |  Excact result: |" << endl;
    ofile << "|       |               |              |                 |" << endl;

    //Looping throug different iteration -values of N
    for(int N = 1; N <= n ; N++){

        //----------------------------------------------------------------------------------------------------------------------
        // Gaussian-Legendre quadrature
        double *x = new double [N];
        double *w = new double [N];

        gauleg(a, b, x, w, N);
        double int_gauss = 0.0;
        double iterations = 0.0;
        for ( int i = 0;  i < N; i++){
            for (int j = 0; j < N; j++){
                for (int k = 0; k < N; k++){
                    for (int l = 0; l < N; l++){
                        for (int p = 0; p < N; p++){
                            for (int q = 0; q < N; q++){
                                int_gauss+=w[i]*w[j]*w[k]*w[l]*w[p]*w[q]*Cartesian_Function(x[i],x[j],x[k],x[l],x[p],x[q]);
                                iterations++;
                            }
                        }
                    }
                }
            }
        }

        //----------------------------------------------------------------------------------------------------------------------

        // Gaussian-Laguerre quadrature

        double *R_Gauss_Laguerre = new double [N+1];
        double *Weights_Gauss_Laguerre = new double [N+1];

        double *Theta = new double[N];
        double *Weights_Theta = new double [N];

        double *Phi = new double [N];
        double *Weights_Phi = new double [N];

        double alf = 0.0;

        gauss_laguerre(R_Gauss_Laguerre, Weights_Gauss_Laguerre, N, alf);
        gauleg(0, PI, Theta, Weights_Theta, N);
        gauleg(0, 2*PI, Phi, Weights_Phi, N);

        double int_gausslag = 0.;
        double iterations_Laguerre = 0.0;
        for ( int i = 1;  i <= N; i++){                      // r1
            for (int j = 1; j <= N; j++){                    // r2
                for (int k = 0; k < N; k++){                 // Kun de to første som skal gå fra 1 til n?
                    for (int l = 0; l < N; l++){
                        for (int p = 0; p < N; p++){
                            for (int q = 0; q < N; q++){
                                int_gausslag += Weights_Gauss_Laguerre[i]*Weights_Gauss_Laguerre[j]*Weights_Theta[k]*Weights_Theta[l]*Weights_Phi[p]*Weights_Phi[q]*Spherical_Function(R_Gauss_Laguerre[i],R_Gauss_Laguerre[j],Theta[k], Theta[l], Phi[p], Phi[q]); //Spherical_Coordinates()*Weights for de forskjellige koordinatene//*sin(xgl[i]);
                                iterations_Laguerre++;
                            }
                        }
                    }
                }
            }
        }
        //----------------------------------------------------------------------------------------------------------------------
        /*
         cout << "------------------------------------------------------------------------------ " << endl;
         cout << "Weigths Laguerre" << endl;
         for (int i = 1; i<N-1; i++ ){
            cout << Weights_Gauss_Laguerre[i] << endl;
        }
         cout << "------------------------------------------------------------------------------ " << endl;
         cout << "Weigths Theta" << endl;
        for (int i = 1; i<N-1; i++ ){
            cout << Weights_Theta[i] << endl;
        }
         cout << "------------------------------------------------------------------------------ " << endl;
         cout << "Weigths Phi" << endl;
        for (int i = 1; i<N-1; i++ ){
            cout << Weights_Phi[i] << endl;
        }
         cout << "------------------------------------------------------------------------------ " << endl;
         cout << "Weights Legendre" << endl;
         for (int i = 1; i<N-1; i++ ){
            cout << w[i] << endl;
        }
         cout << "------------------------------------------------------------------------------ " << endl;


        //Limits of [-3, 3] yields 9.59929050878060E-07
        //Limits of [-4, 4] yields 9.40499789667483E-10
        //Limits of [-5, 5] yields 9.21463782719653E-13
        //Limits of [-6, 6] yields 9.02813070446519E-16

        cout << "Iterations: " << iterations << endl;
        cout  << setiosflags(ios::showpoint | ios::uppercase);
        cout << "------------------------------------------------------------------------------ " << endl;
        cout << "Gaussian-Legendre quadrature = "<< setw(20) << setprecision(15)  << int_gauss << endl;
        cout << "Gaussian-Laguerre quadrature = "<< setw(20) << setprecision(15)  << int_gausslag << endl;
        cout << "------------------------------------------------------------------------------ " << endl;
        cout << "Excact result = "<< Exact_value << endl;
        cout << "------------------------------------------------------------------------------ " << endl;
        */

        //----------------------------------------------------------------------------------------------------------------------
        //----------------------------------------------------------------------------------------------------------------------
        // Open file and write results to file:

        //For interpreting where in the loop we are in terminal
        cout << N << endl;
        ofile << setw(5) << N;
        ofile << setw(18) << setprecision(8) << int_gauss;
        ofile << setw(15) << setprecision(8) << int_gausslag;
        ofile << setw(15) << setprecision(8) << Exact_value << endl;


        delete [] x;
        delete [] w;
        delete [] R_Gauss_Laguerre;
        delete [] Weights_Gauss_Laguerre;
        delete [] Theta;
        delete [] Weights_Theta;
        delete [] Phi;
        delete [] Weights_Phi;
    }
    ofile.close();
    //----------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------------


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
