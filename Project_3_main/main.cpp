#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>  //setiosflags, setw, setprecision
#include <string>

#include "repulsion_function.h"
#include "integration_loops.h"

using namespace std;
ofstream ofile;


int main(int argc, char* argv[])
{
    // Assigning constant PI and calculating exact value of integral.
    const double PI = atan(1.0)*4;
    double Exact_value = (5*pow(PI,2))/(256);

    //----------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------------
    cout << "-----------------------------------------------------------------------------"<<endl;
    integration_limit();
    cout << "-----------------------------------------------------------------------------"<<endl;
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

    // Declaration of command line arguments
    if (argc <= 1){
        cout << "Error occured: " << argv[0] <<endl;
        cout <<  "Need do provide a commandline argument for which version to run, [Legendre, Laguerre, MonteCarlo], or make a file by passing inn [Project_3]" << endl;
        exit(1);
    }
    //----------------------------------------------------------------------------------------------------------------------
    // If test for running both Legendre and laguere, and producing file.
    //----------------------------------------------------------------------------------------------------------------------
    else if (argc == 2){
        string outfilename;
        outfilename = argv[1];
        cout << "------------------------------------------------------------------------------ " << endl;
        cout << "Running Legedendre and Laguerre"<< endl;
        cout  << setiosflags(ios::showpoint | ios::uppercase);
        cout << "------------------------------------------------------------------------------ " << endl;
        // Open file and write results to file:
        ofile.open(outfilename);
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        ofile << "|   N:  |  Legendre:    | Relative error Legendre: |   Laguerre:  |  Relative error Laguerre: |  Excact result: |" << endl;
        ofile << "|       |               |                          |              |                           |                 |" << endl;

        //Looping through different iteration -values of N
        for(int N = 1; N <= n ; N++){

            double int_legendre, Legendre_time;
            Gauss_legendre(a, b, N, int_legendre, Legendre_time);

            double int_laguerre,Laguerre_time;
            Gauss_laguerre(N,PI, int_laguerre, Laguerre_time);

            //----------------------------------------------------------------------------------------------------------------------
            //----------------------------------------------------------------------------------------------------------------------
            //Calculating the relative error
            double relative_error_leg = (Exact_value - int_legendre)/Exact_value;
            double relative_error_lag = (Exact_value - int_laguerre)/Exact_value;
            //----------------------------------------------------------------------------------------------------------------------
            //----------------------------------------------------------------------------------------------------------------------
            cout << N << endl;
            //Writing results in file.
            ofile << setw(5) << N;
            ofile << setw(18) << setprecision(8) << int_legendre;
            ofile << setw(20) << setprecision(8) << relative_error_leg;
            ofile << setw(22) << setprecision(8) << int_laguerre;
            ofile << setw(22) << setprecision(8) << relative_error_lag;
            ofile << setw(22) << setprecision(8) << Exact_value << endl;

        }
        //Loop ended
        ofile.close();}
    //----------------------------------------------------------------------------------------------------------------------
    // If test for either Legendre or Laguerre
    //----------------------------------------------------------------------------------------------------------------------
    else if (argc == 3){
        string outfilename;
        string method;
        outfilename = argv[1];
        method = argv[2];
        if (method == "Legendre"){
            double int_legendre, Legendre_time;
            Gauss_legendre(a, b, n,int_legendre, Legendre_time);

            cout << "------------------------------------------------------------------------------ " << endl;
            cout << "Running Gaussian-Legendre quadrature"<< endl;
            cout  << setiosflags(ios::showpoint | ios::uppercase);
            cout << "------------------------------------------------------------------------------ " << endl;
            cout << "Gaussian-Legendre quadrature = "<< setw(20) << setprecision(15)  << int_legendre << endl;
            cout << "------------------------------------------------------------------------------ " << endl;
            cout << "Excact result = "<< Exact_value << endl;
            cout << "------------------------------------------------------------------------------ " << endl;
        }
        else if (method == "Laguerre"){
            double int_laguerre,Laguerre_time;
            Gauss_laguerre(n,PI, int_laguerre, Laguerre_time);

            cout << "------------------------------------------------------------------------------ " << endl;
            cout << "Running Gauss-Laguerre Quadrature" << endl;
            cout << setiosflags(ios::showpoint | ios::uppercase);
            cout << "------------------------------------------------------------------------------ " << endl;
            cout << "Gaussian-Laguerre quadrature = "<< setw(20) << setprecision(15)  << int_laguerre << endl;
            cout << "------------------------------------------------------------------------------ " << endl;
            cout << "Excact result = "<< Exact_value << endl;
            cout << "------------------------------------------------------------------------------ " << endl;
        }
    }

}  // end of main program

