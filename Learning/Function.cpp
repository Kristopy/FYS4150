// function example
#include<iostream>
#include<fstream>
#include<iomanip>
#include<math.h>
#include<iomanip>
#include<time.h>
#include<algorithm>
#include<vector>
#include <armadillo>

using namespace arma;
using namespace std;

void matrix (int n)
{
    double h = 2.0;

    mat A = zeros<mat>(n,n);

    //The diagonal elements set to 2.
    for (int i=0;i<n;i++)
    {
        A(i,i) = 2.0;
    }

    //The non-diagonal elements set to -1
    for (int i=1;i<n;i++)
    {
        A(i,i-1) = -1.0/h*h;

        A(i-1,i) = -1.0/h*h;
    }


    mat B = zeros<mat>(n,n);
    B.fill(2);

    cout<< A <<endl;
}

void orthogonality (int n)
{
   mat A = zeros<mat>(n,n);

   for (int i=0;i<n;i++)
   {

       A(i,i) = 1.0;

   }
   cout << A << endl;

}

void analytic_solution (int n)
{
    double h;
    double p_max = 1;
    double p_min = 0;
    double p;
    double lambda;
    int d = 2;
    int a = -1;

    h = (p_max-p_min)/(n);

    for (int i = 0; i <= n; i++)
    {
        p = p_min + i*h;
        //cout << "p values:" << p << endl;
    }


    for (int i = 1; i <= n; i++)
    {
        lambda = d/(h*h) + 2*a/(h*h)*cos((i*M_PI)/(n + 1));
        cout << "eigenvalues: "<< "i = "<< i << " :" <<lambda << endl;
    }
}

void jacobi_max (mat A, int n)
{
   double max;
   for (int i = 0; i < n; ++i)
   {
       //cout << A(i,i) << endl;
       for ( int j = 0; j < n; ++j)
       {
           double aij = fabs(A(i,j));
           if ( aij > max)
           {
              max = aij;
              int k = i;
              int l = j;
              cout << l << k << endl;
           }
       }
   }
}
