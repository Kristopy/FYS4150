#include <iostream>
#include <armadillo>
#include <string>

using namespace std;
using namespace arma;

mat Identity(int n){
    mat A = zeros<mat>(n,n);    // Constructs a nxn matrix filled with zeroes
    for (int i = 0; i<n; i++){
        A(i,i) = 1.0;           // Fills the diagonal elements with ones.
    }
    return A;
}

//-----------------------------------------------------------------------------------------------

mat Toeplitz(int n, double rho_max, string method, double wr){
    // Overwriting the value of rho_max given as command line argument
    // for Buckling if it's different from 1.0
    if (rho_max != 1.0 && method == "buckling"){
        cout << "You inserted rho_max =/= 1.0 for this method. Setting it to rho_max = 1.0\n";
        rho_max = 1.0;
    }

    mat A = zeros<mat>(n,n);   // Constructs a nxn matrix filled with zeroes
    vec rho(n+1);              // Declaring a vector with given size

    rho(0) = 0.0;              // Start point
    rho(n) = rho_max;          // End point (Infinity for Qdot1 and Qdot2)

    // Calculating the step size
    double h = (rho(n) - rho(0))/ (double) n;

    // Calculating the values for rho
    for (int i = 1; i < n; i++){
        rho(i) = rho(0) + i*h;
    }

    double Diagonal = 2.0/(h*h);       // Assigning the default value of diagonal
    double NonDiagonal = -1.0/(h*h);   // Assigning the default value of non-diagonal

    if (method == "buckling" ) {
        cout << "Running Buckling" << endl;

        // Indexing the missing values at start and end
        // that won't be assigned in the for loop
        A(0,0) = Diagonal; A(0,1) = NonDiagonal;
        A(n-1,n-1) = Diagonal; A(n-1, n-2) = NonDiagonal;

        for (int i = 1; i < n-1; i++){
            A(i,i) = Diagonal;
            A(i,i-1) = NonDiagonal;
            A(i,i+1) = NonDiagonal;
        }
    }
    else if(method == "qdot1") {
        cout << "Running Qdot1" << endl;

        // Indexing the missing values at start and end
        // that won't be assigned in the for loop
        A(0,0) = Diagonal + rho(1)*rho(1); A(0,1) = NonDiagonal;
        A(n-1,n-1) = Diagonal + rho(n)*rho(n); A(n-1, n-2) = NonDiagonal;

        for (int i = 1; i < n-1; i++){
            A(i,i) = Diagonal + rho(i+1)*rho(i+1);   // Assigning different value to the diagonal element
            A(i,i-1) = NonDiagonal;
            A(i,i+1) = NonDiagonal;
        }
    }

    else if (method == "qdot2"){
        cout << "Running Qdot2 \n";
        if (wr == 0.0){
            cout << "Can't run Qdot2 with frequency wr = 0.0\n";
            exit(1);
        }
        cout << "Frequency wr " << wr << endl;

        // Frequency variabel
        double wr_2 = wr*wr;

        //Indexing the missing values at start and end
        //Since Diagonal contains one more value
        A(0,0) = Diagonal + wr_2*rho(1)*rho(1)+1.0/rho(1); A(0,1) = NonDiagonal;
        A(n-1,n-1) = Diagonal + wr_2*rho(n)*rho(n) + 1.0/rho(n); A(n-1, n-2) = NonDiagonal;

        for (int i = 1; i < n-1; i++){
            A(i,i) = Diagonal + wr_2*rho(i+1)*rho(i+1) + 1.0/rho(i+1);  // Assigning different value to the diagonal element
            A(i,i-1) = NonDiagonal;
            A(i,i+1) = NonDiagonal;

        }
     }
    else {
        cout << "-------------------------------------------------------------------\n";
        cout << "Error: Third command line argument is required to choose the method used." << endl;
        cout << "-------------------------------------------------------------------\n";
     exit(1);
     }
    return A;    // Returning matrix A
}
