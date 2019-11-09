#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>  //setiosflags, setw, setprecision
#include <string>
//#include <mpi.h>
#include "CoordinateSystem.h"
#include "GaussianQuadrature.h"
#include "MonteCarlo.h"

using namespace std;


int main(){
    // Assigning constant PI and calculating exact value of integral.
    const double PI = atan(1.0)*4;
    int n;
    double a,b;

    int Method;
    cout << "What would you like to do?" << endl;
    cout << "Values listed below give output for one value of 'n' in the terminal\n";
    cout << "[1] Legendre, [2] Laguerre, [3] MonteCarlo\n\n";
    cout << "Remaining options write the output to file up to the given value of 'n'\n";
    cout << "[4] Legendre & Laguerre up to N, [5] MonteCarlo up to... \nMethod: ";
    cin >> Method;

    if (Method == 1){
        // Legendre
        cout << "Read in the number of integration points" << endl;
        cin >> n;
        cout << "-----------------------------------------------------------------------------"<<endl;
        integration_limit();
        cout << "-----------------------------------------------------------------------------"<<endl;
        cout << "Read in integration limits" << endl;
        cin >> a >> b;

        bool WriteToFile = false;
        double int_legendre, Legendre_time;
        GQ_legendre(a, b, n, int_legendre, Legendre_time, WriteToFile);

    }
    else if (Method == 2){
        // Laguerre
        cout << "Read in the number of integration points" << endl;
        cin >> n;

        bool WriteToFile = false;
        double int_laguerre, Laguerre_time;
        GQ_laguerre(n,PI, int_laguerre, Laguerre_time, WriteToFile);

    }
    else if (Method == 3){
        cout << "Do you want to run [1] Brute Force or [2] Importance sampling? " << endl;
        int Number;
        cin >> Number;
        bool WriteToFile = false;

        double Int_MonteCarlo = 0;
        double Time = 0;
        double Standard_Dev = 0;

        switch(Number){
            case 1:
                cout << "Read in the number of integration points" << endl;
                cin >> n;
                cout << "-----------------------------------------------------------------------------"<<endl;
                integration_limit();
                cout << "-----------------------------------------------------------------------------"<<endl;
                cout << "Insert integration limits a and b" << endl;
                cin >> a >> b;
                MonteCarloBF(a,b,n, Int_MonteCarlo, Time, Standard_Dev, WriteToFile);
                break;

            case 2:
                cout << "Read in the number of integration points" << endl;
                cin >> n;
                MonteCarlo_IS(n, Int_MonteCarlo, Time, Standard_Dev, WriteToFile);
                break;
         }
    }
    else if (Method == 4){
        cout << "Read in the number of integration points" << endl;
        cin >> n;
        cout << "-----------------------------------------------------------------------------"<<endl;
        integration_limit();
        cout << "-----------------------------------------------------------------------------"<<endl;
        cout << "Insert integration limits a and b" << endl;
        cin >> a >> b;
        string outfilename = "Project_3_GQ.txt"; // DEN HER MÅ FIKSES
        OutputToFile(a, b, n, outfilename);
    }
    else if (Method == 5){
        cout << "Number of sampling points \n";
        cin >> n;
        cout << "-----------------------------------------------------------------------------"<<endl;
        integration_limit();
        cout << "-----------------------------------------------------------------------------"<<endl;
        cout << "Integration limits\n";
        cin >> a >> b;
        string outfilename = "Project_3_MC.txt"; // Den her må fikses
        OutputToFileMC(a, b, n, outfilename);
    }
}




