#include <vector>
#include <stdio.h>
#include <random>
#include <cmath>
#include <algorithm>
using namespace std;

double avgX(vector<double> );

double avg_X_Sqd(vector<double> );

double avg_X_four(vector<double> );

double standard_Deviation(double , double ,double  );

double error_Bars(vector<double> );

double lattice_Hamiltonian(vector<vector<double> > ,unsigned int,double,double);

double lattice_Action(vector<double> ,unsigned int);

double lattice_KineticEnergy(vector<double>,unsigned int);

double Harmonic_hamiltonian(double ,double ,double,double);

double Harmonic_action(double,double);

double Anarmonic_hamiltonian(double ,double ,double,double,double);

double Anarmonic_action(double,double);

double kinetic_Energy(double);