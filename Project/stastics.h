#include <vector>
#include <stdio.h>
#include <random>
#include <cmath>
using namespace std;

double avgX(vector<double> );

void moving_avg_X_Sqd(vector<vector<double> >,vector<vector<double> > &,unsigned int,unsigned int );

double avg_X_Sqd(vector<double> );

double standard_Deviation(double , double ,double  );

double error_Bars(vector<double> );

void autocorrelation_Time(vector<vector<double> > ,unsigned int ,unsigned int );

double lattice_Hamiltonian(vector<vector<double> > ,unsigned int);

double lattice_Action(vector<double> ,unsigned int);

double lattice_KineticEnergy(vector<double>,unsigned int);

double hamiltonian(double ,double ,double);

double action(double,double);

double kinetic_Energy(double);