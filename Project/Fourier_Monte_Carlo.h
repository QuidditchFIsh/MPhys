#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <iostream>
#include <string>
#include <ctime>
#include <vector>
#include <cmath>

using namespace std;

void F_lattice_Evolution(vector<vector<double> > &,unsigned int ,double t_step,unsigned int );

void F_hmcAlgorithm(double ,double ,double ,double ,double ,vector<double> &);

double F_hamiltonian(vector<vector<double> > ,unsigned int ,double );