#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <iostream>
#include <string>
#include <ctime>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

void lattice_Evolution(vector<vector<double> > &,unsigned int ,double ,unsigned int );

int hmcAlgorithm(unsigned int ,double ,vector<vector<double> > &,vector<vector<double> > & );

double hamiltonian(double ,double ,double ,double );

double lattice_Hamiltonian(vector<vector<double> > ,unsigned int ,double );