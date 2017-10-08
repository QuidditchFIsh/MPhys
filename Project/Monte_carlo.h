#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <iostream>
#include <string>
#include <ctime>
#include <vector>
#include <cmath>

using namespace std;

void lattice_Evolution(vector<vector<double> > &,unsigned int ,double t_step,unsigned int );

void hmcAlgorithm(double ,double ,double ,vector<double> &);

double hamiltonian(vector<vector<double> > ,unsigned int ,double );