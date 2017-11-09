#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <iostream>
#include <string>
#include <ctime>
#include <vector>
#include <cmath>
#include <algorithm>
#include "stastics.h"

using namespace std;

void lattice_Evolution(vector<vector<double> > &,unsigned int ,double ,unsigned int ,double,double);

double hmcAlgorithm_Harmonic(unsigned int ,double ,vector<vector<double> > &,vector<vector<double> > & ,vector<double> &,double,unsigned int,double);

double hmcAlgorithm_Anharmonic(unsigned int ,double ,vector<vector<double> > &,vector<vector<double> > & ,vector<double> &,double,double,unsigned int,double);

double hmcAlgorithm_Anharmonic_uncoupled(unsigned int ,double ,vector<vector<double> > &,vector<vector<double> > & ,vector<double> &,double,double,unsigned int);



