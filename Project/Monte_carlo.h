#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <string>


void leapFrog(double p,double q,double t_step,double p_new,double q_new);

void hmcAlgorithm(double t_step,int iter);

double hamiltonian(double p,double q);