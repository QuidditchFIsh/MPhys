/*
	Author: Aneirin John Baker 
	Date: 23/09/17
	Description: Simulation of a Quantum Harmonic/Anharmonic Oscillator using Hamiltonian Monte Carlo 
		techniques to advance the simulation for my Masters Project.This is the main class where all of 
		the main functions will be called from and the final results collected and outputted from. Here
		I shall be using natural units such that h_bar = c = 1. Beginning with a 1d system. 
*/
#include "Main.h"

int main(){

	printf("Beginning Simulation Initalising System\n");
	clock_t t1,t2;

	//Number of iterations of the HMC algorithm to be performed, and number of times the algoirthm is 
	//going to loop

	unsigned int iterations = 1000,length = 1000;
	//unsigned int iterations = 20,length = 10;

	double t_step=0.05;

	vector<double> v2(length,0);
	vector<vector<double> >lattice(iterations,v2);

	printf("Started Simulation with:\n %d Oscillators\n Iterating %d times at a time step of %f\n",length,iterations,t_step);
	
	t1=clock();
	lattice_Evolution(lattice,length,t_step,iterations);
	t2=clock();
	
	float seconds =((float)t2-(float)t1)/(CLOCKS_PER_SEC);
	printf("Simulation Completed in %f seconds\n",seconds);

}
