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

	//Number of iterations of the HMC algorithm to be performed, and number of times the algoirthm is 
	//going to loop
	int iterations = 1,iter=100;

	double q=0.5,p=2,t_step=0.05;



	//Initalise the vector .	MAYBE IT WOULD BE QUICKER TO USE A MALLOC AND 1D ARRAY
	vector<double> v1(iterations,0);
	vector<vector<double> >v2(iter,v1);
	vector<vector<vector<double> > >results(2,v2);

	//create file to store the data into 

	FILE * output;
	output = fopen("HMC_Results.dat","w");



	for(int i=0;i<iterations;i++)
	{
		//begin each simulation at p=0 and q=0

		hmcAlgorithm(t_step,iter,p,q,results);

		//write the data to the stastics array 

	}
	for(unsigned int i=0;i < 1 ;i++)
		for(unsigned int j=0;j<100;j++)
			{
				fprintf(output,"%f ",j*0.05);
				fprintf(output,"%f ",results[i][j][0]);
				fprintf(output,"%f\n",results[i][j][1]);
			}


	//after all of the simulation sare completed then perform all of the stastical analysis on it
	





	



}