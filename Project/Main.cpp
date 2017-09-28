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


	int iterations = 1000,length = 1000;

	double t_step=0.1;




	//Initalise the vector .	MAYBE IT WOULD BE QUICKER TO USE A MALLOC AND 1D ARRAY
	vector<double> v1(2,0);
	vector<vector<double> >v2(length,v1);
	vector<vector<vector<double> > >lattice(iterations,v2);


	//create file to store the data into 

	FILE * output;
	output = fopen("HMC_Results.dat","w");




	//begin each simulation at p=0 and q=0

	lattice_Evolution(t_step,iterations,length,lattice);

	//write the data to the stastics array perhaps!!

	#if 1
	avg_X_Sqd(lattice,iterations,length);

	//for(unsigned int i=0;i < 1 ;i++)
		for(unsigned int j=0;j<length;j++)
			{
				fprintf(output,"%d ",j);
				//fprintf(output,"%f ",lattice[i][j][0]);//p
				//fprintf(output,"%f, ",lattice[499][j][0]);//q
				fprintf(output,"%f \n",lattice[0][j][1]);//<x^2>
				//fprintf(output,"%f\n",lattice[i][j][3]);
			}
	#endif

	//after all of the simulation sare completed then perform all of the stastical analysis on it
	





	



}