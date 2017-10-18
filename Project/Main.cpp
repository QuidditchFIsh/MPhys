/*
	Author: Aneirin John Baker 
	Date: 23/09/17
	Description: Simulation of a Quantum Harmonic/Anharmonic Oscillator using Hamiltonian Monte Carlo 
		techniques to advance the simulation for my Masters Project.This is the main class where all of 
		the main functions will be called from and the final results collected and outputted from. Here
		I shall be using natural units such that h_bar = c = 1. Beginning with a 1d system. 
*/
#include "Main.h"
#define Fouier_Flip 1
//1-No Fourier Acceleration 
//0-Fourier accelerated algorithm


int main(){

	printf("Beginning Simulation Initalising System\n");
	clock_t t1,t2,t3;

	//Number of iterations of the HMC algorithm to be performed, and number of times the algoirthm is 
	//going to loop

	unsigned int iterations = 10000,length = 1000;
	//unsigned int iterations = 20,length = 10;

	double t_step=0.1;

	//Initalise the vector .	MAYBE IT WOULD BE QUICKER TO USE A MALLOC AND 1D ARRAY
	vector<double> v2(length,0);
	vector<vector<double> >lattice(iterations,v2);

	//0-<q>,1-<q^2>,2-<error>,3-autocorrelation
	vector<double> v1(4,0);
	vector<vector<double> >stats_Data(length,v1);

	//create file to store the data into 

	FILE * output;
	output = fopen("HMC_Results.dat","w");


	//begin each simulation at p=0 and q=0
	printf("Started Simulatio with:\n %d Oscillators\n Iterating %d times at a time step of %f\n",length,iterations,t_step);
	
	t1=clock();
	//#if Fouier_Flip
	lattice_Evolution(lattice,length,t_step,iterations);
	//#endif
	//#if !Fouier_Flip
	//F_lattice_Evolution(lattice,length,t_step,iterations);
	//#endif
	t2=clock();
	
	float seconds =((float)t2-(float)t1)/(CLOCKS_PER_SEC);
	printf("Simulation Completed in %f seconds\n",seconds);
	printf("Begingin Stats Calculations\n");

	//stats calculations go here
//	printf("welp");
	//for(unsigned int i=0;i<iterations/5;i++)
	//{
		//stats_Data[i][1]=avg_X_Sqd(lattice[i]);
		moving_avg_X_Sqd(lattice,stats_Data,iterations,length);
		//fprintf(output,"%d %f\n",i*5,stats_Data[i][1]);
	//}
	//moving_avg_X_Sqd(lattice,stats_Data,iterations,length);
	//autocorrelation_Time(stats_Data,iterations,length);
	//printf("welp");



	t3=clock();
	seconds =((float)t3-(float)t2)/(CLOCKS_PER_SEC);
	printf("Statistics computation Completed in %f seconds\n",seconds);

	#if 1

		//for(unsigned int j=0;j<iterations/5;j++)
			//{
				for(unsigned int i=0;i<length;i++)
				{
				//fprintf(output,"%d "j);
				//fprintf(output,"%f ",lattice[i][j][0]);//p
				fprintf(output,"%f\n",lattice[iterations/5 -1][i]);//q
				//fprintf(output,"%f \n",lattice[0][j][1]);//<x^2>
				//fprintf(output,"%f\n",lattice[i][j][3]);
				}
				//fprintf(output,"\n");
			//}
	#endif

	//after all of the simulation sare completed then perform all of the stastical analysis on it
	





	



}