/*
	Author: Aneirin John Baker 
	Date: 23/09/17
	Description: The script is where the HMC algorithm will take place. The intergration methods will be 
		housed here and will be executed here. 
*/
#include <stdlib.h>
#include <stdio.h>

void leapFrog(double p,double q,double t_step,double p_new,double q_new){
	// creation of variables
	double p_New_Half =0;

	// leaf frog method for a harmonic oscaillator so the differentialtion does not need to be 
	//approximated 

	//Leapfrog Method
	p_Half_New = p - (0.5 * t_step * q);

	q_new = q + (t_step * p_Half_New);

	p_new = p_Half_New - (0.5*t_step*q_new);


}

void hmcAlgorithm(double t_step,int iter)
{
	// HMC algorithm executed here
	//pick out new random momentum (CHECK THAT THIS RANDOM NUMBER GENERATOR IS USABLE MAY NEED TO USE A GAUSIAN)

	double p_new,q_new,p,q,p_propose;


	for(int i=0;i<iter;i++)
	{
		p_propose = ((double) rand()/(RAND_MAX));

		//now update the p and q with leapforg

		leapFrog(p,q,t_step,.p_new,q_new);

		if(exp(hamiltonian(p,q) - hamiltonian(p_new,q_new)) < 1)
		{
			//accept the new p and q
			//if rjectted then the p and q stay the same.
			p = p_new;
			q = q_new;

		}
		
	}


}

double hamiltonian(double p,double q){
	//calculate the Hamiltonian of the system

	double H=0;

	H = (p * p * 0.5 *(1/m)) + (0.5 * p * p);

	return H;

}

double gausianRandom(){


}


