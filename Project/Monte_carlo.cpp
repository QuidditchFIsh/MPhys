/*
	Author: Aneirin John Baker 
	Date: 23/09/17
	Description: The script is where the HMC algorithm will take place. The intergration methods will be 
		housed here and will be executed here. 
*/
#include "Monte_carlo.h"

void leapFrog(double p,double q,double t_step,double p_new,double q_new){
	// creation of variables
	double p_New_Half =0;

	// leaf frog method for a harmonic oscaillator so the differentialtion does not need to be 
	//approximated 

	//Leapfrog Method
	p_New_Half = p - (0.5 * t_step * q);

	q_new = q + (t_step * p_New_Half);

	p_new = p_New_Half - (0.5*t_step*q_new);


}

void hmcAlgorithm(double t_step,int iter,double p,double q)
{
	// HMC algorithm executed here
	//pick out new random momentum (CHECK THAT THIS RANDOM NUMBER GENERATOR IS USABLE MAY NEED TO USE A GAUSIAN)

	double p_new,q_new,p_propose;

	//iter is the number of monte carlo updates which will be perfomred. 
	for(int i=0;i<iter;i++)
	{
		p_propose = ((double) rand()/(RAND_MAX));

		//now update the p and q with leapforg

		leapFrog(p,q,t_step,p_new,q_new);

		if(exp(hamiltonian(p,q) - hamiltonian(p_new,q_new)) < 1)
		{
			//accept the new p and q
			//if rjectted then the p and q stay the same.
			p = p_new;
			q = q_new;

		}
		printf("%d %d\n",p,q );
		
	}


}

double hamiltonian(double p,double q){
	//calculate the Hamiltonian of the system

	double H=0;

	//Set m=1;

	H = (p * p * 0.5 ) + (0.5 * p * p);

	return H;

}

double gausianRandom(){


}


