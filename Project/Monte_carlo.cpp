/*
	Author: Aneirin John Baker 
	Date: 23/09/17
	Description: The script is where the HMC algorithm will take place. The intergration methods will be 
		housed here and will be executed here. 
*/
#include "Monte_carlo.h"
#define flip 1


void lattice_Evolution(double t_step,unsigned int iter ,unsigned int length,vector<vector<double> >&lattice)
{
	//initalise the variables
	double p=1,H_Old;
	int acceptance =0;

	//pick out new random P
	default_random_engine generator;
 	normal_distribution<double> distribution(0.0,1.0);

 	vector<double> v1(length,0);
 	vector<vector<double> >store_State(2,v1);

 	for(int i=0;i<length;i++)
 	{
 		//storeing the q's and generating p's
		store_State[i][0]=lattice[i][0];
		store_State[i][1]=distribution(generator);
	}

	//update the lattice
	H_Old = Hamiltonian();

	for(unsigned int i=0;i<length;i++)
	{

		lattice[i][0]= hmcAlgorithm(t_step,p,lattice[i][0]);

	}
	H_New = Hamiltonian();
	if(exp(H_Old - H_New) < 1 )
	{
			//accept the new configuration
		acceptance++;

	}
	else
	{
			//Keep the old one.
		for(unsigned int k =0;k<length;k++)
		{
			lattics[k][j] = store_State[k][0];
		}
	}

	for(unsigned int j=1;j<iter;j++)
	{
		H_Old = Hamiltonian();
		for(unsigned int i=0;i<length;i++)
		{
			p=distribution(generator);
			store_State[i][0]=lattice[i][j-1];

			lattice[i][j]= hmcAlgorithm(t_step,p,lattice[i][j-1]);

		}
		H_New = Hamiltonian();
		if(exp(H_Old - H_New) < 1 )
		{
			//accept the new configuration
			acceptance++;

		}
		else
		{
			//Keep the old one.
			lattics[k][j] = store_State[k];
		}
	}
}
void latticd_Evolution(unsigned int length,unsigned int iterations)
{
	//rewriting algorithm!!!
	//vctor to hold the saved p's and q's. 0-q,1-p
	vector<double> v1(length,0);
	vector<vector<double> >saved_State(2,v1);

	double H_Old=0,H_New=0;
	int acceptance=0;
	//Random Number generator.
	//pick out new random P
	default_random_engine generator;
 	normal_distribution<double> distribution(0.0,1.0);

 	//store the first p's and q's
 	for(unsigned int i = 0; i < length;i++)
 	{
 		//initalise the p's for the state. The first q's will be 0 and evolved randomly by the HMC algorithm.
 		saved_State[i][1] = distribution(generator);
 		saved_State[i][0] = 0;
 	}

 	H_Old = Hamiltonian();
 	for(unsigned int i = 0; i < length ; i++)
 	{
 		lattice[i][0] = hmc();
 	}
 	H_New = Hamiltonian();

 	//Now accept or reject the new state with Metroplis update 
 	if(exp(H_Old - H_New) < 1)
 	{
 		//accept the state 
 		acceptance++;
 	}
 	else
 	{
 		//reject the state and reset the state back to the old one
 		for(unsigned int i = 0; i < length; i++)
 		{
 			//set the lattice to the old state
 			lattics[j][i] = saved_State[i][0]; 
 		}
 	}

 	for(unsigned int i = 1; i < iterations; i++)
 	{
 		for(unsigned int j = 0; j < length;j++)
 		{
 			//initalise the p's for the state. The first q's will be 0 and evolved randomly by the HMC algorithm.
 			saved_State[i][1] = distribution(generator);
 			saved_State[i][j] = lattice[][];
 		}	

 		H_Old = Hamiltonian();
 		for(unsigned int i = 0; i < length ; i++)
 		{
 			lattice[i][j]= hmcAlgorithm(t_step,p,lattice[i][j-1]);	
 		}

 		H_New = Hamiltonian();

 		//Now accept or reject the new state with Metroplis update 
 		if(exp(H_Old - H_New) < 1)
 		{
 			//accept the state 
 			acceptance++;
 		}
 		else
 		{
 			//reject the state and reset the state back to the old one
 			for(unsigned int i = 0; i < length; i++)
 			{
 				//set the lattice to the old state
 				lattice[j][i] = saved_State[i][0];
 			}
 		}
	}	
}
double hmcAlgorithm(double t_step,double p_rand,double q_old)
{
	// HMC algorithm for a single oscillator

	//iter is the number of monte carlo updates which will be perfomred. 
		double q,p=0,p_old;
		unsigned int steps=15;
	//anharmonic algorithm 
#if !flip
		p = p_rand - (0.5 * t_step * (q + (4*q*q*q)));

		for(unsigned int j=0;j<steps;j++)
		{
			q = q + (t_step * p);

			if(j != steps) {p = p - (0.5 * t_step * (q + (4*q*q*q)));}
		}

		p = p - (0.5 * t_step * (q + (4*q*q*q)));

#endif
//harmonic algorithm 
#if flip
		p = p_rand - (0.5 * t_step * q );

		for(unsigned int j=0;j<steps;j++)
		{
			q = q + (t_step * p);

			if(j != steps) {p = p - (0.5 * t_step * q);}
		}

		p = p - (0.5 * t_step * q);

#endif
return q;

}

double hamiltonian(vector<double> q,vector<double> p,unsigned int length,unsigned int iter,double t_step){
	//calculate the Hamiltonian of the whole lattice

	double H=0;

	//Set m=1;
	//loop over entire lattice
	for(unsigned int i = 0;i < length;i++)
	{
		//anharmonic
		#if !flip
			H = (pow((q[iter][i] + q[iter][i+1]),2))/(pow(t_step),2) + (pow(p[iter][i],2) * 0.5) + (pow(q[iter][i],2)* 0.5) + pow(q[iter][i],4);
		#endif

		//harmonic
		#if flip
			H = (pow((q[iter][i] + q[iter][i+1]),2))/(pow(t_step),2) + (pow(p[iter][i],2) * 0.5) + (pow(q[iter][i],2)* 0.5);
		#endif
	}
	return H ;
}



