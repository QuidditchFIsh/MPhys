/*
	Author: Aneirin John Baker 
	Date: 23/09/17
	Description: The script is where the HMC algorithm will take place. The intergration methods will be 
		housed here and will be executed here. 
*/
#include "Monte_carlo.h"

void lattice_Evolution(double t_step,unsigned int iter ,unsigned int length,vector<vector<double> >&lattice)
{
	//initalise the variables
	double p=1;

	//pick out new random P
	default_random_engine generator;
 	normal_distribution<double> distribution(0.0,1.0);

	//update the lattice

	for(unsigned int i=0;i<length;i++)
	{
		p=distribution(generator);
		//printf("%f\n",p);
		lattice[0][i]= hmcAlgorithm(t_step,p,lattice[0][i]);
		//printf("%f %f \n",p,lattice[0][i][1]);
	}

#if 1
	for(unsigned int j=1;j<iter;j++)
	{
		for(unsigned int i=0;i<length;i++)
		{
			//printf("%d %d \n",j,i);
			p=distribution(generator);
			//printf("%f\n",p);
			lattice[j][i]= hmcAlgorithm(t_step,p,lattice[j-1][i]);
			//printf("%f %f \n",p,lattice[0][i][1]);
		}
	}
#endif


}
double hmcAlgorithm(double t_step,double p_rand,double q_old)
{
	// HMC algorithm executed here

	//iter is the number of monte carlo updates which will be perfomred. 
		double q,p=0,p_old;
		unsigned int steps=25;

		q = q_old;
		p_old = p;



//anharmonic algorithm 
#if 0
		p = p_rand - (0.5 * t_step * (q + (4*q*q*q)));

		for(unsigned int j=0;j<steps;j++)
		{
			q = q + (t_step * p);

			if(j != steps) {p = p - (0.5 * t_step * (q + (4*q*q*q)));}
		}

		p = p - (0.5 * t_step * (q + (4*q*q*q)));

#endif
//harmonic algorithm 
#if 1
		p = p_rand - (0.5 * t_step * q );

		for(unsigned int j=0;j<steps;j++)
		{
			q = q + (t_step * p);

			if(j != steps) {p = p - (0.5 * t_step * q);}
		}

		p = p - (0.5 * t_step * q);

#endif


		//printf("%f\n",exp(hamiltonian(p,q)) - exp(hamiltonian(p_new,q_new)));

		if(1 < exp(-hamiltonian(p_old,q_old)) +(hamiltonian(p,q)))
		{
			//accept the new p and q
			//if rjectted then the p and q stay the same.
			//printf("accept %f\n",q);
			return q;

		}
		else
		{
			//printf("reject %f\n",q_old);
			return q_old;
		}
		
	}


double hamiltonian(double p,double q){
	//calculate the Hamiltonian of the system

	double H=0;

	//Set m=1;
	//anharmonic
#if 0
	H = (pow(p,2) * 0.5) + (pow(q,2)* 0.5) + pow(q,4);
#endif

	//harmonic
#if 1
	H = (pow(p,2) * 0.5) + (pow(q,2)* 0.5);
#endif

	return H ;

}



