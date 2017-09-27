/*
	Author: Aneirin John Baker 
	Date: 23/09/17
	Description: The script is where the HMC algorithm will take place. The intergration methods will be 
		housed here and will be executed here. 
*/
#include "Monte_carlo.h"



void leapFrog(double p,double q,double t_step,double& p_new,double& q_new){
	// creation of variables
	double p_New_Half =0;

	// leaf frog method for a harmonic oscaillator so the differentialtion does not need to be 
	//approximated 

	//Leapfrog Method
	p_New_Half = p - (0.5 * t_step * q);

	q_new = q + (t_step * p_New_Half);

	p_new = p_New_Half - (0.5*t_step*q_new);
	//printf("%f %f %f %f\n",p_new,q_new,p,q);


}

void hmcAlgorithm(double t_step,int iter,double p,double q,vector<vector<vector<double> > > &results)
{
	// HMC algorithm executed here
	//pick out new random momentum (CHECK THAT THIS RANDOM NUMBER GENERATOR IS USABLE MAY NEED TO USE A GAUSIAN)

	double p_new,q_new,p_propose;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0, 1);
#if 1
	leapFrog(p,q,t_step,p_new,q_new);

		//printf("%f\n",exp(hamiltonian(p,q)) - exp(hamiltonian(p_new,q_new)));

	if(exp(hamiltonian(p,q)) - exp(hamiltonian(p_new,q_new)) < 1)
	{ x
		//accept the new p and q
		//if rjectted then the p and q stay the same.
		p = p_new;
		q = q_new;
		//printf("welp");

	}
	results[0][0][0]=p;
	results[0][0][1]=q;
#endif
	//iter is the number of monte carlo updates which will be perfomred. 
	for(unsigned int i=0;i<iter;i++)
	{
		p_propose = dis(gen);

		//now update the p and q with leapforg
		leapFrog(p_propose,q,t_step,p_new,q_new);

		//printf("%f\n",exp(hamiltonian(p,q)) - exp(hamiltonian(p_new,q_new)));

		if(exp(hamiltonian(p_propose,q)) - exp(hamiltonian(p_new,q_new)) < 1)
		{
			//accept the new p and q
			//if rjectted then the p and q stay the same.
			p = p_new;
			q = q_new;
		//	printf("%f\n",q);
		//	printf("welp\n");

		}
		printf("%f %f\n",p,q);

		results[0][i][0]=p;
		results[0][i][1]=q;
	}

}

double hamiltonian(double p,double q){
	//calculate the Hamiltonian of the system

	double H=0;

	//Set m=1;

	H = (p * p ) + (q * q);

	return (H * 0.5);

}



