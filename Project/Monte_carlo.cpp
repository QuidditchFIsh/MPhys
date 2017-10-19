/*
	Author: Aneirin John Baker 
	Date: 23/09/17
	Description: The script is where the HMC algorithm will take place. The intergration methods will be 
		housed here and will be executed here. 
*/
#include "Monte_carlo.h"
#define Oscillator_flip 1

void lattice_Evolution(vector<vector<double> > &lattice,unsigned int length,double t_step,unsigned int iterations)
{

	FILE * out;
	out = fopen("HMC_LeapFrog_H","w");
	// p-0,q-1
	vector<double> v1(2,0);
	vector<vector<double> > State(length,v1);

	vector<double> v(2,0);
	vector<vector<double> >temp_State(length,v);

	vector<double> H_store(201,0);
	H_store[0]=0;

	default_random_engine generator(random_device{}());
 	normal_distribution<double> distribution(0.0,1.0);

 	double acceptance =0;

 	for(unsigned int i=0;i<length;i++)
 	{
 		State[i][0]=0;
 		State[i][1]=0;
 		temp_State[i][0]=0;
 		temp_State[i][1]=0;
 	}

 	unsigned int result_no=0;
 	//run main algorithm
 	for(unsigned int i = 0; i<iterations;i++)
 	{
 		default_random_engine generator(random_device{}());
 		for(unsigned int j = 0; j<length;j++)
 		{
 			State[j][0] = distribution(generator);
 		}
 		acceptance += hmcAlgorithm(length,t_step,State,temp_State,H_store);

 		if(i % 5 == 0)
 		{
 			//copy results into results array
 			for(unsigned int k = 0;k<length;k++)
 			{
 				lattice[result_no][k] = State[k][1];

 			}
 			result_no++;
 		}	
 	}
 	printf("The Acceptance Rate is %f percent \n",(acceptance*100)/(double) iterations);
 	for(unsigned int i=0;i<201;i++)
 	{
 		fprintf(out,"%f \n",H_store[i]/(double)iterations);
 	}


}
int hmcAlgorithm(unsigned int length,double t_step,vector<vector<double> > &old_state,vector<vector<double> > &temp_State,vector<double> &H_store)
{

	double min=0;
	unsigned int steps = 15;

	double H_old=0,H_new=0,H_inter=0;

	H_old=lattice_Hamiltonian(old_state,length);

	//half step in the p
	//	printf("################\n");
	temp_State[0][0] = old_state[0][0] -  (0.5*t_step * (old_state[0][1] - (old_state[1][1]+old_state[length-1][1]-(2*old_state[0][1]))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[j][0] = old_state[j][0] - (0.5*t_step * (old_state[j][1] - (old_state[j+1][1]+old_state[j-1][1]-(2*old_state[j][1]))));
		temp_State[j][1] = old_state[j][1];
		//printf("%f %f\n",temp_State[j][0],temp_State[j][1]);
	}
	temp_State[length-1][0] = old_state[length-1][0] - (0.5*t_step * (old_state[length-1][1] - (old_state[0][1]+old_state[length-2][1]-(2*old_state[length-1][1]))));
			//printf("-------------------\n");



	//full step in p and q for n steps
	for(unsigned int i = 0;i<steps;i++)
	{
		//update all q's
		for(unsigned int j = 0;j<length;j++)
		{
			temp_State[j][1] = temp_State[j][1] + (t_step * temp_State[j][0]);
		}

#if 1
//a full step for when running the algorithm normally
		if(i != steps-1)
		{
		temp_State[0][0] = temp_State[0][0] -  (t_step * (temp_State[0][1] - ((temp_State[1][1]+temp_State[length-1][1]-(2*temp_State[0][1])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[j][0] = temp_State[j][0] -  (t_step * (temp_State[j][1] - ((temp_State[j+1][1]+temp_State[j-1][1]-(2*temp_State[j][1])))));
		}

		temp_State[length-1][0] = temp_State[length-1][0] - (t_step * (temp_State[length-1][1] - ((temp_State[0][1]+temp_State[length-2][1]-(2*temp_State[length-1][1])))));
		}
#endif 

#if 0
//two half steps for running when checking the algorithm still is equivalent to one full step only with a calcuation of the hamiltonian in the middle
		temp_State[0][0] = temp_State[0][0] -  (0.5 * t_step * (temp_State[0][1] - ((temp_State[1][1]+temp_State[length-1][1]-(2*temp_State[0][1])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[j][0] = temp_State[j][0] -  (0.5 * t_step * (temp_State[j][1] - ((temp_State[j+1][1]+temp_State[j-1][1]-(2*temp_State[j][1])))));
		}

		temp_State[length-1][0] = temp_State[length-1][0] - (0.5 * t_step * (temp_State[length-1][1] - ((temp_State[0][1]+temp_State[length-2][1]-(2*temp_State[length-1][1])))));
		//calcuate the hamiltonian here.
		H_store[i+1] +=lattice_Hamiltonian(temp_State,length)-H_old;


		temp_State[0][0] = temp_State[0][0] -  (0.5 * t_step * (temp_State[0][1] - ((temp_State[1][1]+temp_State[length-1][1]-(2*temp_State[0][1])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[j][0] = temp_State[j][0] -  (0.5 * t_step * (temp_State[j][1] - ((temp_State[j+1][1]+temp_State[j-1][1]-(2*temp_State[j][1])))));
		}

		temp_State[length-1][0] = temp_State[length-1][0] - (0.5 * t_step * (temp_State[length-1][1] - ((temp_State[0][1]+temp_State[length-2][1]-(2*temp_State[length-1][1])))));

#endif



	}
	//half step in the p
	temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * (temp_State[0][1] - (temp_State[1][1]+temp_State[length-1][1]-(2*temp_State[0][1]))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[j][0] = temp_State[j][0] - (0.5*t_step * (temp_State[j][1] - (temp_State[j+1][1]+temp_State[j-1][1]-(2*temp_State[j][1]))));
	}
	temp_State[length-1][0] = temp_State[length-1][0] - (0.5*t_step * (temp_State[length-1][1] - (temp_State[0][1]+temp_State[length-2][1]-(2*temp_State[length-1][1]))));
	for(unsigned int j=0;j<length;j++)
	{
	//	printf("%f %f\n",temp_State[j][0],temp_State[j][1]);
	}
	//printf("################\n");
	
	H_new = lattice_Hamiltonian(temp_State,length);

	//metroplis update
	double r = ((double) rand() / (RAND_MAX));

	min = (1 < exp(H_old - H_new)) ? 1 : exp(H_old - H_new);
	//printf("Hamiltonians: %f\n",exp(H_old-H_new));
	if(r < min)
	{
		//accept
		for(unsigned int i = 0;i<length;i++)
		{
			old_state[i][1] = temp_State[i][1];

		}
		return 1;
	}
	return 0;


}

double lattice_Hamiltonian(vector<vector<double> > state,unsigned int length )
{
	double H=0;
	//loop for all sites which are not effected by periodic BC's
	for(unsigned int i=0;i<length-1;i++)
	{
		H += hamiltonian(state[i][0],state[i][1],state[i+1][1]);
	}
	//Periodic BC sites
	H += hamiltonian(state[length-1][0],state[length-1][1],state[0][1]);

	return H;

}

double hamiltonian(double p,double q,double q_plus )
{
	double h=0;

	h = (p*p*0.5) + (pow((q_plus - q),2)*0.5) + (0.5*q*q);

	return h;
}






