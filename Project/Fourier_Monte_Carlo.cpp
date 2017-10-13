/*
	Author: Aneirin John Baker 
	Date: 23/09/17
	Description: The script is where the HMC algorithm will take place. The intergration methods will be 
		housed here and will be executed here. 
*/
#include "Fourier_Monte_Carlo.h"
#define Oscillator_flip 1
//1 = harmonic oscillator, 0 = anharmonic oscillator


void F_lattice_Evolution(vector<vector<double> > &results,unsigned int length,double t_step,unsigned int iterations)
{
	//Create two arrays to store the tempory data in. The results will be stored in the results array passed in and recrorded every 5
	// p-0,q-1
	printf("%d %d \n",results.size(),results[0].size());
	vector<double> v1 (2,0);
	vector<vector<double> > old_State(length,v1);

	vector<double> v2 (2,0);
	vector<vector<double> > new_State(length,v2);

	unsigned int acceptance =0,result_no=0;

	double H_new=0,H_old=0;


	//initalise the random number generator as a gausian with mean 0 and std 1
	default_random_engine generator(random_device{}());
 	normal_distribution<double> distribution(0.0,1.0);

 	//PERHAPS LOOK AT HOW EFFICENT IT WOULD BE TO PUT IF STATEMENTS INTO THE FOR LOOP TO IRADICATE THIS 


 	////////////////MAY NOT NEED THIS//////////////////////////
 	for(unsigned int i = 0; i < length;i++)
 	{
 		//create an array of random numbers 
 		old_State[i][0] = 0;
 		old_State[i][1] = 0;
 	}
 	///////////////////////////////////////////////////////////
 	//main loop
 	for(unsigned int i = 0; i < iterations ;i++)
 	{
 		default_random_engine generator(random_device{}());
 		//create new array of random p's 
 		for(unsigned int j = 0; j<length;j++)
 		{
 			new_State[j][0] = distribution(generator);
 		}

 		 H_old =F_hamiltonian(new_State,length,t_step);
 		
 		for(unsigned int j = 0;j < length;j++)
 		{
 			//printf("%f \n",new_State[j][0]);
 			F_hmcAlgorithm(t_step,new_State[j][0],old_State[(j-1) % length][1],old_State[j][1],old_State[(j+1) % length][1],new_State[j]);
 		}
 
 		H_new= F_hamiltonian(new_State,length,t_step);
 		//printf("%f\n",H_old-H_new);

 		if(exp(H_old - H_new) < 1)
 		{
 			acceptance++;
 			for(unsigned int k=0;k<length;k++)
 			{
 				old_State[k][1] = new_State[k][1];
 			}
 			H_old = H_new;
 		}
 		if(i % 5 == 0)
 		{
 			//copy results into results array
 			for(unsigned int k = 0;k<length;k++)
 			{
 				//printf("%d %f %d\n",k,old_State[k][1],result_no);
 				results[result_no][k] = old_State[k][1];
 				//results[result_no][k] = 0;
 
 			}
 			result_no++;
 		}
 	}
 	printf("Acceptance rate is :%f %d/%d\n",(double)acceptance/iterations,acceptance,iterations);





}
void F_hmcAlgorithm(double t_step,double p_rand,double q_minus,double q,double q_plus,vector<double> &state)
{
	// HMC algorithm for a single oscillator

	//iter is the number of monte carlo updates which will be perfomred. 
		double p=0;
		unsigned int steps=15;
	//anharmonic algorithm 
#if !Oscillator_flip
		p = p_rand - (0.5 * t_step * (q + (4*q*q*q)));

		for(unsigned int j=0;j<steps;j++)
		{
			q = q + (t_step * p);

			if(j != steps-1) {p = p - (0.5 * t_step * (q_minus+q_plus-(2*q) + 2*q+ (4*q*q*q)));}
		}

		p = p - (0.5 * t_step * (q + (4*q*q*q)));

#endif
//harmonic algorithm 
#if Oscillator_flip
		p = p_rand - (0.5 * t_step * q );

		for(unsigned int j=0;j<steps;j++)
		{
			q = q + (t_step * p);

			if(j != steps-1) {p = p - (0.5 * t_step * (q_minus+q_plus-(2*q))* 2*q);}
		}

		p = p - (0.5 * t_step * q);

#endif
		state[0] = p;
		state[1] = q;
}

double F_hamiltonian(vector<vector<double> > state,unsigned int length,double t_step)
{
	//calculate the Hamiltonian of the whole lattice

	double H=0;

	//Set m=1;
	//loop over entire lattice
	for(unsigned int i = 0;i < length;i++)
	{
		//anharmonic
		#if !Oscillator_flip

			H += (pow((-state[i][1] + state[(i+1) % length][1]),2))/(0.5*pow(t_step,2)) + (pow(state[i][0],2) * 0.5) + (pow(state[i][1],2)* 0.5) + (pow(state[i][1],4));
		#endif

		//harmonic
		#if Oscillator_flip
			H += (pow((-state[i][1] + state[(i+1) % length][1]),2))/(0.5*pow(t_step,2)) + (pow(state[i][0],2) * 0.5) + (pow(state[i][1],2)* 0.5);
		#endif
	}
	return H ;
}




