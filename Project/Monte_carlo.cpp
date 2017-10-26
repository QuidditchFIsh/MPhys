/*
	Author: Aneirin John Baker 
	Date: 23/09/17
	Description: The script is where the HMC algorithm will take place. The intergration methods will be 
		housed here and will be executed here. 
*/
#include "Monte_carlo.h"
#define Oscillator_flip 1
//1 = Harmonic
//0 = Anharmonic

void lattice_Evolution(vector<vector<double> > &lattice,unsigned int length,double t_step,unsigned int iterations)
{

	FILE * out;
	out = fopen("HMC_LeapFrog_H","w");

	FILE * output_stats;
	output_stats = fopen("HMC_Stats.dat","w");

	FILE * output_X;
	output_X = fopen("HMC_X.dat","w");

	// p-0,q-1
	vector<double> v1(length,0);
	vector<vector<double> > State(2,v1);

	vector<double> v(length,0);
	vector<vector<double> >temp_State(2,v);

	vector<double> v5(length,0);
	vector<vector<double> >Energy_save(3,v);

	vector<double> square_state(length,0);

	vector<double> H_store(201,0);
	H_store[0]=0;

	default_random_engine generator(random_device{}());
 	normal_distribution<double> distribution(0.0,1.0);

 	double acceptance =0,delta_H_Average=0,temp=0,temp1=0,temp2=0,error_x2=0,error_x=0,ACT=0;

 	for(unsigned int i=0;i<length;i++)
 	{
 		State[0][i]=0;
 		State[1][i]=0;
 		temp_State[0][i]=0;
 		temp_State[1][i]=0;
 	}

 	unsigned int result_no=0;
 	//run main algorithm
 	for(unsigned int i = 0; i<iterations;i++)
 	{
 		default_random_engine generator(random_device{}());
 		for(unsigned int j = 0; j<length;j++)
 		{
 			State[0][j] = distribution(generator);
 		}

 #if Oscillator_flip 
 		acceptance += hmcAlgorithm_Harmonic(length,t_step,State,temp_State,H_store);
 #endif

 #if !Oscillator_flip
 		acceptance += hmcAlgorithm_Anharmonic(length,t_step,State,temp_State,H_store);
 #endif
 		
 			//copy results into results array
 			for(unsigned int k = 0;k<length;k++)
 			{
// 				lattice[result_no][k] = State[1][k];
 				square_state[k] = State[1][k] * State[1][k];
 				//printf("%f\n",square_state[k]);
 			}

 			temp1 = avgX(square_state);
 			temp2 = avg_X_Sqd(square_state);
 			error_x = standard_Deviation(temp2,temp1,length);
 			temp1=avgX(State[1]);
 			temp2=avg_X_Sqd(State[1]);
 			error_x2 = standard_Deviation(temp2,temp1,length);

 			ACT = autocorrelation_Time(State[1],State[1],length);

 			fprintf(output_stats,"%d %f %f %f %f %f %f %f %f\n",i,temp1,error_x,temp2,error_x2,lattice_Action(State[1],length),lattice_KineticEnergy(State[0],length),temp,ACT);
 			

 			result_no++;
 		
 	}
 	
 	printf("%f\n",Energy_gap_calc(Energy_save,length));
 	printf("The aacceptance is %f percent \n",(acceptance*100)/(double) iterations);
 	
 	for(unsigned int i=0;i<length;i++)
 	{
 		fprintf(output_X,"%f\n",State[1][i]);
 	}
 	for(unsigned int i=0;i<201;i++)
 	{
 		fprintf(out,"%f \n",H_store[i]/(double)iterations);
 	}
 	

}
double hmcAlgorithm_Harmonic(unsigned int length,double t_step,vector<vector<double> > &old_state,vector<vector<double> > &temp_State,vector<double> &H_store)
{

	double min=0,mu=9;
	unsigned int steps = 10;

	double H_old=0,H_new=0,H_inter=0;

	H_old=lattice_Hamiltonian(old_state,length);

	//half step in the p

	temp_State[0][0] = old_state[0][0] -  (0.5*t_step * ((mu*old_state[1][0]) - (old_state[1][1]+old_state[1][length-1]-(2*old_state[1][0]))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = old_state[0][j] - (0.5*t_step * ((mu*old_state[1][j]) - (old_state[1][j+1]+old_state[1][j-1]-(2*old_state[1][j]))));
		temp_State[1][j] = old_state[1][j];
		//printf("%f %f\n",temp_State[j][0],temp_State[j][1]);
	}
	temp_State[0][length-1] = old_state[0][length-1] - (0.5*t_step * ((mu*old_state[1][length-1]) - (old_state[1][0]+old_state[1][length-2]-(2*old_state[1][length-1]))));




	//full step in p and q for n steps
	for(unsigned int i = 0;i<steps;i++)
	{
		//update all q's
		for(unsigned int j = 0;j<length;j++)
		{
			temp_State[1][j] = temp_State[1][j] + (t_step * temp_State[0][j]);
		}

#if 1
//a full step for when running the algorithm normally
		if(i != steps-1)
		{
		temp_State[0][0] = temp_State[0][0] -  (t_step * ((mu*temp_State[1][0]) - ((temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[0][j] = temp_State[0][j] -  (t_step * ((mu*temp_State[1][j]) - ((temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
		}

		temp_State[0][length-1] = temp_State[0][length-1] - (t_step * ((mu*temp_State[1][length-1]) - ((temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1])))));
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
	temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * ((mu*temp_State[1][0]) - (temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0]))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = temp_State[0][j] - (0.5*t_step * ((mu*temp_State[1][j]) - (temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j]))));
	}
	temp_State[0][length-1] = temp_State[0][length-1] - (0.5*t_step * ((mu*temp_State[1][length-1]) - (temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1]))));
	
	for(unsigned int j=0;j<length;j++)
	{
	//	printf("%f %f\n",temp_State[j][0],temp_State[j][1]);
	}

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
			old_state[1][i] = temp_State[1][i];

		}
		//return H_old - H_new;
		return 1;
	}
	//return H_old - H_new;
	return 0;

}

double hmcAlgorithm_Anharmonic(unsigned int length,double t_step,vector<vector<double> > &old_state,vector<vector<double> > &temp_State,vector<double> &H_store)
{

	double min=0,f=1,mu = -20,lamba=0.1;
	mu = mu * 2;
	lamba = lamba * 4;
	unsigned int steps = 20;

	double H_old=0,H_new=0,H_inter=0;

	H_old=lattice_Hamiltonian(old_state,length);

	//half step in the p

//old potential ie q^4
	temp_State[0][0] = old_state[0][0] -  (0.5*t_step * (lamba*pow(old_state[1][0],3) + (mu* old_state[1][0]) - (old_state[1][1]+old_state[1][length-1]-(2*old_state[1][0]))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = old_state[0][j] - (0.5*t_step * (lamba*pow(old_state[1][j],3) + (mu*old_state[1][j]) - (old_state[1][j+1]+old_state[1][j-1]-(2*old_state[1][j]))));
		temp_State[1][j] = old_state[1][j];
		//printf("%f %f\n",temp_State[j][0],temp_State[j][1]);
	}
	temp_State[0][length-1] = old_state[0][length-1] - (0.5*t_step * (lamba*pow(old_state[1][length-1],3) + (mu *old_state[1][length-1]) - (old_state[1][0]+old_state[1][length-2]-(2*old_state[1][length-1]))));
/*
	//modified potential V = (q^2 - f^2)^2
	temp_State[0][0] = old_state[0][0] -  (0.5*t_step * ((old_state[1][0]*(pow(old_state[1][0],2)-f)) - (old_state[1][1]+old_state[1][length-1]-(2*old_state[1][0]))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = old_state[0][j] - (0.5*t_step * ((old_state[1][j]*(pow(old_state[1][j],2)-f))- (old_state[1][j+1]+old_state[1][j-1]-(2*old_state[1][j]))));
		temp_State[1][j] = old_state[1][j];
		//printf("%f %f\n",temp_State[j][0],temp_State[j][1]);
	}
	temp_State[0][length-1] = old_state[0][length-1] - (0.5*t_step * ((old_state[1][length-1]*(pow(old_state[1][length-1],2)-f))- (old_state[1][0]+old_state[1][length-2]-(2*old_state[1][length-1]))));
*/


	//full step in p and q for n steps
	for(unsigned int i = 0;i<steps;i++)
	{
		//update all q's
		for(unsigned int j = 0;j<length;j++)
		{
			temp_State[1][j] = temp_State[1][j] + (t_step * temp_State[0][j]);
		}

#if 1
//a full step for when running the algorithm normally
		
		if(i != steps-1)

		{
		/*
		//potential as V=(x^2 - f^2)^2
		temp_State[0][0] = temp_State[0][0] -  (t_step * ((temp_State[1][0]*(pow(temp_State[1][0],2)-f)) - ((temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[0][j] = temp_State[0][j] -  (t_step * ((temp_State[1][j]*(pow(temp_State[1][j],2)-f))  - ((temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
		}

		temp_State[0][length-1] = temp_State[0][length-1] - (t_step * ((temp_State[1][length-1]*(pow(temp_State[1][length-1],2)-f))  - ((temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1])))));
		*/
		temp_State[0][0] = temp_State[0][0] -  (t_step * (lamba*pow(temp_State[1][0],3) + (mu*temp_State[1][0]) - ((temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[0][j] = temp_State[0][j] -  (t_step * (lamba*pow(temp_State[1][j],3) + (mu*temp_State[1][j]) - ((temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
		}

		temp_State[0][length-1] = temp_State[0][length-1] - (t_step *(lamba*pow(temp_State[1][length-1],3) + (mu*temp_State[1][length-1]) - ((temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1])))));

		}

#endif 

#if 0
//two half steps for running when checking the algorithm still is equivalent to one full step only with a calcuation of the hamiltonian in the middle
		temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * (pow(temp_State[1][0],3) + temp_State[1][0] - ((temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[0][j] = temp_State[0][j] -  (0.5*t_step * (pow(temp_State[1][j],3) + temp_State[1][j] - ((temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
		}

		temp_State[0][length-1] = temp_State[0][length-1] - (0.5*t_step *(pow(temp_State[1][length-1],3) + (temp_State[1][length-1] - ((temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1]))))));

		//calcuate the hamiltonian here.
		H_store[i+1] +=lattice_Hamiltonian(temp_State,length)-H_old;


		temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * (pow(temp_State[1][0],3) + temp_State[1][0] - ((temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[0][j] = temp_State[0][j] -  (0.5*t_step * (pow(temp_State[1][j],3) + temp_State[1][j] - ((temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
		}

		temp_State[0][length-1] = temp_State[0][length-1] - (0.5 * t_step *(pow(temp_State[1][length-1],3) + t_step * (temp_State[1][length-1] - ((temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1]))))));

#endif

	}
	//half step in the p
	temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * (lamba*pow(temp_State[1][0],3) + (mu *temp_State[1][0]) - (temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0]))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = temp_State[0][j] - (0.5*t_step * (lamba*pow(temp_State[1][j],3) + (mu *temp_State[1][j]) - (temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j]))));
		temp_State[1][j] = temp_State[1][j];
		//printf("%f %f\n",temp_State[j][0],temp_State[j][1]);
	}
	temp_State[0][length-1] = temp_State[0][length-1] - (0.5*t_step * (lamba*pow(temp_State[1][length-1],3) + (mu*temp_State[1][length-1]) - (temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1]))));
//Potential at V=(x^2-f^2)^2
	// temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * ((temp_State[1][0]*(pow(temp_State[1][0],2)-f)) - (temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0]))));
	// for(unsigned int j = 1;j<length-1;j++)
	// {
	// 	temp_State[0][j] = temp_State[0][j] - (0.5*t_step * ((temp_State[1][j]*(pow(temp_State[1][j],2)-f)) - (temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j]))));
	// }
	// temp_State[0][length-1] = temp_State[0][length-1] - (0.5*t_step * ((temp_State[1][length-1]*(pow(temp_State[1][length-1],2)-f)) - (temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1]))));
	
	for(unsigned int j=0;j<length;j++)
	{
	//	printf("%f %f\n",temp_State[j][0],temp_State[j][1]);
	}

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
			old_state[1][i] = temp_State[1][i];

		}
		//return H_old - H_new;
		return 1;
	}
	//return H_old - H_new;
	return 0;


}










